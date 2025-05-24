//
// Created by XTAD on 2025/5/29.
//

#ifndef AHRS_H
#define AHRS_H

#include <vector>
#include <fmt/base.h>
#include "common/angle.h"
#include "common/types.h"
#include "Eigen/Geometry"

using Eigen::Matrix3d;
using Eigen::MatrixXd;
using Eigen::Vector3d;
using Eigen::VectorXd;

class AHRS {
private:
    Matrix3d cbh; /* 旋转矩阵: b -> h */
    Matrix3d chn; /* 旋转矩阵: h -> n */
    Matrix3d cbn; /* 旋转矩阵: b -> n */

    double           cur_gravity{ 9.81 };                   /* 初始的重力加速度值 */
    Vector3d         gyro_bias_estimate = Vector3d::Zero(); /* 陀螺仪偏差的积分估计 */
    double           Kp{ 1 };                               /*  比例增益 (可调) */
    constexpr double Ki{ 0.001 };                           /* 积分增益 (可调) */
    bool             alignment_ready{ false };              /* 是否已经完成对准 */

public:
    /* === 误差角对准 === */
    /* ================================================================================================================= */
    /**
     * @brief 水平面对准
     * @param imu_data 范围内IMU数据
     * @return 对准是否成功
     */
    bool HorizontalCorrect(const std::vector<IMU>& imu_data) {
        using namespace Eigen;
        double sum_acc_x = 0.0, sum_acc_y = 0.0, sum_acc_z = 0.0;
        for (const auto& data : imu_data) {
            sum_acc_x += data.accel_raw.x();
            sum_acc_y += data.accel_raw.y();
            sum_acc_z += data.accel_raw.z();
        }
        const double avg_acc_x = sum_acc_x / imu_data.size();
        const double avg_acc_y = sum_acc_y / imu_data.size();
        const double avg_acc_z = sum_acc_z / imu_data.size();

        cur_gravity = sqrt(avg_acc_x * avg_acc_x + avg_acc_y * avg_acc_y + avg_acc_z * avg_acc_z);

        const double     roll  = asin(avg_acc_y / cur_gravity);
        const double     pitch = asin(avg_acc_x / cur_gravity);
        constexpr double yaw   = 0.0;

        const Matrix3d rz = AngleAxisd(yaw, Vector3d::UnitZ()).toRotationMatrix();
        const Matrix3d ry = AngleAxisd(pitch, Vector3d::UnitY()).toRotationMatrix();
        const Matrix3d rx = AngleAxisd(roll, Vector3d::UnitX()).toRotationMatrix();
        cbh               = rx * ry * rz;

        if (avg_acc_z > 0) {
            fmt::println("avg_acc_z > 0, rotate coordinate system");
            const Matrix3d rzz = AngleAxisd(0, Vector3d::UnitZ()).toRotationMatrix();
            const Matrix3d ryy = AngleAxisd(0, Vector3d::UnitY()).toRotationMatrix();
            const Matrix3d rxx = AngleAxisd(180 * D2R, Vector3d::UnitX()).toRotationMatrix();
            cbh                = rxx * ryy * rzz;
        }
        return true;
    }

    /**
     * @brief 航向角对准
     * @param imu_data 范围内IMU数据
     * @param last_point 上一个GNSS信息
     * @param cur_point 当前GNSS信息
     * @return 航向角是否校准成功
     */
    bool HeadingCorrect(const std::vector<IMU>& imu_data,
                        const GNSS&             last_point, const GNSS& cur_point) {
        using namespace Eigen;
        if (alignment_ready) return true;

        double sum_acc_x = 0.0, sum_acc_y = 0.0, sum_acc_z = 0.0;
        for (const auto& data : imu_data) {
            sum_acc_x += data.accel_raw.x();
            sum_acc_y += data.accel_raw.y();
            sum_acc_z += data.accel_raw.z();
        }
        const double avg_acc_x = sum_acc_x / imu_data.size();
        const double avg_acc_y = sum_acc_y / imu_data.size();
        const double avg_acc_z = sum_acc_z / imu_data.size();

        Vector3d acc{ avg_acc_x, avg_acc_y, avg_acc_z };
        acc = cbh * acc;
        if (acc.x() <= 1 && acc.y() <= 1) return false;

        const double accel_n = (cur_point.vel - last_point.vel) * cos(last_point.heading * D2R) / 1;
        const double accel_e = (cur_point.vel - last_point.vel) * sin(last_point.heading * D2R) / 1;
        const double sin_phi = (-accel_n * acc.y() + accel_e * acc.x()) / (pow(acc.x(), 2) + pow(acc.y(), 2));
        const double cos_phi = (accel_n * acc.x() + accel_e * acc.y()) / (pow(acc.x(), 2) + pow(acc.y(), 2));

        chn << cos_phi, -sin_phi, 0,
                sin_phi, cos_phi, 0,
                0, 0, 1;
        cbn             = chn * cbh;
        alignment_ready = true;
        return true;
    }

    /* ================================================================================================================= */

    /* === Mahony滤波算法 === */
    /* ================================================================================================================= */
    /**
     * @brief Mahony互补滤波
     * @param cbn_current 当前姿态 C_b^n(t_k)
     * @param acc_b_raw 当前加速度计读数 (m/s^2)
     * @param gyro_b_raw 当前陀螺仪读数 (rad/s)
     * @param dt 时间间隔 (s)
     */
    Vector3d ComplementaryFilter(Matrix3d&       cbn_current,
                                 const Vector3d& acc_b_raw,
                                 const Vector3d& gyro_b_raw,
                                 const double    dt) {
        /**
         * 1. 加速度数据预处理 (归一化)
         */
        Vector3d     acc_b_normalized = acc_b_raw;
        const double acc_norm         = acc_b_raw.norm();
        if (acc_norm > 0.01) {
            // 避免除以零，并只在有显著加速度时进行校正
            acc_b_normalized /= acc_norm;
        } else {
            // 如果加速度接近零，则不使用加速度计进行校正（或根据情况处理）
            // 这里简单处理，后续陀螺仪积分仍会进行
            // 或者可以认为此时重力矢量未知，跳过误差计算和校正角速度部分
            // return; // 或仅用陀螺仪更新
        }

        /**
         * 2. 动态更新比例增益项 Kp
         *    不考虑车辆垂直起飞的情况, 此时合外力一定朝下
         */
        Kp = acc_norm / cur_gravity;
        if (Kp > 1) {
            Kp = 1 / Kp;
        }

        /**
         * 3. 从当前姿态 cbn_current 估计体坐标系下的重力矢量方向
         *    假设导航坐标系是NED，重力矢量 g_n = [0, 0, 1]^T (归一化后)
         *    加速度计测量得到比力矢量应为 a_n = [0, 0, -1]^T (归一化后)
         */
        const Vector3d gravity_b_estimated = cbn_current.transpose() * Vector3d(0, 0, -1.0);

        /**
         * 4. 计算加速度测量值与估计重力方向之间的误差 (叉积)
         *    这个误差向量在体坐标系下, 表示了使估计重力转向测量重力的旋转轴和幅度信息
         */
        const Vector3d error_body = acc_b_normalized.cross(gravity_b_estimated);

        /**
         * 5. 更新陀螺仪偏差估计（积分项）
         */
        gyro_bias_estimate += Ki * error_body * dt;
        // fmt::println("---------------------------------------------------------------------------------------------------");
        // fmt::println("acc_b_normalized [{}, {}, {}]^T", acc_b_normalized.x(), acc_b_normalized.y(), acc_b_normalized.z());
        // fmt::println("gravity_b_estimated [{}, {}, {}]^T", gravity_b_estimated.x(), gravity_b_estimated.y(),
        //      gravity_b_estimated.z());
        // fmt::println("error_body [{}, {}, {}]^T", error_body.x(), error_body.y(), error_body.z());
        // fmt::println("gyro_bias_estimate [{}, {}, {}]^T", gyro_bias_estimate.x(), gyro_bias_estimate.y(),
        //      gyro_bias_estimate.z());
        // fmt::println("---------------------------------------------------------------------------------------------------");

        /**
         * 6. 计算校正后的陀螺仪角速度
         *    首先减去估计的偏差，然后加上比例反馈项
         */
        const Vector3d gyro_b_corrected_for_bias = gyro_b_raw - gyro_bias_estimate;
        const Vector3d gyro_b_fused              = gyro_b_corrected_for_bias + Kp * error_body;

        /**
         * 7. 使用融合后的角速度更新 cbn
         */
        const Vector3d delta_alpha = gyro_b_fused * dt;
        Matrix3d       dR          = Matrix3d::Identity();
        const double   angle       = delta_alpha.norm();
        if (angle > 1e-9) {
            // 避免除以零
            const Vector3d axis = delta_alpha / angle;
            dR                  = Eigen::AngleAxisd(angle, axis).toRotationMatrix();
        }
        cbn_current = cbn_current * dR;

        /**
         * 8. 正交化 cbn_current
         */
        const Eigen::JacobiSVD svd(cbn_current, Eigen::ComputeFullU | Eigen::ComputeFullV);
        cbn_current = svd.matrixU() * svd.matrixV().transpose();
        return gyro_b_fused;
    }

    /* ================================================================================================================= */

    /* === AHRS姿态航向角解算 === */
    /* ================================================================================================================= */
    /**
     * @brief 未进行校准误差角的预测
     * @param imu_data 范围内IMU数据
     * @param last_point 上一个GNSS数据
     * @return 预测结果
     */
    GNSS PredictWithoutAlignment(const std::vector<IMU>& imu_data,
                                 const GNSS&             last_point) {
        using namespace Eigen;
        constexpr double dt      = 0.02;
        double           lng     = last_point.llh.x();
        double           lat     = last_point.llh.y();
        double           alt     = last_point.llh.z();
        double           speed   = last_point.vel;
        double           heading = last_point.heading;

        for (auto iter = imu_data.begin(); iter != imu_data.end(); ++iter) {
            Vector3d gyro_b{ iter->gyro_raw.x(), iter->gyro_raw.y(), iter->gyro_raw.z() };
            Vector3d acc_b{ iter->accel_raw.x(), iter->accel_raw.y(), iter->accel_raw.z() };

            if (alignment_ready) {
                gyro_b = ComplementaryFilter(cbn, acc_b, gyro_b, dt);
            }
            // if (alignment_ready) {
            //     const Vector3d delta_alpha = gyro_b * dt;
            //     Matrix3d       dR          = Matrix3d::Identity();
            //     const double   angle       = delta_alpha.norm();
            //     if (angle > 1e-9) {
            //         // 避免除以零
            //         const Vector3d axis = delta_alpha / angle;
            //         dR                  = AngleAxisd(angle, axis).toRotationMatrix();
            //     }
            //     cbn = cbn * dR;
            // }

            const double cur_vt = iter->wheel_speed;
            const double vx     = cur_vt * sin(heading * D2R);
            const double vy     = cur_vt * cos(heading * D2R);
            const double x_dis  = vx * dt;
            const double y_dis  = vy * dt;

            lng += x_dis / (111000.0 * cos(lat));
            lat += y_dis / 111000.0;

            Vector3d actual_gyro;
            if (alignment_ready)
                actual_gyro = cbn * gyro_b;
            else
                actual_gyro = cbh * gyro_b;

            heading += actual_gyro.z() * R2D * dt;

            if (heading < 0) heading += 360;
            if (heading > 360) heading -= 360;
        }

        speed = imu_data.back().wheel_speed;
        // LOGI("vx = {} m/s, vy = {} m/s and speed = {} km/h", vx, vy, speed * 3.6);

        return GNSS{ last_point.time + 1000, { lng, lat, alt }, speed, heading, true };
    }

    /**
     * @brief 已进行校准误差角后的预测
     * @param imu_data 范围内IMU数据
     * @param last_point 上一次GNSS信息
     * @return 预测后结果
     */
    GNSS PredictWithAlignment(const std::vector<IMU>& imu_data,
                              const GNSS&             last_point) {
        using namespace Eigen;
        constexpr double dt  = 0.02;
        double           lng = last_point.llh.x();
        double           lat = last_point.llh.y();
        double           alt = last_point.llh.z();
        // double           speed = last_point.vel;
        double speed   = imu_data.front().wheel_speed;
        double heading = last_point.heading;

        if (speed < 0.2) {
            const double g2 = pow(imu_data.begin()->accel_raw.x(), 2) +
                              pow(imu_data.begin()->accel_raw.y(), 2) +
                              pow(imu_data.begin()->accel_raw.z(), 2);
            cur_gravity = sqrt(g2);
        }

        double vn = speed * cos(heading * D2R);
        double ve = speed * sin(heading * D2R);

        for (auto iter = imu_data.begin(); iter != imu_data.end(); ++iter) {
            Vector3d acc_b{ iter->accel_raw.x(), iter->accel_raw.y(), iter->accel_raw.z() };
            Vector3d gyro_b{ iter->gyro_raw.x(), iter->gyro_raw.y(), iter->gyro_raw.z() };

            gyro_b = ComplementaryFilter(cbn, acc_b, gyro_b, dt);
            // if (alignment_ready) {
            //     const Vector3d delta_alpha = gyro_b * dt;
            //     Matrix3d       dR          = Matrix3d::Identity();
            //     const double   angle       = delta_alpha.norm();
            //     if (angle > 1e-9) {
            //         // 避免除以零
            //         const Vector3d axis = delta_alpha / angle;
            //         dR                  = AngleAxisd(angle, axis).toRotationMatrix();
            //     }
            //     cbn = cbn * dR;
            // }

            Vector3d acc_n  = cbn * acc_b;
            Vector3d gyro_n = cbn * gyro_b;

            const double n_dis = vn * dt + 0.5 * acc_n.x() * pow(dt, 2);
            const double e_dis = ve * dt + 0.5 * acc_n.y() * pow(dt, 2);

            vn += acc_n.x() * dt;
            ve += acc_n.y() * dt;

            lng += e_dis / (111000.0 * cos(lat));
            lat += n_dis / 111000.0;

            heading += gyro_n.z() * R2D * dt;

            if (heading < 0) heading += 360;
            if (heading > 360) heading -= 360;
        }

        speed = sqrt(ve * ve + vn * vn);
        fmt::println("ve = {} m/s, vn = {} m/s and speed = {} km/h", ve, vn, speed * 3.6);

        return GNSS{ last_point.time + 1000, { lng, lat, alt }, speed, heading, true };
    }

    GNSS PredictCoord(const std::vector<IMU>& imu_data,
                      const GNSS&             last_point) {
        if (alignment_ready) return PredictWithAlignment(imu_data, last_point);
        return PredictWithoutAlignment(imu_data, last_point);
    }

    /* ================================================================================================================= */
};

#endif //AHRS_H
