//
// Created by XTAD on 2025/5/23.
//

#ifndef LOOSELY_COUPLED_EKF_H
#define LOOSELY_COUPLED_EKF_H

#include "Eigen/Core"
#include "common/rotation.h"

using Eigen::Vector2d;
using Eigen::Vector3d;
using Eigen::Matrix3d;
using Eigen::VectorXd;
using Eigen::MatrixXd;

// --- 常量定义 ---
const double DEG_TO_RAD            = M_PI / 180.0;     // 角度转弧度的转换因子
const double RAD_TO_DEG            = 180.0 / M_PI;     // 弧度转角度的转换因子
const double EARTH_RADIUS_EQ       = 6378137.0;        // WGS84标准下的地球赤道半径 (米)
const double EARTH_ECCENTRICITY_SQ = 0.00669437999014; // WGS84标准下的地球第一偏心率的平方
const double OMEGA_IE              = 7.292115e-5;      // 地球自转角速率 (弧度/秒)

class GnssImsEskf {
private:
    VectorXd dx_;     /* 误差状态向量 */
    MatrixXd P_;      /* 协方差矩阵 */
    MatrixXd Q_cont_; /* 连续时间过程噪声功率谱密度矩阵 */
    MatrixXd F_c_;    /* 连续时间状态转移矩阵 (雅可比矩阵) */
    MatrixXd R_gnss_; /* GNSS测量噪声协方差矩阵 */

    // 工具函数: 计算给定纬度处的子午圈曲率半径 (R_M) 和卯酉圈曲率半径 (R_N)
    Vector3d radii(double lat_rad) {
        // lat_rad 为纬度，单位弧度
        double sin_lat = sin(lat_rad);
        double den_sq  = 1.0 - EARTH_ECCENTRICITY_SQ * sin_lat * sin_lat; // (1 - e^2 * sin^2(lat))
        if (den_sq < 1e-6) den_sq = 1e-6;                                 // 避免除数为零或过小
        double den = sqrt(den_sq);                                        // sqrt(1 - e^2 * sin^2(lat))

        // R_M = a * (1 - e^2) / (1 - e^2 * sin^2(lat))^(3/2)
        double R_M = EARTH_RADIUS_EQ * (1.0 - EARTH_ECCENTRICITY_SQ) / (den_sq * den); // 子午圈曲率半径
        // R_N = a / sqrt(1 - e^2 * sin^2(lat))
        double R_N = EARTH_RADIUS_EQ / den;    // 卯酉圈曲率半径
        return Eigen::Vector3d(R_M, R_N, 0.0); // Z分量在此处未使用，设为0
    }

public:
    /** --- 误差状态向量 dx_ (11x1) ---
     * 0: dp_N (北向位置误差, 米)
     * 1: dp_E (东向位置误差, 米)
     * 2: dv_N (北向速度误差, 米/秒)
     * 3: dv_E (东向速度误差, 米/秒)
     * 4: dPsi_D (航向平台误差，绕D轴，即NED坐标系下的Z轴, 弧度)
     * 5: dBa_x (x轴加速度计零偏误差, 米/秒^2)
     * 6: dBa_y (y轴加速度计零偏误差, 米/秒^2)
     * 7: dBa_z (z轴加速度计零偏误差, 米/秒^2)
     * 8: dBg_x (x轴陀螺仪零偏误差, 弧度/秒)
     * 9: dBg_y (y轴陀螺仪零偏误差, 弧度/秒)
     * 10: dBg_z (z轴陀螺仪零偏误差, 弧度/秒)
     */
    GnssImsEskf() : dx_(11),
                    P_(11, 11),
                    Q_cont_(11, 11),
                    F_c_(11, 11),
                    R_gnss_(4, 4) {
        dx_.setZero();    // 初始化误差状态向量为零
        P_.setIdentity(); // 初始化协方差矩阵为单位阵
        P_ *= 0.1;        // 设置一个初始的协方差值 (示例值, 需要根据实际情况调整)

        /** --- 设置连续时间过程噪声的功率谱密度矩阵 Q_cont_ ---
         * 这些是功率谱密度值, 需要非常仔细地调整!
         * 功率谱密度的单位:
         * 速度噪声: (米/秒^2)^2 / Hz (来自加速度计噪声)
         * 姿态噪声: (弧度/秒)^2 / Hz (来自陀螺仪噪声)
         * 零偏随机游走: (米/秒^2)^2/秒 / Hz 或 (弧度/秒)^2/秒 /Hz (更准确地说是 (米/秒^2)^2*秒 或 (弧度/秒)^2*秒 对于积分后的零偏)
         * EKF中的Q矩阵是 Q_k = integral(F*G*Q_c*G'*F')dt，对于离散白噪声近似为 G*Q_c*G'*dt
         * 或者简单地 Q_k_diag = sigma_psd_value * dt
         */
        Q_cont_.setZero();                                 // 初始化为零矩阵
        const double pos_process_noise_psd = 0.0;          // 位置过程噪声功率谱密度 (通常位置误差由速度误差驱动，可设为0)
        const double vel_process_noise_psd = pow(1e-2, 2); // 速度过程噪声功率谱密度 ((米/秒^2)/sqrt(Hz) -> (米/秒^2)^2/Hz，对应加速度计噪声)
        // 航向姿态过程噪声功率谱密度 ((弧度/秒)/sqrt(Hz) -> (弧度/秒)^2/Hz，对应陀螺仪噪声)
        const double att_psi_process_noise_psd = pow(1e-4 * DEG_TO_RAD, 2);
        // 加速度计零偏随机游走功率谱密度 ((米/秒^3)/sqrt(Hz) 或 (米/秒^2)/sqrt(秒) -> (米/秒^2)^2/秒)
        const double accel_bias_psd = pow(1e-5, 2);
        // 陀螺仪零偏随机游走功率谱密度 ((弧度/秒^2)/sqrt(Hz) 或 (弧度/秒)/sqrt(秒) -> (弧度/秒)^2/秒)
        const double gyro_bias_psd = pow(1e-6 * DEG_TO_RAD, 2);

        // 对于 dp_N, dp_E (如果速度误差能恰当模型化位置不确定性，通常设为0)
        Q_cont_(0, 0) = pos_process_noise_psd;
        Q_cont_(1, 1) = pos_process_noise_psd;

        // 对于 dv_N, dv_E (由加速度计噪声驱动)
        Q_cont_(2, 2) = vel_process_noise_psd;
        Q_cont_(3, 3) = vel_process_noise_psd;

        // 对于 dPsi_D (由Z轴陀螺仪噪声驱动)
        Q_cont_(4, 4) = att_psi_process_noise_psd;

        // 对于加速度计零偏 (随机游走模型)
        Q_cont_(5, 5) = accel_bias_psd;
        Q_cont_(6, 6) = accel_bias_psd;
        Q_cont_(7, 7) = accel_bias_psd;

        // 对于陀螺仪零偏 (随机游走模型)
        Q_cont_(8, 8)   = gyro_bias_psd;
        Q_cont_(9, 9)   = gyro_bias_psd;
        Q_cont_(10, 10) = gyro_bias_psd;

        /** --- GNSS测量噪声协方差矩阵 R_gnss_ (4x4) ---
         * 顺序: N (米), E (米), Heading (弧度), Ground Speed (米/秒)
         * 需要根据您的GNSS接收机规格进行调整
         */
        R_gnss_       = Eigen::Matrix<double, 4, 4>::Identity(); // 初始化为单位阵
        R_gnss_(0, 0) = pow(2.0, 2);                             // 北向位置方差 (米^2)
        R_gnss_(1, 1) = pow(2.0, 2);                             // 东向位置方差 (米^2)
        R_gnss_(2, 2) = pow(0.5 * DEG_TO_RAD, 2);                // 航向角方差 (弧度^2)
        R_gnss_(3, 3) = pow(0.5, 2);                             // 地速方差 (米/秒)^2
    }

    // EKF 初始化函数
    void initialize(const MatrixXd& initial_P) {
        P_ = initial_P; // 设置初始协方差矩阵
        dx_.setZero();  // 误差状态从零开始
    }

    // EKF 预测步骤
    void predict(const Vector3d& accel_b,              /* b系下的比力测量值 (米/秒^2) */
                 const Vector3d& gyro_b,               /* b系下的角速度测量值 (弧度/秒) */
                 const Matrix3d& C_b_n,                /* 从b系到n系的方向余弦矩阵 (DCM) */
                 const Vector3d& current_vel_ned,      /* 当前n系下的速度 (米/秒) */
                 const Vector3d& current_pos_llh,      /* 当前位置 (纬度, 经度, 高程) (弧度, 弧度, 米) */
                 const Vector3d& current_accel_bias_b, /* 当前加速度计零偏估计值 (b系) */
                 const Vector3d& current_gyro_bias_b,  /* 当前陀螺仪零偏估计值 (b系) */
                 const double    dt) /* 时间间隔 (秒) */ {
        const double Lat   = current_pos_llh.x(); // 当前纬度 (弧度)
        const double h     = current_pos_llh.z(); // 当前高程 (米)
        Vector3d     R_vec = radii(Lat);          // 计算曲率半径
        const double R_M   = R_vec.x();           // 子午圈曲率半径
        const double R_N   = R_vec.y();           // 卯酉圈曲率半径

        // 校正后的比力 (减去零偏)
        const Vector3d f_b_corrected = accel_b - current_accel_bias_b;
        Vector3d       f_n           = C_b_n * f_b_corrected; // n系下的比力

        /** --- 构建连续时间状态转移矩阵 F_c_ ---
         * 详细推导请参考惯性导航误差方程 (例如：秦永元《惯性导航》或 Groves《Principles of GNSS, Inertial, and Multisensor Integrated Navigation Systems》)
         */
        F_c_.setZero(); // 初始化为零矩阵

        // 位置误差微分方程: dp_dot = dv
        F_c_(0, 2) = 1.0; // dp_N_dot = dv_N
        F_c_(1, 3) = 1.0; // dp_E_dot = dv_E

        // 速度误差微分方程: dv_dot = - (2*omega_ie + omega_en) x dv + f_n x dAtt - C_b_n * dBa + ...
        // dv_N_dot 对各状态的偏导
        // F_c_(2,0) = ... ; // dv_N_dot / dp_N (包含重力异常项，此处简化)
        if (R_N + h > 1e-6) {
            // 避免除零
            // dv_N_dot / dv_E (哥氏力 + 运载体速率引起的项)
            F_c_(2, 3) = -(2.0 * OMEGA_IE * sin(Lat) + current_vel_ned.y() * tan(Lat) / (R_N + h));
        } else {
            F_c_(2, 3) = -2.0 * OMEGA_IE * sin(Lat);
        }
        F_c_(2, 4) = f_n.y(); // dv_N_dot / dPsi_D ( = f_E * dPsi_D，其中 f_E 是东向比力)

        // dv_E_dot 对各状态的偏导
        F_c_(3, 2) = 2.0 * OMEGA_IE * sin(Lat); // dv_E_dot / dv_N
        if (R_N + h > 1e-6) {
            F_c_(3, 2) += current_vel_ned.y() * tan(Lat) / (R_N + h);
        }
        if (R_M + h > 1e-6) {
            F_c_(3, 2) += current_vel_ned.x() / (R_M + h); // dv_E_dot / dv_N (部分项来自运载体速率)
        }
        F_c_(3, 4) = -f_n.x(); // dv_E_dot / dPsi_D ( = -f_N * dPsi_D，其中 f_N 是北向比力)

        // 航向误差微分方程: dPsi_D_dot = ...
        // dPsi_D_dot 对各状态的偏导
        if (R_M + h > 1e-6) {
            // dPsi_D_dot / dp_N
            // F_c_(4,0) = (OMEGA_IE * sin(Lat) + current_vel_ned.y() * tan(Lat) / (R_N + h)) / (R_M + h); // 此项较复杂，通常较小
        }
        if (R_M + h > 1e-6) {
            F_c_(4, 2) = -1.0 / (R_M + h); // dPsi_D_dot / dv_N
        }
        if (R_N + h > 1e-6) {
            F_c_(4, 3) = tan(Lat) / (R_N + h); // dPsi_D_dot / dv_E
        }
        // dPsi_D_dot / dPsi_D (运载体速率引起的项, psi_dot 中的一部分)
        // F_c_(4,4) = -(OMEGA_IE * cos(Lat) + current_vel_ned.x() / (R_N + h)) * tan(Lat); // 此项复杂，依赖于 dPhi, dTheta，对于仅有dPsi的简化模型，此项可能为0或被忽略
        // 简化：dPsi_D_dot 由于 dPsi_D 本身的影响 (通常对航向误差自身较小)

        // 速度误差对加速度计零偏误差的偏导: dv_dot / dBa = -C_b_n
        F_c_.block<1, 3>(2, 5) = -C_b_n.row(0); // d(dv_N)/dBa
        F_c_.block<1, 3>(3, 5) = -C_b_n.row(1); // d(dv_E)/dBa

        // 航向误差对陀螺仪零偏误差的偏导: dPsi_D_dot / dBg_z = -C_b_n(2,z_axis_in_body)
        // 假设陀螺仪的Z轴与载体的Z轴（指向天）对齐，则 gyro_b.z() 主要影响航向。
        // d(dPsi_D)/dBg，这里dBg是b系下的误差，需要通过C_b_n转换到n系下的角速度误差，再影响dPsi_D
        // dPsi_D_dot = -omega_z_error_n = - (C_b_n * dBg_b)_z
        // 取C_b_n的第三行与dBg_b做点积
        F_c_.block<1, 3>(4, 8) = -C_b_n.row(2); // d(dPsi_D)/dBg (即 -[C_b_n(2,0), C_b_n(2,1), C_b_n(2,2)])

        // 零偏误差模型为随机游走, F_c 中对应项为0 (d(Bias)/dt = w_bias, w_bias是白噪声)

        // --- 离散时间状态转移矩阵 F_k = I + F_c * dt (一阶近似) ---
        Eigen::Matrix<double, 11, 11> F_k = Eigen::Matrix<double, 11, 11>::Identity() + F_c_ * dt;

        // --- 协方差矩阵传播 P_k|k-1 = F_k * P_k-1|k-1 * F_k^T + Q_k ---
        // Q_k 是离散时间过程噪声协方差矩阵
        Eigen::Matrix<double, 11, 11> Q_k = Q_cont_ * dt; // 简化形式: 假设 G (噪声输入矩阵) 为单位阵，且 F_c 变化缓慢
        // 更精确的形式: Q_k = Phi * G * Q_c_cont * G^T * Phi^T * dt (近似)
        // 或 Q_k = integral_0^dt (F(tau) * G * Q_c_cont * G^T * F(tau)^T) dtau
        // 对于零偏的简单随机游走模型, Q_k_bias = sigma_bias_psd * dt
        // 对于由白噪声驱动的速度/姿态, Q_k_vel = sigma_accel_psd * dt
        P_ = F_k * P_ * F_k.transpose() + Q_k;

        // 误差状态 dx_ 在校正后通常假设为零均值,
        // 所以 dx_ = F_k * dx_ (预测均值) 通常是 dx_ = 0。
        // 预测步骤主要是为了传播协方差矩阵 P_。
    }

    // EKF 更新步骤 (使用GNSS数据)
    bool updateGNSS(const Vector3d& ins_pos_llh,      /* INS 纬度 (弧度), 经度 (弧度), 高程 (米) */
                    const Vector3d& ins_vel_ned,      /* INS 北向速度, 东向速度, 地向速度 (米/秒) */
                    const double    ins_heading_rad,  /* INS 航向角/偏航角 (弧度, -PI 到 PI 或 0 到 2PI) */
                    const Vector3d& gnss_pos_llh,     /* GNSS 纬度 (弧度), 经度 (弧度), 高程 (米) */
                    const double    gnss_heading_rad, /* GNSS 航向角 (弧度) */
                    const double    gnss_ground_speed_ms) /* GNSS 地速 (米/秒) */ {
        // 1. 计算测量新息 y = z - h(x_hat_minus)
        // z = [pos_N_GNSS, pos_E_GNSS, Heading_GNSS, Speed_GNSS]^T
        Eigen::Matrix<double, 4, 1> y; // 新息向量

        const double Lat_ins   = ins_pos_llh.x(); // INS 纬度
        const double Lon_ins   = ins_pos_llh.y(); // INS 经度
        const double h_ins     = ins_pos_llh.z(); // INS 高程
        Vector3d     R_vec_ins = radii(Lat_ins);  // INS位置处的曲率半径
        const double R_M_ins   = R_vec_ins.x();   // 子午圈曲率半径
        const double R_N_ins   = R_vec_ins.y();   // 卯酉圈曲率半径

        // 位置新息 (将LLA差异转换为NED下的米制差异)
        y(0) = (gnss_pos_llh.x() - Lat_ins) * (R_M_ins + h_ins); // dN (米) = (Lat_GNSS - Lat_INS) * (R_M + h)
        y(1) = (gnss_pos_llh.y() - Lon_ins) * (R_N_ins + h_ins) * cos(Lat_ins);
        // dE (米) = (Lon_GNSS - Lon_INS) * (R_N + h) * cos(Lat)

        // 航向新息 (确保角度范围一致, 例如，归一化到 +/- PI)
        double heading_diff = gnss_heading_rad - ins_heading_rad;
        while (heading_diff > M_PI) heading_diff -= 2.0 * M_PI;
        while (heading_diff < -M_PI) heading_diff += 2.0 * M_PI;
        y(2) = heading_diff; // dPsi (弧度)

        // 地速新息
        const double ins_ground_speed = sqrt(pow(ins_vel_ned.x(), 2) + pow(ins_vel_ned.y(), 2)); // INS计算的地速
        y(3)                          = gnss_ground_speed_ms - ins_ground_speed;                 // dVg (米/秒)

        // 2. 构建观测矩阵 H_ (4x11)
        Eigen::Matrix<double, 4, 11> H_;
        H_.setZero(); // 初始化为零矩阵

        // 位置观测对应的H矩阵行 (dp_N, dp_E)
        H_(0, 0) = 1.0; // d(meas_N) / dp_N_error
        H_(1, 1) = 1.0; // d(meas_E) / dp_E_error

        // 航向观测对应的H矩阵行 (dPsi_D)
        H_(2, 4) = 1.0; // d(meas_Head) / dPsi_D_error

        // 地速观测对应的H矩阵行 (dv_N, dv_E)
        // d(Vg_meas - Vg_ins) / d(state_error)
        // Vg_ins = sqrt(v_N^2 + v_E^2)
        // d(Vg_ins) / dv_N_error = v_N / Vg_ins
        // d(Vg_ins) / dv_E_error = v_E / Vg_ins
        // 所以 H 元素是 -d(Vg_ins)/dv_error, 但由于新息是 GNSS - INS，所以 H 应该是 +d(Vg_ins)/dv_error
        // 或者，如果状态误差 dx 定义为 x_true - x_est，则 h(x_est + dx) approx h(x_est) + H dx
        // z - h(x_est) = H dx.  这里的 H 是 d(measurement_model)/d(true_state) @ x_est
        // measurement_model_speed = sqrt( (vN_true)^2 + (vE_true)^2 )
        // vN_true = vN_est + dvN_error, vE_true = vE_est + dvE_error
        // d(measurement_model_speed)/dvN_error = vN_est / Vg_est
        if (ins_ground_speed > 1e-3) {
            // 避免地速过低时除零
            H_(3, 2) = ins_vel_ned.x() / ins_ground_speed; // d(Vg_gnss - Vg_ins) / dv_N_error = - d(Vg_ins)/dv_N_error
            // 如果 y = z_gnss - (h_ins(x_nom) + H_ins * dx_nom_err), H_ins = d(h_ins)/dx_nom_err
            // 这里 y = z_gnss - h_ins(x_nom), H 对应于 z_gnss = h_true(x_nom + dx)
            // 所以 H(3,2) = d(Vg_true)/dvN_error = vN_nom / Vg_nom
            H_(3, 3) = ins_vel_ned.y() / ins_ground_speed; // d(Vg_true)/dvE_error = vE_nom / Vg_nom
        } else {
            H_(3, 2) = 0.0; // 速度过低时可不使用速度更新，或特殊处理
            H_(3, 3) = 0.0;
        }

        // 3. 计算卡尔曼增益 K_
        Eigen::Matrix<double, 11, 4> K_;
        Eigen::Matrix<double, 4, 4>  S = H_ * P_ * H_.transpose() + R_gnss_; // 新息协方差
        K_                             = P_ * H_.transpose() * S.inverse();  // 卡尔曼增益

        // 4. 更新误差状态估计 dx_
        dx_ = K_ * y;

        // 5. 更新协方差矩阵 P_ (使用Joseph形式，更稳定)
        Eigen::Matrix<double, 11, 11> I_KH = Eigen::Matrix<double, 11, 11>::Identity() - K_ * H_;
        P_                                 = I_KH * P_ * I_KH.transpose() + K_ * R_gnss_ * K_.transpose();
        // P_ = (Eigen::Matrix<double, 11, 11>::Identity() - K_ * H_) * P_; // 简单形式，数值稳定性较差

        return true;
    }

    // 将误差状态校正应用到主导航状态
    void applyCorrection(Vector3d& pos_LLH,     /* 纬度, 经度, 高程 (弧度, 弧度, 米) */
                         Vector3d& vel_NED,     /* 北向, 东向, 地向速度 (米/秒) */
                         Matrix3d& C_b_n,       /* 从b系到n系的方向余弦矩阵 */
                         Vector3d& accelBias_b, /* b系下的加速度计零偏 (米/秒^2) */
                         Vector3d& gyroBias_b) /* b系下的陀螺仪零偏 (弧度/秒) */ {
        /** --- 校正位置 (纬度, 经度，由 dp_N, dp_E 校正) ---
         * 注意: dx_[0] 是 dp_N (北向位置误差), dx_[1] 是 dp_E (东向位置误差)
         * 假设误差状态 dx_ 定义为 x_estimated - x_true
         * 则 x_true = x_estimated - dx_
         * 校正后的估计值 x_corrected (更接近x_true) = x_estimated - dx_
         */
        double   current_lat = pos_LLH.x();
        double   current_h   = pos_LLH.z();
        Vector3d R_vec       = radii(current_lat);
        double   R_M         = R_vec.x();
        double   R_N         = R_vec.y();

        pos_LLH.x() -= dx_(0) / (R_M + current_h); // Lat_corrected = Lat_old - dp_N / (R_M+h)
        pos_LLH.y() -= dx_(1) / ((R_N + current_h) * cos(current_lat));
        // Lon_corrected = Lon_old - dp_E / ((R_N+h)*cos(Lat))

        // --- 校正速度 (北向, 东向速度，由 dv_N, dv_E 校正) ---
        vel_NED.x() -= dx_(2); // vN_corrected = vN_old - dvN_error
        vel_NED.y() -= dx_(3); // vE_corrected = vE_old - dvE_error
        // vel_NED.z() (地向速度) 在此EKF状态中不进行校正

        // --- 校正姿态 (航向角，由 dPsi_D 校正) ---
        // Psi_corrected = Psi_old - dPsi_error
        // C_b_n_corrected = (I - [dAtt_err_n x]) * C_b_n_old
        // dAtt_err_n = [0, 0, dPsi_D]^T (NED坐标系下的姿态误差向量)
        const Vector3d att_err_ned(0, 0, dx_(4));                                        // 姿态误差向量，仅包含航向误差
        const Matrix3d C_err_matrix = Matrix3d::Identity() - Rotation::SkewSymmetric(att_err_ned); // (I - [phi_x]) 误差旋转矩阵
        // C_b_n_new = C_n_true_n_est * C_b_n_est
        // 如果 dx_[4] 是 Psi_est - Psi_true, 则 Psi_true = Psi_est - dx_[4]
        // 旋转应该是从估计姿态 n_est 到真实姿态 n_true
        // C_n_true_n_est 表示从 n_est 系到 n_true 系的旋转。
        // C_b_n_true = C_n_true_n_est * C_b_n_est
        C_b_n = C_err_matrix * C_b_n; // C_b_n_corrected = (I - skew(att_err_ned)) * C_b_n_estimated

        // 重新正交化 C_b_n 矩阵, 消除累积的计算误差
        const Eigen::JacobiSVD<Eigen::Matrix3d> svd(C_b_n, Eigen::ComputeFullU | Eigen::ComputeFullV);
        C_b_n = svd.matrixU() * svd.matrixV().transpose();

        // --- 校正零偏 ---
        // 假设 dx_ 中的零偏误差定义为: dBias_error = Bias_estimated - Bias_true
        // 则 Bias_true = Bias_estimated - dBias_error
        // 校正后的零偏 Bias_corrected = Bias_estimated - dBias_error
        accelBias_b -= dx_.segment<3>(5); // accelBias_corrected = accelBias_old - dBa_error
        gyroBias_b -= dx_.segment<3>(8);  // gyroBias_corrected = gyroBias_old - dBg_error

        // 校正后重置误差状态 dx_ 为零
        dx_.setZero();
    }

    // 获取当前误差状态向量
    VectorXd getStateError() const { return dx_; }
    // 获取当前协方差矩阵
    MatrixXd getCovariance() const { return P_; }
};

#endif //LOOSELY_COUPLED_EKF_H
