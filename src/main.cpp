#include "fmt/core.h"

#include "ekf-gins/arbitrarily_mounted_gins.h"
#include "ekf-gins/eskf-gins.h"
#include "common/earth.h"
#include <iostream>


// --- 主函数示例 (仅为演示如何使用此类) ---
int main(int argc, char* argv[]) {
    fmt::println("任意安装车载GNSS/INS失准角估计算法:");
    GnssImsEskf ekf; // 创建EKF对象

    // 使用一些值初始化协方差矩阵 P
    MatrixXd initial_P = MatrixXd::Identity(11,11) * 0.1;
    // 如果已知，可以设置特定的初始不确定性，例如：
    initial_P(0,0) = pow(10.0,2); initial_P(1,1) = pow(10.0,2); // 位置误差 (米^2)
    initial_P(2,2) = pow(1.0,2);  initial_P(3,3) = pow(1.0,2);  // 速度误差 (米/秒)^2
    initial_P(4,4) = pow(5.0*DEG_TO_RAD,2);                     // 航向误差 (弧度^2)
    for(int i=5; i<8; ++i) initial_P(i,i) = pow(1e-2,2);        // 加速度计零偏误差 (米/秒^2)^2
    for(int i=8; i<11; ++i) initial_P(i,i) = pow(1e-3*DEG_TO_RAD,2); // 陀螺仪零偏误差 (弧度/秒)^2
    ekf.initialize(initial_P); // 初始化EKF

    // --- 用于预测的伪数据 ---
    Vector3d accel_b(0.1, 0.05, -9.8); // b系下的比力 (米/秒^2)
    Vector3d gyro_b(0.001, 0.002, 0.0005); // b系下的角速度 (弧度/秒)
    Matrix3d C_b_n = Matrix3d::Identity(); // 方向余弦矩阵示例
    Vector3d current_vel_ned(10.0, 1.0, 0.0); // n系下的速度 (米/秒)
    Vector3d current_pos_llh(34.0*DEG_TO_RAD, -118.0*DEG_TO_RAD, 100.0); // 位置 (纬度,经度,高程) (弧度,弧度,米)
    Vector3d current_accel_bias = Vector3d::Zero(); // 当前加速度计零偏
    Vector3d current_gyro_bias = Vector3d::Zero();  // 当前陀螺仪零偏
    double dt = 0.01; // 时间间隔 (秒)

    // --- 预测循环 (示例) ---
    for (int i=0; i<100; ++i) { // 模拟1秒
        // 从您的INS解算中更新 current_C_b_n, vel_ned, pos_llh
        ekf.predict(accel_b, gyro_b, C_b_n, current_vel_ned, current_pos_llh, current_accel_bias, current_gyro_bias, dt);
    }

    // --- 用于GNSS更新的伪数据 ---
    Vector3d gnss_pos_llh = current_pos_llh; // 假设GNSS与INS位置暂时完美匹配
    gnss_pos_llh.x() += 0.5 / Earth::MeridianPrimeVerticalRadius(current_pos_llh.x()).x(); // 人为增加0.5米的北向误差
    double gnss_heading_rad = 0.0 * DEG_TO_RAD; // GNSS航向角示例
    double gnss_ground_speed_ms = 10.05; // GNSS地速示例

    // ins_heading_rad 应来自您的INS当前的航向估计值
    // 例如，从 C_b_n 矩阵中提取航向角 (Psi)
    double ins_yaw_estimate = atan2(C_b_n(1,0), C_b_n(0,0)); // 假设是标准的 C_b_n 定义

    ekf.updateGNSS(current_pos_llh, current_vel_ned, ins_yaw_estimate,
    gnss_pos_llh, gnss_heading_rad, gnss_ground_speed_ms); // 执行GNSS更新

    VectorXd final_dx = ekf.getStateError(); // 获取最终的误差状态
    MatrixXd final_P = ekf.getCovariance(); // 获取最终的协方差矩阵
    std::cout << "最终误差状态 dx:\n" << final_dx << std::endl;
    std::cout << "最终协方差 P (对角线):\n" << final_P.diagonal().transpose() << std::endl;

    ekf.applyCorrection(current_pos_llh, current_vel_ned, C_b_n, current_accel_bias, current_gyro_bias); // 应用校正
    std::cout << "校正后的加速度计零偏:\n" << current_accel_bias << std::endl;

    return 0;
}