//
// Created by XTAD on 2025/5/23.
//

#include "arbitrarily_mounted_gins.h"
#include "arbitrarily_mounted_gins.h"

ArbitrarilyMountedGINS::ArbitrarilyMountedGINS(double imu_sampling_interval)
    : imu_dt_(imu_sampling_interval) {
    /**
     * 初始化运动识别器参数(基于论文, 但阈值可能需针对20ms调整)
     * 论文静止参数: N=8, window_time=1.6s (for 0.2s interval)
     * 论文转弯参数: N=5, window_time=1.0s (for 0.2s interval)
     * 假设保持1.6s窗口时间用于静止检测, 对于20ms采样率 n_still = 1.6/0.020 = 80
     * 假设保持1.0s窗口时间用于转弯检测, 对于20ms采样率 n_turn = 1.0/0.020 = 50
     */
    double N_still = 80.0; /* 基于1.6s窗口 @ 20ms采样 */
    double N_turn  = 50.0; /* 基于1.0s窗口 @ 20ms采样 */

    /**
     * 阈值需要根据20ms采样间隔下的增量大小重新标定!
     * 以下为示例, 实际需要实验确定
     */
    double epsilon_acc_still_increment  = 0.01 / 10.0; /* 论文值 / (0.2s/0.020s)*/
    double epsilon_gyro_still_increment = 0.0004 / 10.0;
    double mu_turn_gyro_increment_sq    = 0.0005 / (10.0 * 10.0); /* G_ij^2, 所以分母平方 */

    motion_recognizer_ = std::make_unique<MotionRecognizer>{
        N_still, N_turn,
        epsilon_acc_still_increment, epsilon_gyro_still_increment,
        mu_turn_gyro_increment_sq, // 注意这是 G_ij^2 的阈值
        imu_dt_
    };

    // --- 初始化EKF ---
    ekf_ = std::make_unique<LooselyCoupledEKF>();
}
