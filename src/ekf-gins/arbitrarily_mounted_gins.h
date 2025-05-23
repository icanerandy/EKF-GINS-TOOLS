//
// Created by XTAD on 2025/5/23.
//

#ifndef ARBITRARILY_MOUNTED_GINS_H
#define ARBITRARILY_MOUNTED_GINS_H

#include "common/types.h"
#include "motion_recognizer.h"
#include "loosely_coupled_ekf.h"

// --- 核心GNSS/INS集成与失准角估计类 ---
class ArbitrarilyMountedGINS {
private:
    double imu_dt_; /* IMU采样间隔 */
    bool   is_initialized_{false};
    bool   gnss_ekf_initialized_{false};

    IMU  last_imu_data_;  /* 保存最新的IMU数据 */
    GNSS last_gnss_data_; /* 保存最新的GNSS数据 */

    NavigationState current_nav_state_;
    Vector3d        current_accel_bias_;
    Vector3d        current_gyro_bias_;

    Misalignment misalignment_;
    Matrix3d     cmm1; /* 传感器(m)到水平传感器(m1) */

    std::unique_ptr<LooselyCoupledEKF> ekf_;
    std::unique_ptr<MotionRecognizer>  motion_recognizer_;

    // 用于粗略航向失准估计的状态
    double       psi_m1_; /* 水平对准的传感器m1系的估计航向 */
    int          valid_psi_m1_estimation_count_;
    const int    kPsiM1CountThreshHold_{ 15 };          /* cite: 136 */
    const double kPsiM1PosErrorThreshHoldAlpha_{ 0.5 }; /* Eq 13的alpha [cite: 132] (需调整) */
    const double kPsiM1IntegrationIntervalT{ 1.0 };     /* Eq 13的T [cite: 132] */

public:
    explicit ArbitrarilyMountedGINS(double imu_sampling_interval);
};


#endif //ARBITRARILY_MOUNTED_GINS_H
