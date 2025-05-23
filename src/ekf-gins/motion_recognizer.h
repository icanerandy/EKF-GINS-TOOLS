//
// Created by XTAD on 2025/5/23.
//

#ifndef MOTION_RECOGNIZER_H
#define MOTION_RECOGNIZER_H

#include "Eigen/Geometry"

using Eigen::Vector3d;
using Eigen::VectorXd;

// --- 运动模式识别类 ---
class MotionRecognizer {
public:
    MotionRecognizer(double N_still, double                     N_turn,
        double epsilon_acc_still_increment, double epsilon_gyro_still_increment,
        double mu_turn_gyro_increment_sq, double   sampling_interval);

    // 输入原始加速度(m/s^2)和角速度(rad/s)
    void Update(const Vector3d& accel_raw, const Vector3d& gyro_raw);

    bool IsStill() const;

    bool IsTurning() const;

    bool IsStraight() const;

private:
    double N_still_, N_turn_;
    double alpha_still_, beta_still_;
    double alpha_turn_, beta_turn_;
    double epsilon_acc_still_sq_, epsilon_gyro_still_sq_; /* 标准差阈值 */
    double mu_turn_gyro_increment_sq_;                    /* 转弯的陀螺仪平方和阈值 */
    double sampling_interval_;

    VectorXd means_still_;     /* 6x1 vector for E[acc_inc], E[gyro_inc] */
    VectorXd std_deviation_still_; /* 6x1 vector for std_dev(acc_inc), std_dev(gyro_inc) */
    bool     still_initialized_{ false };

    double    current_TD_{ 0.0 }; /* EMA of sum_j(G_ij^2) */
    bool      turn_initialized_{ false };
    int       successive_turn_count_{ 0 };
    const int TURN_COUNT_INDICATOR_{ 5 }; /* [cite: 86] */
};

#endif //MOTION_RECOGNIZER_H
