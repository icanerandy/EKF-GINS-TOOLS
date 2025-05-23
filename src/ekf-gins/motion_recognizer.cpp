//
// Created by XTAD on 2025/5/23.
//

#include "motion_recognizer.h"

MotionRecognizer::MotionRecognizer(const double N_still, const double                     N_turn,
                                   const double epsilon_acc_still_increment, const double epsilon_gyro_still_increment,
                                   const double mu_turn_gyro_increment_sq, const double   sampling_interval)
    : N_still_(N_still),
      N_turn_(N_turn),
      epsilon_acc_still_sq_(epsilon_acc_still_increment),
      epsilon_gyro_still_sq_(epsilon_gyro_still_increment),
      mu_turn_gyro_increment_sq_(mu_turn_gyro_increment_sq), // 这是 G_ij_inc^2 的阈值
      sampling_interval_(sampling_interval),
      means_still_(6),
      std_deviation_still_(6) {
    if (N_still_ > 1) {
        alpha_still_ = 1.0 / N_still_;
        beta_still_  = (N_still_ - 1.0) / N_still_;
    } else {
        alpha_still_ = 1.0;
        beta_still_  = 0.0;
    }

    if (N_turn_ > 1) {
        alpha_turn_ = 1.0 / N_turn_;
        beta_turn_  = (N_turn_ - 1.0) / N_turn_;
    } else {
        alpha_turn_ = 1.0;
        beta_turn_  = 0.0;
    }

    means_still_.setZero();
    std_deviation_still_.setZero();
}

void MotionRecognizer::Update(const Vector3d& accel_raw, const Vector3d& gyro_raw) {
    // 1. 计算增量
    const Vector3d accel_inc = accel_raw * sampling_interval_;
    const Vector3d gyro_inc  = gyro_raw * sampling_interval_;

    VectorXd current_increments(6);
    current_increments << accel_inc, gyro_inc;

    // 2. 更新静止检测的EMA统计量 [cite: 60, 61] (标准差)
    if (!still_initialized_) {
        means_still_ = current_increments;
        std_deviation_still_.setZero(); // 初始标准差为0
        still_initialized_ = true;
    } else {
        for (int i = 0; i < 6; ++i) {
            const double mean_prev = means_still_[i];
            const double dev_prev  = std_deviation_still_[i];
            const double val       = current_increments[i];

            means_still_[i]         = beta_still_ * mean_prev + alpha_still_ * val;
            std_deviation_still_[i] = beta_still_ * dev_prev + alpha_still_ * (val - means_still_[i]);
            if (std_deviation_still_[i] < 0) std_deviation_still_[i] = 0; // 防负
        }
    }

    // 3. 更新转弯检测的EMA统计量 (TD_i) [cite: 80]
    // TD_i = EMA of sum(G_gyro_inc_j^2)
    const double sum_sq_gyro_inc = gyro_inc.squaredNorm();
    if (!turn_initialized_) {
        current_TD_       = sum_sq_gyro_inc; // 初始值
        turn_initialized_ = true;
    } else {
        current_TD_ = beta_turn_ * current_TD_ + alpha_turn_ * sum_sq_gyro_inc;
    }

    // 更新连续转弯计数 [cite: 86]
    if (current_TD_ > mu_turn_gyro_increment_sq_) {
        successive_turn_count_++;
    } else {
        successive_turn_count_ = 0; // 不满足则清零
    }

    printf("accel: %f, %f, %f\n", std_deviation_still_[0], std_deviation_still_[1], std_deviation_still_[2]);
    printf("gyro: %f, %f, %f\n", std_deviation_still_[3], std_deviation_still_[4], std_deviation_still_[5]);
    printf("current_TD_: %f\n", current_TD_);
}

bool MotionRecognizer::IsStill() const {
    if (!still_initialized_) return false;
    bool accel_still = true;
    for (int i = 0; i < 3; ++i) {
        // 加速度计增量方差
        if (std_deviation_still_[i] >= epsilon_acc_still_sq_) {
            // 比较方差与阈值的平方
            accel_still = false;
            break;
        }
    }
    if (!accel_still) return false;

    for (int i = 3; i < 6; ++i) {
        // 陀螺仪增量方差
        if (std_deviation_still_[i] >= epsilon_gyro_still_sq_) {
            return false;
        }
    }
    return true;
}

bool MotionRecognizer::IsTurning() const {
    if (!turn_initialized_) return false;
    // 论文中提到用连续转弯计数来消除错误检测 [cite: 86]
    return successive_turn_count_ >= TURN_COUNT_INDICATOR_ && (current_TD_ > mu_turn_gyro_increment_sq_);
}

bool MotionRecognizer::IsStraight() const {
    // 直线运动: 不是静止, 也不是(已确认的)转弯
    return !IsStill() && !IsTurning();
}
