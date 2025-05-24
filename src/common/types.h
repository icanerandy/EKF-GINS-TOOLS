//
// Created by XTAD on 2025/5/23.
//

#ifndef TYPES_H
#define TYPES_H

#include "Eigen/Geometry"

using Eigen::Vector2d;
using Eigen::Vector3d;
using Eigen::VectorXd;
using Eigen::Matrix3d;
using Eigen::MatrixXd;

/**
 * IMU原始数据(加速度 m/s^2, 角速度 rad/s)
 */
struct IMU {
    double   time;
    Vector3d accel_raw; /* x, y, z 轴加速度(m/s^2) - 传感器(m)坐标系下 */
    Vector3d gyro_raw;  /* x, y, z 轴角速度(rad/s) - 传感器(m)坐标系下 */
};

/**
 * GNSS数据
 */
struct GNSS {
    double   time;
    Vector3d llh;      /* 维度(rad), 经度(rad), 高程(m) */
    double   vel;     /* 对地速度(m/s) */
    double   heading; /* GNSS航向角(rad), 通常来自速度矢量, 运动时才可靠 */
    bool     valid;   /* GNSS是否有效 */
};

/**
 * 导航状态
 */
struct NavigationState {
    double   time;
    Vector2d ll;       /* 维度(rad), 经度(rad) */
    double   vel;      /* 对地速度(m/s) */
    Vector3d attitude; /* 滚转(rad), 俯仰(rad), 航向(rad) - 欧拉角表示车体(b)姿态 */
    Matrix3d cbn;      /* 从车体(b)坐标系到导航(n)坐标系的转换矩阵 */
};

#endif //TYPES_H
