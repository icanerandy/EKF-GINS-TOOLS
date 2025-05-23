//
// Created by XTAD on 2025/5/23.
//

#ifndef ROTATION_H
#define ROTATION_H

#include "Eigen/Geometry"
#include "fmt/core.h"

using Eigen::Matrix3d;
using Eigen::Quaterniond;
using Eigen::Vector3d;

class Rotation {
public:
    static Quaterniond Matrix2Quaternion(const Matrix3d& matrix) {
        return Quaterniond(matrix);
    }

    static Matrix3d Quaternion2Matrix(const Quaterniond& quaternion) {
        return quaternion.toRotationMatrix();
    }

    // ZYX旋转顺序, 前右下的IMU, 输出RPY
    static Vector3d Matrix2Euler(const Matrix3d& dcm) {
        Vector3d euler;

        euler.y() = atan(-dcm(2, 0) / sqrt(dcm(2, 1) * dcm(2, 1) + dcm(2, 2) * dcm(2, 2)));

        if (dcm(2, 0) <= -0.999) {
            euler.x() = 0;
            euler.z() = atan2(dcm(1, 2) - dcm(0, 1), (dcm(0, 2) + dcm(1, 1)));
            fmt::println("[WARNING] Rotation::matrix2euler: Singular Euler Angle! Set the roll angle to 0!");
        } else if (dcm(2, 0) >= 0.999) {
            euler.x() = 0;
            euler.z() = M_PI + atan2((dcm(1, 2) + dcm(0, 1)), (dcm(0, 2) - dcm(1, 1)));
            fmt::println("[WARNING] Rotation::matrix2euler: Singular Euler Angle! Set the roll angle to 0!");
        } else {
            euler.x() = atan2(dcm(2, 1), dcm(2, 2));
            euler.z() = atan2(dcm(1, 0), dcm(0, 0));
        }

        // heading 0~2PI
        if (euler.z() < 0) {
            euler.z() = M_PI * 2 + euler.z();
        }

        return euler;
    }

    static Vector3d Quaternion2Euler(const Quaterniond& quaternion) {
        return Matrix2Euler(quaternion.toRotationMatrix());
    }

    static Quaterniond Rotvec2Quaternion(const Vector3d& rotvec) {
        const double   angle = rotvec.norm();
        const Vector3d vec   = rotvec.normalized();
        return Quaterniond(Eigen::AngleAxisd(angle, vec));
    }

    static Vector3d Quaternion2Vector(const Quaterniond& quaternion) {
        Eigen::AngleAxisd axisd(quaternion);
        return axisd.angle() * axisd.axis();
    }

    // RPY --> C_b^n, ZYX顺序
    static Matrix3d Euler2Matrix(const Vector3d& euler) {
        return Matrix3d(Eigen::AngleAxisd(euler.z(), Vector3d::UnitZ()) *
                        Eigen::AngleAxisd(euler.y(), Vector3d::UnitY()) *
                        Eigen::AngleAxisd(euler.x(), Vector3d::UnitX()));
    }

    static Quaterniond Euler2Quaternion(const Vector3d& euler) {
        return Quaterniond(Eigen::AngleAxisd(euler.z(), Vector3d::UnitZ()) *
                           Eigen::AngleAxisd(euler.y(), Vector3d::UnitY()) *
                           Eigen::AngleAxisd(euler.x(), Vector3d::UnitX()));
    }

    static Matrix3d SkewSymmetric(const Vector3d& vector) {
        Matrix3d mat;
        mat << 0, -vector(2), vector(1),
                vector(2), 0, -vector(0),
                -vector(1), vector(0), 0;
        return mat;
    }

    static Eigen::Matrix4d QuaternionLeft(const Quaterniond& q) {
        Eigen::Matrix4d ans;
        ans(0, 0)             = q.w();
        ans.block<1, 3>(0, 1) = -q.vec().transpose();
        ans.block<3, 1>(1, 0) = q.vec();
        ans.block<3, 3>(1, 1) = q.w() * Matrix3d::Identity() + SkewSymmetric(q.vec());
        return ans;
    }

    static Eigen::Matrix4d QuaternionRight(const Quaterniond& p) {
        Eigen::Matrix4d ans;
        ans(0, 0)             = p.w();
        ans.block<1, 3>(0, 1) = -p.vec().transpose();
        ans.block<3, 1>(1, 0) = p.vec();
        ans.block<3, 3>(1, 1) = p.w() * Matrix3d::Identity() - SkewSymmetric(p.vec());
        return ans;
    }
};

#endif //ROTATION_H
