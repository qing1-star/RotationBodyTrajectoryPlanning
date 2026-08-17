#pragma once

#include <RotationBodyTrajectoryPlanning/Core/PlanningTypes.h>

#include <Eigen/Core>
#include <Eigen/Geometry>

namespace smrobot::spray::rotationbody
{
    constexpr double metersPerMillimeter = 0.001;
    constexpr double millimetersPerMeter = 1000.0;
    constexpr double pi = 3.141592653589793238462643383279502884;

    double millimetersToMeters(double millimeters) noexcept;
    double metersToMillimeters(double meters) noexcept;
    double degreesToRadians(double degrees) noexcept;
    double radiansToDegrees(double radians) noexcept;

    Eigen::Vector3d signedAxisVector(SignedAxis axis) noexcept;
    bool arePerpendicular(SignedAxis first, SignedAxis second, double tolerance = 1.0e-12) noexcept;
    bool isFiniteTransform(const Eigen::Isometry3d& transform) noexcept;

    Eigen::Matrix3d rotationFromRollPitchYaw(const Eigen::Vector3d& rollPitchYawRadians) noexcept;
    Eigen::Isometry3d makeTransform(const TransformComponents& components) noexcept;
    PlanningResult<TransformComponents> transformComponents(const Eigen::Isometry3d& transform);
    Eigen::Matrix3d rotationAligningVectorToVector(
        const Eigen::Vector3d& from,
        const Eigen::Vector3d& to) noexcept;
}
