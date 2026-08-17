#include <RotationBodyTrajectoryPlanning/Core/TransformUtils.h>

#include <algorithm>
#include <cmath>

namespace smrobot::spray::rotationbody
{
    double millimetersToMeters(double millimeters) noexcept
    {
        return millimeters * metersPerMillimeter;
    }

    double metersToMillimeters(double meters) noexcept
    {
        return meters * millimetersPerMeter;
    }

    double degreesToRadians(double degrees) noexcept
    {
        return degrees * pi / 180.0;
    }

    double radiansToDegrees(double radians) noexcept
    {
        return radians * 180.0 / pi;
    }

    Eigen::Vector3d signedAxisVector(SignedAxis axis) noexcept
    {
        switch (axis)
        {
        case SignedAxis::PositiveX:
            return Eigen::Vector3d::UnitX();
        case SignedAxis::NegativeX:
            return -Eigen::Vector3d::UnitX();
        case SignedAxis::PositiveY:
            return Eigen::Vector3d::UnitY();
        case SignedAxis::NegativeY:
            return -Eigen::Vector3d::UnitY();
        case SignedAxis::PositiveZ:
            return Eigen::Vector3d::UnitZ();
        case SignedAxis::NegativeZ:
            return -Eigen::Vector3d::UnitZ();
        }
        return Eigen::Vector3d::UnitZ();
    }

    bool arePerpendicular(SignedAxis first, SignedAxis second, double tolerance) noexcept
    {
        return std::abs(signedAxisVector(first).dot(signedAxisVector(second))) <= tolerance;
    }

    bool isFiniteTransform(const Eigen::Isometry3d& transform) noexcept
    {
        if (!transform.matrix().allFinite())
        {
            return false;
        }
        const Eigen::Matrix3d rotation = transform.linear();
        return std::abs(rotation.determinant() - 1.0) <= 1.0e-8 &&
            (rotation.transpose() * rotation - Eigen::Matrix3d::Identity()).norm() <= 1.0e-8;
    }

    Eigen::Matrix3d rotationFromRollPitchYaw(const Eigen::Vector3d& rollPitchYawRadians) noexcept
    {
        const Eigen::AngleAxisd roll(rollPitchYawRadians.x(), Eigen::Vector3d::UnitX());
        const Eigen::AngleAxisd pitch(rollPitchYawRadians.y(), Eigen::Vector3d::UnitY());
        const Eigen::AngleAxisd yaw(rollPitchYawRadians.z(), Eigen::Vector3d::UnitZ());
        return (yaw * pitch * roll).toRotationMatrix();
    }

    Eigen::Isometry3d makeTransform(const TransformComponents& components) noexcept
    {
        Eigen::Isometry3d result = Eigen::Isometry3d::Identity();
        result.linear() = rotationFromRollPitchYaw(components.rollPitchYawRadians);
        result.translation() = components.translationMeters;
        return result;
    }

    PlanningResult<TransformComponents> transformComponents(const Eigen::Isometry3d& transform)
    {
        if (!isFiniteTransform(transform))
        {
            return PlanningResult<TransformComponents>::failure(
                PlanningErrorCode::NonFiniteGeometry,
                "Transform must contain a finite rigid rotation and translation.");
        }

        TransformComponents result;
        result.translationMeters = transform.translation();
        const Eigen::Matrix3d& rotation = transform.linear();
        const double sinPitch = std::clamp(-rotation(2, 0), -1.0, 1.0);
        const double pitch = std::asin(sinPitch);
        const double cosPitch = std::cos(pitch);
        double roll = 0.0;
        double yaw = 0.0;
        if (std::abs(cosPitch) > 1.0e-10)
        {
            roll = std::atan2(rotation(2, 1), rotation(2, 2));
            yaw = std::atan2(rotation(1, 0), rotation(0, 0));
        }
        else
        {
            roll = 0.0;
            yaw = std::atan2(-rotation(0, 1), rotation(1, 1));
        }
        result.rollPitchYawRadians = { roll, pitch, yaw };
        return PlanningResult<TransformComponents>::success(result);
    }

    Eigen::Matrix3d rotationAligningVectorToVector(
        const Eigen::Vector3d& from,
        const Eigen::Vector3d& to) noexcept
    {
        if (!from.allFinite() || !to.allFinite() || from.norm() <= 1.0e-15 || to.norm() <= 1.0e-15)
        {
            return Eigen::Matrix3d::Identity();
        }
        const Eigen::Vector3d source = from.normalized();
        const Eigen::Vector3d target = to.normalized();
        const double dot = std::clamp(source.dot(target), -1.0, 1.0);
        if (dot >= 1.0 - 1.0e-12)
        {
            return Eigen::Matrix3d::Identity();
        }
        if (dot <= -1.0 + 1.0e-12)
        {
            const Eigen::Vector3d helper =
                std::abs(source.x()) < 0.8 ? Eigen::Vector3d::UnitX() : Eigen::Vector3d::UnitY();
            return Eigen::AngleAxisd(pi, source.cross(helper).normalized()).toRotationMatrix();
        }
        const Eigen::Vector3d axis = source.cross(target).normalized();
        return Eigen::AngleAxisd(std::acos(dot), axis).toRotationMatrix();
    }
}
