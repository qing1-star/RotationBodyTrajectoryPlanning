#include <RotationBodyTrajectoryPlanning/Alignment/ModelTransformOperations.h>

#include <RotationBodyTrajectoryPlanning/Core/TransformUtils.h>

#include <Eigen/Geometry>

#include <algorithm>
#include <cmath>
#include <limits>

namespace smrobot::spray::rotationbody
{
    namespace
    {
        PlanningResult<Eigen::Isometry3d> composePlanningDelta(
            const Eigen::Isometry3d& reference,
            const TransformComponents& delta)
        {
            if (!isFiniteTransform(reference))
            {
                return PlanningResult<Eigen::Isometry3d>::failure(
                    PlanningErrorCode::InvalidArgument,
                    "The reference model transform must be a finite rigid transform.");
            }
            if (!delta.translationMeters.allFinite() || !delta.rollPitchYawRadians.allFinite())
            {
                return PlanningResult<Eigen::Isometry3d>::failure(
                    PlanningErrorCode::InvalidArgument,
                    "Model translation and rotation increments must be finite.");
            }

            const Eigen::Isometry3d result = makeTransform(delta) * reference;
            if (!isFiniteTransform(result))
            {
                return PlanningResult<Eigen::Isometry3d>::failure(
                    PlanningErrorCode::AlignmentFailed,
                    "The model transform increment did not produce a finite rigid transform.");
            }
            return PlanningResult<Eigen::Isometry3d>::success(result);
        }
    }

    PlanningResult<Eigen::Isometry3d> ModelTransformOperations::fromBaseline(
        const Eigen::Isometry3d& automaticBaseline,
        const TransformComponents& delta)
    {
        return composePlanningDelta(automaticBaseline, delta);
    }

    PlanningResult<Eigen::Isometry3d> ModelTransformOperations::applyIncrement(
        const Eigen::Isometry3d& current,
        const TransformComponents& delta)
    {
        return composePlanningDelta(current, delta);
    }

    PlanningResult<Eigen::Isometry3d> ModelTransformOperations::flip(
        const TriangleMesh& mesh,
        const Eigen::Isometry3d& current,
        const Eigen::Vector3d& directedRotaryAxisInMesh,
        const Eigen::Vector3d& pointOnRotaryAxisInMesh)
    {
        const PlanningResult<void> validation = mesh.validate();
        if (!validation)
        {
            return PlanningResult<Eigen::Isometry3d>::failure(
                validation.error.code, validation.error.message);
        }
        if (!isFiniteTransform(current))
        {
            return PlanningResult<Eigen::Isometry3d>::failure(
                PlanningErrorCode::InvalidArgument,
                "The current model transform must be a finite rigid transform.");
        }
        if (!directedRotaryAxisInMesh.allFinite() ||
            directedRotaryAxisInMesh.norm() <= std::numeric_limits<double>::epsilon() ||
            !pointOnRotaryAxisInMesh.allFinite())
        {
            return PlanningResult<Eigen::Isometry3d>::failure(
                PlanningErrorCode::InvalidArgument,
                "The rotary-axis direction and point must be finite, and the direction must be non-zero.");
        }

        Eigen::Vector3d currentAxisInPlanning =
            current.linear() * directedRotaryAxisInMesh.normalized();
        const Eigen::Vector3d currentAxisPointInPlanning =
            current * pointOnRotaryAxisInMesh;
        currentAxisInPlanning.normalize();

        double maximumAxial = -std::numeric_limits<double>::infinity();
        for (const Eigen::Vector3d& position : mesh.positions)
        {
            const Eigen::Vector3d currentPosition = current * position;
            maximumAxial = std::max(
                maximumAxial,
                (currentPosition - currentAxisPointInPlanning).dot(currentAxisInPlanning));
        }
        if (!std::isfinite(maximumAxial))
        {
            return PlanningResult<Eigen::Isometry3d>::failure(
                PlanningErrorCode::NonFiniteGeometry,
                "The model has no finite opposite axial end for flipping.");
        }

        const Eigen::Vector3d newBottomInCurrentPlanning =
            currentAxisPointInPlanning + currentAxisInPlanning * maximumAxial;
        const Eigen::Vector3d newBottomInMesh = current.inverse() * newBottomInCurrentPlanning;

        // Map the old bottom-to-top direction to -Z. Consequently the reversed
        // direction from the new bottom to the old bottom is exactly +Z.
        const Eigen::Matrix3d planningNormalization = rotationAligningVectorToVector(
            currentAxisInPlanning, -Eigen::Vector3d::UnitZ());
        Eigen::Isometry3d result = Eigen::Isometry3d::Identity();
        result.linear() = planningNormalization * current.linear();
        result.translation() = -result.linear() * newBottomInMesh;
        if (!isFiniteTransform(result))
        {
            return PlanningResult<Eigen::Isometry3d>::failure(
                PlanningErrorCode::AlignmentFailed,
                "Flipping the model did not produce a finite rigid transform.");
        }

        const PlanningResult<Eigen::AlignedBox3d> flippedBounds = mesh.bounds(result);
        if (!flippedBounds)
        {
            return PlanningResult<Eigen::Isometry3d>::failure(
                flippedBounds.error.code, flippedBounds.error.message);
        }
        const double zTolerance = std::max(mesh.scale() * 1.0e-10, 1.0e-12);
        if (std::abs(flippedBounds.value.min().z()) > zTolerance ||
            (result * newBottomInMesh).norm() > zTolerance)
        {
            return PlanningResult<Eigen::Isometry3d>::failure(
                PlanningErrorCode::AlignmentFailed,
                "The flipped model could not satisfy the bottom-axis origin and Z-minimum constraints.");
        }
        return PlanningResult<Eigen::Isometry3d>::success(result);
    }

    PlanningResult<Eigen::Isometry3d> ModelTransformOperations::flip(
        const TriangleMesh& mesh,
        const Eigen::Isometry3d& current)
    {
        const Eigen::Vector3d assumedAxisInMesh =
            current.linear().transpose() * Eigen::Vector3d::UnitZ();
        const Eigen::Vector3d assumedBottomInMesh = current.inverse() * Eigen::Vector3d::Zero();
        return flip(mesh, current, assumedAxisInMesh, assumedBottomInMesh);
    }

    PlanningResult<Eigen::Isometry3d> ModelTransformOperations::reset(
        const Eigen::Isometry3d& automaticBaseline)
    {
        if (!isFiniteTransform(automaticBaseline))
        {
            return PlanningResult<Eigen::Isometry3d>::failure(
                PlanningErrorCode::InvalidArgument,
                "The automatic alignment baseline must be a finite rigid transform.");
        }
        return PlanningResult<Eigen::Isometry3d>::success(automaticBaseline);
    }
}
