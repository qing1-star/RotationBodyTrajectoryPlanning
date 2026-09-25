#include <RotationBodyTrajectoryPlanning/TrajectoryPlanning/TrajectoryPlanner.h>

#include <RotationBodyTrajectoryPlanning/TrajectoryPlanning/RelativeHelicalTrajectoryBuilder.h>

#include <algorithm>
#include <cmath>
#include <limits>

namespace smrobot::spray::rotationbody
{
    namespace
    {
        constexpr double epsilon = 1.0e-12;

        bool finiteParameters(const TrajectoryGenerationParameters& parameters) noexcept
        {
            return std::isfinite(parameters.sprayDistanceMeters) &&
                std::isfinite(parameters.tiltRadians) &&
                std::isfinite(parameters.speedMetersPerSecond) &&
                std::isfinite(parameters.startExtensionMeters) &&
                std::isfinite(parameters.endExtensionMeters) &&
                std::isfinite(parameters.positionerRpm);
        }
    }

    Eigen::Matrix3d TrajectoryPlanner::baseFromToolAtZeroTilt()
    {
        Eigen::Matrix3d orientation;
        orientation.col(0) = -Eigen::Vector3d::UnitZ();
        orientation.col(1) = Eigen::Vector3d::UnitX();
        orientation.col(2) = -Eigen::Vector3d::UnitY();
        return orientation;
    }

    Eigen::Matrix3d TrajectoryPlanner::levelSprayAxisAroundLocalY(
        const Eigen::Matrix3d& baseFromTool)
    {
        const double toolXVertical = baseFromTool.col(0).z();
        const double sprayAxisVertical = baseFromTool.col(2).z();
        double correctionRadians = std::atan2(
            -sprayAxisVertical,
            toolXVertical);
        constexpr double halfPi = 0.5 * 3.14159265358979323846;
        constexpr double pi = 3.14159265358979323846;
        if(correctionRadians > halfPi) {
            correctionRadians -= pi;
        } else if(correctionRadians < -halfPi) {
            correctionRadians += pi;
        }
        return baseFromTool *
            Eigen::AngleAxisd(
                correctionRadians,
                Eigen::Vector3d::UnitY()).toRotationMatrix();
    }

    PlanningResult<std::size_t> TrajectoryPlanner::suggestedPointCount(
        double pathLengthMeters)
    {
        if(!std::isfinite(pathLengthMeters) || pathLengthMeters <= 0.0) {
            return PlanningResult<std::size_t>::failure(
                PlanningErrorCode::InvalidArgument,
                "Trajectory length must be finite and positive.");
        }
        const double segmentCount = std::ceil(pathLengthMeters / automaticSpacingMeters);
        if(!std::isfinite(segmentCount) ||
            segmentCount >= static_cast<double>(std::numeric_limits<std::size_t>::max())) {
            return PlanningResult<std::size_t>::failure(
                PlanningErrorCode::InvalidArgument,
                "The automatic 1 mm sampling cannot be represented by the point-count type.");
        }
        return PlanningResult<std::size_t>::success(
            std::max(minimumPointCount, static_cast<std::size_t>(segmentCount) + 1));
    }

    PlanningResult<PlannedTrajectory> TrajectoryPlanner::generate(
        const SprayBoundary& boundary,
        const TrajectoryGenerationParameters& parameters)
    {
        return generate(boundary, parameters, Eigen::Isometry3d::Identity());
    }

    PlanningResult<PlannedTrajectory> TrajectoryPlanner::generate(
        const SprayBoundary& boundary,
        const TrajectoryGenerationParameters& parameters,
        const Eigen::Isometry3d& baseFromPlanning)
    {
        if(!finiteParameters(parameters) ||
            parameters.sprayDistanceMeters < minimumSprayDistanceMeters ||
            parameters.sprayDistanceMeters > maximumSprayDistanceMeters ||
            parameters.speedMetersPerSecond <= 0.0 ||
            parameters.startExtensionMeters < 0.0 ||
            parameters.endExtensionMeters < 0.0 ||
            std::abs(parameters.tiltRadians) > maximumAbsoluteTiltRadians ||
            !baseFromPlanning.matrix().allFinite()) {
            return PlanningResult<PlannedTrajectory>::failure(
                PlanningErrorCode::InvalidArgument,
                "Trajectory parameters are outside their supported finite ranges.");
        }
        if(!std::isfinite(boundary.minimumZ) || !std::isfinite(boundary.maximumZ) ||
            boundary.maximumZ <= boundary.minimumZ ||
            !std::isfinite(boundary.outerY(boundary.minimumZ)) ||
            !std::isfinite(boundary.outerY(boundary.maximumZ))) {
            return PlanningResult<PlannedTrajectory>::failure(
                PlanningErrorCode::InsufficientRegionData,
                "A valid confirmed spray boundary is required to generate a trajectory.");
        }

        const Eigen::Vector3d canonicalEdgeA(
            0.0,
            boundary.outerY(boundary.maximumZ),
            boundary.maximumZ);
        const Eigen::Vector3d canonicalEdgeB(
            0.0,
            boundary.outerY(boundary.minimumZ),
            boundary.minimumZ);

        const Eigen::Vector3d rawTangent = canonicalEdgeB - canonicalEdgeA;
        const double rawLength = rawTangent.norm();
        if(!std::isfinite(rawLength) || rawLength <= epsilon) {
            return PlanningResult<PlannedTrajectory>::failure(
                PlanningErrorCode::DegenerateGeometry,
                "The confirmed boundary outer edge is degenerate.");
        }
        Eigen::Vector3d edgeA = canonicalEdgeA;
        Eigen::Vector3d edgeB = canonicalEdgeB;
        if(parameters.reversed) {
            std::swap(edgeA, edgeB);
        }
        const Eigen::Vector3d tangent = (edgeB - edgeA) / rawLength;

        PlannedTrajectory trajectory;
        trajectory.parameters = parameters;
        trajectory.sourceBoundary = boundary;
        trajectory.targetSurfaceStart = edgeA - parameters.startExtensionMeters * tangent;
        trajectory.targetSurfaceEnd = edgeB + parameters.endExtensionMeters * tangent;

        // At zero tilt, tool X/Y/Z point along base -Z/+X/-Y. Apply the
        // requested tilt around local tool Y, then express that fixed base
        // orientation in planning coordinates for every trajectory point.
        const Eigen::Matrix3d baseFromTiltedTool =
            baseFromToolAtZeroTilt() *
            Eigen::AngleAxisd(
                parameters.tiltRadians,
                Eigen::Vector3d::UnitY()).toRotationMatrix();
        const Eigen::Matrix3d planningFromTool =
            baseFromPlanning.linear().transpose() * baseFromTiltedTool;
        const Eigen::Vector3d sprayDirection = planningFromTool.col(2);
        const Eigen::Vector3d startPosition = trajectory.targetSurfaceStart -
            parameters.sprayDistanceMeters * sprayDirection;
        const Eigen::Vector3d endPosition = trajectory.targetSurfaceEnd -
            parameters.sprayDistanceMeters * sprayDirection;
        const Eigen::Vector3d movement = endPosition - startPosition;
        const double pathLength = movement.norm();
        if(!std::isfinite(pathLength) || pathLength <= epsilon) {
            return PlanningResult<PlannedTrajectory>::failure(
                PlanningErrorCode::DegenerateGeometry,
                "The generated trajectory start and end positions coincide.");
        }

        std::size_t pointCount = parameters.pointCount;
        if(pointCount == 0) {
            PlanningResult<std::size_t> count = suggestedPointCount(pathLength);
            if(!count) {
                return PlanningResult<PlannedTrajectory>::failure(
                    count.error.code,
                    count.error.message);
            }
            pointCount = count.value;
            trajectory.parameters.pointCount = pointCount;
        }
        if(pointCount < minimumPointCount) {
            return PlanningResult<PlannedTrajectory>::failure(
                PlanningErrorCode::InvalidArgument,
                "Trajectory point count must be at least 2.");
        }

        const double duration = pathLength / parameters.speedMetersPerSecond;
        trajectory.linearPoints.reserve(pointCount);
        for(std::size_t index = 0; index < pointCount; ++index) {
            const double ratio = static_cast<double>(index) /
                static_cast<double>(pointCount - 1);
            TrajectoryPosePoint point;
            point.timeSeconds = ratio * duration;
            point.planningFromTool.linear() = planningFromTool;
            point.planningFromTool.translation() =
                (1.0 - ratio) * startPosition + ratio * endPosition;
            trajectory.linearPoints.push_back(std::move(point));
        }
        trajectory.metrics.pathLengthMeters = pathLength;
        trajectory.metrics.durationSeconds = duration;
        trajectory.metrics.minimumPointIntervalSeconds =
            duration / static_cast<double>(pointCount - 1);
        trajectory.metrics.maximumPointIntervalSeconds =
            trajectory.metrics.minimumPointIntervalSeconds;

        PlanningResult<std::vector<TrajectoryPosePoint>> helical =
            RelativeHelicalTrajectoryBuilder::build(
                trajectory.linearPoints,
                parameters.positionerRpm);
        if(!helical) {
            return PlanningResult<PlannedTrajectory>::failure(
                helical.error.code,
                helical.error.message);
        }
        trajectory.relativeHelicalPoints = std::move(helical.value);
        return PlanningResult<PlannedTrajectory>::success(std::move(trajectory));
    }
}
