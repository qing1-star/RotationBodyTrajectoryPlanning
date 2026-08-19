#include <RotationBodyTrajectoryPlanning/TrajectoryPlanning/TrajectoryPlanner.h>

#include <RotationBodyTrajectoryPlanning/TrajectoryPlanning/RelativeHelicalTrajectoryBuilder.h>

#include <algorithm>
#include <cmath>
#include <limits>

namespace smrobot::spray::rotationbody
{
    namespace
    {
        constexpr double pi = 3.14159265358979323846;
        constexpr double epsilon = 1.0e-12;

        Eigen::Vector3d inwardOuterNormal(
            const Eigen::Vector3d& outerEdgeTangent)
        {
            // The confirmed outer edge lies in the planning YZ plane.  Its
            // in-plane normal is chosen toward the inner (minimum-Y) side of
            // the spray boundary, so zero tilt is always normal to the edge.
            Eigen::Vector3d normal(
                0.0,
                -outerEdgeTangent.z(),
                outerEdgeTangent.y());
            if(normal.norm() <= epsilon) {
                return Eigen::Vector3d::Zero();
            }
            normal.normalize();
            if(normal.y() > 0.0) {
                normal = -normal;
            }
            return normal;
        }

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
            segmentCount > static_cast<double>(maximumPointCount - 1)) {
            return PlanningResult<std::size_t>::failure(
                PlanningErrorCode::InvalidArgument,
                "The automatic 1 mm sampling would exceed the supported point count.");
        }
        return PlanningResult<std::size_t>::success(
            std::max(minimumPointCount, static_cast<std::size_t>(segmentCount) + 1));
    }

    PlanningResult<PlannedTrajectory> TrajectoryPlanner::generate(
        const SprayBoundary& boundary,
        const TrajectoryGenerationParameters& parameters)
    {
        if(!finiteParameters(parameters) || parameters.sprayDistanceMeters <= 0.0 ||
            parameters.speedMetersPerSecond <= 0.0 ||
            parameters.startExtensionMeters < 0.0 ||
            parameters.endExtensionMeters < 0.0 ||
            std::abs(parameters.tiltRadians) >= 80.0 * pi / 180.0) {
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
        const Eigen::Vector3d canonicalTangent = rawTangent / rawLength;
        const Eigen::Vector3d zeroTiltSprayDirection =
            inwardOuterNormal(canonicalTangent);
        if(!zeroTiltSprayDirection.allFinite() ||
            zeroTiltSprayDirection.norm() <= epsilon) {
            return PlanningResult<PlannedTrajectory>::failure(
                PlanningErrorCode::DegenerateGeometry,
                "The confirmed boundary outer edge has no valid in-plane normal.");
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

        // Tilt is measured from the outer-edge normal and rotates around the
        // trajectory frame's local Y axis.  The canonical edge orientation is
        // used here so swapping A/B changes travel direction only, not the
        // physical spray direction or the meaning of a positive tilt angle.
        const Eigen::Vector3d canonicalToolY =
            zeroTiltSprayDirection.cross(canonicalTangent).normalized();
        const Eigen::Vector3d sprayDirection =
            Eigen::AngleAxisd(parameters.tiltRadians, canonicalToolY) *
            zeroTiltSprayDirection;
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

        // Reversing a pass changes only its travel direction.  Keep the tool
        // frame tied to the canonical boundary direction so a return pass
        // does not introduce a 180-degree roll around the spray axis.
        Eigen::Vector3d toolX = canonicalTangent -
            canonicalTangent.dot(sprayDirection) * sprayDirection;
        if(toolX.norm() <= epsilon) {
            return PlanningResult<PlannedTrajectory>::failure(
                PlanningErrorCode::DegenerateGeometry,
                "The travel direction is parallel to the spray direction.");
        }
        toolX.normalize();
        Eigen::Vector3d toolY = sprayDirection.cross(toolX);
        if(toolY.norm() <= epsilon) {
            return PlanningResult<PlannedTrajectory>::failure(
                PlanningErrorCode::DegenerateGeometry,
                "A right-handed trajectory frame cannot be constructed.");
        }
        toolY.normalize();
        toolX = toolY.cross(sprayDirection).normalized();

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
        if(pointCount < minimumPointCount || pointCount > maximumPointCount) {
            return PlanningResult<PlannedTrajectory>::failure(
                PlanningErrorCode::InvalidArgument,
                "Trajectory point count must be between 2 and 10001.");
        }

        const double duration = pathLength / parameters.speedMetersPerSecond;
        trajectory.linearPoints.reserve(pointCount);
        for(std::size_t index = 0; index < pointCount; ++index) {
            const double ratio = static_cast<double>(index) /
                static_cast<double>(pointCount - 1);
            TrajectoryPosePoint point;
            point.timeSeconds = ratio * duration;
            point.planningFromTool.linear().col(0) = toolX;
            point.planningFromTool.linear().col(1) = toolY;
            point.planningFromTool.linear().col(2) = sprayDirection;
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
