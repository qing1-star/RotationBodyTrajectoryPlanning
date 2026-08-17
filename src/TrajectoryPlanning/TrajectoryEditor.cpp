#include <RotationBodyTrajectoryPlanning/TrajectoryPlanning/TrajectoryEditor.h>

#include <RotationBodyTrajectoryPlanning/TrajectoryPlanning/RelativeHelicalTrajectoryBuilder.h>

#include <Eigen/Geometry>

#include <algorithm>
#include <cmath>
#include <limits>
#include <set>

namespace smrobot::spray::rotationbody
{
    namespace
    {
        Eigen::Matrix3d deltaRotation(const Eigen::Vector3d& rollPitchYaw)
        {
            return (
                Eigen::AngleAxisd(rollPitchYaw.z(), Eigen::Vector3d::UnitZ()) *
                Eigen::AngleAxisd(rollPitchYaw.y(), Eigen::Vector3d::UnitY()) *
                Eigen::AngleAxisd(rollPitchYaw.x(), Eigen::Vector3d::UnitX()))
                .toRotationMatrix();
        }

        TrajectoryPosePoint interpolatePoint(
            const TrajectoryPosePoint& first,
            const TrajectoryPosePoint& second,
            double timeSeconds)
        {
            const double duration = second.timeSeconds - first.timeSeconds;
            const double ratio = (timeSeconds - first.timeSeconds) / duration;
            Eigen::Quaterniond firstRotation(first.planningFromTool.linear());
            Eigen::Quaterniond secondRotation(second.planningFromTool.linear());
            TrajectoryPosePoint point;
            point.timeSeconds = timeSeconds;
            point.interpolated = true;
            point.planningFromTool.translation() =
                (1.0 - ratio) * first.planningFromTool.translation() +
                ratio * second.planningFromTool.translation();
            point.planningFromTool.linear() =
                firstRotation.slerp(ratio, secondRotation).normalized().toRotationMatrix();
            return point;
        }
    }

    PlanningResult<void> TrajectoryEditor::interpolateRange(
        PlannedTrajectory& trajectory,
        const std::vector<std::size_t>& selectedIndices,
        double intervalSeconds)
    {
        if(trajectory.linearPoints.size() < 2 || selectedIndices.size() < 2 ||
            !std::isfinite(intervalSeconds) || intervalSeconds <= 0.0) {
            return PlanningResult<void>::failure(
                PlanningErrorCode::InvalidArgument,
                "Select at least two trajectory points and enter a positive interpolation interval.");
        }
        const auto bounds = std::minmax_element(selectedIndices.begin(), selectedIndices.end());
        const std::size_t firstIndex = *bounds.first;
        const std::size_t lastIndex = *bounds.second;
        if(firstIndex >= trajectory.linearPoints.size() ||
            lastIndex >= trajectory.linearPoints.size() || firstIndex >= lastIndex) {
            return PlanningResult<void>::failure(
                PlanningErrorCode::InvalidArgument,
                "The selected interpolation range is invalid.");
        }

        std::vector<TrajectoryPosePoint> expanded;
        expanded.reserve(trajectory.linearPoints.size() * 2);
        expanded.insert(
            expanded.end(),
            trajectory.linearPoints.begin(),
            trajectory.linearPoints.begin() + static_cast<std::ptrdiff_t>(firstIndex));
        for(std::size_t index = firstIndex; index < lastIndex; ++index) {
            const TrajectoryPosePoint& first = trajectory.linearPoints[index];
            const TrajectoryPosePoint& second = trajectory.linearPoints[index + 1];
            expanded.push_back(first);
            const double duration = second.timeSeconds - first.timeSeconds;
            if(duration <= 0.0) {
                return PlanningResult<void>::failure(
                    PlanningErrorCode::InvalidArgument,
                    "Trajectory timestamps must increase inside the interpolation range.");
            }
            for(double time = first.timeSeconds + intervalSeconds;
                time < second.timeSeconds - 1.0e-12;
                time += intervalSeconds) {
                expanded.push_back(interpolatePoint(first, second, time));
            }
        }
        expanded.push_back(trajectory.linearPoints[lastIndex]);
        expanded.insert(
            expanded.end(),
            trajectory.linearPoints.begin() + static_cast<std::ptrdiff_t>(lastIndex + 1),
            trajectory.linearPoints.end());
        trajectory.linearPoints = std::move(expanded);
        trajectory.parameters.pointCount = trajectory.linearPoints.size();
        return rebuildDerived(trajectory);
    }

    PlanningResult<void> TrajectoryEditor::transformSelected(
        PlannedTrajectory& trajectory,
        const std::vector<std::size_t>& selectedIndices,
        const TransformComponents& delta)
    {
        if(selectedIndices.empty() || !delta.translationMeters.allFinite() ||
            !delta.rollPitchYawRadians.allFinite()) {
            return PlanningResult<void>::failure(
                PlanningErrorCode::InvalidArgument,
                "Select trajectory points and enter a finite pose adjustment.");
        }
        const Eigen::Matrix3d rotation = deltaRotation(delta.rollPitchYawRadians);
        std::set<std::size_t> uniqueIndices(selectedIndices.begin(), selectedIndices.end());
        for(std::size_t index : uniqueIndices) {
            if(index >= trajectory.linearPoints.size()) {
                return PlanningResult<void>::failure(
                    PlanningErrorCode::InvalidArgument,
                    "A selected trajectory point index is out of range.");
            }
            TrajectoryPosePoint& point = trajectory.linearPoints[index];
            point.planningFromTool.translation() += delta.translationMeters;
            point.planningFromTool.linear() = rotation * point.planningFromTool.linear();
        }
        return rebuildDerived(trajectory);
    }

    PlanningResult<void> TrajectoryEditor::swapDirection(PlannedTrajectory& trajectory)
    {
        if(trajectory.linearPoints.size() < 2) {
            return PlanningResult<void>::failure(
                PlanningErrorCode::InvalidArgument,
                "Generate a trajectory before swapping its direction.");
        }
        const double duration = trajectory.linearPoints.back().timeSeconds;
        std::reverse(trajectory.linearPoints.begin(), trajectory.linearPoints.end());
        for(TrajectoryPosePoint& point : trajectory.linearPoints) {
            point.timeSeconds = duration - point.timeSeconds;
            point.planningFromTool.linear().col(0) *= -1.0;
            point.planningFromTool.linear().col(1) *= -1.0;
        }
        std::swap(trajectory.targetSurfaceStart, trajectory.targetSurfaceEnd);
        trajectory.parameters.reversed = !trajectory.parameters.reversed;
        return rebuildDerived(trajectory);
    }

    PlanningResult<void> TrajectoryEditor::rebuildDerived(
        PlannedTrajectory& trajectory,
        double startOffsetSeconds)
    {
        if(trajectory.linearPoints.size() < 2) {
            return PlanningResult<void>::failure(
                PlanningErrorCode::InvalidArgument,
                "A trajectory requires at least two points.");
        }
        PlanningResult<std::vector<TrajectoryPosePoint>> derived =
            RelativeHelicalTrajectoryBuilder::build(
                trajectory.linearPoints,
                trajectory.parameters.positionerRpm,
                startOffsetSeconds);
        if(!derived) {
            return PlanningResult<void>::failure(derived.error.code, derived.error.message);
        }
        trajectory.relativeHelicalPoints = std::move(derived.value);
        trajectory.parameters.pointCount = trajectory.linearPoints.size();
        trajectory.metrics = calculateMetrics(trajectory.linearPoints);
        return PlanningResult<void>::success();
    }

    TrajectoryMetrics TrajectoryEditor::calculateMetrics(
        const std::vector<TrajectoryPosePoint>& points)
    {
        TrajectoryMetrics metrics;
        if(points.size() < 2) {
            return metrics;
        }
        metrics.minimumPointIntervalSeconds = std::numeric_limits<double>::infinity();
        for(std::size_t index = 1; index < points.size(); ++index) {
            metrics.pathLengthMeters += (
                points[index].planningFromTool.translation() -
                points[index - 1].planningFromTool.translation()).norm();
            const double interval = points[index].timeSeconds - points[index - 1].timeSeconds;
            metrics.minimumPointIntervalSeconds =
                std::min(metrics.minimumPointIntervalSeconds, interval);
            metrics.maximumPointIntervalSeconds =
                std::max(metrics.maximumPointIntervalSeconds, interval);
        }
        metrics.durationSeconds = points.back().timeSeconds - points.front().timeSeconds;
        if(!std::isfinite(metrics.minimumPointIntervalSeconds)) {
            metrics.minimumPointIntervalSeconds = 0.0;
        }
        return metrics;
    }
}
