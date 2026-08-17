#pragma once

#include <RotationBodyTrajectoryPlanning/TrajectoryPlanning/TrajectoryTypes.h>

#include <vector>

namespace smrobot::spray::rotationbody
{
    class TrajectoryEditor
    {
    public:
        static PlanningResult<void> interpolateRange(
            PlannedTrajectory& trajectory,
            const std::vector<std::size_t>& selectedIndices,
            double intervalSeconds);
        static PlanningResult<void> transformSelected(
            PlannedTrajectory& trajectory,
            const std::vector<std::size_t>& selectedIndices,
            const TransformComponents& delta);
        static PlanningResult<void> swapDirection(PlannedTrajectory& trajectory);
        static PlanningResult<void> rebuildDerived(
            PlannedTrajectory& trajectory,
            double startOffsetSeconds = 0.0);
        static TrajectoryMetrics calculateMetrics(
            const std::vector<TrajectoryPosePoint>& points);
    };
}
