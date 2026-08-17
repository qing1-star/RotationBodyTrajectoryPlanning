#pragma once

#include <RotationBodyTrajectoryPlanning/TrajectoryPlanning/TrajectoryTypes.h>

#include <string>

namespace smrobot::spray::rotationbody
{
    class TrajectoryGroupEditor
    {
    public:
        static PlanningResult<std::string> addOrUpdate(
            TrajectoryGroup& group,
            const PlannedTrajectory& trajectory,
            const std::string& editingPassId = {});
        static PlanningResult<void> remove(
            TrajectoryGroup& group,
            const std::string& passId);
        static PlanningResult<void> setVisible(
            TrajectoryGroup& group,
            const std::string& passId,
            bool visible);
        static PlanningResult<void> setTransitionAfter(
            TrajectoryGroup& group,
            const std::string& passId,
            double seconds,
            double zeroIntervalPositionToleranceMeters = 1.0e-4);
        static const TrajectoryPass* find(
            const TrajectoryGroup& group,
            const std::string& passId) noexcept;
        static PlanningResult<void> refreshSchedule(TrajectoryGroup& group);
        static PlanningResult<void> validate(
            const TrajectoryGroup& group,
            double zeroIntervalPositionToleranceMeters = 1.0e-4);
    };
}
