#pragma once

#include <RotationBodyTrajectoryPlanning/TrajectoryPlanning/TrajectoryTypes.h>

namespace smrobot::spray::rotationbody
{
    class RelativeHelicalTrajectoryBuilder
    {
    public:
        static PlanningResult<std::vector<TrajectoryPosePoint>> build(
            const std::vector<TrajectoryPosePoint>& linearPoints,
            double positionerRpm,
            double startOffsetSeconds = 0.0);
    };
}
