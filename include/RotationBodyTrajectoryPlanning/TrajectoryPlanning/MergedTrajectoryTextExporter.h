#pragma once

#include <RotationBodyTrajectoryPlanning/TrajectoryPlanning/TrajectoryTypes.h>

#include <string>

namespace smrobot::spray::rotationbody
{
    class MergedTrajectoryTextExporter
    {
    public:
        // Formats all saved passes as base-frame 4x4 poses in millimeters,
        // followed by each pass's group-relative timestamp in seconds.
        static PlanningResult<std::string> format(
            const PublishedTrajectoryPlan& plan);
    };
}
