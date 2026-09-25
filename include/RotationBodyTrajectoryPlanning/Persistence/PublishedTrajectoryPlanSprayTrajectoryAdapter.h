#pragma once

#include <RotationBodyTrajectoryPlanning/TrajectoryPlanning/TrajectoryTypes.h>
#include <SprayTrajectoryCore/SprayTrajectory.h>

#include <optional>
#include <string>

namespace smrobot::spray::rotationbody
{
    class PublishedTrajectoryPlanSprayTrajectoryAdapter
    {
    public:
        // Converts the saved, positioner-compensated helical passes into the
        // time-stamped base-frame trajectory consumed by coating prediction.
        static std::optional<spraytrajectory::SprayTrajectory> convert(
            const PublishedTrajectoryPlan& plan,
            std::string* errorMessage = nullptr);
    };
}
