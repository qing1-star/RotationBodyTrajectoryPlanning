#pragma once

#include <RotationBodyTrajectoryPlanning/TrajectoryPlanning/TrajectoryTypes.h>

#include <optional>
#include <string>

namespace smrobot::spray::rotationbody
{
    class PublishedTrajectoryPlanStore
    {
    public:
        static std::optional<PublishedTrajectoryPlan> read(
            const std::string& serializedPayload,
            std::string* errorMessage = nullptr);
    };
}
