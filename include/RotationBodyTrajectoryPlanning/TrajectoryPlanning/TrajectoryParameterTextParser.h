#pragma once

#include <RotationBodyTrajectoryPlanning/TrajectoryPlanning/TrajectoryTypes.h>

#include <string_view>
#include <vector>

namespace smrobot::spray::rotationbody
{
    class TrajectoryParameterTextParser
    {
    public:
        static PlanningResult<std::vector<TrajectoryGenerationParameters>> parse(
            std::string_view text);
    };
}
