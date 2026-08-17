#pragma once

#include <RotationBodyTrajectoryPlanning/TrajectoryPlanning/TrajectoryTypes.h>

namespace smrobot::spray::rotationbody
{
    class TrajectoryPlanner
    {
    public:
        static constexpr double automaticSpacingMeters = 0.001;
        static constexpr std::size_t minimumPointCount = 2;
        static constexpr std::size_t maximumPointCount = 10001;

        static PlanningResult<std::size_t> suggestedPointCount(
            double pathLengthMeters);
        static PlanningResult<PlannedTrajectory> generate(
            const SprayBoundary& boundary,
            const TrajectoryGenerationParameters& parameters);
    };
}
