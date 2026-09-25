#pragma once

#include <RotationBodyTrajectoryPlanning/TrajectoryPlanning/TrajectoryTypes.h>

namespace smrobot::spray::rotationbody
{
    class TrajectoryPlanner
    {
    public:
        static constexpr double automaticSpacingMeters = 0.001;
        static constexpr std::size_t minimumPointCount = 2;
        static constexpr double minimumSprayDistanceMeters = -0.1;
        static constexpr double maximumSprayDistanceMeters = 9.999;
        static constexpr double maximumAbsoluteTiltRadians =
            999.0 * 3.14159265358979323846 / 180.0;

        static Eigen::Matrix3d baseFromToolAtZeroTilt();
        static Eigen::Matrix3d levelSprayAxisAroundLocalY(
            const Eigen::Matrix3d& baseFromTool);

        static PlanningResult<std::size_t> suggestedPointCount(
            double pathLengthMeters);
        static PlanningResult<PlannedTrajectory> generate(
            const SprayBoundary& boundary,
            const TrajectoryGenerationParameters& parameters);
        static PlanningResult<PlannedTrajectory> generate(
            const SprayBoundary& boundary,
            const TrajectoryGenerationParameters& parameters,
            const Eigen::Isometry3d& baseFromPlanning);
    };
}
