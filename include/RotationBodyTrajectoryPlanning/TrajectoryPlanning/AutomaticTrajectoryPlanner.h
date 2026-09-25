#pragma once

#include <RotationBodyTrajectoryPlanning/TrajectoryPlanning/TrajectoryTypes.h>

#include <vector>

namespace smrobot::spray::rotationbody
{
    enum class AutomaticTrajectoryMode
    {
        Dual,
        Triple
    };

    struct AutomaticTrajectoryPlan
    {
        std::vector<TrajectoryGenerationParameters> trajectories;
        double upperWallNormalAngleRadians{ 0.0 };
        double lowerWallNormalAngleRadians{ 0.0 };
        double toothTopNormalAngleRadians{ 0.0 };
        double toothBottomNormalAngleRadians{ 0.0 };
        double dualTiltLimitRadians{ 0.0 };
    };

    class AutomaticTrajectoryPlanner
    {
    public:
        static PlanningResult<AutomaticTrajectoryPlan> plan(
            const SectionContour& contour,
            const RegionAssignment& regions,
            AutomaticTrajectoryMode mode);
    };
}
