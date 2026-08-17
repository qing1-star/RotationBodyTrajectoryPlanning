#pragma once

#include <RotationBodyTrajectoryPlanning/Core/PlanningTypes.h>

namespace smrobot::spray::rotationbody
{
    class SprayBoundaryBuilder
    {
    public:
        static PlanningResult<SprayBoundary> build(
            const SectionContour& contour,
            const RegionAssignment& assignment,
            BoundaryMode mode);
    };
}
