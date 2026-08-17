#pragma once

#include <RotationBodyTrajectoryPlanning/Core/PlanningTypes.h>

#include <cstddef>

namespace smrobot::spray::rotationbody
{
    struct ToothRecognitionOptions
    {
        std::size_t minimumResampleCount{ 96 };
        std::size_t maximumResampleCount{ 2048 };
        double smoothingWindowFraction{ 0.015 };
        double minimumRadialProminenceFraction{ 0.025 };
        double axialSurfaceThreshold{ 0.62 };
        double minimumConfidence{ 0.28 };
        double transitionExtentFraction{ 0.035 };
    };

    class ToothRegionRecognizer
    {
    public:
        static PlanningResult<RegionAssignment> recognize(
            const SectionContour& contour,
            const ToothRecognitionOptions& options = {});
    };
}
