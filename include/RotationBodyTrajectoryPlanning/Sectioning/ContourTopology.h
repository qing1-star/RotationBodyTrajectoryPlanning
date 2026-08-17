#pragma once

#include <RotationBodyTrajectoryPlanning/Core/PlanningTypes.h>

#include <Eigen/Core>

#include <cstddef>
#include <vector>

namespace smrobot::spray::rotationbody
{
    struct YzSegment
    {
        Eigen::Vector2d first = Eigen::Vector2d::Zero();
        Eigen::Vector2d second = Eigen::Vector2d::Zero();
    };

    struct ContourTopologyResult
    {
        std::vector<SectionContour> contours;
        std::size_t inputSegmentCount{ 0 };
        std::size_t uniqueSegmentCount{ 0 };
        std::size_t discardedDegenerateSegmentCount{ 0 };
        std::size_t discardedDuplicateSegmentCount{ 0 };
        std::size_t branchNodeCount{ 0 };
    };

    class ContourTopology
    {
    public:
        static PlanningResult<ContourTopologyResult> buildContours(
            const std::vector<YzSegment>& segments,
            double weldToleranceMeters,
            double planeX = 0.0);

        static PlanningResult<std::size_t> selectTargetOuterContour(
            const std::vector<SectionContour>& contours,
            double toleranceMeters);
    };
}
