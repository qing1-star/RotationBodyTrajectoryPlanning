#pragma once

#include <RotationBodyTrajectoryPlanning/Core/PlanningTypes.h>
#include <RotationBodyTrajectoryPlanning/Core/TriangleMesh.h>

#include <Eigen/Geometry>

#include <cstddef>
#include <vector>

namespace smrobot::spray::rotationbody
{
    struct YzSectionOptions
    {
        double relativePlaneTolerance{ 1.0e-9 };
        double relativeWeldTolerance{ 5.0e-9 };
        double minimumToleranceMeters{ 1.0e-12 };
    };

    struct YzSectionResult
    {
        SectionContour targetContour;
        std::vector<SectionContour> candidateContours;
        double meshScaleMeters{ 0.0 };
        double planeToleranceMeters{ 0.0 };
        double weldToleranceMeters{ 0.0 };
        std::size_t processedTriangleCount{ 0 };
        std::size_t coplanarTriangleCount{ 0 };
        std::size_t rawIntersectionSegmentCount{ 0 };
        std::size_t positiveYSegmentCount{ 0 };
        std::size_t discardedSegmentCount{ 0 };
        std::size_t branchNodeCount{ 0 };
    };

    class YzSectionExtractor
    {
    public:
        static PlanningResult<YzSectionResult> extract(
            const TriangleMesh& mesh,
            const Eigen::Isometry3d& planningFromMesh,
            const YzSectionOptions& options = {});

        static PlanningResult<SectionContour> extractTargetContour(
            const TriangleMesh& mesh,
            const Eigen::Isometry3d& planningFromMesh,
            const YzSectionOptions& options = {});
    };
}
