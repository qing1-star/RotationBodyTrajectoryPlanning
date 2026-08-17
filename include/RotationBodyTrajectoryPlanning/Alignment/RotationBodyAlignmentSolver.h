#pragma once

#include <RotationBodyTrajectoryPlanning/Core/PlanningTypes.h>
#include <RotationBodyTrajectoryPlanning/Core/TriangleMesh.h>

namespace smrobot::spray::rotationbody
{
    class RotationBodyAlignmentSolver
    {
    public:
        // Estimates the rotary axis from all three PCA directions, selects the
        // end farther from the minimum-diameter station as the planning bottom,
        // and maps that directed axis to planning +Z.
        static PlanningResult<AlignmentResult> solve(const TriangleMesh& mesh);
    };
}
