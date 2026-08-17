#pragma once

#include <RotationBodyTrajectoryPlanning/Core/PlanningTypes.h>
#include <RotationBodyTrajectoryPlanning/Core/TriangleMesh.h>

namespace smrobot::spray::rotationbody
{
    class SimulationBlockPlacementSolver
    {
    public:
        // rotationAxisInMesh maps to planning +Z and toothOutwardAxisInMesh
        // maps to planning +Y. The partial block is then placed at Xmid=0,
        // Zmin=0 and Ymax=motherMaximumDiameterMeters/2.
        static PlanningResult<AlignmentResult> solve(
            const TriangleMesh& mesh,
            SignedAxis rotationAxisInMesh,
            SignedAxis toothOutwardAxisInMesh,
            double motherMaximumDiameterMeters);
    };
}
