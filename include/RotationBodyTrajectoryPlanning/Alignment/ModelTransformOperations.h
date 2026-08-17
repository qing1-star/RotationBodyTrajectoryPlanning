#pragma once

#include <RotationBodyTrajectoryPlanning/Core/PlanningTypes.h>
#include <RotationBodyTrajectoryPlanning/Core/TriangleMesh.h>

namespace smrobot::spray::rotationbody
{
    class ModelTransformOperations
    {
    public:
        // Rebuilds an editable transform from the automatic alignment baseline.
        // The delta is expressed in the planning frame and uses Rz * Ry * Rx.
        static PlanningResult<Eigen::Isometry3d> fromBaseline(
            const Eigen::Isometry3d& automaticBaseline,
            const TransformComponents& delta);

        // Applies one planning-frame increment to an existing editable transform.
        static PlanningResult<Eigen::Isometry3d> applyIncrement(
            const Eigen::Isometry3d& current,
            const TransformComponents& delta);

        // Exchanges the two axial ends around the supplied mesh-space baseline
        // axis. The direction must point from the current bottom toward the
        // opposite end, and pointOnRotaryAxisInMesh may be the old bottom center
        // or any other point on the same axis line. The current planning axis is
        // derived through current.linear(). The new bottom-to-top direction is
        // normalized to planning +Z, its axis center is placed at the planning
        // origin, and the transformed mesh therefore has Zmin=0. This is a
        // re-normalization command; an arbitrary manual tilt is intentionally
        // removed while the current around-axis phase is retained where possible.
        static PlanningResult<Eigen::Isometry3d> flip(
            const TriangleMesh& mesh,
            const Eigen::Isometry3d& current,
            const Eigen::Vector3d& directedRotaryAxisInMesh,
            const Eigen::Vector3d& pointOnRotaryAxisInMesh);

        // Compatibility overload. It is valid only when the current rotary axis
        // is planning +Z and the current bottom axis center is the origin.
        static PlanningResult<Eigen::Isometry3d> flip(
            const TriangleMesh& mesh,
            const Eigen::Isometry3d& current);

        // Reset means returning to the automatic result, never to Identity.
        static PlanningResult<Eigen::Isometry3d> reset(
            const Eigen::Isometry3d& automaticBaseline);
    };
}
