#include <RotationBodyTrajectoryPlanning/Alignment/SimulationBlockPlacementSolver.h>

#include <RotationBodyTrajectoryPlanning/Core/TransformUtils.h>

#include <algorithm>
#include <cmath>
#include <sstream>

namespace smrobot::spray::rotationbody
{
    namespace
    {
        bool isKnownSignedAxis(SignedAxis axis) noexcept
        {
            switch (axis)
            {
            case SignedAxis::PositiveX:
            case SignedAxis::NegativeX:
            case SignedAxis::PositiveY:
            case SignedAxis::NegativeY:
            case SignedAxis::PositiveZ:
            case SignedAxis::NegativeZ:
                return true;
            }
            return false;
        }
    }

    PlanningResult<AlignmentResult> SimulationBlockPlacementSolver::solve(
        const TriangleMesh& mesh,
        SignedAxis rotationAxisInMesh,
        SignedAxis toothOutwardAxisInMesh,
        double motherMaximumDiameterMeters)
    {
        const PlanningResult<void> validation = mesh.validate();
        if (!validation)
        {
            return PlanningResult<AlignmentResult>::failure(
                validation.error.code, validation.error.message);
        }
        if (!std::isfinite(motherMaximumDiameterMeters) || motherMaximumDiameterMeters <= 0.0)
        {
            return PlanningResult<AlignmentResult>::failure(
                PlanningErrorCode::InvalidArgument,
                "The mother-part maximum diameter must be finite and greater than zero.");
        }
        if (!isKnownSignedAxis(rotationAxisInMesh) || !isKnownSignedAxis(toothOutwardAxisInMesh))
        {
            return PlanningResult<AlignmentResult>::failure(
                PlanningErrorCode::InvalidAxisSelection,
                "The simulation-block axis selection is not a known signed Cartesian axis.");
        }
        if (!arePerpendicular(rotationAxisInMesh, toothOutwardAxisInMesh))
        {
            return PlanningResult<AlignmentResult>::failure(
                PlanningErrorCode::InvalidAxisSelection,
                "The original rotary axis and tooth-outward axis must be perpendicular.");
        }

        const Eigen::Vector3d planningZInMesh = signedAxisVector(rotationAxisInMesh);
        const Eigen::Vector3d planningYInMesh = signedAxisVector(toothOutwardAxisInMesh);
        const Eigen::Vector3d planningXInMesh =
            planningYInMesh.cross(planningZInMesh).normalized();

        Eigen::Matrix3d meshBasis;
        meshBasis.col(0) = planningXInMesh;
        meshBasis.col(1) = planningYInMesh;
        meshBasis.col(2) = planningZInMesh;
        if (!meshBasis.allFinite() || std::abs(meshBasis.determinant() - 1.0) > 1.0e-12)
        {
            return PlanningResult<AlignmentResult>::failure(
                PlanningErrorCode::InvalidAxisSelection,
                "The selected signed axes do not form a right-handed planning frame.");
        }

        Eigen::Isometry3d planningFromMesh = Eigen::Isometry3d::Identity();
        planningFromMesh.linear() = meshBasis.transpose();
        const PlanningResult<Eigen::AlignedBox3d> rotatedBounds = mesh.bounds(planningFromMesh);
        if (!rotatedBounds)
        {
            return PlanningResult<AlignmentResult>::failure(
                rotatedBounds.error.code, rotatedBounds.error.message);
        }

        planningFromMesh.translation().x() =
            -0.5 * (rotatedBounds.value.min().x() + rotatedBounds.value.max().x());
        planningFromMesh.translation().y() =
            0.5 * motherMaximumDiameterMeters - rotatedBounds.value.max().y();
        planningFromMesh.translation().z() = -rotatedBounds.value.min().z();
        if (!isFiniteTransform(planningFromMesh))
        {
            return PlanningResult<AlignmentResult>::failure(
                PlanningErrorCode::AlignmentFailed,
                "Simulation-block placement produced a non-finite rigid transform.");
        }

        const PlanningResult<Eigen::AlignedBox3d> finalBounds = mesh.bounds(planningFromMesh);
        if (!finalBounds)
        {
            return PlanningResult<AlignmentResult>::failure(
                finalBounds.error.code, finalBounds.error.message);
        }

        AlignmentResult alignment;
        alignment.planningFromMesh = planningFromMesh;
        alignment.automaticBaseline = planningFromMesh;
        alignment.bottomAxisCenterInMesh = planningFromMesh.inverse() * Eigen::Vector3d::Zero();
        alignment.lowConfidence = false;
        alignment.statistics.vertexCount = mesh.vertexCount();
        alignment.statistics.triangleCount = mesh.triangleCount();
        alignment.statistics.boundsMinimum = finalBounds.value.min();
        alignment.statistics.boundsMaximum = finalBounds.value.max();
        alignment.statistics.estimatedAxisInMesh = planningZInMesh;
        alignment.statistics.heightMeters = std::max(
            0.0, finalBounds.value.max().z() - finalBounds.value.min().z());
        alignment.statistics.maximumDiameterMeters = motherMaximumDiameterMeters;
        alignment.statistics.minimumDiameterMeters = 0.0;
        alignment.statistics.axisConfidence = 1.0;

        PlanningResult<AlignmentResult> output =
            PlanningResult<AlignmentResult>::success(std::move(alignment));
        output.diagnostics.push_back(
            "Simulation block placed with rotary axis at +Z, tooth-outward direction at +Y, X midpoint at zero and Z minimum at zero.");
        const double outwardSpan = finalBounds.value.max().y() - finalBounds.value.min().y();
        if (outwardSpan > motherMaximumDiameterMeters)
        {
            std::ostringstream warning;
            warning << "The block's Y span (" << outwardSpan
                    << " m) exceeds the supplied mother-part diameter ("
                    << motherMaximumDiameterMeters << " m); verify the diameter and axis selections.";
            output.diagnostics.push_back(warning.str());
        }
        return output;
    }
}
