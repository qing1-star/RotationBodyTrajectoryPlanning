#pragma once

#include <RotationBodyTrajectoryPlanning/Core/PlanningTypes.h>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <cstddef>
#include <cstdint>
#include <vector>

namespace smrobot::spray::rotationbody
{
    class TriangleMesh
    {
    public:
        using Triangle = Eigen::Vector3i;

        std::vector<Eigen::Vector3d> positions;
        std::vector<Triangle> triangles;

        bool empty() const noexcept;
        std::size_t vertexCount() const noexcept;
        std::size_t triangleCount() const noexcept;

        PlanningResult<void> validate(double relativeAreaTolerance = 1.0e-14) const;
        PlanningResult<Eigen::AlignedBox3d> bounds(
            const Eigen::Isometry3d& transform = Eigen::Isometry3d::Identity()) const;
        double scale() const;
        std::uint64_t stableFingerprint() const noexcept;
        PlanningResult<TriangleMesh> transformed(const Eigen::Isometry3d& transform) const;
    };
}
