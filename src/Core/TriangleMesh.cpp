#include <RotationBodyTrajectoryPlanning/Core/TriangleMesh.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>

namespace smrobot::spray::rotationbody
{
    namespace
    {
        constexpr std::uint64_t fnvOffsetBasis = 14695981039346656037ull;
        constexpr std::uint64_t fnvPrime = 1099511628211ull;

        void hashByte(std::uint64_t& hash, std::uint8_t value) noexcept
        {
            hash ^= static_cast<std::uint64_t>(value);
            hash *= fnvPrime;
        }

        template<typename T>
        void hashValue(std::uint64_t& hash, T value) noexcept
        {
            const auto* bytes = reinterpret_cast<const std::uint8_t*>(&value);
            for (std::size_t i = 0; i < sizeof(T); ++i)
            {
                hashByte(hash, bytes[i]);
            }
        }

        double canonicalZero(double value) noexcept
        {
            return value == 0.0 ? 0.0 : value;
        }
    }

    bool TriangleMesh::empty() const noexcept
    {
        return positions.empty() || triangles.empty();
    }

    std::size_t TriangleMesh::vertexCount() const noexcept
    {
        return positions.size();
    }

    std::size_t TriangleMesh::triangleCount() const noexcept
    {
        return triangles.size();
    }

    PlanningResult<void> TriangleMesh::validate(double relativeAreaTolerance) const
    {
        if (empty())
        {
            return PlanningResult<void>::failure(
                PlanningErrorCode::EmptyMesh,
                "Triangle mesh must contain positions and triangles.");
        }
        if (!std::isfinite(relativeAreaTolerance) || relativeAreaTolerance < 0.0)
        {
            return PlanningResult<void>::failure(
                PlanningErrorCode::InvalidArgument,
                "Relative area tolerance must be finite and non-negative.");
        }
        for (const Eigen::Vector3d& position : positions)
        {
            if (!position.allFinite())
            {
                return PlanningResult<void>::failure(
                    PlanningErrorCode::NonFiniteGeometry,
                    "Triangle mesh contains a non-finite position.");
            }
        }

        Eigen::AlignedBox3d meshBounds;
        for (const Eigen::Vector3d& position : positions)
        {
            meshBounds.extend(position);
        }
        const double diagonal = meshBounds.diagonal().norm();
        if (!std::isfinite(diagonal) || diagonal <= std::numeric_limits<double>::epsilon())
        {
            return PlanningResult<void>::failure(
                PlanningErrorCode::DegenerateGeometry,
                "Triangle mesh has no measurable extent.");
        }
        const double areaToleranceSquared =
            std::pow(std::max(diagonal * diagonal * relativeAreaTolerance, 1.0e-30), 2.0);

        for (const Triangle& triangle : triangles)
        {
            for (int corner = 0; corner < 3; ++corner)
            {
                if (triangle[corner] < 0 ||
                    static_cast<std::size_t>(triangle[corner]) >= positions.size())
                {
                    return PlanningResult<void>::failure(
                        PlanningErrorCode::InvalidTriangleIndex,
                        "Triangle mesh contains an out-of-range vertex index.");
                }
            }
            if (triangle[0] == triangle[1] || triangle[1] == triangle[2] || triangle[2] == triangle[0])
            {
                return PlanningResult<void>::failure(
                    PlanningErrorCode::DegenerateGeometry,
                    "Triangle mesh contains a triangle with repeated indices.");
            }

            const Eigen::Vector3d cross =
                (positions[triangle[1]] - positions[triangle[0]]).cross(
                    positions[triangle[2]] - positions[triangle[0]]);
            if (!cross.allFinite() || cross.squaredNorm() <= areaToleranceSquared)
            {
                return PlanningResult<void>::failure(
                    PlanningErrorCode::DegenerateGeometry,
                    "Triangle mesh contains a zero-area triangle.");
            }
        }
        return PlanningResult<void>::success();
    }

    PlanningResult<Eigen::AlignedBox3d> TriangleMesh::bounds(const Eigen::Isometry3d& transform) const
    {
        if (positions.empty())
        {
            return PlanningResult<Eigen::AlignedBox3d>::failure(
                PlanningErrorCode::EmptyMesh,
                "Cannot compute bounds for an empty mesh.");
        }
        if (!transform.matrix().allFinite())
        {
            return PlanningResult<Eigen::AlignedBox3d>::failure(
                PlanningErrorCode::NonFiniteGeometry,
                "Cannot compute bounds with a non-finite transform.");
        }

        Eigen::AlignedBox3d result;
        for (const Eigen::Vector3d& position : positions)
        {
            if (!position.allFinite())
            {
                return PlanningResult<Eigen::AlignedBox3d>::failure(
                    PlanningErrorCode::NonFiniteGeometry,
                    "Cannot compute bounds for non-finite geometry.");
            }
            result.extend(transform * position);
        }
        return PlanningResult<Eigen::AlignedBox3d>::success(result);
    }

    double TriangleMesh::scale() const
    {
        const auto result = bounds();
        return result ? result.value.diagonal().norm() : 0.0;
    }

    std::uint64_t TriangleMesh::stableFingerprint() const noexcept
    {
        std::uint64_t hash = fnvOffsetBasis;
        hashValue(hash, static_cast<std::uint64_t>(positions.size()));
        hashValue(hash, static_cast<std::uint64_t>(triangles.size()));
        for (const Eigen::Vector3d& position : positions)
        {
            for (int component = 0; component < 3; ++component)
            {
                const double value = canonicalZero(position[component]);
                std::uint64_t bits = 0;
                static_assert(sizeof(bits) == sizeof(value), "Unexpected double width.");
                std::memcpy(&bits, &value, sizeof(bits));
                hashValue(hash, bits);
            }
        }
        for (const Triangle& triangle : triangles)
        {
            for (int corner = 0; corner < 3; ++corner)
            {
                hashValue(hash, static_cast<std::int32_t>(triangle[corner]));
            }
        }
        return hash;
    }

    PlanningResult<TriangleMesh> TriangleMesh::transformed(const Eigen::Isometry3d& transform) const
    {
        if (!transform.matrix().allFinite())
        {
            return PlanningResult<TriangleMesh>::failure(
                PlanningErrorCode::NonFiniteGeometry,
                "Cannot transform a mesh with a non-finite transform.");
        }
        TriangleMesh result = *this;
        for (Eigen::Vector3d& position : result.positions)
        {
            position = transform * position;
        }
        const auto validation = result.validate();
        if (!validation)
        {
            return PlanningResult<TriangleMesh>::failure(validation.error.code, validation.error.message);
        }
        return PlanningResult<TriangleMesh>::success(std::move(result));
    }
}
