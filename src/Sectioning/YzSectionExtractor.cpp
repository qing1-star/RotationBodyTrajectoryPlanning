#include <RotationBodyTrajectoryPlanning/Sectioning/YzSectionExtractor.h>

#include <RotationBodyTrajectoryPlanning/Core/TransformUtils.h>
#include <RotationBodyTrajectoryPlanning/Sectioning/ContourTopology.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace smrobot::spray::rotationbody
{
    namespace
    {
        struct QuantizedYz
        {
            std::int64_t y{ 0 };
            std::int64_t z{ 0 };

            bool operator<(const QuantizedYz& other) const noexcept
            {
                return y < other.y || (y == other.y && z < other.z);
            }
        };

        struct QuantizedSegment
        {
            QuantizedYz first;
            QuantizedYz second;

            bool operator<(const QuantizedSegment& other) const noexcept
            {
                if (first < other.first)
                {
                    return true;
                }
                if (other.first < first)
                {
                    return false;
                }
                return second < other.second;
            }
        };

        struct CoplanarEdge
        {
            YzSegment segment;
            std::size_t occurrenceCount{ 0 };
        };

        bool quantizeCoordinate(double value, double tolerance, std::int64_t& output) noexcept
        {
            const long double scaled = std::round(
                static_cast<long double>(value) / static_cast<long double>(tolerance));
            const auto minimum = static_cast<long double>(std::numeric_limits<std::int64_t>::min() + 1);
            const auto maximum = static_cast<long double>(std::numeric_limits<std::int64_t>::max() - 1);
            if (!std::isfinite(scaled) || scaled < minimum || scaled > maximum)
            {
                return false;
            }
            output = static_cast<std::int64_t>(scaled);
            return true;
        }

        bool makeSegmentKey(
            const YzSegment& segment,
            double tolerance,
            QuantizedSegment& output) noexcept
        {
            QuantizedYz first;
            QuantizedYz second;
            if (!quantizeCoordinate(segment.first.x(), tolerance, first.y)
                || !quantizeCoordinate(segment.first.y(), tolerance, first.z)
                || !quantizeCoordinate(segment.second.x(), tolerance, second.y)
                || !quantizeCoordinate(segment.second.y(), tolerance, second.z))
            {
                return false;
            }
            if (second < first)
            {
                std::swap(first, second);
            }
            output = { first, second };
            return true;
        }

        int planeSide(double signedDistance, double tolerance) noexcept
        {
            if (signedDistance > tolerance)
            {
                return 1;
            }
            if (signedDistance < -tolerance)
            {
                return -1;
            }
            return 0;
        }

        bool appendUniqueIntersection(
            std::vector<Eigen::Vector3d>& intersections,
            const Eigen::Vector3d& point,
            double tolerance)
        {
            const double toleranceSquared = tolerance * tolerance;
            const auto duplicate = std::find_if(intersections.begin(), intersections.end(),
                [&](const Eigen::Vector3d& existing)
                {
                    return (existing - point).squaredNorm() <= toleranceSquared;
                });
            if (duplicate != intersections.end())
            {
                return false;
            }
            intersections.push_back(point);
            return true;
        }

        bool clipToPositiveY(YzSegment& segment, double tolerance) noexcept
        {
            const bool firstOutside = segment.first.x() < -tolerance;
            const bool secondOutside = segment.second.x() < -tolerance;
            if (firstOutside && secondOutside)
            {
                return false;
            }

            if (firstOutside != secondOutside)
            {
                const Eigen::Vector2d delta = segment.second - segment.first;
                if (std::abs(delta.x()) <= std::numeric_limits<double>::epsilon())
                {
                    return false;
                }
                Eigen::Vector2d crossing = segment.first
                    + (-segment.first.x() / delta.x()) * delta;
                crossing.x() = 0.0;
                if (firstOutside)
                {
                    segment.first = crossing;
                }
                else
                {
                    segment.second = crossing;
                }
            }

            if (std::abs(segment.first.x()) <= tolerance)
            {
                segment.first.x() = 0.0;
            }
            if (std::abs(segment.second.x()) <= tolerance)
            {
                segment.second.x() = 0.0;
            }
            return (segment.second - segment.first).squaredNorm() > tolerance * tolerance;
        }

        template<typename T>
        PlanningResult<T> failureWithDiagnostics(
            PlanningErrorCode code,
            std::string message,
            std::vector<std::string> diagnostics = {})
        {
            PlanningResult<T> result = PlanningResult<T>::failure(code, std::move(message));
            result.diagnostics = std::move(diagnostics);
            return result;
        }
    }

    PlanningResult<YzSectionResult> YzSectionExtractor::extract(
        const TriangleMesh& mesh,
        const Eigen::Isometry3d& planningFromMesh,
        const YzSectionOptions& options)
    {
        if (!std::isfinite(options.relativePlaneTolerance)
            || !std::isfinite(options.relativeWeldTolerance)
            || !std::isfinite(options.minimumToleranceMeters)
            || options.relativePlaneTolerance <= 0.0
            || options.relativePlaneTolerance > 1.0e-4
            || options.relativeWeldTolerance <= 0.0
            || options.relativeWeldTolerance > 1.0e-3
            || options.minimumToleranceMeters < 0.0)
        {
            return PlanningResult<YzSectionResult>::failure(
                PlanningErrorCode::InvalidArgument,
                "Section tolerances must be finite, non-negative, and within their supported relative range.");
        }
        if (!isFiniteTransform(planningFromMesh))
        {
            return PlanningResult<YzSectionResult>::failure(
                PlanningErrorCode::InvalidArgument,
                "The planning-from-mesh transform contains non-finite values.");
        }

        const PlanningResult<void> validation = mesh.validate();
        if (!validation)
        {
            return failureWithDiagnostics<YzSectionResult>(
                validation.error.code,
                validation.error.message,
                validation.diagnostics);
        }

        std::vector<Eigen::Vector3d> transformedPositions;
        transformedPositions.reserve(mesh.positions.size());
        Eigen::Vector3d minimum = Eigen::Vector3d::Constant(
            std::numeric_limits<double>::infinity());
        Eigen::Vector3d maximum = Eigen::Vector3d::Constant(
            -std::numeric_limits<double>::infinity());
        for (const auto& position : mesh.positions)
        {
            const Eigen::Vector3d transformed = planningFromMesh * position;
            if (!transformed.allFinite())
            {
                return PlanningResult<YzSectionResult>::failure(
                    PlanningErrorCode::NonFiniteGeometry,
                    "Transforming the mesh produced a non-finite vertex.");
            }
            transformedPositions.push_back(transformed);
            minimum = minimum.cwiseMin(transformed);
            maximum = maximum.cwiseMax(transformed);
        }

        const double meshScale = (maximum - minimum).norm();
        if (!std::isfinite(meshScale) || meshScale <= std::numeric_limits<double>::min())
        {
            return PlanningResult<YzSectionResult>::failure(
                PlanningErrorCode::DegenerateGeometry,
                "The transformed mesh has no measurable three-dimensional scale.");
        }

        const double coordinateMagnitude = std::max(
            minimum.cwiseAbs().maxCoeff(),
            maximum.cwiseAbs().maxCoeff());
        const double machineTolerance = 128.0 * std::numeric_limits<double>::epsilon()
            * std::max(meshScale, coordinateMagnitude);
        const double absoluteFloor = std::min(
            options.minimumToleranceMeters,
            meshScale * 1.0e-6);
        const double relativePlaneTolerance = std::clamp(
            options.relativePlaneTolerance,
            1.0e-12,
            1.0e-4);
        double planeTolerance = std::max({
            meshScale * relativePlaneTolerance,
            machineTolerance,
            absoluteFloor });
        planeTolerance = std::min(
            planeTolerance,
            std::max(meshScale * 1.0e-4, machineTolerance));
        const double weldTolerance = std::max({
            meshScale * options.relativeWeldTolerance,
            4.0 * planeTolerance,
            2.0 * machineTolerance });

        YzSectionResult output;
        output.meshScaleMeters = meshScale;
        output.planeToleranceMeters = planeTolerance;
        output.weldToleranceMeters = weldTolerance;

        std::vector<YzSegment> rawSegments;
        rawSegments.reserve(mesh.triangles.size());
        std::map<QuantizedSegment, CoplanarEdge> coplanarEdges;

        for (const auto& triangleIndices : mesh.triangles)
        {
            ++output.processedTriangleCount;
            std::array<Eigen::Vector3d, 3> triangle{
                transformedPositions[static_cast<std::size_t>(triangleIndices.x())],
                transformedPositions[static_cast<std::size_t>(triangleIndices.y())],
                transformedPositions[static_cast<std::size_t>(triangleIndices.z())] };
            std::array<double, 3> signedDistances{
                triangle[0].x(), triangle[1].x(), triangle[2].x() };
            std::array<int, 3> sides{
                planeSide(signedDistances[0], planeTolerance),
                planeSide(signedDistances[1], planeTolerance),
                planeSide(signedDistances[2], planeTolerance) };

            if (sides[0] == 0 && sides[1] == 0 && sides[2] == 0)
            {
                ++output.coplanarTriangleCount;
                for (std::size_t edgeIndex = 0; edgeIndex < 3; ++edgeIndex)
                {
                    const std::size_t next = (edgeIndex + 1) % 3;
                    YzSegment segment{
                        Eigen::Vector2d(triangle[edgeIndex].y(), triangle[edgeIndex].z()),
                        Eigen::Vector2d(triangle[next].y(), triangle[next].z()) };
                    if ((segment.second - segment.first).squaredNorm()
                        <= weldTolerance * weldTolerance)
                    {
                        ++output.discardedSegmentCount;
                        continue;
                    }
                    QuantizedSegment key;
                    if (!makeSegmentKey(segment, weldTolerance, key))
                    {
                        return PlanningResult<YzSectionResult>::failure(
                            PlanningErrorCode::DegenerateGeometry,
                            "Coplanar section coordinates exceed the safe quantization range.");
                    }
                    auto& edge = coplanarEdges[key];
                    if (edge.occurrenceCount == 0)
                    {
                        edge.segment = segment;
                    }
                    ++edge.occurrenceCount;
                }
                continue;
            }

            std::vector<Eigen::Vector3d> intersections;
            intersections.reserve(3);
            for (std::size_t vertexIndex = 0; vertexIndex < 3; ++vertexIndex)
            {
                if (sides[vertexIndex] == 0)
                {
                    Eigen::Vector3d point = triangle[vertexIndex];
                    point.x() = 0.0;
                    appendUniqueIntersection(intersections, point, planeTolerance);
                }
            }
            for (std::size_t edgeIndex = 0; edgeIndex < 3; ++edgeIndex)
            {
                const std::size_t next = (edgeIndex + 1) % 3;
                if (sides[edgeIndex] * sides[next] >= 0)
                {
                    continue;
                }
                const double denominator = signedDistances[edgeIndex] - signedDistances[next];
                if (std::abs(denominator) <= std::numeric_limits<double>::min())
                {
                    continue;
                }
                const double interpolation = signedDistances[edgeIndex] / denominator;
                Eigen::Vector3d point = triangle[edgeIndex]
                    + interpolation * (triangle[next] - triangle[edgeIndex]);
                point.x() = 0.0;
                appendUniqueIntersection(intersections, point, planeTolerance);
            }

            if (intersections.size() < 2)
            {
                continue;
            }

            std::size_t firstIndex = 0;
            std::size_t secondIndex = 1;
            double longestDistanceSquared = 0.0;
            for (std::size_t first = 0; first < intersections.size(); ++first)
            {
                for (std::size_t second = first + 1; second < intersections.size(); ++second)
                {
                    const double distanceSquared =
                        (intersections[first] - intersections[second]).squaredNorm();
                    if (distanceSquared > longestDistanceSquared)
                    {
                        longestDistanceSquared = distanceSquared;
                        firstIndex = first;
                        secondIndex = second;
                    }
                }
            }
            if (longestDistanceSquared <= weldTolerance * weldTolerance)
            {
                ++output.discardedSegmentCount;
                continue;
            }
            rawSegments.push_back({
                Eigen::Vector2d(intersections[firstIndex].y(), intersections[firstIndex].z()),
                Eigen::Vector2d(intersections[secondIndex].y(), intersections[secondIndex].z()) });
        }

        std::size_t cancelledCoplanarEdges = 0;
        for (const auto& entry : coplanarEdges)
        {
            if (entry.second.occurrenceCount % 2 == 1)
            {
                rawSegments.push_back(entry.second.segment);
            }
            else
            {
                ++cancelledCoplanarEdges;
            }
        }
        output.rawIntersectionSegmentCount = rawSegments.size();

        std::vector<YzSegment> positiveYSegments;
        positiveYSegments.reserve(rawSegments.size());
        for (YzSegment segment : rawSegments)
        {
            if (clipToPositiveY(segment, weldTolerance))
            {
                positiveYSegments.push_back(std::move(segment));
            }
            else
            {
                ++output.discardedSegmentCount;
            }
        }
        output.positiveYSegmentCount = positiveYSegments.size();
        if (positiveYSegments.empty())
        {
            return PlanningResult<YzSectionResult>::failure(
                PlanningErrorCode::NoContour,
                "The X=0 plane produced no non-degenerate section segment in Y >= 0.");
        }

        PlanningResult<ContourTopologyResult> topology = ContourTopology::buildContours(
            positiveYSegments,
            weldTolerance,
            0.0);
        if (!topology)
        {
            return failureWithDiagnostics<YzSectionResult>(
                topology.error.code,
                topology.error.message,
                topology.diagnostics);
        }
        output.discardedSegmentCount += topology.value.discardedDegenerateSegmentCount
            + topology.value.discardedDuplicateSegmentCount;
        output.branchNodeCount = topology.value.branchNodeCount;
        output.candidateContours = std::move(topology.value.contours);

        PlanningResult<std::size_t> target = ContourTopology::selectTargetOuterContour(
            output.candidateContours,
            weldTolerance);
        if (!target)
        {
            std::vector<std::string> diagnostics = std::move(topology.diagnostics);
            diagnostics.insert(
                diagnostics.end(),
                target.diagnostics.begin(),
                target.diagnostics.end());
            return failureWithDiagnostics<YzSectionResult>(
                target.error.code,
                target.error.message,
                std::move(diagnostics));
        }

        output.targetContour = output.candidateContours[target.value];
        std::ostringstream summary;
        summary << std::scientific << std::setprecision(3)
            << "YZ section: scale=" << meshScale << " m, plane tolerance="
            << planeTolerance << " m, weld tolerance=" << weldTolerance << " m; "
            << output.processedTriangleCount << " triangles, "
            << output.rawIntersectionSegmentCount << " raw segments, "
            << output.positiveYSegmentCount << " positive-Y segments, "
            << output.candidateContours.size() << " contours.";
        output.targetContour.diagnostics.push_back(summary.str());
        if (cancelledCoplanarEdges > 0)
        {
            std::ostringstream coplanarSummary;
            coplanarSummary << "Removed " << cancelledCoplanarEdges
                << " even-occurrence internal edges from coplanar triangles.";
            output.targetContour.diagnostics.push_back(coplanarSummary.str());
        }
        output.targetContour.diagnostics.insert(
            output.targetContour.diagnostics.end(),
            target.diagnostics.begin(),
            target.diagnostics.end());

        PlanningResult<YzSectionResult> result =
            PlanningResult<YzSectionResult>::success(std::move(output));
        result.diagnostics = std::move(topology.diagnostics);
        result.diagnostics.insert(
            result.diagnostics.end(),
            target.diagnostics.begin(),
            target.diagnostics.end());
        result.diagnostics.push_back(summary.str());
        return result;
    }

    PlanningResult<SectionContour> YzSectionExtractor::extractTargetContour(
        const TriangleMesh& mesh,
        const Eigen::Isometry3d& planningFromMesh,
        const YzSectionOptions& options)
    {
        PlanningResult<YzSectionResult> section = extract(mesh, planningFromMesh, options);
        PlanningResult<SectionContour> result;
        result.error = section.error;
        result.diagnostics = std::move(section.diagnostics);
        if (section)
        {
            result.value = std::move(section.value.targetContour);
        }
        return result;
    }
}
