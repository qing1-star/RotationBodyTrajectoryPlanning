#include <RotationBodyTrajectoryPlanning/RegionPlanning/ToothRegionRecognizer.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numeric>
#include <optional>
#include <sstream>
#include <utility>

namespace smrobot::spray::rotationbody
{
    namespace
    {
        struct SampledContour
        {
            std::vector<Eigen::Vector2d> points;
            std::vector<Eigen::Vector2d> smoothed;
            double totalLength{ 0.0 };
            bool closed{ false };

            std::size_t segmentCount() const noexcept
            {
                if (points.size() < 2)
                {
                    return 0;
                }
                return closed ? points.size() : points.size() - 1;
            }
        };

        std::size_t wrappedIndex(long long index, std::size_t size) noexcept
        {
            const long long signedSize = static_cast<long long>(size);
            long long result = index % signedSize;
            if (result < 0)
            {
                result += signedSize;
            }
            return static_cast<std::size_t>(result);
        }

        template <typename T>
        void rotateLeft(std::vector<T>& values, std::size_t offset)
        {
            if (values.empty())
            {
                return;
            }
            offset %= values.size();
            std::rotate(
                values.begin(),
                values.begin() + static_cast<std::ptrdiff_t>(offset),
                values.end());
        }

        template <typename T>
        void undoRotateLeft(std::vector<T>& values, std::size_t offset)
        {
            if (values.empty())
            {
                return;
            }
            offset %= values.size();
            if (offset != 0)
            {
                std::rotate(
                    values.begin(),
                    values.end() - static_cast<std::ptrdiff_t>(offset),
                    values.end());
            }
        }

        std::size_t canonicalizeClosedContour(SectionContour& contour)
        {
            if (!contour.closed || contour.pointsYz.empty())
            {
                return 0;
            }
            if (std::any_of(
                    contour.pointsYz.begin(),
                    contour.pointsYz.end(),
                    [](const Eigen::Vector2d& point)
                    {
                        return !point.allFinite();
                    }))
            {
                return 0;
            }
            const auto canonical = std::min_element(
                contour.pointsYz.begin(),
                contour.pointsYz.end(),
                [](const Eigen::Vector2d& first, const Eigen::Vector2d& second)
                {
                    if (first.y() != second.y())
                    {
                        return first.y() < second.y();
                    }
                    return first.x() < second.x();
                });
            const std::size_t offset = static_cast<std::size_t>(
                std::distance(contour.pointsYz.begin(), canonical));
            rotateLeft(contour.pointsYz, offset);
            if (contour.points3d.size() == contour.pointsYz.size())
            {
                rotateLeft(contour.points3d, offset);
            }
            contour.cumulativeArcLength.clear();
            return offset;
        }

        Eigen::Vector2d interpolateAtDistance(
            const SectionContour& contour,
            const std::vector<double>& cumulative,
            double totalLength,
            double distance)
        {
            if (contour.closed)
            {
                distance = std::fmod(distance, totalLength);
                if (distance < 0.0)
                {
                    distance += totalLength;
                }
            }
            else
            {
                distance = std::clamp(distance, 0.0, totalLength);
            }

            const auto upper = std::upper_bound(cumulative.begin(), cumulative.end(), distance);
            std::size_t segment = upper == cumulative.begin()
                ? 0
                : static_cast<std::size_t>(std::distance(cumulative.begin(), upper) - 1);
            const std::size_t segmentCount = contour.segmentCount();
            segment = std::min(segment, segmentCount - 1);
            const std::size_t end = (segment + 1) % contour.pointsYz.size();
            const double segmentStart = cumulative[segment];
            const double segmentEnd = cumulative[segment + 1];
            const double denominator = segmentEnd - segmentStart;
            const double factor = denominator > 0.0
                ? std::clamp((distance - segmentStart) / denominator, 0.0, 1.0)
                : 0.0;
            return contour.pointsYz[segment] +
                factor * (contour.pointsYz[end] - contour.pointsYz[segment]);
        }

        PlanningResult<SampledContour> resample(
            const SectionContour& contour,
            const ToothRecognitionOptions& options)
        {
            if (contour.segmentCount() < 3)
            {
                return PlanningResult<SampledContour>::failure(
                    PlanningErrorCode::NoContour,
                    "Tooth recognition requires at least three contour segments.");
            }
            if (options.minimumResampleCount < 8 ||
                options.maximumResampleCount < options.minimumResampleCount ||
                !std::isfinite(options.smoothingWindowFraction) ||
                options.smoothingWindowFraction < 0.0 ||
                !std::isfinite(options.minimumRadialProminenceFraction) ||
                options.minimumRadialProminenceFraction <= 0.0 ||
                !std::isfinite(options.axialSurfaceThreshold) ||
                options.axialSurfaceThreshold <= 0.5 || options.axialSurfaceThreshold >= 1.0 ||
                !std::isfinite(options.minimumConfidence) ||
                options.minimumConfidence < 0.0 || options.minimumConfidence > 1.0 ||
                !std::isfinite(options.transitionExtentFraction) ||
                options.transitionExtentFraction < 0.0)
            {
                return PlanningResult<SampledContour>::failure(
                    PlanningErrorCode::InvalidArgument,
                    "Tooth recognition options are outside their supported ranges.");
            }
            for (const Eigen::Vector2d& point : contour.pointsYz)
            {
                if (!point.allFinite())
                {
                    return PlanningResult<SampledContour>::failure(
                        PlanningErrorCode::NonFiniteGeometry,
                        "Section contour contains a non-finite point.");
                }
            }

            std::vector<double> cumulative(contour.segmentCount() + 1, 0.0);
            std::vector<double> originalSegmentLengths;
            originalSegmentLengths.reserve(contour.segmentCount());
            for (std::size_t segment = 0; segment < contour.segmentCount(); ++segment)
            {
                const Eigen::Vector2d& start = contour.pointsYz[segment];
                const Eigen::Vector2d& end =
                    contour.pointsYz[(segment + 1) % contour.pointsYz.size()];
                const double length = (end - start).norm();
                if (!std::isfinite(length) || length <= 0.0)
                {
                    return PlanningResult<SampledContour>::failure(
                        PlanningErrorCode::DegenerateGeometry,
                        "Section contour contains a zero-length segment.");
                }
                originalSegmentLengths.push_back(length);
                cumulative[segment + 1] = cumulative[segment] + length;
            }
            const double totalLength = cumulative.back();
            if (!std::isfinite(totalLength) || totalLength <= 0.0)
            {
                return PlanningResult<SampledContour>::failure(
                    PlanningErrorCode::DegenerateGeometry,
                    "Section contour has no measurable length.");
            }

            std::sort(originalSegmentLengths.begin(), originalSegmentLengths.end());
            const std::size_t middleSegment = originalSegmentLengths.size() / 2;
            const double medianOriginalSegmentLength =
                originalSegmentLengths.size() % 2 == 0
                    ? 0.5 * (originalSegmentLengths[middleSegment - 1] +
                        originalSegmentLengths[middleSegment])
                    : originalSegmentLengths[middleSegment];
            const std::size_t resolutionCount = static_cast<std::size_t>(
                std::ceil(totalLength /
                    std::max(0.4 * medianOriginalSegmentLength, 1.0e-15)));
            const std::size_t desiredCount = std::max({
                options.minimumResampleCount,
                contour.pointsYz.size() * 8,
                resolutionCount });
            const std::size_t sampleCount = std::clamp(
                desiredCount,
                options.minimumResampleCount,
                options.maximumResampleCount);

            SampledContour result;
            result.totalLength = totalLength;
            result.closed = contour.closed;
            result.points.reserve(sampleCount);
            for (std::size_t sample = 0; sample < sampleCount; ++sample)
            {
                const double denominator = contour.closed
                    ? static_cast<double>(sampleCount)
                    : static_cast<double>(sampleCount - 1);
                const double distance = totalLength * static_cast<double>(sample) / denominator;
                result.points.push_back(interpolateAtDistance(contour, cumulative, totalLength, distance));
            }

            result.smoothed.resize(sampleCount);
            const double sampleSpacing = totalLength /
                static_cast<double>(contour.closed ? sampleCount : sampleCount - 1);
            const double requestedSmoothingLength =
                options.smoothingWindowFraction * totalLength;
            const double smoothingLength = sampleCount >= 512
                ? std::min(
                    requestedSmoothingLength,
                    medianOriginalSegmentLength)
                : requestedSmoothingLength;
            const std::size_t halfWindow = std::max<std::size_t>(
                1,
                static_cast<std::size_t>(
                    std::ceil(smoothingLength / sampleSpacing)));
            for (std::size_t sample = 0; sample < sampleCount; ++sample)
            {
                Eigen::Vector2d weighted = Eigen::Vector2d::Zero();
                double weightSum = 0.0;
                for (long long offset = -static_cast<long long>(halfWindow);
                     offset <= static_cast<long long>(halfWindow);
                     ++offset)
                {
                    long long candidate = static_cast<long long>(sample) + offset;
                    if (!contour.closed &&
                        (candidate < 0 || candidate >= static_cast<long long>(sampleCount)))
                    {
                        continue;
                    }
                    const std::size_t index = contour.closed
                        ? wrappedIndex(candidate, sampleCount)
                        : static_cast<std::size_t>(candidate);
                    const double weight =
                        static_cast<double>(halfWindow + 1) - std::abs(static_cast<double>(offset));
                    weighted += weight * result.points[index];
                    weightSum += weight;
                }
                result.smoothed[sample] = weighted / weightSum;
            }
            return PlanningResult<SampledContour>::success(std::move(result));
        }

        bool isToothLabel(RegionLabel label) noexcept
        {
            return label == RegionLabel::ToothTop ||
                label == RegionLabel::ToothWall ||
                label == RegionLabel::ToothBottom;
        }

        std::size_t stableClosedClassificationSeam(
            const std::vector<RegionLabel>& labels,
            bool closed)
        {
            if (!closed || labels.empty())
            {
                return 0;
            }
            std::vector<std::size_t> toothIndices;
            toothIndices.reserve(labels.size());
            for (std::size_t index = 0; index < labels.size(); ++index)
            {
                if (isToothLabel(labels[index]))
                {
                    toothIndices.push_back(index);
                }
            }
            if (toothIndices.empty() || toothIndices.size() == labels.size())
            {
                return 0;
            }

            std::size_t largestGap = 0;
            std::size_t seam = 0;
            for (std::size_t position = 0; position < toothIndices.size(); ++position)
            {
                const std::size_t current = toothIndices[position];
                const std::size_t next = position + 1 < toothIndices.size()
                    ? toothIndices[position + 1]
                    : toothIndices.front() + labels.size();
                const std::size_t gap = next - current - 1;
                if (gap > largestGap)
                {
                    largestGap = gap;
                    seam = (current + 1) % labels.size();
                }
            }
            return seam;
        }

        bool isOuterToothBoundary(RegionLabel label) noexcept
        {
            return label == RegionLabel::ToothWall ||
                label == RegionLabel::ToothBottom;
        }

        struct ToothCluster
        {
            std::size_t start{ 0 };
            std::size_t end{ 0 };
            bool coversClosedContour{ false };
        };

        struct TransitionDetectionSummary
        {
            std::size_t clusterCount{ 0 };
            std::size_t qualifiedBoundaryCount{ 0 };
            std::size_t markedSegmentCount{ 0 };
            std::size_t outerTopRunCount{ 0 };
        };

        std::vector<ToothCluster> buildToothClusters(
            const std::vector<RegionLabel>& labels,
            bool closed,
            std::size_t maximumInternalGap)
        {
            std::vector<std::size_t> toothIndices;
            toothIndices.reserve(labels.size());
            for (std::size_t index = 0; index < labels.size(); ++index)
            {
                if (isToothLabel(labels[index]))
                {
                    toothIndices.push_back(index);
                }
            }
            if (toothIndices.empty())
            {
                return {};
            }

            if (!closed)
            {
                std::vector<ToothCluster> clusters;
                std::size_t start = toothIndices.front();
                std::size_t previous = start;
                for (std::size_t position = 1; position < toothIndices.size(); ++position)
                {
                    const std::size_t current = toothIndices[position];
                    if (current - previous - 1 > maximumInternalGap)
                    {
                        clusters.push_back({ start, previous, false });
                        start = current;
                    }
                    previous = current;
                }
                clusters.push_back({ start, previous, false });
                return clusters;
            }

            std::size_t largestGap = 0;
            std::size_t largestGapAfter = 0;
            for (std::size_t position = 0; position < toothIndices.size(); ++position)
            {
                const std::size_t current = toothIndices[position];
                const std::size_t next = position + 1 < toothIndices.size()
                    ? toothIndices[position + 1]
                    : toothIndices.front() + labels.size();
                const std::size_t gap = next - current - 1;
                if (gap > largestGap)
                {
                    largestGap = gap;
                    largestGapAfter = position;
                }
            }

            if (largestGap <= maximumInternalGap)
            {
                return { { toothIndices.front(), toothIndices.front() + labels.size() - 1, true } };
            }

            std::vector<std::size_t> ordered;
            ordered.reserve(toothIndices.size());
            const std::size_t firstPosition = (largestGapAfter + 1) % toothIndices.size();
            for (std::size_t offset = 0; offset < toothIndices.size(); ++offset)
            {
                const std::size_t position = (firstPosition + offset) % toothIndices.size();
                std::size_t index = toothIndices[position];
                if (!ordered.empty() && index <= ordered.back())
                {
                    index += labels.size();
                }
                ordered.push_back(index);
            }

            std::vector<ToothCluster> clusters;
            std::size_t start = ordered.front();
            std::size_t previous = start;
            for (std::size_t position = 1; position < ordered.size(); ++position)
            {
                const std::size_t current = ordered[position];
                if (current - previous - 1 > maximumInternalGap)
                {
                    clusters.push_back({ start, previous, false });
                    start = current;
                }
                previous = current;
            }
            clusters.push_back({ start, previous, false });
            return clusters;
        }

        double tangentAngle(
            const Eigen::Vector2d& first,
            const Eigen::Vector2d& second) noexcept
        {
            return std::acos(std::clamp(first.dot(second), -1.0, 1.0));
        }

        std::vector<Eigen::Vector2d> segmentTangents(const SampledContour& sampled)
        {
            std::vector<Eigen::Vector2d> tangents(sampled.segmentCount(), Eigen::Vector2d::Zero());
            for (std::size_t segment = 0; segment < sampled.segmentCount(); ++segment)
            {
                const std::size_t end = (segment + 1) % sampled.smoothed.size();
                const Eigen::Vector2d delta = sampled.smoothed[end] - sampled.smoothed[segment];
                const double length = delta.norm();
                if (length > 1.0e-15)
                {
                    tangents[segment] = delta / length;
                }
            }
            return tangents;
        }

        std::vector<double> localTangentChanges(
            const std::vector<Eigen::Vector2d>& tangents,
            bool closed)
        {
            std::vector<double> changes(tangents.size(), 0.0);
            for (std::size_t index = 0; index < tangents.size(); ++index)
            {
                if (tangents[index].squaredNorm() <= 0.0)
                {
                    changes[index] = std::numeric_limits<double>::infinity();
                    continue;
                }
                if (closed || index > 0)
                {
                    const std::size_t previous = index > 0 ? index - 1 : tangents.size() - 1;
                    if (tangents[previous].squaredNorm() > 0.0)
                    {
                        changes[index] = std::max(
                            changes[index],
                            tangentAngle(tangents[index], tangents[previous]));
                    }
                }
                if (closed || index + 1 < tangents.size())
                {
                    const std::size_t next = index + 1 < tangents.size() ? index + 1 : 0;
                    if (tangents[next].squaredNorm() > 0.0)
                    {
                        changes[index] = std::max(
                            changes[index],
                            tangentAngle(tangents[index], tangents[next]));
                    }
                }
            }
            return changes;
        }

        bool stableTangentRun(
            const std::vector<std::size_t>& outwardIndices,
            std::size_t start,
            std::size_t count,
            const std::vector<Eigen::Vector2d>& tangents,
            const std::vector<double>& localChanges,
            double changeThreshold,
            Eigen::Vector2d* averageTangent = nullptr)
        {
            if (start + count > outwardIndices.size() || count == 0)
            {
                return false;
            }

            const Eigen::Vector2d reference = tangents[outwardIndices[start]];
            if (reference.squaredNorm() <= 0.0)
            {
                return false;
            }
            Eigen::Vector2d sum = Eigen::Vector2d::Zero();
            for (std::size_t offset = 0; offset < count; ++offset)
            {
                const std::size_t index = outwardIndices[start + offset];
                if (tangents[index].squaredNorm() <= 0.0 ||
                    localChanges[index] > changeThreshold ||
                    tangentAngle(reference, tangents[index]) > changeThreshold)
                {
                    return false;
                }
                sum += tangents[index];
            }
            if (sum.norm() <= 1.0e-15)
            {
                return false;
            }
            if (averageTangent)
            {
                *averageTangent = sum.normalized();
            }
            return true;
        }

        std::optional<std::size_t> findOuterBoundary(
            const ToothCluster& cluster,
            const std::vector<RegionLabel>& labels,
            bool leftSide)
        {
            if (leftSide)
            {
                for (std::size_t position = cluster.start; position <= cluster.end; ++position)
                {
                    const std::size_t index = position % labels.size();
                    if (isOuterToothBoundary(labels[index]))
                    {
                        return position;
                    }
                }
                return std::nullopt;
            }

            for (std::size_t position = cluster.end + 1; position-- > cluster.start;)
            {
                const std::size_t index = position % labels.size();
                if (isOuterToothBoundary(labels[index]))
                {
                    return position;
                }
            }
            return std::nullopt;
        }

        std::vector<std::size_t> collectOutwardUnclassified(
            const ToothCluster& cluster,
            const std::vector<RegionLabel>& labels,
            bool closed,
            bool leftSide,
            std::size_t maximumSearch)
        {
            std::vector<std::size_t> indices;
            indices.reserve(maximumSearch);
            for (std::size_t offset = 1; offset <= maximumSearch; ++offset)
            {
                std::size_t index = 0;
                if (leftSide)
                {
                    if (!closed && cluster.start < offset)
                    {
                        break;
                    }
                    const long long position = static_cast<long long>(cluster.start) -
                        static_cast<long long>(offset);
                    index = closed
                        ? wrappedIndex(position, labels.size())
                        : static_cast<std::size_t>(position);
                }
                else
                {
                    const std::size_t position = cluster.end + offset;
                    if (!closed && position >= labels.size())
                    {
                        break;
                    }
                    index = closed ? position % labels.size() : position;
                }

                if (labels[index] != RegionLabel::Unclassified)
                {
                    break;
                }
                indices.push_back(index);
            }
            return indices;
        }

        std::size_t markOneTransitionSide(
            std::vector<RegionLabel>& labels,
            std::vector<double>& confidence,
            const ToothCluster& cluster,
            const SampledContour& sampled,
            const std::vector<Eigen::Vector2d>& tangents,
            const std::vector<double>& localChanges,
            bool closed,
            bool leftSide,
            std::size_t transitionExtent,
            std::size_t stableBodyBuffer,
            double axialEdgeTolerance,
            double minimumConfidence)
        {
            const std::optional<std::size_t> boundary =
                findOuterBoundary(cluster, labels, leftSide);
            if (!boundary)
            {
                return 0;
            }
            const std::size_t boundaryIndex = *boundary % labels.size();
            const double boundaryConfidence = confidence[boundaryIndex];
            if (boundaryConfidence < minimumConfidence)
            {
                return 0;
            }

            const std::size_t maximumSearch = std::min(
                labels.size(),
                std::max<std::size_t>(12, 3 * transitionExtent));
            const std::vector<std::size_t> outward = collectOutwardUnclassified(
                cluster,
                labels,
                closed,
                leftSide,
                maximumSearch);
            constexpr std::size_t stableSampleCount = 3;
            constexpr double tangentChangeThreshold = 0.16;
            if (outward.size() <= stableSampleCount)
            {
                return 0;
            }

            std::optional<std::size_t> farBodyStart;
            Eigen::Vector2d bodyTangent = Eigen::Vector2d::Zero();
            for (std::size_t start = outward.size() - stableSampleCount + 1; start-- > 0;)
            {
                Eigen::Vector2d candidateTangent;
                if (stableTangentRun(
                        outward,
                        start,
                        stableSampleCount,
                        tangents,
                        localChanges,
                        tangentChangeThreshold,
                        &candidateTangent))
                {
                    farBodyStart = start;
                    bodyTangent = candidateTangent;
                    break;
                }
            }
            if (!farBodyStart)
            {
                return 0;
            }

            std::optional<std::size_t> bodyStart;
            for (std::size_t start = 0; start <= *farBodyStart; ++start)
            {
                Eigen::Vector2d candidateTangent;
                if (!stableTangentRun(
                        outward,
                        start,
                        stableSampleCount,
                        tangents,
                        localChanges,
                        tangentChangeThreshold,
                        &candidateTangent))
                {
                    continue;
                }
                if (tangentAngle(candidateTangent, bodyTangent) <= tangentChangeThreshold)
                {
                    bodyStart = start;
                    break;
                }
            }
            if (!bodyStart || *bodyStart == 0)
            {
                return 0;
            }

            double maximumChange = 0.0;
            for (std::size_t offset = 0; offset < *bodyStart; ++offset)
            {
                const std::size_t index = outward[offset];
                maximumChange = std::max(
                    maximumChange,
                    std::max(
                        tangentAngle(tangents[index], bodyTangent),
                        localChanges[index]));
            }
            if (maximumChange < tangentChangeThreshold)
            {
                return 0;
            }

            const double geometricConfidence = std::clamp(
                maximumChange / (3.0 * tangentChangeThreshold),
                0.0,
                1.0);
            const double rangeConfidence = std::min(
                boundaryConfidence,
                0.42 + 0.48 * geometricConfidence);
            std::size_t marked = 0;
            const std::size_t transitionCandidateCount =
                *bodyStart > stableBodyBuffer
                    ? *bodyStart - stableBodyBuffer
                    : *bodyStart;
            const std::size_t maximumMarked = std::min(
                transitionCandidateCount,
                transitionExtent);
            Eigen::AlignedBox2d bounds;
            for (const Eigen::Vector2d& point : sampled.points)
            {
                bounds.extend(point);
            }
            for (std::size_t offset = 0; offset < maximumMarked; ++offset)
            {
                const double decay = 1.0 - 0.18 * static_cast<double>(offset) /
                    static_cast<double>(std::max<std::size_t>(1, maximumMarked));
                const double candidateConfidence = rangeConfidence * decay;
                const std::size_t index = outward[offset];
                const std::size_t next = (index + 1) % sampled.smoothed.size();
                const double startZ = sampled.points[index].y();
                const double endZ = sampled.points[next].y();
                const bool onAxialEndFace =
                    ((std::abs(startZ - bounds.min().y()) <= axialEdgeTolerance &&
                         std::abs(endZ - bounds.min().y()) <= axialEdgeTolerance) ||
                        (std::abs(startZ - bounds.max().y()) <= axialEdgeTolerance &&
                            std::abs(endZ - bounds.max().y()) <= axialEdgeTolerance)) &&
                    std::abs(tangents[index].x()) >= 0.65;
                if (onAxialEndFace)
                {
                    break;
                }
                if (candidateConfidence >= minimumConfidence &&
                    labels[index] == RegionLabel::Unclassified)
                {
                    labels[index] = RegionLabel::Transition;
                    confidence[index] = candidateConfidence;
                    ++marked;
                }
            }
            return marked;
        }

        TransitionDetectionSummary markTransitionRanges(
            std::vector<RegionLabel>& labels,
            std::vector<double>& confidence,
            const SampledContour& sampled,
            std::size_t transitionExtent,
            double minimumProminence,
            double minimumConfidence,
            double axialSurfaceThreshold)
        {
            TransitionDetectionSummary summary;
            if (labels.empty())
            {
                return summary;
            }

            const std::size_t maximumInternalGap = std::max<std::size_t>(
                2,
                transitionExtent / 2);
            const std::vector<ToothCluster> clusters = buildToothClusters(
                labels,
                sampled.closed,
                maximumInternalGap);

            struct TopRun
            {
                std::size_t start{ 0 };
                std::size_t end{ 0 };
                double meanY{ 0.0 };
                double meanZ{ 0.0 };
            };
            struct CandidateGroup
            {
                std::size_t coreStart{ 0 };
                std::size_t coreEnd{ 0 };
                std::vector<TopRun> topRuns;
                bool repeatedPattern{ false };
            };

            Eigen::AlignedBox2d bounds;
            for (const Eigen::Vector2d& point : sampled.smoothed)
            {
                bounds.extend(point);
            }
            const double radialRange = bounds.max().x() - bounds.min().x();
            constexpr std::size_t minimumTopRunSamples = 2;

            const std::size_t microWindow = std::max<std::size_t>(
                2,
                static_cast<std::size_t>(std::ceil(
                    0.0025 * static_cast<double>(labels.size()))));
            const auto isTopSeedSample = [&](std::size_t position)
            {
                if (labels[position] == RegionLabel::ToothTop)
                {
                    return true;
                }
                const std::size_t next = (position + 1) % sampled.smoothed.size();
                const Eigen::Vector2d direction =
                    sampled.smoothed[next] - sampled.smoothed[position];
                const double length = direction.norm();
                if (length <= 1.0e-15 ||
                    std::abs(direction.y()) / length < axialSurfaceThreshold)
                {
                    return false;
                }

                double localMinimumY = std::numeric_limits<double>::infinity();
                double localMaximumY = -std::numeric_limits<double>::infinity();
                for (long long offset = -static_cast<long long>(microWindow);
                     offset <= static_cast<long long>(microWindow);
                     ++offset)
                {
                    const long long candidate =
                        static_cast<long long>(position) + offset;
                    if (!sampled.closed &&
                        (candidate < 0 ||
                            candidate >= static_cast<long long>(sampled.smoothed.size())))
                    {
                        continue;
                    }
                    const std::size_t index = sampled.closed
                        ? wrappedIndex(candidate, sampled.smoothed.size())
                        : static_cast<std::size_t>(candidate);
                    localMinimumY = std::min(
                        localMinimumY,
                        sampled.smoothed[index].x());
                    localMaximumY = std::max(
                        localMaximumY,
                        sampled.smoothed[index].x());
                }
                const double localSpan = localMaximumY - localMinimumY;
                if (localSpan < 0.45 * minimumProminence)
                {
                    return false;
                }
                const double middleY = 0.5 *
                    (sampled.smoothed[position].x() + sampled.smoothed[next].x());
                return middleY >= localMinimumY + 0.62 * localSpan;
            };

            const auto makeTopRun = [&](std::size_t start, std::size_t end)
            {
                double sumY = 0.0;
                double sumZ = 0.0;
                std::size_t count = 0;
                for (std::size_t position = start; position <= end; ++position)
                {
                    const std::size_t index = position % labels.size();
                    const std::size_t next = (index + 1) % sampled.smoothed.size();
                    sumY += 0.5 *
                        (sampled.smoothed[index].x() + sampled.smoothed[next].x());
                    sumZ += 0.5 *
                        (sampled.smoothed[index].y() + sampled.smoothed[next].y());
                    ++count;
                }
                return TopRun{
                    start,
                    end,
                    sumY / static_cast<double>(count),
                    sumZ / static_cast<double>(count) };
            };

            std::vector<TopRun> allTopRuns;
            for (std::size_t position = 0; position < labels.size();)
            {
                if (!isTopSeedSample(position))
                {
                    ++position;
                    continue;
                }
                const std::size_t start = position;
                std::size_t end = position;
                while (end + 1 < labels.size() &&
                    isTopSeedSample(end + 1))
                {
                    ++end;
                }
                if (end - start + 1 >= minimumTopRunSamples)
                {
                    allTopRuns.push_back(makeTopRun(start, end));
                }
                position = end + 1;
            }
            if (sampled.closed && allTopRuns.size() > 1 &&
                isTopSeedSample(0) &&
                isTopSeedSample(labels.size() - 1))
            {
                const std::size_t wrapAnchor = allTopRuns.back().start;
                const TopRun merged = makeTopRun(
                    wrapAnchor,
                    allTopRuns.front().end + labels.size());
                allTopRuns.erase(allTopRuns.begin());
                allTopRuns.pop_back();
                for (TopRun& run : allTopRuns)
                {
                    if (run.start < wrapAnchor)
                    {
                        run.start += labels.size();
                        run.end += labels.size();
                    }
                }
                allTopRuns.push_back(merged);
                std::sort(allTopRuns.begin(), allTopRuns.end(),
                    [](const TopRun& first, const TopRun& second)
                    {
                        return first.start < second.start;
                    });
            }

            const auto median = [](std::vector<double> values)
            {
                std::sort(values.begin(), values.end());
                const std::size_t middle = values.size() / 2;
                return values.size() % 2 == 0
                    ? 0.5 * (values[middle - 1] + values[middle])
                    : values[middle];
            };

            // A closed section also contains the flange, inner wall, and end faces.
            // Direction alone makes them look like tooth tops and walls. A real tooth
            // group instead has at least three narrow top runs at a common radius and
            // a repeating axial pitch. Select that pattern before retaining labels.
            std::vector<TopRun> radialOrder = allTopRuns;
            std::sort(radialOrder.begin(), radialOrder.end(),
                [](const TopRun& first, const TopRun& second)
                {
                    return first.meanY > second.meanY;
                });
            const double radialLayerTolerance = std::max(
                0.15 * radialRange,
                16.0 * std::numeric_limits<double>::epsilon() *
                    std::max(1.0, bounds.max().cwiseAbs().maxCoeff()));
            std::vector<std::vector<TopRun>> radialLayers;
            for (const TopRun& run : radialOrder)
            {
                if (radialLayers.empty() ||
                    std::abs(radialLayers.back().front().meanY - run.meanY) >
                        radialLayerTolerance)
                {
                    radialLayers.push_back({});
                }
                radialLayers.back().push_back(run);
            }

            std::vector<TopRun> repeatedGroup;
            double bestScore = -std::numeric_limits<double>::infinity();
            for (std::vector<TopRun> layer : radialLayers)
            {
                if (layer.size() < 3)
                {
                    continue;
                }
                std::sort(layer.begin(), layer.end(),
                    [](const TopRun& first, const TopRun& second)
                    {
                        return first.meanZ < second.meanZ;
                    });
                for (std::size_t first = 0; first + 2 < layer.size(); ++first)
                {
                    for (std::size_t last = first + 2; last < layer.size(); ++last)
                    {
                        std::vector<TopRun> candidate(
                            layer.begin() + static_cast<std::ptrdiff_t>(first),
                            layer.begin() + static_cast<std::ptrdiff_t>(last + 1));
                        std::vector<double> axialPitches;
                        for (std::size_t index = 1; index < candidate.size(); ++index)
                        {
                            const double pitch =
                                candidate[index].meanZ - candidate[index - 1].meanZ;
                            if (pitch > 0.0)
                            {
                                axialPitches.push_back(pitch);
                            }
                        }
                        if (axialPitches.size() + 1 != candidate.size())
                        {
                            continue;
                        }
                        const double medianPitch = median(axialPitches);
                        if (!std::isfinite(medianPitch) || medianPitch <= 0.0)
                        {
                            continue;
                        }
                        double maximumPitchDeviation = 0.0;
                        for (double pitch : axialPitches)
                        {
                            maximumPitchDeviation = std::max(
                                maximumPitchDeviation,
                                std::abs(pitch - medianPitch) / medianPitch);
                        }
                        if (maximumPitchDeviation > 0.20)
                        {
                            continue;
                        }

                        std::vector<double> runWidths;
                        for (const TopRun& run : candidate)
                        {
                            runWidths.push_back(static_cast<double>(
                                run.end - run.start + 1));
                        }
                        std::vector<TopRun> contourOrder = candidate;
                        std::sort(contourOrder.begin(), contourOrder.end(),
                            [](const TopRun& firstRun, const TopRun& secondRun)
                            {
                                return firstRun.start < secondRun.start;
                            });
                        std::vector<double> samplePitches;
                        for (std::size_t index = 1;
                             index < contourOrder.size();
                             ++index)
                        {
                            const double firstCenter = 0.5 * static_cast<double>(
                                contourOrder[index - 1].start +
                                contourOrder[index - 1].end);
                            const double secondCenter = 0.5 * static_cast<double>(
                                contourOrder[index].start + contourOrder[index].end);
                            if (secondCenter > firstCenter)
                            {
                                samplePitches.push_back(secondCenter - firstCenter);
                            }
                        }
                        if (samplePitches.size() + 1 != contourOrder.size())
                        {
                            continue;
                        }
                        const double medianSamplePitch = median(samplePitches);
                        double maximumSamplePitchDeviation = 0.0;
                        for (double samplePitch : samplePitches)
                        {
                            maximumSamplePitchDeviation = std::max(
                                maximumSamplePitchDeviation,
                                std::abs(samplePitch - medianSamplePitch) /
                                    medianSamplePitch);
                        }
                        if (maximumSamplePitchDeviation > 0.45)
                        {
                            continue;
                        }
                        const double widthToPitch = median(runWidths) /
                            std::max(medianSamplePitch, 1.0);
                        if (widthToPitch > 0.60)
                        {
                            continue;
                        }

                        std::vector<double> toothDepths;
                        for (std::size_t index = 1;
                             index < contourOrder.size();
                             ++index)
                        {
                            const TopRun& previousRun = contourOrder[index - 1];
                            const TopRun& nextRun = contourOrder[index];
                            if (previousRun.end + 1 >= nextRun.start)
                            {
                                continue;
                            }
                            double valleyY = std::numeric_limits<double>::infinity();
                            for (std::size_t position = previousRun.end + 1;
                                 position < nextRun.start;
                                 ++position)
                            {
                                const std::size_t sample = position % labels.size();
                                const std::size_t nextSample =
                                    (sample + 1) % sampled.smoothed.size();
                                valleyY = std::min(
                                    valleyY,
                                    0.5 * (sampled.smoothed[sample].x() +
                                        sampled.smoothed[nextSample].x()));
                            }
                            const double depth =
                                std::min(previousRun.meanY, nextRun.meanY) - valleyY;
                            if (std::isfinite(depth) && depth > 0.5 * minimumProminence)
                            {
                                toothDepths.push_back(depth);
                            }
                        }
                        if (toothDepths.size() + 1 != contourOrder.size())
                        {
                            continue;
                        }
                        const double medianToothDepth = median(toothDepths);
                        double maximumDepthDeviation = 0.0;
                        for (double depth : toothDepths)
                        {
                            maximumDepthDeviation = std::max(
                                maximumDepthDeviation,
                                std::abs(depth - medianToothDepth) /
                                    medianToothDepth);
                        }
                        if (maximumDepthDeviation > 0.65)
                        {
                            continue;
                        }

                        const double score =
                            100.0 * static_cast<double>(candidate.size())
                            - 20.0 * maximumPitchDeviation
                            - 20.0 * maximumSamplePitchDeviation
                            - 10.0 * maximumDepthDeviation
                            - widthToPitch;
                        if (score > bestScore)
                        {
                            bestScore = score;
                            repeatedGroup = std::move(contourOrder);
                        }
                    }
                }
            }

            std::vector<CandidateGroup> candidateGroups;
            if (!repeatedGroup.empty())
            {
                summary.outerTopRunCount = repeatedGroup.size();
                if (sampled.closed)
                {
                    for (TopRun& run : repeatedGroup)
                    {
                        const std::size_t runLength = run.end - run.start + 1;
                        run.start %= labels.size();
                        run.end = run.start + runLength - 1;
                    }
                    std::sort(repeatedGroup.begin(), repeatedGroup.end(),
                        [](const TopRun& first, const TopRun& second)
                        {
                            return first.start < second.start;
                        });

                    std::size_t largestGapAfter = 0;
                    std::size_t largestGap = 0;
                    for (std::size_t index = 0; index < repeatedGroup.size(); ++index)
                    {
                        const TopRun& current = repeatedGroup[index];
                        const TopRun& next = repeatedGroup[
                            (index + 1) % repeatedGroup.size()];
                        std::size_t nextStart = next.start;
                        if (index + 1 == repeatedGroup.size())
                        {
                            nextStart += labels.size();
                        }
                        const std::size_t gap = nextStart > current.end
                            ? nextStart - current.end - 1
                            : 0;
                        if (gap > largestGap)
                        {
                            largestGap = gap;
                            largestGapAfter = index;
                        }
                    }

                    const std::size_t firstRun =
                        (largestGapAfter + 1) % repeatedGroup.size();
                    std::vector<TopRun> unfoldedRuns;
                    unfoldedRuns.reserve(repeatedGroup.size());
                    for (std::size_t offset = 0;
                         offset < repeatedGroup.size();
                         ++offset)
                    {
                        TopRun run = repeatedGroup[
                            (firstRun + offset) % repeatedGroup.size()];
                        if (!unfoldedRuns.empty())
                        {
                            while (run.start <= unfoldedRuns.back().end)
                            {
                                run.start += labels.size();
                                run.end += labels.size();
                            }
                        }
                        unfoldedRuns.push_back(run);
                    }
                    repeatedGroup = std::move(unfoldedRuns);
                }

                std::size_t coreStart = repeatedGroup.front().start;
                std::size_t coreEnd = repeatedGroup.back().end;
                if (sampled.closed)
                {
                    coreStart += labels.size();
                    coreEnd += labels.size();
                    for (TopRun& run : repeatedGroup)
                    {
                        run.start += labels.size();
                        run.end += labels.size();
                    }
                }

                candidateGroups.push_back({
                    coreStart,
                    coreEnd,
                    std::move(repeatedGroup),
                    true });
            }
            else
            {
                // A one-tooth simulation block has no repeat pattern. Keep a wider
                // outside-radius fallback, but still scope it to the preliminary
                // geometric clusters so the rest of the contour remains editable as
                // Unclassified rather than being silently accepted as spray area.
                constexpr double fallbackOuterBandFraction = 0.30;
                const double fallbackMinimumY =
                    bounds.max().x() - fallbackOuterBandFraction * radialRange;
                for (const ToothCluster& cluster : clusters)
                {
                    std::vector<TopRun> topRuns;
                    for (const TopRun& run : allTopRuns)
                    {
                        if (run.meanY < fallbackMinimumY)
                        {
                            continue;
                        }
                        const std::size_t runLength = run.end - run.start + 1;
                        for (long long shift = sampled.closed ? -2 : 0;
                             shift <= (sampled.closed ? 2 : 0);
                             ++shift)
                        {
                            const long long shiftedStart =
                                static_cast<long long>(run.start) +
                                shift * static_cast<long long>(labels.size());
                            const long long shiftedEnd = shiftedStart +
                                static_cast<long long>(runLength) - 1;
                            if (shiftedStart >= 0 &&
                                shiftedStart >= static_cast<long long>(cluster.start) &&
                                shiftedEnd <= static_cast<long long>(cluster.end))
                            {
                                TopRun aligned = run;
                                aligned.start = static_cast<std::size_t>(shiftedStart);
                                aligned.end = static_cast<std::size_t>(shiftedEnd);
                                topRuns.push_back(aligned);
                                break;
                            }
                        }
                    }
                    if (!topRuns.empty())
                    {
                        std::sort(topRuns.begin(), topRuns.end(),
                            [](const TopRun& first, const TopRun& second)
                            {
                                return first.start < second.start;
                            });
                        std::vector<double> fallbackTopWidths;
                        fallbackTopWidths.reserve(topRuns.size());
                        for (const TopRun& run : topRuns)
                        {
                            fallbackTopWidths.push_back(static_cast<double>(
                                run.end - run.start + 1));
                        }
                        const std::size_t boundaryPadding = std::max<std::size_t>(
                            4,
                            static_cast<std::size_t>(std::ceil(
                                3.0 * median(fallbackTopWidths))));
                        const std::size_t coreStart =
                            topRuns.front().start > boundaryPadding
                                ? std::max(
                                    cluster.start,
                                    topRuns.front().start - boundaryPadding)
                                : cluster.start;
                        const std::size_t coreEnd = std::min(
                            cluster.end,
                            topRuns.back().end + boundaryPadding);
                        summary.outerTopRunCount += topRuns.size();
                        candidateGroups.push_back({
                            coreStart,
                            coreEnd,
                            std::move(topRuns),
                            false });
                    }
                }
            }

            summary.clusterCount = candidateGroups.size();
            std::vector<bool> coreMask(labels.size(), false);
            const std::vector<RegionLabel> preliminaryLabels = labels;
            const std::vector<double> preliminaryConfidence = confidence;
            const std::vector<Eigen::Vector2d> tangents = segmentTangents(sampled);
            const auto segmentMiddleY = [&](std::size_t position)
            {
                const std::size_t index = position % labels.size();
                const std::size_t next = (index + 1) % sampled.smoothed.size();
                return 0.5 *
                    (sampled.smoothed[index].x() + sampled.smoothed[next].x());
            };
            const auto assignSurface = [&](
                std::size_t position,
                RegionLabel label,
                double normalStrength)
            {
                const std::size_t index = position % labels.size();
                labels[index] = label;
                confidence[index] = std::clamp(
                    std::max(minimumConfidence, 0.55 + 0.40 * normalStrength),
                    0.0,
                    1.0);
            };

            for (CandidateGroup& group : candidateGroups)
            {
                std::sort(group.topRuns.begin(), group.topRuns.end(),
                    [](const TopRun& first, const TopRun& second)
                    {
                        return first.start < second.start;
                    });

                std::vector<double> samplePitches;
                for (std::size_t runIndex = 1;
                     runIndex < group.topRuns.size();
                     ++runIndex)
                {
                    const TopRun& previousRun = group.topRuns[runIndex - 1];
                    const TopRun& nextRun = group.topRuns[runIndex];
                    const double previousCenter = 0.5 * static_cast<double>(
                        previousRun.start + previousRun.end);
                    const double nextCenter = 0.5 * static_cast<double>(
                        nextRun.start + nextRun.end);
                    if (nextCenter > previousCenter)
                    {
                        samplePitches.push_back(nextCenter - previousCenter);
                    }
                }
                const double medianSamplePitch = samplePitches.empty()
                    ? 2.0 * static_cast<double>(transitionExtent)
                    : median(samplePitches);

                // Use the repeated valleys as the body-side radius of the teeth.
                // Starting at each outer top, walk only until its wall first reaches
                // that radius. This protects the complete outer wall without letting
                // a fixed half-pitch spill into the transition or the stable body.
                std::vector<double> valleyLevels;
                for (std::size_t runIndex = 1;
                     runIndex < group.topRuns.size();
                     ++runIndex)
                {
                    const TopRun& previousRun = group.topRuns[runIndex - 1];
                    const TopRun& nextRun = group.topRuns[runIndex];
                    if (previousRun.end + 1 >= nextRun.start)
                    {
                        continue;
                    }
                    double valleyY = std::numeric_limits<double>::infinity();
                    for (std::size_t position = previousRun.end + 1;
                         position < nextRun.start;
                         ++position)
                    {
                        valleyY = std::min(valleyY, segmentMiddleY(position));
                    }
                    if (std::isfinite(valleyY))
                    {
                        valleyLevels.push_back(valleyY);
                    }
                }
                if (!group.topRuns.empty() && !valleyLevels.empty())
                {
                    const double bodySideY = median(valleyLevels);
                    std::vector<double> topLevels;
                    topLevels.reserve(group.topRuns.size());
                    for (const TopRun& run : group.topRuns)
                    {
                        double maximumY = -std::numeric_limits<double>::infinity();
                        for (std::size_t position = run.start;
                             position <= run.end;
                             ++position)
                        {
                            maximumY = std::max(
                                maximumY,
                                segmentMiddleY(position));
                        }
                        topLevels.push_back(maximumY);
                    }
                    const double toothDepth = std::max(
                        0.0,
                        median(topLevels) - bodySideY);
                    const double localGeometryTolerance = std::max(
                        64.0 * std::numeric_limits<double>::epsilon() *
                            std::max(1.0, bounds.max().cwiseAbs().maxCoeff()),
                        0.04 * toothDepth);
                    const double wallFootY = bodySideY + localGeometryTolerance;
                    const std::size_t maximumWallSearch = std::max<std::size_t>(
                        3,
                        static_cast<std::size_t>(std::ceil(
                            0.75 * medianSamplePitch)));

                    const auto wallBoundary = [&](bool leftSide)
                    {
                        std::size_t boundary = leftSide
                            ? group.topRuns.front().start
                            : group.topRuns.back().end;
                        bool foundRadialWall = false;
                        for (std::size_t offset = 1;
                             offset <= maximumWallSearch;
                             ++offset)
                        {
                            const long long position = leftSide
                                ? static_cast<long long>(
                                    group.topRuns.front().start) -
                                    static_cast<long long>(offset)
                                : static_cast<long long>(
                                    group.topRuns.back().end) +
                                    static_cast<long long>(offset);
                            if (!sampled.closed &&
                                (position < 0 ||
                                    position >= static_cast<long long>(labels.size())))
                            {
                                break;
                            }
                            const std::size_t index = sampled.closed
                                ? wrappedIndex(position, labels.size())
                                : static_cast<std::size_t>(position);
                            const std::size_t next =
                                (index + 1) % sampled.smoothed.size();
                            const double axialNormalStrength =
                                std::abs(tangents[index].x());
                            if (foundRadialWall && axialNormalStrength < 0.45)
                            {
                                break;
                            }
                            boundary = static_cast<std::size_t>(position);
                            foundRadialWall = foundRadialWall ||
                                axialNormalStrength >= 0.58;
                            const double endpointMinimumY = std::min(
                                sampled.smoothed[index].x(),
                                sampled.smoothed[next].x());
                            if (endpointMinimumY <= wallFootY)
                            {
                                break;
                            }
                        }
                        return boundary;
                    };

                    group.coreStart = std::min(
                        group.coreStart,
                        wallBoundary(true));
                    group.coreEnd = std::max(
                        group.coreEnd,
                        wallBoundary(false));
                }
                if (!sampled.closed)
                {
                    group.coreStart = std::min(group.coreStart, labels.size() - 1);
                    group.coreEnd = std::min(group.coreEnd, labels.size() - 1);
                }

                for (std::size_t position = group.coreStart;
                     position <= group.coreEnd;
                     ++position)
                {
                    const std::size_t index = position % labels.size();
                    labels[index] = RegionLabel::Unclassified;
                    confidence[index] = 0.0;
                    coreMask[index] = true;
                }

                // Confirm each tooth top by both its local radial maximum and its
                // radial-facing normal. Preliminary high sloping walls therefore no
                // longer remain ToothTop merely because their tangent is partly axial.
                for (const TopRun& run : group.topRuns)
                {
                    double maximumY = -std::numeric_limits<double>::infinity();
                    std::size_t maximumPosition = run.start;
                    for (std::size_t position = run.start;
                         position <= run.end;
                         ++position)
                    {
                        const double middleY = segmentMiddleY(position);
                        if (middleY > maximumY)
                        {
                            maximumY = middleY;
                            maximumPosition = position;
                        }
                    }

                    const std::size_t neighbourhood = std::max<std::size_t>(
                        2,
                        static_cast<std::size_t>(std::ceil(
                            0.42 * medianSamplePitch)));
                    double nearbyMinimumY = maximumY;
                    for (long long offset = -static_cast<long long>(neighbourhood);
                         offset <= static_cast<long long>(neighbourhood);
                         ++offset)
                    {
                        const long long candidate = static_cast<long long>(
                            maximumPosition) + offset;
                        if (!sampled.closed &&
                            (candidate < 0 ||
                                candidate >= static_cast<long long>(labels.size())))
                        {
                            continue;
                        }
                        nearbyMinimumY = std::min(
                            nearbyMinimumY,
                            segmentMiddleY(sampled.closed
                                ? wrappedIndex(candidate, labels.size())
                                : static_cast<std::size_t>(candidate)));
                    }
                    const double toothDepth = std::max(0.0, maximumY - nearbyMinimumY);
                    const double topFloor = maximumY - std::max(
                        radialLayerTolerance,
                        0.12 * toothDepth);

                    bool assignedTop = false;
                    for (std::size_t position = run.start;
                         position <= run.end;
                         ++position)
                    {
                        const std::size_t index = position % labels.size();
                        const double radialNormalStrength =
                            std::abs(tangents[index].y());
                        if (segmentMiddleY(position) >= topFloor &&
                            radialNormalStrength >= 0.50)
                        {
                            assignSurface(
                                position,
                                RegionLabel::ToothTop,
                                radialNormalStrength);
                            assignedTop = true;
                        }
                    }
                    if (!assignedTop)
                    {
                        const std::size_t index = maximumPosition % labels.size();
                        assignSurface(
                            maximumPosition,
                            RegionLabel::ToothTop,
                            std::abs(tangents[index].y()));
                    }
                }

                // Once two neighbouring tooth tops establish a period, position is
                // more reliable than a single noisy normal. The whole interval is a
                // continuous wall-bottom-wall sequence; only the low radial-facing
                // valley may be ToothBottom and every remaining piece is ToothWall.
                for (std::size_t runIndex = 1;
                     runIndex < group.topRuns.size();
                     ++runIndex)
                {
                    const TopRun& previousRun = group.topRuns[runIndex - 1];
                    const TopRun& nextRun = group.topRuns[runIndex];
                    if (previousRun.end + 1 >= nextRun.start)
                    {
                        continue;
                    }

                    double topReferenceY = std::numeric_limits<double>::infinity();
                    for (const TopRun* run : { &previousRun, &nextRun })
                    {
                        double runMaximumY =
                            -std::numeric_limits<double>::infinity();
                        for (std::size_t position = run->start;
                             position <= run->end;
                             ++position)
                        {
                            runMaximumY = std::max(
                                runMaximumY,
                                segmentMiddleY(position));
                        }
                        topReferenceY = std::min(topReferenceY, runMaximumY);
                    }
                    double valleyY = std::numeric_limits<double>::infinity();
                    for (std::size_t position = previousRun.end + 1;
                         position < nextRun.start;
                         ++position)
                    {
                        valleyY = std::min(valleyY, segmentMiddleY(position));
                    }
                    const double valleyDepth = topReferenceY - valleyY;
                    const double bottomMaximumY = valleyY +
                        0.18 * std::max(0.0, valleyDepth);
                    for (std::size_t position = previousRun.end + 1;
                         position < nextRun.start;
                         ++position)
                    {
                        const std::size_t index = position % labels.size();
                        const double radialNormalStrength =
                            std::abs(tangents[index].y());
                        const double axialNormalStrength =
                            std::abs(tangents[index].x());
                        const bool isValleySurface =
                            std::isfinite(valleyDepth) &&
                            valleyDepth > 0.5 * minimumProminence &&
                            segmentMiddleY(position) <= bottomMaximumY &&
                            radialNormalStrength >= axialSurfaceThreshold;
                        assignSurface(
                            position,
                            isValleySurface
                                ? RegionLabel::ToothBottom
                                : RegionLabel::ToothWall,
                            isValleySurface
                                ? radialNormalStrength
                                : axialNormalStrength);
                    }
                }

                // The two outer half-periods are the exposed walls of the boundary
                // teeth. Fill every non-top sample so short or slanted wall pieces do
                // not disappear and cannot later be relabelled as Transition.
                for (std::size_t position = group.coreStart;
                     position <= group.coreEnd;
                     ++position)
                {
                    const std::size_t index = position % labels.size();
                    if (labels[index] == RegionLabel::Unclassified)
                    {
                        if (group.repeatedPattern)
                        {
                            assignSurface(
                                position,
                                RegionLabel::ToothWall,
                                std::abs(tangents[index].x()));
                        }
                        else if (isToothLabel(preliminaryLabels[index]))
                        {
                            labels[index] = preliminaryLabels[index];
                            confidence[index] = preliminaryConfidence[index];
                        }
                    }
                }
            }

            for (std::size_t index = 0; index < labels.size(); ++index)
            {
                if (!coreMask[index])
                {
                    labels[index] = RegionLabel::Unclassified;
                    confidence[index] = 0.0;
                }
            }

            const std::vector<double> localChanges =
                localTangentChanges(tangents, sampled.closed);
            const double axialRange = bounds.max().y() - bounds.min().y();
            const double axialEdgeTolerance = std::max(
                1.0e-6 * axialRange,
                64.0 * std::numeric_limits<double>::epsilon() *
                    std::max(1.0, bounds.max().cwiseAbs().maxCoeff()));
            for (const CandidateGroup& group : candidateGroups)
            {
                const ToothCluster selectedCluster{
                    group.coreStart,
                    group.coreEnd,
                    group.coreEnd - group.coreStart + 1 >= labels.size() };
                const std::size_t leftMarked = markOneTransitionSide(
                    labels,
                    confidence,
                    selectedCluster,
                    sampled,
                    tangents,
                    localChanges,
                    sampled.closed,
                    true,
                    transitionExtent,
                    group.topRuns.size() >= 3 ? 0 : std::size_t{ 3 },
                    axialEdgeTolerance,
                    minimumConfidence);
                const std::size_t rightMarked = markOneTransitionSide(
                    labels,
                    confidence,
                    selectedCluster,
                    sampled,
                    tangents,
                    localChanges,
                    sampled.closed,
                    false,
                    transitionExtent,
                    group.topRuns.size() >= 3 ? 0 : std::size_t{ 3 },
                    axialEdgeTolerance,
                    minimumConfidence);
                summary.markedSegmentCount += leftMarked + rightMarked;
                summary.qualifiedBoundaryCount += leftMarked > 0 ? 1 : 0;
                summary.qualifiedBoundaryCount += rightMarked > 0 ? 1 : 0;
            }
            return summary;
        }

        const char* labelName(RegionLabel label) noexcept
        {
            switch (label)
            {
            case RegionLabel::ToothTop:
                return "ToothTop";
            case RegionLabel::ToothWall:
                return "ToothWall";
            case RegionLabel::ToothBottom:
                return "ToothBottom";
            case RegionLabel::Transition:
                return "Transition";
            case RegionLabel::Unclassified:
                return "Unclassified";
            }
            return "Unclassified";
        }
    }

    PlanningResult<RegionAssignment> ToothRegionRecognizer::recognize(
        const SectionContour& contour,
        const ToothRecognitionOptions& options)
    {
        SectionContour workingContour = contour;
        const std::size_t contourRotation = canonicalizeClosedContour(workingContour);
        auto sampledResult = resample(workingContour, options);
        if (!sampledResult)
        {
            return PlanningResult<RegionAssignment>::failure(
                sampledResult.error.code,
                sampledResult.error.message);
        }
        SampledContour sampled = std::move(sampledResult.value);
        const std::size_t sampledSegments = sampled.segmentCount();

        Eigen::AlignedBox2d bounds;
        for (const Eigen::Vector2d& point : sampled.smoothed)
        {
            bounds.extend(point);
        }
        const double radialRange = bounds.max().x() - bounds.min().x();
        const double axialRange = bounds.max().y() - bounds.min().y();
        const double geometryScale = std::max(radialRange, axialRange);
        if (!std::isfinite(geometryScale) || geometryScale <= 0.0)
        {
            return PlanningResult<RegionAssignment>::failure(
                PlanningErrorCode::DegenerateGeometry,
                "Section contour has no measurable YZ extent.");
        }
        const double prominenceScale = sampledSegments >= 512
            ? std::max(
                0.35 * std::min(radialRange, axialRange),
                0.005 * geometryScale)
            : geometryScale;
        const double minimumProminence = std::max(
            prominenceScale * options.minimumRadialProminenceFraction,
            1.0e-12);
        const std::size_t localWindow = std::max<std::size_t>(
            4,
            static_cast<std::size_t>(std::ceil(
                (sampledSegments >= 512 ? 0.012 : 0.055) *
                    static_cast<double>(sampledSegments))));

        std::vector<RegionLabel> sampledLabels(sampledSegments, RegionLabel::Unclassified);
        std::vector<double> sampledConfidence(sampledSegments, 0.0);
        for (std::size_t segment = 0; segment < sampledSegments; ++segment)
        {
            const std::size_t endIndex = (segment + 1) % sampled.smoothed.size();
            const Eigen::Vector2d direction = sampled.smoothed[endIndex] - sampled.smoothed[segment];
            const double length = direction.norm();
            if (length <= 1.0e-15)
            {
                continue;
            }
            const double axiality = std::abs(direction.y()) / length;
            const double radiality = std::abs(direction.x()) / length;
            const double middleY = 0.5 *
                (sampled.smoothed[segment].x() + sampled.smoothed[endIndex].x());

            double localMinimumY = std::numeric_limits<double>::infinity();
            double localMaximumY = -std::numeric_limits<double>::infinity();
            for (long long offset = -static_cast<long long>(localWindow);
                 offset <= static_cast<long long>(localWindow);
                 ++offset)
            {
                long long candidate = static_cast<long long>(segment) + offset;
                if (!sampled.closed &&
                    (candidate < 0 || candidate >= static_cast<long long>(sampled.smoothed.size())))
                {
                    continue;
                }
                const std::size_t index = sampled.closed
                    ? wrappedIndex(candidate, sampled.smoothed.size())
                    : static_cast<std::size_t>(candidate);
                localMinimumY = std::min(localMinimumY, sampled.smoothed[index].x());
                localMaximumY = std::max(localMaximumY, sampled.smoothed[index].x());
            }
            const double localSpan = localMaximumY - localMinimumY;
            if (localSpan < minimumProminence)
            {
                continue;
            }
            const double normalizedHeight = std::clamp(
                (middleY - localMinimumY) / std::max(localSpan, 1.0e-15),
                0.0,
                1.0);
            const double prominenceConfidence = std::clamp(
                localSpan / std::max(4.0 * minimumProminence, 1.0e-15),
                0.0,
                1.0);

            RegionLabel label = RegionLabel::Unclassified;
            double confidence = 0.0;
            if (axiality >= options.axialSurfaceThreshold && normalizedHeight >= 0.57)
            {
                label = RegionLabel::ToothTop;
                confidence = axiality * prominenceConfidence *
                    std::clamp((normalizedHeight - 0.48) / 0.52, 0.0, 1.0);
            }
            else if (axiality >= options.axialSurfaceThreshold && normalizedHeight <= 0.43)
            {
                label = RegionLabel::ToothBottom;
                confidence = axiality * prominenceConfidence *
                    std::clamp((0.52 - normalizedHeight) / 0.52, 0.0, 1.0);
            }
            else if (radiality >= 0.52)
            {
                label = RegionLabel::ToothWall;
                confidence = radiality * prominenceConfidence;
            }
            if (confidence >= options.minimumConfidence)
            {
                sampledLabels[segment] = label;
                sampledConfidence[segment] = std::clamp(confidence, 0.0, 1.0);
            }
        }

        const std::size_t sampledRotation = stableClosedClassificationSeam(
            sampledLabels,
            sampled.closed);
        if (sampledRotation != 0)
        {
            rotateLeft(sampled.points, sampledRotation);
            rotateLeft(sampled.smoothed, sampledRotation);
            rotateLeft(sampledLabels, sampledRotation);
            rotateLeft(sampledConfidence, sampledRotation);
        }

        const std::size_t transitionExtent = std::max<std::size_t>(
            2,
            static_cast<std::size_t>(
                std::ceil(options.transitionExtentFraction * static_cast<double>(sampledSegments))));
        const TransitionDetectionSummary transitionSummary = markTransitionRanges(
            sampledLabels,
            sampledConfidence,
            sampled,
            transitionExtent,
            minimumProminence,
            options.minimumConfidence,
            options.axialSurfaceThreshold);

        if (sampledRotation != 0)
        {
            undoRotateLeft(sampledLabels, sampledRotation);
            undoRotateLeft(sampledConfidence, sampledRotation);
        }

        std::vector<double> originalCumulative(workingContour.segmentCount() + 1, 0.0);
        for (std::size_t segment = 0; segment < workingContour.segmentCount(); ++segment)
        {
            originalCumulative[segment + 1] = originalCumulative[segment] +
                (workingContour.pointsYz[(segment + 1) % workingContour.pointsYz.size()] -
                    workingContour.pointsYz[segment]).norm();
        }

        RegionAssignment result;
        result.segmentLabels.resize(workingContour.segmentCount(), RegionLabel::Unclassified);
        result.segmentConfidence.resize(workingContour.segmentCount(), 0.0);
        const double sampleStep = sampled.totalLength /
            static_cast<double>(sampledSegments);
        constexpr std::array<RegionLabel, 5> labelsByIndex{
            RegionLabel::ToothTop,
            RegionLabel::ToothWall,
            RegionLabel::ToothBottom,
            RegionLabel::Transition,
            RegionLabel::Unclassified };
        const auto labelIndex = [](RegionLabel label)
        {
            switch (label)
            {
            case RegionLabel::ToothTop:
                return std::size_t{ 0 };
            case RegionLabel::ToothWall:
                return std::size_t{ 1 };
            case RegionLabel::ToothBottom:
                return std::size_t{ 2 };
            case RegionLabel::Transition:
                return std::size_t{ 3 };
            case RegionLabel::Unclassified:
                return std::size_t{ 4 };
            }
            return std::size_t{ 4 };
        };
        for (std::size_t segment = 0; segment < workingContour.segmentCount(); ++segment)
        {
            const double segmentStart = originalCumulative[segment];
            const double segmentEnd = originalCumulative[segment + 1];
            const std::size_t firstSample = std::min(
                sampledSegments - 1,
                static_cast<std::size_t>(std::floor(
                    segmentStart / sampleStep)));
            const std::size_t sampleEndExclusive = std::max<std::size_t>(
                1,
                static_cast<std::size_t>(std::ceil(
                    segmentEnd / sampleStep)));
            const std::size_t lastSample = std::min(
                sampledSegments - 1,
                sampleEndExclusive - 1);
            std::array<double, 5> coverage{};
            std::array<double, 5> maximumConfidence{};
            for (std::size_t sample = firstSample;
                 sample <= lastSample;
                 ++sample)
            {
                const double sampleStart = sampleStep * static_cast<double>(sample);
                const double sampleEnd = sampleStep * static_cast<double>(sample + 1);
                const double overlap = std::max(
                    0.0,
                    std::min(segmentEnd, sampleEnd) -
                        std::max(segmentStart, sampleStart));
                if (overlap <= 0.0)
                {
                    continue;
                }
                const std::size_t index = labelIndex(sampledLabels[sample]);
                coverage[index] += overlap;
                maximumConfidence[index] = std::max(
                    maximumConfidence[index],
                    sampledConfidence[sample]);
            }

            std::size_t selected = 4;
            const double segmentLength = segmentEnd - segmentStart;
            const double coverageTolerance = 32.0 *
                std::numeric_limits<double>::epsilon() *
                std::max(1.0, sampled.totalLength);
            for (std::size_t candidate = 0; candidate < 4; ++candidate)
            {
                if (coverage[candidate] >
                    0.5 * segmentLength + coverageTolerance)
                {
                    selected = candidate;
                    break;
                }
            }
            result.segmentLabels[segment] = labelsByIndex[selected];
            result.segmentConfidence[segment] = maximumConfidence[selected];
        }

        // The strict-majority remap can place a long original CAD segment across
        // several resampled candidates. Keep radial-facing surfaces only when the
        // original segment normal agrees; otherwise the segment is one of the
        // adjacent sloping/radial walls and must remain ToothWall.
        for (std::size_t segment = 0;
             segment < result.segmentLabels.size();
             ++segment)
        {
            const RegionLabel label = result.segmentLabels[segment];
            if (label != RegionLabel::ToothTop &&
                label != RegionLabel::ToothBottom)
            {
                continue;
            }
            const Eigen::Vector2d direction =
                workingContour.pointsYz[
                    (segment + 1) % workingContour.pointsYz.size()] -
                workingContour.pointsYz[segment];
            const double length = direction.norm();
            if (length <= 1.0e-15)
            {
                result.segmentLabels[segment] = RegionLabel::Unclassified;
                result.segmentConfidence[segment] = 0.0;
                continue;
            }
            const double radialNormalStrength =
                std::abs(direction.y()) / length;
            const double axialNormalStrength =
                std::abs(direction.x()) / length;
            const bool normalDisagrees = label == RegionLabel::ToothTop
                ? radialNormalStrength < 0.58 ||
                    radialNormalStrength <= axialNormalStrength
                : radialNormalStrength < options.axialSurfaceThreshold;
            if (normalDisagrees)
            {
                result.segmentLabels[segment] = RegionLabel::ToothWall;
                result.segmentConfidence[segment] = std::clamp(
                    std::max(
                        options.minimumConfidence,
                        0.55 + 0.40 * axialNormalStrength),
                    0.0,
                    1.0);
            }
        }

        // A coarse original segment at either tooth-group edge can cover both the
        // protected outer wall and the first transition samples. Strict-majority
        // remapping then assigns the whole segment to Transition. Continue the
        // adjacent wall only while its tangent direction remains coherent with the
        // boundary wall, so a gently sampled transition is not consumed piecemeal.
        constexpr double maximumWallContinuationAngle = 0.35;
        const auto originalTangent = [&](std::size_t segment) -> Eigen::Vector2d
        {
            const Eigen::Vector2d direction =
                workingContour.pointsYz[
                    (segment + 1) % workingContour.pointsYz.size()] -
                workingContour.pointsYz[segment];
            const double length = direction.norm();
            if (length <= 1.0e-15)
            {
                return Eigen::Vector2d::Zero();
            }
            return Eigen::Vector2d(direction / length);
        };
        const auto recoverWallContinuation = [&](
            std::size_t wallSegment,
            std::size_t firstTransition,
            bool increasing)
        {
            const Eigen::Vector2d wallTangent = originalTangent(wallSegment);
            if (wallTangent.squaredNorm() <= 0.0)
            {
                return;
            }
            std::size_t segment = firstTransition;
            for (std::size_t visited = 0;
                 visited < result.segmentLabels.size();
                 ++visited)
            {
                if (result.segmentLabels[segment] != RegionLabel::Transition)
                {
                    break;
                }
                const Eigen::Vector2d candidateTangent = originalTangent(segment);
                const double wallNormalStrength = std::abs(candidateTangent.x());
                if (candidateTangent.squaredNorm() <= 0.0 ||
                    wallNormalStrength < options.axialSurfaceThreshold ||
                    tangentAngle(wallTangent, candidateTangent) >
                        maximumWallContinuationAngle)
                {
                    break;
                }

                result.segmentLabels[segment] = RegionLabel::ToothWall;
                result.segmentConfidence[segment] = std::clamp(
                    std::max(
                        options.minimumConfidence,
                        0.55 + 0.40 * wallNormalStrength),
                    0.0,
                    1.0);
                if (increasing)
                {
                    if (!workingContour.closed &&
                        segment + 1 >= result.segmentLabels.size())
                    {
                        break;
                    }
                    segment = (segment + 1) % result.segmentLabels.size();
                }
                else
                {
                    if (!workingContour.closed && segment == 0)
                    {
                        break;
                    }
                    segment = segment == 0
                        ? result.segmentLabels.size() - 1
                        : segment - 1;
                }
            }
        };
        const std::vector<RegionLabel> remappedLabels = result.segmentLabels;
        for (std::size_t segment = 0; segment < remappedLabels.size(); ++segment)
        {
            if (remappedLabels[segment] != RegionLabel::ToothWall)
            {
                continue;
            }
            if ((workingContour.closed || segment + 1 < remappedLabels.size()) &&
                remappedLabels[(segment + 1) % remappedLabels.size()] ==
                    RegionLabel::Transition)
            {
                recoverWallContinuation(
                    segment,
                    (segment + 1) % remappedLabels.size(),
                    true);
            }
            if ((workingContour.closed || segment > 0) &&
                remappedLabels[segment == 0
                    ? remappedLabels.size() - 1
                    : segment - 1] == RegionLabel::Transition)
            {
                recoverWallContinuation(
                    segment,
                    segment == 0 ? remappedLabels.size() - 1 : segment - 1,
                    false);
            }
        }

        // Long original CAD segments may still receive no strict-majority label
        // even though they bridge two surfaces inside the same confirmed tooth
        // period. Restore only short, interior, wall-facing gaps; the large
        // Unclassified body interval and both transition boundaries stay intact.
        for (std::size_t runStart = 0;
             runStart < result.segmentLabels.size();)
        {
            if (result.segmentLabels[runStart] != RegionLabel::Unclassified)
            {
                ++runStart;
                continue;
            }
            std::size_t runEnd = runStart;
            while (runEnd + 1 < result.segmentLabels.size() &&
                result.segmentLabels[runEnd + 1] == RegionLabel::Unclassified)
            {
                ++runEnd;
            }
            const std::size_t runLength = runEnd - runStart + 1;
            if (runLength <= 2 && runStart > 0 &&
                runEnd + 1 < result.segmentLabels.size() &&
                isToothLabel(result.segmentLabels[runStart - 1]) &&
                isToothLabel(result.segmentLabels[runEnd + 1]))
            {
                bool wallFacing = true;
                for (std::size_t segment = runStart; segment <= runEnd; ++segment)
                {
                    const Eigen::Vector2d tangent = originalTangent(segment);
                    wallFacing = wallFacing && tangent.squaredNorm() > 0.0 &&
                        std::abs(tangent.x()) >= 0.45;
                }
                if (wallFacing)
                {
                    for (std::size_t segment = runStart; segment <= runEnd; ++segment)
                    {
                        const double wallNormalStrength =
                            std::abs(originalTangent(segment).x());
                        result.segmentLabels[segment] = RegionLabel::ToothWall;
                        result.segmentConfidence[segment] = std::clamp(
                            std::max(
                                options.minimumConfidence,
                                0.55 + 0.40 * wallNormalStrength),
                            0.0,
                            1.0);
                    }
                }
            }
            runStart = runEnd + 1;
        }

        if (contourRotation != 0)
        {
            undoRotateLeft(result.segmentLabels, contourRotation);
            undoRotateLeft(result.segmentConfidence, contourRotation);
        }

        {
            std::ostringstream diagnostic;
            diagnostic << "Candidate tooth groups: " << transitionSummary.clusterCount
                << ", selected top runs: " << transitionSummary.outerTopRunCount
                << ", qualified outer boundaries: "
                << transitionSummary.qualifiedBoundaryCount
                << ", transition samples: "
                << transitionSummary.markedSegmentCount;
            result.diagnostics.push_back(diagnostic.str());
        }

        for (RegionLabel label : {
                 RegionLabel::ToothTop,
                 RegionLabel::ToothWall,
                 RegionLabel::ToothBottom,
                 RegionLabel::Transition,
                 RegionLabel::Unclassified })
        {
            const std::size_t count = static_cast<std::size_t>(std::count(
                result.segmentLabels.begin(),
                result.segmentLabels.end(),
                label));
            std::ostringstream diagnostic;
            diagnostic << labelName(label) << " segments: " << count;
            result.diagnostics.push_back(diagnostic.str());
        }

        PlanningResult<RegionAssignment> output =
            PlanningResult<RegionAssignment>::success(std::move(result));
        output.diagnostics = output.value.diagnostics;
        return output;
    }
}
