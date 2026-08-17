#include <RotationBodyTrajectoryPlanning/Sectioning/ContourTopology.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>

namespace smrobot::spray::rotationbody
{
    namespace
    {
        struct QuantizedPoint
        {
            std::int64_t y{ 0 };
            std::int64_t z{ 0 };

            bool operator==(const QuantizedPoint& other) const noexcept
            {
                return y == other.y && z == other.z;
            }
        };

        struct QuantizedPointHash
        {
            std::size_t operator()(const QuantizedPoint& point) const noexcept
            {
                const auto first = std::hash<std::int64_t>{}(point.y);
                const auto second = std::hash<std::int64_t>{}(point.z);
                return first ^ (second + 0x9e3779b9U + (first << 6U) + (first >> 2U));
            }
        };

        struct EdgeKey
        {
            std::size_t first{ 0 };
            std::size_t second{ 0 };

            bool operator==(const EdgeKey& other) const noexcept
            {
                return first == other.first && second == other.second;
            }
        };

        struct EdgeKeyHash
        {
            std::size_t operator()(const EdgeKey& edge) const noexcept
            {
                const auto first = std::hash<std::size_t>{}(edge.first);
                const auto second = std::hash<std::size_t>{}(edge.second);
                return first ^ (second + 0x9e3779b9U + (first << 6U) + (first >> 2U));
            }
        };

        struct Node
        {
            Eigen::Vector2d representative = Eigen::Vector2d::Zero();
            Eigen::Vector2d accumulated = Eigen::Vector2d::Zero();
            std::size_t sampleCount{ 0 };
            std::vector<std::size_t> incidentEdges;

            Eigen::Vector2d position() const noexcept
            {
                return sampleCount == 0
                    ? representative
                    : accumulated / static_cast<double>(sampleCount);
            }
        };

        struct Edge
        {
            std::size_t firstNode{ 0 };
            std::size_t secondNode{ 0 };
        };

        bool quantizeCoordinate(double value, double tolerance, std::int64_t& output) noexcept
        {
            const long double scaled = std::floor(
                static_cast<long double>(value) / static_cast<long double>(tolerance));
            const auto minimum = static_cast<long double>(std::numeric_limits<std::int64_t>::min() + 2);
            const auto maximum = static_cast<long double>(std::numeric_limits<std::int64_t>::max() - 2);
            if (!std::isfinite(scaled) || scaled < minimum || scaled > maximum)
            {
                return false;
            }
            output = static_cast<std::int64_t>(scaled);
            return true;
        }

        bool quantizePoint(
            const Eigen::Vector2d& point,
            double tolerance,
            QuantizedPoint& output) noexcept
        {
            return quantizeCoordinate(point.x(), tolerance, output.y)
                && quantizeCoordinate(point.y(), tolerance, output.z);
        }

        double squaredDistance(const Eigen::Vector2d& first, const Eigen::Vector2d& second) noexcept
        {
            return (first - second).squaredNorm();
        }

        double signedArea(const std::vector<Eigen::Vector2d>& points) noexcept
        {
            if (points.size() < 3)
            {
                return 0.0;
            }

            double twiceArea = 0.0;
            for (std::size_t index = 0; index < points.size(); ++index)
            {
                const auto next = (index + 1) % points.size();
                twiceArea += points[index].x() * points[next].y()
                    - points[next].x() * points[index].y();
            }
            return 0.5 * twiceArea;
        }

        double contourLength(const SectionContour& contour) noexcept
        {
            double length = 0.0;
            for (std::size_t index = 1; index < contour.pointsYz.size(); ++index)
            {
                length += (contour.pointsYz[index] - contour.pointsYz[index - 1]).norm();
            }
            if (contour.closed && contour.pointsYz.size() > 2)
            {
                length += (contour.pointsYz.front() - contour.pointsYz.back()).norm();
            }
            return length;
        }

        bool pointOrder(const Eigen::Vector2d& first, const Eigen::Vector2d& second) noexcept
        {
            if (first.y() != second.y())
            {
                return first.y() < second.y();
            }
            return first.x() < second.x();
        }

        void normalizePointOrder(std::vector<Eigen::Vector2d>& points, bool closed)
        {
            if (points.empty())
            {
                return;
            }

            if (!closed)
            {
                if (pointOrder(points.back(), points.front()))
                {
                    std::reverse(points.begin(), points.end());
                }
                return;
            }

            const auto start = std::min_element(points.begin(), points.end(), pointOrder);
            std::rotate(points.begin(), start, points.end());
            if (signedArea(points) < 0.0 && points.size() > 2)
            {
                std::reverse(points.begin() + 1, points.end());
            }
        }

        SectionContour makeContour(
            const std::vector<std::size_t>& nodePath,
            bool closed,
            const std::vector<Node>& nodes,
            double tolerance,
            double planeX)
        {
            SectionContour contour;
            contour.toleranceMeters = tolerance;
            contour.closed = closed;
            contour.pointsYz.reserve(nodePath.size());

            const double toleranceSquared = tolerance * tolerance;
            for (const auto nodeIndex : nodePath)
            {
                const Eigen::Vector2d point = nodes[nodeIndex].position();
                if (contour.pointsYz.empty()
                    || squaredDistance(contour.pointsYz.back(), point) > toleranceSquared)
                {
                    contour.pointsYz.push_back(point);
                }
            }

            if (contour.pointsYz.size() > 2
                && squaredDistance(contour.pointsYz.front(), contour.pointsYz.back()) <= toleranceSquared)
            {
                contour.pointsYz.pop_back();
                contour.closed = true;
            }
            if (contour.closed && contour.pointsYz.size() < 3)
            {
                contour.closed = false;
            }

            normalizePointOrder(contour.pointsYz, contour.closed);
            contour.points3d.reserve(contour.pointsYz.size());
            contour.cumulativeArcLength.reserve(contour.pointsYz.size());

            double length = 0.0;
            for (std::size_t index = 0; index < contour.pointsYz.size(); ++index)
            {
                if (index > 0)
                {
                    length += (contour.pointsYz[index] - contour.pointsYz[index - 1]).norm();
                }
                contour.cumulativeArcLength.push_back(length);
                contour.points3d.emplace_back(
                    planeX,
                    contour.pointsYz[index].x(),
                    contour.pointsYz[index].y());
            }

            std::ostringstream diagnostic;
            diagnostic << (contour.closed ? "Closed" : "Open")
                << " contour assembled with " << contour.segmentCount() << " segments.";
            contour.diagnostics.push_back(diagnostic.str());
            return contour;
        }

        struct ContourMetric
        {
            std::size_t index{ 0 };
            double minimumY{ std::numeric_limits<double>::infinity() };
            double maximumY{ -std::numeric_limits<double>::infinity() };
            double minimumZ{ std::numeric_limits<double>::infinity() };
            double maximumZ{ -std::numeric_limits<double>::infinity() };
            double meanY{ 0.0 };
            double length{ 0.0 };
            double area{ 0.0 };
            bool eligible{ false };
        };

        ContourMetric measureContour(
            const SectionContour& contour,
            std::size_t index,
            double tolerance) noexcept
        {
            ContourMetric metric;
            metric.index = index;
            if (contour.pointsYz.size() < 2)
            {
                return metric;
            }

            bool hasNegativePoint = false;
            for (const auto& point : contour.pointsYz)
            {
                if (!point.allFinite())
                {
                    return metric;
                }
                metric.minimumY = std::min(metric.minimumY, point.x());
                metric.maximumY = std::max(metric.maximumY, point.x());
                metric.minimumZ = std::min(metric.minimumZ, point.y());
                metric.maximumZ = std::max(metric.maximumZ, point.y());
                metric.meanY += point.x();
                hasNegativePoint = hasNegativePoint || point.x() < -tolerance;
            }
            metric.meanY /= static_cast<double>(contour.pointsYz.size());
            metric.length = contourLength(contour);
            metric.area = contour.closed ? std::abs(signedArea(contour.pointsYz)) : 0.0;
            metric.eligible = !hasNegativePoint
                && metric.maximumY > tolerance
                && metric.length > tolerance;
            return metric;
        }

        double normalized(double value, double maximum) noexcept
        {
            return maximum > 0.0 ? std::clamp(value / maximum, 0.0, 1.0) : 0.0;
        }
    }

    PlanningResult<ContourTopologyResult> ContourTopology::buildContours(
        const std::vector<YzSegment>& segments,
        double weldToleranceMeters,
        double planeX)
    {
        if (!std::isfinite(weldToleranceMeters) || weldToleranceMeters <= 0.0
            || !std::isfinite(planeX))
        {
            return PlanningResult<ContourTopologyResult>::failure(
                PlanningErrorCode::InvalidArgument,
                "Contour weld tolerance and plane position must be finite and positive.");
        }
        if (segments.empty())
        {
            return PlanningResult<ContourTopologyResult>::failure(
                PlanningErrorCode::NoContour,
                "No section segments were supplied for contour assembly.");
        }

        ContourTopologyResult output;
        output.inputSegmentCount = segments.size();

        std::vector<Node> nodes;
        std::vector<Edge> edges;
        std::unordered_map<QuantizedPoint, std::vector<std::size_t>, QuantizedPointHash> buckets;
        std::unordered_set<EdgeKey, EdgeKeyHash> uniqueEdges;
        const double toleranceSquared = weldToleranceMeters * weldToleranceMeters;

        auto findOrCreateNode = [&](const Eigen::Vector2d& point, std::size_t& nodeIndex) -> bool
        {
            QuantizedPoint center;
            if (!quantizePoint(point, weldToleranceMeters, center))
            {
                return false;
            }

            double bestDistance = std::numeric_limits<double>::infinity();
            std::size_t bestNode = std::numeric_limits<std::size_t>::max();
            for (std::int64_t yOffset = -1; yOffset <= 1; ++yOffset)
            {
                for (std::int64_t zOffset = -1; zOffset <= 1; ++zOffset)
                {
                    const QuantizedPoint neighbour{ center.y + yOffset, center.z + zOffset };
                    const auto bucket = buckets.find(neighbour);
                    if (bucket == buckets.end())
                    {
                        continue;
                    }
                    for (const auto candidate : bucket->second)
                    {
                        const double distance = squaredDistance(nodes[candidate].representative, point);
                        if (distance <= toleranceSquared
                            && (distance < bestDistance
                                || (distance == bestDistance && candidate < bestNode)))
                        {
                            bestDistance = distance;
                            bestNode = candidate;
                        }
                    }
                }
            }

            if (bestNode != std::numeric_limits<std::size_t>::max())
            {
                nodes[bestNode].accumulated += point;
                ++nodes[bestNode].sampleCount;
                nodeIndex = bestNode;
                return true;
            }

            Node node;
            node.representative = point;
            node.accumulated = point;
            node.sampleCount = 1;
            nodeIndex = nodes.size();
            nodes.push_back(std::move(node));
            buckets[center].push_back(nodeIndex);
            return true;
        };

        for (const auto& segment : segments)
        {
            if (!segment.first.allFinite() || !segment.second.allFinite())
            {
                return PlanningResult<ContourTopologyResult>::failure(
                    PlanningErrorCode::NonFiniteGeometry,
                    "A section segment contains a non-finite endpoint.");
            }
            if (squaredDistance(segment.first, segment.second) <= toleranceSquared)
            {
                ++output.discardedDegenerateSegmentCount;
                continue;
            }

            std::size_t firstNode = 0;
            std::size_t secondNode = 0;
            if (!findOrCreateNode(segment.first, firstNode)
                || !findOrCreateNode(segment.second, secondNode))
            {
                return PlanningResult<ContourTopologyResult>::failure(
                    PlanningErrorCode::DegenerateGeometry,
                    "Section coordinates exceed the safe quantization range.");
            }
            if (firstNode == secondNode)
            {
                ++output.discardedDegenerateSegmentCount;
                continue;
            }

            const EdgeKey edgeKey{
                std::min(firstNode, secondNode),
                std::max(firstNode, secondNode) };
            if (!uniqueEdges.insert(edgeKey).second)
            {
                ++output.discardedDuplicateSegmentCount;
                continue;
            }
            edges.push_back({ firstNode, secondNode });
        }

        if (edges.empty())
        {
            return PlanningResult<ContourTopologyResult>::failure(
                PlanningErrorCode::NoContour,
                "All section segments were duplicate or degenerate.");
        }

        for (std::size_t edgeIndex = 0; edgeIndex < edges.size(); ++edgeIndex)
        {
            nodes[edges[edgeIndex].firstNode].incidentEdges.push_back(edgeIndex);
            nodes[edges[edgeIndex].secondNode].incidentEdges.push_back(edgeIndex);
        }
        for (auto& node : nodes)
        {
            std::sort(node.incidentEdges.begin(), node.incidentEdges.end());
            if (node.incidentEdges.size() > 2)
            {
                ++output.branchNodeCount;
            }
        }

        output.uniqueSegmentCount = edges.size();
        std::vector<bool> edgeUsed(edges.size(), false);

        auto oppositeNode = [&](std::size_t edgeIndex, std::size_t nodeIndex) noexcept
        {
            return edges[edgeIndex].firstNode == nodeIndex
                ? edges[edgeIndex].secondNode
                : edges[edgeIndex].firstNode;
        };

        auto trace = [&](std::size_t startNode, std::size_t firstEdge)
        {
            std::vector<std::size_t> path;
            path.push_back(startNode);
            std::size_t currentNode = startNode;
            std::size_t currentEdge = firstEdge;
            bool closed = false;

            while (!edgeUsed[currentEdge])
            {
                edgeUsed[currentEdge] = true;
                const std::size_t nextNode = oppositeNode(currentEdge, currentNode);
                path.push_back(nextNode);
                if (nextNode == startNode)
                {
                    path.pop_back();
                    closed = true;
                    break;
                }
                if (nodes[nextNode].incidentEdges.size() != 2)
                {
                    break;
                }

                std::size_t nextEdge = std::numeric_limits<std::size_t>::max();
                for (const auto candidate : nodes[nextNode].incidentEdges)
                {
                    if (!edgeUsed[candidate])
                    {
                        nextEdge = candidate;
                        break;
                    }
                }
                if (nextEdge == std::numeric_limits<std::size_t>::max())
                {
                    break;
                }
                currentNode = nextNode;
                currentEdge = nextEdge;
            }

            return std::make_pair(std::move(path), closed);
        };

        for (std::size_t nodeIndex = 0; nodeIndex < nodes.size(); ++nodeIndex)
        {
            if (nodes[nodeIndex].incidentEdges.size() == 2)
            {
                continue;
            }
            for (const auto edgeIndex : nodes[nodeIndex].incidentEdges)
            {
                if (edgeUsed[edgeIndex])
                {
                    continue;
                }
                auto traced = trace(nodeIndex, edgeIndex);
                SectionContour contour = makeContour(
                    traced.first,
                    traced.second,
                    nodes,
                    weldToleranceMeters,
                    planeX);
                if (contour.segmentCount() > 0
                    && contourLength(contour) > weldToleranceMeters)
                {
                    output.contours.push_back(std::move(contour));
                }
            }
        }

        for (std::size_t edgeIndex = 0; edgeIndex < edges.size(); ++edgeIndex)
        {
            if (edgeUsed[edgeIndex])
            {
                continue;
            }
            auto traced = trace(edges[edgeIndex].firstNode, edgeIndex);
            SectionContour contour = makeContour(
                traced.first,
                traced.second,
                nodes,
                weldToleranceMeters,
                planeX);
            if (contour.segmentCount() > 0
                && contourLength(contour) > weldToleranceMeters)
            {
                output.contours.push_back(std::move(contour));
            }
        }

        if (output.contours.empty())
        {
            return PlanningResult<ContourTopologyResult>::failure(
                PlanningErrorCode::NoContour,
                "Section segments could not be assembled into a usable contour.");
        }

        std::stable_sort(output.contours.begin(), output.contours.end(),
            [](const SectionContour& first, const SectionContour& second)
            {
                if (first.closed != second.closed)
                {
                    return first.closed;
                }
                const double firstLength = contourLength(first);
                const double secondLength = contourLength(second);
                if (firstLength != secondLength)
                {
                    return firstLength > secondLength;
                }
                const auto firstMinimum = *std::min_element(
                    first.pointsYz.begin(), first.pointsYz.end(), pointOrder);
                const auto secondMinimum = *std::min_element(
                    second.pointsYz.begin(), second.pointsYz.end(), pointOrder);
                return pointOrder(firstMinimum, secondMinimum);
            });

        PlanningResult<ContourTopologyResult> result =
            PlanningResult<ContourTopologyResult>::success(std::move(output));
        std::ostringstream summary;
        summary << "Contour graph: " << result.value.inputSegmentCount << " input, "
            << result.value.uniqueSegmentCount << " unique, "
            << result.value.contours.size() << " contours, "
            << result.value.branchNodeCount << " branch nodes.";
        result.diagnostics.push_back(summary.str());
        return result;
    }

    PlanningResult<std::size_t> ContourTopology::selectTargetOuterContour(
        const std::vector<SectionContour>& contours,
        double toleranceMeters)
    {
        if (!std::isfinite(toleranceMeters) || toleranceMeters <= 0.0)
        {
            return PlanningResult<std::size_t>::failure(
                PlanningErrorCode::InvalidArgument,
                "Target contour tolerance must be finite and positive.");
        }
        if (contours.empty())
        {
            return PlanningResult<std::size_t>::failure(
                PlanningErrorCode::NoContour,
                "No contour candidates are available for outer-contour selection.");
        }

        std::vector<ContourMetric> metrics;
        metrics.reserve(contours.size());
        double globalMinimumY = std::numeric_limits<double>::infinity();
        double globalMaximumY = -std::numeric_limits<double>::infinity();
        double globalMinimumZ = std::numeric_limits<double>::infinity();
        double globalMaximumZ = -std::numeric_limits<double>::infinity();
        double maximumLength = 0.0;
        double maximumArea = 0.0;

        for (std::size_t index = 0; index < contours.size(); ++index)
        {
            ContourMetric metric = measureContour(contours[index], index, toleranceMeters);
            if (metric.eligible)
            {
                globalMinimumY = std::min(globalMinimumY, metric.minimumY);
                globalMaximumY = std::max(globalMaximumY, metric.maximumY);
                globalMinimumZ = std::min(globalMinimumZ, metric.minimumZ);
                globalMaximumZ = std::max(globalMaximumZ, metric.maximumZ);
                maximumLength = std::max(maximumLength, metric.length);
                maximumArea = std::max(maximumArea, metric.area);
            }
            metrics.push_back(metric);
        }

        const auto eligibleCount = std::count_if(metrics.begin(), metrics.end(),
            [](const ContourMetric& metric) { return metric.eligible; });
        if (eligibleCount == 0)
        {
            return PlanningResult<std::size_t>::failure(
                PlanningErrorCode::NoContour,
                "No contour has a non-degenerate extent in the positive-Y half-plane.");
        }

        constexpr std::size_t envelopeBinCount = 48;
        const double zSpan = std::max(globalMaximumZ - globalMinimumZ, toleranceMeters);
        const double ySpan = std::max(globalMaximumY - globalMinimumY, toleranceMeters);
        const double unavailable = -std::numeric_limits<double>::infinity();
        std::vector<std::array<double, envelopeBinCount>> envelopes(metrics.size());
        std::array<double, envelopeBinCount> globalEnvelope;
        globalEnvelope.fill(unavailable);

        auto binCenter = [&](std::size_t bin) noexcept
        {
            return globalMinimumZ
                + (static_cast<double>(bin) + 0.5) * zSpan
                    / static_cast<double>(envelopeBinCount);
        };

        for (std::size_t metricIndex = 0; metricIndex < metrics.size(); ++metricIndex)
        {
            envelopes[metricIndex].fill(unavailable);
            if (!metrics[metricIndex].eligible)
            {
                continue;
            }

            const auto& contour = contours[metricIndex];
            const std::size_t segmentCount = contour.segmentCount();
            for (std::size_t segmentIndex = 0; segmentIndex < segmentCount; ++segmentIndex)
            {
                const auto next = segmentIndex + 1 < contour.pointsYz.size()
                    ? segmentIndex + 1
                    : 0;
                const Eigen::Vector2d& first = contour.pointsYz[segmentIndex];
                const Eigen::Vector2d& second = contour.pointsYz[next];
                const double minimumZ = std::min(first.y(), second.y()) - toleranceMeters;
                const double maximumZ = std::max(first.y(), second.y()) + toleranceMeters;
                const double deltaZ = second.y() - first.y();

                for (std::size_t bin = 0; bin < envelopeBinCount; ++bin)
                {
                    const double z = binCenter(bin);
                    if (z < minimumZ || z > maximumZ)
                    {
                        continue;
                    }

                    double y = std::max(first.x(), second.x());
                    if (std::abs(deltaZ) > toleranceMeters)
                    {
                        const double interpolation = std::clamp(
                            (z - first.y()) / deltaZ,
                            0.0,
                            1.0);
                        y = first.x() + interpolation * (second.x() - first.x());
                    }
                    envelopes[metricIndex][bin] = std::max(envelopes[metricIndex][bin], y);
                    globalEnvelope[bin] = std::max(globalEnvelope[bin], y);
                }
            }
        }

        const std::size_t globalEnvelopeBins = static_cast<std::size_t>(std::count_if(
            globalEnvelope.begin(), globalEnvelope.end(),
            [&](double value) { return value != unavailable; }));

        double bestScore = -std::numeric_limits<double>::infinity();
        std::size_t bestIndex = 0;
        std::vector<std::string> scoreDiagnostics;
        for (std::size_t metricIndex = 0; metricIndex < metrics.size(); ++metricIndex)
        {
            const auto& metric = metrics[metricIndex];
            if (!metric.eligible)
            {
                continue;
            }

            std::size_t envelopeMatches = 0;
            for (std::size_t bin = 0; bin < envelopeBinCount; ++bin)
            {
                if (envelopes[metricIndex][bin] != unavailable
                    && globalEnvelope[bin] != unavailable
                    && globalEnvelope[bin] - envelopes[metricIndex][bin]
                        <= 4.0 * toleranceMeters)
                {
                    ++envelopeMatches;
                }
            }

            std::size_t dominatedContours = 0;
            for (const auto& other : metrics)
            {
                if (!other.eligible || other.index == metric.index)
                {
                    continue;
                }
                const bool spansOther = metric.minimumZ <= other.minimumZ + toleranceMeters
                    && metric.maximumZ >= other.maximumZ - toleranceMeters;
                if (spansOther && metric.maximumY >= other.maximumY - toleranceMeters)
                {
                    ++dominatedContours;
                }
            }

            const double envelopeCoverage = globalEnvelopeBins == 0
                ? 0.0
                : static_cast<double>(envelopeMatches)
                    / static_cast<double>(globalEnvelopeBins);
            const double zCoverage = normalized(metric.maximumZ - metric.minimumZ, zSpan);
            const double maximumYScore = std::clamp(
                (metric.maximumY - globalMinimumY) / ySpan,
                0.0,
                1.0);
            const double meanYScore = std::clamp(
                (metric.meanY - globalMinimumY) / ySpan,
                0.0,
                1.0);
            const double containment = eligibleCount <= 1
                ? 1.0
                : static_cast<double>(dominatedContours)
                    / static_cast<double>(eligibleCount - 1);
            const double score = 10.0 * envelopeCoverage
                + 2.5 * zCoverage
                + 1.5 * maximumYScore
                + meanYScore
                + 0.75 * normalized(metric.length, maximumLength)
                + 0.5 * normalized(metric.area, maximumArea)
                + 0.35 * containment
                + (contours[metricIndex].closed ? 0.15 : 0.0);

            std::ostringstream detail;
            detail << std::fixed << std::setprecision(4)
                << "Contour " << metric.index << " outer score=" << score
                << ", envelope=" << envelopeCoverage
                << ", zCoverage=" << zCoverage
                << ", meanY=" << metric.meanY
                << ", area=" << metric.area
                << ", length=" << metric.length
                << ", closed=" << (contours[metricIndex].closed ? "true" : "false") << '.';
            scoreDiagnostics.push_back(detail.str());

            if (score > bestScore
                || (score == bestScore && metric.index < bestIndex))
            {
                bestScore = score;
                bestIndex = metric.index;
            }
        }

        PlanningResult<std::size_t> result = PlanningResult<std::size_t>::success(bestIndex);
        result.diagnostics = std::move(scoreDiagnostics);
        std::ostringstream selected;
        selected << "Selected contour " << bestIndex
            << " as the positive-Y outer contour with score "
            << std::fixed << std::setprecision(4) << bestScore << '.';
        result.diagnostics.push_back(selected.str());
        return result;
    }
}
