#include <RotationBodyTrajectoryPlanning/RegionPlanning/SprayBoundaryBuilder.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <vector>

namespace smrobot::spray::rotationbody
{
    namespace
    {
        bool isSprayRegion(RegionLabel label) noexcept
        {
            return label == RegionLabel::ToothTop ||
                label == RegionLabel::ToothWall ||
                label == RegionLabel::ToothBottom ||
                label == RegionLabel::Transition;
        }

        void appendUnique(
            std::vector<Eigen::Vector2d>& points,
            const Eigen::Vector2d& point,
            double tolerance)
        {
            const auto duplicate = std::find_if(
                points.begin(),
                points.end(),
                [&](const Eigen::Vector2d& candidate)
                {
                    return (candidate - point).norm() <= tolerance;
                });
            if (duplicate == points.end())
            {
                points.push_back(point);
            }
        }

        bool weightedLineFit(
            const std::vector<Eigen::Vector2d>& points,
            const std::vector<double>& weights,
            double& slope,
            double& intercept) noexcept
        {
            double weightSum = 0.0;
            double weightedZ = 0.0;
            double weightedY = 0.0;
            for (std::size_t index = 0; index < points.size(); ++index)
            {
                weightSum += weights[index];
                weightedZ += weights[index] * points[index].y();
                weightedY += weights[index] * points[index].x();
            }
            if (weightSum <= 0.0)
            {
                return false;
            }
            const double meanZ = weightedZ / weightSum;
            const double meanY = weightedY / weightSum;
            double covariance = 0.0;
            double variance = 0.0;
            for (std::size_t index = 0; index < points.size(); ++index)
            {
                const double centeredZ = points[index].y() - meanZ;
                covariance += weights[index] * centeredZ * (points[index].x() - meanY);
                variance += weights[index] * centeredZ * centeredZ;
            }
            if (variance <= std::numeric_limits<double>::epsilon() * weightSum)
            {
                return false;
            }
            slope = covariance / variance;
            intercept = meanY - slope * meanZ;
            return std::isfinite(slope) && std::isfinite(intercept);
        }

        PlanningResult<std::pair<double, double>> fitOuterEnvelope(
            const std::vector<Eigen::Vector2d>& topPoints,
            double scale)
        {
            if (topPoints.size() < 3)
            {
                return PlanningResult<std::pair<double, double>>::failure(
                    PlanningErrorCode::InsufficientToothTopData,
                    "At least three distinct tooth-top points are required for envelope fitting.");
            }
            double minimumZ = std::numeric_limits<double>::infinity();
            double maximumZ = -std::numeric_limits<double>::infinity();
            for (const Eigen::Vector2d& point : topPoints)
            {
                minimumZ = std::min(minimumZ, point.y());
                maximumZ = std::max(maximumZ, point.y());
            }
            if (maximumZ - minimumZ <= std::max(scale * 1.0e-8, 1.0e-12))
            {
                return PlanningResult<std::pair<double, double>>::failure(
                    PlanningErrorCode::InsufficientToothTopData,
                    "Tooth-top points do not span enough Z distance for envelope fitting.");
            }

            std::vector<double> weights(topPoints.size(), 1.0);
            double slope = 0.0;
            double intercept = 0.0;
            for (int iteration = 0; iteration < 8; ++iteration)
            {
                if (!weightedLineFit(topPoints, weights, slope, intercept))
                {
                    return PlanningResult<std::pair<double, double>>::failure(
                        PlanningErrorCode::SingularFit,
                        "Tooth-top envelope fit is singular.");
                }
                std::vector<double> absoluteResiduals;
                absoluteResiduals.reserve(topPoints.size());
                for (const Eigen::Vector2d& point : topPoints)
                {
                    absoluteResiduals.push_back(
                        std::abs(point.x() - (slope * point.y() + intercept)));
                }
                const auto median = absoluteResiduals.begin() +
                    static_cast<std::ptrdiff_t>(absoluteResiduals.size() / 2);
                std::nth_element(absoluteResiduals.begin(), median, absoluteResiduals.end());
                const double robustScale = std::max(*median * 1.4826, scale * 1.0e-9);
                const double huberThreshold = 1.5 * robustScale;
                for (std::size_t index = 0; index < topPoints.size(); ++index)
                {
                    const double residual = std::abs(
                        topPoints[index].x() - (slope * topPoints[index].y() + intercept));
                    weights[index] = residual <= huberThreshold
                        ? 1.0
                        : huberThreshold / residual;
                }
            }

            double outwardShift = 0.0;
            for (const Eigen::Vector2d& point : topPoints)
            {
                outwardShift = std::max(
                    outwardShift,
                    point.x() - (slope * point.y() + intercept));
            }
            intercept += outwardShift;
            return PlanningResult<std::pair<double, double>>::success({ slope, intercept });
        }
    }

    PlanningResult<SprayBoundary> SprayBoundaryBuilder::build(
        const SectionContour& contour,
        const RegionAssignment& assignment,
        BoundaryMode mode)
    {
        if (contour.segmentCount() == 0)
        {
            return PlanningResult<SprayBoundary>::failure(
                PlanningErrorCode::NoContour,
                "A valid section contour is required to build the spray boundary.");
        }
        if (!assignment.matches(contour))
        {
            return PlanningResult<SprayBoundary>::failure(
                PlanningErrorCode::InvalidArgument,
                "Region labels do not match the section contour segment count.");
        }

        double minimumY = std::numeric_limits<double>::infinity();
        double maximumY = -std::numeric_limits<double>::infinity();
        double minimumZ = std::numeric_limits<double>::infinity();
        double maximumZ = -std::numeric_limits<double>::infinity();
        double maximumTopY = -std::numeric_limits<double>::infinity();
        std::vector<Eigen::Vector2d> topPoints;
        const double scale = std::max(
            contour.toleranceMeters,
            [&]()
            {
                Eigen::AlignedBox2d contourBounds;
                for (const Eigen::Vector2d& point : contour.pointsYz)
                {
                    contourBounds.extend(point);
                }
                return contourBounds.diagonal().norm();
            }());
        const double deduplicationTolerance = std::max(scale * 1.0e-10, 1.0e-12);
        std::size_t spraySegmentCount = 0;
        for (std::size_t segment = 0; segment < contour.segmentCount(); ++segment)
        {
            const RegionLabel label = assignment.segmentLabels[segment];
            if (!isSprayRegion(label))
            {
                continue;
            }
            ++spraySegmentCount;
            const Eigen::Vector2d points[] = {
                contour.pointsYz[segment],
                contour.pointsYz[(segment + 1) % contour.pointsYz.size()]
            };
            for (const Eigen::Vector2d& point : points)
            {
                if (!point.allFinite())
                {
                    return PlanningResult<SprayBoundary>::failure(
                        PlanningErrorCode::NonFiniteGeometry,
                        "Spray region contains a non-finite contour point.");
                }
                minimumY = std::min(minimumY, point.x());
                maximumY = std::max(maximumY, point.x());
                minimumZ = std::min(minimumZ, point.y());
                maximumZ = std::max(maximumZ, point.y());
                if (label == RegionLabel::ToothTop)
                {
                    maximumTopY = std::max(maximumTopY, point.x());
                    appendUnique(topPoints, point, deduplicationTolerance);
                }
            }
        }
        if (spraySegmentCount == 0)
        {
            return PlanningResult<SprayBoundary>::failure(
                PlanningErrorCode::InsufficientRegionData,
                "No classified spray region is available for boundary generation.");
        }
        if (topPoints.empty() || !std::isfinite(maximumTopY))
        {
            return PlanningResult<SprayBoundary>::failure(
                PlanningErrorCode::InsufficientToothTopData,
                "At least one ToothTop segment is required for boundary generation.");
        }
        if (maximumZ - minimumZ <= deduplicationTolerance)
        {
            return PlanningResult<SprayBoundary>::failure(
                PlanningErrorCode::DegenerateGeometry,
                "Spray regions do not span a non-zero Z range.");
        }

        SprayBoundary boundary;
        boundary.mode = mode;
        boundary.minimumY = minimumY;
        boundary.maximumY = maximumY;
        boundary.minimumZ = minimumZ;
        boundary.maximumZ = maximumZ;

        if (mode == BoundaryMode::MaximumToothTopY)
        {
            boundary.maximumY = maximumTopY;
            boundary.outerLineSlopeYPerZ = 0.0;
            boundary.outerLineInterceptY = maximumTopY;
        }
        else
        {
            auto fit = fitOuterEnvelope(topPoints, scale);
            if (!fit)
            {
                return PlanningResult<SprayBoundary>::failure(fit.error.code, fit.error.message);
            }
            boundary.outerLineSlopeYPerZ = fit.value.first;
            boundary.outerLineInterceptY = fit.value.second;
            boundary.maximumY = std::max(
                boundary.outerY(minimumZ),
                boundary.outerY(maximumZ));
            if (boundary.outerY(minimumZ) < minimumY - deduplicationTolerance ||
                boundary.outerY(maximumZ) < minimumY - deduplicationTolerance)
            {
                return PlanningResult<SprayBoundary>::failure(
                    PlanningErrorCode::SingularFit,
                    "Fitted tooth-top envelope crosses the inner spray boundary.");
            }
        }

        boundary.polygonYz = {
            { minimumY, minimumZ },
            { boundary.outerY(minimumZ), minimumZ },
            { boundary.outerY(maximumZ), maximumZ },
            { minimumY, maximumZ }
        };
        return PlanningResult<SprayBoundary>::success(std::move(boundary));
    }
}
