#include <RotationBodyTrajectoryPlanning/TrajectoryPlanning/AutomaticTrajectoryPlanner.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>

namespace smrobot::spray::rotationbody
{
    namespace
    {
        constexpr double epsilon = 1.0e-12;

        struct RegionAccumulator
        {
            Eigen::Vector2d weightedNormal = Eigen::Vector2d::Zero();
            Eigen::Vector2d weightedMidpoint = Eigen::Vector2d::Zero();
            double totalLength{ 0.0 };

            void add(
                const Eigen::Vector2d& normal,
                const Eigen::Vector2d& midpoint,
                double length)
            {
                weightedNormal += length * normal;
                weightedMidpoint += length * midpoint;
                totalLength += length;
            }
        };

        PlanningResult<double> meanNormalAngle(
            const RegionAccumulator& accumulator,
            const char* missingMessage)
        {
            if(!std::isfinite(accumulator.totalLength) ||
                accumulator.totalLength <= epsilon ||
                !accumulator.weightedNormal.allFinite() ||
                accumulator.weightedNormal.norm() <=
                    epsilon * accumulator.totalLength) {
                return PlanningResult<double>::failure(
                    PlanningErrorCode::InsufficientRegionData,
                    missingMessage);
            }
            return PlanningResult<double>::success(std::atan2(
                accumulator.weightedNormal.y(),
                accumulator.weightedNormal.x()));
        }

        Eigen::Vector2d closestPointOnSegment(
            const Eigen::Vector2d& point,
            const Eigen::Vector2d& start,
            const Eigen::Vector2d& end)
        {
            const Eigen::Vector2d delta = end - start;
            const double lengthSquared = delta.squaredNorm();
            if(lengthSquared <= epsilon) {
                return start;
            }
            const double ratio = std::clamp(
                (point - start).dot(delta) / lengthSquared,
                0.0,
                1.0);
            return start + ratio * delta;
        }

        std::vector<Eigen::Vector2d> toothBottomRunCenters(
            const SectionContour& contour,
            const RegionAssignment& regions)
        {
            const std::size_t segmentCount = contour.segmentCount();
            std::vector<Eigen::Vector2d> result;
            std::size_t firstNonBottom = segmentCount;
            for(std::size_t segment = 0; segment < segmentCount; ++segment) {
                if(regions.segmentLabels[segment] != RegionLabel::ToothBottom) {
                    firstNonBottom = segment;
                    break;
                }
            }

            RegionAccumulator run;
            const std::size_t start = firstNonBottom == segmentCount
                ? 0
                : (firstNonBottom + 1) % segmentCount;
            for(std::size_t offset = 0; offset < segmentCount; ++offset) {
                const std::size_t segment = (start + offset) % segmentCount;
                if(regions.segmentLabels[segment] == RegionLabel::ToothBottom) {
                    const Eigen::Vector2d& first = contour.pointsYz[segment];
                    const Eigen::Vector2d& second =
                        contour.pointsYz[(segment + 1) % contour.pointsYz.size()];
                    const double length = (second - first).norm();
                    run.add(Eigen::Vector2d::Zero(), 0.5 * (first + second), length);
                } else if(run.totalLength > epsilon) {
                    result.push_back(run.weightedMidpoint / run.totalLength);
                    run = {};
                }
            }
            if(run.totalLength > epsilon) {
                result.push_back(run.weightedMidpoint / run.totalLength);
            }
            return result;
        }

        TrajectoryGenerationParameters parameters(double tiltRadians)
        {
            TrajectoryGenerationParameters result;
            result.sprayDistanceMeters = 0.110;
            result.tiltRadians = tiltRadians;
            result.speedMetersPerSecond = 0.004;
            result.startExtensionMeters = 0.015;
            result.endExtensionMeters = 0.015;
            result.pointCount = 0;
            result.positionerRpm = 65.0;
            result.reversed = false;
            return result;
        }
    }

    PlanningResult<AutomaticTrajectoryPlan> AutomaticTrajectoryPlanner::plan(
        const SectionContour& contour,
        const RegionAssignment& regions,
        AutomaticTrajectoryMode mode)
    {
        if(!contour.closed || contour.segmentCount() < 3 ||
            !regions.matches(contour)) {
            return PlanningResult<AutomaticTrajectoryPlan>::failure(
                PlanningErrorCode::InvalidArgument,
                "Automatic trajectory planning requires a closed classified section contour.");
        }

        double signedDoubleArea = 0.0;
        for(std::size_t index = 0; index < contour.pointsYz.size(); ++index) {
            const Eigen::Vector2d& start = contour.pointsYz[index];
            const Eigen::Vector2d& end =
                contour.pointsYz[(index + 1) % contour.pointsYz.size()];
            if(!start.allFinite() || !end.allFinite()) {
                return PlanningResult<AutomaticTrajectoryPlan>::failure(
                    PlanningErrorCode::NonFiniteGeometry,
                    "The section contour contains non-finite points.");
            }
            signedDoubleArea += start.x() * end.y() - end.x() * start.y();
        }
        if(!std::isfinite(signedDoubleArea) || std::abs(signedDoubleArea) <= epsilon) {
            return PlanningResult<AutomaticTrajectoryPlan>::failure(
                PlanningErrorCode::DegenerateGeometry,
                "The section contour orientation cannot be determined.");
        }
        const bool counterClockwise = signedDoubleArea > 0.0;

        RegionAccumulator upperWall;
        RegionAccumulator lowerWall;
        RegionAccumulator toothTop;
        RegionAccumulator toothBottom;
        for(std::size_t segment = 0; segment < contour.segmentCount(); ++segment) {
            const Eigen::Vector2d& start = contour.pointsYz[segment];
            const Eigen::Vector2d& end =
                contour.pointsYz[(segment + 1) % contour.pointsYz.size()];
            const Eigen::Vector2d delta = end - start;
            const double length = delta.norm();
            if(!std::isfinite(length) || length <= epsilon) {
                return PlanningResult<AutomaticTrajectoryPlan>::failure(
                    PlanningErrorCode::DegenerateGeometry,
                    "The section contour contains a zero-length segment.");
            }
            const Eigen::Vector2d tangent = delta / length;
            const Eigen::Vector2d outwardNormal = counterClockwise
                ? Eigen::Vector2d(tangent.y(), -tangent.x())
                : Eigen::Vector2d(-tangent.y(), tangent.x());
            const Eigen::Vector2d midpoint = 0.5 * (start + end);

            switch(regions.segmentLabels[segment]) {
            case RegionLabel::ToothTop:
                toothTop.add(outwardNormal, midpoint, length);
                break;
            case RegionLabel::ToothBottom:
                toothBottom.add(outwardNormal, midpoint, length);
                break;
            case RegionLabel::ToothWall:
                if(outwardNormal.y() > epsilon) {
                    upperWall.add(outwardNormal, midpoint, length);
                } else if(outwardNormal.y() < -epsilon) {
                    lowerWall.add(outwardNormal, midpoint, length);
                }
                break;
            case RegionLabel::Unclassified:
            case RegionLabel::Transition:
                break;
            }
        }

        const PlanningResult<double> upperAngle = meanNormalAngle(
            upperWall,
            "No upper tooth-wall segments with outward normal toward +Z were found.");
        const PlanningResult<double> lowerAngle = meanNormalAngle(
            lowerWall,
            "No lower tooth-wall segments with outward normal toward -Z were found.");
        const PlanningResult<double> topAngle = meanNormalAngle(
            toothTop,
            "No usable tooth-top segments were found.");
        const PlanningResult<double> bottomAngle = meanNormalAngle(
            toothBottom,
            "No usable tooth-bottom segments were found.");
        for(const PlanningResult<double>* result : {
            &upperAngle, &lowerAngle, &topAngle, &bottomAngle }) {
            if(!*result) {
                return PlanningResult<AutomaticTrajectoryPlan>::failure(
                    result->error.code,
                    result->error.message);
            }
        }

        const std::vector<Eigen::Vector2d> bottomCenters =
            toothBottomRunCenters(contour, regions);
        double dualTiltLimit = std::numeric_limits<double>::infinity();
        for(const Eigen::Vector2d& bottomCenter : bottomCenters) {
            Eigen::Vector2d closestTopPoint = Eigen::Vector2d::Zero();
            double closestDistanceSquared = std::numeric_limits<double>::infinity();
            for(std::size_t segment = 0; segment < contour.segmentCount(); ++segment) {
                if(regions.segmentLabels[segment] != RegionLabel::ToothTop) {
                    continue;
                }
                const Eigen::Vector2d& start = contour.pointsYz[segment];
                const Eigen::Vector2d& end =
                    contour.pointsYz[(segment + 1) % contour.pointsYz.size()];
                const Eigen::Vector2d candidate =
                    closestPointOnSegment(bottomCenter, start, end);
                const double distanceSquared =
                    (candidate - bottomCenter).squaredNorm();
                if(distanceSquared < closestDistanceSquared) {
                    closestDistanceSquared = distanceSquared;
                    closestTopPoint = candidate;
                }
            }
            const Eigen::Vector2d limitVector = closestTopPoint - bottomCenter;
            if(limitVector.allFinite() && limitVector.norm() > epsilon) {
                dualTiltLimit = std::min(
                    dualTiltLimit,
                    std::atan2(
                        std::abs(limitVector.x()),
                        std::abs(limitVector.y())));
            }
        }
        if(!std::isfinite(dualTiltLimit)) {
            return PlanningResult<AutomaticTrajectoryPlan>::failure(
                PlanningErrorCode::DegenerateGeometry,
                "No distinct tooth-top and tooth-bottom reference points were found.");
        }

        AutomaticTrajectoryPlan result;
        result.upperWallNormalAngleRadians = upperAngle.value;
        result.lowerWallNormalAngleRadians = lowerAngle.value;
        result.toothTopNormalAngleRadians = topAngle.value;
        result.toothBottomNormalAngleRadians = bottomAngle.value;
        result.dualTiltLimitRadians = dualTiltLimit;

        if(mode == AutomaticTrajectoryMode::Dual) {
            const double upperTilt = 0.4 * upperAngle.value +
                0.3 * topAngle.value + 0.3 * bottomAngle.value;
            const double lowerTilt = 0.4 * lowerAngle.value +
                0.3 * topAngle.value + 0.3 * bottomAngle.value;
            result.trajectories = {
                parameters(std::clamp(upperTilt, -dualTiltLimit, dualTiltLimit)),
                parameters(std::clamp(lowerTilt, -dualTiltLimit, dualTiltLimit))
            };
        } else {
            result.trajectories = {
                parameters(0.5 * upperAngle.value + 0.5 * topAngle.value),
                parameters(0.5 * lowerAngle.value + 0.5 * topAngle.value),
                parameters(bottomAngle.value)
            };
        }
        return PlanningResult<AutomaticTrajectoryPlan>::success(std::move(result));
    }
}
