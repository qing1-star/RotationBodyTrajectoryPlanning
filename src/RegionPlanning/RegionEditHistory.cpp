#include <RotationBodyTrajectoryPlanning/RegionPlanning/RegionEditHistory.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>

namespace smrobot::spray::rotationbody
{
    namespace
    {
        double cross2d(const Eigen::Vector2d& first, const Eigen::Vector2d& second) noexcept
        {
            return first.x() * second.y() - first.y() * second.x();
        }

        bool pointInside(const Eigen::Vector2d& point, const YzRectangle& rectangle, double tolerance) noexcept
        {
            return point.x() >= rectangle.minimum.x() - tolerance &&
                point.x() <= rectangle.maximum.x() + tolerance &&
                point.y() >= rectangle.minimum.y() - tolerance &&
                point.y() <= rectangle.maximum.y() + tolerance;
        }

        bool segmentsIntersect(
            const Eigen::Vector2d& firstStart,
            const Eigen::Vector2d& firstEnd,
            const Eigen::Vector2d& secondStart,
            const Eigen::Vector2d& secondEnd,
            double coordinateTolerance) noexcept
        {
            const Eigen::Vector2d firstDirection = firstEnd - firstStart;
            const Eigen::Vector2d secondDirection = secondEnd - secondStart;
            const double denominator = cross2d(firstDirection, secondDirection);
            const Eigen::Vector2d delta = secondStart - firstStart;
            const double lengthScale = std::max({
                firstDirection.norm(),
                secondDirection.norm(),
                delta.norm(),
                coordinateTolerance,
                std::numeric_limits<double>::min() });
            const double areaTolerance = std::max(
                coordinateTolerance * lengthScale,
                64.0 * std::numeric_limits<double>::epsilon() * lengthScale * lengthScale);
            if (std::abs(denominator) <= areaTolerance)
            {
                if (std::abs(cross2d(delta, firstDirection)) > areaTolerance)
                {
                    return false;
                }
                const int axis = std::abs(firstDirection.x()) >= std::abs(firstDirection.y()) ? 0 : 1;
                const double firstMinimum = std::min(firstStart[axis], firstEnd[axis]);
                const double firstMaximum = std::max(firstStart[axis], firstEnd[axis]);
                const double secondMinimum = std::min(secondStart[axis], secondEnd[axis]);
                const double secondMaximum = std::max(secondStart[axis], secondEnd[axis]);
                return std::max(firstMinimum, secondMinimum) <=
                    std::min(firstMaximum, secondMaximum) + coordinateTolerance;
            }
            const double firstParameter = cross2d(delta, secondDirection) / denominator;
            const double secondParameter = cross2d(delta, firstDirection) / denominator;
            const double firstParameterTolerance = std::max(
                64.0 * std::numeric_limits<double>::epsilon(),
                coordinateTolerance /
                    std::max(firstDirection.norm(), coordinateTolerance));
            const double secondParameterTolerance = std::max(
                64.0 * std::numeric_limits<double>::epsilon(),
                coordinateTolerance /
                    std::max(secondDirection.norm(), coordinateTolerance));
            return firstParameter >= -firstParameterTolerance &&
                firstParameter <= 1.0 + firstParameterTolerance &&
                secondParameter >= -secondParameterTolerance &&
                secondParameter <= 1.0 + secondParameterTolerance;
        }

        bool segmentIntersectsRectangle(
            const Eigen::Vector2d& start,
            const Eigen::Vector2d& end,
            const YzRectangle& rectangle,
            double tolerance) noexcept
        {
            if (pointInside(start, rectangle, tolerance) || pointInside(end, rectangle, tolerance))
            {
                return true;
            }
            if (std::max(start.x(), end.x()) < rectangle.minimum.x() - tolerance ||
                std::min(start.x(), end.x()) > rectangle.maximum.x() + tolerance ||
                std::max(start.y(), end.y()) < rectangle.minimum.y() - tolerance ||
                std::min(start.y(), end.y()) > rectangle.maximum.y() + tolerance)
            {
                return false;
            }
            const Eigen::Vector2d lowerLeft = rectangle.minimum;
            const Eigen::Vector2d upperRight = rectangle.maximum;
            const Eigen::Vector2d lowerRight(upperRight.x(), lowerLeft.y());
            const Eigen::Vector2d upperLeft(lowerLeft.x(), upperRight.y());
            return segmentsIntersect(start, end, lowerLeft, lowerRight, tolerance) ||
                segmentsIntersect(start, end, lowerRight, upperRight, tolerance) ||
                segmentsIntersect(start, end, upperRight, upperLeft, tolerance) ||
                segmentsIntersect(start, end, upperLeft, lowerLeft, tolerance);
        }

        Eigen::Vector2d segmentEnd(const SectionContour& contour, std::size_t segmentIndex)
        {
            return contour.pointsYz[(segmentIndex + 1) % contour.pointsYz.size()];
        }

        double coordinateTolerance(
            const SectionContour& contour,
            const YzRectangle& rectangle) noexcept
        {
            Eigen::AlignedBox2d bounds;
            for (const Eigen::Vector2d& point : contour.pointsYz)
            {
                bounds.extend(point);
            }
            bounds.extend(rectangle.minimum);
            bounds.extend(rectangle.maximum);
            const double extent = bounds.diagonal().norm();
            const double magnitude = std::max(
                bounds.min().cwiseAbs().maxCoeff(),
                bounds.max().cwiseAbs().maxCoeff());
            const double machineTolerance =
                128.0 * std::numeric_limits<double>::epsilon() *
                std::max({ extent, magnitude, std::numeric_limits<double>::min() });
            return std::max(
                std::max(0.0, contour.toleranceMeters),
                machineTolerance);
        }
    }

    RegionEditHistory::RegionEditHistory(
        SectionContour contour,
        RegionAssignment automaticAssignment)
    {
        resetAutomatic(std::move(contour), std::move(automaticAssignment));
    }

    PlanningResult<void> RegionEditHistory::resetAutomatic(
        SectionContour contour,
        RegionAssignment automaticAssignment)
    {
        if (contour.segmentCount() == 0)
        {
            return PlanningResult<void>::failure(
                PlanningErrorCode::NoContour,
                "Region edit history requires a contour with at least one segment.");
        }
        if (!automaticAssignment.matches(contour))
        {
            return PlanningResult<void>::failure(
                PlanningErrorCode::InvalidArgument,
                "Automatic region labels do not match the contour segment count.");
        }
        if (!automaticAssignment.segmentConfidence.empty() &&
            automaticAssignment.segmentConfidence.size() != contour.segmentCount())
        {
            return PlanningResult<void>::failure(
                PlanningErrorCode::InvalidArgument,
                "Automatic region confidence values do not match the contour segment count.");
        }
        m_contour = std::move(contour);
        m_automatic = std::move(automaticAssignment);
        m_commands.clear();
        m_appliedCommandCount = 0;
        return PlanningResult<void>::success();
    }

    PlanningResult<void> RegionEditHistory::applyRectangle(
        const YzRectangle& rectangle,
        RegionLabel label)
    {
        if (m_contour.segmentCount() == 0 || !m_automatic.matches(m_contour))
        {
            return PlanningResult<void>::failure(
                PlanningErrorCode::InsufficientRegionData,
                "Automatic regions must be initialized before applying an edit.");
        }
        if (!rectangle.isFinite())
        {
            return PlanningResult<void>::failure(
                PlanningErrorCode::NonFiniteGeometry,
                "Region edit rectangle must be finite.");
        }
        const YzRectangle normalized = rectangle.normalized();
        const double tolerance = coordinateTolerance(m_contour, normalized);
        if ((normalized.maximum - normalized.minimum).minCoeff() <= tolerance)
        {
            return PlanningResult<void>::failure(
                PlanningErrorCode::InvalidArgument,
                "Region edit rectangle must have a non-zero area.");
        }

        bool intersects = false;
        for (std::size_t segment = 0; segment < m_contour.segmentCount(); ++segment)
        {
            if (segmentIntersectsRectangle(
                    m_contour.pointsYz[segment],
                    segmentEnd(m_contour, segment),
                    normalized,
                    tolerance))
            {
                intersects = true;
                break;
            }
        }
        if (!intersects)
        {
            return PlanningResult<void>::failure(
                PlanningErrorCode::InvalidArgument,
                "Region edit rectangle does not intersect the section contour.");
        }

        m_commands.resize(m_appliedCommandCount);
        m_commands.push_back({ normalized, label });
        m_appliedCommandCount = m_commands.size();
        return PlanningResult<void>::success();
    }

    bool RegionEditHistory::canUndo() const noexcept
    {
        return m_appliedCommandCount > 0;
    }

    bool RegionEditHistory::canRedo() const noexcept
    {
        return m_appliedCommandCount < m_commands.size();
    }

    bool RegionEditHistory::undo() noexcept
    {
        if (!canUndo())
        {
            return false;
        }
        --m_appliedCommandCount;
        return true;
    }

    bool RegionEditHistory::redo() noexcept
    {
        if (!canRedo())
        {
            return false;
        }
        ++m_appliedCommandCount;
        return true;
    }

    void RegionEditHistory::restoreAutomatic() noexcept
    {
        m_commands.clear();
        m_appliedCommandCount = 0;
    }

    RegionAssignment RegionEditHistory::resolved() const
    {
        RegionAssignment result = m_automatic;
        if (!result.matches(m_contour))
        {
            return {};
        }
        for (std::size_t commandIndex = 0; commandIndex < m_appliedCommandCount; ++commandIndex)
        {
            const RegionOverrideCommand& command = m_commands[commandIndex];
            const double tolerance = coordinateTolerance(m_contour, command.rectangle);
            for (std::size_t segment = 0; segment < m_contour.segmentCount(); ++segment)
            {
                if (segmentIntersectsRectangle(
                        m_contour.pointsYz[segment],
                        segmentEnd(m_contour, segment),
                        command.rectangle,
                        tolerance))
                {
                    result.segmentLabels[segment] = command.label;
                    if (result.segmentConfidence.size() == result.segmentLabels.size())
                    {
                        result.segmentConfidence[segment] = 1.0;
                    }
                }
            }
        }
        return result;
    }

    const RegionAssignment& RegionEditHistory::automatic() const noexcept
    {
        return m_automatic;
    }

    const std::vector<RegionOverrideCommand>& RegionEditHistory::commands() const noexcept
    {
        return m_commands;
    }

    std::size_t RegionEditHistory::appliedCommandCount() const noexcept
    {
        return m_appliedCommandCount;
    }
}
