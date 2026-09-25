#include <RotationBodyTrajectoryPlanning/TrajectoryPlanning/TrajectoryGroupEditor.h>

#include <RotationBodyTrajectoryPlanning/TrajectoryPlanning/TrajectoryEditor.h>

#include <algorithm>
#include <cmath>
#include <iterator>
#include <set>

namespace smrobot::spray::rotationbody
{
    namespace
    {
        TrajectoryPass* findMutable(TrajectoryGroup& group, const std::string& passId)
        {
            const auto found = std::find_if(
                group.passes.begin(),
                group.passes.end(),
                [&](const TrajectoryPass& pass) { return pass.id == passId; });
            return found == group.passes.end() ? nullptr : &*found;
        }

        std::string nextId(const TrajectoryGroup& group)
        {
            int suffix = 1;
            for(;; ++suffix) {
                const std::string candidate = "trajectory-" + std::to_string(suffix);
                if(TrajectoryGroupEditor::find(group, candidate) == nullptr) {
                    return candidate;
                }
            }
        }
    }

    PlanningResult<std::string> TrajectoryGroupEditor::addOrUpdate(
        TrajectoryGroup& group,
        const PlannedTrajectory& trajectory,
        const std::string& editingPassId)
    {
        if(!trajectory.hasValidPoints()) {
            return PlanningResult<std::string>::failure(
                PlanningErrorCode::InvalidArgument,
                "Generate a valid linear and relative helical trajectory first.");
        }
        std::string id = editingPassId;
        if(!id.empty()) {
            TrajectoryPass* existing = findMutable(group, id);
            if(existing == nullptr) {
                return PlanningResult<std::string>::failure(
                    PlanningErrorCode::InvalidArgument,
                    "The trajectory selected for update no longer exists.");
            }
            existing->trajectory = trajectory;
        } else {
            id = nextId(group);
            TrajectoryPass pass;
            pass.id = id;
            pass.order = static_cast<int>(group.passes.size()) + 1;
            pass.trajectory = trajectory;
            group.passes.push_back(std::move(pass));
        }
        PlanningResult<void> refreshed = refreshSchedule(group);
        if(!refreshed) {
            return PlanningResult<std::string>::failure(
                refreshed.error.code,
                refreshed.error.message);
        }
        return PlanningResult<std::string>::success(std::move(id));
    }

    PlanningResult<void> TrajectoryGroupEditor::remove(
        TrajectoryGroup& group,
        const std::string& passId)
    {
        const auto found = std::find_if(
            group.passes.begin(),
            group.passes.end(),
            [&](const TrajectoryPass& pass) { return pass.id == passId; });
        if(found == group.passes.end()) {
            return PlanningResult<void>::failure(
                PlanningErrorCode::InvalidArgument,
                "The trajectory to remove does not exist.");
        }
        group.passes.erase(found);
        for(std::size_t index = 0; index < group.passes.size(); ++index) {
            group.passes[index].order = static_cast<int>(index) + 1;
        }
        return refreshSchedule(group);
    }

    PlanningResult<void> TrajectoryGroupEditor::setVisible(
        TrajectoryGroup& group,
        const std::string& passId,
        bool visible)
    {
        TrajectoryPass* pass = findMutable(group, passId);
        if(pass == nullptr) {
            return PlanningResult<void>::failure(
                PlanningErrorCode::InvalidArgument,
                "The trajectory visibility target does not exist.");
        }
        pass->visible = visible;
        return PlanningResult<void>::success();
    }

    PlanningResult<void> TrajectoryGroupEditor::setTransitionAfter(
        TrajectoryGroup& group,
        const std::string& passId,
        double seconds,
        double zeroIntervalPositionToleranceMeters)
    {
        if(!std::isfinite(seconds) || seconds < 0.0 ||
            !std::isfinite(zeroIntervalPositionToleranceMeters) ||
            zeroIntervalPositionToleranceMeters < 0.0) {
            return PlanningResult<void>::failure(
                PlanningErrorCode::InvalidArgument,
                "Trajectory transition time must be finite and non-negative.");
        }
        const auto found = std::find_if(
            group.passes.begin(),
            group.passes.end(),
            [&](const TrajectoryPass& pass) { return pass.id == passId; });
        if(found == group.passes.end()) {
            return PlanningResult<void>::failure(
                PlanningErrorCode::InvalidArgument,
                "The trajectory transition target does not exist.");
        }

        const std::size_t passIndex = static_cast<std::size_t>(
            std::distance(group.passes.begin(), found));
        if(seconds == 0.0 && passIndex + 1 < group.passes.size()) {
            const TrajectoryPass& nextPass = group.passes[passIndex + 1];
            if(!found->trajectory.hasValidPoints() || !nextPass.trajectory.hasValidPoints()) {
                return PlanningResult<void>::failure(
                    PlanningErrorCode::InvalidArgument,
                    "A zero transition requires valid adjacent trajectories.");
            }
            const Eigen::Vector3d previousEnd =
                found->trajectory.linearPoints.back().planningFromTool.translation();
            const Eigen::Vector3d nextStart =
                nextPass.trajectory.linearPoints.front().planningFromTool.translation();
            if((previousEnd - nextStart).norm() > zeroIntervalPositionToleranceMeters) {
                return PlanningResult<void>::failure(
                    PlanningErrorCode::InvalidArgument,
                    "A zero transition requires coincident positions between adjacent passes.");
            }
        }
        found->transitionAfterSeconds = seconds;
        return refreshSchedule(group);
    }

    const TrajectoryPass* TrajectoryGroupEditor::find(
        const TrajectoryGroup& group,
        const std::string& passId) noexcept
    {
        const auto found = std::find_if(
            group.passes.begin(),
            group.passes.end(),
            [&](const TrajectoryPass& pass) { return pass.id == passId; });
        return found == group.passes.end() ? nullptr : &*found;
    }

    PlanningResult<void> TrajectoryGroupEditor::refreshSchedule(TrajectoryGroup& group)
    {
        double offset = 0.0;
        for(std::size_t index = 0; index < group.passes.size(); ++index) {
            TrajectoryPass& pass = group.passes[index];
            pass.order = static_cast<int>(index) + 1;
            pass.startOffsetSeconds = offset;
            PlanningResult<void> rebuilt =
                TrajectoryEditor::rebuildDerived(pass.trajectory, offset);
            if(!rebuilt) {
                return rebuilt;
            }
            offset += pass.trajectory.metrics.durationSeconds +
                pass.transitionAfterSeconds;
        }
        return PlanningResult<void>::success();
    }

    double TrajectoryGroupEditor::cyclePeriodSeconds(const TrajectoryGroup& group) noexcept
    {
        if(group.passes.empty()) return 0.0;
        const TrajectoryPass& last = group.passes.back();
        return last.startOffsetSeconds + last.trajectory.metrics.durationSeconds +
            std::max(last.transitionAfterSeconds, 1.0e-3);
    }

    PlanningResult<void> TrajectoryGroupEditor::validate(
        const TrajectoryGroup& group,
        double zeroIntervalPositionToleranceMeters)
    {
        if(group.cycleCount < 1 || group.cycleCount > 100) {
            return PlanningResult<void>::failure(
                PlanningErrorCode::InvalidArgument,
                "Trajectory group cycle count must be between 1 and 100.");
        }
        std::set<std::string> ids;
        for(std::size_t index = 0; index < group.passes.size(); ++index) {
            const TrajectoryPass& pass = group.passes[index];
            if(pass.id.empty() || !ids.insert(pass.id).second ||
                pass.order != static_cast<int>(index) + 1 ||
                !std::isfinite(pass.startOffsetSeconds) || pass.startOffsetSeconds < 0.0 ||
                !std::isfinite(pass.transitionAfterSeconds) ||
                pass.transitionAfterSeconds < 0.0 || !pass.trajectory.hasValidPoints()) {
                return PlanningResult<void>::failure(
                    PlanningErrorCode::InvalidArgument,
                    "The trajectory group contains invalid pass data.");
            }
            for(std::size_t pointIndex = 0;
                pointIndex < pass.trajectory.linearPoints.size();
                ++pointIndex) {
                const TrajectoryPosePoint& linear =
                    pass.trajectory.linearPoints[pointIndex];
                const TrajectoryPosePoint& helical =
                    pass.trajectory.relativeHelicalPoints[pointIndex];
                if(!std::isfinite(linear.timeSeconds) ||
                    !linear.planningFromTool.matrix().allFinite() ||
                    !helical.planningFromTool.matrix().allFinite() ||
                    std::abs(linear.timeSeconds - helical.timeSeconds) > 1.0e-12) {
                    return PlanningResult<void>::failure(
                        PlanningErrorCode::InvalidArgument,
                        "Linear and relative helical trajectory points are not aligned.");
                }
            }
            if(index + 1 < group.passes.size() && pass.transitionAfterSeconds == 0.0) {
                const Eigen::Vector3d previousEnd =
                    pass.trajectory.linearPoints.back().planningFromTool.translation();
                const Eigen::Vector3d nextStart =
                    group.passes[index + 1].trajectory.linearPoints.front()
                        .planningFromTool.translation();
                if((previousEnd - nextStart).norm() >
                    zeroIntervalPositionToleranceMeters) {
                    return PlanningResult<void>::failure(
                        PlanningErrorCode::InvalidArgument,
                        "A zero transition requires coincident positions between adjacent passes.");
                }
            }
        }
        return PlanningResult<void>::success();
    }
}
