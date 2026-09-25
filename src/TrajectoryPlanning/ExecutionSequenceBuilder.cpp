#include <RotationBodyTrajectoryPlanning/TrajectoryPlanning/ExecutionSequenceBuilder.h>

#include <RotationBodyTrajectoryPlanning/TrajectoryPlanning/TrajectoryGroupEditor.h>
#include <RotationBodyTrajectoryPlanning/TrajectoryPlanning/TrajectoryPlanner.h>

#include <Eigen/Geometry>

#include <algorithm>
#include <cmath>
#include <optional>

namespace smrobot::spray::rotationbody
{
    namespace
    {
        constexpr double kSafetyAngularSpeedRadiansPerSecond =
            100.0 * 3.14159265358979323846 / 180.0;
        constexpr double kPoseTolerance = 1.0e-9;
        constexpr double kReturnPathTolerance = 1.0e-7;

        const TrajectoryPass* sequencePass(
            const PublishedTrajectoryPlan& plan,
            const RapidSequenceEntry& entry)
        {
            return entry.kind == RapidSequenceEntryKind::Trajectory
                ? TrajectoryGroupEditor::find(plan.group, entry.trajectoryPassId)
                : nullptr;
        }

        double transferDuration(
            const Eigen::Isometry3d& from,
            const Eigen::Isometry3d& to,
            double linearSpeedMetersPerSecond)
        {
            const double linearSeconds =
                (to.translation() - from.translation()).norm() /
                linearSpeedMetersPerSecond;
            const Eigen::AngleAxisd rotation(from.linear().transpose() * to.linear());
            const double angularSeconds =
                std::abs(rotation.angle()) / kSafetyAngularSpeedRadiansPerSecond;
            return std::max(linearSeconds, angularSeconds);
        }

        bool samePose(
            const Eigen::Isometry3d& lhs,
            const Eigen::Isometry3d& rhs)
        {
            return lhs.translation().isApprox(rhs.translation(), kPoseTolerance) &&
                lhs.linear().isApprox(rhs.linear(), kPoseTolerance);
        }

        bool isCoincidentReversePass(
            const TrajectoryPass& previous,
            const TrajectoryPass& candidate)
        {
            if(!previous.trajectory.hasValidPoints() ||
                !candidate.trajectory.hasValidPoints()) {
                return false;
            }
            const TrajectoryPosePoint& previousStart =
                previous.trajectory.linearPoints.front();
            const TrajectoryPosePoint& previousEnd =
                previous.trajectory.linearPoints.back();
            const TrajectoryPosePoint& candidateStart =
                candidate.trajectory.linearPoints.front();
            const TrajectoryPosePoint& candidateEnd =
                candidate.trajectory.linearPoints.back();
            const bool endpointsReversed =
                (previousEnd.planningFromTool.translation() -
                    candidateStart.planningFromTool.translation()).norm() <=
                    kReturnPathTolerance &&
                (previousStart.planningFromTool.translation() -
                    candidateEnd.planningFromTool.translation()).norm() <=
                    kReturnPathTolerance;
            const bool sprayAxesMatch =
                previousEnd.planningFromTool.linear().col(2).isApprox(
                    candidateStart.planningFromTool.linear().col(2),
                    kReturnPathTolerance) &&
                previousStart.planningFromTool.linear().col(2).isApprox(
                    candidateEnd.planningFromTool.linear().col(2),
                    kReturnPathTolerance);
            return endpointsReversed && sprayAxesMatch;
        }

        Eigen::Matrix3d orientationWithoutTilt(
            const PublishedTrajectoryPlan& plan,
            const TrajectoryPosePoint& point)
        {
            return TrajectoryPlanner::levelSprayAxisAroundLocalY(
                (plan.baseFromPlanning * point.planningFromTool).linear());
        }

        PlanningResult<Eigen::Matrix3d> safetyOrientation(
            const PublishedTrajectoryPlan& plan,
            std::size_t safetyIndex)
        {
            for(std::size_t index = safetyIndex + 1;
                index < plan.executionSequence.size();
                ++index) {
                if(const TrajectoryPass* pass =
                    sequencePass(plan, plan.executionSequence[index])) {
                    return PlanningResult<Eigen::Matrix3d>::success(
                        orientationWithoutTilt(
                            plan, pass->trajectory.linearPoints.front()));
                }
            }
            for(std::size_t index = safetyIndex; index-- > 0;) {
                if(const TrajectoryPass* pass =
                    sequencePass(plan, plan.executionSequence[index])) {
                    return PlanningResult<Eigen::Matrix3d>::success(
                        orientationWithoutTilt(
                            plan, pass->trajectory.linearPoints.back()));
                }
            }
            return PlanningResult<Eigen::Matrix3d>::failure(
                PlanningErrorCode::InvalidArgument,
                "A safety point requires an adjacent trajectory from which to remove tilt.");
        }

    }

    PlanningResult<TimedExecutionTargets> ExecutionSequenceBuilder::build(
        const PublishedTrajectoryPlan& plan,
        const Eigen::Isometry3d& initialBaseFromTool)
    {
        if(plan.executionSequence.empty() ||
            !initialBaseFromTool.matrix().allFinite() ||
            !plan.baseFromPlanning.matrix().allFinite() ||
            !plan.safetyPositionBaseMeters.allFinite() ||
            !std::isfinite(plan.safetySpeedMetersPerSecond) ||
            plan.safetySpeedMetersPerSecond <= 0.0) {
            return PlanningResult<TimedExecutionTargets>::failure(
                PlanningErrorCode::InvalidArgument,
                "The execution sequence or safety motion settings are invalid.");
        }
        const PlanningResult<void> validation =
            TrajectoryGroupEditor::validate(plan.group);
        if(!validation) {
            return PlanningResult<TimedExecutionTargets>::failure(
                validation.error.code,
                validation.error.message);
        }

        TimedExecutionTargets targets;
        targets.reserve(plan.executionSequence.size() * 2 + 1);
        TimedExecutionTarget initial;
        initial.baseFromTool = initialBaseFromTool;
        targets.push_back(initial);

        double timeSeconds = 0.0;
        Eigen::Isometry3d current = initialBaseFromTool;
        const TrajectoryPass* previousTrajectoryPass = nullptr;
        const auto appendTransfer = [&](
            const Eigen::Isometry3d& pose,
            TimedExecutionTargetKind kind,
            const std::string& passId,
            TimedExecutionTargets& output,
            double& time,
            Eigen::Isometry3d& previous) {
            if(samePose(previous, pose)) {
                return;
            }
            time += transferDuration(
                previous, pose, plan.safetySpeedMetersPerSecond);
            TimedExecutionTarget target;
            target.timeSeconds = time;
            target.baseFromTool = pose;
            target.kind = kind;
            target.trajectoryPassId = passId;
            output.push_back(std::move(target));
            previous = pose;
        };

        for(std::size_t sequenceIndex = 0;
            sequenceIndex < plan.executionSequence.size();
            ++sequenceIndex) {
            const RapidSequenceEntry& entry = plan.executionSequence[sequenceIndex];
            if(entry.kind == RapidSequenceEntryKind::SafetyPoint) {
                const PlanningResult<Eigen::Matrix3d> orientation =
                    safetyOrientation(plan, sequenceIndex);
                if(!orientation) {
                    return PlanningResult<TimedExecutionTargets>::failure(
                        orientation.error.code,
                        orientation.error.message);
                }
                Eigen::Isometry3d safetyPose = Eigen::Isometry3d::Identity();
                safetyPose.translation() = plan.safetyPositionBaseMeters;
                safetyPose.linear() = orientation.value;
                appendTransfer(
                    safetyPose,
                    TimedExecutionTargetKind::SafetyPoint,
                    {},
                    targets,
                    timeSeconds,
                    current);
                previousTrajectoryPass = nullptr;
                continue;
            }

            const TrajectoryPass* pass = sequencePass(plan, entry);
            if(pass == nullptr || !pass->trajectory.hasValidPoints()) {
                return PlanningResult<TimedExecutionTargets>::failure(
                    PlanningErrorCode::InvalidArgument,
                    "The execution sequence references a missing trajectory.");
            }
            const Eigen::Isometry3d start = plan.baseFromPlanning *
                pass->trajectory.linearPoints.front().planningFromTool;
            const std::optional<Eigen::Matrix3d> returnOrientation =
                previousTrajectoryPass != nullptr &&
                    isCoincidentReversePass(*previousTrajectoryPass, *pass)
                ? std::optional<Eigen::Matrix3d>(current.linear())
                : std::nullopt;
            Eigen::Isometry3d continuousStart = start;
            if(returnOrientation) {
                continuousStart.linear() = *returnOrientation;
            }
            appendTransfer(
                continuousStart,
                TimedExecutionTargetKind::Trajectory,
                entry.trajectoryPassId,
                targets,
                timeSeconds,
                current);
            const double passStartTime = timeSeconds;
            for(std::size_t pointIndex = 1;
                pointIndex < pass->trajectory.linearPoints.size();
                ++pointIndex) {
                const TrajectoryPosePoint& point =
                    pass->trajectory.linearPoints[pointIndex];
                TimedExecutionTarget target;
                target.timeSeconds = passStartTime + point.timeSeconds;
                target.baseFromTool = plan.baseFromPlanning * point.planningFromTool;
                if(returnOrientation) {
                    target.baseFromTool.linear() = *returnOrientation;
                }
                target.kind = TimedExecutionTargetKind::Trajectory;
                target.trajectoryPassId = entry.trajectoryPassId;
                targets.push_back(std::move(target));
            }
            current = plan.baseFromPlanning *
                pass->trajectory.linearPoints.back().planningFromTool;
            if(returnOrientation) {
                current.linear() = *returnOrientation;
            }
            timeSeconds = passStartTime + pass->trajectory.metrics.durationSeconds;
            previousTrajectoryPass = pass;

            if(pass->transitionAfterSeconds > 0.0) {
                timeSeconds += pass->transitionAfterSeconds;
                TimedExecutionTarget hold;
                hold.timeSeconds = timeSeconds;
                hold.baseFromTool = current;
                hold.kind = TimedExecutionTargetKind::Trajectory;
                hold.trajectoryPassId = entry.trajectoryPassId;
                targets.push_back(std::move(hold));
            }
        }

        return PlanningResult<TimedExecutionTargets>::success(std::move(targets));
    }
}
