#include <RotationBodyTrajectoryPlanning/Persistence/PublishedTrajectoryPlanSprayTrajectoryAdapter.h>

#include <RotationBodyTrajectoryPlanning/TrajectoryPlanning/TrajectoryGroupEditor.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace smrobot::spray::rotationbody
{
    namespace
    {
        constexpr const char* kCoatingProcessId = "paper_gaussian";

        bool finiteNonNegative(double value)
        {
            return std::isfinite(value) && value >= 0.0;
        }

        spraytrajectory::SprayPathPoint makePoint(
            const Eigen::Isometry3d& baseFromTool,
            double timestampSeconds,
            double sprayDistanceMeters,
            int passIndex,
            bool sprayEnabled)
        {
            spraytrajectory::SprayPathPoint point;
            point.time = timestampSeconds;
            point.tcpPose = baseFromTool;
            point.sprayEnabled = sprayEnabled;
            point.processId = kCoatingProcessId;
            point.targetDistance = sprayDistanceMeters;
            point.targetNormal = baseFromTool.linear().col(2).normalized();
            point.workpieceRegionId = passIndex;
            return point;
        }

        void appendTransitionGuards(
            spraytrajectory::SpraySegment& segment,
            const Eigen::Isometry3d& finalPose,
            double finalTimestampSeconds,
            double nextPassStartSeconds)
        {
            const double gapSeconds = nextPassStartSeconds - finalTimestampSeconds;
            if(!(gapSeconds > 0.0) || !std::isfinite(gapSeconds)) {
                return;
            }

            // The sampling API flattens segments. Keep the gun off through a
            // scheduling gap so the gap cannot become an artificial deposit.
            const double epsilon = std::min(1.0e-6, gapSeconds * 0.25);
            segment.points.push_back(makePoint(
                finalPose,
                finalTimestampSeconds + epsilon,
                0.0,
                segment.passIndex,
                false));
            segment.points.push_back(makePoint(
                finalPose,
                nextPassStartSeconds - epsilon,
                0.0,
                segment.passIndex,
                false));
        }
    }

    std::optional<spraytrajectory::SprayTrajectory>
    PublishedTrajectoryPlanSprayTrajectoryAdapter::convert(
        const PublishedTrajectoryPlan& plan,
        std::string* errorMessage)
    {
        const auto fail = [errorMessage](const std::string& message)
            -> std::optional<spraytrajectory::SprayTrajectory> {
            if(errorMessage != nullptr) {
                *errorMessage = message;
            }
            return std::nullopt;
        };

        if(plan.objectId.empty() || !plan.baseFromPlanning.matrix().allFinite()) {
            return fail("The saved trajectory plan has an invalid workpiece or base transform.");
        }
        const PlanningResult<void> validation = TrajectoryGroupEditor::validate(plan.group);
        if(!validation) {
            return fail(validation.error.message);
        }

        spraytrajectory::SprayTrajectory result;
        result.name = "Saved rotation-body helical trajectory";
        const double cyclePeriod = TrajectoryGroupEditor::cyclePeriodSeconds(plan.group);
        if(!finiteNonNegative(cyclePeriod)) {
            return fail("The saved trajectory group has an invalid cycle duration.");
        }
        result.segments.reserve(plan.group.passes.size() * plan.group.cycleCount);
        for(std::size_t cycle = 0; cycle < plan.group.cycleCount; ++cycle) {
          for(std::size_t passIndex = 0; passIndex < plan.group.passes.size(); ++passIndex) {
            const TrajectoryPass& pass = plan.group.passes[passIndex];
            if(pass.trajectory.relativeHelicalPoints.size() < 2 ||
                !finiteNonNegative(pass.startOffsetSeconds) ||
                !finiteNonNegative(pass.trajectory.parameters.sprayDistanceMeters)) {
                return fail("The saved trajectory plan contains an invalid helical pass.");
            }

            spraytrajectory::SpraySegment segment;
            segment.processId = kCoatingProcessId;
            segment.sprayEnabled = true;
            segment.passIndex = pass.order;
            segment.points.reserve(pass.trajectory.relativeHelicalPoints.size() + 2);
            for(const TrajectoryPosePoint& helicalPoint : pass.trajectory.relativeHelicalPoints) {
                const double timestampSeconds = cycle * cyclePeriod +
                    pass.startOffsetSeconds + helicalPoint.timeSeconds;
                const Eigen::Isometry3d baseFromTool =
                    plan.baseFromPlanning * helicalPoint.planningFromTool;
                if(!finiteNonNegative(timestampSeconds) ||
                    !baseFromTool.matrix().allFinite()) {
                    return fail("The saved helical trajectory contains an invalid matrix or timestamp.");
                }
                segment.points.push_back(makePoint(
                    baseFromTool,
                    timestampSeconds,
                    pass.trajectory.parameters.sprayDistanceMeters,
                    pass.order,
                    true));
            }

            if(passIndex + 1 < plan.group.passes.size() ||
                cycle + 1 < plan.group.cycleCount) {
                const bool nextCycle = passIndex + 1 == plan.group.passes.size();
                const TrajectoryPass& nextPass = plan.group.passes[
                    nextCycle ? 0 : passIndex + 1];
                appendTransitionGuards(
                    segment,
                    segment.points.back().tcpPose,
                    segment.points.back().time,
                    (cycle + (nextCycle ? 1 : 0)) * cyclePeriod +
                        nextPass.startOffsetSeconds);
            }
            result.segments.push_back(std::move(segment));
          }
        }
        if(result.empty()) {
            return fail("The saved trajectory plan contains no helical spray points.");
        }
        if(errorMessage != nullptr) {
            errorMessage->clear();
        }
        return result;
    }
}
