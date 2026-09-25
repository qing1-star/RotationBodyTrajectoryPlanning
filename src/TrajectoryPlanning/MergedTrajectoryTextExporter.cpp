#include <RotationBodyTrajectoryPlanning/TrajectoryPlanning/MergedTrajectoryTextExporter.h>

#include <RotationBodyTrajectoryPlanning/TrajectoryPlanning/TrajectoryGroupEditor.h>

#include <cmath>
#include <iomanip>
#include <locale>
#include <sstream>

namespace smrobot::spray::rotationbody
{
    namespace
    {
        constexpr double kOutputPrecision = 1.0e-9;

        double outputValue(double value)
        {
            return std::abs(value) < kOutputPrecision ? 0.0 : value;
        }

        void appendPose(
            std::ostringstream& stream,
            const Eigen::Isometry3d& baseFromTool,
            double timestampSeconds)
        {
            const Eigen::Matrix4d matrix = baseFromTool.matrix();
            stream << std::fixed << std::setprecision(6);
            for(int row = 0; row < 4; ++row) {
                stream << "  ";
                for(int column = 0; column < 4; ++column) {
                    double value = matrix(row, column);
                    if(column == 3 && row < 3) {
                        value *= 1000.0;
                    }
                    if(column > 0) {
                        stream << ", ";
                    }
                    stream << outputValue(value);
                }
                stream << '\n';
            }
            stream << outputValue(timestampSeconds) << '\n';
        }
    }

    PlanningResult<std::string> MergedTrajectoryTextExporter::format(
        const PublishedTrajectoryPlan& plan)
    {
        if(plan.group.passes.empty() || !plan.baseFromPlanning.matrix().allFinite()) {
            return PlanningResult<std::string>::failure(
                PlanningErrorCode::InvalidArgument,
                "The trajectory group or base-frame transform is invalid.");
        }
        const PlanningResult<void> validation =
            TrajectoryGroupEditor::validate(plan.group);
        if(!validation) {
            return PlanningResult<std::string>::failure(
                validation.error.code,
                validation.error.message);
        }

        std::ostringstream stream;
        stream.imbue(std::locale::classic());
        bool firstPoint = true;
        const double cyclePeriod = TrajectoryGroupEditor::cyclePeriodSeconds(plan.group);
        if(!std::isfinite(cyclePeriod)) {
            return PlanningResult<std::string>::failure(
                PlanningErrorCode::InvalidArgument,
                "The trajectory group has an invalid cycle duration.");
        }
        for(std::size_t cycle = 0; cycle < plan.group.cycleCount; ++cycle) {
          for(const TrajectoryPass& pass : plan.group.passes) {
            for(const TrajectoryPosePoint& point : pass.trajectory.relativeHelicalPoints) {
                if(!firstPoint) {
                    stream << '\n';
                }
                const double timestampSeconds = cycle * cyclePeriod +
                    pass.startOffsetSeconds + point.timeSeconds;
                if(!std::isfinite(timestampSeconds)) {
                    return PlanningResult<std::string>::failure(
                        PlanningErrorCode::InvalidArgument,
                        "The repeated trajectory timestamp is invalid.");
                }
                appendPose(
                    stream,
                    plan.baseFromPlanning * point.planningFromTool,
                    timestampSeconds);
                firstPoint = false;
            }
          }
        }
        if(firstPoint) {
            return PlanningResult<std::string>::failure(
                PlanningErrorCode::InvalidArgument,
                "The trajectory group contains no exportable helical poses.");
        }
        return PlanningResult<std::string>::success(stream.str());
    }
}
