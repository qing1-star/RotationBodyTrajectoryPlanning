#include <RotationBodyTrajectoryPlanning/TrajectoryPlanning/RelativeHelicalTrajectoryBuilder.h>

#include <Eigen/Geometry>

#include <cmath>

namespace smrobot::spray::rotationbody
{
    PlanningResult<std::vector<TrajectoryPosePoint>>
    RelativeHelicalTrajectoryBuilder::build(
        const std::vector<TrajectoryPosePoint>& linearPoints,
        double positionerRpm,
        double startOffsetSeconds)
    {
        if(linearPoints.empty() || !std::isfinite(positionerRpm) ||
            !std::isfinite(startOffsetSeconds) || startOffsetSeconds < 0.0) {
            return PlanningResult<std::vector<TrajectoryPosePoint>>::failure(
                PlanningErrorCode::InvalidArgument,
                "Relative helical trajectory input is invalid.");
        }
        constexpr double pi = 3.14159265358979323846;
        const double radiansPerSecond = positionerRpm * 2.0 * pi / 60.0;
        std::vector<TrajectoryPosePoint> result;
        result.reserve(linearPoints.size());
        double previousTime = -1.0;
        for(const TrajectoryPosePoint& source : linearPoints) {
            if(!std::isfinite(source.timeSeconds) || source.timeSeconds < 0.0 ||
                source.timeSeconds < previousTime ||
                !source.planningFromTool.matrix().allFinite()) {
                return PlanningResult<std::vector<TrajectoryPosePoint>>::failure(
                    PlanningErrorCode::InvalidArgument,
                    "Linear trajectory contains invalid or decreasing timestamps.");
            }
            const double angle =
                radiansPerSecond * (startOffsetSeconds + source.timeSeconds);
            Eigen::Isometry3d rotation = Eigen::Isometry3d::Identity();
            rotation.linear() = Eigen::AngleAxisd(angle, Eigen::Vector3d::UnitZ()).toRotationMatrix();
            TrajectoryPosePoint derived = source;
            derived.planningFromTool = rotation * source.planningFromTool;
            result.push_back(std::move(derived));
            previousTime = source.timeSeconds;
        }
        return PlanningResult<std::vector<TrajectoryPosePoint>>::success(std::move(result));
    }
}
