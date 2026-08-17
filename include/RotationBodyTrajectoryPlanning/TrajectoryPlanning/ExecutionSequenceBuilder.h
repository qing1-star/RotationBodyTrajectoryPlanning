#pragma once

#include <RotationBodyTrajectoryPlanning/TrajectoryPlanning/TrajectoryTypes.h>

#include <Eigen/Geometry>
#include <Eigen/StdVector>

#include <string>
#include <vector>

namespace smrobot::spray::rotationbody
{
    enum class TimedExecutionTargetKind
    {
        InitialPose,
        SafetyPoint,
        Trajectory
    };

    struct TimedExecutionTarget
    {
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        double timeSeconds{ 0.0 };
        Eigen::Isometry3d baseFromTool = Eigen::Isometry3d::Identity();
        TimedExecutionTargetKind kind{ TimedExecutionTargetKind::InitialPose };
        std::string trajectoryPassId;
    };

    using TimedExecutionTargets = std::vector<
        TimedExecutionTarget,
        Eigen::aligned_allocator<TimedExecutionTarget>>;

    class ExecutionSequenceBuilder
    {
    public:
        static PlanningResult<TimedExecutionTargets> build(
            const PublishedTrajectoryPlan& plan,
            const Eigen::Isometry3d& initialBaseFromTool);
    };
}
