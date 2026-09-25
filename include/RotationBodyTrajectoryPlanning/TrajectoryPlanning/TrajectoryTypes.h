#pragma once

#include <RotationBodyTrajectoryPlanning/Core/PlanningTypes.h>

#include <Eigen/Geometry>

#include <cstddef>
#include <optional>
#include <string>
#include <vector>

namespace smrobot::spray::rotationbody
{
    enum class RapidSequenceEntryKind
    {
        SafetyPoint,
        Trajectory
    };

    struct RapidSequenceEntry
    {
        RapidSequenceEntryKind kind{ RapidSequenceEntryKind::SafetyPoint };
        std::string trajectoryPassId;
    };

    struct TrajectoryGenerationParameters
    {
        double sprayDistanceMeters{ 0.1 };
        double tiltRadians{ 0.0 };
        double speedMetersPerSecond{ 0.006 };
        double startExtensionMeters{ 0.0 };
        double endExtensionMeters{ 0.0 };
        std::size_t pointCount{ 0 };
        double positionerRpm{ 0.0 };
        bool reversed{ false };
    };

    struct TrajectoryPosePoint
    {
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        double timeSeconds{ 0.0 };
        Eigen::Isometry3d planningFromTool = Eigen::Isometry3d::Identity();
        bool interpolated{ false };
    };

    struct TrajectoryMetrics
    {
        double pathLengthMeters{ 0.0 };
        double durationSeconds{ 0.0 };
        double minimumPointIntervalSeconds{ 0.0 };
        double maximumPointIntervalSeconds{ 0.0 };
    };

    struct PlannedTrajectory
    {
        TrajectoryGenerationParameters parameters;
        SprayBoundary sourceBoundary;
        Eigen::Vector3d targetSurfaceStart = Eigen::Vector3d::Zero();
        Eigen::Vector3d targetSurfaceEnd = Eigen::Vector3d::Zero();
        std::vector<TrajectoryPosePoint> linearPoints;
        std::vector<TrajectoryPosePoint> relativeHelicalPoints;
        TrajectoryMetrics metrics;

        bool hasValidPoints() const noexcept
        {
            return linearPoints.size() >= 2 &&
                relativeHelicalPoints.size() == linearPoints.size();
        }
    };

    struct TrajectoryPass
    {
        std::string id;
        int order{ 0 };
        bool visible{ true };
        double startOffsetSeconds{ 0.0 };
        double transitionAfterSeconds{ 0.0 };
        PlannedTrajectory trajectory;
    };

    struct TrajectoryGroup
    {
        std::vector<TrajectoryPass> passes;
        std::size_t cycleCount{ 1 };
    };

    enum class TrajectoryDisplayMode
    {
        Linear,
        RelativeHelical
    };

    struct TrajectoryWorkspace
    {
        TrajectoryGenerationParameters parameters;
        std::optional<PlannedTrajectory> currentTrajectory;
        std::string editingPassId;
        TrajectoryGroup group;
        TrajectoryDisplayMode displayMode{ TrajectoryDisplayMode::Linear };
    };

    struct PublishedTrajectoryPlan
    {
        int schemaVersion{ 1 };
        std::string objectId;
        Eigen::Isometry3d baseFromPlanning = Eigen::Isometry3d::Identity();
        // These optional fields were added for consumers that need to
        // evaluate predicted coating statistics per planned tooth region.
        // They remain optional so schema-v1 plans stay usable by legacy
        // trajectory consumers.
        Eigen::Isometry3d planningFromMesh = Eigen::Isometry3d::Identity();
        std::optional<SectionContour> section;
        std::optional<RegionAssignment> regions;
        TrajectoryGroup group;
        Eigen::Vector3d safetyPositionBaseMeters = Eigen::Vector3d::Zero();
        double safetySpeedMetersPerSecond{ 0.2 };
        std::vector<RapidSequenceEntry> executionSequence;
    };
}
