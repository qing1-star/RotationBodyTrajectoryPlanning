#pragma once

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <cstddef>
#include <string>
#include <utility>
#include <vector>

namespace smrobot::spray::rotationbody
{
    enum class PlanningErrorCode
    {
        None,
        InvalidArgument,
        EmptyMesh,
        NonFiniteGeometry,
        InvalidTriangleIndex,
        DegenerateGeometry,
        InvalidAxisSelection,
        AlignmentFailed,
        SectionFailed,
        NoContour,
        InsufficientRegionData,
        InsufficientToothTopData,
        SingularFit
    };

    struct PlanningError
    {
        PlanningErrorCode code{ PlanningErrorCode::None };
        std::string message;

        explicit operator bool() const noexcept
        {
            return code != PlanningErrorCode::None;
        }
    };

    template<typename T>
    struct PlanningResult
    {
        T value{};
        PlanningError error;
        std::vector<std::string> diagnostics;

        bool ok() const noexcept
        {
            return !static_cast<bool>(error);
        }

        explicit operator bool() const noexcept
        {
            return ok();
        }

        static PlanningResult success(T result)
        {
            PlanningResult output;
            output.value = std::move(result);
            return output;
        }

        static PlanningResult failure(PlanningErrorCode code, std::string message)
        {
            PlanningResult output;
            output.error = { code, std::move(message) };
            return output;
        }
    };

    template<>
    struct PlanningResult<void>
    {
        PlanningError error;
        std::vector<std::string> diagnostics;

        bool ok() const noexcept
        {
            return !static_cast<bool>(error);
        }

        explicit operator bool() const noexcept
        {
            return ok();
        }

        static PlanningResult success()
        {
            return {};
        }

        static PlanningResult failure(PlanningErrorCode code, std::string message)
        {
            PlanningResult output;
            output.error = { code, std::move(message) };
            return output;
        }
    };

    enum class PlanningObjectType
    {
        CompletePart,
        SimulationBlock
    };

    enum class SignedAxis
    {
        PositiveX,
        NegativeX,
        PositiveY,
        NegativeY,
        PositiveZ,
        NegativeZ
    };

    enum class PlanningStage
    {
        NoModel,
        ModelLoaded,
        FrameConfirmed,
        SectionReady,
        RegionsReady,
        SprayBoundaryConfirmed
    };

    enum class RegionLabel
    {
        Unclassified,
        ToothTop,
        ToothWall,
        ToothBottom,
        Transition
    };

    enum class BoundaryMode
    {
        MaximumToothTopY,
        ToothTopEnvelope
    };

    struct TransformComponents
    {
        Eigen::Vector3d translationMeters = Eigen::Vector3d::Zero();
        Eigen::Vector3d rollPitchYawRadians = Eigen::Vector3d::Zero();
    };

    struct ModelStatistics
    {
        std::size_t vertexCount{ 0 };
        std::size_t triangleCount{ 0 };
        Eigen::Vector3d boundsMinimum = Eigen::Vector3d::Zero();
        Eigen::Vector3d boundsMaximum = Eigen::Vector3d::Zero();
        Eigen::Vector3d estimatedAxisInMesh = Eigen::Vector3d::UnitZ();
        double heightMeters{ 0.0 };
        double maximumDiameterMeters{ 0.0 };
        double minimumDiameterMeters{ 0.0 };
        double axisConfidence{ 0.0 };
    };

    struct AlignmentResult
    {
        Eigen::Isometry3d planningFromMesh = Eigen::Isometry3d::Identity();
        Eigen::Isometry3d automaticBaseline = Eigen::Isometry3d::Identity();
        ModelStatistics statistics;
        Eigen::Vector3d bottomAxisCenterInMesh = Eigen::Vector3d::Zero();
        bool lowConfidence{ false };
    };

    struct SectionContour
    {
        std::vector<Eigen::Vector2d> pointsYz;
        std::vector<Eigen::Vector3d> points3d;
        std::vector<double> cumulativeArcLength;
        bool closed{ false };
        double toleranceMeters{ 0.0 };
        std::vector<std::string> diagnostics;

        std::size_t segmentCount() const noexcept
        {
            if (pointsYz.size() < 2)
            {
                return 0;
            }
            return closed ? pointsYz.size() : pointsYz.size() - 1;
        }
    };

    struct RegionAssignment
    {
        std::vector<RegionLabel> segmentLabels;
        std::vector<double> segmentConfidence;
        std::vector<std::string> diagnostics;

        bool matches(const SectionContour& contour) const noexcept
        {
            return segmentLabels.size() == contour.segmentCount();
        }
    };

    struct YzRectangle
    {
        Eigen::Vector2d minimum = Eigen::Vector2d::Zero();
        Eigen::Vector2d maximum = Eigen::Vector2d::Zero();

        YzRectangle normalized() const noexcept
        {
            YzRectangle result;
            result.minimum = minimum.cwiseMin(maximum);
            result.maximum = minimum.cwiseMax(maximum);
            return result;
        }

        bool isFinite() const noexcept
        {
            return minimum.allFinite() && maximum.allFinite();
        }
    };

    struct RegionOverrideCommand
    {
        YzRectangle rectangle;
        RegionLabel label{ RegionLabel::Unclassified };
    };

    struct SprayBoundary
    {
        BoundaryMode mode{ BoundaryMode::MaximumToothTopY };
        std::vector<Eigen::Vector2d> polygonYz;
        double minimumY{ 0.0 };
        double maximumY{ 0.0 };
        double minimumZ{ 0.0 };
        double maximumZ{ 0.0 };
        double outerLineSlopeYPerZ{ 0.0 };
        double outerLineInterceptY{ 0.0 };

        double outerY(double z) const noexcept
        {
            return outerLineSlopeYPerZ * z + outerLineInterceptY;
        }
    };
}
