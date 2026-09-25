#include <RotationBodyTrajectoryPlanning/Persistence/PublishedTrajectoryPlanStore.h>
#include <RotationBodyTrajectoryPlanning/Persistence/PublishedTrajectoryPlanContract.h>

#include <nlohmann/json.hpp>

#include <cmath>
#include <stdexcept>

namespace smrobot::spray::rotationbody
{
    namespace
    {
        using Json = nlohmann::json;

        Eigen::Vector2d vector2(const Json& value)
        {
            return { value.at(0).get<double>(), value.at(1).get<double>() };
        }

        Eigen::Vector3d vector3(const Json& value)
        {
            Eigen::Vector3d result(value.at(0).get<double>(), value.at(1).get<double>(), value.at(2).get<double>());
            if(!result.allFinite()) throw std::runtime_error("Vector contains a non-finite value.");
            return result;
        }

        Eigen::Isometry3d transform(const Json& value)
        {
            Eigen::Matrix4d matrix;
            for(int row = 0; row < 4; ++row)
                for(int column = 0; column < 4; ++column)
                    matrix(row, column) = value.at(row * 4 + column).get<double>();
            Eigen::Isometry3d result = Eigen::Isometry3d::Identity();
            result.matrix() = matrix;
            return result;
        }

        RegionLabel regionLabel(const std::string& value)
        {
            if(value == "unclassified") return RegionLabel::Unclassified;
            if(value == "toothTop") return RegionLabel::ToothTop;
            if(value == "toothWall") return RegionLabel::ToothWall;
            if(value == "toothBottom") return RegionLabel::ToothBottom;
            if(value == "transition") return RegionLabel::Transition;
            throw std::runtime_error("Unknown region label.");
        }

        SectionContour section(const Json& value)
        {
            SectionContour result;
            for(const Json& point : value.at("pointsYz")) {
                result.pointsYz.push_back(vector2(point));
            }
            for(const Json& point : value.at("points3d")) {
                result.points3d.push_back(vector3(point));
            }
            result.cumulativeArcLength = value.at("cumulativeArcLength")
                .get<std::vector<double>>();
            result.closed = value.at("closed").get<bool>();
            result.toleranceMeters = value.at("toleranceMeters").get<double>();
            result.diagnostics = value.value("diagnostics", std::vector<std::string>{});
            if(result.pointsYz.size() < 2 || result.points3d.size() != result.pointsYz.size()) {
                throw std::runtime_error("Stored section contour dimensions are invalid.");
            }
            return result;
        }

        RegionAssignment regions(const Json& value)
        {
            RegionAssignment result;
            for(const Json& label : value.at("labels")) {
                result.segmentLabels.push_back(regionLabel(label.get<std::string>()));
            }
            result.segmentConfidence = value.value(
                "confidence", std::vector<double>{});
            result.diagnostics = value.value("diagnostics", std::vector<std::string>{});
            return result;
        }

        TrajectoryPosePoint point(const Json& value)
        {
            TrajectoryPosePoint result;
            result.timeSeconds = value.at("timeSeconds").get<double>();
            result.planningFromTool = transform(value.at("planningFromTool"));
            result.interpolated = value.value("interpolated", false);
            return result;
        }

        TrajectoryGenerationParameters parameters(const Json& value)
        {
            TrajectoryGenerationParameters result;
            result.sprayDistanceMeters = value.at("sprayDistanceMeters").get<double>();
            result.tiltRadians = value.at("tiltRadians").get<double>();
            result.speedMetersPerSecond = value.at("speedMetersPerSecond").get<double>();
            result.startExtensionMeters = value.at("startExtensionMeters").get<double>();
            result.endExtensionMeters = value.at("endExtensionMeters").get<double>();
            result.pointCount = value.at("pointCount").get<std::size_t>();
            result.positionerRpm = value.at("positionerRpm").get<double>();
            result.reversed = value.at("reversed").get<bool>();
            return result;
        }

        PlannedTrajectory trajectory(const Json& value)
        {
            PlannedTrajectory result;
            result.parameters = parameters(value.at("parameters"));
            const Json& boundary = value.at("sourceBoundary");
            result.sourceBoundary.mode = boundary.at("mode").get<std::string>() == "toothTopEnvelope"
                ? BoundaryMode::ToothTopEnvelope : BoundaryMode::MaximumToothTopY;
            for(const Json& item : boundary.at("polygonYz")) result.sourceBoundary.polygonYz.push_back(vector2(item));
            result.sourceBoundary.minimumY = boundary.at("minimumY").get<double>();
            result.sourceBoundary.maximumY = boundary.at("maximumY").get<double>();
            result.sourceBoundary.minimumZ = boundary.at("minimumZ").get<double>();
            result.sourceBoundary.maximumZ = boundary.at("maximumZ").get<double>();
            result.sourceBoundary.outerLineSlopeYPerZ = boundary.at("outerLineSlopeYPerZ").get<double>();
            result.sourceBoundary.outerLineInterceptY = boundary.at("outerLineInterceptY").get<double>();
            result.targetSurfaceStart = vector3(value.at("targetSurfaceStart"));
            result.targetSurfaceEnd = vector3(value.at("targetSurfaceEnd"));
            for(const Json& item : value.at("linearPoints")) result.linearPoints.push_back(point(item));
            for(const Json& item : value.at("relativeHelicalPoints")) result.relativeHelicalPoints.push_back(point(item));
            const Json& metrics = value.at("metrics");
            result.metrics.pathLengthMeters = metrics.at("pathLengthMeters").get<double>();
            result.metrics.durationSeconds = metrics.at("durationSeconds").get<double>();
            result.metrics.minimumPointIntervalSeconds = metrics.at("minimumPointIntervalSeconds").get<double>();
            result.metrics.maximumPointIntervalSeconds = metrics.at("maximumPointIntervalSeconds").get<double>();
            return result;
        }
    }

    std::optional<PublishedTrajectoryPlan> PublishedTrajectoryPlanStore::read(
        const std::string& serializedPayload,
        std::string* errorMessage)
    {
        try {
            const Json root = Json::parse(serializedPayload);
            PublishedTrajectoryPlan result;
            result.schemaVersion = root.value("schemaVersion", 1);
            if(result.schemaVersion < kPublishedTrajectoryPlanMinimumSchemaVersion ||
                result.schemaVersion > kPublishedTrajectoryPlanSchemaVersion) {
                throw std::runtime_error("The published trajectory plan uses an unsupported schema version.");
            }
            result.objectId = root.at("objectId").get<std::string>();
            result.baseFromPlanning = transform(root.at("baseFromPlanning"));
            if(root.contains("planningFromMesh")) {
                result.planningFromMesh = transform(root.at("planningFromMesh"));
            }
            if(root.contains("section")) {
                result.section = section(root.at("section"));
            }
            if(root.contains("regions")) {
                result.regions = regions(root.at("regions"));
                if(!result.section || !result.regions->matches(*result.section)) {
                    throw std::runtime_error("Stored region assignment does not match the section contour.");
                }
            }
            result.group.cycleCount = root.at("group").value("cycleCount", std::size_t{ 1 });
            if(result.group.cycleCount < 1 || result.group.cycleCount > 100) {
                throw std::runtime_error("The saved trajectory group has an invalid cycle count.");
            }
            for(const Json& item : root.at("group").at("passes")) {
                TrajectoryPass pass;
                pass.id = item.at("id").get<std::string>();
                pass.order = item.at("order").get<int>();
                pass.visible = item.at("visible").get<bool>();
                pass.startOffsetSeconds = item.at("startOffsetSeconds").get<double>();
                pass.transitionAfterSeconds = item.at("transitionAfterSeconds").get<double>();
                pass.trajectory = trajectory(item.at("trajectory"));
                result.group.passes.push_back(std::move(pass));
            }
            if(root.contains("safetyPositionBaseMeters")) result.safetyPositionBaseMeters = vector3(root.at("safetyPositionBaseMeters"));
            result.safetySpeedMetersPerSecond = root.value("safetySpeedMetersPerSecond", 0.2);
            if(root.contains("executionSequence")) {
                for(const Json& item : root.at("executionSequence")) {
                    RapidSequenceEntry entry;
                    entry.kind = item.value("kind", "safetyPoint") == "trajectory"
                        ? RapidSequenceEntryKind::Trajectory : RapidSequenceEntryKind::SafetyPoint;
                    entry.trajectoryPassId = item.value("trajectoryPassId", std::string());
                    result.executionSequence.push_back(std::move(entry));
                }
            }
            if(result.executionSequence.empty()) {
                result.executionSequence.push_back({ RapidSequenceEntryKind::SafetyPoint, {} });
                for(const TrajectoryPass& pass : result.group.passes)
                    result.executionSequence.push_back({ RapidSequenceEntryKind::Trajectory, pass.id });
                result.executionSequence.push_back({ RapidSequenceEntryKind::SafetyPoint, {} });
            }
            return result;
        } catch(const std::exception& exception) {
            if(errorMessage != nullptr) *errorMessage = exception.what();
            return std::nullopt;
        }
    }
}
