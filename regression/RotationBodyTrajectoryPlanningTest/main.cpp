#include <RotationBodyTrajectoryPlanning/Core/PlanningTypes.h>
#include <RotationBodyTrajectoryPlanning/Core/TransformUtils.h>
#include <RotationBodyTrajectoryPlanning/Core/TriangleMesh.h>
#include <RotationBodyTrajectoryPlanning/Alignment/ModelTransformOperations.h>
#include <RotationBodyTrajectoryPlanning/Alignment/RotationBodyAlignmentSolver.h>
#include <RotationBodyTrajectoryPlanning/Alignment/SimulationBlockPlacementSolver.h>
#include <RotationBodyTrajectoryPlanning/RegionPlanning/RegionEditHistory.h>
#include <RotationBodyTrajectoryPlanning/RegionPlanning/SprayBoundaryBuilder.h>
#include <RotationBodyTrajectoryPlanning/RegionPlanning/ToothRegionRecognizer.h>
#include <RotationBodyTrajectoryPlanning/Sectioning/ContourTopology.h>
#include <RotationBodyTrajectoryPlanning/Sectioning/YzSectionExtractor.h>
#include <RotationBodyTrajectoryPlanning/TrajectoryPlanning/TrajectoryEditor.h>
#include <RotationBodyTrajectoryPlanning/TrajectoryPlanning/ExecutionSequenceBuilder.h>
#include <RotationBodyTrajectoryPlanning/TrajectoryPlanning/TrajectoryGroupEditor.h>
#include <RotationBodyTrajectoryPlanning/TrajectoryPlanning/TrajectoryPlanner.h>
#include <RotationBodyTrajectoryPlanning/ABBTranslation/RapidModuleGenerator.h>
#include <RotationBodyTrajectoryPlanning/Calibration/WorkpieceCalibration.h>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace rotationbody = smrobot::spray::rotationbody;

namespace
{
    int failureCount = 0;

    void expect(bool condition, const std::string& message)
    {
        if (!condition)
        {
            ++failureCount;
            std::cerr << "FAILED: " << message << '\n';
        }
    }

    bool near(double first, double second, double tolerance = 1.0e-9)
    {
        return std::abs(first - second) <= tolerance;
    }

    bool nearVector(
        const Eigen::Vector3d& first,
        const Eigen::Vector3d& second,
        double tolerance = 1.0e-9)
    {
        return (first - second).norm() <= tolerance;
    }

    rotationbody::TriangleMesh makeTetrahedron()
    {
        rotationbody::TriangleMesh mesh;
        mesh.positions = {
            { 0.0, 0.0, 0.0 },
            { 1.0, 0.0, 0.0 },
            { 0.0, 1.0, 0.0 },
            { 0.0, 0.0, 1.0 }
        };
        mesh.triangles = {
            { 0, 2, 1 },
            { 0, 1, 3 },
            { 1, 2, 3 },
            { 2, 0, 3 }
        };
        return mesh;
    }

    rotationbody::TriangleMesh makeBox(
        const Eigen::Vector3d& minimum,
        const Eigen::Vector3d& maximum)
    {
        rotationbody::TriangleMesh mesh;
        mesh.positions = {
            { minimum.x(), minimum.y(), minimum.z() },
            { maximum.x(), minimum.y(), minimum.z() },
            { maximum.x(), maximum.y(), minimum.z() },
            { minimum.x(), maximum.y(), minimum.z() },
            { minimum.x(), minimum.y(), maximum.z() },
            { maximum.x(), minimum.y(), maximum.z() },
            { maximum.x(), maximum.y(), maximum.z() },
            { minimum.x(), maximum.y(), maximum.z() }
        };
        mesh.triangles = {
            { 0, 2, 1 }, { 0, 3, 2 },
            { 4, 5, 6 }, { 4, 6, 7 },
            { 0, 1, 5 }, { 0, 5, 4 },
            { 1, 2, 6 }, { 1, 6, 5 },
            { 2, 3, 7 }, { 2, 7, 6 },
            { 3, 0, 4 }, { 3, 4, 7 }
        };
        return mesh;
    }

    rotationbody::TriangleMesh makeRevolvedProfile(
        const std::vector<double>& axialStations,
        const std::vector<double>& radii,
        std::size_t sectorCount = 48)
    {
        rotationbody::TriangleMesh mesh;
        if (axialStations.size() != radii.size() || axialStations.size() < 2 || sectorCount < 3)
        {
            return mesh;
        }
        for (std::size_t station = 0; station < axialStations.size(); ++station)
        {
            for (std::size_t sector = 0; sector < sectorCount; ++sector)
            {
                const double angle = 2.0 * rotationbody::pi *
                    static_cast<double>(sector) / static_cast<double>(sectorCount);
                mesh.positions.push_back({
                    radii[station] * std::cos(angle),
                    radii[station] * std::sin(angle),
                    axialStations[station]
                });
            }
        }
        for (std::size_t station = 0; station + 1 < axialStations.size(); ++station)
        {
            for (std::size_t sector = 0; sector < sectorCount; ++sector)
            {
                const std::size_t next = (sector + 1) % sectorCount;
                const int lower = static_cast<int>(station * sectorCount + sector);
                const int lowerNext = static_cast<int>(station * sectorCount + next);
                const int upper = static_cast<int>((station + 1) * sectorCount + sector);
                const int upperNext = static_cast<int>((station + 1) * sectorCount + next);
                mesh.triangles.push_back({ lower, lowerNext, upperNext });
                mesh.triangles.push_back({ lower, upperNext, upper });
            }
        }
        const int lowerCenter = static_cast<int>(mesh.positions.size());
        mesh.positions.push_back({ 0.0, 0.0, axialStations.front() });
        const int upperCenter = static_cast<int>(mesh.positions.size());
        mesh.positions.push_back({ 0.0, 0.0, axialStations.back() });
        const int upperOffset = static_cast<int>((axialStations.size() - 1) * sectorCount);
        for (std::size_t sector = 0; sector < sectorCount; ++sector)
        {
            const int current = static_cast<int>(sector);
            const int next = static_cast<int>((sector + 1) % sectorCount);
            mesh.triangles.push_back({ lowerCenter, next, current });
            mesh.triangles.push_back({ upperCenter, upperOffset + current, upperOffset + next });
        }
        return mesh;
    }

    void testTriangleMeshContract()
    {
        rotationbody::TriangleMesh mesh = makeTetrahedron();
        expect(mesh.validate().ok(), "finite indexed tetrahedron is valid");
        expect(mesh.vertexCount() == 4, "mesh reports vertex count");
        expect(mesh.triangleCount() == 4, "mesh reports triangle count");
        expect(near(mesh.scale(), std::sqrt(3.0)), "mesh scale is bounds diagonal");

        const std::uint64_t fingerprint = mesh.stableFingerprint();
        expect(fingerprint != 0, "mesh fingerprint is non-zero");
        expect(fingerprint == mesh.stableFingerprint(), "mesh fingerprint is stable");

        rotationbody::TransformComponents components;
        components.translationMeters = { 2.0, -3.0, 4.0 };
        const Eigen::Isometry3d transform = rotationbody::makeTransform(components);
        const auto transformed = mesh.transformed(transform);
        expect(transformed.ok(), "mesh accepts a finite rigid transform");
        if (transformed)
        {
            expect(nearVector(
                transformed.value.positions.front(),
                components.translationMeters),
                "mesh transform updates positions");
            expect(transformed.value.stableFingerprint() != fingerprint,
                "mesh fingerprint changes with geometry");
        }

        mesh.triangles.front()[0] = 99;
        expect(
            mesh.validate().error.code == rotationbody::PlanningErrorCode::InvalidTriangleIndex,
            "out-of-range triangle index is rejected");
        mesh = makeTetrahedron();
        mesh.positions.front().x() = std::numeric_limits<double>::quiet_NaN();
        expect(
            mesh.validate().error.code == rotationbody::PlanningErrorCode::NonFiniteGeometry,
            "NaN position is rejected");
    }

    void testTransformContract()
    {
        expect(near(rotationbody::millimetersToMeters(1250.0), 1.25),
            "millimeters convert to meters");
        expect(near(rotationbody::metersToMillimeters(1.25), 1250.0),
            "meters convert to millimeters");
        expect(near(rotationbody::radiansToDegrees(rotationbody::degreesToRadians(37.0)), 37.0),
            "degree and radian conversion round trips");
        expect(rotationbody::arePerpendicular(
            rotationbody::SignedAxis::NegativeX,
            rotationbody::SignedAxis::PositiveZ),
            "signed orthogonal axes are perpendicular");
        expect(!rotationbody::arePerpendicular(
            rotationbody::SignedAxis::PositiveY,
            rotationbody::SignedAxis::NegativeY),
            "parallel signed axes are not perpendicular");

        rotationbody::TransformComponents components;
        components.translationMeters = { 0.2, -0.4, 1.1 };
        components.rollPitchYawRadians = {
            rotationbody::degreesToRadians(90.0),
            rotationbody::degreesToRadians(20.0),
            rotationbody::degreesToRadians(-30.0)
        };
        const Eigen::Isometry3d transform = rotationbody::makeTransform(components);
        const Eigen::Matrix3d expected =
            (Eigen::AngleAxisd(components.rollPitchYawRadians.z(), Eigen::Vector3d::UnitZ()) *
             Eigen::AngleAxisd(components.rollPitchYawRadians.y(), Eigen::Vector3d::UnitY()) *
             Eigen::AngleAxisd(components.rollPitchYawRadians.x(), Eigen::Vector3d::UnitX()))
                .toRotationMatrix();
        expect((transform.linear() - expected).norm() <= 1.0e-12,
            "Euler composition is Rz * Ry * Rx");
        const auto roundTrip = rotationbody::transformComponents(transform);
        expect(roundTrip.ok(), "rigid transform decomposes into UI components");
        if (roundTrip)
        {
            const Eigen::Isometry3d rebuilt = rotationbody::makeTransform(roundTrip.value);
            expect((rebuilt.matrix() - transform.matrix()).norm() <= 1.0e-10,
                "transform component conversion round trips");
        }
    }

    void testCompletePartAlignment()
    {
        rotationbody::TriangleMesh rotary = makeRevolvedProfile(
            { 0.0, 0.55, 1.10, 1.70, 2.35, 3.0 },
            { 0.72, 0.42, 0.66, 0.70, 0.74, 0.78 });
        expect(rotary.validate().ok(), "synthetic rotary mesh is valid");

        rotationbody::TransformComponents rawPose;
        rawPose.translationMeters = { 2.4, -1.2, 0.7 };
        rawPose.rollPitchYawRadians = { 0.42, -0.31, 0.73 };
        const Eigen::Isometry3d rawFromModel = rotationbody::makeTransform(rawPose);
        const auto transformedMesh = rotary.transformed(rawFromModel);
        expect(transformedMesh.ok(), "arbitrarily posed rotary mesh is valid");
        if (!transformedMesh)
        {
            return;
        }

        const auto alignment = rotationbody::RotationBodyAlignmentSolver::solve(
            transformedMesh.value);
        expect(alignment.ok(), "complete part alignment succeeds");
        if (!alignment)
        {
            return;
        }
        expect(rotationbody::isFiniteTransform(alignment.value.planningFromMesh),
            "complete part alignment returns a finite rigid transform");
        const auto plannedBounds = transformedMesh.value.bounds(
            alignment.value.planningFromMesh);
        expect(plannedBounds.ok(), "aligned complete part has finite bounds");
        if (plannedBounds)
        {
            expect(near(plannedBounds.value.min().z(), 0.0, 1.0e-8),
                "selected bottom is placed on Z=0");
            expect(near(plannedBounds.value.max().z(), 3.0, 2.0e-2),
                "aligned height follows rotary axis");
        }
        const Eigen::Vector3d expectedDirectedAxis =
            -(rawFromModel.linear() * Eigen::Vector3d::UnitZ());
        expect(
            alignment.value.statistics.estimatedAxisInMesh.dot(expectedDirectedAxis) > 0.95,
            "minimum-diameter station selects the farther high end as bottom");
        expect(alignment.value.statistics.maximumDiameterMeters >
            alignment.value.statistics.minimumDiameterMeters,
            "complete part statistics preserve diameter contrast");

        rotationbody::TriangleMesh cylinder = makeRevolvedProfile(
            { 0.0, 0.6, 1.2, 1.8, 2.4 },
            { 0.5, 0.5, 0.5, 0.5, 0.5 });
        const auto cylinderAlignment = rotationbody::RotationBodyAlignmentSolver::solve(cylinder);
        expect(cylinderAlignment.ok(), "equal-diameter cylinder still returns an editable candidate");
        if (cylinderAlignment)
        {
            expect(cylinderAlignment.value.lowConfidence,
                "equal-diameter cylinder reports ambiguous bottom selection");
        }
    }

    void testSimulationBlockPlacementAndTransformOperations()
    {
        const rotationbody::TriangleMesh block = makeBox(
            { -0.10, -0.20, -0.30 },
            { 0.30, 0.40, 0.50 });
        const auto placement = rotationbody::SimulationBlockPlacementSolver::solve(
            block,
            rotationbody::SignedAxis::PositiveX,
            rotationbody::SignedAxis::PositiveZ,
            2.0);
        expect(placement.ok(), "simulation block placement accepts perpendicular axes");
        if (!placement)
        {
            return;
        }
        const auto bounds = block.bounds(placement.value.planningFromMesh);
        expect(bounds.ok(), "placed simulation block has finite bounds");
        if (bounds)
        {
            expect(near(0.5 * (bounds.value.min().x() + bounds.value.max().x()), 0.0),
                "simulation block thickness midpoint is X=0");
            expect(near(bounds.value.min().z(), 0.0),
                "simulation block minimum Z is zero");
            expect(near(bounds.value.max().y(), 1.0),
                "simulation block outer tooth is mother diameter divided by two");
        }
        expect(near(placement.value.planningFromMesh.linear().determinant(), 1.0),
            "simulation block basis is right handed");

        const auto invalidAxes = rotationbody::SimulationBlockPlacementSolver::solve(
            block,
            rotationbody::SignedAxis::PositiveY,
            rotationbody::SignedAxis::NegativeY,
            2.0);
        expect(invalidAxes.error.code == rotationbody::PlanningErrorCode::InvalidAxisSelection,
            "parallel simulation block axes are rejected");

        rotationbody::TransformComponents delta;
        delta.translationMeters = { 0.01, -0.02, 0.03 };
        delta.rollPitchYawRadians = { 0.1, 0.2, -0.15 };
        const auto edited = rotationbody::ModelTransformOperations::fromBaseline(
            placement.value.automaticBaseline,
            delta);
        expect(edited.ok(), "manual transform delta applies from automatic baseline");
        const auto reset = rotationbody::ModelTransformOperations::reset(
            placement.value.automaticBaseline);
        expect(reset.ok() &&
            (reset.value.matrix() - placement.value.automaticBaseline.matrix()).norm() <= 1.0e-12,
            "reset returns to automatic baseline rather than identity");
        const auto flipped = rotationbody::ModelTransformOperations::flip(
            block,
            placement.value.planningFromMesh);
        expect(flipped.ok(), "model can be flipped after automatic placement");
        if (flipped)
        {
            const auto flippedBounds = block.bounds(flipped.value);
            expect(flippedBounds.ok() && near(flippedBounds.value.min().z(), 0.0),
                "flipped model is re-seated on Z=0");
        }

        rotationbody::TransformComponents arbitraryTilt;
        arbitraryTilt.translationMeters = { 0.17, -0.08, 0.23 };
        arbitraryTilt.rollPitchYawRadians = { 0.41, -0.36, 0.27 };
        const auto tilted = rotationbody::ModelTransformOperations::applyIncrement(
            placement.value.planningFromMesh,
            arbitraryTilt);
        expect(tilted.ok(), "arbitrary manual tilt applies before flip");
        if (tilted)
        {
            const Eigen::Vector3d axisInMesh =
                placement.value.statistics.estimatedAxisInMesh.normalized();
            expect(std::abs((tilted.value.linear() * axisInMesh).z()) < 0.95,
                "rotated flip regression starts from a genuinely tilted rotary axis");
            const auto normalizedFlip = rotationbody::ModelTransformOperations::flip(
                block,
                tilted.value,
                axisInMesh,
                placement.value.bottomAxisCenterInMesh);
            expect(normalizedFlip.ok(), "flip accepts a genuinely tilted current rotary axis");
            if (normalizedFlip)
            {
                expect((normalizedFlip.value.linear() * axisInMesh +
                    Eigen::Vector3d::UnitZ()).norm() <= 1.0e-10,
                    "flip maps the old bottom-to-top axis to -Z");
                double oppositeEndDistance = -std::numeric_limits<double>::infinity();
                for (const Eigen::Vector3d& position : block.positions)
                {
                    oppositeEndDistance = std::max(
                        oppositeEndDistance,
                        (position - placement.value.bottomAxisCenterInMesh).dot(axisInMesh));
                }
                const Eigen::Vector3d newBottomInMesh =
                    placement.value.bottomAxisCenterInMesh +
                    oppositeEndDistance * axisInMesh;
                expect((normalizedFlip.value * newBottomInMesh).norm() <= 1.0e-10,
                    "flip places the opposite axial end center at the planning origin");
                const auto normalizedBounds = block.bounds(normalizedFlip.value);
                expect(normalizedBounds.ok() && near(normalizedBounds.value.min().z(), 0.0),
                    "normalized flip restores Zmin=0 after arbitrary rotation");
            }
        }
    }

    void testYzSectioning()
    {
        const rotationbody::TriangleMesh box = makeBox(
            { -1.0, -2.0, -3.0 },
            { 1.0, 2.0, 3.0 });
        const auto section = rotationbody::YzSectionExtractor::extract(
            box,
            Eigen::Isometry3d::Identity());
        expect(section.ok(), "box intersects the planning YZ plane");
        if (section)
        {
            const rotationbody::SectionContour& contour = section.value.targetContour;
            expect(contour.segmentCount() >= 3, "positive-Y box section is an ordered contour");
            expect(contour.points3d.size() == contour.pointsYz.size(),
                "section emits matching YZ and 3D points");
            expect(contour.cumulativeArcLength.size() == contour.pointsYz.size(),
                "section emits cumulative arc length at every point");
            for (std::size_t index = 0; index < contour.pointsYz.size(); ++index)
            {
                expect(contour.pointsYz[index].x() >= -contour.toleranceMeters,
                    "section target contains only Y>=0 points");
                expect(near(contour.points3d[index].x(), 0.0, contour.toleranceMeters),
                    "section 3D intersection lies on X=0");
                if (index > 0)
                {
                    expect((contour.pointsYz[index] - contour.pointsYz[index - 1]).norm() >
                        contour.toleranceMeters,
                        "section has no consecutive degenerate segment");
                }
            }
        }

        rotationbody::TransformComponents rawPose;
        rawPose.translationMeters = { 5.0, 0.4, -0.7 };
        rawPose.rollPitchYawRadians = { 0.0, 0.0, 0.0 };
        const auto movedBox = box.transformed(rotationbody::makeTransform(rawPose));
        expect(movedBox.ok(), "translated section test mesh is valid");
        if (movedBox)
        {
            const auto movedSection = rotationbody::YzSectionExtractor::extract(
                movedBox.value,
                rotationbody::makeTransform(rawPose).inverse());
            expect(movedSection.ok(), "planningFromMesh transform is applied before sectioning");
            if (section && movedSection)
            {
                expect(near(
                    movedSection.value.targetContour.pointsYz.front().x(),
                    section.value.targetContour.pointsYz.front().x(),
                    1.0e-7),
                    "transformed section matches direct section geometry");
            }
        }

        std::vector<rotationbody::YzSegment> topologySegments = {
            { { 0.0, 0.0 }, { 1.0, 0.0 } },
            { { 1.0, 0.0 }, { 1.0, 1.0 } },
            { { 1.0, 1.0 }, { 0.0, 1.0 } },
            { { 0.0, 1.0 }, { 0.0, 0.0 } },
            { { 1.0, 0.0 }, { 0.0, 0.0 } }
        };
        const auto topology = rotationbody::ContourTopology::buildContours(
            topologySegments,
            1.0e-8);
        expect(topology.ok(), "contour topology accepts unordered duplicate segments");
        if (topology)
        {
            expect(topology.value.discardedDuplicateSegmentCount == 1,
                "contour topology removes reversed duplicate segment");
            expect(!topology.value.contours.empty() && topology.value.contours.front().closed,
                "contour topology stitches a closed loop");
        }
    }

    rotationbody::SectionContour makeRegionContour()
    {
        rotationbody::SectionContour contour;
        contour.pointsYz = {
            { 1.0, 0.00 },
            { 1.0, 0.18 },
            { 1.03, 0.23 },
            { 1.0, 0.28 },
            { 1.22, 0.28 },
            { 1.22, 0.36 },
            { 1.0, 0.36 },
            { 1.0, 0.42 },
            { 1.24, 0.42 },
            { 1.24, 0.50 },
            { 1.0, 0.50 },
            { 1.0, 0.56 },
            { 1.21, 0.56 },
            { 1.21, 0.64 },
            { 1.0, 0.64 },
            { 1.03, 0.69 },
            { 1.0, 0.74 },
            { 1.0, 0.92 }
        };
        contour.points3d.reserve(contour.pointsYz.size());
        for (const Eigen::Vector2d& point : contour.pointsYz)
        {
            contour.points3d.push_back({ 0.0, point.x(), point.y() });
        }
        contour.closed = false;
        contour.toleranceMeters = 1.0e-8;
        return contour;
    }

    rotationbody::SectionContour makeNonUniformMultiClusterContour()
    {
        rotationbody::SectionContour contour;
        contour.pointsYz = {
            { 1.00, 0.00 }, { 1.00, 0.04 }, { 1.00, 0.17 }, { 1.00, 0.28 },
            { 1.02, 0.32 }, { 1.05, 0.35 }, { 1.00, 0.39 },
            { 1.22, 0.39 }, { 1.22, 0.49 }, { 1.00, 0.49 }, { 1.00, 0.57 },
            { 1.05, 0.61 }, { 1.02, 0.65 }, { 1.00, 0.70 },
            { 1.00, 0.82 }, { 1.00, 1.03 },
            { 1.03, 1.07 }, { 1.06, 1.11 }, { 1.00, 1.15 },
            { 1.27, 1.15 }, { 1.27, 1.24 }, { 1.00, 1.24 }, { 1.00, 1.33 },
            { 1.06, 1.37 }, { 1.02, 1.42 }, { 1.00, 1.48 },
            { 1.00, 1.61 }, { 1.00, 1.86 }
        };
        for (const Eigen::Vector2d& point : contour.pointsYz)
        {
            contour.points3d.push_back({ 0.0, point.x(), point.y() });
        }
        contour.closed = false;
        contour.toleranceMeters = 1.0e-10;
        return contour;
    }

    rotationbody::SectionContour makeClosedRegionContour()
    {
        rotationbody::SectionContour contour = makeRegionContour();
        contour.pointsYz.push_back({ 0.72, 0.92 });
        contour.pointsYz.push_back({ 0.72, 0.00 });
        contour.points3d.clear();
        for (const Eigen::Vector2d& point : contour.pointsYz)
        {
            contour.points3d.push_back({ 0.0, point.x(), point.y() });
        }
        contour.closed = true;
        return contour;
    }

    rotationbody::SectionContour rotateClosedContour(
        rotationbody::SectionContour contour,
        std::size_t offset)
    {
        offset %= contour.pointsYz.size();
        std::rotate(
            contour.pointsYz.begin(),
            contour.pointsYz.begin() + static_cast<std::ptrdiff_t>(offset),
            contour.pointsYz.end());
        std::rotate(
            contour.points3d.begin(),
            contour.points3d.begin() + static_cast<std::ptrdiff_t>(offset),
            contour.points3d.end());
        contour.cumulativeArcLength.clear();
        return contour;
    }

    std::array<std::size_t, 5> regionLabelCounts(
        const rotationbody::RegionAssignment& assignment)
    {
        std::array<std::size_t, 5> counts{};
        for (rotationbody::RegionLabel label : assignment.segmentLabels)
        {
            ++counts[static_cast<std::size_t>(label)];
        }
        return counts;
    }

    rotationbody::SectionContour makeSharpShoulderContour()
    {
        rotationbody::SectionContour contour;
        contour.pointsYz = {
            { 1.00, 0.00 }, { 1.00, 0.20 },
            { 1.22, 0.20 }, { 1.22, 0.27 }, { 1.00, 0.27 }, { 1.00, 0.33 },
            { 1.22, 0.33 }, { 1.22, 0.40 }, { 1.00, 0.40 }, { 1.00, 0.46 },
            { 1.22, 0.46 }, { 1.22, 0.53 }, { 1.00, 0.53 },
            { 1.00, 0.82 }
        };
        for (const Eigen::Vector2d& point : contour.pointsYz)
        {
            contour.points3d.push_back({ 0.0, point.x(), point.y() });
        }
        contour.closed = false;
        contour.toleranceMeters = 1.0e-10;
        return contour;
    }

    void testAutomaticRegionRecognition()
    {
        const rotationbody::SectionContour contour = makeRegionContour();
        const auto recognition = rotationbody::ToothRegionRecognizer::recognize(contour);
        expect(recognition.ok(), "tooth profile recognition succeeds");
        if (!recognition)
        {
            return;
        }
        expect(recognition.value.matches(contour), "automatic labels match contour segments");
        expect(
            recognition.value.segmentLabels.front() == rotationbody::RegionLabel::Unclassified &&
                recognition.value.segmentLabels.back() == rotationbody::RegionLabel::Unclassified,
            "coarse straight body segments are not promoted by a small transition overlap");

        rotationbody::ToothRecognitionOptions coarseShoulderOptions;
        coarseShoulderOptions.minimumResampleCount = 64;
        coarseShoulderOptions.maximumResampleCount = 64;
        coarseShoulderOptions.smoothingWindowFraction = 0.025;
        const auto coarseShoulderRecognition =
            rotationbody::ToothRegionRecognizer::recognize(
                contour,
                coarseShoulderOptions);
        expect(coarseShoulderRecognition.ok(),
            "coarse shoulder recognition succeeds");
        if (coarseShoulderRecognition)
        {
            expect(
                coarseShoulderRecognition.value.segmentLabels[1] ==
                    rotationbody::RegionLabel::Unclassified,
                "a mixed coarse shoulder remains Unclassified without a strict majority");
        }
        for (rotationbody::RegionLabel required : {
                 rotationbody::RegionLabel::ToothTop,
                 rotationbody::RegionLabel::ToothWall,
                 rotationbody::RegionLabel::ToothBottom,
                 rotationbody::RegionLabel::Transition })
        {
            expect(
                std::find(
                    recognition.value.segmentLabels.begin(),
                    recognition.value.segmentLabels.end(),
                    required) != recognition.value.segmentLabels.end(),
                "synthetic teeth contain every automatic spray category");
        }
        for (std::size_t wallSegment : { 3u, 5u, 7u, 9u, 11u, 13u })
        {
            expect(
                recognition.value.segmentLabels[wallSegment] ==
                    rotationbody::RegionLabel::ToothWall,
                "each tooth period keeps both complete radial walls");
        }
        for (std::size_t topSegment : { 4u, 8u, 12u })
        {
            expect(
                recognition.value.segmentLabels[topSegment] ==
                    rotationbody::RegionLabel::ToothTop,
                "only the local radial maxima remain ToothTop");
        }
        for (std::size_t bottomSegment : { 6u, 10u })
        {
            expect(
                recognition.value.segmentLabels[bottomSegment] ==
                    rotationbody::RegionLabel::ToothBottom,
                "only the axial-facing valley floor remains ToothBottom");
        }

        rotationbody::SectionContour reversed = contour;
        std::reverse(reversed.pointsYz.begin(), reversed.pointsYz.end());
        std::reverse(reversed.points3d.begin(), reversed.points3d.end());
        const auto reverseRecognition = rotationbody::ToothRegionRecognizer::recognize(reversed);
        expect(reverseRecognition.ok(), "recognition accepts reversed point order");
        if (reverseRecognition)
        {
            for (rotationbody::RegionLabel required : {
                     rotationbody::RegionLabel::ToothTop,
                     rotationbody::RegionLabel::ToothWall,
                     rotationbody::RegionLabel::ToothBottom })
            {
                expect(
                    std::find(
                        reverseRecognition.value.segmentLabels.begin(),
                        reverseRecognition.value.segmentLabels.end(),
                        required) != reverseRecognition.value.segmentLabels.end(),
                    "reversed teeth preserve geometric categories");
            }
        }

        const rotationbody::SectionContour multiCluster =
            makeNonUniformMultiClusterContour();
        const auto multiRecognition =
            rotationbody::ToothRegionRecognizer::recognize(multiCluster);
        expect(multiRecognition.ok(), "non-uniform multi-cluster recognition succeeds");
        if (multiRecognition)
        {
            const auto& labels = multiRecognition.value.segmentLabels;
            const auto isTooth = [](rotationbody::RegionLabel label)
            {
                return label == rotationbody::RegionLabel::ToothTop ||
                    label == rotationbody::RegionLabel::ToothWall ||
                    label == rotationbody::RegionLabel::ToothBottom;
            };

            std::size_t toothRunCount = 0;
            std::size_t nonToothGap = 2;
            bool inToothRun = false;
            for (rotationbody::RegionLabel label : labels)
            {
                if (isTooth(label))
                {
                    if (!inToothRun && nonToothGap >= 2)
                    {
                        ++toothRunCount;
                    }
                    inToothRun = true;
                    nonToothGap = 0;
                }
                else
                {
                    ++nonToothGap;
                    if (nonToothGap >= 2)
                    {
                        inToothRun = false;
                    }
                }
            }
            expect(toothRunCount >= 2,
                "non-uniform profile retains two distinct tooth clusters");

            for (std::size_t index = 14; index <= 15; ++index)
            {
                expect(labels[index] != rotationbody::RegionLabel::Transition,
                    "stable body gap between tooth clusters is not marked Transition");
            }

            std::size_t transitionRunCount = 0;
            bool inTransition = false;
            for (rotationbody::RegionLabel label : labels)
            {
                if (label == rotationbody::RegionLabel::Transition && !inTransition)
                {
                    ++transitionRunCount;
                    inTransition = true;
                }
                else if (label != rotationbody::RegionLabel::Transition)
                {
                    inTransition = false;
                }
            }
            expect(transitionRunCount >= 2,
                "multi-cluster transitions stay localized at cluster boundaries");
        }

        rotationbody::ToothRecognitionOptions mixedClassOptions;
        mixedClassOptions.minimumResampleCount = 128;
        mixedClassOptions.maximumResampleCount = 128;
        mixedClassOptions.smoothingWindowFraction = 0.04;
        const auto mixedClassRecognition =
            rotationbody::ToothRegionRecognizer::recognize(
                multiCluster,
                mixedClassOptions);
        expect(mixedClassRecognition.ok(),
            "mixed-class coarse recognition succeeds");
        if (mixedClassRecognition)
        {
            expect(
                mixedClassRecognition.value.segmentLabels[10] ==
                    rotationbody::RegionLabel::Unclassified,
                "a multi-class coarse segment remains Unclassified without a strict majority");
        }

        const rotationbody::SectionContour closed = makeClosedRegionContour();
        constexpr std::size_t cyclicOffset = 7;
        const rotationbody::SectionContour shifted =
            rotateClosedContour(closed, cyclicOffset);
        const auto closedRecognition =
            rotationbody::ToothRegionRecognizer::recognize(closed);
        const auto shiftedRecognition =
            rotationbody::ToothRegionRecognizer::recognize(shifted);
        expect(closedRecognition.ok() && shiftedRecognition.ok(),
            "closed tooth recognition accepts cyclically equivalent contours");
        if (closedRecognition && shiftedRecognition)
        {
            expect(
                regionLabelCounts(closedRecognition.value) ==
                    regionLabelCounts(shiftedRecognition.value),
                "closed contour label statistics are invariant under cyclic shifts");
            bool labelsFollowGeometry = true;
            bool confidenceFollowsGeometry = true;
            for (std::size_t shiftedIndex = 0;
                 shiftedIndex < shifted.segmentCount();
                 ++shiftedIndex)
            {
                const std::size_t originalIndex =
                    (shiftedIndex + cyclicOffset) % closed.segmentCount();
                labelsFollowGeometry = labelsFollowGeometry &&
                    shiftedRecognition.value.segmentLabels[shiftedIndex] ==
                        closedRecognition.value.segmentLabels[originalIndex];
                confidenceFollowsGeometry = confidenceFollowsGeometry && near(
                    shiftedRecognition.value.segmentConfidence[shiftedIndex],
                    closedRecognition.value.segmentConfidence[originalIndex],
                    1.0e-10);
            }
            expect(labelsFollowGeometry,
                "closed contour labels remain attached to the same geometric segments");
            expect(confidenceFollowsGeometry,
                "closed contour confidence remains attached to the same geometric segments");
        }

        const rotationbody::SectionContour sharpShoulder =
            makeSharpShoulderContour();
        const auto sharpRecognition =
            rotationbody::ToothRegionRecognizer::recognize(sharpShoulder);
        expect(sharpRecognition.ok(), "sharp-shoulder tooth recognition succeeds");
        if (sharpRecognition)
        {
            expect(
                std::count(
                    sharpRecognition.value.segmentLabels.begin(),
                    sharpRecognition.value.segmentLabels.end(),
                    rotationbody::RegionLabel::Transition) == 0,
                "straight body beside a sharp tooth shoulder is not invented as Transition");
        }
    }

    void testRegionEditHistory()
    {
        rotationbody::SectionContour contour;
        contour.pointsYz = {
            { 1.0, 0.0 },
            { 2.0, 0.0 },
            { 2.0, 1.0 },
            { 1.0, 1.0 }
        };
        contour.closed = false;
        contour.toleranceMeters = 1.0e-9;
        rotationbody::RegionAssignment automatic;
        automatic.segmentLabels.assign(contour.segmentCount(), rotationbody::RegionLabel::ToothTop);
        automatic.segmentConfidence.assign(contour.segmentCount(), 0.8);

        rotationbody::RegionEditHistory history(contour, automatic);
        rotationbody::YzRectangle first;
        first.minimum = { 1.3, -0.1 };
        first.maximum = { 2.1, 0.2 };
        expect(history.applyRectangle(first, rotationbody::RegionLabel::ToothWall).ok(),
            "first rectangle edit intersects contour");
        rotationbody::YzRectangle second = first;
        second.minimum.x() = 1.5;
        expect(history.applyRectangle(second, rotationbody::RegionLabel::Unclassified).ok(),
            "overlapping later rectangle edit succeeds");
        expect(history.resolved().segmentLabels.front() == rotationbody::RegionLabel::Unclassified,
            "later rectangle wins on overlapping segments");
        expect(history.undo(), "region edit can undo");
        expect(history.resolved().segmentLabels.front() == rotationbody::RegionLabel::ToothWall,
            "undo restores previous rectangle label");
        expect(history.redo(), "region edit can redo");
        expect(history.resolved().segmentLabels.front() == rotationbody::RegionLabel::Unclassified,
            "redo reapplies later label");
        history.restoreAutomatic();
        expect(history.resolved().segmentLabels.front() == rotationbody::RegionLabel::ToothTop,
            "restoreAutomatic clears all overrides");
        expect(!history.canUndo() && !history.canRedo(),
            "restoreAutomatic clears edit history cursor");
    }

    void testRegionEditScaleInvariance()
    {
        for (double scale : { 1.0e-9, 1.0, 1.0e6 })
        {
            rotationbody::SectionContour contour;
            contour.pointsYz = {
                { 0.0, 0.0 },
                { scale, scale }
            };
            contour.closed = false;
            contour.toleranceMeters = scale * 1.0e-12;
            rotationbody::RegionAssignment automatic;
            automatic.segmentLabels = { rotationbody::RegionLabel::ToothTop };
            automatic.segmentConfidence = { 0.8 };

            rotationbody::RegionEditHistory history(contour, automatic);
            rotationbody::YzRectangle crossing;
            crossing.minimum = { 0.45 * scale, 0.40 * scale };
            crossing.maximum = { 0.55 * scale, 0.60 * scale };
            const auto edit = history.applyRectangle(
                crossing,
                rotationbody::RegionLabel::Unclassified);
            expect(edit.ok(), "rectangle intersection is invariant across coordinate scales");
            if (edit)
            {
                expect(history.resolved().segmentLabels.front() ==
                    rotationbody::RegionLabel::Unclassified,
                    "scaled rectangle overrides a segment through its interior");
            }

            rotationbody::RegionEditHistory missHistory(contour, automatic);
            rotationbody::YzRectangle separate;
            separate.minimum = { 0.45 * scale, 0.70 * scale };
            separate.maximum = { 0.55 * scale, 0.80 * scale };
            expect(!missHistory.applyRectangle(
                    separate,
                    rotationbody::RegionLabel::Unclassified).ok(),
                "scaled rectangle does not hit a separated segment");
        }
    }

    void testSprayBoundaries()
    {
        rotationbody::SectionContour contour;
        contour.pointsYz = {
            { 9.0, -4.0 },
            { 1.0, 0.0 },
            { 2.0, 0.0 },
            { 2.1, 1.0 },
            { 2.2, 2.0 },
            { 1.0, 2.0 }
        };
        contour.closed = false;
        contour.toleranceMeters = 1.0e-10;
        rotationbody::RegionAssignment assignment;
        assignment.segmentLabels = {
            rotationbody::RegionLabel::Unclassified,
            rotationbody::RegionLabel::Transition,
            rotationbody::RegionLabel::ToothTop,
            rotationbody::RegionLabel::ToothTop,
            rotationbody::RegionLabel::ToothBottom
        };

        const auto maximumY = rotationbody::SprayBoundaryBuilder::build(
            contour,
            assignment,
            rotationbody::BoundaryMode::MaximumToothTopY);
        expect(maximumY.ok(), "maximum tooth-top boundary succeeds");
        if (maximumY)
        {
            expect(near(maximumY.value.minimumY, 1.0),
                "unclassified radial outlier is excluded from minimum Y");
            expect(near(maximumY.value.maximumY, 2.2),
                "maximum Y is exactly the maximum tooth-top Y");
            expect(near(maximumY.value.minimumZ, 0.0) && near(maximumY.value.maximumZ, 2.0),
                "boundary uses tight classified Z bounds without padding");
            expect(near(maximumY.value.polygonYz[1].x(), 2.2),
                "straight outer edge contains no one-percent expansion");
        }

        const auto envelope = rotationbody::SprayBoundaryBuilder::build(
            contour,
            assignment,
            rotationbody::BoundaryMode::ToothTopEnvelope);
        expect(envelope.ok(), "tooth-top envelope boundary succeeds");
        if (envelope)
        {
            expect(near(envelope.value.outerLineSlopeYPerZ, 0.1, 1.0e-8),
                "envelope recovers slanted tooth-top slope");
            expect(near(envelope.value.outerY(0.0), 2.0, 1.0e-8),
                "envelope lower endpoint is tight");
            expect(near(envelope.value.outerY(2.0), 2.2, 1.0e-8),
                "envelope upper endpoint is tight");
        }

        assignment.segmentLabels[2] = rotationbody::RegionLabel::Unclassified;
        assignment.segmentLabels[3] = rotationbody::RegionLabel::Unclassified;
        const auto missingTop = rotationbody::SprayBoundaryBuilder::build(
            contour,
            assignment,
            rotationbody::BoundaryMode::MaximumToothTopY);
        expect(
            missingTop.error.code == rotationbody::PlanningErrorCode::InsufficientToothTopData,
            "boundary rejects assignments without tooth-top data");
    }

    void testWorkpieceCalibration()
    {
        const Eigen::Vector3d expectedAxis =
            Eigen::Vector3d(0.25, -0.35, 0.902).normalized();
        const Eigen::Vector3d expectedAxisPoint(0.2, -0.1, 0.3);
        Eigen::Vector3d basisFirst = Eigen::Vector3d::UnitZ().cross(expectedAxis);
        if(basisFirst.norm() < 1.0e-9) {
            basisFirst = Eigen::Vector3d::UnitX().cross(expectedAxis);
        }
        basisFirst.normalize();
        const Eigen::Vector3d basisSecond = expectedAxis.cross(basisFirst).normalized();

        rotationbody::CalibrationPointList cylinderPoints;
        constexpr double radius = 0.05;
        for(std::size_t index = 0; index < 12; ++index)
        {
            const double angle = 2.0 * rotationbody::pi *
                static_cast<double>(index % 6) / 6.0;
            const double axial = index < 6 ? -0.08 : 0.09;
            cylinderPoints.push_back(
                expectedAxisPoint + axial * expectedAxis +
                radius * (std::cos(angle) * basisFirst +
                    std::sin(angle) * basisSecond));
        }
        const auto cylinder =
            rotationbody::WorkpieceCalibrationSolver::fitCylinder3d(cylinderPoints);
        expect(cylinder.ok(), "twelve side touch points fit a spatial cylinder");
        if(cylinder)
        {
            expect(std::abs(cylinder.value.axisDirectionBase.dot(expectedAxis)) > 0.9999,
                "spatial cylinder fit recovers the tilted rotary axis");
            expect(std::abs(cylinder.value.radiusMeters - radius) < 1.0e-6,
                "spatial cylinder fit recovers the cylinder radius");
            const Eigen::Vector3d axisOffset =
                cylinder.value.axisPointBaseMeters - expectedAxisPoint;
            expect((axisOffset - axisOffset.dot(expectedAxis) * expectedAxis).norm() < 1.0e-6,
                "spatial cylinder fit recovers the cylinder axis line");
        }

        rotationbody::CalibrationPointList circlePoints;
        const Eigen::Vector3d expectedCenter(0.4, -0.25, 0.7);
        for(std::size_t index = 0; index < 6; ++index)
        {
            const double angle = 2.0 * rotationbody::pi *
                static_cast<double>(index) / 6.0;
            circlePoints.emplace_back(
                expectedCenter.x() + 0.12 * std::cos(angle),
                expectedCenter.y() + 0.12 * std::sin(angle),
                expectedCenter.z() + (index % 2 == 0 ? 0.00002 : -0.00002));
        }
        const auto circle =
            rotationbody::WorkpieceCalibrationSolver::fitCircle2d(circlePoints);
        expect(circle.ok(), "six coplanar touch points fit a horizontal circle");
        if(circle)
        {
            expect((circle.value.axisPointBaseMeters - expectedCenter).norm() < 1.0e-6,
                "horizontal circle fit recovers its center");
            expect(std::abs(circle.value.radiusMeters - 0.12) < 1.0e-9,
                "horizontal circle fit recovers its radius");
        }

        rotationbody::CalibrationPointList nonPlanar = circlePoints;
        nonPlanar.back().z() += 0.001;
        expect(
            !rotationbody::WorkpieceCalibrationSolver::fitCircle2d(nonPlanar),
            "horizontal circle fit rejects touch points outside the Z tolerance");

        rotationbody::CalibrationAxisFit axisFit;
        axisFit.mode = rotationbody::CalibrationMode::Cylinder3d;
        axisFit.axisPointBaseMeters = Eigen::Vector3d(0.5, -0.2, 1.0);
        axisFit.axisDirectionBase = Eigen::Vector3d::UnitZ();
        axisFit.radiusMeters = 0.1;
        const auto frame =
            rotationbody::WorkpieceCalibrationSolver::computeWorkpieceFrame(
                axisFit,
                Eigen::Vector3d(0.6, -0.2, 1.3),
                0.2,
                Eigen::Vector3d(0.5, -0.2, 1.1),
                Eigen::Vector3d(0.5, 0.0, 1.1));
        expect(frame.ok(), "axis, tooth-top point, Y baseline and height build a frame");
        if(frame)
        {
            expect(nearVector(
                frame.value.baseFromPlanning.translation(),
                Eigen::Vector3d(0.5, -0.2, 1.1)),
                "calibration places the workpiece origin at the projected bottom center");
            expect(frame.value.baseFromPlanning.linear().isApprox(
                Eigen::Matrix3d::Identity(), 1.0e-9),
                "calibration constructs X/Y/Z with the platform Rz-Ry-Rx convention");
            expect(nearVector(
                frame.value.baseFromPlanningComponents.translationMeters,
                Eigen::Vector3d(0.5, -0.2, 1.1)),
                "calibration exposes components for the base-coordinate pose fields");
        }

        rotationbody::CalibrationAxisFit tiltedAxisFit;
        tiltedAxisFit.mode = rotationbody::CalibrationMode::Cylinder3d;
        const Eigen::Vector3d tiltedZ =
            Eigen::Vector3d(0.3, -0.4, 0.866025403784).normalized();
        Eigen::Vector3d tiltedYSeed(-0.2, 0.9, 0.35);
        Eigen::Vector3d tiltedY =
            tiltedYSeed - tiltedYSeed.dot(tiltedZ) * tiltedZ;
        tiltedY.normalize();
        const Eigen::Vector3d tiltedX = tiltedY.cross(tiltedZ).normalized();
        tiltedY = tiltedZ.cross(tiltedX).normalized();
        const Eigen::Vector3d tiltedOrigin(0.42, -0.18, 0.76);
        const Eigen::Vector3d tiltedTopCenter = tiltedOrigin + 0.24 * tiltedZ;
        tiltedAxisFit.axisPointBaseMeters = tiltedOrigin + 0.08 * tiltedZ;
        tiltedAxisFit.axisDirectionBase = tiltedZ;
        tiltedAxisFit.radiusMeters = 0.08;
        const auto tiltedFrame =
            rotationbody::WorkpieceCalibrationSolver::computeWorkpieceFrame(
                tiltedAxisFit,
                tiltedTopCenter + 0.08 * tiltedX,
                0.24,
                tiltedOrigin - 0.06 * tiltedY - 0.01 * tiltedZ,
                tiltedOrigin + 0.06 * tiltedY + 0.02 * tiltedZ);
        expect(tiltedFrame.ok(), "tilted axis and nontrivial Y baseline build a frame");
        if(tiltedFrame)
        {
            Eigen::Matrix3d expectedRotation;
            expectedRotation.col(0) = tiltedX;
            expectedRotation.col(1) = tiltedY;
            expectedRotation.col(2) = tiltedZ;
            expect(tiltedFrame.value.baseFromPlanning.translation().isApprox(
                tiltedOrigin, 1.0e-9),
                "tilted calibration preserves the projected bottom-center origin");
            expect(tiltedFrame.value.baseFromPlanning.linear().isApprox(
                expectedRotation, 1.0e-9),
                "tilted calibration constructs the expected workpiece axes");
            expect(rotationbody::makeTransform(
                tiltedFrame.value.baseFromPlanningComponents).matrix().isApprox(
                    tiltedFrame.value.baseFromPlanning.matrix(), 1.0e-9),
                "nonzero calibration Rx/Ry/Rz round-trip with Rz-Ry-Rx ordering");
        }
        expect(
            !rotationbody::WorkpieceCalibrationSolver::computeWorkpieceFrame(
                axisFit,
                Eigen::Vector3d(0.6, -0.2, 1.3),
                0.2,
                Eigen::Vector3d::Zero(),
                Eigen::Vector3d(0.0, 0.0, 0.0005)),
            "calibration rejects a degenerate Y-direction baseline");
    }

    rotationbody::SprayBoundary makeTrajectoryBoundary()
    {
        rotationbody::SprayBoundary boundary;
        boundary.mode = rotationbody::BoundaryMode::MaximumToothTopY;
        boundary.minimumY = 0.01;
        boundary.maximumY = 0.02;
        boundary.minimumZ = 0.0;
        boundary.maximumZ = 0.01;
        boundary.outerLineSlopeYPerZ = 0.0;
        boundary.outerLineInterceptY = 0.02;
        boundary.polygonYz = {
            { 0.01, 0.0 },
            { 0.02, 0.0 },
            { 0.02, 0.01 },
            { 0.01, 0.01 }
        };
        return boundary;
    }

    rotationbody::SprayBoundary makeSlopedTrajectoryBoundary()
    {
        rotationbody::SprayBoundary boundary = makeTrajectoryBoundary();
        boundary.outerLineSlopeYPerZ = 0.5;
        boundary.outerLineInterceptY = 0.02;
        boundary.maximumY = boundary.outerY(boundary.maximumZ);
        boundary.polygonYz = {
            { 0.01, 0.0 },
            { boundary.outerY(boundary.minimumZ), boundary.minimumZ },
            { boundary.outerY(boundary.maximumZ), boundary.maximumZ },
            { 0.01, boundary.maximumZ }
        };
        boundary.mode = rotationbody::BoundaryMode::ToothTopEnvelope;
        return boundary;
    }

    void testTrajectoryPlanningAndEditing()
    {
        rotationbody::TrajectoryGenerationParameters parameters;
        parameters.sprayDistanceMeters = 0.05;
        parameters.speedMetersPerSecond = 0.01;
        parameters.positionerRpm = 6.0;
        const auto generated = rotationbody::TrajectoryPlanner::generate(
            makeTrajectoryBoundary(),
            parameters);
        expect(generated.ok(), "trajectory generation accepts a confirmed boundary");
        if(!generated)
        {
            return;
        }

        rotationbody::PlannedTrajectory trajectory = generated.value;
        expect(trajectory.linearPoints.size() == 11,
            "automatic trajectory sampling uses at most one millimeter spacing");
        expect(trajectory.relativeHelicalPoints.size() == trajectory.linearPoints.size(),
            "linear and relative helical trajectories have one-to-one points");
        expect(near(trajectory.linearPoints.front().planningFromTool.translation().z(), 0.01),
            "the larger planning-local Z endpoint is the initial A point");
        expect(nearVector(
            trajectory.linearPoints.front().planningFromTool.linear().col(2),
            -Eigen::Vector3d::UnitY()),
            "zero tilt points the tool Z spray axis normal to a horizontal outer edge");

        rotationbody::TrajectoryGenerationParameters tiltedParameters = parameters;
        const auto slopedZero = rotationbody::TrajectoryPlanner::generate(
            makeSlopedTrajectoryBoundary(),
            tiltedParameters);
        expect(slopedZero.ok(),
            "sloped outer edge accepts a zero relative-normal tilt");
        tiltedParameters.tiltRadians = 30.0 * std::acos(-1.0) / 180.0;
        const auto tilted = rotationbody::TrajectoryPlanner::generate(
            makeSlopedTrajectoryBoundary(),
            tiltedParameters);
        expect(tilted.ok(),
            "sloped outer edge accepts a local-Y trajectory tilt");
        if(tilted)
        {
            const Eigen::Vector3d edgeTangent =
                (Eigen::Vector3d(0.0, 0.02, 0.0) -
                    Eigen::Vector3d(0.0, 0.025, 0.01)).normalized();
            Eigen::Vector3d expectedNormal(
                0.0,
                -edgeTangent.z(),
                edgeTangent.y());
            expectedNormal.normalize();
            if(expectedNormal.y() > 0.0) expectedNormal = -expectedNormal;
            const Eigen::Vector3d expectedY =
                expectedNormal.cross(edgeTangent).normalized();
            const Eigen::Vector3d expectedSpray =
                Eigen::AngleAxisd(tiltedParameters.tiltRadians, expectedY) *
                expectedNormal;
            if(slopedZero)
            {
                const Eigen::Vector3d zeroSpray =
                    slopedZero.value.linearPoints.front().planningFromTool.linear().col(2);
                expect(near(zeroSpray.dot(edgeTangent), 0.0, 1.0e-9) &&
                    zeroSpray.y() < 0.0,
                    "zero tilt is perpendicular to the fitted edge and points into the workpiece");
                expect(near(
                    std::acos(std::clamp(zeroSpray.dot(expectedSpray), -1.0, 1.0)),
                    tiltedParameters.tiltRadians,
                    1.0e-9) && expectedSpray.z() < zeroSpray.z(),
                    "positive tilt rotates the spray Z axis downward around local Y");
            }
            expect(nearVector(
                tilted.value.linearPoints.front().planningFromTool.linear().col(2),
                expectedSpray,
                1.0e-9),
                "tilt is measured from the fitted outer-edge normal");

            tiltedParameters.reversed = true;
            const auto reversedTilted = rotationbody::TrajectoryPlanner::generate(
                makeSlopedTrajectoryBoundary(),
                tiltedParameters);
            expect(reversedTilted.ok(),
                "reversing a tilted sloped trajectory remains valid");
            if(reversedTilted)
            {
                expect(nearVector(
                    reversedTilted.value.linearPoints.front().planningFromTool.linear().col(2),
                    expectedSpray,
                    1.0e-9),
                    "swapping A/B preserves the physical spray direction");
                expect(
                    reversedTilted.value.linearPoints.front().planningFromTool.linear().isApprox(
                        tilted.value.linearPoints.back().planningFromTool.linear(),
                        1.0e-9) &&
                    reversedTilted.value.linearPoints.back().planningFromTool.linear().isApprox(
                        tilted.value.linearPoints.front().planningFromTool.linear(),
                        1.0e-9),
                    "reversing generation preserves the complete tool orientation");
            }
        }
        bool timestampsMatch = true;
        for(std::size_t index = 0; index < trajectory.linearPoints.size(); ++index)
        {
            timestampsMatch = timestampsMatch && near(
                trajectory.linearPoints[index].timeSeconds,
                trajectory.relativeHelicalPoints[index].timeSeconds);
        }
        expect(timestampsMatch,
            "relative helical points preserve every linear trajectory timestamp");

        const Eigen::Matrix3d forwardStartOrientation =
            trajectory.linearPoints.front().planningFromTool.linear();
        const Eigen::Matrix3d forwardEndOrientation =
            trajectory.linearPoints.back().planningFromTool.linear();
        expect(rotationbody::TrajectoryEditor::swapDirection(trajectory).ok(),
            "trajectory direction can be swapped");
        expect(trajectory.parameters.reversed &&
            near(trajectory.linearPoints.front().planningFromTool.translation().z(), 0.0) &&
            trajectory.linearPoints.front().planningFromTool.linear().isApprox(
                forwardEndOrientation, 1.0e-9) &&
            trajectory.linearPoints.back().planningFromTool.linear().isApprox(
                forwardStartOrientation, 1.0e-9),
            "swapping exchanges A and B without changing tool orientation");

        const std::size_t previousCount = trajectory.linearPoints.size();
        expect(rotationbody::TrajectoryEditor::interpolateRange(
            trajectory,
            { 2, 5 },
            0.01).ok(),
            "a Shift/Ctrl selection range can be time-interpolated");
        expect(trajectory.linearPoints.size() > previousCount &&
            trajectory.relativeHelicalPoints.size() == trajectory.linearPoints.size(),
            "interpolation densifies both trajectories with aligned indices");

        const Eigen::Vector3d untouched =
            trajectory.linearPoints.front().planningFromTool.translation();
        const Eigen::Vector3d selectedBefore =
            trajectory.linearPoints[1].planningFromTool.translation();
        rotationbody::TransformComponents delta;
        delta.translationMeters = { 0.001, 0.0, 0.0 };
        delta.rollPitchYawRadians = { 0.0, 0.0, 0.1 };
        expect(rotationbody::TrajectoryEditor::transformSelected(
            trajectory,
            { 1 },
            delta).ok(),
            "selected trajectory poses accept local pose adjustment");
        expect(nearVector(
            trajectory.linearPoints.front().planningFromTool.translation(),
            untouched) && nearVector(
                trajectory.linearPoints[1].planningFromTool.translation(),
                selectedBefore + delta.translationMeters),
            "pose adjustment affects only explicitly selected indices");
    }

    void testTrajectoryGroupAndRapidTranslation()
    {
        rotationbody::TrajectoryGenerationParameters parameters;
        parameters.sprayDistanceMeters = 0.05;
        parameters.speedMetersPerSecond = 0.01;
        const auto generated = rotationbody::TrajectoryPlanner::generate(
            makeTrajectoryBoundary(),
            parameters);
        expect(generated.ok(), "RAPID fixture trajectory can be generated");
        if(!generated)
        {
            return;
        }

        rotationbody::TrajectoryGroup group;
        const auto saved = rotationbody::TrajectoryGroupEditor::addOrUpdate(
            group,
            generated.value);
        expect(saved.ok() && rotationbody::TrajectoryGroupEditor::validate(group).ok(),
            "a finalized trajectory can be saved to a valid trajectory group");
        if(!saved)
        {
            return;
        }

        rotationbody::PublishedTrajectoryPlan plan;
        plan.objectId = "workpiece-1";
        plan.baseFromPlanning.translation() = Eigen::Vector3d(1.0, 2.0, 3.0);
        plan.group = group;
        rotationbody::RapidExportSettings settings;
        settings.safetyPositionBaseMeters = { 1.2, 2.3, 3.4 };
        std::vector<rotationbody::RapidSequenceEntry> sequence(2);
        sequence[0].kind = rotationbody::RapidSequenceEntryKind::SafetyPoint;
        sequence[1].kind = rotationbody::RapidSequenceEntryKind::Trajectory;
        sequence[1].trajectoryPassId = saved.value;
        const auto rapid = rotationbody::RapidModuleGenerator::generate(
            plan,
            settings,
            sequence);
        expect(rapid.ok(), "ABB RAPID translation accepts safety and trajectory entries");
        if(rapid)
        {
            expect(rapid.value.code.find("MODULE SprayRotation") != std::string::npos &&
                rapid.value.code.find("ConfJ \\Off;") != std::string::npos &&
                rapid.value.code.find("ConfL \\Off;") != std::string::npos &&
                rapid.value.code.find("MoveJ pSafe01,vSafeCustom,fine,penqiang;") !=
                    std::string::npos &&
                rapid.value.code.find("MoveJ pPass01Start") != std::string::npos &&
                rapid.value.code.find("MoveL pPass01End") != std::string::npos &&
                rapid.value.code.find("PERS num nTableRPM") != std::string::npos &&
                rapid.value.code.find("PERS num nSprayTimes:=15") != std::string::npos &&
                rapid.value.code.find("FOR i FROM 1 TO nSprayTimes") != std::string::npos &&
                rapid.value.code.find("! Pass1, tilt=0.000 deg, D=50.000 mm") !=
                    std::string::npos,
                "RAPID output disables configuration monitoring and uses the approved motion pattern");
            expect(rapid.value.code.find("1000") != std::string::npos &&
                rapid.value.code.find("\\WObj:=wobj0") == std::string::npos &&
                rapid.value.code.find(
                    "[9E+09,90.6655,-0.000945636,9E+09,9E+09,9E+09]") !=
                    std::string::npos &&
                rapid.value.code.find(
                    "[9E+09,90.6666,-0.000918239,9E+09,9E+09,9E+09]") !=
                    std::string::npos,
                "RAPID poses use base coordinates, implicit wobj0 and fixed external axes");
            const std::size_t sprayOnce = rapid.value.code.find("PROC SprayOnce()");
            const std::size_t firstSafety = rapid.value.code.find("MoveJ pSafe01");
            expect(sprayOnce != std::string::npos && firstSafety > sprayOnce,
                "every repetition starts from the first safety point inside SprayOnce");
            expect(rapid.value.previewSteps.size() == 3 &&
                rapid.value.previewSteps[0].instruction == "MoveJ" &&
                rapid.value.previewSteps[0].sourceKind ==
                    rotationbody::RapidSequenceEntryKind::SafetyPoint &&
                rapid.value.previewSteps[1].instruction == "MoveJ" &&
                rapid.value.previewSteps[2].instruction == "MoveL",
                "RAPID translation exposes safety, trajectory start and end preview steps");
            expect(rapid.value.code.find("Return") == std::string::npos,
                "RAPID translation does not add an implicit return to each trajectory");
            expect(rapid.value.previewSteps[0].baseFromTool.translation().isApprox(
                    settings.safetyPositionBaseMeters) &&
                rapid.value.previewSteps[1].baseFromTool.translation().isApprox(
                    plan.baseFromPlanning *
                    generated.value.linearPoints.front().planningFromTool.translation()),
                "RAPID preview poses use base coordinates for viewport instruction stepping");
        }

        plan.safetyPositionBaseMeters = settings.safetyPositionBaseMeters;
        plan.safetySpeedMetersPerSecond = 0.05;
        plan.executionSequence = sequence;
        Eigen::Isometry3d initialBaseFromTool = Eigen::Isometry3d::Identity();
        initialBaseFromTool.translation() =
            settings.safetyPositionBaseMeters + Eigen::Vector3d(0.0, 0.0, 0.1);
        initialBaseFromTool.linear() =
            (plan.baseFromPlanning *
                generated.value.linearPoints.front().planningFromTool).linear();
        const auto timedExecution = rotationbody::ExecutionSequenceBuilder::build(
            plan,
            initialBaseFromTool);
        expect(timedExecution.ok() && timedExecution.value.size() >= 3,
            "execution timing accepts a safety-to-trajectory sequence");
        if(timedExecution && timedExecution.value.size() >= 3) {
            const Eigen::Vector3d trajectoryStart =
                (plan.baseFromPlanning *
                    generated.value.linearPoints.front().planningFromTool).translation();
            const double expectedTransferSeconds =
                (trajectoryStart - settings.safetyPositionBaseMeters).norm() /
                plan.safetySpeedMetersPerSecond;
            const double actualTransferSeconds =
                timedExecution.value[2].timeSeconds -
                timedExecution.value[1].timeSeconds;
            expect(std::abs(actualTransferSeconds - expectedTransferSeconds) < 1.0e-6,
                "execution timing derives safety transfer duration from distance and configured speed");
            expect(timedExecution.value.front().timeSeconds == 0.0 &&
                timedExecution.value.front().kind ==
                    rotationbody::TimedExecutionTargetKind::InitialPose,
                "execution timing preserves the current robot pose at time zero");
        }

        rotationbody::TrajectoryGroup returnGroup;
        rotationbody::TrajectoryPass firstPass;
        firstPass.id = "trajectory-1";
        firstPass.order = 1;
        firstPass.trajectory = generated.value;
        rotationbody::TrajectoryPass returnPass;
        returnPass.id = "trajectory-2";
        returnPass.order = 2;
        returnPass.trajectory = generated.value;
        expect(rotationbody::TrajectoryEditor::swapDirection(returnPass.trajectory).ok(),
            "reverse trajectory fixture can be generated");
        returnGroup.passes = { firstPass, returnPass };
        expect(rotationbody::TrajectoryGroupEditor::refreshSchedule(returnGroup).ok(),
            "reversed trajectory pair has a valid schedule");
        rotationbody::PublishedTrajectoryPlan returnPlan = plan;
        returnPlan.group = returnGroup;
        std::vector<rotationbody::RapidSequenceEntry> returnSequence(2);
        returnSequence[0].kind = rotationbody::RapidSequenceEntryKind::Trajectory;
        returnSequence[0].trajectoryPassId = firstPass.id;
        returnSequence[1].kind = rotationbody::RapidSequenceEntryKind::Trajectory;
        returnSequence[1].trajectoryPassId = returnPass.id;
        const auto returnRapid = rotationbody::RapidModuleGenerator::generate(
            returnPlan, settings, returnSequence);
        expect(returnRapid.ok(), "ABB translation accepts a coincident reverse pass");
        if(returnRapid && returnRapid.value.previewSteps.size() == 4) {
            const auto& forwardEnd = returnRapid.value.previewSteps[1].baseFromTool;
            const auto& reverseStart = returnRapid.value.previewSteps[2].baseFromTool;
            const auto& reverseEnd = returnRapid.value.previewSteps[3].baseFromTool;
            expect(forwardEnd.translation().isApprox(reverseStart.translation(), 1.0e-9) &&
                forwardEnd.linear().isApprox(reverseStart.linear(), 1.0e-9) &&
                forwardEnd.linear().isApprox(reverseEnd.linear(), 1.0e-9),
                "coincident reverse pass inherits the previous tool orientation");
        }

        // Projects saved before return-pass orientation continuity was fixed
        // contain X/Y axes flipped by 180 degrees around the spray Z axis.
        // Motion Planning must normalize that legacy representation too.
        rotationbody::TrajectoryPass legacyReturnPass = returnPass;
        for(auto& point : legacyReturnPass.trajectory.linearPoints) {
            point.planningFromTool.linear().col(0) *= -1.0;
            point.planningFromTool.linear().col(1) *= -1.0;
        }
        rotationbody::PublishedTrajectoryPlan legacyReturnPlan = returnPlan;
        legacyReturnPlan.group.passes = { firstPass, legacyReturnPass };
        legacyReturnPlan.executionSequence = returnSequence;
        legacyReturnPlan.safetyPositionBaseMeters = settings.safetyPositionBaseMeters;
        legacyReturnPlan.safetySpeedMetersPerSecond = 0.05;
        const Eigen::Isometry3d returnStart = legacyReturnPlan.baseFromPlanning *
            firstPass.trajectory.linearPoints.front().planningFromTool;
        const Eigen::Matrix3d expectedReturnOrientation =
            (legacyReturnPlan.baseFromPlanning *
                firstPass.trajectory.linearPoints.back().planningFromTool).linear();
        const auto legacyTimedReturn = rotationbody::ExecutionSequenceBuilder::build(
            legacyReturnPlan,
            returnStart);
        bool foundLegacyReturnTarget = false;
        bool legacyReturnOrientationContinuous = true;
        if(legacyTimedReturn) {
            for(const auto& target : legacyTimedReturn.value) {
                if(target.trajectoryPassId != legacyReturnPass.id) continue;
                foundLegacyReturnTarget = true;
                legacyReturnOrientationContinuous =
                    legacyReturnOrientationContinuous &&
                    target.baseFromTool.linear().isApprox(
                        expectedReturnOrientation, 1.0e-9);
            }
        }
        expect(legacyTimedReturn.ok() && foundLegacyReturnTarget &&
            legacyReturnOrientationContinuous,
            "motion execution removes the legacy 180-degree return-pass rotation");

        rotationbody::TrajectoryPass thirdPass = firstPass;
        thirdPass.id = "trajectory-3";
        thirdPass.order = 3;
        for(auto& point : thirdPass.trajectory.linearPoints) {
            point.planningFromTool.translation().z() -= 0.1;
        }
        rotationbody::TrajectoryPass fourthPass = thirdPass;
        fourthPass.id = "trajectory-4";
        fourthPass.order = 4;
        expect(rotationbody::TrajectoryEditor::swapDirection(fourthPass.trajectory).ok(),
            "second reverse trajectory fixture can be generated");
        returnPass.transitionAfterSeconds = 0.1;
        rotationbody::TrajectoryGroup completeGroup;
        completeGroup.passes = { firstPass, returnPass, thirdPass, fourthPass };
        expect(rotationbody::TrajectoryGroupEditor::refreshSchedule(completeGroup).ok(),
            "two reversed trajectory pairs have a valid schedule");
        rotationbody::PublishedTrajectoryPlan completePlan = plan;
        completePlan.group = completeGroup;
        std::vector<rotationbody::RapidSequenceEntry> completeSequence(6);
        completeSequence[0].kind = rotationbody::RapidSequenceEntryKind::SafetyPoint;
        completeSequence[1].kind = rotationbody::RapidSequenceEntryKind::Trajectory;
        completeSequence[1].trajectoryPassId = firstPass.id;
        completeSequence[2].kind = rotationbody::RapidSequenceEntryKind::Trajectory;
        completeSequence[2].trajectoryPassId = returnPass.id;
        completeSequence[3].kind = rotationbody::RapidSequenceEntryKind::SafetyPoint;
        completeSequence[4].kind = rotationbody::RapidSequenceEntryKind::Trajectory;
        completeSequence[4].trajectoryPassId = thirdPass.id;
        completeSequence[5].kind = rotationbody::RapidSequenceEntryKind::Trajectory;
        completeSequence[5].trajectoryPassId = fourthPass.id;
        const auto completeRapid = rotationbody::RapidModuleGenerator::generate(
            completePlan, settings, completeSequence);
        expect(completeRapid.ok(),
            "ABB translation accepts safety-1-2-safety-3-4 repetition order");
        if(completeRapid && completeRapid.value.previewSteps.size() == 10) {
            const std::string& code = completeRapid.value.code;
            const std::size_t safety1 = code.find("MoveJ pSafe01");
            const std::size_t pass1 = code.find("MoveJ pPass01Start");
            const std::size_t pass2 = code.find("MoveJ pPass02Start");
            const std::size_t safety2 = code.find("MoveJ pSafe02");
            const std::size_t pass3 = code.find("MoveJ pPass03Start");
            const std::size_t pass4 = code.find("MoveJ pPass04Start");
            expect(safety1 < pass1 && pass1 < pass2 && pass2 < safety2 &&
                safety2 < pass3 && pass3 < pass4,
                "SprayOnce preserves safety-1-2-safety-3-4 instruction order");
            expect(completeRapid.value.previewSteps[2].baseFromTool.linear().isApprox(
                    completeRapid.value.previewSteps[3].baseFromTool.linear(), 1.0e-9) &&
                completeRapid.value.previewSteps[2].baseFromTool.linear().isApprox(
                    completeRapid.value.previewSteps[4].baseFromTool.linear(), 1.0e-9) &&
                completeRapid.value.previewSteps[7].baseFromTool.linear().isApprox(
                    completeRapid.value.previewSteps[8].baseFromTool.linear(), 1.0e-9) &&
                completeRapid.value.previewSteps[7].baseFromTool.linear().isApprox(
                    completeRapid.value.previewSteps[9].baseFromTool.linear(), 1.0e-9),
                "both reversed passes return along the line without changing tool orientation");
        }
        returnPlan.safetyPositionBaseMeters = settings.safetyPositionBaseMeters;
        returnPlan.safetySpeedMetersPerSecond = 0.05;
        returnPlan.executionSequence.resize(3);
        returnPlan.executionSequence[0] = returnSequence[0];
        returnPlan.executionSequence[1].kind =
            rotationbody::RapidSequenceEntryKind::SafetyPoint;
        returnPlan.executionSequence[2] = returnSequence[1];
        const Eigen::Isometry3d returnInitial = returnPlan.baseFromPlanning *
            firstPass.trajectory.linearPoints.front().planningFromTool;
        const auto timedReturn = rotationbody::ExecutionSequenceBuilder::build(
            returnPlan,
            returnInitial);
        expect(timedReturn.ok(),
            "execution timing accepts a trajectory-safety-trajectory sequence");
        if(timedReturn) {
            const auto safetyTarget = std::find_if(
                timedReturn.value.begin(),
                timedReturn.value.end(),
                [](const auto& target) {
                    return target.kind ==
                        rotationbody::TimedExecutionTargetKind::SafetyPoint;
                });
            expect(safetyTarget != timedReturn.value.end() &&
                safetyTarget != timedReturn.value.begin() &&
                safetyTarget + 1 != timedReturn.value.end() &&
                safetyTarget->timeSeconds > (safetyTarget - 1)->timeSeconds &&
                (safetyTarget + 1)->timeSeconds > safetyTarget->timeSeconds,
                "middle safety transfers have positive durations on both sides");
        }

        rotationbody::TrajectoryPass independentPass = returnPass;
        independentPass.id = "trajectory-3";
        independentPass.trajectory = generated.value;
        independentPass.trajectory.linearPoints.front().planningFromTool.translation().x()
            += 0.02;
        independentPass.trajectory.linearPoints.back().planningFromTool.translation().x()
            += 0.02;
        const Eigen::Matrix3d independentRotation =
            Eigen::AngleAxisd(0.35, Eigen::Vector3d::UnitZ()).toRotationMatrix();
        for(auto& point : independentPass.trajectory.linearPoints) {
            point.planningFromTool.linear() =
                independentRotation * point.planningFromTool.linear();
        }
        firstPass.transitionAfterSeconds = 0.1;
        rotationbody::TrajectoryGroup independentGroup;
        independentGroup.passes = { firstPass, independentPass };
        expect(rotationbody::TrajectoryGroupEditor::refreshSchedule(independentGroup).ok(),
            "independent trajectory pair has a valid schedule");
        rotationbody::PublishedTrajectoryPlan independentPlan = plan;
        independentPlan.group = independentGroup;
        returnSequence[1].trajectoryPassId = independentPass.id;
        const auto independentRapid = rotationbody::RapidModuleGenerator::generate(
            independentPlan, settings, returnSequence);
        expect(independentRapid.ok(), "ABB translation accepts independent passes");
        if(independentRapid && independentRapid.value.previewSteps.size() == 4) {
            expect(!independentRapid.value.previewSteps[1].baseFromTool.linear().isApprox(
                    independentRapid.value.previewSteps[2].baseFromTool.linear(), 1.0e-9),
                "independent pass keeps its own orientation");
        }
    }
}

int main(int argc, char** argv)
{
    if(argc > 1 && std::string(argv[1]) == "--region-only")
    {
        testAutomaticRegionRecognition();
        if(failureCount != 0)
        {
            std::cerr << failureCount << " region regression checks failed.\n";
            return EXIT_FAILURE;
        }
        std::cout << "Region recognition regression checks passed.\n";
        return EXIT_SUCCESS;
    }
    if(argc > 1 && std::string(argv[1]) == "--calibration-only")
    {
        testWorkpieceCalibration();
        if(failureCount != 0)
        {
            std::cerr << failureCount << " calibration regression checks failed.\n";
            return EXIT_FAILURE;
        }
        std::cout << "Workpiece calibration regression checks passed.\n";
        return EXIT_SUCCESS;
    }
    if(argc > 1 && std::string(argv[1]) == "--trajectory-only")
    {
        testTrajectoryPlanningAndEditing();
        testTrajectoryGroupAndRapidTranslation();
        if (failureCount != 0)
        {
            std::cerr << failureCount << " trajectory regression checks failed.\n";
            return EXIT_FAILURE;
        }
        std::cout << "Trajectory planning regression checks passed.\n";
        return EXIT_SUCCESS;
    }

    testTriangleMeshContract();
    testTransformContract();
    testCompletePartAlignment();
    testSimulationBlockPlacementAndTransformOperations();
    testYzSectioning();
    testAutomaticRegionRecognition();
    testRegionEditHistory();
    testRegionEditScaleInvariance();
    testSprayBoundaries();
    testWorkpieceCalibration();
    testTrajectoryPlanningAndEditing();
    testTrajectoryGroupAndRapidTranslation();

    if (failureCount != 0)
    {
        std::cerr << failureCount << " regression checks failed.\n";
        return EXIT_FAILURE;
    }
    std::cout << "RotationBodyTrajectoryPlanning regression checks passed.\n";
    return EXIT_SUCCESS;
}
