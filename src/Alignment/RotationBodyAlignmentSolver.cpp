#include <RotationBodyTrajectoryPlanning/Alignment/RotationBodyAlignmentSolver.h>

#include <RotationBodyTrajectoryPlanning/Core/TransformUtils.h>

#include <Eigen/Eigenvalues>

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <limits>
#include <numeric>
#include <sstream>
#include <vector>

namespace smrobot::spray::rotationbody
{
    namespace
    {
        constexpr int angularSectorCount = 24;
        constexpr double geometryEpsilon = 1.0e-12;

        struct ProjectedPoint
        {
            double axial{ 0.0 };
            double x{ 0.0 };
            double y{ 0.0 };
        };

        struct ExtentBin
        {
            std::size_t count{ 0 };
            double minimumX{ std::numeric_limits<double>::infinity() };
            double maximumX{ -std::numeric_limits<double>::infinity() };
            double minimumY{ std::numeric_limits<double>::infinity() };
            double maximumY{ -std::numeric_limits<double>::infinity() };
        };

        struct RadiusBin
        {
            std::size_t count{ 0 };
            double axialSum{ 0.0 };
            std::vector<double> radii;
            std::array<double, angularSectorCount> sectorMaximum{};

            RadiusBin()
            {
                sectorMaximum.fill(-1.0);
            }
        };

        struct ProfileSample
        {
            double axial{ 0.0 };
            double outerRadius{ 0.0 };
            double smoothedRadius{ 0.0 };
        };

        struct AxisCandidate
        {
            bool valid{ false };
            int pcaIndex{ -1 };
            Eigen::Vector3d axis = Eigen::Vector3d::UnitZ();
            Eigen::Vector3d axisPoint = Eigen::Vector3d::Zero();
            double minimumAxial{ 0.0 };
            double maximumAxial{ 0.0 };
            double minimumDiameterAxial{ 0.0 };
            double minimumOuterRadius{ 0.0 };
            double maximumOuterRadius{ 0.0 };
            double profileContrast{ 0.0 };
            double angularCoverage{ 0.0 };
            double circularityError{ 1.0 };
            double centerDrift{ 1.0 };
            double score{ std::numeric_limits<double>::infinity() };
        };

        double percentile(std::vector<double> values, double fraction)
        {
            if (values.empty())
            {
                return 0.0;
            }
            fraction = std::clamp(fraction, 0.0, 1.0);
            const std::size_t index = static_cast<std::size_t>(
                std::floor(fraction * static_cast<double>(values.size() - 1)));
            std::nth_element(values.begin(), values.begin() + index, values.end());
            return values[index];
        }

        double median(std::vector<double> values)
        {
            return percentile(std::move(values), 0.5);
        }

        int axialBinIndex(double axial, double minimum, double extent, int binCount)
        {
            if (extent <= geometryEpsilon)
            {
                return 0;
            }
            const double normalized = std::clamp((axial - minimum) / extent, 0.0, 1.0);
            return std::min(binCount - 1, static_cast<int>(normalized * binCount));
        }

        Eigen::Vector3d deterministicDirection(Eigen::Vector3d axis)
        {
            Eigen::Index dominant = 0;
            axis.cwiseAbs().maxCoeff(&dominant);
            if (axis[dominant] < 0.0)
            {
                axis = -axis;
            }
            return axis.normalized();
        }

        AxisCandidate evaluateAxis(
            const TriangleMesh& mesh,
            const Eigen::Vector3d& centroid,
            const Eigen::Vector3d& candidateAxis,
            int pcaIndex,
            double covariancePairMismatch)
        {
            AxisCandidate result;
            result.pcaIndex = pcaIndex;
            result.axis = deterministicDirection(candidateAxis);

            Eigen::Vector3d basisX = result.axis.unitOrthogonal().normalized();
            Eigen::Vector3d basisY = result.axis.cross(basisX).normalized();

            std::vector<ProjectedPoint> projected;
            projected.reserve(mesh.positions.size());
            result.minimumAxial = std::numeric_limits<double>::infinity();
            result.maximumAxial = -std::numeric_limits<double>::infinity();
            for (const Eigen::Vector3d& position : mesh.positions)
            {
                const Eigen::Vector3d centered = position - centroid;
                ProjectedPoint point;
                point.axial = centered.dot(result.axis);
                point.x = centered.dot(basisX);
                point.y = centered.dot(basisY);
                projected.push_back(point);
                result.minimumAxial = std::min(result.minimumAxial, point.axial);
                result.maximumAxial = std::max(result.maximumAxial, point.axial);
            }

            const double axialExtent = result.maximumAxial - result.minimumAxial;
            const double meshScale = mesh.scale();
            if (!std::isfinite(axialExtent) || !std::isfinite(meshScale) ||
                axialExtent <= meshScale * 1.0e-9)
            {
                return result;
            }

            const int binCount = std::clamp(
                static_cast<int>(std::sqrt(static_cast<double>(projected.size()))), 16, 128);
            std::vector<ExtentBin> extentBins(static_cast<std::size_t>(binCount));
            for (const ProjectedPoint& point : projected)
            {
                ExtentBin& bin = extentBins[static_cast<std::size_t>(axialBinIndex(
                    point.axial, result.minimumAxial, axialExtent, binCount))];
                ++bin.count;
                bin.minimumX = std::min(bin.minimumX, point.x);
                bin.maximumX = std::max(bin.maximumX, point.x);
                bin.minimumY = std::min(bin.minimumY, point.y);
                bin.maximumY = std::max(bin.maximumY, point.y);
            }

            std::vector<double> centerXs;
            std::vector<double> centerYs;
            const double usefulSpan = std::max(meshScale * 1.0e-7, geometryEpsilon);
            for (const ExtentBin& bin : extentBins)
            {
                if (bin.count >= 4 && bin.maximumX - bin.minimumX > usefulSpan &&
                    bin.maximumY - bin.minimumY > usefulSpan)
                {
                    centerXs.push_back(0.5 * (bin.minimumX + bin.maximumX));
                    centerYs.push_back(0.5 * (bin.minimumY + bin.maximumY));
                }
            }
            const double centerX = centerXs.empty() ? 0.0 : median(centerXs);
            const double centerY = centerYs.empty() ? 0.0 : median(centerYs);
            result.axisPoint = centroid + basisX * centerX + basisY * centerY;

            std::vector<double> centerOffsets;
            centerOffsets.reserve(centerXs.size());
            for (std::size_t i = 0; i < centerXs.size(); ++i)
            {
                centerOffsets.push_back(std::hypot(centerXs[i] - centerX, centerYs[i] - centerY));
            }

            std::vector<RadiusBin> radiusBins(static_cast<std::size_t>(binCount));
            double maximumRadius = 0.0;
            for (ProjectedPoint& point : projected)
            {
                point.x -= centerX;
                point.y -= centerY;
                const double radius = std::hypot(point.x, point.y);
                maximumRadius = std::max(maximumRadius, radius);
                RadiusBin& bin = radiusBins[static_cast<std::size_t>(axialBinIndex(
                    point.axial, result.minimumAxial, axialExtent, binCount))];
                ++bin.count;
                bin.axialSum += point.axial;
                bin.radii.push_back(radius);
                double angle = std::atan2(point.y, point.x);
                if (angle < 0.0)
                {
                    angle += 2.0 * pi;
                }
                const int sector = std::min(
                    angularSectorCount - 1,
                    static_cast<int>(angle * angularSectorCount / (2.0 * pi)));
                bin.sectorMaximum[static_cast<std::size_t>(sector)] = std::max(
                    bin.sectorMaximum[static_cast<std::size_t>(sector)], radius);
            }
            if (!std::isfinite(maximumRadius) || maximumRadius <= usefulSpan)
            {
                return result;
            }

            std::vector<double> circularityErrors;
            std::vector<double> angularCoverages;
            std::vector<ProfileSample> profile;
            std::size_t nonEmptyBinCount = 0;
            std::size_t qualifiedBinCount = 0;
            const std::size_t minimumBinPopulation = std::max<std::size_t>(
                4, projected.size() / static_cast<std::size_t>(binCount * 100));
            for (const RadiusBin& bin : radiusBins)
            {
                if (bin.count == 0)
                {
                    continue;
                }
                ++nonEmptyBinCount;
                if (bin.count < minimumBinPopulation || bin.radii.size() < 4)
                {
                    continue;
                }

                std::vector<double> sectorRadii;
                for (double sectorRadius : bin.sectorMaximum)
                {
                    if (sectorRadius >= 0.0)
                    {
                        sectorRadii.push_back(sectorRadius);
                    }
                }
                const double coverage = static_cast<double>(sectorRadii.size()) /
                    static_cast<double>(angularSectorCount);
                if (sectorRadii.size() >= 3)
                {
                    const double sectorMedian = median(sectorRadii);
                    std::vector<double> deviations;
                    deviations.reserve(sectorRadii.size());
                    for (double radius : sectorRadii)
                    {
                        deviations.push_back(std::abs(radius - sectorMedian));
                    }
                    const double robustDeviation = 1.4826 * median(std::move(deviations));
                    circularityErrors.push_back(robustDeviation /
                        std::max(sectorMedian, usefulSpan));
                    angularCoverages.push_back(coverage);
                    ++qualifiedBinCount;
                }

                ProfileSample sample;
                sample.axial = bin.axialSum / static_cast<double>(bin.count);
                sample.outerRadius = percentile(bin.radii, 0.92);
                sample.smoothedRadius = sample.outerRadius;
                if (sample.outerRadius > usefulSpan)
                {
                    profile.push_back(sample);
                }
            }
            if (profile.empty() || qualifiedBinCount == 0)
            {
                return result;
            }

            std::sort(profile.begin(), profile.end(), [](const ProfileSample& first, const ProfileSample& second)
            {
                return first.axial < second.axial;
            });
            for (std::size_t i = 0; i < profile.size(); ++i)
            {
                const std::size_t begin = i > 2 ? i - 2 : 0;
                const std::size_t end = std::min(profile.size(), i + 3);
                std::vector<double> neighborhood;
                neighborhood.reserve(end - begin);
                for (std::size_t j = begin; j < end; ++j)
                {
                    neighborhood.push_back(profile[j].outerRadius);
                }
                profile[i].smoothedRadius = median(std::move(neighborhood));
            }

            std::size_t searchBegin = 0;
            std::size_t searchEnd = profile.size();
            if (profile.size() >= 7)
            {
                searchBegin = 1;
                searchEnd = profile.size() - 1;
            }
            result.minimumOuterRadius = std::numeric_limits<double>::infinity();
            result.maximumOuterRadius = 0.0;
            for (std::size_t i = searchBegin; i < searchEnd; ++i)
            {
                result.minimumOuterRadius = std::min(
                    result.minimumOuterRadius, profile[i].smoothedRadius);
                result.maximumOuterRadius = std::max(
                    result.maximumOuterRadius, profile[i].smoothedRadius);
            }
            if (!std::isfinite(result.minimumOuterRadius) || result.maximumOuterRadius <= usefulSpan)
            {
                return result;
            }

            const double nearMinimumTolerance = std::max(
                usefulSpan,
                0.015 * (result.maximumOuterRadius - result.minimumOuterRadius));
            double minimumStationSum = 0.0;
            std::size_t minimumStationCount = 0;
            for (std::size_t i = searchBegin; i < searchEnd; ++i)
            {
                if (profile[i].smoothedRadius <= result.minimumOuterRadius + nearMinimumTolerance)
                {
                    minimumStationSum += profile[i].axial;
                    ++minimumStationCount;
                }
            }
            result.minimumDiameterAxial = minimumStationCount > 0
                ? minimumStationSum / static_cast<double>(minimumStationCount)
                : 0.5 * (result.minimumAxial + result.maximumAxial);

            result.circularityError = median(std::move(circularityErrors));
            result.angularCoverage = median(std::move(angularCoverages));
            result.centerDrift = centerOffsets.empty()
                ? 0.0
                : median(std::move(centerOffsets)) / std::max(result.maximumOuterRadius, usefulSpan);
            const double binCoverage = nonEmptyBinCount == 0
                ? 0.0
                : static_cast<double>(qualifiedBinCount) / static_cast<double>(nonEmptyBinCount);
            result.profileContrast = std::clamp(
                (result.maximumOuterRadius - result.minimumOuterRadius) /
                    result.maximumOuterRadius,
                0.0,
                1.0);
            result.score =
                0.45 * std::min(result.circularityError, 2.0) +
                0.20 * (1.0 - result.angularCoverage) +
                0.15 * std::min(result.centerDrift, 2.0) +
                0.10 * (1.0 - std::clamp(binCoverage, 0.0, 1.0)) +
                0.10 * std::min(covariancePairMismatch, 2.0);
            result.valid = std::isfinite(result.score);
            return result;
        }

        std::string candidateDiagnostic(const AxisCandidate& candidate)
        {
            std::ostringstream stream;
            stream << std::fixed << std::setprecision(4)
                   << "PCA axis " << candidate.pcaIndex
                   << ": score=" << candidate.score
                   << ", circularity=" << candidate.circularityError
                   << ", angularCoverage=" << candidate.angularCoverage
                   << ", centerDrift=" << candidate.centerDrift;
            return stream.str();
        }
    }

    PlanningResult<AlignmentResult> RotationBodyAlignmentSolver::solve(const TriangleMesh& mesh)
    {
        const PlanningResult<void> validation = mesh.validate();
        if (!validation)
        {
            return PlanningResult<AlignmentResult>::failure(
                validation.error.code, validation.error.message);
        }
        if (mesh.positions.size() < 6)
        {
            return PlanningResult<AlignmentResult>::failure(
                PlanningErrorCode::AlignmentFailed,
                "At least six finite mesh positions are required to estimate a rotary axis.");
        }

        Eigen::Vector3d centroid = Eigen::Vector3d::Zero();
        for (const Eigen::Vector3d& position : mesh.positions)
        {
            centroid += position;
        }
        centroid /= static_cast<double>(mesh.positions.size());

        Eigen::Matrix3d covariance = Eigen::Matrix3d::Zero();
        for (const Eigen::Vector3d& position : mesh.positions)
        {
            const Eigen::Vector3d centered = position - centroid;
            covariance.noalias() += centered * centered.transpose();
        }
        covariance /= static_cast<double>(mesh.positions.size());

        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigenSolver(covariance);
        if (eigenSolver.info() != Eigen::Success || !eigenSolver.eigenvalues().allFinite() ||
            !eigenSolver.eigenvectors().allFinite())
        {
            return PlanningResult<AlignmentResult>::failure(
                PlanningErrorCode::AlignmentFailed,
                "The covariance eigensystem could not be solved for this mesh.");
        }

        std::vector<AxisCandidate> candidates;
        candidates.reserve(3);
        const Eigen::Vector3d eigenvalues = eigenSolver.eigenvalues();
        for (int candidateIndex = 0; candidateIndex < 3; ++candidateIndex)
        {
            int firstOther = 0;
            int secondOther = 1;
            if (candidateIndex == 0)
            {
                firstOther = 1;
                secondOther = 2;
            }
            else if (candidateIndex == 1)
            {
                firstOther = 0;
                secondOther = 2;
            }
            const double pairScale = std::max({
                std::abs(eigenvalues[firstOther]),
                std::abs(eigenvalues[secondOther]),
                geometryEpsilon });
            const double pairMismatch = std::abs(
                eigenvalues[firstOther] - eigenvalues[secondOther]) / pairScale;
            AxisCandidate candidate = evaluateAxis(
                mesh,
                centroid,
                eigenSolver.eigenvectors().col(candidateIndex),
                candidateIndex,
                pairMismatch);
            if (candidate.valid)
            {
                candidates.push_back(std::move(candidate));
            }
        }
        if (candidates.empty())
        {
            return PlanningResult<AlignmentResult>::failure(
                PlanningErrorCode::AlignmentFailed,
                "No PCA direction produced a measurable rotary profile.");
        }

        std::sort(candidates.begin(), candidates.end(), [](const AxisCandidate& first, const AxisCandidate& second)
        {
            if (first.score != second.score)
            {
                return first.score < second.score;
            }
            return first.pcaIndex < second.pcaIndex;
        });
        const AxisCandidate& best = candidates.front();

        const double distanceToLow = std::abs(best.minimumDiameterAxial - best.minimumAxial);
        const double distanceToHigh = std::abs(best.maximumAxial - best.minimumDiameterAxial);
        const bool highEndIsBottom = distanceToHigh > distanceToLow;
        const double bottomAxial = highEndIsBottom ? best.maximumAxial : best.minimumAxial;
        const Eigen::Vector3d directedAxis = highEndIsBottom ? -best.axis : best.axis;
        const Eigen::Vector3d bottomAxisCenter = best.axisPoint + best.axis * bottomAxial;

        Eigen::Isometry3d planningFromMesh = Eigen::Isometry3d::Identity();
        planningFromMesh.linear() = rotationAligningVectorToVector(
            directedAxis, Eigen::Vector3d::UnitZ());
        planningFromMesh.translation() = -planningFromMesh.linear() * bottomAxisCenter;
        if (!isFiniteTransform(planningFromMesh))
        {
            return PlanningResult<AlignmentResult>::failure(
                PlanningErrorCode::AlignmentFailed,
                "Axis alignment produced a non-finite rigid transform.");
        }

        const PlanningResult<Eigen::AlignedBox3d> bounds = mesh.bounds(planningFromMesh);
        if (!bounds)
        {
            return PlanningResult<AlignmentResult>::failure(
                bounds.error.code, bounds.error.message);
        }
        double maximumRadius = 0.0;
        for (const Eigen::Vector3d& position : mesh.positions)
        {
            const Eigen::Vector3d planned = planningFromMesh * position;
            maximumRadius = std::max(maximumRadius, std::hypot(planned.x(), planned.y()));
        }

        const double symmetryQuality = std::exp(-3.0 * std::max(0.0, best.score));
        double candidateSeparation = 1.0;
        if (candidates.size() > 1)
        {
            candidateSeparation = std::clamp(
                (candidates[1].score - best.score) /
                    std::max(candidates[1].score, 0.05),
                0.0,
                1.0);
        }
        const double endSelectionQuality = std::clamp(best.profileContrast / 0.08, 0.0, 1.0);
        const double confidence = std::clamp(
            symmetryQuality *
                (0.20 + 0.80 * candidateSeparation) *
                (0.35 + 0.65 * endSelectionQuality),
            0.0,
            1.0);

        AlignmentResult alignment;
        alignment.planningFromMesh = planningFromMesh;
        alignment.automaticBaseline = planningFromMesh;
        alignment.bottomAxisCenterInMesh = bottomAxisCenter;
        alignment.lowConfidence = confidence < 0.55;
        alignment.statistics.vertexCount = mesh.vertexCount();
        alignment.statistics.triangleCount = mesh.triangleCount();
        alignment.statistics.boundsMinimum = bounds.value.min();
        alignment.statistics.boundsMaximum = bounds.value.max();
        alignment.statistics.estimatedAxisInMesh = directedAxis;
        alignment.statistics.heightMeters = std::max(
            0.0, bounds.value.max().z() - bounds.value.min().z());
        alignment.statistics.maximumDiameterMeters = 2.0 * maximumRadius;
        alignment.statistics.minimumDiameterMeters = 2.0 * best.minimumOuterRadius;
        alignment.statistics.axisConfidence = confidence;

        PlanningResult<AlignmentResult> output =
            PlanningResult<AlignmentResult>::success(std::move(alignment));
        for (const AxisCandidate& candidate : candidates)
        {
            output.diagnostics.push_back(candidateDiagnostic(candidate));
        }
        std::ostringstream stationDiagnostic;
        stationDiagnostic << std::fixed << std::setprecision(6)
                          << "Minimum-diameter station=" << best.minimumDiameterAxial
                          << " m; selected " << (highEndIsBottom ? "high" : "low")
                          << " axial end as planning bottom; confidence=" << confidence << '.';
        output.diagnostics.push_back(stationDiagnostic.str());
        if (output.value.lowConfidence)
        {
            output.diagnostics.push_back(
                "Alignment is usable but ambiguous; verify the rotary axis and bottom before confirming the workpiece frame.");
        }
        return output;
    }
}
