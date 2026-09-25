#include <RotationBodyTrajectoryPlanning/TrajectoryPlanning/TrajectoryParameterTextParser.h>

#include <RotationBodyTrajectoryPlanning/TrajectoryPlanning/TrajectoryPlanner.h>

#include <cmath>
#include <exception>
#include <locale>
#include <sstream>
#include <string>
#include <utility>

namespace smrobot::spray::rotationbody
{
    PlanningResult<std::vector<TrajectoryGenerationParameters>>
    TrajectoryParameterTextParser::parse(std::string_view text)
    {
        std::istringstream stream{ std::string(text) };
        stream.imbue(std::locale::classic());
        std::vector<double> values;
        std::string token;
        while(stream >> token) {
            std::size_t consumed = 0;
            double value = 0.0;
            try {
                value = std::stod(token, &consumed);
            } catch(const std::exception&) {
                return PlanningResult<std::vector<TrajectoryGenerationParameters>>::failure(
                    PlanningErrorCode::InvalidArgument,
                    "The trajectory parameter TXT contains a non-numeric token.");
            }
            if(consumed != token.size() || !std::isfinite(value)) {
                return PlanningResult<std::vector<TrajectoryGenerationParameters>>::failure(
                    PlanningErrorCode::InvalidArgument,
                    "The trajectory parameter TXT contains an invalid numeric value.");
            }
            values.push_back(value);
        }
        if(values.empty() || values.size() % 6 != 0) {
            return PlanningResult<std::vector<TrajectoryGenerationParameters>>::failure(
                PlanningErrorCode::InvalidArgument,
                "The trajectory parameter TXT must contain six numbers per trajectory.");
        }

        constexpr double pi = 3.14159265358979323846;
        std::vector<TrajectoryGenerationParameters> result;
        result.reserve(values.size() / 6);
        for(std::size_t offset = 0; offset < values.size(); offset += 6) {
            TrajectoryGenerationParameters parameters;
            parameters.sprayDistanceMeters = values[offset] / 1000.0;
            parameters.tiltRadians = values[offset + 1] * pi / 180.0;
            parameters.speedMetersPerSecond = values[offset + 2] / 1000.0;
            parameters.startExtensionMeters = values[offset + 3] / 1000.0;
            parameters.endExtensionMeters = values[offset + 4] / 1000.0;
            parameters.pointCount = 0;
            parameters.positionerRpm = values[offset + 5];
            parameters.reversed = false;
            if(parameters.sprayDistanceMeters <
                    TrajectoryPlanner::minimumSprayDistanceMeters ||
                parameters.sprayDistanceMeters >
                    TrajectoryPlanner::maximumSprayDistanceMeters ||
                std::abs(parameters.tiltRadians) >
                    TrajectoryPlanner::maximumAbsoluteTiltRadians ||
                parameters.speedMetersPerSecond <= 0.0 ||
                parameters.startExtensionMeters < 0.0 ||
                parameters.endExtensionMeters < 0.0) {
                return PlanningResult<std::vector<TrajectoryGenerationParameters>>::failure(
                    PlanningErrorCode::InvalidArgument,
                    "A trajectory parameter set is outside the supported ranges.");
            }
            result.push_back(parameters);
        }
        return PlanningResult<std::vector<TrajectoryGenerationParameters>>::success(
            std::move(result));
    }
}
