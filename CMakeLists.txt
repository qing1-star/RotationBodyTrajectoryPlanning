#############################################################
# CMake requirement (From Modern CMake)
#############################################################
cmake_minimum_required(VERSION 3.20)

get_filename_component( TARGET_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME )
message( STATUS "\n-------------  Config ${PACKAGE_NAME}::${TARGET_NAME}: ${CMAKE_CURRENT_SOURCE_DIR}  -------------" )

set(SOURCES
    src/Core/TriangleMesh.cpp
    src/Core/TransformUtils.cpp
    src/Alignment/RotationBodyAlignmentSolver.cpp
    src/Alignment/SimulationBlockPlacementSolver.cpp
    src/Alignment/ModelTransformOperations.cpp
    src/Sectioning/ContourTopology.cpp
    src/Sectioning/YzSectionExtractor.cpp
    src/RegionPlanning/ToothRegionRecognizer.cpp
    src/RegionPlanning/RegionEditHistory.cpp
    src/RegionPlanning/SprayBoundaryBuilder.cpp
    src/TrajectoryPlanning/TrajectoryPlanner.cpp
    src/TrajectoryPlanning/TrajectoryEditor.cpp
    src/TrajectoryPlanning/RelativeHelicalTrajectoryBuilder.cpp
    src/TrajectoryPlanning/TrajectoryGroupEditor.cpp
    src/TrajectoryPlanning/ExecutionSequenceBuilder.cpp
    src/ABBTranslation/RapidModuleGenerator.cpp
    src/Calibration/WorkpieceCalibration.cpp
    src/Persistence/PublishedTrajectoryPlanStore.cpp
)

set(INTERFACE_HEADERS
    include/RotationBodyTrajectoryPlanning/Core/PlanningTypes.h
    include/RotationBodyTrajectoryPlanning/Core/TriangleMesh.h
    include/RotationBodyTrajectoryPlanning/Core/TransformUtils.h
    include/RotationBodyTrajectoryPlanning/Alignment/RotationBodyAlignmentSolver.h
    include/RotationBodyTrajectoryPlanning/Alignment/SimulationBlockPlacementSolver.h
    include/RotationBodyTrajectoryPlanning/Alignment/ModelTransformOperations.h
    include/RotationBodyTrajectoryPlanning/Sectioning/ContourTopology.h
    include/RotationBodyTrajectoryPlanning/Sectioning/YzSectionExtractor.h
    include/RotationBodyTrajectoryPlanning/RegionPlanning/ToothRegionRecognizer.h
    include/RotationBodyTrajectoryPlanning/RegionPlanning/RegionEditHistory.h
    include/RotationBodyTrajectoryPlanning/RegionPlanning/SprayBoundaryBuilder.h
    include/RotationBodyTrajectoryPlanning/TrajectoryPlanning/TrajectoryTypes.h
    include/RotationBodyTrajectoryPlanning/TrajectoryPlanning/TrajectoryPlanner.h
    include/RotationBodyTrajectoryPlanning/TrajectoryPlanning/TrajectoryEditor.h
    include/RotationBodyTrajectoryPlanning/TrajectoryPlanning/RelativeHelicalTrajectoryBuilder.h
    include/RotationBodyTrajectoryPlanning/TrajectoryPlanning/TrajectoryGroupEditor.h
    include/RotationBodyTrajectoryPlanning/TrajectoryPlanning/ExecutionSequenceBuilder.h
    include/RotationBodyTrajectoryPlanning/ABBTranslation/RapidModuleGenerator.h
    include/RotationBodyTrajectoryPlanning/Calibration/WorkpieceCalibration.h
    include/RotationBodyTrajectoryPlanning/Persistence/PublishedTrajectoryPlanStore.h
)

cmake_path(SET abs_include_dir NORMALIZE "${CMAKE_CURRENT_SOURCE_DIR}/include")

add_library(${TARGET_NAME} STATIC ${SOURCES} ${INTERFACE_HEADERS})
add_library(${PACKAGE_NAME}::${TARGET_NAME} ALIAS ${TARGET_NAME})

target_compile_features(${TARGET_NAME} PUBLIC cxx_std_17)
target_include_directories(${TARGET_NAME}
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${PACKAGE_NAME}/include>
    PRIVATE
        $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
)

target_link_libraries(${TARGET_NAME}
    PUBLIC ${${PACKAGE_NAME}${TARGET_NAME}_RequiredLibsPublic}
    PRIVATE ${${PACKAGE_NAME}${TARGET_NAME}_RequiredLibsPrivate}
)

set_target_properties(${TARGET_NAME} PROPERTIES
    FOLDER ${PACKAGE_NAME}
    DEBUG_POSTFIX _shared_${CMAKE_DEBUG_POSTFIX}
    RELEASE_POSTFIX _shared_${CMAKE_RELEASE_POSTFIX}
)

if(BuildTests)
    add_subdirectory(regression/RotationBodyTrajectoryPlanningTest)
endif()

if(ExportPacakge${PACKAGE_NAME}_${TARGET_NAME})
    function_InstallTarget(
        TARGET_NAME ${TARGET_NAME}
        PACKAGE_NAME ${PACKAGE_NAME}
        PACKAGE_INTERFACE_HEADERS ${abs_include_dir}
    )
endif()

message( STATUS "-------------  Config ${PACKAGE_NAME}::${TARGET_NAME}: ${CMAKE_CURRENT_SOURCE_DIR}  -------------\n" )
