cmake_minimum_required(VERSION 3.20)

add_executable(RotationBodyTrajectoryPlanningTest main.cpp)
target_link_libraries(RotationBodyTrajectoryPlanningTest
    PRIVATE SMRobotSpray::RotationBodyTrajectoryPlanning
)
target_compile_features(RotationBodyTrajectoryPlanningTest PRIVATE cxx_std_17)
set_target_properties(RotationBodyTrajectoryPlanningTest PROPERTIES
    FOLDER tests/SMRobotSpray
)
add_test(
    NAME RotationBodyTrajectoryPlanningTest
    COMMAND $<TARGET_FILE:RotationBodyTrajectoryPlanningTest>
)
