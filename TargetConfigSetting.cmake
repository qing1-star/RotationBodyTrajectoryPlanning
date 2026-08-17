##########################################################################
# CMake requirement
##########################################################################
cmake_minimum_required( VERSION 3.20 )

set( ${TARGET_NAME}_RequiredLibsPublic
    Eigen3::Eigen
)

set( ${TARGET_NAME}_RequiredLibsPrivate
    nlohmann_json::nlohmann_json
)
