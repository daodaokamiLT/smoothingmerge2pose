project(smoothingMergePose)

cmake_minimum_required(VERSION 3.2)

set(CMAKE_BUILD_TYPE Release)

find_package(Ceres REQUIRED)
include_directories(${Ceres_INCLUDE_DIRS}
    "/usr/include/eigen3"
    ${OpenCV_INCLUDE_DIRS})

add_library(svvfusion src/core/svvfusion.cpp)
add_library(vinsfusion src/vins_fusion/vinsfusion.cpp)

target_link_libraries(svvfusion ${CERES_LIBRARIES})
target_link_libraries(vinsfusion ${CERES_LIBRARIES})

add_executable(test_svvfusion test_fusion.cpp)
target_link_libraries(test_svvfusion svvfusion)

add_executable(test_vinsfusion test_vinsfusion.cpp)
target_link_libraries(test_vinsfusion vinsfusion)