cmake_minimum_required(VERSION 3.31)
project(MisalignmentAngleEstimation)

set(CMAKE_CXX_STANDARD 20)

# 添加源文件
file(GLOB_RECURSE SRC_FILES
        ${CMAKE_SOURCE_DIR}/src/*.cpp
)

add_executable(${PROJECT_NAME} ${SRC_FILES})

target_compile_definitions(
        ${PROJECT_NAME} PRIVATE
        _USE_MATH_DEFINES
)

target_include_directories(
        ${PROJECT_NAME} PRIVATE
        ${CMAKE_SOURCE_DIR}/src
)

# 将库文件链接到目标文件
target_link_libraries(
        ${PROJECT_NAME} PUBLIC
        range-v3::range-v3
        fmt::fmt
        spdlog::spdlog
        Eigen3::Eigen
        #        nlohmann_json::nlohmann_json
        #        tomlplusplus::tomlplusplus
        #        GeographicLib::GeographicLib
)

# conan 管理的三方库
find_package(range-v3 REQUIRED)
find_package(fmt REQUIRED)
find_package(spdlog REQUIRED)
find_package(Eigen3 REQUIRED)
#find_package(nlohmann_json REQUIRED)
#find_package(tomlplusplus REQUIRED)
#find_package(geographiclib REQUIRED)
