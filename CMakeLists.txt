cmake_minimum_required(VERSION 3.16)
project(ProjectHelium LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(Qt6 REQUIRED COMPONENTS Widgets DataVisualization)

add_executable(scatter Numerical_calculation.cpp)
target_link_libraries(scatter PRIVATE Qt6::Widgets Qt6::DataVisualization)