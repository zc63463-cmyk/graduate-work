#
# CMakeLists.txt for GMRES implementation
# ================================
#
# Based on Saad & Schultz (1986)
# Date: 2026-04-25
#

cmake_minimum_required(VERSION 3.10)
project(GMRES_Implementation)

# Set C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add executable
add_executable(gmres_example main.cpp)

# Set output directory
set_target_properties(gmres_example
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

# Print configuration
message(STATUS "Project: ${PROJECT_NAME}")
message(STATUS "C++ Standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "Output: ${CMAKE_BINARY_DIR}/bin/gmres_example")
