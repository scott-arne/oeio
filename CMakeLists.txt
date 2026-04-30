cmake_minimum_required(VERSION 3.16)

if(POLICY CMP0144)
    cmake_policy(SET CMP0144 NEW)
endif()

project(oeio VERSION 0.2.2 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Options
option(OEIO_BUILD_TESTS "Build tests" ON)
option(OEIO_BUILD_PYTHON "Build Python bindings" ON)
option(OEIO_UNIVERSAL2 "Build universal2 binary for macOS" OFF)
option(OEIO_USE_STABLE_ABI "Use Python stable ABI" ON)
option(OEIO_BUILD_BENCHMARKS "Build benchmarks" OFF)

if(APPLE AND OEIO_UNIVERSAL2)
    set(CMAKE_OSX_ARCHITECTURES "x86_64;arm64" CACHE STRING "Build architectures" FORCE)
elseif(APPLE)
    unset(CMAKE_OSX_ARCHITECTURES CACHE)
endif()

if(APPLE AND NOT CMAKE_OSX_DEPLOYMENT_TARGET)
    set(CMAKE_OSX_DEPLOYMENT_TARGET "11.0" CACHE STRING "Minimum macOS version")
endif()

# Dependencies
include(FetchContent)

FetchContent_Declare(
    cmake_openeye
    GIT_REPOSITORY https://github.com/scott-arne/cmake-openeye.git
    GIT_TAG v1.2.0
    GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(cmake_openeye)
list(APPEND CMAKE_MODULE_PATH "${cmake_openeye_SOURCE_DIR}")

if(OEIO_BUILD_PYTHON)
    find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
    find_package(SWIG 4.0)
endif()

if(OPENEYE_USE_SHARED)
    find_package(OpenEyePython REQUIRED)
endif()

find_package(OpenEye REQUIRED)

# Library
set(OEIO_SOURCES
    src/format_registry.cpp
    src/oechem_handler.cpp
    src/mol_range.cpp
    src/read.cpp
    src/write.cpp
    src/pipeline.cpp
)

add_library(oeio ${OEIO_SOURCES})

set_target_properties(oeio PROPERTIES POSITION_INDEPENDENT_CODE ON)

target_include_directories(oeio
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

target_link_libraries(oeio
    PUBLIC
        OpenEye::OEChem
        OpenEye::OEGrid
)

# Tests
if(OEIO_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

# Python bindings
if(OEIO_BUILD_PYTHON AND SWIG_FOUND)
    add_subdirectory(swig)
endif()

# Benchmarks
if(OEIO_BUILD_BENCHMARKS)
    add_subdirectory(benchmarks)
endif()

# Install
install(TARGETS oeio
    EXPORT oeio-targets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)

install(DIRECTORY include/oeio DESTINATION include)

install(EXPORT oeio-targets
    FILE OEIOTargets.cmake
    NAMESPACE OEIO::
    DESTINATION lib/cmake/oeio
)
