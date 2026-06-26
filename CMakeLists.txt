cmake_minimum_required(VERSION 3.24)

# ----------------------------------------------------------------------------
# Real-Time Options Pricing Engine
#
# C++20, closed-form Black-Scholes-Merton + analytical Greeks, CRR binomial
# tree for American options, and a cache-friendly Structure-of-Arrays batch
# pricing API. Numerical stability is a primary deliverable.
# ----------------------------------------------------------------------------
project(options_pricing_engine
        VERSION 1.0.0
        DESCRIPTION "Real-time options pricing engine with analytical Greeks"
        LANGUAGES CXX)

# --- Standard --------------------------------------------------------------
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Default to Release so the optimized hot path is what we build and benchmark.
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# --- Warning / optimization flags applied to OUR targets only --------------
# (Fetched/installed dependencies bring their own flags; their headers are
#  consumed as SYSTEM includes via imported targets so they do not pollute
#  our zero-warnings build.)
add_library(pricing_flags INTERFACE)
if(MSVC)
    target_compile_options(pricing_flags INTERFACE
        /W4 /WX /permissive- /Zc:__cplusplus /EHsc
        $<$<CONFIG:Release>:/O2>)
else()
    target_compile_options(pricing_flags INTERFACE
        -Wall -Wextra -Wpedantic -Werror
        $<$<CONFIG:Release>:-O3>
        $<$<CONFIG:Release>:-march=native>)
endif()

# ----------------------------------------------------------------------------
# Dependencies: Catch2 v3 (tests) and Google Benchmark (benchmarks).
#
# Declared via FetchContent per project requirements so a fresh checkout needs
# no manual install. FIND_PACKAGE_ARGS makes FetchContent first try a locally
# installed package (find_package) and only clone from upstream if none is
# present -- this keeps configure fast and works offline when the packages are
# already available on the machine.
# ----------------------------------------------------------------------------
include(FetchContent)

set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "" FORCE)
set(BENCHMARK_ENABLE_INSTALL OFF CACHE BOOL "" FORCE)
set(BENCHMARK_ENABLE_GTEST_TESTS OFF CACHE BOOL "" FORCE)
set(CATCH_INSTALL_DOCS OFF CACHE BOOL "" FORCE)

FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG        v3.5.2
    FIND_PACKAGE_ARGS 3 CONFIG)

FetchContent_Declare(
    benchmark
    GIT_REPOSITORY https://github.com/google/benchmark.git
    GIT_TAG        v1.8.3
    FIND_PACKAGE_ARGS CONFIG)

FetchContent_MakeAvailable(Catch2 benchmark)

# The Catch.cmake helper (catch_discover_tests) ships in the package/clone
# "extras" directory. find_package(Catch2) already puts it on the module path;
# for the fetched case add it explicitly.
if(catch2_SOURCE_DIR)
    list(APPEND CMAKE_MODULE_PATH "${catch2_SOURCE_DIR}/extras")
endif()
include(CTest)
include(Catch)

# ----------------------------------------------------------------------------
# Core library
# ----------------------------------------------------------------------------
add_library(pricing STATIC
    src/numerics.cpp
    src/black_scholes.cpp
    src/binomial_tree.cpp
    src/option_chain.cpp)
target_include_directories(pricing PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(pricing PRIVATE pricing_flags)
set_target_properties(pricing PROPERTIES POSITION_INDEPENDENT_CODE ON)

# ----------------------------------------------------------------------------
# Real-time demo executable
# ----------------------------------------------------------------------------
add_executable(pricing_demo src/main.cpp)
target_link_libraries(pricing_demo PRIVATE pricing pricing_flags)

# ----------------------------------------------------------------------------
# Tests
# ----------------------------------------------------------------------------
enable_testing()
add_executable(unit_tests
    tests/test_numerics.cpp
    tests/test_pricing_accuracy.cpp
    tests/test_greeks_finite_diff.cpp
    tests/test_numerical_stability.cpp
    tests/test_performance.cpp)
target_link_libraries(unit_tests PRIVATE pricing pricing_flags Catch2::Catch2WithMain)
catch_discover_tests(unit_tests)

# ----------------------------------------------------------------------------
# Benchmarks
# ----------------------------------------------------------------------------
add_executable(pricing_benchmarks benchmarks/benchmark_main.cpp)
target_link_libraries(pricing_benchmarks PRIVATE pricing pricing_flags benchmark::benchmark)
