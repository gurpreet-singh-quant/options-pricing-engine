#!/usr/bin/env bash
#
# Configure (Release), build, and run the Google Benchmark suite.
# Usage:  scripts/run_benchmarks.sh [build-dir]
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
BUILD_DIR="${1:-${PROJECT_DIR}/build}"

echo ">> Configuring (Release) in ${BUILD_DIR}"
cmake -S "${PROJECT_DIR}" -B "${BUILD_DIR}" -G Ninja -DCMAKE_BUILD_TYPE=Release

echo ">> Building benchmark target"
cmake --build "${BUILD_DIR}" --target pricing_benchmarks

echo ">> Running benchmarks"
"${BUILD_DIR}/pricing_benchmarks" \
    --benchmark_report_aggregates_only=true \
    --benchmark_time_unit=us \
    "$@"
