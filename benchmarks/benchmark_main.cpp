/// \file benchmark_main.cpp
/// \brief Google Benchmark suite: single option, batch of 50, batch of 500.
///
/// Each benchmark is repeated many times; we report the built-in median plus a
/// custom p99 statistic so latency tails are visible. Run with, e.g.:
///   ./pricing_benchmarks --benchmark_report_aggregates_only=true

#include "pricing/black_scholes.hpp"
#include "pricing/binomial_tree.hpp"
#include "pricing/option_chain.hpp"

#include <benchmark/benchmark.h>

#include <algorithm>
#include <cmath>
#include <vector>

using namespace pricing;

namespace {

// p99 over the per-repetition measurements (same time unit Benchmark reports).
double p99(const std::vector<double>& v) {
    if (v.empty()) return 0.0;
    std::vector<double> s = v;
    std::sort(s.begin(), s.end());
    const std::size_t idx =
        static_cast<std::size_t>(0.99 * static_cast<double>(s.size() - 1) + 0.5);
    return s[idx];
}

OptionChain make_chain(std::size_t count) {
    OptionChain chain;
    chain.reserve(count);
    for (std::size_t i = 0; i < count; ++i) {
        const double m = 0.7 + 0.6 * (static_cast<double>(i % 25) / 24.0); // 0.7..1.3
        const double T = 0.05 + 2.0 * (static_cast<double>(i % 7) / 6.0);  // ~0.05..2.05
        chain.add(100.0, 100.0 * m, T, 0.03, 0.20, 0.01,
                  (i % 2 == 0) ? OptionType::Call : OptionType::Put);
    }
    return chain;
}

} // namespace

// --- Single option: price + all 5 Greeks (the <1 microsecond target) --------
static void BM_SingleOption(benchmark::State& state) {
    double S = 100.0;
    for (auto _ : state) {
        S = (S > 130.0) ? 70.0 : S + 0.01; // vary inputs so nothing is hoisted
        PriceGreeks g = black_scholes(S, 100.0, 1.0, 0.03, 0.20, 0.01, OptionType::Call);
        benchmark::DoNotOptimize(g);
    }
    state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_SingleOption)
    ->Repetitions(50)
    ->ComputeStatistics("p99", p99)
    ->ReportAggregatesOnly(true);

// --- Batch of 50 (the explicit <5 ms chain target) --------------------------
static void BM_Batch50(benchmark::State& state) {
    OptionChain chain = make_chain(50);
    double bump = 0.0;
    for (auto _ : state) {
        bump += 0.01;
        auto spots = chain.spots();
        std::fill(spots.begin(), spots.end(), 100.0 + std::sin(bump));
        const auto& r = chain.price_european_batch();
        // Observe actual computed outputs (not just the buffer pointer) via
        // mutable locals -- the const-ref DoNotOptimize overload is deprecated
        // precisely because it can permit the optimizer to weaken the work --
        // and clobber memory so writes to every result array are observable.
        double obs_p = r.price[0];
        double obs_d = r.delta[0];
        benchmark::DoNotOptimize(obs_p);
        benchmark::DoNotOptimize(obs_d);
        benchmark::ClobberMemory();
    }
    state.SetItemsProcessed(state.iterations() * 50);
}
BENCHMARK(BM_Batch50)
    ->Repetitions(50)
    ->ComputeStatistics("p99", p99)
    ->ReportAggregatesOnly(true);

// --- Batch of 500 -----------------------------------------------------------
static void BM_Batch500(benchmark::State& state) {
    OptionChain chain = make_chain(500);
    double bump = 0.0;
    for (auto _ : state) {
        bump += 0.01;
        auto spots = chain.spots();
        std::fill(spots.begin(), spots.end(), 100.0 + std::sin(bump));
        const auto& r = chain.price_european_batch();
        double obs_p = r.price[0];
        double obs_d = r.delta[0];
        benchmark::DoNotOptimize(obs_p);
        benchmark::DoNotOptimize(obs_d);
        benchmark::ClobberMemory();
    }
    state.SetItemsProcessed(state.iterations() * 500);
}
BENCHMARK(BM_Batch500)
    ->Repetitions(50)
    ->ComputeStatistics("p99", p99)
    ->ReportAggregatesOnly(true);

// --- American option via CRR tree (reference cost of the tree path) ---------
static void BM_AmericanTree500Steps(benchmark::State& state) {
    BinomialTree tree(500);
    Option opt;
    opt.S = 100; opt.K = 100; opt.T = 1.0; opt.r = 0.05; opt.sigma = 0.25;
    opt.q = 0.02; opt.type = OptionType::Put; opt.style = ExerciseStyle::American;
    for (auto _ : state) {
        double v = tree.price(opt);
        benchmark::DoNotOptimize(v);
    }
}
BENCHMARK(BM_AmericanTree500Steps)
    ->Repetitions(20)
    ->ComputeStatistics("p99", p99)
    ->ReportAggregatesOnly(true);

BENCHMARK_MAIN();
