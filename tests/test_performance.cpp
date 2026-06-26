/// \file test_performance.cpp
/// \brief Performance guardrail: pricing a 50-option chain (price + all five
///        Greeks) must complete in under 5 ms. This is a Catch2 test, so a
///        regression fails `ctest` (and therefore any CI gate that runs it),
///        not the compile step itself.

#include "pricing/option_chain.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <vector>

using namespace pricing;
using Clock = std::chrono::steady_clock;

namespace {

OptionChain make_chain_50() {
    OptionChain chain;
    chain.reserve(50);
    const double expiries[] = {1.0 / 12, 3.0 / 12, 6.0 / 12, 1.0, 2.0};
    const double moneyness[] = {0.80, 0.90, 0.95, 1.00, 1.05,
                                1.10, 1.20, 1.30, 1.40, 1.50};
    int n = 0;
    for (double T : expiries)
        for (double m : moneyness) {
            chain.add(100.0, 100.0 * m, T, 0.03, 0.20, 0.01,
                      (n % 2 == 0) ? OptionType::Call : OptionType::Put);
            ++n;
        }
    return chain;
}

} // namespace

TEST_CASE("Batch of 50: price + 5 Greeks under 5 ms", "[performance]") {
    OptionChain chain = make_chain_50();
    REQUIRE(chain.size() == 50);

    // Warm up (populate caches, page in).
    volatile double sink = 0.0;
    for (int i = 0; i < 50; ++i) {
        const auto& r = chain.price_european_batch();
        sink += r.price[0];
    }

    // Measure the median over many repetitions.
    constexpr int kReps = 500;
    std::vector<double> times_us;
    times_us.reserve(kReps);
    for (int i = 0; i < kReps; ++i) {
        const auto t0 = Clock::now();
        const auto& r = chain.price_european_batch();
        const auto t1 = Clock::now();
        sink += r.price[0]; // defeat dead-code elimination
        times_us.push_back(std::chrono::duration<double, std::micro>(t1 - t0).count());
    }
    (void)sink;

    std::sort(times_us.begin(), times_us.end());
    const double median_us = times_us[times_us.size() / 2];
    const double p99_us = times_us[static_cast<std::size_t>(0.99 * times_us.size())];

    std::cout << "[performance] batch-of-50 price+Greeks: median "
              << median_us << " us, p99 " << p99_us << " us (limit 5000 us)\n";

    // Hard gate: median must be under 5 ms.
    REQUIRE(median_us < 5000.0);

    // Sanity: every result is populated and finite (value-correctness vs. the
    // scalar pricer is covered by the batch test in test_pricing_accuracy).
    const auto& res = chain.results();
    REQUIRE(res.size() == 50);
    for (std::size_t i = 0; i < res.size(); ++i) {
        REQUIRE(std::isfinite(res.price[i]));
        REQUIRE(std::isfinite(res.delta[i]));
        REQUIRE(std::isfinite(res.gamma[i]));
        REQUIRE(std::isfinite(res.vega[i]));
        REQUIRE(std::isfinite(res.theta[i]));
        REQUIRE(std::isfinite(res.rho[i]));
    }
}
