/// \file main.cpp
/// \brief Real-time market-data feed simulation.
///
/// Simulates a live underlying price feed: every tick (default ~1 ms apart) the
/// underlying takes a geometric-Brownian-motion step, all options in a ~50-name
/// chain are re-priced (price + all five Greeks) via the SoA batch API, and the
/// wall-clock latency of that recomputation is logged to stdout and to a CSV
/// file. The program runs a fixed number of ticks and then exits, printing
/// summary latency statistics.
///
/// Note on tick spacing: std::this_thread::sleep_for(1ms) is subject to OS
/// scheduler granularity (notably coarse on Windows), so real spacing may
/// exceed 1 ms. This only affects the *cadence* of the feed; the latency we
/// report is the isolated recomputation time, measured with steady_clock.

#include "pricing/option_chain.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <thread>
#include <vector>

namespace {

using pricing::ExerciseStyle;
using pricing::OptionChain;
using pricing::OptionType;
using Clock = std::chrono::steady_clock;

/// Build a ~50-option chain centred on \p spot0: a ladder of strikes across
/// several expiries, alternating calls and puts.
OptionChain build_chain(double spot0) {
    OptionChain chain;
    const double r = 0.03;     // 3% risk-free
    const double q = 0.01;     // 1% dividend yield
    const double sigma = 0.20; // 20% vol
    const double expiries[] = {1.0 / 12, 3.0 / 12, 6.0 / 12, 1.0, 2.0};
    const double moneyness[] = {0.80, 0.90, 0.95, 1.00, 1.05, 1.10, 1.20, 1.30,
                                1.40, 1.50};
    chain.reserve(50);
    int n = 0;
    for (double T : expiries) {
        for (double m : moneyness) {
            const double K = spot0 * m;
            const OptionType ty = (n % 2 == 0) ? OptionType::Call : OptionType::Put;
            chain.add(spot0, K, T, r, sigma, q, ty, ExerciseStyle::European);
            ++n;
        }
    }
    return chain; // exactly 5 * 10 = 50 options
}

double percentile(std::vector<double>& v, double pct) {
    if (v.empty()) return 0.0;
    std::sort(v.begin(), v.end());
    const double idx = (pct / 100.0) * static_cast<double>(v.size() - 1);
    const std::size_t lo = static_cast<std::size_t>(std::floor(idx));
    const std::size_t hi = static_cast<std::size_t>(std::ceil(idx));
    if (lo == hi) return v[lo];
    const double frac = idx - static_cast<double>(lo);
    return v[lo] * (1.0 - frac) + v[hi] * frac;
}

} // namespace

int main(int argc, char** argv) {
    const long num_ticks = (argc > 1) ? std::strtol(argv[1], nullptr, 10) : 1000;
    const long interval_ms = (argc > 2) ? std::strtol(argv[2], nullptr, 10) : 1;

    if (num_ticks <= 0) {
        std::cerr << "usage: " << argv[0]
                  << " [num_ticks>0] [interval_ms>=0]\n";
        return 1;
    }

    double spot = 100.0;
    OptionChain chain = build_chain(spot);
    const std::size_t chain_size = chain.size();

    // GBM tick parameters (per-tick, not annualized; small for a smooth feed).
    std::mt19937_64 rng(12345ULL); // fixed seed -> reproducible run
    std::normal_distribution<double> z(0.0, 1.0);
    const double mu_step = 0.0;
    const double sig_step = 0.0008; // ~0.08% std move per tick

    std::ofstream csv("latency_log.csv");
    if (!csv) {
        std::cerr << "error: could not open latency_log.csv for writing\n";
        return 1;
    }
    csv << "tick,underlying,latency_us\n";
    csv << std::fixed << std::setprecision(4);

    std::vector<double> latencies_us;
    latencies_us.reserve(static_cast<std::size_t>(num_ticks));

    std::cout << "Real-time options pricing feed\n"
              << "  chain size : " << chain_size << " options\n"
              << "  ticks      : " << num_ticks << "\n"
              << "  interval   : " << interval_ms << " ms (subject to OS granularity)\n"
              << "  output CSV : latency_log.csv\n\n";

    for (long tick = 0; tick < num_ticks; ++tick) {
        // GBM step for the shared underlying.
        spot *= std::exp((mu_step - 0.5 * sig_step * sig_step) + sig_step * z(rng));

        // Push the new spot into every option (in-place, no allocation).
        auto spots = chain.spots();
        std::fill(spots.begin(), spots.end(), spot);

        // --- timed recomputation: price + all 5 Greeks for the whole chain ---
        const auto t0 = Clock::now();
        chain.price_european_batch();
        const auto t1 = Clock::now();

        const double us =
            std::chrono::duration<double, std::micro>(t1 - t0).count();
        latencies_us.push_back(us);
        csv << tick << ',' << spot << ',' << us << '\n';

        if (tick % 100 == 0) {
            const auto& res = chain.results();
            std::cout << "tick " << std::setw(5) << tick << " | S=" << std::fixed
                      << std::setprecision(2) << std::setw(8) << spot
                      << " | recompute " << std::setprecision(2) << std::setw(7)
                      << us << " us"
                      << " | px[0]=" << std::setprecision(4) << res.price[0]
                      << " delta[0]=" << res.delta[0] << '\n';
        }

        if (interval_ms > 0) {
            std::this_thread::sleep_for(std::chrono::milliseconds(interval_ms));
        }
    }

    csv.close();

    // Summary statistics.
    std::vector<double> sorted = latencies_us;
    const double mn = *std::min_element(sorted.begin(), sorted.end());
    const double mx = *std::max_element(sorted.begin(), sorted.end());
    const double med = percentile(sorted, 50.0);
    const double p99 = percentile(sorted, 99.0);
    double sum = 0.0;
    for (double x : latencies_us) sum += x;
    const double mean = sum / static_cast<double>(latencies_us.size());

    std::cout << "\nLatency over " << num_ticks << " ticks (recompute of "
              << chain_size << "-option chain, price + 5 Greeks):\n"
              << std::fixed << std::setprecision(3) << "  min    : " << mn << " us\n"
              << "  mean   : " << mean << " us\n"
              << "  median : " << med << " us\n"
              << "  p99    : " << p99 << " us\n"
              << "  max    : " << mx << " us\n"
              << "\nFull per-tick latencies written to latency_log.csv\n";

    return 0;
}
