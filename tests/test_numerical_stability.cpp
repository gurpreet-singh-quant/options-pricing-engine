/// \file test_numerical_stability.cpp
/// \brief Edge-case stability: T->0, sigma->0, S=0, deep ITM/OTM. Asserts no
///        NaN/Inf is ever produced and the documented limiting values hold.

#include "pricing/binomial_tree.hpp"
#include "pricing/black_scholes.hpp"
#include "pricing/numerics.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>

using Catch::Approx;
using namespace pricing;

namespace {

void require_all_finite(const PriceGreeks& g) {
    REQUIRE(std::isfinite(g.price));
    REQUIRE(std::isfinite(g.delta));
    REQUIRE(std::isfinite(g.gamma));
    REQUIRE(std::isfinite(g.vega));
    REQUIRE(std::isfinite(g.theta));
    REQUIRE(std::isfinite(g.rho));
}

} // namespace

TEST_CASE("T = 0 returns intrinsic value, finite Greeks", "[stability][expiry]") {
    // In-the-money call: price == S - K.
    {
        const PriceGreeks g = black_scholes(120, 100, 0.0, 0.05, 0.2, 0.0, OptionType::Call);
        require_all_finite(g);
        REQUIRE(g.price == Approx(20.0).margin(1e-12));
        REQUIRE(g.delta == Approx(1.0).margin(1e-12)); // ITM call delta -> 1
    }
    // Out-of-the-money call: price == 0.
    {
        const PriceGreeks g = black_scholes(80, 100, 0.0, 0.05, 0.2, 0.0, OptionType::Call);
        require_all_finite(g);
        REQUIRE(g.price == 0.0);
        REQUIRE(g.delta == 0.0);
    }
    // In-the-money put: price == K - S.
    {
        const PriceGreeks g = black_scholes(80, 100, 0.0, 0.05, 0.2, 0.0, OptionType::Put);
        require_all_finite(g);
        REQUIRE(g.price == Approx(20.0).margin(1e-12));
        REQUIRE(g.delta == Approx(-1.0).margin(1e-12));
    }
    // Exactly at the money: intrinsic 0, half-delta convention.
    {
        const PriceGreeks g = black_scholes(100, 100, 0.0, 0.05, 0.2, 0.0, OptionType::Call);
        require_all_finite(g);
        REQUIRE(g.price == 0.0);
        REQUIRE(g.delta == Approx(0.5).margin(1e-12));
    }
}

TEST_CASE("sigma = 0 gives the deterministic forward value", "[stability][zero-vol]") {
    const double S = 100, K = 95, T = 1.0, r = 0.05, q = 0.02;
    const double F = S * std::exp((r - q) * T);
    const double disc = std::exp(-r * T);

    const PriceGreeks call = black_scholes(S, K, T, r, 0.0, q, OptionType::Call);
    require_all_finite(call);
    REQUIRE(call.price == Approx(disc * std::max(F - K, 0.0)).margin(1e-12));
    REQUIRE(call.gamma == 0.0);
    REQUIRE(call.vega == 0.0);

    const PriceGreeks put = black_scholes(S, 120, T, r, 0.0, q, OptionType::Put);
    require_all_finite(put);
    const double F2 = S * std::exp((r - q) * T);
    REQUIRE(put.price == Approx(disc * std::max(120.0 - F2, 0.0)).margin(1e-12));
}

TEST_CASE("S = 0 is handled with no NaN/Inf", "[stability][zero-spot]") {
    const double K = 100, T = 1.0, r = 0.05, sigma = 0.2, q = 0.01;
    const PriceGreeks call = black_scholes(0.0, K, T, r, sigma, q, OptionType::Call);
    require_all_finite(call);
    REQUIRE(call.price == 0.0); // worthless underlying -> worthless call

    const PriceGreeks put = black_scholes(0.0, K, T, r, sigma, q, OptionType::Put);
    require_all_finite(put);
    // Put receives K at expiry for certain: value == K * e^{-rT}.
    REQUIRE(put.price == Approx(K * std::exp(-r * T)).margin(1e-12));
}

TEST_CASE("Deep ITM / OTM produce no overflow, sane limits", "[stability][deep]") {
    // Deep ITM call: ~ forward intrinsic, delta ~ e^{-qT}, no inf.
    {
        const double S = 1.0e6, K = 1.0, T = 1.0, r = 0.03, sigma = 0.2, q = 0.01;
        const PriceGreeks g = black_scholes(S, K, T, r, sigma, q, OptionType::Call);
        require_all_finite(g);
        REQUIRE(g.delta == Approx(std::exp(-q * T)).margin(1e-9));
        REQUIRE(g.price > 0.0);
    }
    // Deep OTM call: ~ 0, all Greeks ~ 0, finite.
    {
        const double S = 1.0, K = 1.0e6, T = 1.0, r = 0.03, sigma = 0.2, q = 0.01;
        const PriceGreeks g = black_scholes(S, K, T, r, sigma, q, OptionType::Call);
        require_all_finite(g);
        REQUIRE(g.price == Approx(0.0).margin(1e-12));
        REQUIRE(g.delta == Approx(0.0).margin(1e-12));
    }
    // Extremely large volatility: still finite; call -> S e^{-qT} bound.
    {
        const double S = 100, K = 100, T = 1.0, r = 0.03, sigma = 50.0, q = 0.0;
        const PriceGreeks g = black_scholes(S, K, T, r, sigma, q, OptionType::Call);
        require_all_finite(g);
        REQUIRE(g.price <= S + 1e-9);
        REQUIRE(g.price >= 0.0);
    }
}

TEST_CASE("Stress grid: never NaN/Inf across extreme parameter combinations",
          "[stability][grid]") {
    const double Ss[] = {0.0, 1e-6, 1.0, 100.0, 1e6};
    const double Ks[] = {1e-6, 1.0, 100.0, 1e6};
    const double Ts[] = {0.0, 1e-10, 1e-6, 1.0, 30.0};
    const double sigmas[] = {0.0, 1e-8, 0.2, 5.0, 100.0};
    const double rs[] = {-0.05, 0.0, 0.10};
    const double qs[] = {0.0, 0.03};

    for (double S : Ss)
        for (double K : Ks)
            for (double T : Ts)
                for (double sig : sigmas)
                    for (double r : rs)
                        for (double q : qs)
                            for (OptionType ty : {OptionType::Call, OptionType::Put}) {
                                const PriceGreeks g = black_scholes(S, K, T, r, sig, q, ty);
                                require_all_finite(g);
                                REQUIRE(g.price >= -1e-9); // never materially negative
                            }
}

TEST_CASE("Binomial tree stays finite on degenerate inputs", "[stability][tree]") {
    BinomialTree tree(200);
    Option o;
    o.type = OptionType::Put;
    o.style = ExerciseStyle::American;

    for (double S : {0.0, 1.0, 100.0})
        for (double T : {0.0, 1e-8, 1.0})
            for (double sig : {0.0, 0.2, 5.0}) {
                o.S = S; o.K = 100.0; o.T = T; o.r = 0.05; o.sigma = sig; o.q = 0.0;
                const double p = tree.price(o);
                REQUIRE(std::isfinite(p));
                REQUIRE(p >= -1e-9);
                const PriceGreeks g = tree.price_greeks(o);
                require_all_finite(g);
            }
}
