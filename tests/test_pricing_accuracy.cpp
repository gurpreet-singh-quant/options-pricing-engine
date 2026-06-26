/// \file test_pricing_accuracy.cpp
/// \brief Black-Scholes price accuracy against known reference values, plus
///        internal consistency checks (put-call parity, tree-vs-closed-form).

#include "pricing/binomial_tree.hpp"
#include "pricing/black_scholes.hpp"
#include "pricing/option_chain.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <iterator>

using Catch::Approx;
using namespace pricing;

TEST_CASE("European call/put match textbook reference", "[accuracy][bsm]") {
    // Hull, "Options, Futures, and Other Derivatives": S=K=100, T=1, r=5%,
    // sigma=20%, q=0  ->  call = 10.450583572, put = 5.573526022.
    const double S = 100, K = 100, T = 1.0, r = 0.05, sigma = 0.20, q = 0.0;

    const double call = black_scholes_price(S, K, T, r, sigma, q, OptionType::Call);
    const double put = black_scholes_price(S, K, T, r, sigma, q, OptionType::Put);

    REQUIRE(call == Approx(10.450583572).margin(1e-6));
    REQUIRE(put == Approx(5.573526022).margin(1e-6));
}

TEST_CASE("Reference Greeks for the ATM example", "[accuracy][greeks]") {
    const double S = 100, K = 100, T = 1.0, r = 0.05, sigma = 0.20, q = 0.0;
    const PriceGreeks c = black_scholes(S, K, T, r, sigma, q, OptionType::Call);

    // d1 = 0.35, d2 = 0.15. (Reference values from textbook formulas; the
    // finite-difference suite is the rigorous per-Greek check.)
    REQUIRE(c.delta == Approx(0.6368306512).margin(1e-7)); // N(d1)
    REQUIRE(c.gamma == Approx(0.0187620173).margin(1e-7));
    REQUIRE(c.vega == Approx(37.52403).margin(1e-3));   // per 1.0 vol
    REQUIRE(c.theta == Approx(-6.414028).margin(1e-3)); // per year
    REQUIRE(c.rho == Approx(53.232482).margin(1e-3));   // per 1.0 rate
}

TEST_CASE("Put-call parity holds to machine precision", "[accuracy][parity]") {
    struct Case { double S, K, T, r, sigma, q; };
    const Case cases[] = {
        {100, 100, 1.0, 0.05, 0.20, 0.00},
        {100, 90, 0.5, 0.03, 0.35, 0.02},
        {42, 50, 2.0, 0.01, 0.15, 0.00},
        {7.5, 5.0, 0.25, 0.08, 0.50, 0.04},
    };
    for (const auto& c : cases) {
        const double call = black_scholes_price(c.S, c.K, c.T, c.r, c.sigma, c.q, OptionType::Call);
        const double put = black_scholes_price(c.S, c.K, c.T, c.r, c.sigma, c.q, OptionType::Put);
        // C - P = S e^{-qT} - K e^{-rT}
        const double lhs = call - put;
        const double rhs = c.S * std::exp(-c.q * c.T) - c.K * std::exp(-c.r * c.T);
        REQUIRE(lhs == Approx(rhs).margin(1e-10));
    }
}

TEST_CASE("CRR European price converges to closed-form BSM", "[accuracy][tree]") {
    // A European-style option priced on the CRR tree must approach the
    // closed-form value as the step count grows.
    Option opt;
    opt.S = 100; opt.K = 105; opt.T = 1.0; opt.r = 0.05; opt.sigma = 0.25;
    opt.q = 0.02; opt.type = OptionType::Call; opt.style = ExerciseStyle::European;

    const double bsm = black_scholes(opt).price;
    BinomialTree tree(2000);
    const double tree_px = tree.price(opt);
    REQUIRE(tree_px == Approx(bsm).margin(2e-3));
}

TEST_CASE("American >= European; American call w/o dividends == European",
          "[accuracy][tree][american]") {
    BinomialTree tree(1000);

    // American put is worth strictly more than European (early-exercise premium).
    Option am_put;
    am_put.S = 100; am_put.K = 110; am_put.T = 1.0; am_put.r = 0.06;
    am_put.sigma = 0.30; am_put.q = 0.0; am_put.type = OptionType::Put;
    am_put.style = ExerciseStyle::American;
    Option eu_put = am_put; eu_put.style = ExerciseStyle::European;

    const double am = tree.price(am_put);
    const double eu = tree.price(eu_put);
    REQUIRE(am >= eu - 1e-9);
    REQUIRE(am > eu); // genuine early-exercise premium for an ITM put

    // With no dividends, an American call equals the European call.
    Option am_call;
    am_call.S = 100; am_call.K = 95; am_call.T = 1.0; am_call.r = 0.05;
    am_call.sigma = 0.25; am_call.q = 0.0; am_call.type = OptionType::Call;
    am_call.style = ExerciseStyle::American;
    const double am_c = tree.price(am_call);
    const double eu_c = black_scholes(am_call).price;
    REQUIRE(am_c == Approx(eu_c).margin(5e-3));
}

TEST_CASE("Batch SoA API reproduces the scalar pricer exactly",
          "[accuracy][batch]") {
    // Every batch result must equal the single-option closed form bit-for-bit
    // (same kernel), for price AND all five Greeks -- not merely be finite.
    struct Row { double S, K, T, r, sigma, q; OptionType type; };
    const Row rows[] = {
        {100, 100, 1.00, 0.05, 0.20, 0.00, OptionType::Call},
        {100, 120, 0.50, 0.03, 0.35, 0.02, OptionType::Put},
        {42,  40,  2.00, 0.01, 0.15, 0.01, OptionType::Call},
        {100, 100, 0.00, 0.05, 0.20, 0.00, OptionType::Put},   // expiry edge
        {100, 90,  1.00, 0.05, 0.00, 0.02, OptionType::Call},  // zero-vol edge
        {0.0, 100, 1.00, 0.05, 0.20, 0.01, OptionType::Put},   // S = 0 edge
        {1e6, 1.0, 1.00, 0.03, 0.20, 0.01, OptionType::Call},  // deep ITM
    };

    OptionChain chain;
    for (const auto& r : rows)
        chain.add(r.S, r.K, r.T, r.r, r.sigma, r.q, r.type);

    const ChainResults& res = chain.price_european_batch();
    REQUIRE(res.size() == std::size(rows));

    for (std::size_t i = 0; i < std::size(rows); ++i) {
        const auto& r = rows[i];
        const PriceGreeks g =
            black_scholes(r.S, r.K, r.T, r.r, r.sigma, r.q, r.type);
        CAPTURE(i);
        REQUIRE(res.price[i] == g.price);
        REQUIRE(res.delta[i] == g.delta);
        REQUIRE(res.gamma[i] == g.gamma);
        REQUIRE(res.vega[i] == g.vega);
        REQUIRE(res.theta[i] == g.theta);
        REQUIRE(res.rho[i] == g.rho);
    }
}

TEST_CASE("Tree does not silently misprice low effective volatility",
          "[accuracy][tree][stability]") {
    // Regression: for small sigma the risk-neutral up-probability escapes
    // [0,1]; the tree must fall back to the deterministic limit (matching BSM)
    // instead of returning a clamped, materially wrong value.
    BinomialTree tree(500);
    Option opt;
    opt.S = 100; opt.K = 100; opt.r = 0.10; opt.q = 0.0; opt.T = 1.0;
    opt.type = OptionType::Call; opt.style = ExerciseStyle::European;

    for (double sig : {1e-3, 1e-4, 1e-6}) {
        opt.sigma = sig;
        const double tree_px = tree.price(opt);
        const double bsm_px = black_scholes(opt).price;
        CAPTURE(sig);
        REQUIRE(std::isfinite(tree_px));
        REQUIRE(tree_px == Approx(bsm_px).margin(1e-6)); // ~9.52, not clamped ~1.3
    }

    // American put with tiny vol stays finite and respects the intrinsic floor.
    Option am;
    am.S = 90; am.K = 100; am.r = 0.05; am.q = 0.0; am.T = 1.0; am.sigma = 1e-5;
    am.type = OptionType::Put; am.style = ExerciseStyle::American;
    const double am_px = tree.price(am);
    REQUIRE(std::isfinite(am_px));
    REQUIRE(am_px >= intrinsic_value(am.S, am.K, OptionType::Put) - 1e-9);
}
