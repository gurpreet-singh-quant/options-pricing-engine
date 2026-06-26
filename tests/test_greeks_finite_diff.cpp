/// \file test_greeks_finite_diff.cpp
/// \brief Cross-check every analytical Greek against a centered finite-
///        difference of the price function. Analytical and numerical must
///        agree to within 1e-6.
///
/// Bump sizes are chosen per variable to balance truncation error O(h^2)
/// against floating-point roundoff O(eps/h), keeping the central-difference
/// error well below the 1e-6 tolerance for the (moderate-moneyness) scenarios
/// exercised here.

#include "pricing/black_scholes.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

using Catch::Approx;
using namespace pricing;

namespace {

struct Scenario { double S, K, T, r, sigma, q; OptionType type; const char* name; };

double px(const Scenario& s, double S, double K, double T, double r,
          double sigma, double q) {
    return black_scholes_price(S, K, T, r, sigma, q, s.type);
}

} // namespace

TEST_CASE("Analytical Greeks match central finite differences within 1e-6",
          "[greeks][finite-diff]") {
    const Scenario scenarios[] = {
        {100, 100, 1.00, 0.05, 0.20, 0.00, OptionType::Call, "ATM call"},
        {100, 100, 1.00, 0.05, 0.20, 0.00, OptionType::Put,  "ATM put"},
        {100, 90,  0.50, 0.03, 0.35, 0.02, OptionType::Call, "ITM call w/ div"},
        {100, 110, 0.75, 0.04, 0.25, 0.01, OptionType::Put,  "ITM put w/ div"},
        {100, 115, 1.50, 0.06, 0.30, 0.03, OptionType::Call, "OTM call w/ div"},
        {50,  45,  0.25, 0.02, 0.40, 0.00, OptionType::Put,  "OTM put"},
        // Wider moneyness to broaden coverage beyond near-ATM.
        {100, 70,  1.00, 0.05, 0.25, 0.01, OptionType::Call, "deeper ITM call"},
        {100, 140, 1.00, 0.05, 0.25, 0.01, OptionType::Call, "deeper OTM call"},
        {100, 140, 1.00, 0.05, 0.25, 0.01, OptionType::Put,  "deeper ITM put"},
        {100, 70,  1.00, 0.05, 0.25, 0.01, OptionType::Put,  "deeper OTM put"},
    };

    for (const auto& s : scenarios) {
        CAPTURE(s.name);
        const PriceGreeks g = black_scholes(s.S, s.K, s.T, s.r, s.sigma, s.q, s.type);

        const double hS = s.S * 1e-4;
        const double hV = s.sigma * 1e-4;
        const double hT = s.T * 1e-4;
        const double hR = 1e-4;

        // Delta = dP/dS, Gamma = d2P/dS2.
        const double pSup = px(s, s.S + hS, s.K, s.T, s.r, s.sigma, s.q);
        const double pSdn = px(s, s.S - hS, s.K, s.T, s.r, s.sigma, s.q);
        const double p0 = px(s, s.S, s.K, s.T, s.r, s.sigma, s.q);
        const double fd_delta = (pSup - pSdn) / (2.0 * hS);
        const double fd_gamma = (pSup - 2.0 * p0 + pSdn) / (hS * hS);

        // Vega = dP/dsigma.
        const double fd_vega =
            (px(s, s.S, s.K, s.T, s.r, s.sigma + hV, s.q) -
             px(s, s.S, s.K, s.T, s.r, s.sigma - hV, s.q)) / (2.0 * hV);

        // Theta = -dP/dT.
        const double fd_theta =
            -(px(s, s.S, s.K, s.T + hT, s.r, s.sigma, s.q) -
              px(s, s.S, s.K, s.T - hT, s.r, s.sigma, s.q)) / (2.0 * hT);

        // Rho = dP/dr.
        const double fd_rho =
            (px(s, s.S, s.K, s.T, s.r + hR, s.sigma, s.q) -
             px(s, s.S, s.K, s.T, s.r - hR, s.sigma, s.q)) / (2.0 * hR);

        REQUIRE(g.delta == Approx(fd_delta).margin(1e-6));
        REQUIRE(g.gamma == Approx(fd_gamma).margin(1e-6));
        REQUIRE(g.vega == Approx(fd_vega).margin(1e-6));
        REQUIRE(g.theta == Approx(fd_theta).margin(1e-6));
        REQUIRE(g.rho == Approx(fd_rho).margin(1e-6));
    }
}
