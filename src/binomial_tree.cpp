/// \file binomial_tree.cpp
/// \brief CRR binomial tree implementation with reusable buffer and FD Greeks.
///
/// Lattice construction (Cox-Ross-Rubinstein)
/// ------------------------------------------
///     dt   = T / N
///     u    = exp(sigma * sqrt(dt)),  d = 1 / u
///     a    = exp((r - q) * dt)                 (risk-neutral growth per step)
///     p    = (a - d) / (u - d)                 (up probability)
///     disc = exp(-r * dt)
///
/// Terminal spot at node j (0..N) is S * u^(2j - N). We roll back through the
/// lattice; for American options each node takes max(intrinsic, continuation),
/// capturing the early-exercise premium.
///
/// Stability: if the effective total volatility is too small relative to the
/// drift (a leaves the interval (d, u)), the risk-neutral probability p leaves
/// [0, 1] by a wide margin; the lattice can no longer represent a valid measure,
/// so we fall back to the deterministic limit rather than silently clamping to a
/// wrong value. Only last-ulp excursions are absorbed by the clamp. Degenerate
/// inputs (T <= 0, sigma <= 0, S <= 0, N < 1) bypass the lattice entirely and
/// return a finite closed-form fallback.
///
/// Greeks: closed-form Greeks do not exist for American options, so we use
/// CENTERED FINITE DIFFERENCES over re-priced trees -- explicitly an
/// approximation, labeled as such in the header and here.

#include "pricing/binomial_tree.hpp"
#include "pricing/black_scholes.hpp"
#include "pricing/numerics.hpp"

#include <algorithm>
#include <cmath>

namespace pricing {

BinomialTree::BinomialTree(int steps)
    : steps_(steps < 1 ? 1 : steps), buf_(static_cast<std::size_t>(steps_) + 1) {}

double BinomialTree::price_impl(double S, double K, double T, double r,
                                double sigma, double q, OptionType type,
                                ExerciseStyle style) noexcept {
    const int N = steps_;
    const bool american = (style == ExerciseStyle::American);

    // Deterministic-limit fallback shared by every degenerate path. For a
    // European contract the zero-vol value IS the deterministic forward value;
    // an American contract is worth at least that and at least intrinsic-now.
    // (NOTE: for an American option this is only a lower bound on the exact
    // optimal-stopping value in the measure-zero sigma->0 regime, which is an
    // edge case; it is finite, arbitrage-consistent, and matches the BSM
    // pricer for European contracts.)
    auto deterministic = [&]() noexcept {
        const double euro = black_scholes_price(S, K, T, r, sigma, q, type);
        return american ? std::max(intrinsic_value(S, K, type), euro) : euro;
    };

    // --- Degenerate guards: keep everything finite -------------------------
    if (T <= 0.0 || S <= 0.0) {
        return intrinsic_value(S, K, type);
    }
    if (sigma <= 0.0) {
        return deterministic();
    }

    const double dt = T / N;
    const double u = numerics::safe_exp(sigma * std::sqrt(dt));
    const double d = 1.0 / u;
    const double a = numerics::safe_exp((r - q) * dt);
    double p = (a - d) / (u - d);

    // The risk-neutral up-probability must lie in [0, 1]. When the effective
    // total volatility is tiny relative to the drift spacing -- i.e. a falls
    // outside (d, u) -- p escapes [0, 1] by a wide margin and clamping it would
    // silently return a deterministic-up/down lattice value that is materially
    // wrong (and inconsistent with the BSM pricer). Detect that here and fall
    // back to the deterministic limit; only absorb genuine last-ulp excursions
    // with the clamp below.
    constexpr double kPTol = 1e-9;
    if (p < -kPTol || p > 1.0 + kPTol) {
        return deterministic();
    }
    p = std::clamp(p, 0.0, 1.0);

    const double disc = numerics::safe_exp(-r * dt);
    const double pu = disc * p;
    const double pd = disc * (1.0 - p);

    const double logd = std::log(d); // = -sigma*sqrt(dt); hoisted out of loops

    // Terminal payoffs. Spot at node j is S * u^(2j-N) = S * d^N * u^(2j).
    // Build iteratively to avoid pow() per node.
    double spot = S * numerics::safe_exp(static_cast<double>(N) * logd);
    const double u2 = u * u;
    for (int j = 0; j <= N; ++j) {
        buf_[static_cast<std::size_t>(j)] = intrinsic_value(spot, K, type);
        spot *= u2;
    }

    // Backward induction.
    for (int i = N - 1; i >= 0; --i) {
        // Spot at node (i, 0) = S * u^(-i) = S * d^i.
        double node_spot = S * numerics::safe_exp(static_cast<double>(i) * logd);
        for (int j = 0; j <= i; ++j) {
            const std::size_t jj = static_cast<std::size_t>(j);
            const double cont = pu * buf_[jj + 1] + pd * buf_[jj];
            if (american) {
                const double exercise = intrinsic_value(node_spot, K, type);
                buf_[jj] = exercise > cont ? exercise : cont;
            } else {
                buf_[jj] = cont;
            }
            node_spot *= u2;
        }
    }

    return buf_[0];
}

double BinomialTree::price(const Option& opt) noexcept {
    return price_impl(opt.S, opt.K, opt.T, opt.r, opt.sigma, opt.q, opt.type,
                      opt.style);
}

PriceGreeks BinomialTree::price_greeks(const Option& opt) noexcept {
    // Centered finite differences. Bumps are relative where the quantity has a
    // natural scale, and floored so they never collapse to zero.
    const double S = opt.S, K = opt.K, T = opt.T, r = opt.r, sig = opt.sigma,
                 q = opt.q;
    const OptionType ty = opt.type;
    const ExerciseStyle st = opt.style;

    auto P = [&](double s, double t, double sg, double rr) {
        return price_impl(s, K, t, rr, sg, q, ty, st);
    };

    PriceGreeks g{};
    g.price = P(S, T, sig, r);

    // Spot bumps -> delta, gamma.
    const double hS = std::max(1e-4, S * 1e-3);
    const double pUp = P(S + hS, T, sig, r);
    const double pDn = P(S - hS, T, sig, r);
    g.delta = (pUp - pDn) / (2.0 * hS);
    g.gamma = (pUp - 2.0 * g.price + pDn) / (hS * hS);

    // Volatility bump -> vega (per 1.0 of vol).
    const double hV = std::max(1e-5, sig * 1e-3);
    g.vega = (P(S, T, sig + hV, r) - P(S, T, sig - hV, r)) / (2.0 * hV);

    // Time bump -> theta = -dPrice/dT (per year). hT is always strictly
    // positive (floored), so this never divides by zero even at T = 0; near or
    // at expiry we fall back to a one-sided forward difference.
    const double hT = std::max(1e-6, T * 1e-3);
    if (T - hT > 0.0) {
        g.theta = -(P(S, T + hT, sig, r) - P(S, T - hT, sig, r)) / (2.0 * hT);
    } else {
        g.theta = -(P(S, T + hT, sig, r) - g.price) / hT;
    }

    // Rate bump -> rho (per 1.0 of rate).
    const double hR = 1e-5;
    g.rho = (P(S, T, sig, r + hR) - P(S, T, sig, r - hR)) / (2.0 * hR);

    return g;
}

} // namespace pricing
