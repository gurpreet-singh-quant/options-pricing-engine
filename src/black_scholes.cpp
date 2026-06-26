/// \file black_scholes.cpp
/// \brief Black-Scholes-Merton price + analytical Greeks, with explicit,
///        documented limiting behavior at every degenerate edge.
///
/// Numerical stability contract
/// ============================
///
/// d1/d2 without catastrophic cancellation
/// ---------------------------------------
/// The total volatility w = sigma * sqrt(T) is computed once. Then
///
///     d1 = (ln(S/K) + (r - q) * T) / w + 0.5 * w
///     d2 = d1 - w
///
/// d2 is formed by subtracting w directly rather than recomputing
/// (ln(S/K) + (r - q) * T)/w - 0.5*w; this guarantees d1 - d2 == w exactly and
/// avoids subtracting two large, nearly-equal numbers (catastrophic
/// cancellation) when |ln(S/K)| dominates.
///
/// Degenerate / limiting branch
/// ----------------------------
/// When the model degenerates -- expiry (T <= 0), zero volatility (sigma <= 0),
/// vanishing total volatility (w < EPS_TIME), or a non-positive spot/strike --
/// the diffusion term disappears and the underlying is deterministic. We then
/// return the *deterministic forward* value
///
///     F    = S * exp((r - q) * T)         (the forward price; F = S at T = 0)
///     disc = exp(-r * T)                  (disc = 1 at T = 0)
///     price(call) = disc * max(F - K, 0)
///     price(put)  = disc * max(K - F, 0)
///
/// This single expression is correct at BOTH limits and is continuous with the
/// full BSM formula:
///   * As T -> 0 it reduces to the undiscounted intrinsic value max(S-K,0) /
///     max(K-S,0) -- exactly the required limit.
///   * As sigma -> 0 with T > 0 it gives the exact zero-vol (deterministic
///     payoff) value, with no spurious discontinuity at the threshold.
/// Greeks in this branch are the analytic one-sided limits (gamma = vega = 0;
/// delta/theta/rho the deterministic carry terms), all finite -- never NaN/inf.
///
/// Overflow / underflow
/// --------------------
/// All exp() calls go through numerics::safe_exp, which saturates the exponent
/// into the representable band. Deep-ITM/OTM options drive d1/d2 to +/-inf; the
/// CDF (via erfc) saturates cleanly to 1 or 0 and the PDF underflows to +0, so
/// every Greek stays finite.

#include "pricing/black_scholes.hpp"
#include "pricing/numerics.hpp"

#include <cmath>

namespace pricing {

namespace {

using numerics::EPS_TIME;
using numerics::norm_cdf;
using numerics::norm_pdf;
using numerics::safe_exp;

/// Deterministic (zero-diffusion) limit shared by every degenerate case.
/// Teff is assumed already clamped to be >= 0.
PriceGreeks deterministic_limit(double S, double K, double r, double q,
                                double Teff, OptionType type) noexcept {
    // A negative spot is non-physical; floor it at 0 so the forward stays
    // non-negative and the result remains sensible (a negative S would
    // otherwise drive an unbounded, non-physical forward intrinsic).
    if (S < 0.0) S = 0.0;
    const double df_r = safe_exp(-r * Teff);  // 1 at Teff == 0
    const double df_q = safe_exp(-q * Teff);  // 1 at Teff == 0
    const double F = S * safe_exp((r - q) * Teff); // forward; == S at Teff == 0

    PriceGreeks g{}; // zero-initialized: gamma, vega default to 0

    if (type == OptionType::Call) {
        if (F > K) { // in-the-money forward
            g.price = df_r * (F - K);
            g.delta = df_q;                       // d/dS [df_q*S - df_r*K]
            g.theta = q * S * df_q - r * K * df_r; // = -d(price)/dT
            g.rho = Teff * K * df_r;               // d(price)/dr
        } else if (F == K) {
            g.delta = 0.5 * df_q; // kink: average of the one-sided limits
        }
        // OTM: price and all Greeks remain 0.
    } else { // Put
        if (F < K) {
            g.price = df_r * (K - F);
            g.delta = -df_q;
            g.theta = r * K * df_r - q * S * df_q;
            g.rho = -Teff * K * df_r;
        } else if (F == K) {
            g.delta = -0.5 * df_q;
        }
    }
    return g;
}

} // namespace

PriceGreeks black_scholes(double S, double K, double T, double r, double sigma,
                          double q, OptionType type) noexcept {
    const double Teff = T > 0.0 ? T : 0.0;

    // Route every degenerate configuration to the deterministic limit.
    if (T <= 0.0 || sigma <= 0.0 || S <= 0.0 || K <= 0.0 ||
        sigma * std::sqrt(Teff) < EPS_TIME) {
        return deterministic_limit(S, K, r, q, Teff, type);
    }

    // --- Regular BSM regime: S>0, K>0, T>0, sigma>0, w >= EPS_TIME ----------
    const double sqrtT = std::sqrt(T);
    const double w = sigma * sqrtT; // total volatility == d1 - d2 (exact)
    const double df_r = safe_exp(-r * T);
    const double df_q = safe_exp(-q * T);

    const double d1 = (std::log(S / K) + (r - q) * T) / w + 0.5 * w;
    const double d2 = d1 - w; // direct subtraction: no cancellation

    const double Nd1 = norm_cdf(d1);
    const double Nd2 = norm_cdf(d2);
    const double nd1 = norm_pdf(d1);
    // Tail-accurate complements (avoid 1 - Nd which loses digits in the tail).
    const double Nmd1 = norm_cdf(-d1);
    const double Nmd2 = norm_cdf(-d2);

    PriceGreeks g{};

    // Greeks common to calls and puts.
    g.gamma = (df_q * nd1) / (S * w);     // d2 Price / dS2
    g.vega = S * df_q * nd1 * sqrtT;      // per 1.0 of volatility

    const double decay = -(S * df_q * nd1 * sigma) / (2.0 * sqrtT); // shared theta term

    if (type == OptionType::Call) {
        g.price = S * df_q * Nd1 - K * df_r * Nd2;
        g.delta = df_q * Nd1;
        g.theta = decay - r * K * df_r * Nd2 + q * S * df_q * Nd1;
        g.rho = K * T * df_r * Nd2;
    } else {
        g.price = K * df_r * Nmd2 - S * df_q * Nmd1;
        g.delta = -df_q * Nmd1; // == df_q * (Nd1 - 1)
        g.theta = decay + r * K * df_r * Nmd2 - q * S * df_q * Nmd1;
        g.rho = -K * T * df_r * Nmd2;
    }

    return g;
}

PriceGreeks black_scholes(const Option& opt) noexcept {
    return black_scholes(opt.S, opt.K, opt.T, opt.r, opt.sigma, opt.q, opt.type);
}

double black_scholes_price(double S, double K, double T, double r, double sigma,
                           double q, OptionType type) noexcept {
    return black_scholes(S, K, T, r, sigma, q, type).price;
}

} // namespace pricing
