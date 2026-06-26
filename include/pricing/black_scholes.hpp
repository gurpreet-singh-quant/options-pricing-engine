#ifndef PRICING_BLACK_SCHOLES_HPP
#define PRICING_BLACK_SCHOLES_HPP

/// \file black_scholes.hpp
/// \brief Closed-form Black-Scholes-Merton pricer with analytical Greeks.
///
/// All Greeks are derived analytically (not by finite-difference bumping). The
/// implementation guards every edge case (T -> 0, sigma -> 0, S = 0, deep
/// ITM/OTM): for any *finite* input it returns finite results (never NaN or
/// +/-inf). Non-finite inputs (NaN/inf parameters) are not screened in the hot
/// path and propagate as usual -- callers are responsible for passing finite
/// data. See black_scholes.cpp for the limiting-value contract.

#include "pricing/option.hpp"

namespace pricing {

/// \brief Price + all five analytical Greeks for a European option.
///
/// Hot-path friendly: takes scalars by value, allocates nothing, and is
/// `noexcept`. This scalar overload is the kernel the batch API loops over.
///
/// \param S      Underlying spot price.
/// \param K      Strike.
/// \param T      Time to expiry in years.
/// \param r      Continuously-compounded risk-free rate.
/// \param sigma  Annualized volatility.
/// \param q      Continuous dividend yield.
/// \param type   Call or put.
/// \return Price and the five Greeks; finite for any finite input.
[[nodiscard]] PriceGreeks black_scholes(double S, double K, double T,
                                        double r, double sigma, double q,
                                        OptionType type) noexcept;

/// \brief Convenience overload taking an Option. Prices it as European
///        regardless of the contract's exercise style (use the binomial tree
///        for American early-exercise value).
[[nodiscard]] PriceGreeks black_scholes(const Option& opt) noexcept;

/// \brief Just the European price (no Greeks). Handy for finite-difference
///        cross-checks in the test suite.
[[nodiscard]] double black_scholes_price(double S, double K, double T,
                                         double r, double sigma, double q,
                                         OptionType type) noexcept;

} // namespace pricing

#endif // PRICING_BLACK_SCHOLES_HPP
