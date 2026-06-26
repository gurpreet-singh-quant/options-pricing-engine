#ifndef PRICING_NUMERICS_HPP
#define PRICING_NUMERICS_HPP

/// \file numerics.hpp
/// \brief Numerically stable building blocks: standard-normal CDF/PDF, a
///        saturating exp, and small edge-case helpers.
///
/// These primitives are the foundation every pricer depends on, so stability
/// across the *entire* real line (not just the typical |x| < 8 range) is the
/// design goal. See numerics.cpp for the algorithm choices and the rationale
/// behind each one.

namespace pricing::numerics {

/// \brief A small positive threshold below which a quantity (T, sigma, or
///        sigma*sqrt(T)) is treated as effectively zero.
///
/// Chosen well above the smallest representable double so that the reciprocal
/// 1/EPS_TIME stays finite and the d1/d2 denominator never underflows to a
/// denormal. Any total volatility sigma*sqrt(T) below this triggers the
/// closed-form degenerate (intrinsic-value) branch in the pricer.
inline constexpr double EPS_TIME = 1e-12;

/// \brief Standard normal probability density function, phi(x).
///
/// phi(x) = exp(-x^2 / 2) / sqrt(2*pi). Stable on the whole real line: the
/// argument to exp is always <= 0, so the result is in (0, 1/sqrt(2*pi)] and
/// underflows smoothly to +0 for large |x| with no risk of overflow.
[[nodiscard]] double norm_pdf(double x) noexcept;

/// \brief Standard normal cumulative distribution function, Phi(x).
///
/// Computed as Phi(x) = 0.5 * erfc(-x / sqrt(2)) using the C++ standard library
/// erfc. Using erfc rather than (1 + erf)/2 preserves accuracy in the left tail
/// (no catastrophic cancellation), and erfc is engineered to stay accurate as
/// its argument grows, so Phi is monotone and correct across the full real line.
/// Returns values clamped to the closed interval [0, 1].
[[nodiscard]] double norm_cdf(double x) noexcept;

/// \brief exp(x) saturated to avoid overflow to +inf and gradual underflow noise.
///
/// For x above ~709.78 std::exp(x) overflows to +inf; for very negative x it
/// underflows to 0. This helper clamps the exponent into a safe band so callers
/// in deep-ITM / deep-OTM regimes get a large-but-finite or exactly-zero result
/// instead of inf/NaN propagating through later arithmetic.
[[nodiscard]] double safe_exp(double x) noexcept;

/// \brief True if \p x is finite (neither NaN nor +/-inf).
[[nodiscard]] bool is_finite(double x) noexcept;

} // namespace pricing::numerics

#endif // PRICING_NUMERICS_HPP
