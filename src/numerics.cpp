/// \file numerics.cpp
/// \brief Implementation of the stable normal CDF/PDF and saturating exp.
///
/// Algorithm choices and stability rationale
/// ==========================================
///
/// Standard normal CDF, Phi(x)
/// ---------------------------
/// We use the identity
///
///     Phi(x) = 0.5 * erfc(-x / sqrt(2))
///
/// where erfc is std::erfc from <cmath>. This is a deliberate choice over the
/// once-popular Abramowitz & Stegun / Hastings polynomial approximations:
///
///   * Accuracy on the whole real line. std::erfc is implemented (in glibc,
///     UCRT, libc++ math) with near machine precision and is specifically
///     designed to remain accurate as its argument grows large, which is
///     exactly the tail region (deep ITM/OTM) the engine must handle. The
///     A&S 7.1.26 polynomial, by contrast, carries ~1e-7 absolute error and
///     degrades further in the tails.
///
///   * No catastrophic cancellation in the left tail. Computing Phi as
///     (1 + erf(x/sqrt(2)))/2 loses all significant digits when x is very
///     negative, because erf -> -1 and we subtract two nearly-equal numbers.
///     erfc(-x/sqrt(2)) computes the small tail probability *directly*, so
///     the relative error stays bounded.
///
///   * Monotonicity. Because the single erfc call is monotone, Phi is monotone
///     by construction -- no spurious non-monotone wiggles that polynomial fits
///     can introduce near 0.
///
/// std::erfc is part of the C++ standard library, not a third-party math or
/// finance package, so this stays within the "no external libraries" rule
/// while delegating the hard tail accuracy to a well-tested primitive.
///
/// Standard normal PDF, phi(x)
/// ---------------------------
/// phi(x) = exp(-x^2/2) / sqrt(2*pi). The exponent is always <= 0, so this can
/// only underflow toward +0 (which is the correct limit), never overflow.
///
/// safe_exp
/// --------
/// std::exp overflows to +inf for arguments above ~709.78 (DBL_MAX_EXP). In the
/// pricer the only place a positive, potentially large exponent appears is the
/// up-factor of the binomial tree and discounting with negative rates; we clamp
/// the exponent to a band that maps to [smallest normal, ~DBL_MAX] so the result
/// is always finite, letting callers continue arithmetic without special-casing
/// inf. The clamp bounds are intentionally just inside the representable range.

#include "pricing/numerics.hpp"

#include <cmath>

namespace pricing::numerics {

namespace {
// 1/sqrt(2): scales x for the erfc-based CDF.
constexpr double kInvSqrt2 = 0.7071067811865475244008443621048490392848359376884740;
// 1/sqrt(2*pi): normalizing constant for the PDF.
constexpr double kInvSqrt2Pi = 0.3989422804014326779399460599343818684758586311649347;

// Largest exponent for which std::exp stays finite (ln(DBL_MAX) ~ 709.7827).
// We back off slightly so the result is comfortably below +inf.
constexpr double kMaxExpArg = 708.0;
// Below this, exp(x) is +0 to double precision (ln(DBL_MIN_DENORM) ~ -745).
constexpr double kMinExpArg = -708.0;
} // namespace

double norm_pdf(double x) noexcept {
    return kInvSqrt2Pi * std::exp(-0.5 * x * x);
}

double norm_cdf(double x) noexcept {
    const double v = 0.5 * std::erfc(-x * kInvSqrt2);
    // erfc is mathematically in [0, 2], so 0.5*erfc is in [0, 1]; clamp guards
    // against a last-ulp excursion so downstream callers never see e.g. 1+1e-17.
    if (v <= 0.0) return 0.0;
    if (v >= 1.0) return 1.0;
    return v;
}

double safe_exp(double x) noexcept {
    if (x > kMaxExpArg) x = kMaxExpArg;
    else if (x < kMinExpArg) return 0.0;
    return std::exp(x);
}

bool is_finite(double x) noexcept {
    return std::isfinite(x);
}

} // namespace pricing::numerics
