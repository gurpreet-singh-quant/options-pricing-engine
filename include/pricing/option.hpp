#ifndef PRICING_OPTION_HPP
#define PRICING_OPTION_HPP

/// \file option.hpp
/// \brief Core option contract description: the Option struct and its enums.

namespace pricing {

/// \brief Whether an option is a call or a put.
enum class OptionType {
    Call,
    Put,
};

/// \brief Exercise style: European (expiry only) or American (any time).
enum class ExerciseStyle {
    European,
    American,
};

/// \brief A single option contract and its current market parameters.
///
/// All rates and the volatility are expressed in annualized, continuous-
/// compounding terms; \p T is measured in years. The struct is a plain
/// aggregate (trivially copyable) so it can live on the stack with no heap
/// involvement in the pricing hot path.
struct Option {
    double S = 0.0;          ///< Underlying spot price.
    double K = 0.0;          ///< Strike price.
    double T = 0.0;          ///< Time to expiry in years.
    double r = 0.0;          ///< Continuously-compounded risk-free rate.
    double sigma = 0.0;      ///< Annualized volatility (e.g. 0.20 == 20%).
    double q = 0.0;          ///< Continuous dividend yield.
    OptionType type = OptionType::Call;          ///< Call or put.
    ExerciseStyle style = ExerciseStyle::European; ///< European or American.
};

/// \brief Price + the five primary Greeks for one option.
///
/// Conventions (all "raw" / per-unit, no 1% or per-day scaling):
///  - delta : dPrice/dS              (per 1.0 of underlying)
///  - gamma : d2Price/dS2
///  - vega  : dPrice/dSigma          (per 1.0 of volatility, i.e. per 100%)
///  - theta : dPrice/dt = -dPrice/dT (per 1.0 year; negative = decay)
///  - rho   : dPrice/dr              (per 1.0 of rate, i.e. per 100%)
struct PriceGreeks {
    double price = 0.0;
    double delta = 0.0;
    double gamma = 0.0;
    double vega = 0.0;
    double theta = 0.0;
    double rho = 0.0;
};

/// \brief Intrinsic (immediate-exercise) value of an option, ignoring time value.
/// \return max(S-K,0) for a call, max(K-S,0) for a put.
[[nodiscard]] constexpr double intrinsic_value(double S, double K, OptionType type) noexcept {
    if (type == OptionType::Call) {
        return S > K ? S - K : 0.0;
    }
    return K > S ? K - S : 0.0;
}

} // namespace pricing

#endif // PRICING_OPTION_HPP
