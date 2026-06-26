#ifndef PRICING_OPTION_CHAIN_HPP
#define PRICING_OPTION_CHAIN_HPP

/// \file option_chain.hpp
/// \brief Structure-of-Arrays option chain and its batch pricing API.

#include "pricing/option.hpp"

#include <cstddef>
#include <span>
#include <vector>

namespace pricing {

/// \brief Batch pricing results in Structure-of-Arrays form.
///
/// One parallel array per output. results.price[i], results.delta[i], ...
/// correspond to the i-th option in the chain.
struct ChainResults {
    std::vector<double> price;
    std::vector<double> delta;
    std::vector<double> gamma;
    std::vector<double> vega;
    std::vector<double> theta;
    std::vector<double> rho;

    /// Resize all arrays to \p n (no-op if already that size).
    void resize(std::size_t n);
    [[nodiscard]] std::size_t size() const noexcept { return price.size(); }
};

/// \brief An option chain stored Structure-of-Arrays (SoA) for cache-friendly,
///        auto-vectorizable batch pricing.
///
/// Inputs are kept as parallel arrays (one per field) rather than an array of
/// Option structs, so a batch price walks each field linearly through memory.
/// Results are written into a pre-sized ChainResults, so a steady-state pricing
/// call performs NO heap allocation -- the hot path the real-time feed hammers.
///
/// Typical real-time use: build the chain once, then on each market tick mutate
/// the spot array in place via spots() and call price_european_batch().
class OptionChain {
public:
    /// Reserve capacity for \p n options across all field arrays.
    void reserve(std::size_t n);

    /// Append one option (any field values; style recorded but the batch
    /// closed-form API prices European -- see price_european_batch).
    void add(const Option& opt);
    void add(double S, double K, double T, double r, double sigma, double q,
             OptionType type, ExerciseStyle style = ExerciseStyle::European);

    /// Number of options in the chain.
    [[nodiscard]] std::size_t size() const noexcept { return s_.size(); }
    [[nodiscard]] bool empty() const noexcept { return s_.empty(); }
    void clear() noexcept;

    /// Mutable view of the spot array for fast in-place tick updates.
    [[nodiscard]] std::span<double> spots() noexcept { return s_; }
    [[nodiscard]] std::span<const double> spots() const noexcept { return s_; }

    /// Set the spot of option \p i (bounds-unchecked; caller guarantees i<size).
    void set_spot(std::size_t i, double S) noexcept { s_[i] = S; }

    /// \brief Price every option in the chain as European (closed-form BSM),
    ///        computing price + all five Greeks in one pass.
    ///
    /// Writes into the internal results buffer (resized once on first call or
    /// after the chain grows) and returns a const reference to it. No heap
    /// allocation occurs once the results buffer matches the chain size.
    const ChainResults& price_european_batch();

    /// Access the most recently computed results.
    [[nodiscard]] const ChainResults& results() const noexcept { return results_; }

private:
    // Structure-of-Arrays input fields.
    std::vector<double> s_, k_, t_, r_, sigma_, q_;
    std::vector<OptionType> type_;
    std::vector<ExerciseStyle> style_;
    ChainResults results_;
};

} // namespace pricing

#endif // PRICING_OPTION_CHAIN_HPP
