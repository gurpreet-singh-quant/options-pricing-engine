#ifndef PRICING_BINOMIAL_TREE_HPP
#define PRICING_BINOMIAL_TREE_HPP

/// \file binomial_tree.hpp
/// \brief Cox-Ross-Rubinstein binomial tree pricer for American options.

#include "pricing/option.hpp"

#include <vector>

namespace pricing {

/// \brief Cox-Ross-Rubinstein (CRR) binomial tree pricer.
///
/// Handles American early exercise (and European, by skipping the exercise
/// test). The per-level value buffer is allocated once at construction and
/// reused for every pricing call, so repeated pricing performs no heap
/// allocation. A single instance is therefore *not* thread-safe; give each
/// thread its own tree.
///
/// Greeks are computed by centered finite differences over re-priced trees
/// (see price_greeks). This is the standard approach for American options,
/// whose Greeks have no simple closed form, and is labeled as such.
class BinomialTree {
public:
    /// Default step count, balancing accuracy against cost.
    static constexpr int kDefaultSteps = 500;

    /// \param steps Number of time steps (>= 1); defaults to kDefaultSteps.
    explicit BinomialTree(int steps = kDefaultSteps);

    /// \brief Price a single option via the CRR tree.
    /// \return The option value; intrinsic value for degenerate inputs.
    [[nodiscard]] double price(const Option& opt) noexcept;

    /// \brief Price + all five Greeks via centered finite differences.
    ///
    /// Greeks are FINITE-DIFFERENCE estimates (not closed form): each is a
    /// centered difference of re-priced trees with a small bump in the relevant
    /// parameter. Bumps are chosen to dominate the tree's O(1/steps)
    /// discretization noise while keeping truncation error small.
    [[nodiscard]] PriceGreeks price_greeks(const Option& opt) noexcept;

    /// \brief Number of time steps this tree uses.
    [[nodiscard]] int steps() const noexcept { return steps_; }

private:
    // Core CRR valuation on raw scalars; uses (and resizes once) buf_.
    double price_impl(double S, double K, double T, double r, double sigma,
                      double q, OptionType type, ExerciseStyle style) noexcept;

    int steps_;
    std::vector<double> buf_; ///< Reused per-level value array, size steps_+1.
};

} // namespace pricing

#endif // PRICING_BINOMIAL_TREE_HPP
