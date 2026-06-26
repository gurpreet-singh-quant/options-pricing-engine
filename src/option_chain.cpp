/// \file option_chain.cpp
/// \brief SoA option chain batch pricing.
///
/// The batch loop walks the parallel input arrays linearly and writes the
/// parallel output arrays linearly -- a cache-friendly access pattern with no
/// pointer chasing. The per-element kernel is the closed-form BSM routine; its
/// transcendental calls (erfc/exp/log) dominate cost and are not auto-
/// vectorized by the compiler, but the surrounding arithmetic and memory
/// traffic are linear and branch-light, and no allocation occurs in the loop.

#include "pricing/option_chain.hpp"
#include "pricing/black_scholes.hpp"

namespace pricing {

void ChainResults::resize(std::size_t n) {
    if (price.size() == n) return;
    price.resize(n);
    delta.resize(n);
    gamma.resize(n);
    vega.resize(n);
    theta.resize(n);
    rho.resize(n);
}

void OptionChain::reserve(std::size_t n) {
    s_.reserve(n);
    k_.reserve(n);
    t_.reserve(n);
    r_.reserve(n);
    sigma_.reserve(n);
    q_.reserve(n);
    type_.reserve(n);
    style_.reserve(n);
}

void OptionChain::add(double S, double K, double T, double r, double sigma,
                      double q, OptionType type, ExerciseStyle style) {
    s_.push_back(S);
    k_.push_back(K);
    t_.push_back(T);
    r_.push_back(r);
    sigma_.push_back(sigma);
    q_.push_back(q);
    type_.push_back(type);
    style_.push_back(style);
}

void OptionChain::add(const Option& opt) {
    add(opt.S, opt.K, opt.T, opt.r, opt.sigma, opt.q, opt.type, opt.style);
}

void OptionChain::clear() noexcept {
    s_.clear();
    k_.clear();
    t_.clear();
    r_.clear();
    sigma_.clear();
    q_.clear();
    type_.clear();
    style_.clear();
}

const ChainResults& OptionChain::price_european_batch() {
    const std::size_t n = s_.size();
    results_.resize(n); // one-time (or on-growth) allocation; reused thereafter

    // Hoist raw pointers so the compiler sees independent linear streams.
    const double* __restrict S = s_.data();
    const double* __restrict K = k_.data();
    const double* __restrict T = t_.data();
    const double* __restrict R = r_.data();
    const double* __restrict V = sigma_.data();
    const double* __restrict Q = q_.data();
    const OptionType* __restrict TY = type_.data();

    double* __restrict pr = results_.price.data();
    double* __restrict de = results_.delta.data();
    double* __restrict ga = results_.gamma.data();
    double* __restrict ve = results_.vega.data();
    double* __restrict th = results_.theta.data();
    double* __restrict rh = results_.rho.data();

    for (std::size_t i = 0; i < n; ++i) {
        const PriceGreeks g =
            black_scholes(S[i], K[i], T[i], R[i], V[i], Q[i], TY[i]);
        pr[i] = g.price;
        de[i] = g.delta;
        ga[i] = g.gamma;
        ve[i] = g.vega;
        th[i] = g.theta;
        rh[i] = g.rho;
    }

    return results_;
}

} // namespace pricing
