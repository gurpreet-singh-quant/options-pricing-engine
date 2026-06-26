/// \file test_numerics.cpp
/// \brief Stability and correctness of the normal CDF/PDF and safe_exp.

#include "pricing/numerics.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <limits>

using Catch::Approx;
using namespace pricing::numerics;

TEST_CASE("norm_cdf known values", "[numerics][cdf]") {
    REQUIRE(norm_cdf(0.0) == Approx(0.5).margin(1e-12));
    REQUIRE(norm_cdf(1.0) == Approx(0.8413447460685429).margin(1e-12));
    REQUIRE(norm_cdf(-1.0) == Approx(0.15865525393145705).margin(1e-12));
    REQUIRE(norm_cdf(1.959963984540054) == Approx(0.975).margin(1e-12));
    REQUIRE(norm_cdf(-1.959963984540054) == Approx(0.025).margin(1e-12));
}

TEST_CASE("norm_cdf is monotone, bounded, and symmetric", "[numerics][cdf]") {
    double prev = norm_cdf(-50.0);
    for (double x = -50.0; x <= 50.0; x += 0.1) {
        const double v = norm_cdf(x);
        REQUIRE(std::isfinite(v));
        REQUIRE(v >= 0.0);
        REQUIRE(v <= 1.0);
        REQUIRE(v >= prev - 1e-15); // non-decreasing (allow last-ulp)
        prev = v;
        // Reflection symmetry Phi(-x) = 1 - Phi(x).
        REQUIRE(norm_cdf(-x) == Approx(1.0 - v).margin(1e-12));
    }
}

TEST_CASE("norm_cdf is accurate deep in the tails", "[numerics][cdf][stability]") {
    // The whole point of the erfc-based formula is that the tail probability is
    // computed *accurately* (not flushed to 0 or approximated coarsely). These
    // reference values would be unreachable by a naive (1+erf)/2 or a low-order
    // polynomial fit. Tolerances are RELATIVE.
    REQUIRE(norm_cdf(-5.0) == Approx(2.8665157187919391e-07).epsilon(1e-9));
    REQUIRE(norm_cdf(-6.0) == Approx(9.8658764503769814e-10).epsilon(1e-9));
    REQUIRE(norm_cdf(-8.0) == Approx(6.2209605742718204e-16).epsilon(1e-6));
    REQUIRE(norm_cdf(-10.0) == Approx(7.6198530241605260e-24).epsilon(1e-5));
    // Right-tail complement stays accurate too.
    REQUIRE(norm_cdf(6.0) == Approx(1.0 - 9.8658764503769814e-10).epsilon(1e-15));
}

TEST_CASE("norm_cdf saturates cleanly in the tails", "[numerics][cdf][stability]") {
    REQUIRE(norm_cdf(40.0) == Approx(1.0).margin(1e-12));
    REQUIRE(norm_cdf(-40.0) == 0.0);          // underflows to exactly 0
    REQUIRE(std::isfinite(norm_cdf(1e300)));
    REQUIRE(std::isfinite(norm_cdf(-1e300)));
    REQUIRE(norm_cdf(1e300) == 1.0);
    REQUIRE(norm_cdf(-1e300) == 0.0);
}

TEST_CASE("norm_pdf known values and tails", "[numerics][pdf]") {
    REQUIRE(norm_pdf(0.0) == Approx(0.3989422804014327).margin(1e-15));
    REQUIRE(norm_pdf(1.0) == Approx(0.24197072451914337).margin(1e-15));
    REQUIRE(norm_pdf(-1.0) == Approx(norm_pdf(1.0)).margin(1e-15)); // even
    // Deep tail underflows to a finite (zero) value, never inf/NaN.
    REQUIRE(std::isfinite(norm_pdf(1e6)));
    REQUIRE(norm_pdf(1e6) == 0.0);
    REQUIRE(norm_pdf(-1e6) == 0.0);
}

TEST_CASE("safe_exp never overflows and underflows to zero", "[numerics][exp][stability]") {
    REQUIRE(safe_exp(0.0) == Approx(1.0).margin(1e-15));
    REQUIRE(safe_exp(1.0) == Approx(std::exp(1.0)).margin(1e-12));
    // Above the natural overflow threshold: finite, not +inf.
    REQUIRE(std::isfinite(safe_exp(1000.0)));
    REQUIRE(std::isfinite(safe_exp(1e9)));
    REQUIRE(safe_exp(1000.0) > 0.0);
    // Far negative: exactly 0, not a denormal or NaN.
    REQUIRE(safe_exp(-1000.0) == 0.0);
    REQUIRE(safe_exp(-1e9) == 0.0);
}

TEST_CASE("is_finite helper", "[numerics]") {
    REQUIRE(is_finite(1.0));
    REQUIRE_FALSE(is_finite(std::numeric_limits<double>::infinity()));
    REQUIRE_FALSE(is_finite(-std::numeric_limits<double>::infinity()));
    REQUIRE_FALSE(is_finite(std::numeric_limits<double>::quiet_NaN()));
}
