# Real-Time Options Pricing Engine (C++20)

A low-latency options pricing engine that computes option prices and the five
primary Greeks in real time, with **numerical stability at the edges treated as
a first-class deliverable**. It provides:

- A closed-form **Black-Scholes-Merton** pricer for European options.
- **Analytical** Delta, Gamma, Vega, Theta, and Rho (derived in closed form, not
  bumped).
- A **Cox-Ross-Rubinstein binomial tree** for American options (default 500
  steps) with finite-difference Greeks.
- A cache-friendly **Structure-of-Arrays batch API** that prices a whole option
  chain (price + all 5 Greeks) in a single pass with no heap allocation in the
  hot path.
- A **real-time tick-feed simulation** that re-prices a ~50-option chain on every
  market tick and logs recomputation latency to stdout and CSV.

Single option (price + 5 Greeks) costs **~0.14 µs**; a 50-option chain costs
**~8.6 µs** — nearly three orders of magnitude under the 5 ms budget.

---

## Build

Requirements: a C++20 compiler, CMake ≥ 3.24, and a generator (Ninja recommended).
Catch2 v3 and Google Benchmark are pulled in automatically (see
[Dependencies](#dependencies)).

```bash
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

Release builds compile with `-O3 -march=native -Wall -Wextra -Wpedantic -Werror`
(MSVC: `/W4 /WX /O2`). The build is **warning-free**.

### Dependencies

`Catch2` (v3) and `Google Benchmark` are declared with CMake `FetchContent` so a
fresh checkout needs no manual install. The declarations use `FIND_PACKAGE_ARGS`,
so CMake first tries a locally installed copy (`find_package`) and only clones
from upstream if none is found — fast configures and offline builds when the
packages are already present, automatic download otherwise.

## Run the tests

```bash
ctest --test-dir build --output-on-failure
```

22 test cases cover numerics stability (incl. deep-tail CDF accuracy), pricing
accuracy vs. reference values, batch-vs-scalar value equality, analytical-vs-
finite-difference Greek cross-checks, edge-case stability, and a performance gate.

## Run the benchmarks

```bash
cmake --build build --target pricing_benchmarks
./build/pricing_benchmarks --benchmark_report_aggregates_only=true --benchmark_time_unit=us
# or, in one shot:
scripts/run_benchmarks.sh
```

## Run the real-time demo

```bash
./build/pricing_demo [num_ticks] [interval_ms]   # defaults: 1000 ticks, 1 ms
```

Each tick advances the underlying by a GBM step, re-prices the 50-option chain
(price + all 5 Greeks), and logs per-tick recompute latency to stdout and
`latency_log.csv`. A summary (min/mean/median/p99/max) prints at the end.

---

## Numerical stability techniques

Stability is the headline concern; here is exactly what is done and why.

### Standard normal CDF `Φ(x)` — `src/numerics.cpp`
Computed as `Φ(x) = ½·erfc(−x/√2)` using the standard-library `std::erfc`.
- **Tail accuracy on the whole real line.** `erfc` is built to stay accurate as
  its argument grows, which is precisely the deep-ITM/OTM region. Polynomial
  approximations (e.g. Abramowitz–Stegun 7.1.26) carry ~1e-7 error and degrade
  in the tails.
- **No catastrophic cancellation.** Computing `Φ` as `(1+erf)/2` loses all
  significant digits in the left tail (subtracting two near-equal numbers);
  `erfc(−x/√2)` evaluates the small tail probability directly.
- **Monotone by construction**, then clamped to `[0,1]` against last-ulp
  excursions. `std::erfc` is part of the C++ standard library, not a third-party
  math/finance package.

### Standard normal PDF `φ(x)`
`φ(x) = exp(−x²/2)/√(2π)`. The exponent is always ≤ 0, so it can only underflow
smoothly to `+0` for large `|x|` — never overflow.

### `d1`/`d2` without catastrophic cancellation — `src/black_scholes.cpp`
The total volatility `w = σ·√T` is computed once, then
```
d1 = (ln(S/K) + (r − q)·T) / w + ½·w
d2 = d1 − w
```
`d2` is formed by subtracting `w` directly, so `d1 − d2 == σ·√T` exactly and we
never subtract two large near-equal numbers.

### Degenerate / limiting branch
When the model degenerates — `T ≤ 0`, `σ ≤ 0`, vanishing `w < 1e-12`, or a
non-positive spot/strike — the diffusion term disappears and the underlying is
deterministic. We return the **deterministic forward value**:
```
F    = S·e^{(r−q)T}      (forward price; F = S at T = 0)
disc = e^{−rT}           (disc = 1 at T = 0)
price(call) = disc·max(F − K, 0)
price(put)  = disc·max(K − F, 0)
```
This single expression is correct at **both** limits and continuous with the
full BSM formula:
- As `T → 0` it reduces to the **undiscounted intrinsic value** `max(S−K,0)` /
  `max(K−S,0)` — the required limit at expiry.
- As `σ → 0` with `T > 0` it gives the exact zero-vol payoff value, with no
  spurious discontinuity at the threshold.

Greeks in this branch are the analytic one-sided limits (`gamma = vega = 0`;
`delta`/`theta`/`rho` the deterministic carry terms), all finite. The ATM kink
uses the half-delta convention.

### Overflow / underflow guards
All `exp()` calls go through `numerics::safe_exp`, which saturates the exponent
into the representable band (no `+inf` above ~709.78, exact `0` far negative).
Deep-ITM/OTM options push `d1`/`d2` to `±∞`; the CDF saturates cleanly to 1 or 0
and the PDF underflows to `+0`, so **every Greek stays finite**.

### Dedicated stability suite — `tests/test_numerical_stability.cpp`
Explicitly exercises `T = 0`, `σ = 0`, `S = 0`, deep ITM, deep OTM, plus a
6-dimensional stress grid (1,200 parameter combinations) asserting **no NaN/Inf
ever** and the documented limiting values.

---

## Performance

- **No heap allocation in the hot path.** The batch API writes into a results
  buffer that is sized once and reused; per-tick re-pricing allocates nothing.
- **Structure-of-Arrays layout.** Inputs are parallel arrays (one per field), so
  a batch price walks memory linearly with no pointer chasing. (The kernel cost
  is dominated by the transcendental calls — `erfc`/`exp`/`log` — which the
  compiler does not auto-vectorize; the SoA layout's win is the linear,
  cache-friendly, branch-light access pattern rather than SIMD of the math.)
- Scalar BSM kernel is `noexcept`, branch-light, and takes inputs by value.

### Benchmark results

Measured on **Intel Core i5-7300U @ 2.60 GHz** (2 cores / 4 threads; L1 32 KiB,
L2 256 KiB, L3 3 MiB), Windows 10, g++ 15.2.0, `-O3 -march=native`, Release.
Google Benchmark, 50 repetitions; times below are wall-clock ("Time"), with CPU
time ~10–20 % lower. The p99 tails reflect OS scheduling on a 2-core laptop.

| Benchmark                         | Median    | p99       | Throughput      |
|-----------------------------------|-----------|-----------|-----------------|
| Single option (price + 5 Greeks)  | 0.135 µs  | 0.178 µs  | 7.88 M opt/s    |
| Batch of 50 (price + 5 Greeks)    | 8.57 µs   | 12.3 µs   | 6.18 M opt/s    |
| Batch of 500 (price + 5 Greeks)   | 101 µs    | 161 µs    | 5.64 M opt/s    |
| American (CRR, 500 steps), price  | 210 µs    | 229 µs    | —               |

**Targets vs. actual**

| Target                                  | Budget   | Actual (median) | Margin   |
|-----------------------------------------|----------|-----------------|----------|
| Single option, price + 5 Greeks         | < 1 µs   | 0.135 µs        | ~7×      |
| Batch of 50, price + 5 Greeks           | < 5 ms   | 8.57 µs         | ~580×    |

The performance gate in `tests/test_performance.cpp` fails the build if the
batch-of-50 median exceeds 5 ms.

---

## Project layout

```
options-pricing-engine/
├── CMakeLists.txt
├── README.md
├── include/pricing/
│   ├── option.hpp           # Option struct, OptionType/ExerciseStyle enums, PriceGreeks
│   ├── numerics.hpp         # stable Φ/φ, safe_exp
│   ├── black_scholes.hpp    # European price + analytical Greeks
│   ├── binomial_tree.hpp    # CRR tree for American options
│   └── option_chain.hpp     # SoA batch container + batch pricing API
├── src/
│   ├── numerics.cpp
│   ├── black_scholes.cpp
│   ├── binomial_tree.cpp
│   ├── option_chain.cpp
│   └── main.cpp             # real-time tick-feed simulation + latency logging
├── tests/
│   ├── test_numerics.cpp
│   ├── test_pricing_accuracy.cpp
│   ├── test_greeks_finite_diff.cpp
│   ├── test_numerical_stability.cpp
│   └── test_performance.cpp
├── benchmarks/benchmark_main.cpp
└── scripts/run_benchmarks.sh
```

## Greek conventions

All Greeks are "raw" / per-unit (no 1% or per-day rescaling):

| Greek | Definition            | Unit                         |
|-------|-----------------------|------------------------------|
| Delta | `∂Price/∂S`           | per 1.0 of underlying        |
| Gamma | `∂²Price/∂S²`         | —                            |
| Vega  | `∂Price/∂σ`           | per 1.0 of vol (i.e. per 100%) |
| Theta | `∂Price/∂t = −∂Price/∂T` | per 1.0 year (negative = decay) |
| Rho   | `∂Price/∂r`           | per 1.0 of rate (i.e. per 100%) |

American-option Greeks are **finite-difference** estimates over re-priced trees
(centered differences, labeled as such), since no simple closed form exists.
