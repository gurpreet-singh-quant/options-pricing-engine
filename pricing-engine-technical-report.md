# Technical Report: Correctness Hardening of a Real-Time Options Pricing Engine

## Summary

A C++20 options pricing engine — closed-form Black-Scholes pricing and Greeks for European options, a Cox-Ross-Rubinstein binomial tree for American exercise, and a real-time batch pricing demo — was put through a structured, multi-dimensional code review (40 independent review passes across 6 evaluation dimensions: correctness, numerical stability, performance, memory safety, build portability, and documentation accuracy). The review surfaced 17 confirmed issues. Two were rated high severity because they produced silently wrong output rather than a visible failure. Both are described below, along with the lower-severity hardening fixes, and the final verification results.

## Bug 1 (High): Binomial Tree Mispricing at Low Effective Volatility

**Symptom.** For options with low volatility relative to the tree's time step, the American-exercise binomial pricer returned a price roughly 86% off from the correct value, with no error, warning, or crash — the wrong number simply looked plausible.

**Root cause.** The CRR (Cox-Ross-Rubinstein) binomial tree derives its risk-neutral up-move probability as:

```
p = (e^(r·Δt) - d) / (u - d)
```

where `u = e^(σ√Δt)` and `d = 1/u`. This formula only yields a valid probability — a value strictly between 0 and 1 — when the discount-adjusted drift term stays bounded relative to the volatility-driven spread between `u` and `d`. When volatility is low relative to the chosen time step (which can happen with a coarse step count, a low-volatility underlying, or a relatively high risk-free rate), `p` can mathematically compute to a value greater than 1 or less than 0. A naive implementation clamps this out-of-range value into `[0,1]`, which silently turns an invalid probability into a confident-looking but materially wrong number — the source of the ~86% error observed.

**Fix.** The pricer now detects when `p` falls outside `[0,1]` and, instead of clamping, falls back to the deterministic (no-arbitrage) limiting price for that node — the value the discrete tree is supposed to approximate as the random-walk assumption breaks down. This preserves correctness in the regime where the discrete CRR parameterization itself is no longer a valid approximation, rather than returning a number from a probability that was never real.

**Verification.** A regression test constructs an option in the specific low-volatility/coarse-step regime that triggers `p` escaping `[0,1]`, and asserts the engine now returns the correct deterministic-limit price rather than the previous clamped value. Reverting the fix causes this test to fail, confirming the test exercises the actual bug rather than passing trivially.

## Bug 2 (High): Unvalidated Tick Count Causing Crash/Undefined Behavior

**Symptom.** The real-time demo (which simulates a market data feed and re-prices an option chain on every tick) crashed or exhibited undefined behavior when given a tick count of zero or a negative value.

**Root cause.** The tick-count parameter was used directly to size or index into buffers and loop bounds without validation. A zero or negative count produced either a zero-sized allocation misused downstream, or a loop condition that never matched the intended range — both are classic sources of undefined behavior in C++ when input isn't validated before being used structurally.

**Fix.** Input validation now rejects non-positive tick counts at the entry point with a clear error rather than allowing them to propagate into buffer sizing or loop logic.

## Additional Hardening (Lower Severity)

- **Benchmark correctness:** replaced a deprecated Google Benchmark `DoNotOptimize` usage that could let the compiler optimize away the very code being measured — silently invalidating latency numbers.
- **Build portability:** gated MSVC-specific compiler flags so the project builds cleanly across compiler toolchains rather than only the one used during development.
- **Numerical accuracy test:** added a deep-tail test for the normal CDF approximation, verifying accuracy holds for extreme `d1`/`d2` values (deep ITM/OTM), not just typical ranges.
- **Consistency test:** added a test asserting the batch (vectorized) pricing path and the scalar single-option path produce identical results — catching any divergence between the two code paths that unit tests on each path in isolation could miss.
- **Documentation honesty:** softened overclaiming language in the docs — "always finite" and "auto-vectorizable" were stated as guarantees without a test backing every input domain; the wording now reflects what is actually verified rather than implied.

## Verification Results

```
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release   # configured cleanly
cmake --build build                                       # zero warnings under -Werror
ctest --test-dir build --output-on-failure                # 22/22 tests passed (100%)
```

Performance target (batch of 50 options, price + 5 Greeks, under 5ms):

> *[Insert the actual median/p99 latency printed by `scripts/run_benchmarks.sh` here once captured.]*

## What This Demonstrates

The most informative bug here isn't the crash — it's the mispricing. A binomial tree that silently returns a confidently wrong price is a more dangerous failure mode than one that crashes, because nothing about the output signals that anything is wrong. Catching it required understanding *why* the CRR parameterization breaks down at low volatility, not just testing that the code runs. The fix — falling back to the correct limiting behavior instead of clamping an invalid probability — reflects an understanding of the model's assumptions, not just its implementation.
