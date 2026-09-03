# fixed_analytics

Fixed-point mathematical functions which are accurate, deterministic, and guaranteed not to panic.

[![Crates.io](https://img.shields.io/crates/v/fixed_analytics.svg)](https://crates.io/crates/fixed_analytics)
[![CI](https://github.com/GeEom/fixed_analytics/actions/workflows/ci.yml/badge.svg)](https://github.com/GeEom/fixed_analytics/actions/workflows/ci.yml)
[![unsafe forbidden](https://img.shields.io/badge/unsafe-forbidden-success.svg)](https://github.com/rust-secure-code/safety-dance/)
[![MIT licensed](https://img.shields.io/badge/license-MIT-blue.svg)](./LICENSE.txt)
[![codecov](https://codecov.io/gh/GeEom/fixed_analytics/branch/main/graph/badge.svg)](https://codecov.io/gh/GeEom/fixed_analytics)

## Examples

```rust
use fixed::types::I16F16;
use fixed_analytics::{sin, cos, sqrt, ln};

let angle = I16F16::from_num(0.5);
let (s, c) = (sin(angle), cos(angle));

let root = sqrt(I16F16::from_num(2.0)).unwrap();
assert!((root.to_num::<f32>() - 1.414).abs() < 0.001);

let log = ln(I16F16::E).unwrap();
assert!((log.to_num::<f32>() - 1.0).abs() < 0.01);
```

## Installation

Requires Rust 1.95 or later.

```toml
[dependencies]
fixed_analytics = "3.1.0"
```

For `no_std` environments:

```toml
[dependencies]
fixed_analytics = { version = "3.1.0", default-features = false }
```

## Available Functions

### Function Categories

**Total functions** return `T` directly and handle all inputs, possibly with saturation.
**Fallible functions** return `Result<T, Error>` and fail on domain violations.

| Category | Total Functions | Fallible Functions |
|----------|-----------------|-------------------|
| Trigonometric | `sin`, `cos`, `tan`, `sin_cos`, `atan`, `atan2` | `asin`, `acos` |
| Hyperbolic | `sinh`, `cosh`, `tanh`, `sinh_cosh`, `asinh` | `acosh`, `atanh`, `acoth`, `coth` |
| Exponential | `exp`, `pow2` | `ln`, `log2`, `log10`, `pow` |
| Algebraic | — | `sqrt` |

Functions are calculated via polynomial evaluation, CORDIC, and Newton-Raphson techniques. Complete absence of panic is verified at the linker level via the [`no-panic`](https://github.com/dtolnay/no-panic) crate.

### Saturation Behavior

The following total functions saturate, clamping to the representable range near the following thresholds.

| Function | I16F16 Threshold | I32F32 Threshold | Result |
|----------|------------------|------------------|--------|
| `exp` | x ≥ 10.4 | x ≥ 21.5 | `T::MAX` |
| `exp` | x ≤ -11.1 | x ≤ -22.2 | Zero |
| `pow2` | x ≥ 15.0 | x ≥ 31.0 | `T::MAX` |
| `pow2` | x ≤ -16.1 | x ≤ -32.1 | Zero |
| `sinh` | \|x\| ≥ 11.1 | \|x\| ≥ 22.2 | `T::MAX` or `T::MIN` |
| `cosh` | \|x\| ≥ 11.1 | \|x\| ≥ 22.2 | `T::MAX` |
| `tan` | \|x - pole\| < 4e-5 | \|x - pole\| < 5e-10 | `T::MAX` or `T::MIN` |

Where for `tan`, "pole" refers to ±π/2, ±3π/2, ±5π/2, ...

<!-- ACCURACY_START -->
### Accuracy

Relative error statistics measured against MPFR reference implementations. Accuracy regressions are not permitted; every change is benchmarked against the baseline before merging. The file tools/accuracy-bench/baseline.json contains further measurements.

| Function | I16F16 Mean | I16F16 Median | I16F16 P95 | I32F32 Mean | I32F32 Median | I32F32 P95 |
|----------|-------------|---------------|------------|-------------|---------------|------------|
| sin | 6.06e-4 | 8.78e-5 | 1.28e-3 | 1.16e-8 | 1.68e-9 | 2.43e-8 |
| cos | 6.45e-4 | 9.03e-5 | 1.38e-3 | 1.22e-8 | 1.72e-9 | 2.64e-8 |
| tan | 7.20e-5 | 3.57e-5 | 2.20e-4 | 1.28e-9 | 3.98e-10 | 3.03e-9 |
| asin | 1.13e-4 | 2.80e-5 | 3.83e-4 | 3.68e-9 | 5.72e-10 | 6.40e-9 |
| acos | 1.79e-5 | 1.17e-5 | 4.81e-5 | 3.50e-10 | 2.11e-10 | 1.04e-9 |
| atan | 7.79e-6 | 6.49e-6 | 1.82e-5 | 1.88e-10 | 1.51e-10 | 4.38e-10 |
| sinh | 9.80e-5 | 6.23e-5 | 2.79e-4 | 1.52e-9 | 9.64e-10 | 4.29e-9 |
| cosh | 9.40e-5 | 5.75e-5 | 2.77e-4 | 1.44e-9 | 8.90e-10 | 4.25e-9 |
| tanh | 1.60e-5 | 1.32e-5 | 2.56e-5 | 2.25e-10 | 1.22e-10 | 3.90e-10 |
| coth | 6.68e-6 | 3.54e-6 | 1.80e-5 | 1.41e-10 | 1.16e-10 | 2.74e-10 |
| asinh | 2.42e-5 | 1.61e-5 | 5.08e-5 | 6.27e-10 | 5.01e-10 | 1.19e-9 |
| acosh | 1.87e-5 | 1.48e-5 | 4.65e-5 | 5.29e-10 | 4.75e-10 | 1.11e-9 |
| atanh | 2.21e-4 | 3.99e-5 | 3.33e-4 | 3.63e-9 | 7.20e-10 | 6.24e-9 |
| acoth | 1.21e-3 | 6.85e-4 | 4.09e-3 | 1.94e-8 | 1.23e-8 | 6.28e-8 |
| exp | 5.98e-3 | 1.47e-5 | 4.13e-2 | 9.49e-8 | 1.50e-9 | 6.50e-7 |
| ln | 1.18e-5 | 7.61e-6 | 2.09e-5 | 4.35e-10 | 3.78e-10 | 6.47e-10 |
| log2 | 9.98e-6 | 6.72e-6 | 2.01e-5 | 1.92e-10 | 1.31e-10 | 3.90e-10 |
| log10 | 1.25e-5 | 8.96e-6 | 2.36e-5 | 4.06e-10 | 3.57e-10 | 6.39e-10 |
| pow2 | 3.62e-4 | 2.24e-5 | 2.37e-3 | 5.64e-9 | 4.29e-10 | 3.67e-8 |
| pow | 6.90e-4 | 6.83e-5 | 3.13e-3 | 1.09e-8 | 1.23e-9 | 4.86e-8 |
| sqrt | 8.88e-8 | 5.80e-8 | 2.42e-7 | 1.37e-12 | 8.85e-13 | 3.62e-12 |
<!-- ACCURACY_END -->