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

Requires Rust 1.88 or later.

```toml
[dependencies]
fixed_analytics = "2.0.0"
```

For `no_std` environments:

```toml
[dependencies]
fixed_analytics = { version = "2.0.0", default-features = false }
```

## Available Functions

### Function Categories

**Total functions** return `T` directly and handle all inputs, possibly with saturation.
**Fallible functions** return `Result<T, Error>` and fail on domain violations.

| Category | Total Functions | Fallible Functions |
|----------|-----------------|-------------------|
| Trigonometric | `sin`, `cos`, `tan`, `sin_cos`, `atan`, `atan2` | `asin`, `acos` |
| Hyperbolic | `sinh`, `cosh`, `tanh`, `sinh_cosh`, `asinh` | `acosh`, `atanh`, `acoth`, `coth` |
| Exponential | `exp`, `pow2` | `ln`, `log2`, `log10` |
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
| asin | 2.87e-4 | 5.93e-5 | 6.46e-4 | 5.34e-9 | 8.82e-10 | 1.03e-8 |
| acos | 3.61e-5 | 2.18e-5 | 1.14e-4 | 5.37e-10 | 3.19e-10 | 1.71e-9 |
| atan | 2.71e-5 | 2.21e-5 | 6.29e-5 | 3.69e-10 | 2.92e-10 | 8.74e-10 |
| sinh | 9.80e-5 | 6.23e-5 | 2.79e-4 | 1.52e-9 | 9.64e-10 | 4.29e-9 |
| cosh | 9.40e-5 | 5.75e-5 | 2.77e-4 | 1.44e-9 | 8.90e-10 | 4.25e-9 |
| tanh | 1.60e-5 | 1.32e-5 | 2.56e-5 | 2.25e-10 | 1.22e-10 | 3.90e-10 |
| coth | 6.68e-6 | 3.54e-6 | 1.80e-5 | 1.41e-10 | 1.16e-10 | 2.74e-10 |
| asinh | 6.44e-4 | 4.83e-4 | 1.75e-3 | 1.03e-8 | 7.59e-9 | 2.85e-8 |
| acosh | 6.74e-4 | 5.21e-4 | 1.80e-3 | 1.05e-8 | 7.96e-9 | 2.88e-8 |
| atanh | 3.01e-4 | 5.90e-5 | 6.25e-4 | 6.68e-9 | 1.32e-9 | 1.44e-8 |
| acoth | 2.10e-3 | 1.33e-3 | 6.67e-3 | 4.26e-8 | 2.62e-8 | 1.39e-7 |
| exp | 1.14e-2 | 2.32e-5 | 7.88e-2 | 1.91e-7 | 1.73e-9 | 1.30e-6 |
| ln | 1.35e-5 | 8.76e-6 | 2.97e-5 | 4.50e-10 | 3.48e-10 | 9.17e-10 |
| log2 | 1.33e-5 | 8.48e-6 | 2.92e-5 | 3.46e-10 | 2.24e-10 | 7.21e-10 |
| log10 | 1.44e-5 | 9.28e-6 | 3.14e-5 | 4.49e-10 | 3.27e-10 | 9.07e-10 |
| pow2 | 7.21e-4 | 2.58e-5 | 4.72e-3 | 1.13e-8 | 4.93e-10 | 7.43e-8 |
| sqrt | 1.77e-7 | 1.16e-7 | 4.74e-7 | 2.70e-12 | 1.78e-12 | 7.16e-12 |
<!-- ACCURACY_END -->