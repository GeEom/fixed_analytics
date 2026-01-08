# fixed_analytics

Fixed-point trigonometric, hyperbolic, exponential, and algebraic functions via CORDIC.

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
fixed_analytics = "0.3"
```

For `no_std` environments:

```toml
[dependencies]
fixed_analytics = { version = "0.3", default-features = false }
```

## Available Functions

| Category | Total Functions | Fallible Functions |
|----------|-----------------|-------------------|
| Trigonometric | `sin`, `cos`, `tan`, `sin_cos`, `atan`, `atan2` | `asin`, `acos` |
| Hyperbolic | `sinh`, `cosh`, `tanh`, `sinh_cosh`, `asinh` | `acosh`, `atanh`, `acoth`, `coth` |
| Exponential | `exp`, `pow2` | `ln`, `log2`, `log10` |
| Algebraic | — | `sqrt` |

Fallible functions return `Result<T, Error>` and fail on domain violations.
Total functions return `T` directly and handle all inputs.
