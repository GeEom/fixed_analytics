# fixed_analytics

Fixed-point mathematical functions using the CORDIC algorithm. Designed for embedded systems and environments without hardware floating-point support.

[![Crates.io](https://img.shields.io/crates/v/fixed_analytics.svg)](https://crates.io/crates/fixed_analytics)
[![CI](https://github.com/GeEom/fixed_analytics/actions/workflows/ci.yml/badge.svg)](https://github.com/GeEom/fixed_analytics/actions/workflows/ci.yml)
[![Rust](https://github.com/GeEom/fixed_analytics/actions/workflows/rust.yml/badge.svg)](https://github.com/GeEom/fixed_analytics/actions/workflows/rust.yml)
[![unsafe forbidden](https://img.shields.io/badge/unsafe-forbidden-success.svg)](https://github.com/rust-secure-code/safety-dance/)
[![MIT licensed](https://img.shields.io/badge/license-MIT-blue.svg)](./LICENSE.txt)
[![codecov](https://codecov.io/gh/GeEom/fixed_analytics/branch/master/graph/badge.svg)](https://codecov.io/gh/GeEom/fixed_analytics)

## Examples

```rust
use fixed::types::I16F16;
use fixed_analytics::{sin, cos, sqrt, ln};

let angle = I16F16::from_num(0.5);
let (s, c) = (sin(angle), cos(angle));

let root = sqrt(I16F16::from_num(2.0));  // ≈ 1.414

let log = ln(I16F16::from_num(2.718)).unwrap();  // ≈ 1.0
```

## Installation

Requires Rust 1.88 or later.

```toml
[dependencies]
fixed_analytics = "0.1"
```

For `no_std` environments:

```toml
[dependencies]
fixed_analytics = { version = "0.1", default-features = false }
```

## Details

| Category | Functions |
|----------|-----------|
| Trigonometric | `sin`, `cos`, `tan`, `sin_cos`, `asin`, `acos`, `atan`, `atan2` |
| Hyperbolic | `sinh`, `cosh`, `tanh`, `coth`, `sinh_cosh`, `asinh`, `acosh`, `atanh`, `acoth` |
| Exponential | `exp`, `ln`, `log2`, `log10` |
| Algebraic | `sqrt` |
