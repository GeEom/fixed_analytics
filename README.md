# fixed_analytics

Fixed-point mathematical functions which are accurate, deterministic, and guaranteed not to panic.

[![Crates.io](https://img.shields.io/crates/v/fixed_analytics.svg)](https://crates.io/crates/fixed_analytics)
[![CI](https://github.com/GeEom/fixed_analytics/actions/workflows/ci.yml/badge.svg)](https://github.com/GeEom/fixed_analytics/actions/workflows/ci.yml)
[![unsafe forbidden](https://img.shields.io/badge/unsafe-forbidden-success.svg)](https://github.com/rust-secure-code/safety-dance/)
[![MIT licensed](https://img.shields.io/badge/license-MIT-blue.svg)](./LICENSE.txt)
[![codecov](https://codecov.io/gh/GeEom/fixed_analytics/branch/main/graph/badge.svg)](https://codecov.io/gh/GeEom/fixed_analytics)

Every function is generic over the signed types of the [`fixed`](https://crates.io/crates/fixed) crate, from `I8F0` to `I3F125`. None uses floating point or allocates, and all run under `no_std`. The [`no-panic`](https://github.com/dtolnay/no-panic) crate checks at link time that no function can panic.

## Example

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
fixed_analytics = "4.0.0"
```

For `no_std`, turn off the default `std` feature, which only implements `std::error::Error` for `Error`:

```toml
[dependencies]
fixed_analytics = { version = "4.0.0", default-features = false }
```

## Functions

Total functions return `T` for every input. Fallible functions return `Result<T, Error>` and fail only outside their domain.

| Category | Total | Fallible |
|----------|-------|----------|
| Trigonometric | `sin`, `cos`, `tan`, `sin_cos`, `atan`, `atan2` | `asin`, `acos`, `asin_from_gap` |
| Hyperbolic | `sinh`, `cosh`, `tanh`, `sinh_cosh`, `asinh` | `acosh`, `atanh`, `acoth`, `coth` |
| Exponential | `exp`, `pow2`, `expm1`, `exprel`, `exprel2`, `exp_moments` | `ln`, `log2`, `log10`, `pow`, `log1p`, `log1p_ratio`, `log1p_residual`, `log1p_residual_ratio` |
| Turns | `sin_cos_turns`, `atan2_turns`, `div_by_tau`, `mul_by_tau`, `wrap_radians`, `wrap_principal` | — |
| Ratios | `stumpff`, `sin_ratio`, `tan_ratio`, `atan_ratio`, `tanh_ratio` | `atanh_ratio` |
| Algebraic | `sqrt_or_zero`, `hypot`, `norm`, `scale_exponent`, `scale_by` | `sqrt`, `sqrt_product`, `rsqrt`, `normalize` |
| Special | `erf`, `erfc`, `fresnel` | `gamma` |
| Complex (`Complex`) | `exp`, `exprel`, `exprel2`, `faddeeva`, `erf` | `ln`, `pow`, `exprel_weighted`, `e1`, `gamma_upper` |
| Elliptic | `ellip_e`, `jacobi_sn_cn_dn`, `jacobi_am` | `ellip_k`, `ellip_pi`, `ellip_pi_complete`, `ellip_pi_increment`, `jacobi_inverse`, `carlson_rf`, `carlson_rf_roots`, `carlson_rj`, `carlson_rj_roots`, `carlson_rc`, `carlson_rd`, `carlson_rd_roots`, `carlson_rf_rd_roots`, `carlson_rd_rj_roots` |
| Arithmetic | `mul_add`, `mul_div`, `div_or_saturate`, `dot`, and by `Rounding` `mul_rounded`, `div_rounded`, `mul_div_rounded`, `mul_add_rounded`, `dot_rounded` | — |
| Solvers | — | `root_in`, `root_in_newton`, `kepler_universal` |
| Enclosures | `enclose`, `sin_cos_range`, `sin_cos_range_turns`, `exp_range`, `atan_range`, `tanh_range`, `exprel_range`, `exprel2_range`, `exp_moments_range`, `tanh_ratio_range`, `atan_ratio_range`, `tan_ratio_range`, `stumpff_range`, `sqrt_range`, `hypot_range`, `atan2_range` | `log1p_range`, `log1p_ratio_range`, `log1p_residual_ratio_range`, `atanh_ratio_range` |

## Accuracy

Each function evaluates in binary floating point at twice its type's width or more, from its arguments converted exactly, and rounds once. Every result lies within one ulp of the exact value on every layout, the bound `ERROR_ULPS` states. The arithmetic does better. `mul_add`, `mul_div`, `div_or_saturate`, `dot` and `sqrt` round the exact result to nearest, and the `_rounded` forms round it toward −∞ or +∞. Odd and even functions keep their symmetry, and monotone functions stay monotone.

The accuracy tool checks every function against MPFR on 104 layouts from 8 to 128 bits. On the 8- and 16-bit layouts it tries every input of a one-argument function; elsewhere it draws seeded inputs. Up to the working error at a midpoint, it has found the nearest value in all but two functions. `wrap_principal` keeps its result in (−π, π] at a cost of up to 0.98 ulp next to ±π. The real part of `Complex::erf` reaches 0.80 ulp on I80F48. No symmetry or monotonicity check has failed. CI compares every change with the baseline and fails on a regression.

## Saturation

Out of range, a total function returns `T::MAX` or `T::MIN`. A result under half an ulp rounds to zero. Dividing by zero gives the end with the numerator's sign, and 0/0 gives 0. On I16F16, `exp` saturates from 10.40 and `pow2` from 15, `sinh` and `cosh` from ±11.09, and `tan` within 3.1e-5 of a pole. Inside their domains the fallible functions saturate too, `pow` and `gamma` among them.

## Documentation

- [Accuracy](docs/accuracy.md): the bound, how it is measured, and relative error per function.
- [Worst case](docs/worst-case.md): the largest error of every function, on all 104 layouts.
- [Saturation](docs/saturation.md): thresholds and domains for every function.
- [API reference](https://docs.rs/fixed_analytics) and [changelog](CHANGELOG.md).
