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
fixed_analytics = "4.6.0"
```

For `no_std` environments:

```toml
[dependencies]
fixed_analytics = { version = "4.6.0", default-features = false }
```

## Available Functions

### Function Categories

**Total functions** return `T` directly and handle all inputs, possibly with saturation.
**Fallible functions** return `Result<T, Error>` and fail on domain violations.

| Category | Total Functions | Fallible Functions |
|----------|-----------------|-------------------|
| Trigonometric | `sin`, `cos`, `tan`, `sin_cos`, `atan`, `atan2` | `asin`, `acos`, `asin_from_gap` |
| Hyperbolic | `sinh`, `cosh`, `tanh`, `sinh_cosh`, `asinh` | `acosh`, `atanh`, `acoth`, `coth` |
| Exponential | `exp`, `pow2`, `expm1`, `exprel`, `exprel2`, `exp_moments` | `ln`, `log2`, `log10`, `pow`, `log1p`, `log1p_ratio`, `log1p_residual`, `log1p_residual_ratio` |
| Turns | `sin_cos_turns`, `atan2_turns`, `div_by_tau`, `mul_by_tau`, `wrap_radians`, `wrap_principal` | — |
| Ratios | `stumpff`, `tan_ratio`, `atan_ratio`, `tanh_ratio` | `atanh_ratio` |
| Algebraic | `sqrt_or_zero`, `hypot`, `norm`, `scale_exponent`, `scale_by` | `sqrt`, `sqrt_product`, `rsqrt`, `normalize` |
| Special | `erf`, `erfc`, `fresnel` | `gamma` |
| Complex (`Complex`) | `exp`, `exprel`, `exprel2`, `faddeeva`, `erf` | `ln`, `pow`, `exprel_weighted`, `e1`, `gamma_upper` |
| Elliptic | `ellip_e`, `jacobi_sn_cn_dn`, `jacobi_am` | `ellip_k`, `ellip_pi`, `ellip_pi_complete`, `jacobi_inverse`, `carlson_rf`, `carlson_rf_roots`, `carlson_rj`, `carlson_rj_roots`, `carlson_rc`, `carlson_rd` |
| Arithmetic | `mul_add`, `mul_div`, `div_or_saturate`, `dot` | — |
| Solvers | — | `root_in`, `kepler_universal` |
| Enclosures | `enclose`, `sin_cos_range`, `sin_cos_range_turns`, `exp_range`, `atan_range`, `tanh_range` | — |

Each function works in binary floating point at twice its type's width or more, from its arguments converted exactly, and rounds once to the type. Every one is within one ulp of the exact result on every type, `ERROR_ULPS`, and the arithmetic functions within half an ulp. Odd and even functions are exactly so and monotone functions stay monotone. Complete absence of panic is verified at the linker level via the [`no-panic`](https://github.com/dtolnay/no-panic) crate.

### Saturation Behavior

The following total functions saturate, clamping to the representable range near the following thresholds.

| Function | I16F16 Threshold | I32F32 Threshold | Result |
|----------|------------------|------------------|--------|
| `exp` | x ≥ 10.4 | x ≥ 21.5 | `T::MAX` |
| `exp` | x ≤ -11.8 | x ≤ -22.9 | Zero |
| `pow2` | x ≥ 15.0 | x ≥ 31.0 | `T::MAX` |
| `pow2` | x ≤ -17.0 | x ≤ -33.0 | Zero |
| `sinh` | \|x\| ≥ 11.1 | \|x\| ≥ 22.2 | `T::MAX` or `T::MIN` |
| `cosh` | \|x\| ≥ 11.1 | \|x\| ≥ 22.2 | `T::MAX` |
| `tan` | \|x - pole\| < 3.1e-5 | \|x - pole\| < 4.7e-10 | `T::MAX` or `T::MIN` |

Where for `tan`, "pole" refers to ±π/2, ±3π/2, ±5π/2, ...

<!-- ACCURACY_START -->
### Accuracy

Relative error statistics measured against MPFR reference implementations. Accuracy regressions are not permitted; every change is benchmarked against the baseline before merging. The file tools/accuracy-bench/baseline.json contains further measurements.

| Function | I16F16 Mean | I16F16 Median | I16F16 P95 | I32F32 Mean | I32F32 Median | I32F32 P95 |
|----------|-------------|---------------|------------|-------------|---------------|------------|
| sin | 3.11e-5 | 6.12e-6 | 6.18e-5 | 4.65e-10 | 9.35e-11 | 9.59e-10 |
| cos | 3.14e-5 | 6.17e-6 | 6.29e-5 | 4.73e-10 | 9.33e-11 | 9.89e-10 |
| tan | 3.28e-5 | 1.26e-5 | 7.51e-5 | 5.37e-10 | 1.91e-10 | 1.13e-9 |
| asin | 3.89e-5 | 1.12e-5 | 9.52e-5 | 6.56e-10 | 1.72e-10 | 1.50e-9 |
| acos | 8.39e-6 | 3.69e-6 | 2.85e-5 | 1.29e-10 | 5.64e-11 | 4.25e-10 |
| atan | 2.72e-6 | 2.52e-6 | 4.78e-6 | 4.19e-11 | 3.88e-11 | 7.29e-11 |
| sinh | 7.41e-6 | 4.17e-6 | 1.26e-5 | 1.18e-10 | 6.33e-11 | 1.91e-10 |
| cosh | 3.79e-6 | 3.73e-6 | 7.35e-6 | 5.76e-11 | 5.65e-11 | 1.12e-10 |
| tanh | 5.38e-6 | 2.10e-6 | 9.51e-6 | 1.09e-10 | 6.18e-11 | 1.41e-10 |
| coth | 3.53e-6 | 2.24e-6 | 8.26e-6 | 7.05e-11 | 5.99e-11 | 1.28e-10 |
| asinh | 2.85e-6 | 1.37e-6 | 4.96e-6 | 4.43e-11 | 2.08e-11 | 7.56e-11 |
| acosh | 2.29e-6 | 1.36e-6 | 4.80e-6 | 3.59e-11 | 2.09e-11 | 7.45e-11 |
| atanh | 4.22e-5 | 1.40e-5 | 1.01e-4 | 7.11e-10 | 2.13e-10 | 1.56e-9 |
| acoth | 1.92e-4 | 1.40e-4 | 5.49e-4 | 2.94e-9 | 2.16e-9 | 8.40e-9 |
| exp | 5.98e-3 | 9.19e-6 | 4.13e-2 | 9.45e-8 | 1.39e-10 | 6.50e-7 |
| ln | 1.02e-6 | 6.57e-7 | 1.46e-6 | 1.39e-11 | 1.01e-11 | 2.23e-11 |
| log2 | 6.76e-7 | 4.55e-7 | 1.01e-6 | 9.69e-12 | 6.89e-12 | 1.53e-11 |
| log10 | 1.92e-6 | 1.50e-6 | 3.35e-6 | 3.00e-11 | 2.30e-11 | 5.07e-11 |
| pow2 | 3.54e-4 | 5.03e-6 | 2.37e-3 | 5.47e-9 | 7.69e-11 | 3.67e-8 |
| pow | 6.47e-4 | 1.33e-5 | 3.13e-3 | 1.01e-8 | 2.02e-10 | 4.85e-8 |
| sqrt | 8.88e-8 | 5.80e-8 | 2.42e-7 | 1.37e-12 | 8.85e-13 | 3.62e-12 |
<!-- ACCURACY_END -->

<!-- WORST_START -->
### Worst-case error

The largest error seen, in ulps of the type, against MPFR on exact inputs. The sweep covers 104 layouts from 8 to 128 bits, every value of each 8- and 16-bit one and 40000 seeded samples of each wider one; tools/accuracy-bench/worst.txt holds every layout.

| Function | I16F16 | I32F32 | I64F64 |
|----------|------|------|------|
| sin | 0.50 | 0.50 | 0.50 |
| cos | 0.50 | 0.50 | 0.50 |
| tan | 0.50 | 0.50 | 0.50 |
| asin | 0.50 | 0.50 | 0.50 |
| acos | 0.50 | 0.50 | 0.50 |
| atan | 0.50 | 0.50 | 0.50 |
| atan2 | 0.50 | 0.50 | 0.50 |
| sinh | 0.50 | 0.50 | 0.50 |
| cosh | 0.50 | 0.50 | 0.50 |
| tanh | 0.50 | 0.50 | 0.50 |
| coth | 0.50 | 0.50 | 0.50 |
| asinh | 0.50 | 0.50 | 0.50 |
| acosh | 0.50 | 0.50 | 0.50 |
| atanh | 0.50 | 0.50 | 0.50 |
| acoth | 0.50 | 0.50 | 0.50 |
| exp | 0.50 | 0.50 | 0.50 |
| pow2 | 0.50 | 0.50 | 0.50 |
| ln | 0.50 | 0.50 | 0.50 |
| log2 | 0.50 | 0.50 | 0.50 |
| log10 | 0.50 | 0.50 | 0.50 |
| pow | 0.50 | 0.50 | 0.50 |
| sqrt | 0.50 | 0.50 | 0.50 |
| sqrt_or_zero | 0.50 | 0.50 | 0.50 |
| mul_add | 0.50 | 0.50 | 0.50 |
| mul_div | 0.50 | 0.50 | 0.50 |
| div_or_saturate | 0.50 | 0.50 | 0.50 |
| dot3 | 0.50 | 0.50 | 0.50 |
| hypot | 0.50 | 0.50 | 0.50 |
| norm3 | 0.50 | 0.50 | 0.50 |
| normalize3 | 0.50 | 0.50 | 0.50 |
| scale_by | 0.50 | 0.50 | 0.50 |
| sqrt_product | 0.50 | 0.50 | 0.50 |
| asin_from_gap | 0.50 | 0.50 | 0.50 |
| rsqrt | 0.50 | 0.50 | 0.50 |
| expm1 | 0.50 | 0.50 | 0.50 |
| exprel | 0.50 | 0.50 | 0.50 |
| exprel2 | 0.50 | 0.50 | 0.50 |
| exp_moments_e1 | 0.50 | 0.50 | 0.50 |
| exp_moments_e2 | 0.50 | 0.50 | 0.50 |
| exp_moments_e3 | 0.50 | 0.50 | 0.50 |
| exp_moments_e4 | 0.50 | 0.50 | 0.50 |
| exp_moments_e5 | 0.50 | 0.50 | 0.50 |
| log1p | 0.50 | 0.50 | 0.50 |
| log1p_ratio | 0.50 | 0.50 | 0.50 |
| log1p_residual | 0.50 | 0.50 | 0.50 |
| log1p_residual_ratio | 0.50 | 0.50 | 0.50 |
| stumpff_c0 | 0.50 | 0.50 | 0.50 |
| stumpff_c1 | 0.50 | 0.50 | 0.50 |
| stumpff_c2 | 0.50 | 0.50 | 0.50 |
| stumpff_c3 | 0.50 | 0.50 | 0.50 |
| tan_ratio | 0.50 | 0.50 | 0.50 |
| atan_ratio | 0.50 | 0.50 | 0.50 |
| atanh_ratio | 0.50 | 0.50 | 0.50 |
| tanh_ratio | 0.50 | 0.50 | 0.50 |
| ellip_k | 0.50 | 0.50 | 0.50 |
| ellip_k_from_m | 0.50 | 0.50 | 0.50 |
| ellip_e | 0.50 | 0.50 | 0.50 |
| jacobi_sn | 0.50 | 0.50 | 0.50 |
| jacobi_cn | 0.50 | 0.50 | 0.50 |
| jacobi_dn | 0.50 | 0.50 | 0.50 |
| jacobi_am | 0.50 | 0.50 | 0.50 |
| jacobi_inverse | 0.50 | 0.50 | 0.50 |
| carlson_rf | 0.50 | 0.50 | 0.50 |
| carlson_rf_roots | 0.50 | 0.50 | 0.50 |
| carlson_rj | 0.50 | 0.50 | 0.50 |
| carlson_rc | 0.50 | 0.50 | 0.50 |
| carlson_rd | 0.50 | 0.50 | 0.50 |
| ellip_pi | 0.50 | 0.50 | 0.50 |
| ellip_pi_complete | 0.50 | 0.50 | 0.50 |
| sin_turns | 0.50 | 0.50 | 0.50 |
| cos_turns | 0.50 | 0.50 | 0.50 |
| atan2_turns | 0.50 | 0.50 | 0.50 |
| div_by_tau | 0.50 | 0.50 | 0.50 |
| wrap_radians | 0.50 | 0.50 | 0.50 |
| wrap_principal | 0.50 | 0.50 | 0.50 |
| mul_by_tau | 0.50 | 0.50 | 0.50 |
| c_mul_re | 0.50 | 0.50 | 0.50 |
| c_mul_im | 0.50 | 0.50 | 0.50 |
| c_div_re | 0.50 | 0.50 | 0.50 |
| c_div_im | 0.50 | 0.50 | 0.50 |
| c_exp_re | 0.50 | 0.50 | 0.50 |
| c_exp_im | 0.50 | 0.50 | 0.50 |
| c_exprel_re | 0.50 | 0.50 | 0.50 |
| c_exprel_im | 0.50 | 0.50 | 0.50 |
| c_exprel2_re | 0.50 | 0.50 | 0.50 |
| c_exprel2_im | 0.50 | 0.50 | 0.50 |
| faddeeva_re | 0.50 | 0.50 | 0.50 |
| faddeeva_im | 0.50 | 0.50 | 0.50 |
| c_erf_re | 0.50 | 0.50 | 0.50 |
| c_erf_im | 0.50 | 0.50 | 0.50 |
| e1_re | 0.50 | 0.50 | 0.50 |
| e1_im | 0.50 | 0.50 | 0.50 |
| gamma_upper_re | 0.50 | 0.50 | 0.50 |
| gamma_upper_im | 0.50 | 0.50 | 0.50 |
| c_ln_re | 0.50 | 0.50 | 0.50 |
| c_ln_im | 0.50 | 0.50 | 0.50 |
| c_pow_re | 0.50 | 0.50 | 0.50 |
| c_pow_im | 0.50 | 0.50 | 0.50 |
| exprel_weighted_k1_re | 0.50 | 0.50 | 0.50 |
| exprel_weighted_k1_im | 0.50 | 0.50 | 0.49 |
| exprel_weighted_k2_re | 0.50 | 0.50 | 0.49 |
| exprel_weighted_k2_im | 0.50 | 0.50 | 0.50 |
| kepler_universal | 0.50 | 0.50 | 0.50 |
| erf | 0.50 | 0.50 | 0.50 |
| erfc | 0.50 | 0.50 | 0.50 |
| fresnel_c | 0.50 | 0.50 | 0.50 |
| fresnel_s | 0.50 | 0.50 | 0.50 |
| gamma | 0.50 | 0.50 | 0.50 |
<!-- WORST_END -->