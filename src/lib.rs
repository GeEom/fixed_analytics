#![doc = include_str!("../docs/crate.md")]
#![no_std]
#![cfg_attr(docsrs, feature(doc_cfg))]

pub mod bounded;
mod consts;
pub mod error;
mod eval;
pub mod kernel;
pub mod ops;
mod real;
mod series;
pub mod tables;
pub mod traits;

pub use fixed;

/// Every function's worst-case error in ulps, strict, wherever the exact
/// value is in range; out of range a result saturates.
pub const ERROR_ULPS: u32 = 1;

pub use error::{Error, Result};
pub use traits::CordicNumber;

pub use ops::algebraic::{
    hypot, norm, normalize, rsqrt, scale_by, scale_exponent, sqrt, sqrt_product,
};
pub use ops::arith::{div_or_saturate, dot, mul_add, mul_div, sqrt_or_zero};
pub use ops::circular::{acos, asin, asin_from_gap, atan, atan2, cos, sin, sin_cos, tan};
pub use ops::exponential::{
    exp, exp_moments, expm1, exprel, exprel2, ln, log1p, log1p_ratio, log1p_residual,
    log1p_residual_ratio, log2, log10, pow, pow2,
};
pub use ops::hyperbolic::{acosh, acoth, asinh, atanh, cosh, coth, sinh, sinh_cosh, tanh};
pub use ops::interval::{atan_range, enclose, exp_range, sin_cos_range, tanh_range};
pub use ops::ratio::{atan_ratio, atanh_ratio, stumpff, tan_ratio, tanh_ratio};
