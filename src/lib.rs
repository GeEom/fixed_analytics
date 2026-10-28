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

/// Every function's worst-case error, in ulps of the type it computes in.
///
/// A result lies strictly within one ulp of the exact value wherever that
/// value is in range, and is the nearer end of the range otherwise. The
/// [`ops::arith`] functions and [`sqrt`] round an exact value once, so they
/// stay within half an ulp.
pub const ERROR_ULPS: u32 = 1;

pub use error::{Error, Result};
pub use traits::CordicNumber;

pub use ops::algebraic::sqrt;
pub use ops::arith::{div_or_saturate, dot, mul_add, mul_div, sqrt_or_zero};
pub use ops::circular::{acos, asin, atan, atan2, cos, sin, sin_cos, tan};
pub use ops::exponential::{exp, ln, log2, log10, pow, pow2};
pub use ops::hyperbolic::{acosh, acoth, asinh, atanh, cosh, coth, sinh, sinh_cosh, tanh};
