#![doc = include_str!("../docs/crate.md")]
#![no_std]
#![cfg_attr(docsrs, feature(doc_cfg))]

pub mod bounded;
pub mod error;
pub mod kernel;
pub mod ops;
pub mod tables;
pub mod traits;

pub use fixed;

pub use error::{Error, Result};
pub use traits::CordicNumber;

pub use ops::algebraic::sqrt;
pub use ops::circular::{acos, asin, atan, atan2, cos, sin, sin_cos, tan};
pub use ops::exponential::{exp, ln, log2, log10, pow, pow2};
pub use ops::hyperbolic::{acosh, acoth, asinh, atanh, cosh, coth, sinh, sinh_cosh, tanh};
