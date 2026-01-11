//! Fixed-point math via CORDIC. No floating-point ops, `no_std` compatible.
//!
//! # Quick Start
//!
//! ```rust
//! use fixed::types::I16F16;
//! use fixed_analytics::{sin, cos, sqrt, ln};
//!
//! let angle = I16F16::from_num(0.5);
//! let (s, c) = (sin(angle), cos(angle));
//! let root = sqrt(I16F16::from_num(2.0)).unwrap();
//! let log = ln(I16F16::E).unwrap();
//! ```
//!
//! # Precision
//!
//! | Type | Accuracy |
//! |------|----------|
//! | `I16F16` | ~4 decimal digits |
//! | `I32F32` | ~8 decimal digits |
//!
//! # Features
//!
//! - `std` (default): Enables `std::error::Error` impl
//!
//! See [`kernel`] module for algorithm details.

#![no_std]
#![cfg_attr(docsrs, feature(doc_cfg))]

pub mod bounded;
pub mod error;
pub mod kernel;
pub mod ops;
pub mod tables;
pub mod traits;

// Re-export the fixed crate for convenience
pub use fixed;

// Re-export main types
pub use error::{Error, Result};
pub use traits::CordicNumber;

// Re-export all mathematical functions at crate root for convenience
pub use ops::algebraic::sqrt;
pub use ops::circular::{acos, asin, atan, atan2, cos, sin, sin_cos, tan};
pub use ops::exponential::{exp, ln, log2, log10, pow2};
pub use ops::hyperbolic::{acosh, acoth, asinh, atanh, cosh, coth, sinh, sinh_cosh, tanh};
