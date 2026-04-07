//! Fixed-point mathematical functions: accurate, deterministic, and guaranteed not to panic.
//!
//! `no_std` compatible with no floating-point operations.
//!
//! # Quick Start
//!
//! ```rust
//! use fixed::types::I16F16;
//! use fixed_analytics::{sin, cos, sqrt, ln};
//!
//! let angle = I16F16::from_num(0.5);
//! let (s, c) = (sin(angle), cos(angle));
//!
//! let root = sqrt(I16F16::from_num(2.0)).unwrap();
//! assert!((root.to_num::<f32>() - 1.414).abs() < 0.001);
//!
//! let log = ln(I16F16::E).unwrap();
//! assert!((log.to_num::<f32>() - 1.0).abs() < 0.01);
//! ```
//!
//! # Available Functions
//!
//! **Total functions** return `T` directly, saturating on overflow.
//! **Fallible functions** return [`Result<T, Error>`] on domain violations.
//!
//! | Category | Total | Fallible |
//! |--------------|-------|----------|
//! | Trigonometric | [`sin`], [`cos`], [`tan`], [`sin_cos`], [`atan`], [`atan2`] | [`asin`], [`acos`] |
//! | Hyperbolic | [`sinh`], [`cosh`], [`tanh`], [`sinh_cosh`], [`asinh`] | [`acosh`], [`atanh`], [`acoth`], [`coth`] |
//! | Exponential | [`exp`], [`pow2`] | [`ln`], [`log2`], [`log10`] |
//! | Algebraic | — | [`sqrt`] |
//!
//! Functions use polynomial evaluation, CORDIC, and Newton-Raphson techniques.
//! Complete absence of panic is verified at the linker level via the
//! [`no-panic`](https://github.com/dtolnay/no-panic) crate.
//!
//! # Accuracy
//!
//! | Type | Typical Accuracy |
//! |------|------------------|
//! | `I16F16` | ~4 decimal digits |
//! | `I32F32` | ~8 decimal digits |
//!
//! All functions are benchmarked against MPFR reference implementations.
//! Accuracy regressions are not permitted across releases.
//!
//! # Features
//!
//! - **`std`** (default): Enables [`std::error::Error`] impl on [`Error`]
//!
//! See the [`kernel`] module for algorithm details.

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
