//! # `fixed_analytics`
//!
//! Fixed-point mathematical functions using the CORDIC algorithm.
//!
//! This crate provides efficient implementations of trigonometric, hyperbolic,
//! exponential, and algebraic functions for fixed-point numbers. All algorithms
//! use only addition, subtraction, and bit shifts, making them suitable for
//! embedded systems without hardware floating-point support.
//!
//! ## Features
//!
//! - **No floating-point operations**: All computations use fixed-point arithmetic
//! - **`no_std` compatible**: Works on embedded systems without an allocator
//! - **Comprehensive function coverage**: Trig, hyperbolic, exponential, and more
//! - **Compile-time tables**: Lookup tables are embedded in the binary
//! - **Proper error handling**: Domain errors return `Result` types
//!
//! ## Supported Functions
//!
//! | Category | Functions |
//! |----------|-----------|
//! | Circular | [`sin`], [`cos`], [`tan`], [`sin_cos`], [`asin`], [`acos`], [`atan`], [`atan2`] |
//! | Hyperbolic | [`sinh`], [`cosh`], [`tanh`], [`coth`], [`sinh_cosh`], [`asinh`], [`acosh`], [`atanh`], [`acoth`] |
//! | Exponential | [`exp`], [`ln`], [`log2`], [`log10`] |
//! | Algebraic | [`sqrt`] |
//!
//! ## Quick Start
//!
//! ```rust
//! use fixed::types::I16F16;
//! use fixed_analytics::{sin, cos, sqrt, ln};
//!
//! // Compute sin and cos of an angle
//! let angle = I16F16::from_num(0.5); // 0.5 radians
//! let sine = sin(angle);
//! let cosine = cos(angle);
//!
//! // Square root
//! let x = I16F16::from_num(2.0);
//! let root = sqrt(x); // ≈ 1.414
//!
//! // Natural logarithm
//! let y = I16F16::from_num(2.718);
//! let log = ln(y).unwrap(); // ≈ 1.0
//! ```
//!
//! ## Supported Fixed-Point Types
//!
//! This crate works with signed fixed-point types from the [`fixed`] crate:
//!
//! - [`FixedI8<Fract>`](fixed::FixedI8) - 8-bit (limited precision)
//! - [`FixedI16<Fract>`](fixed::FixedI16) - 16-bit
//! - [`FixedI32<Fract>`](fixed::FixedI32) - 32-bit (recommended)
//! - [`FixedI64<Fract>`](fixed::FixedI64) - 64-bit (highest precision)
//!
//! Common type aliases like [`I16F16`](fixed::types::I16F16) and
//! [`I32F32`](fixed::types::I32F32) work well with this library.
//!
//! ## Algorithm Overview
//!
//! CORDIC (Coordinate Rotation Digital Computer) is an iterative algorithm
//! invented by Jack Volder in 1959. It computes trigonometric, hyperbolic,
//! and other functions using only:
//!
//! - Addition and subtraction
//! - Bit shifts (multiplication/division by powers of 2)
//! - Table lookups
//!
//! This makes it ideal for hardware without a hardware multiplier or FPU.
//!
//! ### How CORDIC Works
//!
//! The algorithm rotates a vector through a series of predetermined angles.
//! Each rotation uses the identities:
//!
//! ```text
//! x' = x - σ × y × 2^(-i)
//! y' = y + σ × x × 2^(-i)
//! z' = z - σ × angle[i]
//! ```
//!
//! Where σ = ±1 determines the rotation direction. After n iterations:
//!
//! - **Rotation mode** (z → 0): Computes sin and cos of the initial angle
//! - **Vectoring mode** (y → 0): Computes the angle of the initial vector
//!
//! ## Precision
//!
//! Accuracy depends on the fixed-point type used:
//!
//! | Type | Typical Accuracy |
//! |------|------------------|
//! | `I8F8` | ~2 decimal digits |
//! | `I16F16` | ~4 decimal digits |
//! | `I32F32` | ~8 decimal digits |
//! | `I64F64` | ~15 decimal digits |
//!
//! ## Feature Flags
//!
//! - `std` (default): Enables `std::error::Error` implementation
//! - Without `std`: `#![no_std]` compatible
//!
//! ## References
//!
//! - [CORDIC on Wikipedia](https://en.wikipedia.org/wiki/CORDIC)
//! - Volder, J.E. "The CORDIC Trigonometric Computing Technique" (1959)
//! - Walther, J.S. "A Unified Algorithm for Elementary Functions" (1971)

#![no_std]
#![cfg_attr(docsrs, feature(doc_cfg))]
#![warn(missing_docs)]
#![warn(clippy::pedantic)]
#![allow(clippy::module_name_repetitions)]
#![allow(clippy::similar_names)]

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
