//! High-level mathematical operations built on CORDIC kernels.
//!
//! This module provides user-friendly functions for computing various
//! mathematical operations using the CORDIC algorithm.
//!
//! # Modules
//!
//! - [`circular`]: Trigonometric functions (sin, cos, tan, asin, acos, atan, atan2)
//! - [`hyperbolic`]: Hyperbolic functions (sinh, cosh, tanh, asinh, acosh, atanh, acoth)
//! - [`exponential`]: Exponential and logarithmic functions (exp, ln, log2, log10, pow2)
//! - [`algebraic`]: Algebraic functions (sqrt)

pub mod algebraic;
pub mod circular;
pub mod exponential;
pub mod hyperbolic;

// Re-export all public functions
pub use algebraic::sqrt;
pub use circular::{acos, asin, atan, atan2, cos, sin, sin_cos, tan};
pub use exponential::{exp, ln, log2, log10, pow2};
pub use hyperbolic::{acosh, acoth, asinh, atanh, cosh, coth, sinh, sinh_cosh, tanh};
