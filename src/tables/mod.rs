//! Precomputed lookup tables for CORDIC algorithms and polynomial evaluation.
//!
//! Tables are stored as `i64` values representing signed I1F63 fixed-point
//! numbers (1 sign bit, 63 fractional bits), which are then converted to
//! the target type at runtime.
//!
//! # Table Contents
//!
//! - [`ATAN_TABLE`]: `atan(2^-i)` values for circular CORDIC mode
//! - [`ATANH_TABLE`]: `atanh(2^-i)` values for hyperbolic CORDIC mode
//! - [`chebyshev`]: Minimax polynomial coefficients for sin/cos evaluation

pub mod chebyshev;
pub mod circular;
pub mod hyperbolic;

pub use circular::ATAN_TABLE;
pub use hyperbolic::ATANH_TABLE;
