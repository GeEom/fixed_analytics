//! Precomputed lookup tables for CORDIC algorithms.
//!
//! Tables are stored as `i64` values representing signed I1F63 fixed-point
//! numbers (1 sign bit, 63 fractional bits), which are then converted to
//! the target type at runtime. Some constants (like `HYPERBOLIC_GAIN_INV`)
//! use I2F62 format for values >= 1.
//!
//! # Table Contents
//!
//! - [`ATAN_TABLE`]: `atan(2^-i)` values for circular CORDIC mode
//! - [`ATANH_TABLE`]: `atanh(2^-i)` values for hyperbolic CORDIC mode
//! - [`CIRCULAR_GAIN_INV`]: Inverse of circular gain factor (1/K ≈ 0.6073)
//! - [`HYPERBOLIC_GAIN`]: Hyperbolic gain factor (`K_h` ≈ 0.8282)
//! - [`HYPERBOLIC_GAIN_INV`]: Inverse hyperbolic gain factor (`1/K_h` ≈ 1.2075)

pub mod circular;
pub mod hyperbolic;

pub use circular::{ATAN_TABLE, CIRCULAR_GAIN_INV};
pub use hyperbolic::{ATANH_TABLE, HYPERBOLIC_GAIN, HYPERBOLIC_GAIN_INV};
