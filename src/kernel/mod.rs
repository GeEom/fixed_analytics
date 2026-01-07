//! CORDIC (Coordinate Rotation Digital Computer) kernels.
//!
//! # Algorithm
//!
//! Iteratively rotates vectors using only shifts and adds:
//!
//! ```text
//! x' = x - σ·y·2^(-i)
//! y' = y + σ·x·2^(-i)
//! z' = z - σ·atan(2^-i)
//! ```
//!
//! **Rotation mode** (z→0): computes sin/cos from angle.
//! **Vectoring mode** (y→0): computes atan from coordinates.
//!
//! Hyperbolic mode uses `atanh(2^-i)` tables and requires iteration repeats
//! at indices 4, 13, 40, ... for convergence.
//!
//! | Mode | Rotation (z → 0) | Vectoring (y → 0) |
//! |------|------------------|-------------------|
//! | Circular | sin, cos | atan |
//! | Hyperbolic | sinh, cosh | atanh, ln |
//!
//! Users should call functions in [`crate::ops`] rather than kernels directly.

mod cordic;

pub use crate::kernel::cordic::{circular_rotation, circular_vectoring, cordic_scale_factor};
pub use crate::kernel::cordic::{
    hyperbolic_gain, hyperbolic_gain_inv, hyperbolic_rotation, hyperbolic_vectoring,
};
