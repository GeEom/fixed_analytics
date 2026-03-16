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
//! **Vectoring mode** (y→0): computes atan/atanh from coordinates.
//!
//! Hyperbolic mode uses `atanh(2^-i)` tables and requires iteration repeats
//! at indices 4, 13, 40, ... for convergence.
//!
//! | Mode | Vectoring (y → 0) |
//! |------|-------------------|
//! | Circular | atan |
//! | Hyperbolic | atanh, ln |
//!
//! Users should call functions in [`crate::ops`] rather than kernels directly.

mod cordic;

pub use crate::kernel::cordic::{circular_vectoring, hyperbolic_vectoring};
