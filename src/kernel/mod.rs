//! Core CORDIC iteration kernels.
//!
//! This module contains the fundamental CORDIC algorithm implementations
//! for both circular and hyperbolic modes, in both rotation and vectoring
//! directions.
//!
//! # CORDIC Overview
//!
//! CORDIC (Coordinate Rotation Digital Computer) is an iterative algorithm
//! that computes trigonometric, hyperbolic, and other functions using only
//! addition, subtraction, and bit shifts.
//!
//! ## Modes and Directions
//!
//! | Mode | Rotation (z → 0) | Vectoring (y → 0) |
//! |------|------------------|-------------------|
//! | Circular | sin, cos from angle | atan from (x, y) |
//! | Hyperbolic | sinh, cosh from arg | atanh, ln |
//!
//! ## Usage
//!
//! The kernels are building blocks for higher-level functions. Users should
//! generally use the functions in [`crate::ops`] instead of calling kernels
//! directly.

mod cordic;

pub use crate::kernel::cordic::{circular_gain_inv, circular_rotation, circular_vectoring};
pub use crate::kernel::cordic::{
    hyperbolic_gain, hyperbolic_gain_inv, hyperbolic_rotation, hyperbolic_vectoring,
};
