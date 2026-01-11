//! Core CORDIC iteration implementations.
//!
//! The CORDIC algorithm operates in two modes, each with two directions:
//!
//! | Mode | Rotation (z → 0) | Vectoring (y → 0) |
//! |------|------------------|-------------------|
//! | Circular | sin, cos | atan |
//! | Hyperbolic | sinh, cosh | atanh, ln |
//!
//! # Algorithm
//!
//! Each iteration performs a micro-rotation:
//! ```text
//! x' = x - σ * d * y * 2^(-i)
//! y' = y + σ * x * 2^(-i)
//! z' = z - σ * angle[i]
//! ```
//!
//! Where:
//! - σ = ±1 (direction of rotation)
//! - d = +1 for circular, -1 for hyperbolic, 0 for linear
//! - angle[i] = atan(2^-i) for circular, atanh(2^-i) for hyperbolic

use crate::tables::hyperbolic::needs_repeat;
use crate::tables::{
    ATAN_TABLE, ATANH_TABLE, CIRCULAR_GAIN_INV, HYPERBOLIC_GAIN, HYPERBOLIC_GAIN_INV,
};
use crate::traits::CordicNumber;

/// Table lookup for CORDIC iteration.
///
/// Index is bounded by CORDIC iteration limits:
/// - Circular mode: `min(frac_bits, 62)` → max index 61
/// - Hyperbolic mode: `min(frac_bits, 54)` with `i.saturating_sub(1)` → max index 53
///
/// Since the tables have 64 elements and max index is 61, bounds are always satisfied.
#[inline]
const fn table_lookup(table: &[i64; 64], index: u32) -> i64 {
    #[allow(
        clippy::indexing_slicing,
        reason = "index bounded by CORDIC iteration limits"
    )]
    table[index as usize]
}

/// Returns the CORDIC scale factor (1/K ≈ 0.6073).
///
/// Pre-multiply initial vectors by this to compensate for CORDIC gain.
#[inline]
#[must_use]
pub fn cordic_scale_factor<T: CordicNumber>() -> T {
    T::from_i1f63(CIRCULAR_GAIN_INV)
}

/// Returns the hyperbolic gain factor (`K_h` ≈ 0.8282).
///
/// After hyperbolic CORDIC iterations, results are scaled by `1/K_h`.
/// To compensate, divide by `K_h` (or multiply by `1/K_h`).
#[inline]
#[must_use]
pub fn hyperbolic_gain<T: CordicNumber>() -> T {
    T::from_i1f63(HYPERBOLIC_GAIN)
}

/// Returns the inverse hyperbolic gain factor (`1/K_h` ≈ 1.2075).
///
/// Pre-multiply initial vectors by this to compensate for hyperbolic CORDIC gain.
/// This uses a precomputed constant, avoiding runtime division.
#[inline]
#[must_use]
pub fn hyperbolic_gain_inv<T: CordicNumber>() -> T {
    T::from_i2f62(HYPERBOLIC_GAIN_INV)
}

/// Performs circular CORDIC in rotation mode.
///
/// Given an initial vector (x, y) and angle z, rotates the vector by angle z.
/// After iteration:
/// - x ≈ K * (x₀ * cos(z₀) - y₀ * sin(z₀))
/// - y ≈ K * (y₀ * cos(z₀) + x₀ * sin(z₀))
/// - z ≈ 0
///
/// Where K is the circular gain factor (~1.6468).
///
/// # Arguments
///
/// * `x` - Initial x coordinate
/// * `y` - Initial y coordinate
/// * `z` - Angle to rotate by (in radians)
///
/// # Returns
///
/// Tuple of (x, y, z) after CORDIC iterations.
///
/// # Note
///
/// The input angle should be in the range [-1.74, 1.74] radians for
/// convergence. Use argument reduction for larger angles.
#[must_use]
pub fn circular_rotation<T: CordicNumber>(mut x: T, mut y: T, mut z: T) -> (T, T, T) {
    let zero = T::zero();
    let iterations = T::frac_bits().min(62);

    for i in 0..iterations {
        let angle = T::from_i1f63(table_lookup(&ATAN_TABLE, i));

        if z >= zero {
            let x_new = x.saturating_sub(y >> i);
            y = y.saturating_add(x >> i);
            x = x_new;
            z -= angle;
        } else {
            let x_new = x.saturating_add(y >> i);
            y = y.saturating_sub(x >> i);
            x = x_new;
            z += angle;
        }
    }

    (x, y, z)
}

/// Performs circular CORDIC in vectoring mode.
///
/// Given an initial vector (x, y), rotates it until y ≈ 0.
/// After iteration:
/// - x ≈ K * sqrt(x₀² + y₀²)
/// - y ≈ 0
/// - z ≈ z₀ + atan(y₀/x₀)
///
/// # Arguments
///
/// * `x` - Initial x coordinate (should be positive for standard use)
/// * `y` - Initial y coordinate
/// * `z` - Initial angle accumulator (usually 0)
///
/// # Returns
///
/// Tuple of (x, y, z) after CORDIC iterations.
///
/// # Note
///
/// For computing atan(y/x), pass (1, y/x, 0) or (x, y, 0).
#[must_use]
pub fn circular_vectoring<T: CordicNumber>(mut x: T, mut y: T, mut z: T) -> (T, T, T) {
    let zero = T::zero();
    let iterations = T::frac_bits().min(62);

    for i in 0..iterations {
        let angle = T::from_i1f63(table_lookup(&ATAN_TABLE, i));

        if y < zero {
            // y is negative, rotate counter-clockwise to bring y toward zero
            let x_new = x.saturating_sub(y >> i);
            y = y.saturating_add(x >> i);
            x = x_new;
            z -= angle;
        } else {
            // y is positive or zero, rotate clockwise
            let x_new = x.saturating_add(y >> i);
            y = y.saturating_sub(x >> i);
            x = x_new;
            z += angle;
        }
    }

    (x, y, z)
}

/// Performs hyperbolic CORDIC in rotation mode.
///
/// Given initial values (x, y, z), performs hyperbolic pseudo-rotations
/// to drive z toward zero.
///
/// After iteration:
/// - x ≈ `K_h` * (x₀ * cosh(z₀) + y₀ * sinh(z₀))
/// - y ≈ `K_h` * (y₀ * cosh(z₀) + x₀ * sinh(z₀))
/// - z ≈ 0
///
/// Where `K_h` is the hyperbolic gain factor (~1.2075).
///
/// # Arguments
///
/// * `x` - Initial x value
/// * `y` - Initial y value
/// * `z` - Hyperbolic angle to "rotate" by
///
/// # Returns
///
/// Tuple of (x, y, z) after CORDIC iterations.
///
/// # Note
///
/// - Hyperbolic CORDIC starts at i=1 (not i=0)
/// - Certain iterations must be repeated for convergence
/// - Input z should be in range [-1.12, 1.12] for convergence
#[must_use]
pub fn hyperbolic_rotation<T: CordicNumber>(mut x: T, mut y: T, mut z: T) -> (T, T, T) {
    let zero = T::zero();
    // Use frac_bits iterations, capped at 54 for table bounds.
    let max_iterations = T::frac_bits().min(54);

    let mut i: u32 = 1; // Hyperbolic starts at i=1
    let mut iteration_count: u32 = 0;
    let mut repeated = false;

    while iteration_count < max_iterations && i < 64 {
        let table_index = i.saturating_sub(1);
        let angle = T::from_i1f63(table_lookup(&ATANH_TABLE, table_index));

        if z >= zero {
            // "Rotate" in positive direction
            let x_new = x.saturating_add(y >> i);
            y = y.saturating_add(x >> i);
            x = x_new;
            z -= angle;
        } else {
            // "Rotate" in negative direction
            let x_new = x.saturating_sub(y >> i);
            y = y.saturating_sub(x >> i);
            x = x_new;
            z += angle;
        }

        iteration_count += 1;

        // Handle repetition for convergence
        if needs_repeat(i) && !repeated {
            repeated = true;
            // Don't increment i, repeat this iteration
        } else {
            repeated = false;
            i += 1;
        }
    }

    (x, y, z)
}

/// Performs hyperbolic CORDIC in vectoring mode.
///
/// Drives y toward zero while accumulating the hyperbolic angle.
///
/// After iteration:
/// - x ≈ `K_h` * sqrt(x₀² - y₀²) (for |x| > |y|)
/// - y ≈ 0
/// - z ≈ z₀ + atanh(y₀/x₀)
///
/// # Arguments
///
/// * `x` - Initial x value (should satisfy |x| > |y|)
/// * `y` - Initial y value
/// * `z` - Initial angle accumulator
///
/// # Returns
///
/// Tuple of (x, y, z) after CORDIC iterations.
///
/// # Note
///
/// For computing atanh(v), pass (1, v, 0) where |v| < 1.
/// For computing ln(x), use the identity: ln(x) = 2 * atanh((x-1)/(x+1))
#[must_use]
pub fn hyperbolic_vectoring<T: CordicNumber>(mut x: T, mut y: T, mut z: T) -> (T, T, T) {
    let zero = T::zero();
    // Use at least 24 iterations for better accuracy, even for lower precision types.
    let max_iterations = T::frac_bits().clamp(24, 54);

    let mut i: u32 = 1;
    let mut iteration_count: u32 = 0;
    let mut repeated = false;

    while iteration_count < max_iterations && i < 64 {
        let table_index = i.saturating_sub(1);
        let angle = T::from_i1f63(table_lookup(&ATANH_TABLE, table_index));

        // Hyperbolic pseudo-rotation equations:
        // x' = x + σ*y*2^(-i)
        // y' = y + σ*x*2^(-i)
        // z' = z + σ*angle  (accumulating for vectoring)
        // where σ = -sign(y) to drive y toward zero

        if y < zero {
            // y is negative: σ = +1
            // x' = x + y*2^(-i) [y is negative, so this subtracts magnitude]
            // y' = y + x*2^(-i) [adds positive to make less negative]
            // z' = z - angle    [accumulate negative contribution]
            let x_new = x.saturating_add(y >> i);
            y = y.saturating_add(x >> i);
            x = x_new;
            z -= angle;
        } else {
            // y is positive or zero: σ = -1
            // x' = x - y*2^(-i) [subtracts positive]
            // y' = y - x*2^(-i) [subtracts to decrease toward zero]
            // z' = z + angle    [accumulate positive contribution]
            let x_new = x.saturating_sub(y >> i);
            y = y.saturating_sub(x >> i);
            x = x_new;
            z += angle;
        }

        iteration_count += 1;

        if needs_repeat(i) && !repeated {
            repeated = true;
        } else {
            repeated = false;
            i += 1;
        }
    }

    (x, y, z)
}
