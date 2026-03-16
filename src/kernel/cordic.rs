//! Core CORDIC iteration implementations.
//!
//! CORDIC vectoring mode drives y toward zero while accumulating angles:
//!
//! | Mode | Vectoring (y → 0) |
//! |------|-------------------|
//! | Circular | atan |
//! | Hyperbolic | atanh, ln |
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
use crate::tables::{ATAN_TABLE, ATANH_TABLE};
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
