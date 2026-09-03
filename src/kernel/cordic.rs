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
//!
//! Both kernels stop after about a third of the fractional bits and close
//! the residual angle with one division: for `|w| = |y/x| < 2^-⌈frac/3⌉`,
//! `atan(w)` and `atanh(w)` differ from `w` by under `|w|³/3 < 2^-frac/3`.

use crate::tables::hyperbolic::needs_repeat;
use crate::tables::{ATAN_TABLE, ATANH_TABLE};
use crate::traits::CordicNumber;

/// Index of the last shift stage before the closing division: after the
/// stage with shift `2^-k` the residual angle, and so `|y/x|`, is below
/// `2^-k`. At most 42, so every table lookup is in bounds.
const fn last_stage(frac_bits: u32) -> u32 {
    frac_bits.div_ceil(3)
}

/// Table lookup for CORDIC iteration. The modulo is a no-op that keeps the
/// access provably in bounds for the no-panic check.
#[inline]
const fn table_lookup(table: &[i64; 64], index: u32) -> i64 {
    #[allow(clippy::indexing_slicing, reason = "index reduced modulo table length")]
    table[index as usize % table.len()]
}

/// Converts an I1F63 table constant to `T`, rounding to nearest.
///
/// `CordicNumber::from_i1f63` truncates, which biases every converted
/// angle low by up to one ulp. For circular CORDIC the accumulated z is
/// dominated by this table quantization, so rounding to nearest halves
/// the worst case and removes the systematic bias.
///
/// Hyperbolic CORDIC deliberately keeps the truncating conversion: its
/// datapath shifts carry their own bias, which the truncated table
/// angles partially cancel (measured, not designed — see the accuracy
/// baseline history).
fn from_i1f63_rounded<T: CordicNumber>(bits: i64) -> T {
    let truncated = T::from_i1f63(bits);
    let frac = T::frac_bits();
    if frac >= 63 {
        // Conversion is exact; nothing to round.
        return truncated;
    }
    // The conversion discards the low (63 - frac) bits. If the highest
    // discarded bit is set, the true value is at least half an output
    // ulp above the truncated result, so round up by one ulp.
    let shift = 63 - frac;
    if (bits >> (shift - 1)) & 1 == 1 {
        let ulp = T::one() >> frac;
        truncated.saturating_add(ulp)
    } else {
        truncated
    }
}

/// Performs circular CORDIC in vectoring mode.
///
/// Given an initial vector (x, y) with x > 0, rotates it toward the x axis
/// through `⌈frac/3⌉ + 1` shift stages, then closes the residual angle
/// with one division.
///
/// # Arguments
///
/// * `x` - Initial x coordinate (must be positive)
/// * `y` - Initial y coordinate
/// * `z` - Initial angle accumulator (usually 0)
///
/// # Returns
///
/// Tuple of (x, y, z): z ≈ z₀ + atan(y₀/x₀); x ≈ K·sqrt(x₀² + y₀²) with K
/// the gain of the stages run; y is the residual the shift stages left,
/// whose angle is already in z.
///
/// # Note
///
/// For computing atan(y/x), pass (1, y/x, 0) or (x, y, 0).
#[must_use]
pub fn circular_vectoring<T: CordicNumber>(mut x: T, mut y: T, mut z: T) -> (T, T, T) {
    let zero = T::zero();

    for i in 0..=last_stage(T::frac_bits()) {
        let angle = from_i1f63_rounded::<T>(table_lookup(&ATAN_TABLE, i));

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

    // Close the residual angle: atan(y/x) = y/x to within |y/x|³/3.
    z = z.saturating_add(y.div(x));

    (x, y, z)
}

/// Performs hyperbolic CORDIC in vectoring mode.
///
/// Drives y toward zero through shift stages 1..=⌈frac/3⌉ (repeating 4,
/// 13, 40, … for convergence), then closes the residual angle with one
/// division.
///
/// # Arguments
///
/// * `x` - Initial x value (should satisfy |x| > |y|)
/// * `y` - Initial y value
/// * `z` - Initial angle accumulator
///
/// # Returns
///
/// Tuple of (x, y, z): z ≈ z₀ + atanh(y₀/x₀); x ≈ `K_h`·sqrt(x₀² - y₀²)
/// with `K_h` the gain of the stages run; y is the residual the shift
/// stages left, whose angle is already in z.
///
/// # Note
///
/// For computing atanh(v), pass (1, v, 0) where |v| < 1.
/// For computing ln(x), use the identity: ln(x) = 2 * atanh((x-1)/(x+1))
#[must_use]
pub fn hyperbolic_vectoring<T: CordicNumber>(mut x: T, mut y: T, mut z: T) -> (T, T, T) {
    let zero = T::zero();

    for i in 1..=last_stage(T::frac_bits()) {
        let angle = T::from_i1f63(table_lookup(&ATANH_TABLE, i - 1));
        // Stages 4, 13, 40, … run twice so the sequence converges.
        let passes = if needs_repeat(i) { 2 } else { 1 };

        for _ in 0..passes {
            // Hyperbolic pseudo-rotation equations:
            // x' = x + σ*y*2^(-i)
            // y' = y + σ*x*2^(-i)
            // z' = z + σ*angle  (accumulating for vectoring)
            // where σ = -sign(y) to drive y toward zero
            if y < zero {
                // y is negative: σ = +1
                let x_new = x.saturating_add(y >> i);
                y = y.saturating_add(x >> i);
                x = x_new;
                z -= angle;
            } else {
                // y is positive or zero: σ = -1
                let x_new = x.saturating_sub(y >> i);
                y = y.saturating_sub(x >> i);
                x = x_new;
                z += angle;
            }
        }
    }

    // Close the residual angle: atanh(y/x) = y/x to within |y/x|³/3.
    z = z.saturating_add(y.div(x));

    (x, y, z)
}
