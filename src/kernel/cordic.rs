//! CORDIC vectoring kernels; the [kernel module docs](super) cover the algorithm.

use crate::tables::hyperbolic::needs_repeat;
use crate::tables::{ATAN_TABLE, ATANH_TABLE};
use crate::traits::CordicNumber;

/// Last shift stage before the closing division, at most 42. After stage `k`
/// the residual `|y/x|` is below `2^-k`.
const fn last_stage(frac_bits: u32) -> u32 {
    frac_bits.div_ceil(3)
}

/// Table lookup; the modulo keeps the index in bounds, so the fallback is dead.
#[inline]
fn table_lookup(table: &[i64; 64], index: u32) -> i64 {
    match table.get(index as usize % table.len()) {
        Some(&v) => v,
        None => 0,
    }
}

/// Converts an I1F63 table angle to `T`, rounding to nearest: truncation would
/// bias circular z low. Hyperbolic mode keeps truncation, which partly cancels
/// its shift bias (measured in the accuracy baseline).
fn from_i1f63_rounded<T: CordicNumber>(bits: i64) -> T {
    let truncated = T::from_i1f63(bits);
    let frac = T::frac_bits();
    if frac >= 63 {
        return truncated;
    }
    // Round up when the highest discarded bit is set.
    let shift = 63 - frac;
    if (bits >> (shift - 1)) & 1 == 1 {
        let ulp = T::one() >> frac;
        truncated.saturating_add(ulp)
    } else {
        truncated
    }
}

/// Circular CORDIC in vectoring mode for `x > 0`; `(x, y, 0)` gives atan(y/x).
///
/// Returns z ≈ z₀ + atan(y₀/x₀), x ≈ K·√(x₀² + y₀²) with K the stages' gain,
/// and the residual y, whose angle is already in z.
#[must_use]
pub fn circular_vectoring<T: CordicNumber>(mut x: T, mut y: T, mut z: T) -> (T, T, T) {
    let zero = T::zero();

    for i in 0..=last_stage(T::frac_bits()) {
        let angle = from_i1f63_rounded::<T>(table_lookup(&ATAN_TABLE, i));

        if y < zero {
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

    // Close the residual angle: atan(y/x) = y/x to within |y/x|³/3.
    z = z.saturating_add(y.div(x));

    (x, y, z)
}

/// Hyperbolic CORDIC in vectoring mode for `|x| > |y|`; `(1, v, 0)` gives atanh(v).
///
/// Returns z ≈ z₀ + atanh(y₀/x₀), x ≈ `K_h`·√(x₀² − y₀²) with `K_h` the stages'
/// gain, and the residual y, whose angle is already in z.
#[must_use]
pub fn hyperbolic_vectoring<T: CordicNumber>(mut x: T, mut y: T, mut z: T) -> (T, T, T) {
    let zero = T::zero();

    for i in 1..=last_stage(T::frac_bits()) {
        let angle = T::from_i1f63(table_lookup(&ATANH_TABLE, i - 1));
        // Stages 4, 13, 40, … run twice so the sequence converges.
        let passes = if needs_repeat(i) { 2 } else { 1 };

        for _ in 0..passes {
            // σ = −sign(y) drives y toward zero.
            if y < zero {
                let x_new = x.saturating_add(y >> i);
                y = y.saturating_add(x >> i);
                x = x_new;
                z -= angle;
            } else {
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
