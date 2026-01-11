//! Trigonometric functions via circular CORDIC.

use crate::bounded::{NonNegative, UnitInterval};
use crate::error::{Error, Result};
use crate::kernel::{circular_rotation, circular_vectoring, cordic_scale_factor};
use crate::ops::algebraic::sqrt_nonneg;
use crate::traits::CordicNumber;

/// Sine and cosine. More efficient than separate calls. Accepts any angle.
#[must_use]
pub fn sin_cos<T: CordicNumber>(angle: T) -> (T, T) {
    let pi = T::pi();
    let frac_pi_2 = T::frac_pi_2();
    let zero = T::zero();

    // Reduce angle to [-π, π] range first.
    let mut reduced = angle;
    let two_pi = pi + pi;
    let mut i = 0;
    while reduced > pi && i < 64 {
        reduced -= two_pi;
        i += 1;
    }
    i = 0;
    while reduced < -pi && i < 64 {
        reduced += two_pi;
        i += 1;
    }

    // Further reduce to [-π/2, π/2] and track sign
    let (reduced, negate) = if reduced > frac_pi_2 {
        (reduced - pi, true)
    } else if reduced < -frac_pi_2 {
        (reduced + pi, true)
    } else {
        (reduced, false)
    };

    // Run CORDIC with unit vector scaled by inverse gain
    let inv_gain = cordic_scale_factor();
    let (cos_val, sin_val, _) = circular_rotation(inv_gain, zero, reduced);

    if negate {
        (-sin_val, -cos_val)
    } else {
        (sin_val, cos_val)
    }
}

/// Sine. Accepts any angle (reduced internally).
#[inline]
#[must_use]
pub fn sin<T: CordicNumber>(angle: T) -> T {
    sin_cos(angle).0
}

/// Cosine. Accepts any angle (reduced internally).
#[inline]
#[must_use]
pub fn cos<T: CordicNumber>(angle: T) -> T {
    sin_cos(angle).1
}

/// Tangent. May overflow near ±π/2.
#[must_use]
pub fn tan<T: CordicNumber>(angle: T) -> T {
    let (s, c) = sin_cos(angle);
    s.div(c)
}

/// Arcsine. Domain: `[-1, 1]`. Returns angle in `[-π/2, π/2]`.
///
/// # Errors
/// Returns `DomainError` if `|x| > 1`.
#[must_use = "returns the arcsine result which should be handled"]
pub fn asin<T: CordicNumber>(x: T) -> Result<T> {
    let Some(unit_x) = UnitInterval::new(x) else {
        return Err(Error::domain("asin", "value in range [-1, 1]"));
    };

    // Special cases
    if x == T::one() {
        return Ok(T::frac_pi_2());
    }
    if x == -T::one() {
        return Ok(-T::frac_pi_2());
    }
    if x == T::zero() {
        return Ok(T::zero());
    }

    // Use the identity: asin(x) = atan(x / sqrt(1 - x²))
    // NonNegative::one_minus_square gives 1 - x², which is ≥ 0 since |x| ≤ 1
    let sqrt_term = sqrt_nonneg(NonNegative::one_minus_square(unit_x));

    // Handle case where sqrt_term is very small (x close to ±1)
    if sqrt_term < T::from_i1f63(0x0001_0000_0000_0000) {
        // Very close to ±1, return ±π/2
        return if x.is_positive() {
            Ok(T::frac_pi_2())
        } else {
            Ok(-T::frac_pi_2())
        };
    }

    Ok(atan(x.div(sqrt_term)))
}

/// Arccosine. Domain: `[-1, 1]`. Returns angle in `[0, π]`.
///
/// # Errors
/// Returns `DomainError` if `|x| > 1`.
#[must_use = "returns the arccosine result which should be handled"]
pub fn acos<T: CordicNumber>(x: T) -> Result<T> {
    // acos(x) = π/2 - asin(x)
    asin(x).map(|a| T::frac_pi_2() - a)
}

/// Arctangent. Accepts any value. Returns angle in `(-π/2, π/2)`.
#[must_use]
pub fn atan<T: CordicNumber>(x: T) -> T {
    let zero = T::zero();
    let one = T::one();

    // Special cases
    if x == zero {
        return zero;
    }

    // For |x| > 1, use atan(x) = sign(x) * π/2 - atan(1/x)
    // This keeps the argument in the convergent range
    let abs_x = x.abs();
    if abs_x > one {
        let recip = one.div(x);
        let atan_recip = circular_vectoring(one, recip, zero).2;

        if x.is_positive() {
            T::frac_pi_2() - atan_recip
        } else {
            -T::frac_pi_2() - atan_recip
        }
    } else {
        // |x| <= 1, use CORDIC directly
        circular_vectoring(one, x, zero).2
    }
}

/// Four-quadrant arctangent. Returns angle in `[-π, π]`. Returns 0 for (0, 0).
#[must_use]
pub fn atan2<T: CordicNumber>(y: T, x: T) -> T {
    let zero = T::zero();
    let pi = T::pi();
    let frac_pi_2 = T::frac_pi_2();

    // Handle special cases
    if x == zero {
        return if y.is_negative() {
            -frac_pi_2
        } else if y == zero {
            zero // Undefined, but return 0
        } else {
            frac_pi_2
        };
    }

    if y == zero {
        return if x.is_negative() { pi } else { zero };
    }

    // Compute atan(|y|/|x|) using CORDIC vectoring mode
    // Using absolute values ensures the base angle is always positive
    let (_, _, base_angle) = circular_vectoring(x.abs(), y.abs(), zero);

    // Adjust for quadrant based on signs of original x and y
    match (x.is_negative(), y.is_negative()) {
        // Q1: x positive, y positive -> angle is base_angle
        (false, false) => base_angle,
        // Q4: x positive, y negative -> angle is -base_angle
        (false, true) => -base_angle,
        // Q2: x negative, y positive -> angle is π - base_angle
        (true, false) => pi - base_angle,
        // Q3: x negative, y negative -> angle is -(π - base_angle) = base_angle - π
        (true, true) => base_angle - pi,
    }
}
