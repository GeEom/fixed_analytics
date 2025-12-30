//! Circular (trigonometric) functions.
//!
//! Provides sin, cos, tan, asin, acos, atan, and atan2 using circular CORDIC.
//!
//! # Precision
//!
//! These functions achieve accuracy comparable to the precision of the
//! fixed-point type used. For `I16F16`, expect ~4 decimal digits of accuracy.
//!
//! # Range
//!
//! - `sin`, `cos`, `tan`: Accept any angle, with automatic reduction to [-π, π]
//! - `asin`, `acos`: Domain is [-1, 1]
//! - `atan`: Accepts any value
//! - `atan2`: Accepts any (y, x) pair except (0, 0)

use crate::error::{Error, Result};
use crate::kernel::{circular_gain_inv, circular_rotation, circular_vectoring};
use crate::traits::CordicNumber;

/// Computes the sine and cosine of an angle simultaneously.
///
/// This is more efficient than calling [`sin`] and [`cos`] separately since
/// both values are produced by a single CORDIC iteration.
///
/// # Arguments
///
/// * `angle` - The angle in radians
///
/// # Returns
///
/// A tuple `(sin(angle), cos(angle))`.
///
/// # Examples
///
/// ```
/// use fixed::types::I16F16;
/// use fixed_analytics::sin_cos;
///
/// let angle = I16F16::from_num(1.0); // 1 radian
/// let (s, c) = sin_cos(angle);
/// // s ≈ 0.841, c ≈ 0.540
/// ```
#[must_use]
pub fn sin_cos<T: CordicNumber>(angle: T) -> (T, T) {
    let pi = T::pi();
    let frac_pi_2 = T::frac_pi_2();
    let zero = T::zero();

    // Reduce angle to [-π, π] range first.
    // Guard against infinite loops for extreme values (limit iterations to prevent hangs).
    let mut reduced = angle;
    let two_pi = pi + pi;
    let max_iterations = 64_u32; // More than enough for any representable angle
    let mut iterations = 0_u32;
    while reduced > pi && iterations < max_iterations {
        reduced -= two_pi;
        iterations += 1;
    }
    while reduced < -pi && iterations < max_iterations {
        reduced += two_pi;
        iterations += 1;
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
    let inv_gain = circular_gain_inv();
    let (cos_val, sin_val, _) = circular_rotation(inv_gain, zero, reduced);

    if negate {
        (-sin_val, -cos_val)
    } else {
        (sin_val, cos_val)
    }
}

/// Computes the sine of an angle.
///
/// # Arguments
///
/// * `angle` - The angle in radians
///
/// # Returns
///
/// The sine of the angle, in the range [-1, 1].
///
/// # Examples
///
/// ```
/// use fixed::types::I16F16;
/// use fixed_analytics::sin;
///
/// let angle = I16F16::FRAC_PI_2; // π/2
/// let result = sin(angle);
/// // result ≈ 1.0
/// ```
#[inline]
#[must_use]
pub fn sin<T: CordicNumber>(angle: T) -> T {
    sin_cos(angle).0
}

/// Computes the cosine of an angle.
///
/// # Arguments
///
/// * `angle` - The angle in radians
///
/// # Returns
///
/// The cosine of the angle, in the range [-1, 1].
///
/// # Examples
///
/// ```
/// use fixed::types::I16F16;
/// use fixed_analytics::cos;
///
/// let angle = I16F16::ZERO;
/// let result = cos(angle);
/// // result ≈ 1.0
/// ```
#[inline]
#[must_use]
pub fn cos<T: CordicNumber>(angle: T) -> T {
    sin_cos(angle).1
}

/// Computes the tangent of an angle.
///
/// # Arguments
///
/// * `angle` - The angle in radians
///
/// # Returns
///
/// The tangent of the angle.
///
/// # Note
///
/// May produce very large values or overflow for angles near ±π/2
/// where tangent approaches infinity.
///
/// # Examples
///
/// ```
/// use fixed::types::I16F16;
/// use fixed_analytics::tan;
///
/// let angle = I16F16::FRAC_PI_4; // π/4
/// let result = tan(angle);
/// // result ≈ 1.0
/// ```
#[must_use]
pub fn tan<T: CordicNumber>(angle: T) -> T {
    let (s, c) = sin_cos(angle);
    s.div(c)
}

/// Computes the arcsine (inverse sine) of a value.
///
/// # Arguments
///
/// * `x` - A value in the range [-1, 1]
///
/// # Returns
///
/// The angle whose sine is `x`, in the range [-π/2, π/2].
///
/// # Errors
///
/// Returns [`Error::DomainError`] if `|x| > 1`.
///
/// # Examples
///
/// ```
/// use fixed::types::I16F16;
/// use fixed_analytics::asin;
///
/// let x = I16F16::from_num(0.5);
/// let angle = asin(x).unwrap();
/// // angle ≈ π/6 ≈ 0.524
/// ```
#[must_use = "returns the arcsine result which should be handled"]
pub fn asin<T: CordicNumber>(x: T) -> Result<T> {
    let one = T::one();
    let neg_one = -one;

    if x > one || x < neg_one {
        return Err(Error::DomainError {
            function: "asin",
            expected: "value in range [-1, 1]",
        });
    }

    // Special cases
    if x == one {
        return Ok(T::frac_pi_2());
    }
    if x == neg_one {
        return Ok(-T::frac_pi_2());
    }
    if x == T::zero() {
        return Ok(T::zero());
    }

    // Use the identity: asin(x) = atan(x / sqrt(1 - x²))
    // This gives better accuracy than iterative methods
    let x_sq = x.saturating_mul(x);
    let one_minus_x_sq = one.saturating_sub(x_sq);
    let sqrt_term = crate::ops::algebraic::sqrt(one_minus_x_sq);

    // Handle case where sqrt_term is very small (x close to ±1)
    if sqrt_term < T::from_i64_frac(0x0001_0000_0000_0000) {
        // Very close to ±1, return ±π/2
        return if x.is_positive() {
            Ok(T::frac_pi_2())
        } else {
            Ok(-T::frac_pi_2())
        };
    }

    Ok(atan(x.div(sqrt_term)))
}

/// Computes the arccosine (inverse cosine) of a value.
///
/// # Arguments
///
/// * `x` - A value in the range [-1, 1]
///
/// # Returns
///
/// The angle whose cosine is `x`, in the range [0, π].
///
/// # Errors
///
/// Returns [`Error::DomainError`] if `|x| > 1`.
///
/// # Examples
///
/// ```
/// use fixed::types::I16F16;
/// use fixed_analytics::acos;
///
/// let x = I16F16::from_num(0.5);
/// let angle = acos(x).unwrap();
/// // angle ≈ π/3 ≈ 1.047
/// ```
#[must_use = "returns the arccosine result which should be handled"]
pub fn acos<T: CordicNumber>(x: T) -> Result<T> {
    // acos(x) = π/2 - asin(x)
    asin(x).map(|a| T::frac_pi_2() - a)
}

/// Computes the arctangent (inverse tangent) of a value.
///
/// # Arguments
///
/// * `x` - Any real value
///
/// # Returns
///
/// The angle whose tangent is `x`, in the range (-π/2, π/2).
///
/// # Examples
///
/// ```
/// use fixed::types::I16F16;
/// use fixed_analytics::atan;
///
/// let x = I16F16::from_num(1.0);
/// let angle = atan(x);
/// // angle ≈ π/4 ≈ 0.785
/// ```
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

/// Computes the four-quadrant arctangent of `y/x`.
///
/// This function returns the angle θ in radians between the positive x-axis
/// and the point (x, y), taking into account the signs of both arguments
/// to determine the correct quadrant.
///
/// # Arguments
///
/// * `y` - The y coordinate
/// * `x` - The x coordinate
///
/// # Returns
///
/// An angle in the range [-π, π].
///
/// # Special Cases
///
/// - `atan2(0, 0)` returns `0` (though mathematically undefined)
/// - `atan2(y, 0)` returns `±π/2` depending on the sign of `y`
/// - `atan2(0, x)` returns `0` or `π` depending on the sign of `x`
///
/// # Examples
///
/// ```
/// use fixed::types::I16F16;
/// use fixed_analytics::atan2;
///
/// let y = I16F16::from_num(1.0);
/// let x = I16F16::from_num(1.0);
/// let angle = atan2(y, x);
/// // angle ≈ π/4 ≈ 0.785 (first quadrant)
///
/// let y = I16F16::from_num(1.0);
/// let x = I16F16::from_num(-1.0);
/// let angle = atan2(y, x);
/// // angle ≈ 3π/4 ≈ 2.356 (second quadrant)
/// ```
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
