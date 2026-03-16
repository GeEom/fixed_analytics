//! Trigonometric functions via circular CORDIC.

use crate::bounded::{NonNegative, UnitInterval};
use crate::error::{Error, Result};
use crate::kernel::circular_vectoring;
use crate::ops::algebraic::sqrt_nonneg;
use crate::tables::chebyshev::{COS_Q_HI, COS_Q_LO, SIN_P_HI, SIN_P_LO, horner};
use crate::traits::CordicNumber;

/// Sine and cosine. More efficient than separate calls. Accepts any angle.
#[must_use]
#[cfg_attr(feature = "verify-no-panic", no_panic::no_panic)]
pub fn sin_cos<T: CordicNumber>(angle: T) -> (T, T) {
    let pi = T::pi();
    let frac_pi_2 = T::frac_pi_2();
    let two_pi = pi + pi;

    // Reduce angle to [-π, π] using direct quotient computation.
    // This handles arbitrarily large angles without iteration limits.
    let reduced = if angle > pi || angle < -pi {
        // Compute n = round(angle / 2π), then reduced = angle - n * 2π
        let quotient = angle.div(two_pi);
        let n = quotient.round();
        angle.saturating_sub(n.saturating_mul(two_pi))
    } else {
        angle
    };

    // Clamp to [-π, π] to handle any residual from saturation.
    // This is a safety net; mathematically unnecessary for valid inputs.
    let reduced = if reduced > pi {
        reduced.saturating_sub(two_pi)
    } else if reduced < -pi {
        reduced.saturating_add(two_pi)
    } else {
        reduced
    };

    // Further reduce to [-π/2, π/2] and track sign
    let (reduced, negate) = if reduced > frac_pi_2 {
        (reduced - pi, true)
    } else if reduced < -frac_pi_2 {
        (reduced + pi, true)
    } else {
        (reduced, false)
    };

    // Polynomial evaluation via factored Horner form.
    // To avoid catastrophic cancellation near π/2, reduce to [0, π/4]:
    //   For |x| ∈ [0, π/4]:      sin(x) = sin_poly(x), cos(x) = cos_poly(x)
    //   For |x| ∈ (π/4, π/2]:    sin(x) = cos_poly(π/2-|x|), cos(x) = sin_poly(π/2-|x|)
    let one = T::one();
    let frac_pi_4 = T::frac_pi_4();
    let abs_reduced = reduced.abs();
    let (poly_arg, swapped) = if abs_reduced >= frac_pi_4 {
        (frac_pi_2.saturating_sub(abs_reduced), true)
    } else {
        (abs_reduced, false)
    };
    let u = poly_arg.saturating_mul(poly_arg);

    // Evaluate sin and cos polynomials over [0, π/4] using minimax
    // (Chebyshev) coefficients. Uses multiply-by-constant instead of
    // division, avoiding cumulative rounding error from per-step divides.
    //
    // sin(x) = x + x³·P(x²)   where P = minimax poly of (sin(x)-x)/x³
    // cos(x) = 1 + x²·Q(x²)   where Q = minimax poly of (cos(x)-1)/x²
    let (sp_val, cp_val) = if T::frac_bits() >= 24 {
        // High precision: degree 15 sin, degree 14 cos
        let sp = horner(&SIN_P_HI, u);
        let sin_approx = poly_arg.saturating_add(poly_arg.saturating_mul(u).saturating_mul(sp));
        let cp = horner(&COS_Q_HI, u);
        (sin_approx, one.saturating_add(u.saturating_mul(cp)))
    } else {
        // Low precision: degree 9 sin, degree 8 cos
        let sp = horner(&SIN_P_LO, u);
        let sin_approx = poly_arg.saturating_add(poly_arg.saturating_mul(u).saturating_mul(sp));
        let cp = horner(&COS_Q_LO, u);
        (sin_approx, one.saturating_add(u.saturating_mul(cp)))
    };

    // Map back: if we swapped, sin(x) = cos_poly, cos(x) = sin_poly
    // Also restore sign of sin for negative angles.
    let (sin_unsigned, cos_val) = if swapped {
        (cp_val, sp_val)
    } else {
        (sp_val, cp_val)
    };
    let sin_val = if reduced < T::zero() {
        -sin_unsigned
    } else {
        sin_unsigned
    };

    if negate {
        (-sin_val, -cos_val)
    } else {
        (sin_val, cos_val)
    }
}

/// Sine. Accepts any angle (reduced internally).
#[inline]
#[must_use]
#[cfg_attr(feature = "verify-no-panic", no_panic::no_panic)]
pub fn sin<T: CordicNumber>(angle: T) -> T {
    sin_cos(angle).0
}

/// Cosine. Accepts any angle (reduced internally).
#[inline]
#[must_use]
#[cfg_attr(feature = "verify-no-panic", no_panic::no_panic)]
pub fn cos<T: CordicNumber>(angle: T) -> T {
    sin_cos(angle).1
}

/// Tangent. Returns `sin(angle) / cos(angle)`.
///
/// # Overflow Behavior
///
/// Tangent has poles at ±π/2, ±3π/2, etc. where it approaches ±∞.
/// Since these poles occur at irrational values that cannot be exactly
/// represented in fixed-point, this function will never divide by
/// exactly zero. However, near poles the result may:
///
/// - Saturate to `T::MAX` or `T::MIN` for very small denominators
/// - Produce very large finite values that may overflow in subsequent operations
///
/// The threshold for potential overflow is approximately:
/// - `|angle - π/2| < 2^(-frac_bits/2)` for the nearest pole
///
/// For I16F16, this means angles within ~0.004 radians of π/2 may overflow.
/// For I32F32, within ~0.00003 radians.
///
/// If you need to detect near-pole conditions, check `cos(angle).abs()`
/// against a threshold before calling `tan`.
///
/// # Example
///
/// ```
/// use fixed::types::I16F16;
/// use fixed_analytics::{tan, cos};
///
/// let angle = I16F16::from_num(1.5); // Close to π/2 ≈ 1.571
///
/// // Safe: check cosine magnitude first
/// let c = cos(angle);
/// if c.abs() > I16F16::from_num(0.01) {
///     let t = tan(angle);
///     // Use t safely
/// } else {
///     // Handle near-pole case
/// }
/// ```
#[must_use]
#[cfg_attr(feature = "verify-no-panic", no_panic::no_panic)]
pub fn tan<T: CordicNumber>(angle: T) -> T {
    let (s, c) = sin_cos(angle);
    s.div(c)
}

/// Arcsine. Domain: `[-1, 1]`. Returns angle in `[-π/2, π/2]`.
///
/// # Errors
/// Returns `DomainError` if `|x| > 1`.
#[must_use = "returns the arcsine result which should be handled"]
#[cfg_attr(feature = "verify-no-panic", no_panic::no_panic)]
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
#[cfg_attr(feature = "verify-no-panic", no_panic::no_panic)]
pub fn acos<T: CordicNumber>(x: T) -> Result<T> {
    // acos(x) = π/2 - asin(x)
    asin(x).map(|a| T::frac_pi_2().saturating_sub(a))
}

/// Arctangent. Accepts any value. Returns angle in `(-π/2, π/2)`.
#[must_use]
#[cfg_attr(feature = "verify-no-panic", no_panic::no_panic)]
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
#[cfg_attr(feature = "verify-no-panic", no_panic::no_panic)]
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
