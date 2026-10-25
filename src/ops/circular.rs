//! Trigonometric functions via circular CORDIC.

use crate::bounded::{NonNegative, UnitInterval};
use crate::error::{Error, Result};
use crate::kernel::circular_vectoring;
use crate::ops::algebraic::sqrt_nonneg;
use crate::tables::chebyshev::{COS_Q_HI, COS_Q_LO, SIN_P_HI, SIN_P_LO, horner};
use crate::traits::CordicNumber;

/// `1/(2π)` as I1F63, for the angle reduction multiply in [`sin_cos`].
const FRAC_1_2PI_I1F63: i64 = 0x145F_306D_C9C8_82A5;

/// Sine and cosine of any angle, cheaper than separate calls.
#[must_use]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn sin_cos<T: CordicNumber>(angle: T) -> (T, T) {
    let pi = T::pi();
    let frac_pi_2 = T::frac_pi_2();
    let two_pi = pi + pi;

    // Reduce to [-π, π] by n = round(angle / 2π) at any magnitude. The
    // reciprocal multiply leaves n off by at most one, which the clamp below
    // absorbs bit for bit (docs/sin_cos.md); other types divide.
    let reduced = if angle > pi || angle < -pi {
        let quotient = if T::total_bits() <= 2 * T::frac_bits() {
            angle.saturating_mul(T::from_i1f63(FRAC_1_2PI_I1F63))
        } else {
            angle.div(two_pi)
        };
        let n = quotient.round();
        // n·2π is exact and the difference small, so wrapping recovers it even
        // when n·2π exceeds MAX.
        angle.wrapping_sub(n.wrapping_mul(two_pi))
    } else {
        angle
    };

    let reduced = if reduced > pi {
        reduced.saturating_sub(two_pi)
    } else if reduced < -pi {
        reduced.saturating_add(two_pi)
    } else {
        reduced
    };

    let (reduced, negate) = if reduced > frac_pi_2 {
        (reduced - pi, true)
    } else if reduced < -frac_pi_2 {
        (reduced + pi, true)
    } else {
        (reduced, false)
    };

    // Past π/4 the polynomials take π/2 − |x| and swap roles, which avoids
    // cancellation near π/2.
    let one = T::one();
    let frac_pi_4 = T::frac_pi_4();
    let abs_reduced = reduced.abs();
    let (poly_arg, swapped) = if abs_reduced >= frac_pi_4 {
        (frac_pi_2.saturating_sub(abs_reduced), true)
    } else {
        (abs_reduced, false)
    };
    let u = poly_arg.saturating_mul(poly_arg);

    // Minimax polynomials, multiplying by constants rather than dividing at
    // each step: degree 15 and 14 from 24 fractional bits, else 9 and 8.
    let (sp_val, cp_val) = if T::frac_bits() >= 24 {
        let sp = horner(&SIN_P_HI, u);
        let sin_approx = poly_arg.saturating_add(poly_arg.saturating_mul(u).saturating_mul(sp));
        let cp = horner(&COS_Q_HI, u);
        (sin_approx, one.saturating_add(u.saturating_mul(cp)))
    } else {
        let sp = horner(&SIN_P_LO, u);
        let sin_approx = poly_arg.saturating_add(poly_arg.saturating_mul(u).saturating_mul(sp));
        let cp = horner(&COS_Q_LO, u);
        (sin_approx, one.saturating_add(u.saturating_mul(cp)))
    };

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

/// Sine of any angle.
#[inline]
#[must_use]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn sin<T: CordicNumber>(angle: T) -> T {
    sin_cos(angle).0
}

/// Cosine of any angle.
#[inline]
#[must_use]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn cos<T: CordicNumber>(angle: T) -> T {
    sin_cos(angle).1
}

/// Tangent, as `sin(angle) / cos(angle)`.
///
/// It saturates to `T::MAX` or `T::MIN` within about 4e-5 of a pole for
/// I16F16 and 5e-10 for I32F32. Check the cosine to catch near-pole angles:
///
/// ```
/// use fixed::types::I16F16;
/// use fixed_analytics::{tan, cos};
///
/// let angle = I16F16::from_num(1.5); // Close to π/2 ≈ 1.571
///
/// let c = cos(angle);
/// if c.abs() > I16F16::from_num(0.01) {
///     let t = tan(angle);
/// } else {
///     // Handle near-pole case
/// }
/// ```
#[must_use]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn tan<T: CordicNumber>(angle: T) -> T {
    let (s, c) = sin_cos(angle);
    s.div(c)
}

/// Arcsine. Domain: `[-1, 1]`. Returns an angle in `[-π/2, π/2]`.
/// # Errors
/// Returns `DomainError` if `|x| > 1`.
#[must_use = "returns the arcsine result which should be handled"]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn asin<T: CordicNumber>(x: T) -> Result<T> {
    asin_unguarded(x)
}

/// [`asin`] without the no-panic guard, which [`acos`] calls so the link
/// check sees no guarded callee inside it.
#[inline]
fn asin_unguarded<T: CordicNumber>(x: T) -> Result<T> {
    let Some(unit_x) = UnitInterval::new(x) else {
        return Err(Error::domain("asin", "value in range [-1, 1]"));
    };

    if x == T::one() {
        return Ok(T::frac_pi_2());
    }
    if x == -T::one() {
        return Ok(-T::frac_pi_2());
    }
    if x == T::zero() {
        return Ok(T::zero());
    }

    // asin(x) = atan(x / √(1 − x²)).
    let sqrt_term = sqrt_nonneg(NonNegative::one_minus_square(unit_x));

    // A root below 2^-15 puts x all but at ±1, so the result is ±π/2.
    if sqrt_term < T::from_i1f63(0x0001_0000_0000_0000) {
        return if x.is_positive() {
            Ok(T::frac_pi_2())
        } else {
            Ok(-T::frac_pi_2())
        };
    }

    Ok(atan(x.div(sqrt_term)))
}

/// Arccosine. Domain: `[-1, 1]`. Returns an angle in `[0, π]`.
/// # Errors
/// Returns `DomainError` if `|x| > 1`.
#[must_use = "returns the arccosine result which should be handled"]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn acos<T: CordicNumber>(x: T) -> Result<T> {
    asin_unguarded(x).map(|a| T::frac_pi_2().saturating_sub(a))
}

/// Arctangent of any value. Returns an angle in `(-π/2, π/2)`.
#[must_use]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn atan<T: CordicNumber>(x: T) -> T {
    let zero = T::zero();
    let one = T::one();

    if x == zero {
        return zero;
    }

    // Past 1, atan(x) = sign(x)·π/2 − atan(1/x) keeps the CORDIC convergent.
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
        circular_vectoring(one, x, zero).2
    }
}

/// Four-quadrant arctangent. Returns an angle in `[-π, π]`, and 0 for (0, 0).
#[must_use]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn atan2<T: CordicNumber>(y: T, x: T) -> T {
    let zero = T::zero();
    let pi = T::pi();
    let frac_pi_2 = T::frac_pi_2();

    if x == zero {
        return if y.is_negative() {
            -frac_pi_2
        } else if y == zero {
            zero
        } else {
            frac_pi_2
        };
    }

    if y == zero {
        return if x.is_negative() { pi } else { zero };
    }

    let (_, _, base_angle) = circular_vectoring(x.abs(), y.abs(), zero);

    match (x.is_negative(), y.is_negative()) {
        (false, false) => base_angle,
        (false, true) => -base_angle,
        (true, false) => pi - base_angle,
        (true, true) => base_angle - pi,
    }
}
