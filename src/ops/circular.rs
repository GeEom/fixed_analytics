//! Trigonometric functions.

use crate::bounded::UnitInterval;
use crate::error::{Error, Result};
use crate::eval::{W, atan_ratio, beyond, fixed, pi, sin_cos_octant, tiered, work};
use crate::real::Work;
use crate::traits::CordicNumber;

/// Sine and cosine of `angle`, reduced by 2π exactly, in `as_work`'s working type.
fn sin_cos_work<T: CordicNumber, V: Work>(angle: T, _as_work: fn(T) -> V) -> (V, V) {
    let raw = angle.raw_i128();
    let (oct, w) = V::octant(raw.unsigned_abs(), T::frac_bits(), true);
    let (s, c) = sin_cos_octant(oct, w);
    (if raw < 0 { -s } else { s }, c)
}

/// `√(1 − a²)` for `a` in `[0, 1]`, from `(1 − a)(1 + a)` so it keeps its
/// bits near 1.
fn cofactor<V: Work>(a: V) -> V {
    let one = V::one();
    ((one - a) * (one + a)).sqrt()
}

/// Sine and cosine of any angle, cheaper than separate calls.
#[must_use]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn sin_cos<T: CordicNumber>(angle: T) -> (T, T) {
    if crate::fast::applies::<T>() {
        let raw = angle.raw_i128();
        let (s, c) = crate::fast::sin_cos(raw.unsigned_abs(), T::frac_bits(), true);
        let s = if raw < 0 { -s } else { s };
        return (crate::fast::round(s), crate::fast::round(c));
    }
    tiered!(T, |v| sin_cos_work(angle, v))
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

/// Tangent, as `sin(angle) / cos(angle)` divided before rounding.
///
/// It saturates to `T::MAX` or `T::MIN` only where the tangent itself is out
/// of range. Check the cosine to catch near-pole angles:
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
    tiered!(T, |v| tan_work(angle, v))
}

fn tan_work<T: CordicNumber, V: Work>(angle: T, as_work: fn(T) -> V) -> V {
    let (s, c) = sin_cos_work(angle, as_work);
    if c.is_zero() {
        beyond::<V>().with_sign(s.is_neg())
    } else {
        s / c
    }
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
    if UnitInterval::new(x).is_none() {
        return Err(Error::domain("asin", "value in range [-1, 1]"));
    }
    Ok(tiered!(T, |v| {
        let a = v(x).abs();
        atan_ratio(a, cofactor(a)).with_sign(x.is_negative())
    }))
}

/// `asin(1 − gap)` from the gap to 1, which keeps all of its bits where
/// `1 − gap` would round. Domain: `gap` in `[0, 2]`.
/// # Errors
/// Returns `DomainError` if `gap` is outside `[0, 2]`.
#[must_use = "returns the arcsine result which should be handled"]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn asin_from_gap<T: CordicNumber>(gap: T) -> Result<T> {
    let g = work(gap);
    let two = W::<T>::from_int(2);
    if g.is_neg() || g > two {
        return Err(Error::domain("asin_from_gap", "gap in range [0, 2]"));
    }
    let a = W::<T>::one() - g;
    let s = (g * (two - g)).sqrt();
    Ok(fixed(atan_ratio(a.abs(), s).with_sign(a.is_neg())))
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
    if UnitInterval::new(x).is_none() {
        return Err(Error::domain("acos", "value in range [-1, 1]"));
    }
    Ok(tiered!(T, |v| acos_work(v(x))))
}

fn acos_work<V: Work>(x: V) -> V {
    let a = x.abs();
    let near = atan_ratio(cofactor(a), a);
    if x.is_neg() { pi::<V>() - near } else { near }
}

/// Arctangent of any value. Returns an angle in `[-π/2, π/2]`.
#[must_use]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn atan<T: CordicNumber>(x: T) -> T {
    tiered!(T, |v| atan_work(v(x)))
}

fn atan_work<V: Work>(x: V) -> V {
    let (a, one) = (x.abs(), V::one());
    atan_ratio(a, one).with_sign(x.is_neg())
}

/// Four-quadrant arctangent. Returns an angle in `[-π, π]`: `π` along the
/// negative x-axis, and 0 for (0, 0).
#[must_use]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn atan2<T: CordicNumber>(y: T, x: T) -> T {
    tiered!(T, |v| crate::cx::atan2(v(y), v(x)))
}
