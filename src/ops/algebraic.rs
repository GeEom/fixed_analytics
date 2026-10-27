//! Roots, norms and power-of-two scaling, formed so no square or product underflows.

use crate::bounded::NonNegative;
use crate::error::{Error, Result};
use crate::eval::{W, fixed, work};
use crate::real::Work;
use crate::traits::CordicNumber;

/// Square root, rounded to the nearest representable value. Domain: `x ≥ 0`.
/// # Errors
/// Returns `DomainError` if `x < 0`.
#[must_use = "returns the square root result which should be handled"]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn sqrt<T: CordicNumber>(x: T) -> Result<T> {
    NonNegative::new(x)
        .map(sqrt_nonneg)
        .ok_or_else(|| Error::domain("sqrt", "non-negative value"))
}

/// Square root of a [`NonNegative`] value, with no runtime domain check.
///
/// The result is correctly rounded at every magnitude: raw bits
/// `⌊√(X·2^f) + ½⌋`, as [`CordicNumber::sqrt_round`] computes them.
#[must_use]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn sqrt_nonneg<T: CordicNumber>(x: NonNegative<T>) -> T {
    x.get().sqrt_round()
}

/// `√(x² + y²)`, rounded once, with no overflow or underflow in between.
#[must_use]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn hypot<T: CordicNumber>(x: T, y: T) -> T {
    let sum = work(x).square() + work(y).square();
    fixed(sum.sqrt())
}

/// The Euclidean norm `√(Σ vᵢ²)` of any number of components, rounded once.
#[must_use]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn norm<T: CordicNumber>(v: &[T]) -> T {
    let sum = v
        .iter()
        .fold(W::<T>::zero(), |acc, &x| acc + work(x).square());
    fixed(sum.sqrt())
}

/// `v / |v|`, each component rounded once, the direction exact even a few ulps long.
/// # Errors
/// Returns `DomainError` for the zero vector.
#[must_use = "returns the unit vector which should be handled"]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn normalize<T: CordicNumber, const N: usize>(v: [T; N]) -> Result<[T; N]> {
    let sum = v
        .iter()
        .fold(W::<T>::zero(), |acc, &x| acc + work(x).square());
    if sum.is_zero() {
        return Err(Error::domain("normalize", "non-zero vector"));
    }
    let inv = W::<T>::one() / sum.sqrt();
    Ok(v.map(|x| fixed(work(x) * inv)))
}

/// `⌊log₂|x|⌋`, the power of two of the leading bit, or `None` for zero.
#[must_use]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn scale_exponent<T: CordicNumber>(x: T) -> Option<i32> {
    let mag = x.raw_i128().unsigned_abs();
    if mag == 0 {
        return None;
    }
    let lead = 127 - mag.leading_zeros().cast_signed();
    Some(lead - T::frac_bits().cast_signed())
}

/// `x · 2^e`, exact unless it rounds below an ulp or saturates.
#[must_use]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn scale_by<T: CordicNumber>(x: T, e: i32) -> T {
    fixed(work(x).ldexp(e.clamp(-512, 512)))
}

/// `√(a·b)`, rounded once, with the product formed without underflow.
/// # Errors
/// Returns `DomainError` if `a` and `b` have opposite signs.
#[must_use = "returns the geometric mean which should be handled"]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn sqrt_product<T: CordicNumber>(a: T, b: T) -> Result<T> {
    let (wa, wb) = (work(a), work(b));
    if wa.is_neg() != wb.is_neg() && !wa.is_zero() && !wb.is_zero() {
        return Err(Error::domain("sqrt_product", "arguments of one sign"));
    }
    Ok(fixed((wa * wb).sqrt()))
}

/// Reciprocal square root `1/√x`, rounded once. Domain: `x > 0`.
/// # Errors
/// Returns `DomainError` if `x ≤ 0`.
#[must_use = "returns the reciprocal square root which should be handled"]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn rsqrt<T: CordicNumber>(x: T) -> Result<T> {
    if x <= T::zero() {
        return Err(Error::domain("rsqrt", "positive value"));
    }
    Ok(fixed(W::<T>::one() / work(x).sqrt()))
}
