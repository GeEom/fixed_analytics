#![doc = include_str!("../../docs/hyperbolic.md")]

use crate::bounded::OpenUnitInterval;
use crate::error::{Error, Result};
use crate::eval::{HUGE, W, expm1, fixed, log1p, sinh_cosh as sinh_cosh_work, tiered, work};
use crate::real::Work;
use crate::traits::CordicNumber;

/// `|x|` in the working type, and whether it is past [`HUGE`].
fn magnitude<T: CordicNumber>(x: T) -> (W<T>, bool) {
    let a = work(x).abs();
    let huge = a > W::<T>::from_int(HUGE);
    (a, huge)
}

/// `expm1(2|x|)`, from which tanh and coth follow without cancelling.
fn expm1_double<T: CordicNumber>(a: W<T>) -> W<T> {
    expm1(a.ldexp(1))
}

/// Hyperbolic sine and cosine, cheaper than separate calls. Both saturate
/// together for large `|x|`, as the [module docs](self) describe.
#[must_use]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn sinh_cosh<T: CordicNumber>(x: T) -> (T, T) {
    if magnitude(x).1 {
        let s = if x.is_negative() {
            T::min_value()
        } else {
            T::max_value()
        };
        return (s, T::max_value());
    }
    tiered!(T, |v| {
        let (s, c) = sinh_cosh_work(v(x).abs());
        (s.with_sign(x.is_negative()), c)
    })
}

/// Hyperbolic sine. Saturates to `T::MAX` or `T::MIN` past `±asinh(T::MAX)`;
/// the [module docs](self) list the thresholds.
#[inline]
#[must_use]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn sinh<T: CordicNumber>(x: T) -> T {
    sinh_cosh(x).0
}

/// Hyperbolic cosine, always at least 1. Saturates to `T::MAX` past
/// `acosh(T::MAX)`; the [module docs](self) list the thresholds.
#[inline]
#[must_use]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn cosh<T: CordicNumber>(x: T) -> T {
    sinh_cosh(x).1
}

/// Hyperbolic tangent. Result in `[-1, 1]`.
#[must_use]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn tanh<T: CordicNumber>(x: T) -> T {
    if magnitude(x).1 {
        return if x.is_negative() { -T::one() } else { T::one() };
    }
    tiered!(T, |v| tanh_work(v(x)))
}

fn tanh_work<V: Work>(x: V) -> V {
    let u = expm1(x.abs().ldexp(1));
    (u / (u + V::from_int(2))).with_sign(x.is_neg())
}

/// Hyperbolic cotangent. Domain: `x ≠ 0`.
/// # Errors
/// Returns `DomainError` if `x = 0`.
#[must_use = "returns the hyperbolic cotangent result which should be handled"]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn coth<T: CordicNumber>(x: T) -> Result<T> {
    if x == T::zero() {
        return Err(Error::domain("coth", "non-zero value"));
    }
    let (a, huge) = magnitude(x);
    let one = W::<T>::one();
    let c = if huge {
        one
    } else {
        let u = expm1_double::<T>(a);
        (u + one.ldexp(1)) / u
    };
    Ok(fixed(c.with_sign(x.is_negative())))
}

/// Inverse hyperbolic sine. Accepts any value.
#[must_use]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn asinh<T: CordicNumber>(x: T) -> T {
    tiered!(T, |v| asinh_work(v(x)))
}

/// `asinh(z) = log1p(z + z²/(1 + √(1 + z²)))`, exact in form for any z.
fn asinh_work<V: Work>(x: V) -> V {
    let (z, one) = (x.abs(), V::one());
    let z2 = z.square();
    log1p(z + z2 / (one + (one + z2).sqrt())).with_sign(x.is_neg())
}

/// Inverse hyperbolic cosine. Domain: `x ≥ 1`.
/// # Errors
/// Returns `DomainError` if `x < 1`.
#[must_use = "returns the inverse hyperbolic cosine result which should be handled"]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn acosh<T: CordicNumber>(x: T) -> Result<T> {
    if x < T::one() {
        return Err(Error::domain("acosh", "value >= 1"));
    }
    Ok(tiered!(T, |v| acosh_work(v(x))))
}

/// `acosh(1 + t) = log1p(t + √(t(t + 2)))`, with t exact.
fn acosh_work<V: Work>(x: V) -> V {
    let t = x - V::one();
    log1p(t + (t * (t + V::from_int(2))).sqrt())
}

/// Inverse hyperbolic tangent. Domain: `(-1, 1)`.
/// # Errors
/// Returns `DomainError` if `|x| ≥ 1`.
#[must_use = "returns the inverse hyperbolic tangent result which should be handled"]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn atanh<T: CordicNumber>(x: T) -> Result<T> {
    OpenUnitInterval::new(x)
        .map(atanh_open)
        .ok_or_else(|| Error::domain("atanh", "value in range (-1, 1)"))
}

/// Inverse hyperbolic tangent of an [`OpenUnitInterval`] value, with no
/// runtime domain check.
#[must_use]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn atanh_open<T: CordicNumber>(x: OpenUnitInterval<T>) -> T {
    tiered!(T, |v| atanh_work(v(x.get())))
}

/// `atanh(a) = log1p(2a/(1 − a))/2`, with 1 − a exact.
fn atanh_work<V: Work>(x: V) -> V {
    let a = x.abs();
    log1p((a / (V::one() - a)).ldexp(1))
        .ldexp(-1)
        .with_sign(x.is_neg())
}

/// Inverse hyperbolic cotangent. Domain: `|x| > 1`.
/// # Errors
/// Returns `DomainError` if `|x| ≤ 1`.
#[must_use = "returns the inverse hyperbolic cotangent result which should be handled"]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn acoth<T: CordicNumber>(x: T) -> Result<T> {
    let v = work(x);
    let a = v.abs();
    let one = W::<T>::one();
    if a <= one {
        return Err(Error::domain("acoth", "|value| > 1"));
    }
    // acoth(a) = log1p(2/(a − 1))/2, with a − 1 exact.
    let r = log1p((one / (a - one)).ldexp(1)).ldexp(-1);
    Ok(fixed(r.with_sign(v.is_neg())))
}
