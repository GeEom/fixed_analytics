//! Ratios finite where their parts vanish, to full relative precision near 0.

use crate::bounded::OpenUnitInterval;
use crate::error::{Error, Result};
use crate::eval::{
    W, atan_ratio as atan_of_ratio, expm1, fixed, log1p, sin_cos, stumpff as stumpff_work, work,
};
use crate::real::Work;
use crate::traits::CordicNumber;

/// The Stumpff functions `[c₀, c₁, c₂, c₃]`, `cₖ(z) = Σ (−z)ⁿ/(2n + k)!`, of either sign.
#[must_use]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn stumpff<T: CordicNumber>(z: T) -> [T; 4] {
    stumpff_work(work(z)).map(fixed)
}

/// `(tan x / x, cos x)`, the first 1 at 0 and saturating at the poles.
#[must_use]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn tan_ratio<T: CordicNumber>(x: T) -> (T, T) {
    let v = work(x);
    let (s, c) = sin_cos(v);
    let ratio = if v.is_zero() {
        W::<T>::one()
    } else if c.is_zero() {
        crate::eval::beyond::<W<T>>().with_sign(s.is_neg() != v.is_neg())
    } else {
        s / (v * c)
    };
    (fixed(ratio), fixed(c))
}

/// `atan z / z`, which is 1 at 0.
#[must_use]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn atan_ratio<T: CordicNumber>(z: T) -> T {
    let a = work(z).abs();
    if a.is_zero() {
        return T::one();
    }
    fixed(atan_of_ratio(a, W::<T>::one()) / a)
}

/// `atanh z / z`, which is 1 at 0. Domain: `(-1, 1)`.
/// # Errors
/// Returns `DomainError` if `|z| ≥ 1`.
#[must_use = "returns the ratio which should be handled"]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn atanh_ratio<T: CordicNumber>(z: T) -> Result<T> {
    if OpenUnitInterval::new(z).is_none() {
        return Err(Error::domain("atanh_ratio", "value in range (-1, 1)"));
    }
    let a = work(z).abs();
    if a.is_zero() {
        return Ok(T::one());
    }
    let atanh = log1p((a / (W::<T>::one() - a)).ldexp(1)).ldexp(-1);
    Ok(fixed(atanh / a))
}

/// `tanh x / x`, which is 1 at 0.
#[must_use]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn tanh_ratio<T: CordicNumber>(x: T) -> T {
    let a = work(x).abs();
    if a.is_zero() {
        return T::one();
    }
    if a > W::<T>::from_int(crate::eval::HUGE) {
        return fixed(W::<T>::one() / a);
    }
    let u = expm1(a.ldexp(1));
    fixed(u / ((u + W::<T>::from_int(2)) * a))
}
