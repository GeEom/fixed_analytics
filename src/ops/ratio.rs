//! Ratios finite where their parts vanish, to full relative precision near 0.

use crate::bounded::OpenUnitInterval;
use crate::error::{Error, Result};
use crate::eval::{
    HUGE, atan_ratio as atan_of_ratio, beyond, expm1, log1p, sin_cos, stumpff as stumpff_work,
    tiered,
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
    tiered!(T, |v| stumpff_work(v(z)))
}

/// `(tan x / x, cos x)`, the first 1 at 0 and saturating at the poles.
#[must_use]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn tan_ratio<T: CordicNumber>(x: T) -> (T, T) {
    tiered!(T, |v| tan_ratio_work(v(x)))
}

fn tan_ratio_work<V: Work>(v: V) -> (V, V) {
    let (s, c) = sin_cos(v);
    let ratio = if v.is_zero() {
        V::one()
    } else if c.is_zero() {
        beyond::<V>().with_sign(s.is_neg() != v.is_neg())
    } else {
        s / (v * c)
    };
    (ratio, c)
}

/// `atan z / z`, which is 1 at 0.
#[must_use]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn atan_ratio<T: CordicNumber>(z: T) -> T {
    if z == T::zero() {
        return T::one();
    }
    tiered!(T, |v| atan_ratio_work(v(z)))
}

fn atan_ratio_work<V: Work>(z: V) -> V {
    let a = z.abs();
    atan_of_ratio(a, V::one()) / a
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
    if z == T::zero() {
        return Ok(T::one());
    }
    Ok(tiered!(T, |v| atanh_ratio_work(v(z))))
}

fn atanh_ratio_work<V: Work>(z: V) -> V {
    let a = z.abs();
    log1p((a / (V::one() - a)).ldexp(1)).ldexp(-1) / a
}

/// `tanh x / x`, which is 1 at 0.
#[must_use]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn tanh_ratio<T: CordicNumber>(x: T) -> T {
    if x == T::zero() {
        return T::one();
    }
    tiered!(T, |v| tanh_ratio_work(v(x)))
}

fn tanh_ratio_work<V: Work>(x: V) -> V {
    let a = x.abs();
    if a > V::from_int(HUGE) {
        return V::one() / a;
    }
    let u = expm1(a.ldexp(1));
    u / ((u + V::from_int(2)) * a)
}
