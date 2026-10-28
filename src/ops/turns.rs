//! Angles in turns, and radians wrapped by exact multiples of 2π.

use crate::consts::FRAC_1_2PI;
use crate::eval::{W, fixed, pi, sin_cos_octant, tiered, turn_fraction, work};
use crate::real::Work;
use crate::traits::CordicNumber;

/// Sine and cosine of `t` turns, `t·2π` radians, of any magnitude.
#[must_use]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn sin_cos_turns<T: CordicNumber>(t: T) -> (T, T) {
    if crate::fast::applies::<T>() {
        let raw = t.raw_i128();
        let (s, c) = crate::fast::sin_cos(raw.unsigned_abs(), T::frac_bits(), false);
        let s = if raw < 0 { -s } else { s };
        return (crate::fast::round(s), crate::fast::round(c));
    }
    tiered!(T, |v| turns_work(t, v))
}

fn turns_work<T: CordicNumber, V: Work>(t: T, _as_work: fn(T) -> V) -> (V, V) {
    let raw = t.raw_i128();
    let (oct, w) = V::octant(raw.unsigned_abs(), T::frac_bits(), false);
    let (s, c) = sin_cos_octant(oct, w);
    (if raw < 0 { -s } else { s }, c)
}

/// The four-quadrant arctangent in turns, `atan2(y, x)/2π` in `(−½, ½]`, rounded once.
#[must_use]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn atan2_turns<T: CordicNumber>(y: T, x: T) -> T {
    tiered!(T, |v| atan2_turns_work(v(y), v(x)))
}

fn atan2_turns_work<V: Work>(y: V, x: V) -> V {
    crate::cx::atan2(y, x) * V::constant(&FRAC_1_2PI)
}

/// `x/τ`, radians to turns, rounded once.
#[must_use]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn div_by_tau<T: CordicNumber>(x: T) -> T {
    tiered!(T, |v| by_constant(v(x), &FRAC_1_2PI, 0))
}

fn by_constant<V: Work>(x: V, c: &crate::consts::C, k: i32) -> V {
    x * V::constant(c).ldexp(k)
}

/// `x·τ`, turns to radians, rounded once.
#[must_use]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn mul_by_tau<T: CordicNumber>(x: T) -> T {
    tiered!(T, |v| by_constant(v(x), &crate::consts::PI, 1))
}

/// `r` rounded into `[low, high]`, else to the end nearer it modulo 2π.
fn round_into<T: CordicNumber>(r: W<T>, low: W<T>, high: W<T>) -> T {
    let t = fixed::<T>(r);
    if work(t) >= low && work(t) <= high {
        return t;
    }
    let tau = pi::<W<T>>().ldexp(1);
    let ulp = T::from_raw_i128(1);
    // The values of T just inside each end.
    let first = fixed::<T>(low);
    let first = if work(first) < low {
        first.saturating_add(ulp)
    } else {
        first
    };
    let last = fixed::<T>(high);
    let last = if work(last) > high {
        last.saturating_sub(ulp)
    } else {
        last
    };
    let away = |c: T| {
        let d = (work(c) - r).abs();
        if d > tau.ldexp(-1) { tau - d } else { d }
    };
    if away(first) <= away(last) {
        first
    } else {
        last
    }
}

/// `x` less the multiple of 2π that leaves it in `[0, 2π)`.
#[must_use]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn wrap_radians<T: CordicNumber>(x: T) -> T {
    let tau = pi::<W<T>>().ldexp(1);
    let below = tau - tau.ldexp(-(W::<T>::PREC.cast_signed()));
    round_into::<T>(turn_fraction(x) * tau, W::<T>::zero(), below)
}

/// `x` less the multiple of 2π that leaves it in `(−π, π]`.
#[must_use]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn wrap_principal<T: CordicNumber>(x: T) -> T {
    let one = W::<T>::one();
    let u = turn_fraction(x);
    let u = if u > one.ldexp(-1) { u - one } else { u };
    let half = pi::<W<T>>();
    let above = -half + half.ldexp(-(W::<T>::PREC.cast_signed()));
    round_into::<T>(u * half.ldexp(1), above, half)
}
