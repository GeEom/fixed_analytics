//! Intervals certain to hold the exact value, by [`ERROR_ULPS`](crate::ERROR_ULPS).
//! An end at `T::MIN` or `T::MAX` stands for everything past it.

use crate::eval::{W, pi, turn_fraction, work};
use crate::ops::circular::{atan, sin_cos};
use crate::ops::exponential::exp;
use crate::ops::hyperbolic::tanh;
use crate::ops::turns::sin_cos_turns;
use crate::real::Work;
use crate::traits::CordicNumber;

/// `[y − ulp, y + ulp]`, which holds the exact value of any function giving `y`.
#[must_use]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn enclose<T: CordicNumber>(y: T) -> (T, T) {
    let ulp = T::from_raw_i128(1);
    (y.saturating_sub(ulp), y.saturating_add(ulp))
}

/// `(lo, hi)` ordered.
fn ordered<T: CordicNumber>(a: T, b: T) -> (T, T) {
    if a <= b { (a, b) } else { (b, a) }
}

fn lesser<T: CordicNumber>(a: T, b: T) -> T {
    ordered(a, b).0
}

fn greater<T: CordicNumber>(a: T, b: T) -> T {
    ordered(a, b).1
}

/// `[f(lo) − ulp, f(hi) + ulp]` for a rising `f`, clamped to its range.
fn rising<T: CordicNumber>(f: fn(T) -> T, lo: T, hi: T, floor: T, ceiling: T) -> (T, T) {
    let (lo, hi) = ordered(lo, hi);
    (
        greater(enclose(f(lo)).0, floor),
        lesser(enclose(f(hi)).1, ceiling),
    )
}

/// The range of `e^x` over `[lo, hi]`, rounded outward.
#[must_use]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn exp_range<T: CordicNumber>(lo: T, hi: T) -> (T, T) {
    rising(exp, lo, hi, T::zero(), T::max_value())
}

/// The range of `atan x` over `[lo, hi]`, rounded outward.
#[must_use]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn atan_range<T: CordicNumber>(lo: T, hi: T) -> (T, T) {
    rising(atan, lo, hi, T::min_value(), T::max_value())
}

/// The range of `tanh x` over `[lo, hi]`, rounded outward.
#[must_use]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn tanh_range<T: CordicNumber>(lo: T, hi: T) -> (T, T) {
    rising(tanh, lo, hi, -T::one(), T::one())
}

/// The ranges `((sin_min, sin_max), (cos_min, cos_max))` over `[lo, hi]`,
/// rounded outward; an extreme within the working error of an end counts.
#[must_use]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn sin_cos_range<T: CordicNumber>(lo: T, hi: T) -> ((T, T), (T, T)) {
    let (lo, hi) = ordered(lo, hi);
    let span = (work(hi) - work(lo)) / pi::<W<T>>().ldexp(1);
    range_in_turns(span, turn_fraction(lo), sin_cos(lo), sin_cos(hi))
}

/// [`sin_cos_range`] over `[lo, hi]` in turns, each end reduced exactly.
#[must_use]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn sin_cos_range_turns<T: CordicNumber>(lo: T, hi: T) -> ((T, T), (T, T)) {
    let (lo, hi) = ordered(lo, hi);
    let span = work(hi) - work(lo);
    let start = lo.raw_i128().rem_euclid(T::one().raw_i128());
    let start = W::<T>::from_raw(start, T::frac_bits());
    range_in_turns(span, start, sin_cos_turns(lo), sin_cos_turns(hi))
}

/// The ranges over `span` turns from `start` past a whole one, given the ends' sines and cosines.
fn range_in_turns<T: CordicNumber>(
    span: W<T>,
    start: W<T>,
    (sl, cl): (T, T),
    (sh, ch): (T, T),
) -> ((T, T), (T, T)) {
    let one = T::one();
    if span >= W::<T>::one() {
        return ((-one, one), (-one, one));
    }
    let slack = W::<T>::one().ldexp(16 - W::<T>::PREC.cast_signed());
    let reaches = |quarter: i64| {
        let at = W::<T>::from_int(quarter).ldexp(-2);
        let ahead = at - start;
        let ahead = if ahead.is_neg() {
            ahead + W::<T>::one()
        } else {
            ahead
        };
        ahead <= span + slack || ahead >= W::<T>::one() - slack
    };
    let band = |a: T, b: T, top: bool, bottom: bool| {
        let low = if bottom {
            -one
        } else {
            greater(enclose(lesser(a, b)).0, -one)
        };
        let high = if top {
            one
        } else {
            lesser(enclose(greater(a, b)).1, one)
        };
        (low, high)
    };
    (
        band(sl, sh, reaches(1), reaches(3)),
        band(cl, ch, reaches(0), reaches(2)),
    )
}
