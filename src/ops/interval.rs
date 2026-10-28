//! Intervals certain to hold the exact value, by [`ERROR_ULPS`](crate::ERROR_ULPS).
//! An end at `T::MIN` or `T::MAX` stands for everything past it.

use crate::error::Result;
use crate::eval::{W, pi, stumpff as stumpff_work, turn_fraction, work};
use crate::ops::algebraic::hypot;
use crate::ops::arith::sqrt_or_zero;
use crate::ops::circular::{atan, atan2, sin_cos};
use crate::ops::exponential::{
    exp, exp_moments, exprel, exprel2, log1p, log1p_ratio, log1p_residual_ratio,
};
use crate::ops::hyperbolic::tanh;
use crate::ops::ratio::{atan_ratio, atanh_ratio, stumpff, tan_ratio, tanh_ratio};
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

/// `[a − ulp, b + ulp]` within `[floor, ceiling]`, holding exact values within an ulp of `a ≤ b`.
fn hull<T: CordicNumber>(a: T, b: T, floor: T, ceiling: T) -> (T, T) {
    (greater(enclose(a).0, floor), lesser(enclose(b).1, ceiling))
}

/// `[f(lo) − ulp, f(hi) + ulp]` for a rising `f`, clamped to its range.
fn rising<T: CordicNumber>(f: fn(T) -> T, lo: T, hi: T, floor: T, ceiling: T) -> (T, T) {
    let (lo, hi) = ordered(lo, hi);
    hull(f(lo), f(hi), floor, ceiling)
}

/// The range over `[lo, hi]`, ordered, of an even function monotone in `|x|`, 1 at 0, from its ends' values.
fn even<T: CordicNumber>(
    (lo, hi): (T, T),
    (at_lo, at_hi): (T, T),
    rises: bool,
    clamp: (T, T),
) -> (T, T) {
    let zero = T::zero();
    let far = enclose(
        if lo.raw_i128().unsigned_abs() > hi.raw_i128().unsigned_abs() {
            at_lo
        } else {
            at_hi
        },
    );
    let near = if lo <= zero && hi >= zero {
        (T::one(), T::one())
    } else {
        enclose(if lo > zero { at_lo } else { at_hi })
    };
    let (low, high) = if rises {
        (near.0, far.1)
    } else {
        (far.0, near.1)
    };
    (greater(low, clamp.0), lesser(high, clamp.1))
}

/// The greatest value of `T` at most `v`, then an ulp below it, for a `v` within an ulp.
fn below<T: CordicNumber>(v: W<T>) -> T {
    let (lo, hi) = (T::min_value().raw_i128(), T::max_value().raw_i128());
    T::from_raw_i128(v.to_raw(T::frac_bits(), lo, hi, true)).saturating_sub(T::from_raw_i128(1))
}

fn above<T: CordicNumber>(v: W<T>) -> T {
    T::zero().saturating_sub(below::<T>(-v))
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

/// The range of `e^x` past its first terms, [`exprel`], over `[lo, hi]`, rounded outward.
#[must_use]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn exprel_range<T: CordicNumber>(lo: T, hi: T) -> (T, T) {
    rising(exprel, lo, hi, T::zero(), T::max_value())
}

/// The range of [`exprel2`] over `[lo, hi]`, rounded outward.
#[must_use]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn exprel2_range<T: CordicNumber>(lo: T, hi: T) -> (T, T) {
    rising(exprel2, lo, hi, T::zero(), T::max_value())
}

/// The ranges of [`exp_moments`] over `[lo, hi]`, each falling, rounded outward.
#[must_use]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn exp_moments_range<T: CordicNumber>(lo: T, hi: T) -> [(T, T); 6] {
    let (lo, hi) = ordered(lo, hi);
    let (low, high) = (exp_moments(hi), exp_moments(lo));
    let mut out = [(T::zero(), T::zero()); 6];
    for ((slot, a), b) in out.iter_mut().zip(low).zip(high) {
        *slot = hull(a, b, T::zero(), T::max_value());
    }
    out
}

/// The ranges of [`log1p`] and its ratios over `[lo, hi]`, rounded outward.
fn log1p_family<T: CordicNumber>(
    f: fn(T) -> Result<T>,
    lo: T,
    hi: T,
    rises: bool,
) -> Result<(T, T)> {
    let (lo, hi) = ordered(lo, hi);
    let (a, b) = (f(lo)?, f(hi)?);
    let (min, max) = (T::min_value(), T::max_value());
    Ok(if rises {
        hull(a, b, min, max)
    } else {
        hull(b, a, min, max)
    })
}

/// The range of [`log1p`] over `[lo, hi]`, rounded outward.
/// # Errors
/// Returns `DomainError` unless both ends are above −1.
#[must_use = "returns the range which should be handled"]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn log1p_range<T: CordicNumber>(lo: T, hi: T) -> Result<(T, T)> {
    log1p_family(log1p, lo, hi, true)
}

/// The range of [`log1p_ratio`], which falls, over `[lo, hi]`, rounded outward.
/// # Errors
/// Returns `DomainError` unless both ends are above −1.
#[must_use = "returns the range which should be handled"]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn log1p_ratio_range<T: CordicNumber>(lo: T, hi: T) -> Result<(T, T)> {
    log1p_family(log1p_ratio, lo, hi, false)
}

/// The range of [`log1p_residual_ratio`], which rises, over `[lo, hi]`, rounded outward.
/// # Errors
/// Returns `DomainError` unless both ends are above −1.
#[must_use = "returns the range which should be handled"]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn log1p_residual_ratio_range<T: CordicNumber>(lo: T, hi: T) -> Result<(T, T)> {
    log1p_family(log1p_residual_ratio, lo, hi, true)
}

/// The range of [`tanh_ratio`] over `[lo, hi]`, rounded outward.
#[must_use]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn tanh_ratio_range<T: CordicNumber>(lo: T, hi: T) -> (T, T) {
    let (lo, hi) = ordered(lo, hi);
    even(
        (lo, hi),
        (tanh_ratio(lo), tanh_ratio(hi)),
        false,
        (T::zero(), T::one()),
    )
}

/// The range of [`atan_ratio`] over `[lo, hi]`, rounded outward.
#[must_use]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn atan_ratio_range<T: CordicNumber>(lo: T, hi: T) -> (T, T) {
    let (lo, hi) = ordered(lo, hi);
    even(
        (lo, hi),
        (atan_ratio(lo), atan_ratio(hi)),
        false,
        (T::zero(), T::one()),
    )
}

/// The range of [`atanh_ratio`] over `[lo, hi]`, rounded outward.
/// # Errors
/// Returns `DomainError` unless both ends are in `(-1, 1)`.
#[must_use = "returns the range which should be handled"]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn atanh_ratio_range<T: CordicNumber>(lo: T, hi: T) -> Result<(T, T)> {
    let (lo, hi) = ordered(lo, hi);
    let ends = (atanh_ratio(lo)?, atanh_ratio(hi)?);
    Ok(even((lo, hi), ends, true, (T::one(), T::max_value())))
}

/// The range of `tan x / x` over `[lo, hi]`, rounded outward, everything if a pole may lie in it.
#[must_use]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn tan_ratio_range<T: CordicNumber>(lo: T, hi: T) -> (T, T) {
    let (lo, hi) = ordered(lo, hi);
    let ((at_lo, cos_lo), (at_hi, cos_hi)) = (tan_ratio(lo), tan_ratio(hi));
    let zero = T::zero();
    // Between poles cos keeps its sign, and two poles are π apart.
    if cos_lo == zero
        || cos_hi == zero
        || cos_lo.is_negative() != cos_hi.is_negative()
        || work(hi) - work(lo) >= pi()
    {
        return (T::min_value(), T::max_value());
    }
    let (min, max) = (T::min_value(), T::max_value());
    if lo <= zero && hi >= zero {
        (T::one(), lesser(enclose(greater(at_lo, at_hi)).1, max))
    } else if hi <= zero {
        hull(at_hi, at_lo, min, max)
    } else {
        hull(at_lo, at_hi, min, max)
    }
}

/// The range of `√x`, 0 for `x ≤ 0` like [`sqrt_or_zero`], over `[lo, hi]`, rounded outward.
#[must_use]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn sqrt_range<T: CordicNumber>(lo: T, hi: T) -> (T, T) {
    rising(sqrt_or_zero, lo, hi, T::zero(), T::max_value())
}

/// The end of `[lo, hi]` nearest 0, or 0 inside it, and the farthest.
fn near_far<T: CordicNumber>((a, b): (T, T)) -> (T, T) {
    let (lo, hi) = ordered(a, b);
    let zero = T::zero();
    let near = if lo > zero {
        lo
    } else if hi < zero {
        hi
    } else {
        zero
    };
    let far = if lo.raw_i128().unsigned_abs() > hi.raw_i128().unsigned_abs() {
        lo
    } else {
        hi
    };
    (near, far)
}

/// The range of [`hypot`] over the box `x × y`, rounded outward.
#[must_use]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn hypot_range<T: CordicNumber>(x: (T, T), y: (T, T)) -> (T, T) {
    let ((nx, fx), (ny, fy)) = (near_far(x), near_far(y));
    hull(hypot(nx, ny), hypot(fx, fy), T::zero(), T::max_value())
}

/// The range of [`atan2`] over the box `y × x`, rounded outward: `[−π, π]` when the box
/// holds the origin or crosses the negative x-axis, else its corners' angles.
#[must_use]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn atan2_range<T: CordicNumber>(y: (T, T), x: (T, T)) -> (T, T) {
    let ((ylo, yhi), (xlo, xhi)) = (ordered(y.0, y.1), ordered(x.0, x.1));
    let zero = T::zero();
    let half_turn = above::<T>(pi());
    let spans_zero = ylo <= zero && yhi >= zero;
    if spans_zero && (xlo < zero && ylo < zero || xhi >= zero && xlo <= zero) {
        return (-half_turn, half_turn);
    }
    let corners = [
        atan2(ylo, xlo),
        atan2(ylo, xhi),
        atan2(yhi, xlo),
        atan2(yhi, xhi),
    ];
    let low = corners.iter().copied().fold(T::max_value(), lesser);
    let high = corners.iter().copied().fold(T::min_value(), greater);
    hull(low, high, -half_turn, half_turn)
}

/// The ranges of [`stumpff`]'s `[c₀, c₁, c₂, c₃]` over `[lo, hi]`, rounded outward.
///
/// Each falls to `z = 4`, `c₃` everywhere; past it `|cₖ''| ≤ 2/(k + 4)!` rules an extreme out from
/// the slopes at the ends, or puts it within `2/(k + 4)!·w²/8` of them, `|c₀| ≤ 1`,
/// `|c₁| ≤ 1/√(1 + z)` and `0 ≤ c₂ ≤ 2/(4 + z)`.
#[must_use]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn stumpff_range<T: CordicNumber>(lo: T, hi: T) -> [(T, T); 4] {
    let (lo, hi) = ordered(lo, hi);
    let (at_lo, at_hi) = (stumpff(lo), stumpff(hi));
    let pick = |c: [T; 4], k: usize| c.get(k).copied().unwrap_or_else(T::zero);
    let mut out = [0, 1, 2, 3].map(|k| {
        hull(
            pick(at_hi, k),
            pick(at_lo, k),
            T::min_value(),
            T::max_value(),
        )
    });
    let (four, top) = (W::<T>::from_int(4), work(hi));
    if top <= four {
        return out;
    }
    let past = work(lo) >= four;
    let start = work(lo).max(four);
    let (from, to) = (stumpff_work(start), stumpff_work(top));
    let width = top - start;
    let at = |c: &[W<T>; 4], k: usize| c.get(k).copied().unwrap_or_else(W::<T>::zero);
    // 2z·cₖ' = cₖ₋₁ − k·cₖ, and c₀' = −c₁/2.
    let slope = |c: &[W<T>; 4], z: W<T>, k: usize| match k {
        0 => -at(c, 1).ldexp(-1),
        1 => (at(c, 0) - at(c, 1)) / z.ldexp(1),
        _ => (at(c, 1) - at(c, 2).ldexp(1)) / z.ldexp(1),
    };
    let slack = W::<T>::one().ldexp(16 - W::<T>::PREC.cast_signed());
    for (k, (slot, factorial)) in out.iter_mut().zip([24, 120, 720]).enumerate() {
        let bend = W::<T>::from_int(2).div_int(factorial);
        let first = if past {
            enclose(pick(at_lo, k))
        } else {
            (below::<T>(at(&from, k)), above::<T>(at(&from, k)))
        };
        let last = enclose(pick(at_hi, k));
        let (mut low, mut high) = (lesser(first.0, last.0), greater(first.1, last.1));
        if slope(&from, start, k).abs() + slope(&to, top, k).abs() <= bend * width + slack {
            let spread = above::<T>(bend * width.square().ldexp(-3));
            // c₁'s extremes past 4 are cos s at tan s = s, at least −0.21724, at most 0.1284.
            let fall = |c: i64| W::<T>::from_int(c).div_int(10_000);
            let root = W::<T>::one() / (W::<T>::one() + start).sqrt();
            let (floor, peak) = match k {
                0 => (-W::<T>::one(), W::<T>::one()),
                1 => (-root.min(fall(2173)), root.min(fall(1284))),
                _ => (W::<T>::zero(), W::<T>::from_int(2) / (four + start)),
            };
            low = greater(low.saturating_sub(spread), lesser(low, below::<T>(floor)));
            high = lesser(high.saturating_add(spread), greater(high, above::<T>(peak)));
        }
        // Below 4 it falls from its value at `lo`.
        *slot = (low, if past { high } else { greater(high, slot.1) });
    }
    out
}
