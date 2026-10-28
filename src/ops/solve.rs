//! A bracketed root finder, and the universal-variable Kepler equation.

use crate::error::{Error, Result};
use crate::eval::{HUGE, W, beyond, fixed, stumpff, work};
use crate::real::Work;
use crate::traits::CordicNumber;

/// A bracket `(a, b)` around a root of `f` in `[lo, hi]`, where `f(lo)` and `f(hi)` differ in sign.
///
/// ITP (Oliveira and Takahashi, 2020) is never slower than bisection: from `f(lo)` and `f(hi)`
/// at most `budget` more calls make `a` and `b`, signed like them, adjacent; an exact 0 ends both.
/// # Errors
/// Returns `DomainError` if `f(lo)` and `f(hi)` are non-zero with one sign.
#[must_use = "returns the bracket which should be handled"]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn root_in<T: CordicNumber, F: FnMut(T) -> T>(
    mut func: F,
    lo: T,
    hi: T,
    budget: u32,
) -> Result<(T, T)> {
    let (lo, hi) = if lo <= hi { (lo, hi) } else { (hi, lo) };
    let (mut ya, mut yb) = (func(lo), func(hi));
    let zero = T::zero();
    if ya == zero {
        return Ok((lo, lo));
    }
    if yb == zero {
        return Ok((hi, hi));
    }
    if ya.is_negative() == yb.is_negative() {
        return Err(Error::domain(
            "root_in",
            "f(lo) and f(hi) of opposite signs",
        ));
    }
    let (mut a, mut b) = (lo.raw_i128(), hi.raw_i128());
    let at = |raw: i128| W::<T>::from_raw(raw, 0);
    let width = |from: i128, to: i128| to.abs_diff(from);
    // ITP with κ₁ = 0.2/(b₀ − a₀), κ₂ = 2, n₀ = 1 and ε half an ulp.
    let start = at(b) - at(a);
    let kappa = W::<T>::one().div_int(5) / start;
    let n_max = 129 - width(a, b).leading_zeros().min(128);
    for j in 0..budget {
        let span = width(a, b);
        if span <= 1 {
            break;
        }
        let half = a.saturating_add_unsigned(span / 2);
        let (wa, wb, wh) = (at(a), at(b), at(half));
        let (fa, fb) = (work(ya), work(yb));
        let falsi = wa + (wb - wa) * (fa / (fa - fb));
        let toward = wh - falsi;
        let delta = kappa * (wb - wa).square();
        let truncated = if delta <= toward.abs() {
            falsi + delta.with_sign(toward.is_neg())
        } else {
            wh
        };
        let reach =
            W::<T>::one().ldexp(n_max.saturating_sub(j + 1).cast_signed()) - (wb - wa).ldexp(-1);
        let x = if (truncated - wh).abs() <= reach {
            truncated
        } else {
            wh - reach.with_sign(toward.is_neg())
        };
        let raw = x.to_raw(0, a.saturating_add(1), b.saturating_sub(1), true);
        let point = T::from_raw_i128(raw);
        let y = func(point);
        if y == zero {
            return Ok((point, point));
        }
        if y.is_negative() == ya.is_negative() {
            (a, ya) = (raw, y);
        } else {
            (b, yb) = (raw, y);
        }
    }
    Ok((T::from_raw_i128(a), T::from_raw_i128(b)))
}

fn kepler_residual<V: Work>(chi: V, r0: V, sigma0: V, alpha: V, tau: V) -> (V, V) {
    let x2 = chi.square();
    let z = alpha * x2;
    // Past HUGE² e^√|z| swamps every term, its coefficient on a real orbit of χ's sign.
    let huge = V::from_int(HUGE);
    if z < -huge.square() {
        return (beyond::<V>().with_sign(chi.is_neg()), beyond());
    }
    let [c0, c1, c2, c3] = stumpff(z);
    let beta = V::one() - alpha * r0;
    let f = r0 * chi + sigma0 * x2 * c2 + beta * x2 * chi * c3 - tau;
    let r = r0 * c0 + sigma0 * chi * c1 + x2 * c2;
    (f, r)
}

/// The universal anomaly `χ` after time `t` on a two-body orbit.
///
/// It solves `√μ·t = r₀χ + σ₀χ²c₂(αχ²) + (1 − αr₀)χ³c₃(αχ²)`, rising in `χ`, for `r₀ = |r⃗₀|`,
/// `σ₀ = r⃗₀·v⃗₀/√μ`, `α = 2/r₀ − v₀²/μ` and `τ = √μ·t`, by at most 256 steps of safeguarded Newton.
/// # Errors
/// Returns `DomainError` unless `r₀ > 0` and `r₀(2 − αr₀) ≥ σ₀²`, a real orbit.
#[must_use = "returns the anomaly which should be handled"]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn kepler_universal<T: CordicNumber>(r0: T, sigma0: T, alpha: T, tau: T) -> Result<T> {
    let (r, s, a, t) = (work(r0), work(sigma0), work(alpha), work(tau));
    let two = W::<T>::from_int(2);
    let h2 = W::<T>::dot2(r, two - a * r, -s, s);
    let slack = (r.ldexp(1) + s.square()).ldexp(8 - W::<T>::PREC.cast_signed());
    if !r.is_neg() && !r.is_zero() && h2 >= -slack {
        return Ok(fixed(kepler_work(r, s, a, t)));
    }
    Err(Error::domain(
        "kepler_universal",
        "r0 > 0 and r0 (2 - alpha r0) >= sigma0^2",
    ))
}

fn kepler_work<V: Work>(r0: V, sigma0: V, alpha: V, tau: V) -> V {
    if tau.is_zero() {
        return V::zero();
    }
    let residual = |chi: V| kepler_residual(chi, r0, sigma0, alpha, tau);
    let neg = tau.is_neg();
    let past = |m: V| {
        let (f, _) = residual(m.with_sign(neg));
        f.is_zero() || f.is_neg() == neg
    };
    let start = (tau / r0).abs();
    let (mut near, mut far) = (start, start);
    if past(start) {
        for _ in 0..1200 {
            near = far.ldexp(-1);
            if !past(near) {
                break;
            }
            far = near;
        }
    } else {
        for _ in 0..1200 {
            far = near.ldexp(1);
            if past(far) {
                break;
            }
            near = far;
        }
    }
    let (mut lo, mut hi) = if neg { (-far, -near) } else { (near, far) };
    let mut chi = (lo + hi).ldexp(-1);
    let mut stride = hi - lo;
    for _ in 0..256 {
        let (f, slope) = residual(chi);
        if f.is_zero() {
            break;
        }
        if f.is_neg() {
            lo = chi;
        } else {
            hi = chi;
        }
        // Bisection where Newton leaves the bracket or steps over half its last (rtsafe).
        let step = f / slope;
        let newton = chi - step;
        let slow = step.abs().ldexp(1) > stride.abs();
        let next = if slow || slope.is_neg() || slope.is_zero() || newton <= lo || newton >= hi {
            (lo + hi).ldexp(-1)
        } else {
            newton
        };
        stride = next - chi;
        let settled = V::negligible(&(next - chi), &chi.abs().max(lo.abs()).max(hi.abs()));
        chi = next;
        if settled || next == lo || next == hi {
            break;
        }
    }
    chi
}
