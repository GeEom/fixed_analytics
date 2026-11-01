//! A bracketed root finder, and the universal-variable Kepler equation.

use crate::error::{Error, Result};
use crate::eval::{GUARD, HUGE, W, beyond, ln, stumpff, tiered, work};
use crate::real::{Real, Work};
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
        let raw = x
            .to_raw(0, a, b, true)
            .max(a.saturating_add(1))
            .min(b.saturating_sub(1));
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

/// The bracket of [`root_in`] in at most `budget` more calls, by safeguarded Newton from `guess`.
///
/// `step(x, f(x))` is the distance to the root, such as `|f/f′|` or a Halley step; the probe
/// gives the side. Each step overshoots by a nudge from one ulp, doubled while probes stay on one
/// side, and one that fails to halve the last bisects, as does a `guess` outside `(lo, hi)`.
/// # Errors
/// Returns `DomainError` if `f(lo)` and `f(hi)` are non-zero with one sign.
#[must_use = "returns the bracket which should be handled"]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn root_in_newton<T: CordicNumber, F: FnMut(T) -> T, S: FnMut(T, T) -> T>(
    mut func: F,
    mut step: S,
    guess: T,
    lo: T,
    hi: T,
    budget: u32,
) -> Result<(T, T)> {
    let (lo, hi) = if lo <= hi { (lo, hi) } else { (hi, lo) };
    let (ya, yb) = (func(lo), func(hi));
    let zero = T::zero();
    if ya == zero {
        return Ok((lo, lo));
    }
    if yb == zero {
        return Ok((hi, hi));
    }
    let low = ya.is_negative();
    if yb.is_negative() == low {
        return Err(Error::domain(
            "root_in_newton",
            "f(lo) and f(hi) of opposite signs",
        ));
    }
    let (mut a, mut b, mut at) = (lo.raw_i128(), hi.raw_i128(), guess.raw_i128());
    let (mut last, mut nudge, mut side) = (u128::MAX, 1u128, None);
    for _ in 0..budget {
        let span = b.abs_diff(a);
        if span <= 1 {
            break;
        }
        if at <= a || at >= b {
            at = a.saturating_add_unsigned(span / 2);
        }
        let x = T::from_raw_i128(at);
        let y = func(x);
        if y == zero {
            return Ok((x, x));
        }
        let under = y.is_negative() == low;
        if under {
            a = at;
        } else {
            b = at;
        }
        if b.abs_diff(a) <= 1 {
            break;
        }
        nudge = if side == Some(under) {
            nudge.saturating_mul(2)
        } else {
            1
        };
        side = Some(under);
        let size = step(x, y).raw_i128().unsigned_abs();
        let reach = size.saturating_add(nudge);
        // Nor does a second zero step on one side.
        at = if size > last / 2 || (size == 0 && last == 0 && nudge > 1) {
            a
        } else if under {
            at.saturating_add_unsigned(reach)
        } else {
            at.saturating_sub_unsigned(reach)
        };
        last = size;
    }
    Ok((T::from_raw_i128(a), T::from_raw_i128(b)))
}

#[derive(Clone, Copy)]
struct Orbit<V> {
    r0: V,
    sigma0: V,
    alpha: V,
    tau: V,
    beta: V,
}

impl<V: Work> Orbit<V> {
    fn new(r0: V, sigma0: V, alpha: V, tau: V) -> Self {
        let beta = V::one() - alpha * r0;
        Self {
            r0,
            sigma0,
            alpha,
            tau,
            beta,
        }
    }

    fn to<U: Work>(self) -> Orbit<U> {
        let [r0, sigma0, alpha, tau] = [self.r0, self.sigma0, self.alpha, self.tau].map(convert);
        Orbit::new(r0, sigma0, alpha, tau)
    }
}

fn convert<A: Work, B: Work>(v: A) -> B {
    if v.is_zero() {
        return B::zero();
    }
    let e = v.exponent();
    B::from_raw(v.ldexp(-e).to_raw(126, i128::MIN, i128::MAX, false), 126).ldexp(e)
}

/// `F − τ`, `r = F′`, `σ = r′`, `r″ = 1 − αr` and `χᵏcₖ`; past `HUGE²` a slope of 0.
fn kepler_residual<V: Work>(o: &Orbit<V>, chi: V) -> ([V; 4], [V; 4]) {
    let x2 = chi.square();
    let z = o.alpha * x2;
    // Past HUGE² e^√|z| swamps every term, its coefficient on a real orbit of χ's sign.
    let huge = V::from_int(HUGE);
    if z < -huge.square() {
        let zero = V::zero();
        return (
            [beyond::<V>().with_sign(chi.is_neg()), zero, zero, zero],
            [zero; 4],
        );
    }
    let [c0, c1, c2, c3] = stumpff(z);
    let (g1, g2, g3) = (chi * c1, x2 * c2, x2 * chi * c3);
    let f = o.r0 * chi + o.sigma0 * g2 + o.beta * g3 - o.tau;
    let r = o.r0 * c0 + o.sigma0 * g1 + g2;
    let sigma = o.sigma0 * c0 + o.beta * g1;
    ([f, r, sigma, V::one() - o.alpha * r], [c0, g1, g2, g3])
}

/// The universal anomaly `χ` after time `t` on a two-body orbit.
///
/// It solves `√μ·t = r₀χ + σ₀χ²c₂(αχ²) + (1 − αr₀)χ³c₃(αχ²)`, rising in `χ`, for `r₀ = |r⃗₀|`,
/// `σ₀ = r⃗₀·v⃗₀/√μ`, `α = 2/r₀ − v₀²/μ` and `τ = √μ·t`, by safeguarded Halley steps.
/// # Errors
/// Returns `DomainError` unless `r₀ > 0` and `r₀(2 − αr₀) ≥ σ₀²`, a real orbit.
#[must_use = "returns the anomaly which should be handled"]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn kepler_universal<T: CordicNumber>(r0: T, sigma0: T, alpha: T, tau: T) -> Result<T> {
    if kepler_real(r0, sigma0, alpha) {
        let frac = T::frac_bits();
        return Ok(tiered!(T, |v| kepler_work(
            v(r0),
            v(sigma0),
            v(alpha),
            v(tau),
            v(T::zero()),
            frac
        )
        .0));
    }
    Err(Error::domain(
        "kepler_universal",
        "r0 > 0 and r0 (2 - alpha r0) >= sigma0^2",
    ))
}

fn kepler_real<T: CordicNumber>(r0: T, sigma0: T, alpha: T) -> bool {
    let (r, s, a) = (work(r0), work(sigma0), work(alpha));
    let two = W::<T>::from_int(2);
    let h2 = W::<T>::dot2(r, two - a * r, -s, s);
    let slack = (r.ldexp(1) + s.square()).ldexp(8 - W::<T>::PREC.cast_signed());
    !r.is_neg() && !r.is_zero() && h2 >= -slack
}

/// `χ`, the last point evaluated and its `χᵏcₖ`; a wide `V` starts at a 64-bit root.
fn kepler_work<V: Work>(r0: V, sigma0: V, alpha: V, tau: V, near: V, frac: u32) -> (V, V, [V; 4]) {
    let zero = V::zero();
    if tau.is_zero() {
        return (zero, zero, [V::one(), zero, zero, zero]);
    }
    let orbit = Orbit::new(r0, sigma0, alpha, tau);
    let seeded = !near.is_zero() && near.is_neg() == tau.is_neg();
    let start = if V::PREC <= 64 {
        if seeded { near } else { kepler_start(&orbit) }
    } else {
        let cheap: Orbit<Real<1>> = orbit.to();
        let from = if seeded {
            convert(near)
        } else {
            kepler_start(&cheap)
        };
        convert(kepler_iterate(&cheap, from, i32::MIN).0)
    };
    let (chi, at, powers) = kepler_iterate(&orbit, start, -frac.cast_signed() - GUARD);
    // Within 2^-20 ulp of a tie, the residual at the midpoint picks the side.
    let raw = chi.to_raw(frac, i128::MIN, i128::MAX, false);
    let offset = chi - V::from_raw(raw, frac);
    let half = V::one().ldexp(-frac.cast_signed() - 1);
    if offset.abs() <= half && offset.abs() > half - half.ldexp(-19) {
        let mid = chi - offset + half.with_sign(offset.is_neg());
        let ([res, ..], at_mid) = kepler_residual(&orbit, mid);
        let side = if res.is_zero() {
            zero
        } else {
            half.ldexp(-21).with_sign(!res.is_neg())
        };
        return (mid + side, mid, at_mid);
    }
    (chi, at, powers)
}

/// `|χ|` to start from: a Halley step on the `α = 0` cubic from `|τ|/r₀`, or past `|αχ²| = 1`
/// the long arc's `ατ` or `ln(2|τ|(−α)^{3/2}/(β + σ₀√−α))/√−α`.
fn kepler_start<V: Work>(orbit: &Orbit<V>) -> V {
    let Orbit {
        r0, alpha, beta, ..
    } = *orbit;
    let time = orbit.tau.abs();
    let sigma = if orbit.tau.is_neg() {
        -orbit.sigma0
    } else {
        orbit.sigma0
    };
    let linear = time / r0;
    let slope = r0 + sigma * linear + (beta * linear.square()).ldexp(-1);
    let bend = sigma + beta * linear;
    let excess = linear.square() * (sigma.mul_int(3) + beta * linear).div_int(6);
    let halley = slope.square().ldexp(1) - excess * bend;
    let mut cubic = if slope.is_neg() || slope.is_zero() {
        linear
    } else if halley > slope.square() {
        linear - (excess * slope).ldexp(1) / halley
    } else {
        linear - excess / slope
    };
    if cubic.is_neg() || cubic.is_zero() {
        cubic = linear.ldexp(-1);
    }
    let z = alpha * cubic.square();
    let one = V::one();
    if z > one {
        let mean = alpha * time;
        if alpha * mean.square() > V::from_int(4) {
            return mean;
        }
    } else if z < -one {
        let root = (-alpha).sqrt();
        let lead = beta + sigma * root;
        if !lead.is_neg() && !lead.is_zero() {
            let arg = (time * (-alpha) * root).ldexp(1) / lead;
            if arg > V::from_int(3) {
                return ln(arg) / root;
            }
        }
    }
    cubic
}

/// `χ` by Halley steps in a bracket doubled while open, else split where a step leaves it or
/// fails to halve the last, to an error below `2^tol` or `V`'s noise.
fn kepler_iterate<V: Work>(orbit: &Orbit<V>, start: V, tol: i32) -> (V, V, [V; 4]) {
    let neg = orbit.tau.is_neg();
    let (mut lo, mut hi, mut open) = (V::zero(), V::zero(), true);
    let (mut mag, mut stride) = (start.abs(), beyond::<V>());
    let mut last = (mag, [V::zero(); 4]);
    let quarter = V::one().ldexp(-2);
    for _ in 0..1500 {
        let ([res, slope, sigma, curl], powers) = kepler_residual(orbit, mag.with_sign(neg));
        last = (mag, powers);
        let (res, sigma) = if neg { (-res, -sigma) } else { (res, sigma) };
        if res.is_zero() {
            break;
        }
        if res.is_neg() {
            lo = mag;
        } else {
            (hi, open) = (mag, false);
        }
        let mut next = None;
        if !slope.is_neg() && !slope.is_zero() {
            let newton = res / slope;
            let k1 = sigma / slope;
            let bend = (newton * k1).ldexp(-1);
            let halley = bend < quarter.ldexp(1);
            let step = if halley {
                newton / (V::one() - bend)
            } else {
                newton
            };
            // Newton leaves (σ/2r)δ², Halley ((σ/2r)² − r″/6r)δ³, bounded twice over.
            let size = step.abs();
            let tilt = size * k1.abs();
            let k2 = curl.abs() / slope;
            let order = if halley { tilt.square() } else { tilt };
            let left = size * (order + size.square() * k2);
            let floor = tol.max(mag.exponent() - (V::PREC.cast_signed() - 4));
            if size.is_zero()
                || size.exponent() < floor
                || (tilt < quarter && left.exponent() + 2 < floor)
            {
                mag = mag - step;
                break;
            }
            let to = mag - step;
            let halves = stride.abs() >= size.ldexp(1);
            if halves && to > lo && (open || to < hi) {
                next = Some(to);
            }
        }
        let to = next.unwrap_or_else(|| {
            let gap = hi.exponent() - lo.exponent();
            if open {
                lo.ldexp(1)
            } else if lo.is_zero() {
                hi.ldexp(-1)
            } else if gap > 1 {
                lo.ldexp(gap / 2)
            } else {
                (lo + hi).ldexp(-1)
            }
        });
        stride = to - mag;
        mag = to;
        if mag == lo || mag == hi {
            break;
        }
    }
    let (at, powers) = last;
    (mag.with_sign(neg), at.with_sign(neg), powers)
}

/// [`kepler_universal`]'s `χ` and the Stumpff functions `[c₀, c₁, c₂, c₃]` at `αχ²` for it.
/// # Errors
/// Returns `DomainError` unless `r₀ > 0` and `r₀(2 − αr₀) ≥ σ₀²`, a real orbit.
#[must_use = "returns the anomaly which should be handled"]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn kepler_universal_stumpff<T: CordicNumber>(
    r0: T,
    sigma0: T,
    alpha: T,
    tau: T,
) -> Result<(T, [T; 4])> {
    kepler_with_stumpff(
        "kepler_universal_stumpff",
        r0,
        sigma0,
        alpha,
        tau,
        T::zero(),
    )
}

/// [`kepler_universal_stumpff`] from `near`, an anomaly close by, such as a nearby time's.
///
/// The start changes only the work, not the result; 0 or one of the wrong sign is ignored.
/// # Errors
/// Returns `DomainError` unless `r₀ > 0` and `r₀(2 − αr₀) ≥ σ₀²`, a real orbit.
#[must_use = "returns the anomaly which should be handled"]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn kepler_universal_near<T: CordicNumber>(
    r0: T,
    sigma0: T,
    alpha: T,
    tau: T,
    near: T,
) -> Result<(T, [T; 4])> {
    kepler_with_stumpff("kepler_universal_near", r0, sigma0, alpha, tau, near)
}

fn kepler_with_stumpff<T: CordicNumber>(
    name: &'static str,
    r0: T,
    sigma0: T,
    alpha: T,
    tau: T,
    near: T,
) -> Result<(T, [T; 4])> {
    if kepler_real(r0, sigma0, alpha) {
        let frac = T::frac_bits();
        return Ok(tiered!(T, |v| {
            let (chi, at, g) = kepler_work(v(r0), v(sigma0), v(alpha), v(tau), v(near), frac);
            kepler_stumpff_at::<T, _>(v(alpha), chi, at, g)
        }));
    }
    Err(Error::domain(
        name,
        "r0 > 0 and r0 (2 - alpha r0) >= sigma0^2",
    ))
}

/// `χ` rounded to `T`, and `cₖ(αχ²)` there by Taylor in `gₖ = xᵏcₖ`, `gₖ′ = gₖ₋₁`, else afresh.
fn kepler_stumpff_at<T: CordicNumber, V: Work>(
    alpha: V,
    chi: V,
    at: V,
    powers: [V; 4],
) -> (V, [V; 4]) {
    let frac = T::frac_bits();
    let raw = chi.to_raw(
        frac,
        T::min_value().raw_i128(),
        T::max_value().raw_i128(),
        false,
    );
    let grid = V::from_raw(raw, frac);
    let gap = grid - at;
    let below = |v: V, e: i32| v.is_zero() || v.exponent() < e;
    let third = V::PREC.cast_signed() / 3;
    let close = !grid.is_zero()
        && !at.is_zero()
        && below(gap, at.exponent() - third)
        && below(alpha * gap.square(), -2 * third);
    if !close {
        return (grid, stumpff(alpha * grid.square()));
    }
    let [g0, g1, g2, g3] = powers;
    let half = gap.square().ldexp(-1);
    let (by_gap, by_half) = (alpha * gap, alpha * half);
    let s0 = g0 - by_gap * g1 - by_half * g0;
    let s1 = g1 + gap * g0 - by_half * g1;
    let s2 = g2 + gap * g1 + half * g0;
    let s3 = g3 + gap * g2 + half * g1;
    let inv = V::one() / grid;
    let inv2 = inv.square();
    (grid, [s0, s1 * inv, s2 * inv2, s3 * inv2 * inv])
}
