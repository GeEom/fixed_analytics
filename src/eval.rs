//! Function cores in the working arithmetic, shared by the public functions,
//! which convert in exactly and round once on the way out.

use crate::consts::{
    ATAN_SIXTEENTHS, COS_PI_128THS, EXP2_32NDS, FRAC_1_LN_2, INV_FACT, LN_2, LN_SIXTEENTHS, PI,
    SIN_PI_128THS,
};
use crate::real::Work;
use crate::series::{FACT, ODD};
use crate::traits::CordicNumber;
use crate::traits::sealed::Sealed;

/// The working type of `T`.
pub type W<T> = <T as Sealed>::Work;

/// The narrow working type of `T`, as wide as `T`.
pub type Narrow<T> = <T as Sealed>::Narrow;

/// Bits a narrow result must keep past `T`'s ulp to be trusted.
const GUARD: i32 = 24;

pub fn narrow<T: CordicNumber>(x: T) -> Narrow<T> {
    Narrow::<T>::from_raw(x.raw_i128(), T::frac_bits())
}

/// Whether `T` has a working type narrower than `W<T>` that leaves guard bits,
/// and `extra` more for what the function itself may cancel.
pub fn has_narrow<T: CordicNumber>(extra: i32) -> bool {
    Narrow::<T>::PREC < W::<T>::PREC
        && T::frac_bits().cast_signed() + GUARD + extra < Narrow::<T>::PREC.cast_signed()
}

/// Results in either working type, rounded to `T` all at once.
pub trait Outputs<T> {
    type Rounded;
    fn fits(&self) -> bool;
    fn round(self) -> Self::Rounded;
}

impl<T: CordicNumber, V: Work> Outputs<T> for V {
    type Rounded = T;
    fn fits(&self) -> bool {
        let room = V::PREC.cast_signed() - T::frac_bits().cast_signed() - GUARD;
        self.is_zero() || self.exponent() <= room
    }
    fn round(self) -> T {
        let lo = T::min_value().raw_i128();
        let hi = T::max_value().raw_i128();
        T::from_raw_i128(self.to_raw(T::frac_bits(), lo, hi, false))
    }
}

impl<T: CordicNumber, A: Outputs<T>, B: Outputs<T>> Outputs<T> for (A, B) {
    type Rounded = (A::Rounded, B::Rounded);
    fn fits(&self) -> bool {
        self.0.fits() && self.1.fits()
    }
    fn round(self) -> Self::Rounded {
        (self.0.round(), self.1.round())
    }
}

impl<T: CordicNumber, A: Outputs<T>, const K: usize> Outputs<T> for [A; K] {
    type Rounded = [A::Rounded; K];
    fn fits(&self) -> bool {
        self.iter().all(Outputs::fits)
    }
    fn round(self) -> Self::Rounded {
        self.map(Outputs::round)
    }
}

/// `$body`, arguments converted by `$v`, in `Narrow<T>` where every result fits, else in `W<T>`.
macro_rules! tiered {
    ($t:ty, |$v:ident| $body:expr) => {
        $crate::eval::tiered!($t, 0, |$v| $body)
    };
    ($t:ty, $extra:expr, |$v:ident| $body:expr) => {{
        let mut done = None;
        if $crate::eval::has_narrow::<$t>($extra) {
            let $v = $crate::eval::narrow::<$t>;
            let out = $body;
            if $crate::eval::Outputs::<$t>::fits(&out) {
                done = Some($crate::eval::Outputs::<$t>::round(out));
            }
        }
        match done {
            Some(r) => r,
            None => {
                let $v = $crate::eval::work::<$t>;
                $crate::eval::Outputs::<$t>::round($body)
            }
        }
    }};
}
pub(crate) use tiered;

/// Arguments past this saturate every function of them that grows.
pub const HUGE: i64 = 4096;

pub fn work<T: CordicNumber>(x: T) -> W<T> {
    W::<T>::from_raw(x.raw_i128(), T::frac_bits())
}

/// Rounds to the nearest value of `T`, saturating.
pub fn fixed<T: CordicNumber>(v: W<T>) -> T {
    let lo = T::min_value().raw_i128();
    let hi = T::max_value().raw_i128();
    T::from_raw_i128(v.to_raw(T::frac_bits(), lo, hi, false))
}

pub fn pi<V: Work>() -> V {
    V::constant(&PI)
}

pub fn half_pi<V: Work>() -> V {
    pi::<V>().ldexp(-1)
}

/// The entry for this width, from entries for 64-, 128- and 192-bit mantissas.
pub const fn by_width<V: Work>(n: [u32; 3]) -> u32 {
    match V::PREC {
        0..=64 => n[0],
        65..=128 => n[1],
        _ => n[2],
    }
}

fn entry<V: Work>(table: &[crate::consts::C], i: i64) -> V {
    usize::try_from(i)
        .ok()
        .and_then(|i| table.get(i))
        .map_or_else(V::zero, |c| V::constant(c))
}

/// `expm1(r)` for `|r| ≤ ln2/64`, as `r·Σ r^n/(n+1)!`.
fn expm1_series<V: Work>(r: V) -> V {
    r * V::horner(r, &FACT, 1, 1, by_width::<V>([8, 14, 20]), false)
}

/// `e^x` for `|x| < HUGE`, as `2^(k/32)·e^r`, `|r| ≤ ln2/64`, the power from the table.
pub fn exp<V: Work>(x: V) -> V {
    let steps = (x * V::constant(&FRAC_1_LN_2).ldexp(5)).round_int();
    let rest = x - V::constant(&LN_2).ldexp(-5).mul_int(steps);
    let index = steps.rem_euclid(32);
    let shift = i32::try_from(steps.div_euclid(32)).unwrap_or(0);
    let power = entry::<V>(&EXP2_32NDS, index);
    (power + power * expm1_series(rest)).ldexp(shift)
}

/// `e^x − 1` for `|x| < HUGE`, to full relative precision near 0.
pub fn expm1<V: Work>(x: V) -> V {
    if x.abs() <= V::constant(&LN_2).ldexp(-6) {
        expm1_series(x)
    } else {
        exp(x) - V::one()
    }
}

/// `atanh(s) = s + s³/3 + s⁵/5 + …` for `|s| ≤ 0.023`.
fn atanh_series<V: Work>(s: V) -> V {
    s * V::horner(s.square(), &ODD, 0, 1, by_width::<V>([8, 13, 19]), false)
}

/// `ln x` on `[√½, √2)`: `ln c`, `c = 1 + j/16` nearest, plus `2·atanh((x − c)/(x + c))`.
fn ln_near_one<V: Work>(y: V) -> V {
    let one = V::one();
    let j = (y - one).ldexp(4).round_int().clamp(-5, 7);
    let c = one + V::from_int(j).ldexp(-4);
    entry::<V>(&LN_SIXTEENTHS, j + 5) + atanh_series((y - c) / (y + c)).ldexp(1)
}

/// `x = y·2^k` with `y` in `[√½, √2)`.
fn split_log<V: Work>(x: V) -> (V, i32) {
    let e = x.exponent();
    let y = x.ldexp(-e);
    if y.square() < V::one().ldexp(-1) {
        (y.ldexp(1), e - 1)
    } else {
        (y, e)
    }
}

/// Natural logarithm of `x > 0`.
pub fn ln<V: Work>(x: V) -> V {
    let (y, k) = split_log(x);
    ln_near_one(y) + V::constant(&LN_2).mul_int(i64::from(k))
}

/// Base-2 logarithm of `x > 0`, its integer part exact.
pub fn log2<V: Work>(x: V) -> V {
    let (y, k) = split_log(x);
    ln_near_one(y) * V::constant(&FRAC_1_LN_2) + V::from_int(i64::from(k))
}

/// `ln(1 + u)` for `u > −1`, to full relative precision near 0.
pub fn log1p<V: Work>(u: V) -> V {
    let one = V::one();
    if u.abs() <= one.ldexp(-5) {
        atanh_series(u / (u + V::from_int(2))).ldexp(1)
    } else {
        ln(one + u)
    }
}

/// `atan(num/den)` for `0 ≤ num ≤ den`: `atan c`, `c = j/16` nearest, plus that of
/// `(16num − j·den)/(16den + j·num)`, one division in all.
fn atan_unit<V: Work>(num: V, den: V) -> V {
    let j = V::sixteenths(num, den).clamp(0, 16);
    let t = (num.ldexp(4) - den.mul_int(j)) / (den.ldexp(4) + num.mul_int(j));
    let series = V::horner(t.square(), &ODD, 0, 1, by_width::<V>([8, 14, 21]), true);
    entry::<V>(&ATAN_SIXTEENTHS, j) + t * series
}

/// `atan(num/den)` for `0 ≤ num` and `0 < den`, never forming a ratio above 1.
pub fn atan_ratio<V: Work>(num: V, den: V) -> V {
    if num <= den {
        atan_unit(num, den)
    } else {
        half_pi::<V>() - atan_unit(den, num)
    }
}

/// Sine and cosine of `(π/4)·(oct + w)`, `w` from [`Work::octant`]: the table's
/// `jπ/128` nearest, then series for the rest, `|d| ≤ π/256`.
pub fn sin_cos_octant<V: Work>(oct: u32, frac: V) -> (V, V) {
    let index = frac.ldexp(5).round_int().clamp(0, 32);
    let rest = (frac - V::from_int(index).ldexp(-5)) * pi::<V>().ldexp(-2);
    let rest2 = rest.square();
    let terms = by_width::<V>([4, 7, 11]);
    let (sd, cd) = V::horner2(rest2, &FACT, [1, 0], 2, terms, true);
    let sd = rest * sd;
    let (sj, cj) = (
        entry::<V>(&SIN_PI_128THS, index),
        entry::<V>(&COS_PI_128THS, index),
    );
    let s = sj * cd + cj * sd;
    let c = cj * cd - sj * sd;
    match oct & 7 {
        0 => (s, c),
        1 => (c, s),
        2 => (c, -s),
        3 => (s, -c),
        4 => (-s, -c),
        5 => (-c, -s),
        6 => (-c, s),
        _ => (-s, c),
    }
}

/// A value past `e^HUGE`, which every type saturates from even divided by its values.
pub fn beyond<V: Work>() -> V {
    V::one().ldexp(8192)
}

/// `e^v`, past [`HUGE`] a value every type saturates or rounds to zero from.
pub fn exp_wide<V: Work>(v: V) -> V {
    if v.abs() > V::from_int(HUGE) {
        return if v.is_neg() { V::zero() } else { beyond() };
    }
    exp(v)
}

/// `φₖ(z) = (e^z − Σ_{n<k} zⁿ/n!)/z^k` for `k` in 1 to 3, by its series near 0.
pub fn phi<V: Work>(k: u32, z: V) -> V {
    if z.abs() <= V::one().ldexp(-1) {
        return V::horner(z, &FACT, k, 1, by_width::<V>([18, 30, 39]), false);
    }
    if z > V::from_int(HUGE) {
        return beyond();
    }
    let mut num = if z < V::from_int(-HUGE) {
        -V::one()
    } else {
        expm1(z)
    };
    if k >= 2 {
        num = num - z;
    }
    if k >= 3 {
        num = num - z.square().ldexp(-1);
    }
    let mut den = z;
    for _ in 1..k {
        den = den * z;
    }
    num / den
}

/// `φ₁` to `φ₅` given `e^(−z)`: near 0 `φ₅`'s series and `φₖ = 1/k! + z·φₖ₊₁`, past it
/// `φ₁` from `e^z − 1` and `φₖ₊₁ = (φₖ − 1/k!)/z`, losing at most 11 bits, at `|z| = ½`.
pub fn phis<V: Work>(z: V, decay: V) -> [V; 5] {
    let inv_fact = |k: usize| INV_FACT.get(k).map_or_else(V::zero, V::constant);
    let mut out = [V::zero(); 5];
    if z.abs() <= V::one().ldexp(-1) {
        let mut next = V::horner(z, &FACT, 5, 1, by_width::<V>([18, 30, 39]), false);
        for (k, slot) in out.iter_mut().enumerate().rev() {
            *slot = next;
            next = inv_fact(k) + z * next;
        }
        return out;
    }
    if z > V::from_int(HUGE) {
        return [beyond(); 5];
    }
    // e^z − 1 = 1/e^(−z) − 1, which past |z| = ½ cancels by at most two bits.
    let mut cur = if z < V::from_int(-HUGE) {
        -V::one()
    } else {
        V::one() / decay - V::one()
    } / z;
    for (k, slot) in out.iter_mut().enumerate() {
        *slot = cur;
        cur = (cur - inv_fact(k + 1)) / z;
    }
    out
}

/// `ln(1 + v) − v`, `v > −1`; near 0 `2s³·Σ s^{2k}/(2k + 3) − v²/(2 + v)`, `s = v/(2 + v)`,
/// two terms of order `v²` that do not cancel.
pub fn log1p_minus<V: Work>(v: V) -> V {
    if v.abs() <= V::one().ldexp(-2) {
        let den = v + V::from_int(2);
        let s = v / den;
        let sum = V::horner(s.square(), &ODD, 1, 1, by_width::<V>([13, 24, 36]), false);
        return (s * s.square() * sum).ldexp(1) - v.square() / den;
    }
    log1p(v) - v
}

/// `(ln(1 + u) − u)/u²`, which is −½ at 0.
pub fn log1p_minus_ratio<V: Work>(u: V) -> V {
    if u.abs() <= V::one().ldexp(-2) {
        let den = u + V::from_int(2);
        let s = u / den;
        let sum = V::horner(s.square(), &ODD, 1, 1, by_width::<V>([13, 24, 36]), false);
        return (u * sum).ldexp(1) / (den * den.square()) - V::one() / den;
    }
    log1p_minus(u) / u.square()
}

/// `ln(1 + z)/z`, which is 1 at 0.
pub fn log1p_ratio<V: Work>(z: V) -> V {
    if z.abs() <= V::one().ldexp(-2) {
        let den = z + V::from_int(2);
        let s = z / den;
        let sum = V::horner(s.square(), &ODD, 0, 1, by_width::<V>([13, 24, 36]), false);
        return sum.ldexp(1) / den;
    }
    log1p(z) / z
}

/// The Stumpff functions `cₖ(z) = Σ (−z)ⁿ/(2n + k)!`, k from 0 to 3: the series near
/// 0, past it closed forms in `√|z|`, `c₂` from the half angle so it does not cancel.
pub fn stumpff<V: Work>(z: V) -> [V; 4] {
    if z.abs() <= V::one().ldexp(-2) {
        let n = by_width::<V>([10, 17, 21]);
        return [0, 1, 2, 3].map(|k| V::horner(-z, &FACT, k, 2, n, false));
    }
    let one = V::one();
    let w = z.abs().sqrt();
    if z.is_neg() {
        if w > V::from_int(HUGE) {
            return [beyond(); 4];
        }
        let (half, half_c) = sinh_cosh(w.ldexp(-1));
        let c1 = (half * half_c).ldexp(1) / w;
        [
            one + half.square().ldexp(1),
            c1,
            (half.square() / w.square()).ldexp(1),
            (c1 - one) / -z,
        ]
    } else {
        // sin w = 2·s·c and cos w = (c − s)(c + s), exact in form at its zeros.
        let (half, half_c) = sin_cos(w.ldexp(-1));
        let c1 = (half * half_c).ldexp(1) / w;
        let c0 = (half_c - half) * (half_c + half);
        [c0, c1, (half.square() / z).ldexp(1), (one - c1) / z]
    }
}

/// Sine and cosine of a working value, reduced exactly past π/4.
pub fn sin_cos<V: Work>(w: V) -> (V, V) {
    let a = w.abs();
    let quarter = pi::<V>().ldexp(-2);
    let (oct, frac) = if a <= quarter {
        (0, a / quarter)
    } else {
        a.octant_of()
    };
    let (s, c) = sin_cos_octant(oct, frac);
    (if w.is_neg() { -s } else { s }, c)
}

/// `frac(x/2π)` in `[0, 1)`, from the same exact reduction as `sin` and `cos`.
pub fn turn_fraction<T: CordicNumber>(x: T) -> W<T> {
    let raw = x.raw_i128();
    let (oct, w) = W::<T>::octant(raw.unsigned_abs(), T::frac_bits(), true);
    let within = if oct & 1 == 1 { W::<T>::one() - w } else { w };
    let u = (W::<T>::from_int(i64::from(oct)) + within).ldexp(-3);
    if raw < 0 && !u.is_zero() {
        W::<T>::one() - u
    } else {
        u
    }
}

/// `(sinh x, cosh x)` for `0 ≤ x < HUGE`, from `u = expm1(x)` so neither
/// cancels: `sinh = (u + u/(1 + u))/2` and `cosh = 1 + u²/(2(1 + u))`.
pub fn sinh_cosh<V: Work>(x: V) -> (V, V) {
    let u = expm1(x);
    let one = V::one();
    let e = one + u;
    ((u + u / e).ldexp(-1), one + (u.square() / e).ldexp(-1))
}
