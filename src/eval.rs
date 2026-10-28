//! Function cores in the working arithmetic, shared by the public functions.
//! Each takes and returns working values; the public functions convert in
//! exactly and round once on the way out.

use crate::consts::{
    ATAN_SIXTEENTHS, COS_PI_128THS, EXP2_32NDS, FRAC_1_LN_2, LN_2, LN_SIXTEENTHS, PI, SIN_PI_128THS,
};
use crate::real::Work;
use crate::series::{FACT, ODD};
use crate::traits::CordicNumber;
use crate::traits::sealed::Sealed;

/// The working type of `T`.
pub type W<T> = <T as Sealed>::Work;

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

/// `e^x` for `|x| < HUGE`, as `2^(k/32)·e^r` with `|r| ≤ ln2/64`: the
/// power from the table and a shift, the rest from [`expm1_series`].
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

/// `ln x` for `x` in `[√½, √2)`: `ln(1 + j/16)` from the table, plus
/// `2·atanh((x − c)/(x + c))` with `c = 1 + j/16` nearest `x`.
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

/// `atan z` for `z` in `[0, 1]`: `atan(j/16)` from the table, plus the series
/// for `atan((z − c)/(1 + z·c))` with `c = j/16` nearest `z`.
pub fn atan_unit<V: Work>(z: V) -> V {
    let j = z.ldexp(4).round_int().clamp(0, 16);
    let c = V::from_int(j).ldexp(-4);
    let t = (z - c) / (V::one() + z * c);
    let series = V::horner(t.square(), &ODD, 0, 1, by_width::<V>([8, 14, 21]), true);
    entry::<V>(&ATAN_SIXTEENTHS, j) + t * series
}

/// `atan(num/den)` for `0 ≤ num` and `0 < den`, without forming a ratio
/// above 1.
pub fn atan_ratio<V: Work>(num: V, den: V) -> V {
    if num <= den {
        atan_unit(num / den)
    } else {
        half_pi::<V>() - atan_unit(den / num)
    }
}

/// Sine and cosine of `(π/4)·(oct + w)` for `w` from [`Work::octant`]: the
/// table's angle `jπ/128` nearest `(π/4)·w`, then the two series for the
/// remainder `|d| ≤ π/256`.
pub fn sin_cos_octant<V: Work>(oct: u32, frac: V) -> (V, V) {
    let index = frac.ldexp(5).round_int().clamp(0, 32);
    let rest = (frac - V::from_int(index).ldexp(-5)) * pi::<V>().ldexp(-2);
    let rest2 = rest.square();
    let terms = by_width::<V>([4, 7, 11]);
    let sd = rest * V::horner(rest2, &FACT, 1, 2, terms, true);
    let cd = V::horner(rest2, &FACT, 0, 2, terms, true);
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

/// A value past every type's range, which any of them saturates from.
pub fn beyond<V: Work>() -> V {
    V::one().ldexp(256)
}

/// `φₖ(z) = (e^z − Σ_{n<k} zⁿ/n!)/z^k` for `k` in 1 to 3, which is `1/k!` at
/// 0: by its series near 0, so it keeps relative precision there.
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

/// `ln(1 + v) − v` for `v > −1`. Near 0 it is `2s³·Σ s^{2k}/(2k + 3) −
/// v²/(2 + v)` with `s = v/(2 + v)`, two terms of order `v²` that do not
/// cancel.
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

/// The Stumpff functions `cₖ(z) = Σ (−z)ⁿ/(2n + k)!` for k from 0 to 3, of
/// either sign of `z`: the series near 0, and past it the closed forms in
/// `√|z|`, with `c₂` from the half angle so that it does not cancel.
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
        let (sh, ch) = sinh_cosh(w);
        let (half, _) = sinh_cosh(w.ldexp(-1));
        let c1 = sh / w;
        [
            ch,
            c1,
            (half.square() / w.square()).ldexp(1),
            (c1 - one) / -z,
        ]
    } else {
        let (s, c) = sin_cos(w);
        let (half, _) = sin_cos(w.ldexp(-1));
        let c1 = s / w;
        [c, c1, (half.square() / z).ldexp(1), (one - c1) / z]
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

/// `frac(x/2π)`, the angle as a fraction of a turn in `[0, 1)`, from the
/// same exact reduction as `sin` and `cos`.
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
