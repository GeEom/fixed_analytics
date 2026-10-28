//! Kernels for 128-bit types in Q126, good to 2^−120: 24 guard bits at 96 fractional.

use crate::consts::{
    ATAN_SIXTEENTHS, C, COS_PI_128THS, EXP2_32NDS, FRAC_1_LN_2, LN_2, LN_SIXTEENTHS, PI,
    SIN_PI_128THS,
};
use crate::real::{Real, mul_q126, recip_seed};
use crate::series::{FACT, ODD, q126};
use crate::traits::CordicNumber;

const fn table(src: &[C; 33]) -> [i128; 33] {
    let mut out = [0i128; 33];
    let mut rest: &mut [i128] = &mut out;
    let mut from: &[C] = src;
    while let ([slot, tail @ ..], [c, more @ ..]) = (rest, from) {
        *slot = q126(c);
        rest = tail;
        from = more;
    }
    out
}

const fn table_of<const K: usize>(src: &[C; K]) -> [i128; K] {
    let mut out = [0i128; K];
    let mut rest: &mut [i128] = &mut out;
    let mut from: &[C] = src;
    while let ([slot, tail @ ..], [c, more @ ..]) = (rest, from) {
        *slot = q126(c);
        rest = tail;
        from = more;
    }
    out
}

const fn top(c: &C) -> u128 {
    let (_, _, m) = *c;
    (m[5] as u128) << 64 | m[4] as u128
}

const SIN: [i128; 33] = table(&SIN_PI_128THS);
const EXP2: [i128; 32] = table_of(&EXP2_32NDS);
const ATAN16: [i128; 17] = table_of(&ATAN_SIXTEENTHS);
const LN16: [i128; 13] = table_of(&LN_SIXTEENTHS);
/// `ln2/32·2^132`, rounded; `ln 2·2^128` is its top 128 bits.
const LN2_32NDS: i128 = ((top(&LN_2) >> 1) + (top(&LN_2) & 1)).cast_signed();
/// `32/ln 2·2^40`, from `2^127/ln 2`.
const PER_LN2_32NDS: i128 = (top(&FRAC_1_LN_2) >> 82).cast_signed();
const COS: [i128; 33] = table(&COS_PI_128THS);
/// `π/4·2^126`, from π's top 128 bits, `π·2^126`.
const QUARTER_PI: i128 = {
    let (_, _, m) = PI;
    (((m[5] as u128) << 64 | m[4] as u128) >> 2).cast_signed()
};

/// Whether `T` takes these kernels: 128 bits, at most 96 of them fractional.
pub fn applies<T: CordicNumber>() -> bool {
    T::total_bits() == 128 && T::frac_bits() <= 96
}

/// `|v|·2^(frac − 126)` rounded to nearest even, with `v`'s sign.
pub fn round<T: CordicNumber>(v: i128) -> T {
    let shift = 126 - T::frac_bits();
    let m = v.unsigned_abs();
    let (q, rem, half) = (m >> shift, m & ((1 << shift) - 1), 1u128 << (shift - 1));
    let q = q + u128::from(rem > half || (rem == half && q & 1 == 1));
    let raw = q.cast_signed();
    T::from_raw_i128(if v < 0 { -raw } else { raw })
}

/// `sin` and `cos` of `mag·2^−point` radians, or turns, in Q126.
pub fn sin_cos(mag: u128, point: u32, radians: bool) -> (i128, i128) {
    let (oct, w) = Real::<2>::octant_q127(mag, point, radians);
    let index = ((w + (1 << 121)) >> 122).min(32);
    // (w − index/32)·π/4, from Q127 to Q126; |d| ≤ π/256.
    let rest = w.wrapping_sub(index << 122).cast_signed();
    let d = mul_q126(rest, QUARTER_PI) >> 1;
    let d2 = mul_q126(d, d);
    let coef = |i: u32| FACT.q126.get(i as usize).copied().unwrap_or(0);
    let (mut ps, mut pc) = (0i128, 0i128);
    for k in (0..=7u32).rev() {
        let sign = if k % 2 == 1 { -1 } else { 1 };
        ps = sign * coef(2 * k + 1) + mul_q126(d2, ps);
        pc = sign * coef(2 * k) + mul_q126(d2, pc);
    }
    let sd = mul_q126(d, ps);
    let at = |t: &[i128; 33]| {
        t.get(usize::try_from(index).unwrap_or(0))
            .copied()
            .unwrap_or(0)
    };
    let (sj, cj) = (at(&SIN), at(&COS));
    let s = mul_q126(sj, pc) + mul_q126(cj, sd);
    let c = mul_q126(cj, pc) - mul_q126(sj, sd);
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

/// `e^x = v·2^(q − 126)`, `v` in `[2^126, 2^127)`, for `x = raw·2^−frac` with
/// `|x| ≤ 64` and `8 ≤ frac ≤ 96`: `x = k·ln2/32 + r` exactly, as `r` fits
/// 128 bits whatever the products wrap to, then `2^(j/32)·(1 + expm1 r)`.
pub fn exp(raw: i128, frac: u32) -> Option<(i128, i32)> {
    if raw.unsigned_abs() > 64u128 << frac || !(8..=96).contains(&frac) {
        return None;
    }
    let estimate = (raw << (120 - frac)) >> 80;
    let k = (estimate * PER_LN2_32NDS + (1 << 79)) >> 80;
    let r = raw
        .wrapping_shl(132 - frac)
        .wrapping_sub(k.wrapping_mul(LN2_32NDS))
        >> 6;
    let coef = |i: usize| FACT.q126.get(i).copied().unwrap_or(0);
    let mut p = 0i128;
    for n in (0..=14).rev() {
        p = coef(n + 1) + mul_q126(r, p);
    }
    let base = EXP2
        .get(usize::try_from(k.rem_euclid(32)).unwrap_or(0))
        .copied()
        .unwrap_or(0);
    let value = base + mul_q126(base, mul_q126(r, p));
    Some((value, i32::try_from(k.div_euclid(32)).unwrap_or(0)))
}

/// `v·2^−s` rounded to nearest even, `v ≥ 0`.
pub fn shift_round(v: i128, s: u32) -> i128 {
    if s >= 128 {
        return 0;
    }
    let m = v.unsigned_abs();
    let (q, rem, half) = (m >> s, m & ((1 << s) - 1), 1u128 << (s - 1));
    (q + u128::from(rem > half || (rem == half && q & 1 == 1))).cast_signed()
}

/// `2^252/d` in Q126 for `d` in `[2^125, 2^127)`: a seed, then two steps of `g + g(1 − δg)`.
fn recip(d: i128) -> i128 {
    let z = d.unsigned_abs().leading_zeros();
    let seed = recip_seed(u64::try_from((d.unsigned_abs() << z) >> 64).unwrap_or(u64::MAX));
    let mut g = (seed << (z + 60)).cast_signed();
    for _ in 0..2 {
        g += mul_q126(g, (1i128 << 126) - mul_q126(d, g));
    }
    g
}

fn series(u: i128, table: &[i128], n: usize, alternate: bool) -> i128 {
    let mut p = 0i128;
    for k in (0..=n).rev() {
        let c = table.get(k).copied().unwrap_or(0);
        p = if alternate && k % 2 == 1 { -c } else { c } + mul_q126(u, p);
    }
    p
}

/// `atan2(y, x)` in Q125, neither 0: `atan(j/16)` plus that of `(16a − j·b)/(16b + j·a)`,
/// `a ≤ b` the sizes cut to 122 bits.
pub fn atan2(y: i128, x: i128) -> i128 {
    let (ay, ax) = (y.unsigned_abs(), x.unsigned_abs());
    let (small, large) = if ay <= ax { (ay, ax) } else { (ax, ay) };
    let lz = large.leading_zeros();
    let (small, large) = if lz >= 6 {
        (small << (lz - 6), large << (lz - 6))
    } else {
        (small >> (6 - lz), large >> (6 - lz))
    };
    let j = ((((small >> 57) << 5) / (large >> 57).max(1)) + 1) >> 1;
    let j = j.min(16);
    let (small, large, ji) = (small.cast_signed(), large.cast_signed(), j.cast_signed());
    let ratio = mul_q126((small << 4) - ji * large, recip((large << 4) + ji * small));
    let rest = mul_q126(ratio, series(mul_q126(ratio, ratio), ODD.q126, 13, true));
    let base = (ATAN16
        .get(usize::try_from(j).unwrap_or(0))
        .copied()
        .unwrap_or(0)
        + rest)
        >> 1;
    let pi = (top(&PI) >> 1).cast_signed();
    let angle = if ay > ax { (pi >> 1) - base } else { base };
    let angle = if x < 0 { pi - angle } else { angle };
    if y < 0 { -angle } else { angle }
}

/// `ln x` in Q119, `x = Y·2^k`, `Y` in `[√½, √2)`: `ln c + 2·atanh((Y − c)/(Y + c)) + k·ln 2`.
pub fn ln(raw: i128, frac: u32) -> i128 {
    let lz = raw.unsigned_abs().leading_zeros();
    let mut mant = ((raw.unsigned_abs() << lz) >> 1).cast_signed();
    let mut k = 127 - i64::from(lz) - i64::from(frac);
    if mant > 0x5A82_7999_FCEF_3242_2CBE_C4E8_DB3B_9C48 {
        mant >>= 1;
        k += 1;
    }
    let one = 1i128 << 126;
    let j = ((mant - one + (1 << 121)) >> 122).clamp(-5, 7);
    let centre = one + (j << 122);
    // Halved, so mant + centre stays below 2 in Q126.
    let ratio = mul_q126((mant - centre) >> 1, recip((mant >> 1) + (centre >> 1)));
    let tail = mul_q126(ratio, series(mul_q126(ratio, ratio), ODD.q126, 11, false)) << 1;
    let head = LN16
        .get(usize::try_from(j + 5).unwrap_or(0))
        .copied()
        .unwrap_or(0);
    let ln2 = (top(&LN_2) >> 2).cast_signed();
    let (hi, lo) = (ln2 >> 64, ln2 & ((1 << 64) - 1));
    let k = i128::from(k);
    ((head + tail) >> 7) + ((k * hi) << 57) + ((k * lo) >> 7)
}
