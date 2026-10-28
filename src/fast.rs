//! Kernels for 128-bit types in Q126, good to 2^−120: 24 guard bits at 96 fractional.

use crate::consts::{C, COS_PI_128THS, EXP2_32NDS, FRAC_1_LN_2, LN_2, PI, SIN_PI_128THS};
use crate::real::{Real, mul_q126};
use crate::series::{FACT, q126};
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

const fn table32(src: &[C; 32]) -> [i128; 32] {
    let mut out = [0i128; 32];
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
const EXP2: [i128; 32] = table32(&EXP2_32NDS);
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
