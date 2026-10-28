//! Arithmetic that rounds once.
//!
//! Products, quotients and sums of products are formed exactly in up to 384
//! bits, then rounded to nearest, ties to even. Results beyond the type's
//! range saturate to its nearer end.

use crate::traits::CordicNumber;

const LOW: u128 = u64::MAX as u128;

/// `a · b` as `(high, low)` halves of 256 bits.
const fn wide(a: u128, b: u128) -> (u128, u128) {
    let (ah, al) = (a >> 64, a & LOW);
    let (bh, bl) = (b >> 64, b & LOW);
    let (mid, carry) = (al * bh).overflowing_add(ah * bl);
    let (lo, c2) = (al * bl).overflowing_add(mid << 64);
    let hi = ah * bh + (mid >> 64) + ((carry as u128) << 64) + c2 as u128;
    (hi, lo)
}

/// `(hi·2^128 + lo) / d` and its remainder, for `hi < d`, by two 64-bit
/// digits (Hacker's Delight, divlu).
const fn div_wide(hi: u128, lo: u128, d: u128) -> (u128, u128) {
    let s = d.leading_zeros();
    let d = d << s;
    let (hi, lo) = if s == 0 {
        (hi, lo)
    } else {
        ((hi << s) | (lo >> (128 - s)), lo << s)
    };
    let (dh, dl) = (d >> 64, d & LOW);
    let (lh, ll) = (lo >> 64, lo & LOW);
    let (mut q1, mut r) = div_rem(hi, dh);
    while q1 >> 64 != 0 || q1 * dl > ((r << 64) | lh) {
        q1 -= 1;
        r += dh;
        if r >> 64 != 0 {
            break;
        }
    }
    let mid = ((hi << 64) | lh).wrapping_sub(q1.wrapping_mul(d));
    let (mut q0, mut r0) = div_rem(mid, dh);
    while q0 >> 64 != 0 || q0 * dl > ((r0 << 64) | ll) {
        q0 -= 1;
        r0 += dh;
        if r0 >> 64 != 0 {
            break;
        }
    }
    let rem = ((mid << 64) | ll).wrapping_sub(q0.wrapping_mul(d));
    ((q1 << 64) | q0, rem >> s)
}

/// `a / b` and `a % b`, or zeros for a zero `b`, which callers rule out.
const fn div_rem(a: u128, b: u128) -> (u128, u128) {
    match (a.checked_div(b), a.checked_rem(b)) {
        (Some(q), Some(r)) => (q, r),
        _ => (0, 0),
    }
}

/// The raw value of `T` nearest `±mag`, saturating.
fn saturate<T: CordicNumber>(neg: bool, mag: u128) -> T {
    let raw = if neg {
        0i128.checked_sub_unsigned(mag).unwrap_or(i128::MIN)
    } else {
        i128::try_from(mag).unwrap_or(i128::MAX)
    };
    T::from_raw_i128(raw)
}

/// `(hi, lo) / 2^s` rounded to nearest even, saturating past 128 bits, for
/// `s < 128`.
const fn shift_round(hi: u128, lo: u128, s: u32) -> u128 {
    if s == 0 {
        return if hi == 0 { lo } else { u128::MAX };
    }
    if hi >> s != 0 {
        return u128::MAX;
    }
    let q = (hi << (128 - s)) | (lo >> s);
    let half = (lo >> (s - 1)) & 1 == 1;
    let rest = s > 1 && lo << (129 - s) != 0;
    if half && (rest || q & 1 == 1) {
        q.saturating_add(1)
    } else {
        q
    }
}

/// `a · b + c`, rounded once.
#[must_use]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn mul_add<T: CordicNumber>(a: T, b: T, c: T) -> T {
    let (ra, rb, rc) = (a.raw_i128(), b.raw_i128(), c.raw_i128());
    let (p_hi, p_lo) = wide(ra.unsigned_abs(), rb.unsigned_abs());
    let p_neg = (ra < 0) != (rb < 0);
    let f = T::frac_bits();
    let cm = rc.unsigned_abs();
    let (c_hi, c_lo) = if f == 0 {
        (0, cm)
    } else {
        (cm >> (128 - f), cm << f)
    };
    let (neg, hi, lo) = signed_add(p_neg, p_hi, p_lo, rc < 0, c_hi, c_lo);
    saturate(neg, shift_round(hi, lo, f))
}

/// `±(ah, al) ± (bh, bl)` as a sign and magnitude.
const fn signed_add(
    an: bool,
    ah: u128,
    al: u128,
    bn: bool,
    bh: u128,
    bl: u128,
) -> (bool, u128, u128) {
    if an == bn {
        let (lo, c) = al.overflowing_add(bl);
        return (an, ah.saturating_add(bh).saturating_add(c as u128), lo);
    }
    if ah > bh || (ah == bh && al >= bl) {
        let (lo, b) = al.overflowing_sub(bl);
        (an, ah - bh - b as u128, lo)
    } else {
        let (lo, b) = bl.overflowing_sub(al);
        (bn, bh - ah - b as u128, lo)
    }
}

/// `a · b / c` from the exact 256-bit product, rounded once. A zero `c`
/// saturates by the sign of `a · b`, and gives 0 when that product is 0.
#[must_use]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn mul_div<T: CordicNumber>(a: T, b: T, c: T) -> T {
    let (ra, rb, rc) = (a.raw_i128(), b.raw_i128(), c.raw_i128());
    let neg = ((ra < 0) != (rb < 0)) != (rc < 0);
    let (hi, lo) = wide(ra.unsigned_abs(), rb.unsigned_abs());
    let divisor = rc.unsigned_abs();
    if hi == 0 && lo == 0 {
        return T::zero();
    }
    if divisor == 0 || hi >= divisor {
        return saturate(neg, u128::MAX);
    }
    let (quotient, rem) = div_wide(hi, lo, divisor);
    let twice = rem << 1;
    let up = twice > divisor || (twice == divisor && quotient & 1 == 1);
    saturate(neg, quotient.saturating_add(u128::from(up)))
}

/// `a / b` rounded to nearest. A zero `b` saturates by the sign of `a`, and
/// 0/0 gives 0; it never panics, unlike `fixed`'s division.
#[must_use]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn div_or_saturate<T: CordicNumber>(a: T, b: T) -> T {
    mul_div(a, T::one(), b)
}

/// `x += (hi, lo)` for a 384-bit `x`, most significant word first.
const fn accumulate(x: &mut (u128, u128, u128), hi: u128, lo: u128) {
    let (l, c0) = x.2.overflowing_add(lo);
    let (m1, c1) = x.1.overflowing_add(hi);
    let (m, c2) = m1.overflowing_add(c0 as u128);
    x.2 = l;
    x.1 = m;
    x.0 = x.0.wrapping_add(c1 as u128 + c2 as u128);
}

/// `Σ aᵢ·bᵢ` over the pairs of `a` and `b`, up to the shorter, summed exactly
/// and rounded once.
#[must_use]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn dot<T: CordicNumber>(a: &[T], b: &[T]) -> T {
    let mut pos = (0, 0, 0);
    let mut neg = (0, 0, 0);
    for (&x, &y) in a.iter().zip(b) {
        let (rx, ry) = (x.raw_i128(), y.raw_i128());
        let (hi, lo) = wide(rx.unsigned_abs(), ry.unsigned_abs());
        accumulate(
            if (rx < 0) == (ry < 0) {
                &mut pos
            } else {
                &mut neg
            },
            hi,
            lo,
        );
    }
    let negative = neg > pos;
    let (big, small) = if negative { (neg, pos) } else { (pos, neg) };
    let (lo, b0) = big.2.overflowing_sub(small.2);
    let (m1, b1) = big.1.overflowing_sub(small.1);
    let (mid, b2) = m1.overflowing_sub(u128::from(b0));
    let top = big
        .0
        .wrapping_sub(small.0)
        .wrapping_sub(u128::from(b1 || b2));
    if top != 0 {
        return saturate(negative, u128::MAX);
    }
    saturate(negative, shift_round(mid, lo, T::frac_bits()))
}

/// Square root, 0 for `x ≤ 0`, rounded to nearest like [`sqrt`](crate::sqrt).
#[must_use]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn sqrt_or_zero<T: CordicNumber>(x: T) -> T {
    if x > T::zero() {
        x.sqrt_round()
    } else {
        T::zero()
    }
}
