//! The arithmetic every function evaluates in: binary floating point with an
//! `N`-limb mantissa, rounded to nearest even each step, 1.5 to 8 times a type's width.

use crate::consts::{C, FRAC_1_2PI, FRAC_1_2PI_LONG};
use crate::series::Series;
use core::cmp::Ordering;
use core::fmt::Debug;
use core::ops::{Add, Div, Mul, Neg, Sub};

const SCRATCH: usize = 10;
type Scratch = [u64; SCRATCH];

fn limb(a: &[u64], i: usize) -> u64 {
    a.get(i).copied().unwrap_or(0)
}

fn set(a: &mut [u64], i: usize, v: u64) {
    if let Some(x) = a.get_mut(i) {
        *x = v;
    }
}

fn lo64(x: u128) -> u64 {
    u64::try_from(x & u128::from(u64::MAX)).unwrap_or(0)
}

fn hi64(x: u128) -> u64 {
    u64::try_from(x >> 64).unwrap_or(0)
}

fn u32_of(x: usize) -> u32 {
    u32::try_from(x).unwrap_or(u32::MAX)
}

fn i32_of(x: i64) -> i32 {
    i32::try_from(x.clamp(i64::from(i32::MIN / 4), i64::from(i32::MAX / 4))).unwrap_or(0)
}

fn shl(limbs: &mut [u64], len: usize, shift: u32) {
    let words = (shift / 64) as usize;
    let bits = shift % 64;
    for i in (0..len).rev() {
        let src = if i >= words {
            limb(limbs, i - words)
        } else {
            0
        };
        let below = if i > words {
            limb(limbs, i - words - 1)
        } else {
            0
        };
        let word = if bits == 0 {
            src
        } else {
            (src << bits) | (below >> (64 - bits))
        };
        set(limbs, i, word);
    }
}

/// Shifts right by `s`, returning whether a set bit fell off the end.
fn shr(limbs: &mut [u64], len: usize, shift: u32) -> bool {
    let words = (shift / 64) as usize;
    let bits = shift % 64;
    let mut sticky = (0..words.min(len)).any(|i| limb(limbs, i) != 0);
    if words < len && bits > 0 {
        sticky |= limb(limbs, words) << (64 - bits) != 0;
    }
    for i in 0..len {
        let src = if i + words < len {
            limb(limbs, i + words)
        } else {
            0
        };
        let above = if i + words + 1 < len {
            limb(limbs, i + words + 1)
        } else {
            0
        };
        let word = if bits == 0 {
            src
        } else {
            (src >> bits) | (above << (64 - bits))
        };
        set(limbs, i, word);
    }
    sticky
}

fn leading_zeros(a: &Scratch, len: usize) -> u32 {
    for i in (0..len).rev() {
        let v = limb(a, i);
        if v != 0 {
            return u32_of(len - 1 - i) * 64 + v.leading_zeros();
        }
    }
    u32_of(len) * 64
}

fn add_into(a: &mut Scratch, b: &Scratch, len: usize) -> bool {
    let mut carry = false;
    for i in 0..len {
        let (s1, c1) = limb(a, i).overflowing_add(limb(b, i));
        let (s2, c2) = s1.overflowing_add(u64::from(carry));
        set(a, i, s2);
        carry = c1 || c2;
    }
    carry
}

fn sub_into(a: &mut Scratch, b: &Scratch, len: usize) {
    let mut borrow = false;
    for i in 0..len {
        let (d1, b1) = limb(a, i).overflowing_sub(limb(b, i));
        let (d2, b2) = d1.overflowing_sub(u64::from(borrow));
        set(a, i, d2);
        borrow = b1 || b2;
    }
}

fn bit(a: &Scratch, k: u32) -> bool {
    (limb(a, (k / 64) as usize) >> (k % 64)) & 1 == 1
}

fn any_below(a: &Scratch, k: u32) -> bool {
    let w = (k / 64) as usize;
    let b = k % 64;
    (0..w).any(|i| limb(a, i) != 0) || (b > 0 && limb(a, w) << (64 - b) != 0)
}

/// `|v|·2^62` for a mantissa top limb whose value is `top · 2^(e − 64)`, `|v| < 2`.
fn q62(neg: bool, e: i32, top: u64) -> i64 {
    let sh = 2 - i64::from(e);
    let q = u32::try_from(sh)
        .ok()
        .and_then(|s| top.checked_shr(s))
        .map_or(0, u64::cast_signed);
    if neg { -q } else { q }
}

fn q126(neg: bool, e: i32, top: u128) -> i128 {
    let sh = 2 - i64::from(e);
    let q = u32::try_from(sh)
        .ok()
        .and_then(|s| top.checked_shr(s))
        .map_or(0, u128::cast_signed);
    if neg { -q } else { q }
}

/// `a·b·2^−126` truncated, for values below 2 in magnitude.
fn mul_q126(a: i128, b: i128) -> i128 {
    let (ua, ub) = (a.unsigned_abs(), b.unsigned_abs());
    let (ah, al) = (ua >> 64, ua & u128::from(u64::MAX));
    let (bh, bl) = (ub >> 64, ub & u128::from(u64::MAX));
    let ll = al * bl;
    let (mid, carry) = (al * bh).overflowing_add(ah * bl);
    let (lo, c2) = ll.overflowing_add(mid << 64);
    let hi = ah * bh + (mid >> 64) + (u128::from(carry) << 64) + u128::from(c2);
    let r = ((hi << 2) | (lo >> 126)).cast_signed();
    if (a < 0) == (b < 0) { r } else { -r }
}

/// A signed value with 190 fractional bits in three limbs, below 2 in magnitude.
type Q190 = (bool, [u64; 3]);

fn q190(neg: bool, e: i32, m: &[u64]) -> Q190 {
    let mut ext = [0; SCRATCH];
    for (i, &v) in m.iter().rev().take(3).rev().enumerate() {
        set(&mut ext, i, v);
    }
    let s = 2 - i64::from(e);
    match u32::try_from(s) {
        Ok(s) if s < 192 => {
            shr(&mut ext, 3, s);
        }
        _ => return (false, [0; 3]),
    }
    (neg, [limb(&ext, 0), limb(&ext, 1), limb(&ext, 2)])
}

fn mul_q190(a: Q190, b: Q190) -> Q190 {
    let mut p = [0u64; 6];
    for (i, &x) in a.1.iter().enumerate() {
        let mut carry: u128 = 0;
        for (j, &y) in b.1.iter().enumerate() {
            let cur = u128::from(limb(&p, i + j)) + u128::from(x) * u128::from(y) + carry;
            set(&mut p, i + j, lo64(cur));
            carry = cur >> 64;
        }
        set(&mut p, i + 3, lo64(carry));
    }
    let w = |i: usize| (limb(&p, i) >> 62) | (limb(&p, i + 1) << 2);
    (a.0 != b.0, [w(2), w(3), w(4)])
}

fn split190(m: [u64; 3]) -> (u128, u64) {
    ((u128::from(m[2]) << 64) | u128::from(m[1]), m[0])
}

fn join190((hi, lo): (u128, u64)) -> [u64; 3] {
    [lo, lo64(hi), hi64(hi)]
}

fn add_q190(a: Q190, b: Q190) -> Q190 {
    let ((ah, al), (bh, bl)) = (split190(a.1), split190(b.1));
    if a.0 == b.0 {
        let (lo, carry) = al.overflowing_add(bl);
        return (
            a.0,
            join190((ah.wrapping_add(bh).wrapping_add(u128::from(carry)), lo)),
        );
    }
    let ((xh, xl), (yh, yl), neg) = if (ah, al) < (bh, bl) {
        ((bh, bl), (ah, al), b.0)
    } else {
        ((ah, al), (bh, bl), a.0)
    };
    let (lo, borrow) = xl.overflowing_sub(yl);
    (
        neg,
        join190((xh.wrapping_sub(yh).wrapping_sub(u128::from(borrow)), lo)),
    )
}

/// `2^70/(2i + 1)`, `1/f` in Q62 at the middle of `[i/128, (i + 1)/128)`, for `i` from 64.
const RECIP_SEEDS: [u64; 64] = {
    let mut t = [0u64; 64];
    let mut rest: &mut [u64] = &mut t;
    let mut i = 64u128;
    while let [slot, tail @ ..] = rest {
        *slot = low64((1u128 << 70) / (2 * i + 1));
        rest = tail;
        i += 1;
    }
    t
};

/// `1/√f` in Q62 at the middle of `[i/128, (i + 1)/128)`, for `i` from 16.
const RSQRT_SEEDS: [u64; 112] = {
    let mut t = [0u64; 112];
    let mut rest: &mut [u64] = &mut t;
    let mut i = 16u128;
    while let [slot, tail @ ..] = rest {
        *slot = low64(((1u128 << 124) / (2 * i + 1)).isqrt() << 4);
        rest = tail;
        i += 1;
    }
    t
};

const fn low64(x: u128) -> u64 {
    let b = x.to_le_bytes();
    u64::from_le_bytes([b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7]])
}

/// `2^128/t` to about 2^−57, `t` in `[2^63, 2^64)`: a seed, then four steps of `g ← g(2 − fg)`.
fn recip_seed(t: u64) -> u128 {
    let index = usize::try_from(t >> 57).unwrap_or(64).saturating_sub(64);
    let mut g = RECIP_SEEDS.get(index).copied().unwrap_or(1 << 62);
    for _ in 0..4 {
        let fg = lo64((u128::from(t) * u128::from(g)) >> 64);
        g = lo64((u128::from(g) * u128::from((2u64 << 62).wrapping_sub(fg))) >> 62);
    }
    u128::from(g) << 2
}

/// `2^128/√x` to about 2^−57, `x` in `[2^125, 2^128)`: a seed, then four steps of `g ← g(3 − fg²)/2`.
fn rsqrt_seed(x: u128) -> u128 {
    let s = hi64(x);
    let index = usize::try_from(s >> 57).unwrap_or(16).saturating_sub(16);
    let mut g = RSQRT_SEEDS.get(index).copied().unwrap_or(1 << 62);
    for _ in 0..4 {
        let g2 = hi64(u128::from(g) * u128::from(g));
        let fg2 = u128::from(s) * u128::from(g2);
        let rest = hi64((3u128 << 124).wrapping_sub(fg2));
        g = lo64((u128::from(g) * u128::from(rest)) >> 61);
    }
    u128::from(g) << 2
}

/// `⌊√x⌋`, `x` in `[2^125, 2^128)`: `x·g/2^128`, a Newton correction by `g`, then exact steps.
fn isqrt_top(x: u128) -> u128 {
    let g = lo64(rsqrt_seed(x) >> 2);
    let estimate = (u128::from(hi64(x)) * u128::from(g)) >> 62;
    let error = x
        .cast_signed()
        .wrapping_sub(estimate.wrapping_mul(estimate).cast_signed());
    let step = ((error >> 32) * i128::from(g)) >> 95;
    let mut r = estimate
        .saturating_add_signed(step)
        .min(u128::from(u64::MAX));
    for _ in 0..4 {
        if r * r > x {
            r -= 1;
        } else if (r + 1).checked_mul(r + 1).is_some_and(|v| v <= x) {
            r += 1;
        }
    }
    r
}

type Wide = (u128, u128);

/// `x >> s` and whether a set bit fell off, for any `s`.
const fn wide_shr((hi, lo): Wide, s: u32) -> (Wide, bool) {
    match s {
        0 => ((hi, lo), false),
        1..=127 => (
            (hi >> s, (lo >> s) | (hi << (128 - s))),
            lo << (128 - s) != 0,
        ),
        128 => ((0, hi), lo != 0),
        129..=255 => ((0, hi >> (s - 128)), lo != 0 || hi << (256 - s) != 0),
        _ => ((0, 0), hi != 0 || lo != 0),
    }
}

const fn wide_shl((hi, lo): Wide, s: u32) -> Wide {
    match s {
        0 => (hi, lo),
        1..=127 => ((hi << s) | (lo >> (128 - s)), lo << s),
        128..=255 => (lo << (s - 128), 0),
        _ => (0, 0),
    }
}

const fn wide_leading_zeros((hi, lo): Wide) -> u32 {
    if hi == 0 {
        128 + lo.leading_zeros()
    } else {
        hi.leading_zeros()
    }
}

/// What the functions need from the working arithmetic, for the three widths.
pub trait Work:
    Copy
    + Debug
    + Ord
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Div<Output = Self>
    + Neg<Output = Self>
{
    /// Mantissa bits.
    const PREC: u32;
    fn zero() -> Self;
    fn one() -> Self;
    fn from_int(i: i64) -> Self;
    /// `raw · 2^−frac`, exact for any raw value of the type it serves.
    fn from_raw(raw: i128, frac: u32) -> Self;
    fn constant(c: &C) -> Self;
    /// `self · 2^frac` rounded to nearest even, or toward −∞, clamped to `[lo, hi]`.
    fn to_raw(self, frac: u32, lo: i128, hi: i128, floor: bool) -> i128;
    fn is_zero(&self) -> bool;
    fn is_neg(&self) -> bool;
    /// `⌊log₂|x|⌋ + 1`.
    fn exponent(&self) -> i32;
    #[must_use]
    fn abs(self) -> Self;
    #[must_use]
    fn with_sign(self, neg: bool) -> Self;
    /// `self · 2^k`.
    #[must_use]
    fn ldexp(self, k: i32) -> Self;
    #[must_use]
    fn square(self) -> Self;
    #[must_use]
    fn mul_int(self, k: i64) -> Self;
    /// Division by a positive integer, correctly rounded.
    #[must_use]
    fn div_int(self, k: u64) -> Self;
    #[must_use]
    fn sqrt(self) -> Self;
    /// For `mag · 2^−frac` turns, or radians when `radians`: the octant
    /// `⌊8·frac(t)⌋` and the fraction of it, measured from the octant's near
    /// end to π/4's multiple, so odd octants count back from their top.
    fn octant(mag: u128, frac: u32, radians: bool) -> (u32, Self);
    /// [`octant`](Self::octant) of `|self|` radians, below `2^PREC`.
    fn octant_of(self) -> (u32, Self);
    /// [`octant`](Self::octant) of `|a·b|` radians from the exact product, below `2^300`.
    fn octant_of_product(a: Self, b: Self) -> (u32, Self);
    /// `a·b + c·d` from the exact products, rounded once.
    fn dot2(a: Self, b: Self, c: Self, d: Self) -> Self;
    /// `Σ ±table[first + step·k]·u^k` for k from 0 to `n` by Horner's rule,
    /// the signs alternating when `alternate`, for `|u| < 1` and terms whose
    /// sum stays near 1, in fixed point with 62, 126 or 190 fractional bits.
    fn horner(u: Self, table: &Series, first: u32, step: u32, n: u32, alternate: bool) -> Self;

    fn round_int(self) -> i64 {
        let r = self.to_raw(0, i128::from(i64::MIN), i128::from(i64::MAX), false);
        i64::try_from(r).unwrap_or(0)
    }

    fn floor_int(self) -> i64 {
        let r = self.to_raw(0, i128::from(i64::MIN), i128::from(i64::MAX), true);
        i64::try_from(r).unwrap_or(0)
    }

    /// Whether `term` no longer moves `sum`.
    fn negligible(term: &Self, sum: &Self) -> bool {
        term.is_zero()
            || term.exponent() < sum.exponent().saturating_sub(Self::PREC.cast_signed() + 2)
    }
}

/// A binary floating-point value `m · 2^(e − 64N)`, `m`'s top bit set unless zero.
#[derive(Clone, Copy, Debug)]
pub struct Real<const N: usize> {
    neg: bool,
    e: i32,
    m: [u64; N],
}

impl<const N: usize> Real<N> {
    fn top(&self) -> u64 {
        self.m.last().copied().unwrap_or(0)
    }

    /// Rounds the top `N` of `len` normalized limbs by the rest and `sticky`.
    fn pack(neg: bool, e: i32, ext: &Scratch, len: usize, sticky: bool) -> Self {
        let low = len.saturating_sub(N);
        let mut m = [0u64; N];
        for (i, x) in m.iter_mut().enumerate() {
            *x = limb(ext, low + i);
        }
        let (half, rest) = if low == 0 {
            (false, sticky)
        } else {
            let g = limb(ext, low - 1);
            let below = (0..low - 1).any(|i| limb(ext, i) != 0);
            (g >> 63 == 1, g << 1 != 0 || below || sticky)
        };
        let odd = m.first().is_some_and(|&x| x & 1 == 1);
        let mut e = e;
        if half && (rest || odd) {
            let mut carry = true;
            for x in &mut m {
                let (v, c) = x.overflowing_add(u64::from(carry));
                *x = v;
                carry = c;
            }
            if carry {
                if let Some(t) = m.last_mut() {
                    *t = 1 << 63;
                }
                e = e.saturating_add(1);
            }
        }
        Self { neg, e, m }
    }

    /// The value `X · 2^(e − 64·len)` of the integer `X` in `len` limbs.
    fn normalize(neg: bool, e: i32, ext: &mut Scratch, len: usize, sticky: bool) -> Self {
        let z = leading_zeros(ext, len);
        if z >= u32_of(len) * 64 {
            return <Self as Work>::zero();
        }
        shl(ext, len, z);
        Self::pack(neg, e.saturating_sub(z.cast_signed()), ext, len, sticky)
    }

    /// `mag · 2^(e − 128)`.
    fn from_u128(neg: bool, mag: u128, e: i32) -> Self {
        if N == 1 {
            if mag == 0 {
                return <Self as Work>::zero();
            }
            let z = mag.leading_zeros();
            return Self::round128(neg, e.saturating_sub(z.cast_signed()), mag << z, false);
        }
        if mag == 0 {
            return <Self as Work>::zero();
        }
        let z = mag.leading_zeros();
        let top = mag << z;
        let mut m = [0u64; N];
        set(&mut m, N - 1, hi64(top));
        set(&mut m, N.saturating_sub(2), lo64(top));
        Self {
            neg,
            e: e.saturating_sub(z.cast_signed()),
            m,
        }
    }

    /// One limb: the value `x · 2^(e − 128)` of a normalized `x`, rounded.
    fn round128(neg: bool, e: i32, x: u128, sticky: bool) -> Self {
        let mut m = [0u64; N];
        let mut top = hi64(x);
        let low = lo64(x);
        let mut e = e;
        if low >> 63 == 1 && (low << 1 != 0 || sticky || top & 1 == 1) {
            top = top.wrapping_add(1);
            if top == 0 {
                top = 1 << 63;
                e = e.saturating_add(1);
            }
        }
        if let Some(t) = m.last_mut() {
            *t = top;
        }
        Self { neg, e, m }
    }

    fn add1(a: Self, b: Self) -> Self {
        let ma = u128::from(a.top()) << 64;
        let full = u128::from(b.top()) << 64;
        let d = a.e.abs_diff(b.e);
        let (mb, sticky) = if d < 128 {
            let mb = full >> d;
            (mb, mb << d != full)
        } else {
            (0, true)
        };
        if a.neg == b.neg {
            let (sum, carry) = ma.overflowing_add(mb);
            if carry {
                let x = (sum >> 1) | 1 << 127;
                Self::round128(a.neg, a.e.saturating_add(1), x, sticky || sum & 1 == 1)
            } else {
                Self::round128(a.neg, a.e, sum, sticky)
            }
        } else {
            let diff = ma - mb - u128::from(sticky);
            if diff == 0 {
                return <Self as Work>::zero();
            }
            let z = diff.leading_zeros();
            Self::round128(
                a.neg,
                a.e.saturating_sub(z.cast_signed()),
                diff << z,
                sticky,
            )
        }
    }

    /// `|m · 2^sh|` rounded to nearest even, or away from zero for a negative
    /// value when `floor`, for `sh ≤ 128 − 64N`.
    fn mag_scaled(&self, sh: i64, floor: bool) -> u128 {
        let mut ext = [0; SCRATCH];
        for (i, &v) in self.m.iter().enumerate() {
            set(&mut ext, i, v);
        }
        let len = N + 2;
        if sh >= 0 {
            shl(&mut ext, len, u32::try_from(sh).unwrap_or(0));
            return u128::from(limb(&ext, 0)) | (u128::from(limb(&ext, 1)) << 64);
        }
        let s = u32::try_from(-sh).unwrap_or(u32::MAX).min(64 * u32_of(len));
        let half = bit(&ext, s - 1);
        let rest = any_below(&ext, s - 1);
        shr(&mut ext, len, s);
        let q = u128::from(limb(&ext, 0)) | (u128::from(limb(&ext, 1)) << 64);
        let up = if floor {
            self.neg && (half || rest)
        } else {
            half && (rest || q & 1 == 1)
        };
        q.saturating_add(u128::from(up))
    }

    /// The exact product of the mantissas, in `2N` limbs.
    fn product(&self, o: &Self) -> Scratch {
        let mut p = [0; SCRATCH];
        for i in 0..N {
            let mut carry: u128 = 0;
            for j in 0..N {
                let cur = u128::from(limb(&p, i + j))
                    + u128::from(limb(&self.m, i)) * u128::from(limb(&o.m, j))
                    + carry;
                set(&mut p, i + j, lo64(cur));
                carry = cur >> 64;
            }
            set(&mut p, i + N, lo64(carry));
        }
        p
    }

    /// The octant and remainder of `mag · 2^−point` turns, or radians, `mag` in ≤ 3 limbs.
    fn octant_limbs(mag: &[u64], point: u32, radians: bool) -> (u32, Self) {
        let width = 64 * u32_of(N + 1);
        let mut z = [0; SCRATCH];
        let point = if radians {
            let skip = 6usize.saturating_sub(N + 3);
            for (i, &cv) in FRAC_1_2PI.2.iter().skip(skip).enumerate() {
                let mut carry: u128 = 0;
                for (j, &mv) in mag.iter().enumerate() {
                    let cur = u128::from(limb(&z, i + j)) + u128::from(cv) * u128::from(mv) + carry;
                    set(&mut z, i + j, lo64(cur));
                    carry = cur >> 64;
                }
                set(&mut z, i + mag.len(), lo64(carry));
            }
            386 + point - 64 * u32_of(skip)
        } else {
            for (i, &mv) in mag.iter().enumerate() {
                set(&mut z, i, mv);
            }
            point
        };
        let need = width + 3;
        if point < need {
            shl(&mut z, SCRATCH, need - point);
        } else {
            shr(&mut z, SCRATCH, point - need);
        }
        Self::octant_from(&z)
    }

    /// The octant and remainder of turns shifted to put the octant at bit `64(N + 1)`.
    fn octant_from(z: &[u64]) -> (u32, Self) {
        let oct = u32::try_from(limb(z, N + 1) & 7).unwrap_or(0);
        let mut v = [0; SCRATCH];
        for i in 0..=N {
            set(&mut v, i, limb(z, i));
        }
        if oct & 1 == 1 {
            if v.iter().all(|&x| x == 0) {
                return (oct, <Self as Work>::one());
            }
            let mut full = [0; SCRATCH];
            sub_into(&mut full, &v, N + 1);
            v = full;
        }
        (oct, Self::normalize(false, 0, &mut v, N + 1, false))
    }

    fn wide(&self) -> Wide {
        let (m0, m1, m2) = (limb(&self.m, 0), limb(&self.m, 1), limb(&self.m, 2));
        if N == 2 {
            ((u128::from(m1) << 64) | u128::from(m0), 0)
        } else {
            (
                (u128::from(m2) << 64) | u128::from(m1),
                u128::from(m0) << 64,
            )
        }
    }

    /// `x·2^(e − 256)`, `x`'s top bit set, rounded to nearest even, `sticky` for bits below.
    fn round_wide(neg: bool, e: i32, (hi, lo): Wide, sticky: bool) -> Self {
        let mut m = [0u64; N];
        let mut e = e;
        if N == 2 {
            let mut top = hi;
            if lo >> 127 == 1 && (lo << 1 != 0 || sticky || top & 1 == 1) {
                top = top.wrapping_add(1);
                if top == 0 {
                    top = 1 << 127;
                    e = e.saturating_add(1);
                }
            }
            set(&mut m, 0, lo64(top));
            set(&mut m, 1, hi64(top));
        } else {
            let (mut top, mut low) = (hi, hi64(lo));
            let below = lo64(lo);
            if below >> 63 == 1 && (below << 1 != 0 || sticky || low & 1 == 1) {
                low = low.wrapping_add(1);
                if low == 0 {
                    top = top.wrapping_add(1);
                    if top == 0 {
                        top = 1 << 127;
                        e = e.saturating_add(1);
                    }
                }
            }
            set(&mut m, 0, low);
            set(&mut m, 1, lo64(top));
            set(&mut m, 2, hi64(top));
        }
        Self { neg, e, m }
    }

    /// The correctly rounded sum for two or three limbs, `|big| ≥ |small|`, both non-zero.
    fn add_wide(big: Self, small: Self) -> Self {
        let gap = big.e.abs_diff(small.e);
        let ((ah, al), sticky) = wide_shr(small.wide(), gap);
        let (bh, bl) = big.wide();
        if big.neg == small.neg {
            let (lo, c1) = bl.overflowing_add(al);
            let (hi, c2) = bh.overflowing_add(ah);
            let (hi, c3) = hi.overflowing_add(u128::from(c1));
            if c2 || c3 {
                let shifted = ((hi >> 1) | 1 << 127, (lo >> 1) | (hi << 127));
                return Self::round_wide(
                    big.neg,
                    big.e.saturating_add(1),
                    shifted,
                    sticky || lo & 1 == 1,
                );
            }
            return Self::round_wide(big.neg, big.e, (hi, lo), sticky);
        }
        // big − small, less one unit if bits fell off small, which `sticky` then marks.
        let (lo, b1) = bl.overflowing_sub(al);
        let (lo, b2) = lo.overflowing_sub(u128::from(sticky));
        let hi = bh.wrapping_sub(ah).wrapping_sub(u128::from(b1 || b2));
        if hi == 0 && lo == 0 {
            return <Self as Work>::zero();
        }
        let z = wide_leading_zeros((hi, lo));
        Self::round_wide(
            big.neg,
            big.e.saturating_sub(z.cast_signed()),
            wide_shl((hi, lo), z),
            sticky,
        )
    }

    /// The correctly rounded product for two or three limbs of normalized mantissas.
    fn mul_wide(self, o: Self) -> Option<Self> {
        let neg = self.neg != o.neg;
        let e = self.e.saturating_add(o.e);
        let (acc, sticky) = if N == 2 {
            let (a, b) = (self.wide().0, o.wide().0);
            let (a1, a0, b1, b0) = (
                a >> 64,
                a & u128::from(u64::MAX),
                b >> 64,
                b & u128::from(u64::MAX),
            );
            let (mid, c1) = (a0 * b1).overflowing_add(a1 * b0);
            let (lo, c2) = (a0 * b0).overflowing_add(mid << 64);
            let hi = a1 * b1 + (mid >> 64) + (u128::from(c1) << 64) + u128::from(c2);
            ((hi, lo), false)
        } else {
            let mut p = [0u64; 6];
            for i in 0..3 {
                let x = u128::from(limb(&self.m, i));
                let mut carry: u128 = 0;
                for j in 0..3 {
                    let cur = u128::from(limb(&p, i + j)) + x * u128::from(limb(&o.m, j)) + carry;
                    set(&mut p, i + j, lo64(cur));
                    carry = cur >> 64;
                }
                set(&mut p, i + 3, lo64(carry));
            }
            let pair = |i: usize| (u128::from(limb(&p, i + 1)) << 64) | u128::from(limb(&p, i));
            let top = (pair(4), pair(2));
            let below = pair(0);
            if top.0 >> 127 == 0 {
                let shifted = wide_shl(top, 1);
                return (top.0 >> 126 == 1).then(|| {
                    Self::round_wide(
                        neg,
                        e.saturating_sub(1),
                        (shifted.0, shifted.1 | below >> 127),
                        below << 1 != 0,
                    )
                });
            }
            (top, below != 0)
        };
        if acc.0 >> 127 == 0 {
            return (acc.0 >> 126 == 1)
                .then(|| Self::round_wide(neg, e.saturating_sub(1), wide_shl(acc, 1), sticky));
        }
        Some(Self::round_wide(neg, e, acc, sticky))
    }

    fn mag_cmp(&self, o: &Self) -> Ordering {
        match (self.is_zero(), o.is_zero()) {
            (true, true) => return Ordering::Equal,
            (true, false) => return Ordering::Less,
            (false, true) => return Ordering::Greater,
            (false, false) => {}
        }
        self.e
            .cmp(&o.e)
            .then_with(|| self.top().cmp(&o.top()))
            .then_with(|| self.m.iter().rev().cmp(o.m.iter().rev()))
    }

    fn recip(self) -> Self {
        let top = self.top();
        if top == 0 {
            return self;
        }
        let one = <Self as Work>::one();
        if N == 1 {
            return one / self;
        }
        let b = self.abs();
        let mut y = Self::from_u128(false, recip_seed(top), 64i32.saturating_sub(self.e));
        for _ in 0..2 {
            y = y + y * (one - b * y);
        }
        y.with_sign(self.neg)
    }
}

impl<const N: usize> PartialEq for Real<N> {
    fn eq(&self, o: &Self) -> bool {
        self.cmp(o) == Ordering::Equal
    }
}

impl<const N: usize> Eq for Real<N> {}

impl<const N: usize> PartialOrd for Real<N> {
    fn partial_cmp(&self, o: &Self) -> Option<Ordering> {
        Some(self.cmp(o))
    }
}

impl<const N: usize> Ord for Real<N> {
    fn cmp(&self, o: &Self) -> Ordering {
        let sa = !self.is_zero() && self.neg;
        let sb = !o.is_zero() && o.neg;
        match (sa, sb) {
            (false, true) => Ordering::Greater,
            (true, false) => Ordering::Less,
            (false, false) => self.mag_cmp(o),
            (true, true) => o.mag_cmp(self),
        }
    }
}

impl<const N: usize> Neg for Real<N> {
    type Output = Self;
    fn neg(self) -> Self {
        if self.is_zero() {
            return self;
        }
        Self {
            neg: !self.neg,
            ..self
        }
    }
}

impl<const N: usize> Add for Real<N> {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        if N == 1 {
            let (top_a, top_b) = (self.top(), other.top());
            if top_b == 0 {
                return self;
            }
            if top_a == 0 {
                return other;
            }
            return if (self.e, top_a) < (other.e, top_b) {
                Self::add1(other, self)
            } else {
                Self::add1(self, other)
            };
        }
        if other.is_zero() {
            return self;
        }
        if self.is_zero() {
            return other;
        }
        let (big, small) = if self.mag_cmp(&other) == Ordering::Less {
            (other, self)
        } else {
            (self, other)
        };
        if N <= 3 {
            return Self::add_wide(big, small);
        }
        let len = N + 1;
        let mut acc = [0; SCRATCH];
        let mut addend = [0; SCRATCH];
        for i in 0..N {
            set(&mut acc, i + 1, limb(&big.m, i));
            set(&mut addend, i + 1, limb(&small.m, i));
        }
        let gap = i64::from(big.e) - i64::from(small.e);
        let sticky = if gap >= 64 * i64::from(u32_of(len)) {
            addend = [0; SCRATCH];
            true
        } else {
            shr(&mut addend, len, u32::try_from(gap).unwrap_or(0))
        };
        if big.neg == small.neg {
            if add_into(&mut acc, &addend, len) {
                let lost = shr(&mut acc, len, 1);
                let top = limb(&acc, len - 1) | 1 << 63;
                set(&mut acc, len - 1, top);
                Self::pack(big.neg, big.e.saturating_add(1), &acc, len, sticky || lost)
            } else {
                Self::pack(big.neg, big.e, &acc, len, sticky)
            }
        } else {
            sub_into(&mut acc, &addend, len);
            if sticky {
                let mut one = [0; SCRATCH];
                set(&mut one, 0, 1);
                sub_into(&mut acc, &one, len);
            }
            Self::normalize(big.neg, big.e, &mut acc, len, sticky)
        }
    }
}

impl<const N: usize> Sub for Real<N> {
    type Output = Self;
    fn sub(self, o: Self) -> Self {
        self + -o
    }
}

impl<const N: usize> Mul for Real<N> {
    type Output = Self;
    fn mul(self, o: Self) -> Self {
        if self.is_zero() || o.is_zero() {
            return <Self as Work>::zero();
        }
        let neg = self.neg != o.neg;
        let e = self.e.saturating_add(o.e);
        if N == 1 {
            let p = u128::from(self.top()) * u128::from(o.top());
            return if p >> 127 == 1 {
                Self::round128(neg, e, p, false)
            } else {
                Self::round128(neg, e.saturating_sub(1), p << 1, false)
            };
        }
        if N <= 3
            && let Some(product) = self.mul_wide(o)
        {
            return product;
        }
        let mut p = self.product(&o);
        Self::normalize(neg, e, &mut p, 2 * N, false)
    }
}

impl<const N: usize> Div for Real<N> {
    type Output = Self;
    /// To within a few units in the last place; correctly rounded for one limb.
    fn div(self, o: Self) -> Self {
        if self.is_zero() || o.is_zero() {
            return <Self as Work>::zero();
        }
        if N == 1 {
            let num = u128::from(limb(&self.m, 0)) << 64;
            let den = u128::from(limb(&o.m, 0));
            let q = num / den;
            let mut ext = [0; SCRATCH];
            set(&mut ext, 0, lo64(q));
            set(&mut ext, 1, hi64(q));
            let e = self.e.saturating_sub(o.e).saturating_add(64);
            return Self::normalize(self.neg != o.neg, e, &mut ext, 2, num % den != 0);
        }
        self * o.recip()
    }
}

impl<const N: usize> Work for Real<N> {
    const PREC: u32 = match N {
        1 => 64,
        2 => 128,
        _ => 192,
    };

    fn zero() -> Self {
        Self {
            neg: false,
            e: 0,
            m: [0; N],
        }
    }

    fn one() -> Self {
        Self::from_int(1)
    }

    fn from_int(i: i64) -> Self {
        Self::from_raw(i128::from(i), 0)
    }

    fn from_raw(raw: i128, frac: u32) -> Self {
        Self::from_u128(raw < 0, raw.unsigned_abs(), 128 - frac.cast_signed())
    }

    fn constant(c: &C) -> Self {
        let (neg, e, limbs) = *c;
        let low = 6 - N;
        let mut m = [0u64; N];
        for (i, x) in m.iter_mut().enumerate() {
            *x = limb(&limbs, low + i);
        }
        let g = limb(&limbs, low.wrapping_sub(1));
        let rest = g << 1 != 0 || (0..low.saturating_sub(1)).any(|i| limb(&limbs, i) != 0);
        if g >> 63 == 1 && (rest || m.first().is_some_and(|&x| x & 1 == 1)) {
            let mut ext = [0; SCRATCH];
            for (i, &x) in m.iter().enumerate() {
                set(&mut ext, i + 1, x);
            }
            set(&mut ext, 0, 1 << 63);
            Self::pack(neg, e, &ext, N + 1, true)
        } else {
            Self { neg, e, m }
        }
    }

    fn to_raw(self, frac: u32, lo: i128, hi: i128, floor: bool) -> i128 {
        if self.is_zero() {
            return 0;
        }
        let saturated = if self.neg { lo } else { hi };
        let top = i64::from(self.e) + i64::from(frac);
        if top > 128 {
            return saturated;
        }
        let sh = top - i64::from(Self::PREC);
        let mag = if N == 1 {
            let m = u128::from(self.top());
            match u32::try_from(-sh) {
                Err(_) => m << u32::try_from(sh).unwrap_or(0),
                Ok(s) if s > 64 => u128::from(floor && self.neg),
                Ok(s) => {
                    let dropped = m & ((1u128 << s) - 1);
                    let q = m >> s;
                    let half = s > 0 && dropped >> (s - 1) == 1;
                    let rest = s > 1 && dropped & ((1u128 << (s - 1)) - 1) != 0;
                    let up = if floor {
                        self.neg && dropped != 0
                    } else {
                        half && (rest || q & 1 == 1)
                    };
                    q + u128::from(up)
                }
            }
        } else {
            self.mag_scaled(sh, floor)
        };
        if self.neg {
            0i128.checked_sub_unsigned(mag).map_or(lo, |v| v.max(lo))
        } else {
            i128::try_from(mag).map_or(hi, |v| v.min(hi))
        }
    }

    fn is_zero(&self) -> bool {
        // Normalized, so the top bit is set unless every bit is clear.
        self.top() == 0
    }

    fn is_neg(&self) -> bool {
        self.neg && !self.is_zero()
    }

    fn exponent(&self) -> i32 {
        self.e
    }

    fn abs(self) -> Self {
        Self { neg: false, ..self }
    }

    fn with_sign(self, neg: bool) -> Self {
        if neg { -self.abs() } else { self.abs() }
    }

    fn ldexp(self, k: i32) -> Self {
        if self.is_zero() {
            return self;
        }
        Self {
            e: self.e.saturating_add(k),
            ..self
        }
    }

    fn square(self) -> Self {
        self * self
    }

    fn mul_int(self, k: i64) -> Self {
        if self.is_zero() || k == 0 {
            return Self::zero();
        }
        let ku = k.unsigned_abs();
        let mut p = [0; SCRATCH];
        let mut carry: u128 = 0;
        for i in 0..N {
            let cur = u128::from(limb(&self.m, i)) * u128::from(ku) + carry;
            set(&mut p, i, lo64(cur));
            carry = cur >> 64;
        }
        set(&mut p, N, lo64(carry));
        Self::normalize(
            self.neg != (k < 0),
            self.e.saturating_add(64),
            &mut p,
            N + 1,
            false,
        )
    }

    fn div_int(self, k: u64) -> Self {
        if self.is_zero() || k == 0 {
            return self;
        }
        let mut q = [0; SCRATCH];
        let mut rem: u128 = 0;
        let small = u32::try_from(k).is_ok();
        for i in (0..=N).rev() {
            let digit = if i == 0 { 0 } else { limb(&self.m, i - 1) };
            if small {
                // rem < k < 2^32, so each half fits a 64-bit division.
                let high = (lo64(rem) << 32) | (digit >> 32);
                let low = ((high % k) << 32) | (digit & u64::from(u32::MAX));
                set(&mut q, i, ((high / k) << 32) | (low / k));
                rem = u128::from(low % k);
                continue;
            }
            let cur = (rem << 64) | u128::from(digit);
            set(&mut q, i, lo64(cur / u128::from(k)));
            rem = cur % u128::from(k);
        }
        Self::normalize(self.neg, self.e, &mut q, N + 1, rem != 0)
    }

    fn sqrt(self) -> Self {
        if self.is_zero() {
            return Self::zero();
        }
        let top2 = if N == 1 {
            u128::from(self.top()) << 64
        } else {
            (u128::from(self.top()) << 64) | u128::from(limb(&self.m, N.saturating_sub(2)))
        };
        let ee = i64::from(self.e) - 128;
        let (x, ee) = if ee % 2 == 0 {
            (top2, ee)
        } else {
            (top2 >> 1, ee + 1)
        };
        let half_e = i32_of(ee / 2);
        if N == 1 {
            let r = isqrt_top(x);
            let up = x - r * r > r;
            return Self::from_u128(false, r + u128::from(up), half_e.saturating_add(128));
        }
        let a = self.abs();
        let three = Self::from_int(3);
        let mut y = Self::from_u128(false, rsqrt_seed(x), 0i32.saturating_sub(half_e));
        for _ in 0..2 {
            y = (y * (three - a * y.square())).ldexp(-1);
        }
        a * y
    }

    fn octant(mag: u128, frac: u32, radians: bool) -> (u32, Self) {
        Self::octant_limbs(&[lo64(mag), hi64(mag)], frac, radians)
    }

    fn octant_of(self) -> (u32, Self) {
        let point = Self::PREC
            .cast_signed()
            .saturating_sub(self.e)
            .max(0)
            .cast_unsigned();
        Self::octant_limbs(&self.m, point, true)
    }

    fn octant_of_product(a: Self, b: Self) -> (u32, Self) {
        if a.is_zero() || b.is_zero() {
            return (0, Self::zero());
        }
        let prod = a.product(&b);
        let mut turns = [0u64; 16];
        for (i, &cv) in FRAC_1_2PI_LONG.iter().enumerate() {
            let mut carry: u128 = 0;
            for j in 0..2 * N {
                let cur = u128::from(limb(&turns, i + j))
                    + u128::from(cv) * u128::from(limb(&prod, j))
                    + carry;
                set(&mut turns, i + j, lo64(cur));
                carry = cur >> 64;
            }
            set(&mut turns, i + 2 * N, lo64(carry));
        }
        // |a·b| = prod·2^(e − 128N) and 1/2π = c·2^−578, so |a·b|/2π = turns·2^−point.
        let point = i64::from(128 * u32_of(N) + 578) - i64::from(a.e) - i64::from(b.e);
        let need = i64::from(64 * u32_of(N + 1) + 3);
        let shift = u32::try_from((point - need).unsigned_abs())
            .unwrap_or(u32::MAX)
            .min(1024);
        if point < need {
            shl(&mut turns, 16, shift);
        } else {
            shr(&mut turns, 16, shift);
        }
        Self::octant_from(&turns)
    }

    fn dot2(a: Self, b: Self, c: Self, d: Self) -> Self {
        if a.is_zero() || b.is_zero() {
            return c * d;
        }
        if c.is_zero() || d.is_zero() {
            return a * b;
        }
        let mut first = (a.neg != b.neg, a.e.saturating_add(b.e), a.product(&b));
        let mut second = (c.neg != d.neg, c.e.saturating_add(d.e), c.product(&d));
        if second.1 > first.1 {
            core::mem::swap(&mut first, &mut second);
        }
        let low = 64 * u32_of(SCRATCH - 1 - 2 * N);
        shl(&mut first.2, SCRATCH, low);
        shl(&mut second.2, SCRATCH, low);
        let sticky = shr(&mut second.2, SCRATCH, first.1.abs_diff(second.1).min(1024));
        let top = first.1.saturating_add(64);
        if first.0 == second.0 {
            add_into(&mut first.2, &second.2, SCRATCH);
            return Self::normalize(first.0, top, &mut first.2, SCRATCH, sticky);
        }
        if first.2.iter().rev().cmp(second.2.iter().rev()) == Ordering::Less {
            sub_into(&mut second.2, &first.2, SCRATCH);
            return Self::normalize(second.0, top, &mut second.2, SCRATCH, false);
        }
        sub_into(&mut first.2, &second.2, SCRATCH);
        if sticky {
            // What fell off `second` takes one more unit, and the rest rounds.
            let mut one = [0; SCRATCH];
            set(&mut one, 0, 1);
            sub_into(&mut first.2, &one, SCRATCH);
        }
        Self::normalize(first.0, top, &mut first.2, SCRATCH, sticky)
    }

    fn horner(u: Self, table: &Series, first: u32, step: u32, n: u32, alternate: bool) -> Self {
        let sign = |k: u32| alternate && k % 2 == 1;
        let at = |k: u32| (first + step * k) as usize;
        if N == 1 {
            let uq = q62(u.neg, u.e, u.top());
            let mut p: i64 = 0;
            for k in (0..=n).rev() {
                let c = table.q62.get(at(k)).copied().unwrap_or(0);
                let up = i64::try_from((i128::from(uq) * i128::from(p)) >> 62).unwrap_or(0);
                p = (if sign(k) { -c } else { c }) + up;
            }
            return Self::from_raw(i128::from(p), 62);
        }
        if N == 2 {
            let top = (u128::from(u.top()) << 64) | u128::from(limb(&u.m, 0));
            let uq = q126(u.neg, u.e, top);
            let mut p: i128 = 0;
            for k in (0..=n).rev() {
                let c = table.q126.get(at(k)).copied().unwrap_or(0);
                p = (if sign(k) { -c } else { c }) + mul_q126(uq, p);
            }
            return Self::from_raw(p, 126);
        }
        let uq = q190(u.neg, u.e, &u.m);
        let mut p: Q190 = (false, [0; 3]);
        for k in (0..=n).rev() {
            let c = table.q190.get(at(k)).copied().unwrap_or([0; 3]);
            p = add_q190((sign(k), c), mul_q190(uq, p));
        }
        let mut ext = [0; SCRATCH];
        for (i, &v) in p.1.iter().enumerate() {
            set(&mut ext, i, v);
        }
        Self::normalize(p.0, 2, &mut ext, 3, false)
    }
}

#[cfg(test)]
mod tests {
    use super::{isqrt_top, recip_seed, rsqrt_seed};

    #[test]
    fn seeds_and_roots_are_exact_enough() {
        let mut state = 0x9E37_79B9_7F4A_7C15u64;
        let mut next = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        for _ in 0..200_000 {
            let x = (u128::from(next()) << 64 | u128::from(next())) >> (next() % 3);
            let x = x.max(1 << 125);
            assert_eq!(isqrt_top(x), x.isqrt(), "{x}");
            let t = next() | 1 << 63;
            let inverse = u128::MAX / u128::from(t);
            assert!(recip_seed(t).abs_diff(inverse) < inverse >> 55, "{t}");
            let root = u128::MAX / x.isqrt();
            assert!(rsqrt_seed(x).abs_diff(root) < root >> 55, "{x}");
        }
    }
}
