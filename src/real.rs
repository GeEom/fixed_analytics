//! The arithmetic every function evaluates in: binary floating point with an
//! `N`-limb mantissa, rounded to nearest even after each operation. A type
//! works at twice its width or more and rounds once to itself at the end.

use crate::consts::{C, FRAC_1_2PI};
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

fn shl(limbs: &mut Scratch, len: usize, shift: u32) {
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
fn shr(limbs: &mut Scratch, len: usize, shift: u32) -> bool {
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
    let mut p = [0; SCRATCH];
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

fn add_q190(a: Q190, b: Q190) -> Q190 {
    let (mut x, mut y) = ([0; SCRATCH], [0; SCRATCH]);
    for i in 0..3 {
        set(&mut x, i, limb(&a.1, i));
        set(&mut y, i, limb(&b.1, i));
    }
    let neg = if a.0 == b.0 {
        add_into(&mut x, &y, 3);
        a.0
    } else if a.1.iter().rev().cmp(b.1.iter().rev()) == Ordering::Less {
        sub_into(&mut y, &x, 3);
        x = y;
        b.0
    } else {
        sub_into(&mut x, &y, 3);
        a.0
    };
    (neg, [limb(&x, 0), limb(&x, 1), limb(&x, 2)])
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

/// A binary floating-point value `m · 2^(e − 64N)`, the top bit of `m` set
/// unless the value is zero.
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

    /// Rounds the top `N` of `len` normalized limbs, the rest and `sticky`
    /// deciding the rounding.
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
        let mut ext = [0; SCRATCH];
        set(&mut ext, N, lo64(mag));
        set(&mut ext, N + 1, hi64(mag));
        Self::normalize(neg, e, &mut ext, N + 2, false)
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

    fn mag_cmp(&self, o: &Self) -> Ordering {
        match (self.is_zero(), o.is_zero()) {
            (true, true) => return Ordering::Equal,
            (true, false) => return Ordering::Less,
            (false, true) => return Ordering::Greater,
            (false, false) => {}
        }
        self.e
            .cmp(&o.e)
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
        let mut y = Self::from_u128(
            false,
            u128::MAX / u128::from(top),
            64i32.saturating_sub(self.e),
        );
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
        self.m.iter().all(|&x| x == 0)
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
        for i in (0..=N).rev() {
            let digit = if i == 0 { 0 } else { limb(&self.m, i - 1) };
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
        let r = x.isqrt();
        let half_e = i32_of(ee / 2);
        if N == 1 {
            let up = x - r * r > r;
            return Self::from_u128(false, r + u128::from(up), half_e.saturating_add(128));
        }
        let a = self.abs();
        let three = Self::from_int(3);
        let mut y = Self::from_u128(false, u128::MAX / r.max(1), 0i32.saturating_sub(half_e));
        for _ in 0..2 {
            y = (y * (three - a * y.square())).ldexp(-1);
        }
        a * y
    }

    fn octant(mag: u128, frac: u32, radians: bool) -> (u32, Self) {
        let width = 64 * u32_of(N + 1);
        let mut z = [0; SCRATCH];
        let point = if radians {
            let skip = 6usize.saturating_sub(N + 3);
            for (i, &cv) in FRAC_1_2PI.2.iter().skip(skip).enumerate() {
                let mut carry: u128 = 0;
                for (j, mv) in [lo64(mag), hi64(mag)].into_iter().enumerate() {
                    let cur = u128::from(limb(&z, i + j)) + u128::from(cv) * u128::from(mv) + carry;
                    set(&mut z, i + j, lo64(cur));
                    carry = cur >> 64;
                }
                set(&mut z, i + 2, lo64(carry));
            }
            386 + frac - 64 * u32_of(skip)
        } else {
            set(&mut z, 0, lo64(mag));
            set(&mut z, 1, hi64(mag));
            frac
        };
        let need = width + 3;
        if point < need {
            shl(&mut z, SCRATCH, need - point);
        } else {
            shr(&mut z, SCRATCH, point - need);
        }
        let oct = u32::try_from(limb(&z, N + 1) & 7).unwrap_or(0);
        let mut v = [0; SCRATCH];
        for i in 0..=N {
            set(&mut v, i, limb(&z, i));
        }
        if oct & 1 == 1 {
            if v.iter().all(|&x| x == 0) {
                return (oct, Self::one());
            }
            let mut full = [0; SCRATCH];
            sub_into(&mut full, &v, N + 1);
            v = full;
        }
        (oct, Self::normalize(false, 0, &mut v, N + 1, false))
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
