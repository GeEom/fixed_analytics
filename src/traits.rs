//! Trait definitions for types compatible with CORDIC algorithms.

use core::ops::{Add, AddAssign, Mul, Neg, Shl, Shr, Sub, SubAssign};
use fixed::traits::{Fixed, FixedSigned};
use fixed::types::extra::{IsLessOrEqual, LeEqU128, True, Unsigned};
use fixed::{FixedI8, FixedI16, FixedI32, FixedI64, FixedI128};

pub(crate) mod sealed {
    /// The working arithmetic of a type, and its raw bits; this crate's types
    /// alone implement it.
    pub trait Sealed {
        type Work: crate::real::Work;
        fn raw_i128(self) -> i128;
        fn from_raw_i128(raw: i128) -> Self;
    }
}

/// The operations and constants the algorithms need from a fixed-point type.
///
/// Every signed `fixed` type with three or more integer bits implements it,
/// from [`FixedI8`] to [`FixedI128`] with up
/// to 5, 13, 29, 61 or 125 fractional bits; `I16F16` and `I24F8` qualify. It
/// is sealed: no other type can implement it.
pub trait CordicNumber:
    sealed::Sealed
    + Copy
    + PartialEq
    + PartialOrd
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Neg<Output = Self>
    + AddAssign
    + SubAssign
    + Shl<u32, Output = Self>
    + Shr<u32, Output = Self>
{
    /// Zero.
    fn zero() -> Self;
    /// One.
    fn one() -> Self;
    /// Two.
    #[must_use]
    fn two() -> Self {
        Self::one() + Self::one()
    }
    /// Half.
    #[must_use]
    fn half() -> Self {
        Self::from_num(0.5)
    }
    /// π. Requires ≥2 integer bits.
    fn pi() -> Self;
    /// π/2. Requires ≥1 integer bit.
    fn frac_pi_2() -> Self;
    /// π/4.
    #[must_use]
    fn frac_pi_4() -> Self {
        Self::frac_pi_2() >> 1
    }
    /// Euler's number e.
    fn e() -> Self;
    /// ln(2).
    fn ln_2() -> Self;
    /// ln(10).
    fn ln_10() -> Self;
    /// Absolute value.
    #[must_use]
    fn abs(self) -> Self;
    /// Fractional bits, which set the CORDIC iteration count.
    fn frac_bits() -> u32;
    /// Total bits.
    fn total_bits() -> u32;
    /// Converts a constant in (-1, 1) from raw I1F63 bits, truncating.
    fn from_i1f63(bits: i64) -> Self;
    /// Returns true if negative.
    fn is_negative(self) -> bool;
    /// Returns true if positive.
    fn is_positive(self) -> bool {
        !self.is_negative() && self != Self::zero()
    }
    /// Saturating multiplication.
    #[must_use]
    fn saturating_mul(self, rhs: Self) -> Self;
    /// Saturating addition.
    #[must_use]
    fn saturating_add(self, rhs: Self) -> Self;
    /// Saturating subtraction.
    #[must_use]
    fn saturating_sub(self, rhs: Self) -> Self;
    /// Wrapping multiplication, only for where overflow is impossible and it
    /// matches [`saturating_mul`](Self::saturating_mul) without the check.
    #[must_use]
    fn wrapping_mul(self, rhs: Self) -> Self;
    /// Wrapping addition. Same contract as [`wrapping_mul`](Self::wrapping_mul).
    #[must_use]
    fn wrapping_add(self, rhs: Self) -> Self;
    /// Wrapping subtraction. Same contract as [`wrapping_mul`](Self::wrapping_mul).
    #[must_use]
    fn wrapping_sub(self, rhs: Self) -> Self;
    /// Saturating multiplication by an integer. Exact unless it saturates.
    #[must_use]
    fn mul_int(self, rhs: i32) -> Self;
    /// Division by a positive integer, truncated toward zero. For non-negative
    /// `self` it matches [`div`](Self::div) bit for bit, several times cheaper
    /// on 128-bit types; a zero divisor saturates the same way.
    #[must_use]
    fn div_int(self, divisor: u32) -> Self;
    /// Division.
    #[must_use]
    fn div(self, rhs: Self) -> Self;
    /// `⌊log₂(self)⌋`, or `None` if `self ≤ 0`.
    fn checked_int_log2(self) -> Option<i32>;
    /// Square root, rounded to nearest, of a non-negative `self`; the `fixed`
    /// types return the root of the magnitude.
    #[must_use]
    fn sqrt_round(self) -> Self;
    /// Convert from numeric type.
    fn from_num<N: fixed::traits::ToFixed>(n: N) -> Self;
    /// Maximum value.
    fn max_value() -> Self;
    /// Minimum value.
    fn min_value() -> Self;
    /// Round to the nearest integer, half away from zero, saturating.
    #[must_use]
    fn round(self) -> Self;
    /// Convert to i32, rounding toward −∞ and saturating.
    #[must_use]
    fn to_i32(self) -> i32;
}

// The square root derivations are in docs/sqrt.md.

/// Rounded square root when `X << (FRAC_NBITS + 2)` fits in `$wide`.
macro_rules! sqrt_round_widened {
    ($self:expr, $bits:ty, $wide:ty) => {{
        let n4 = <$wide>::from($self.to_bits().unsigned_abs()) << (Self::FRAC_NBITS + 2);
        Self::from_bits(<$bits>::try_from((n4.isqrt() + 1) >> 1).unwrap_or(<$bits>::MAX))
    }};
}

/// Rounded square root for 128-bit types, through [`sqrt_round_u128`].
macro_rules! sqrt_round_i128 {
    ($self:expr, $bits:ty, $wide:ty) => {{
        let root = sqrt_round_u128($self.to_bits().unsigned_abs(), Self::FRAC_NBITS, |a, b| {
            let quotient = Fixed::checked_div(
                Self::from_bits(a.cast_signed()),
                Self::from_bits(b.cast_signed()),
            );
            quotient.map_or(i128::MAX, Self::to_bits).cast_unsigned()
        });
        Self::from_bits(root.cast_signed())
    }};
}

/// A 256-bit unsigned integer as `(high, low)` limbs.
type U256 = (u128, u128);

/// `a · b` as a 256-bit product.
const fn wide_mul(a: u128, b: u128) -> U256 {
    const MASK: u128 = u64::MAX as u128;
    let (a_hi, a_lo) = (a >> 64, a & MASK);
    let (b_hi, b_lo) = (b >> 64, b & MASK);
    let ll = a_lo * b_lo;
    let lh = a_lo * b_hi;
    let hl = a_hi * b_lo;
    let hh = a_hi * b_hi;
    let (mid, mid_carry) = lh.overflowing_add(hl);
    let (lo, lo_carry) = ll.overflowing_add(mid << 64);
    let hi = hh + (mid >> 64) + ((mid_carry as u128) << 64) + (lo_carry as u128);
    (hi, lo)
}

/// `x << shift` as a 256-bit value, for `0 < shift < 128`.
const fn wide_shl(x: u128, shift: u32) -> U256 {
    (x >> (128 - shift), x << shift)
}

/// `a > b` for 256-bit values.
const fn wide_gt(a: U256, b: U256) -> bool {
    a.0 > b.0 || (a.0 == b.0 && a.1 > b.1)
}

/// Low limb of `a − b`, for `a ≥ b` with a difference below 2^128.
const fn wide_sub_lo(a: U256, b: U256) -> u128 {
    a.1.wrapping_sub(b.1)
}

/// `⌊√N + ½⌋` with `N = x·2^frac` for `x < 2^127`, by one Newton step and a
/// 256-bit remainder. `div_scaled(a, b)` must return `⌊a·2^frac / b⌋`.
fn sqrt_round_u128(x: u128, frac: u32, div_scaled: impl FnOnce(u128, u128) -> u128) -> u128 {
    if x == 0 {
        return 0;
    }
    let lz = x.leading_zeros();
    let shift = lz - ((lz ^ frac) & 1);
    let root_hi = (x << shift).isqrt();

    if shift >= frac {
        let excess = (shift - frac) / 2;
        if excess >= 1 {
            return ((root_hi >> (excess - 1)) + 1) >> 1;
        }
        // N fits in 128 bits since shift = frac ≤ lz.
        let root_sq = root_hi * root_hi;
        let rem = (x << frac) - root_sq;
        return if rem > root_hi { root_hi + 1 } else { root_hi };
    }

    let missing = (frac - shift) / 2;
    let seed = ((root_hi + 1) << missing).min(i128::MAX.cast_unsigned());
    let newton = u128::midpoint(seed, div_scaled(x, seed));

    let n = wide_shl(x, frac);
    let newton_sq = wide_mul(newton, newton);
    let (root, rem) = if wide_gt(newton_sq, n) {
        // N − (newton − 1)² = (2·newton − 1) − (newton² − N).
        (newton - 1, (2 * newton - 1) - wide_sub_lo(newton_sq, n))
    } else {
        (newton, wide_sub_lo(n, newton_sq))
    };
    if rem > root { root + 1 } else { root }
}

/// Implements `CordicNumber` for a `FixedI*` family, bounding `Fract` so the
/// type holds its constants: π needs two integer bits, π/2, π/4 and ln 2 one.
macro_rules! impl_cordic_generic {
    (
        $fixed_type:ident,
        $bits_type:ty,
        $wide_type:ty,
        $sqrt_impl:ident,
        $work:ty,
        $total_bits:expr,
        $max_frac:ty,
        $pi_frac:ty,
        $frac_pi_2:ty,
        $frac_pi_4:ty
    ) => {
        impl<Fract> sealed::Sealed for $fixed_type<Fract>
        where
            Fract: Unsigned
                + IsLessOrEqual<$max_frac, Output = True>
                + IsLessOrEqual<$pi_frac, Output = True>
                + IsLessOrEqual<$frac_pi_2, Output = True>
                + IsLessOrEqual<$frac_pi_4, Output = True>
                + LeEqU128,
        {
            type Work = $work;

            #[inline]
            fn raw_i128(self) -> i128 {
                i128::from(self.to_bits())
            }

            #[inline]
            fn from_raw_i128(raw: i128) -> Self {
                let saturated = if raw < 0 {
                    <$bits_type>::MIN
                } else {
                    <$bits_type>::MAX
                };
                Self::from_bits(<$bits_type>::try_from(raw).unwrap_or(saturated))
            }
        }

        impl<Fract> CordicNumber for $fixed_type<Fract>
        where
            Fract: Unsigned
                + IsLessOrEqual<$max_frac, Output = True>
                + IsLessOrEqual<$pi_frac, Output = True>
                + IsLessOrEqual<$frac_pi_2, Output = True>
                + IsLessOrEqual<$frac_pi_4, Output = True>
                + LeEqU128,
        {
            #[inline]
            fn zero() -> Self {
                Self::ZERO
            }

            #[inline]
            fn one() -> Self {
                Self::ONE
            }

            #[inline]
            fn pi() -> Self {
                Self::PI
            }

            #[inline]
            fn frac_pi_2() -> Self {
                Self::FRAC_PI_2
            }

            #[inline]
            fn e() -> Self {
                Self::E
            }

            #[inline]
            fn ln_2() -> Self {
                Self::LN_2
            }

            #[inline]
            fn ln_10() -> Self {
                Self::LN_10
            }

            #[inline]
            fn abs(self) -> Self {
                FixedSigned::saturating_abs(self)
            }

            #[inline]
            fn frac_bits() -> u32 {
                Self::FRAC_NBITS
            }

            #[inline]
            fn total_bits() -> u32 {
                $total_bits
            }

            #[inline]
            fn from_i1f63(bits: i64) -> Self {
                let shift = 63 - Self::FRAC_NBITS.cast_signed();
                // The trait bounds leave integer bits to spare, so a value in
                // (-1, 1) always fits and the fallback is dead.
                if shift >= 0 {
                    Self::from_bits(<$bits_type>::try_from(bits >> shift).unwrap_or(0))
                } else {
                    // Only 128-bit types get here; widen before the left shift.
                    Self::from_bits(<$bits_type>::try_from(bits).unwrap_or(0) << (-shift))
                }
            }

            #[inline]
            fn is_negative(self) -> bool {
                self < Self::ZERO
            }

            #[inline]
            fn saturating_mul(self, rhs: Self) -> Self {
                Fixed::saturating_mul(self, rhs)
            }

            #[inline]
            fn saturating_add(self, rhs: Self) -> Self {
                Fixed::saturating_add(self, rhs)
            }

            #[inline]
            fn saturating_sub(self, rhs: Self) -> Self {
                Fixed::saturating_sub(self, rhs)
            }

            #[inline]
            fn wrapping_mul(self, rhs: Self) -> Self {
                Fixed::wrapping_mul(self, rhs)
            }

            #[inline]
            fn wrapping_add(self, rhs: Self) -> Self {
                Fixed::wrapping_add(self, rhs)
            }

            #[inline]
            fn wrapping_sub(self, rhs: Self) -> Self {
                Fixed::wrapping_sub(self, rhs)
            }

            #[inline]
            fn mul_int(self, rhs: i32) -> Self {
                match <$bits_type>::try_from(rhs) {
                    Ok(k) => Fixed::saturating_mul_int(self, k),
                    // |rhs| is past the raw range, so only zero escapes saturation.
                    Err(_) => {
                        if self == Self::ZERO {
                            Self::ZERO
                        } else if self.is_negative() != (rhs < 0) {
                            Self::MIN
                        } else {
                            Self::MAX
                        }
                    }
                }
            }

            #[inline]
            fn div_int(self, divisor: u32) -> Self {
                match <$bits_type>::try_from(divisor) {
                    Ok(k) => match Fixed::checked_div_int(self, k) {
                        Some(v) => v,
                        None => {
                            if self.is_negative() {
                                Self::MIN
                            } else {
                                Self::MAX
                            }
                        }
                    },
                    // A divisor past the raw range exceeds |self|: the quotient is 0.
                    Err(_) => Self::ZERO,
                }
            }

            #[inline]
            fn div(self, rhs: Self) -> Self {
                match Fixed::checked_div(self, rhs) {
                    Some(v) => v,
                    None => {
                        if self.is_negative() != rhs.is_negative() {
                            Self::MIN
                        } else {
                            Self::MAX
                        }
                    }
                }
            }

            #[inline]
            fn checked_int_log2(self) -> Option<i32> {
                Fixed::checked_int_log2(self)
            }

            #[inline]
            fn sqrt_round(self) -> Self {
                $sqrt_impl!(self, $bits_type, $wide_type)
            }

            #[inline]
            fn from_num<N: fixed::traits::ToFixed>(n: N) -> Self {
                Self::from_num(n)
            }

            #[inline]
            fn max_value() -> Self {
                Self::MAX
            }

            #[inline]
            fn min_value() -> Self {
                Self::MIN
            }

            #[inline]
            fn round(self) -> Self {
                Fixed::saturating_round(self)
            }

            #[inline]
            fn to_i32(self) -> i32 {
                self.saturating_to_num::<i32>()
            }
        }
    };
}

use crate::real::Real;
use fixed::types::extra::{
    U5, U6, U7, U8, U13, U14, U15, U16, U29, U30, U31, U32, U61, U62, U63, U64, U125, U126, U127,
    U128,
};

// Each type's π bound is one bit stricter than π needs, for headroom.
impl_cordic_generic!(
    FixedI8,
    i8,
    u16,
    sqrt_round_widened,
    Real<1>,
    8,
    U8,
    U5,
    U6,
    U7
);

impl_cordic_generic!(
    FixedI16,
    i16,
    u32,
    sqrt_round_widened,
    Real<1>,
    16,
    U16,
    U13,
    U14,
    U15
);

impl_cordic_generic!(
    FixedI32,
    i32,
    u64,
    sqrt_round_widened,
    Real<1>,
    32,
    U32,
    U29,
    U30,
    U31
);

impl_cordic_generic!(
    FixedI64,
    i64,
    u128,
    sqrt_round_widened,
    Real<2>,
    64,
    U64,
    U61,
    U62,
    U63
);

impl_cordic_generic!(
    FixedI128,
    i128,
    u128,
    sqrt_round_i128,
    Real<3>,
    128,
    U128,
    U125,
    U126,
    U127
);
