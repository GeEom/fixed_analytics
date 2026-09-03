//! Trait definitions for types compatible with CORDIC algorithms.

use core::ops::{Add, AddAssign, Mul, Neg, Shl, Shr, Sub, SubAssign};
use fixed::traits::{Fixed, FixedSigned};
use fixed::types::extra::{IsLessOrEqual, LeEqU128, True, Unsigned};
use fixed::{FixedI8, FixedI16, FixedI32, FixedI64, FixedI128};

/// A number type that can be used with CORDIC-based algorithms.
///
/// This trait abstracts over fixed-point number types, providing the
/// operations and constants necessary for CORDIC computations.
///
/// # Implementors
///
/// This trait is implemented generically for all signed fixed-point types
/// from the `fixed` crate that have sufficient fractional bits to represent
/// the required constants (π, e, etc.):
///
/// - [`FixedI8<Fract>`](fixed::FixedI8) where Fract ≤ 5 (for π to fit)
/// - [`FixedI16<Fract>`](fixed::FixedI16) where Fract ≤ 13
/// - [`FixedI32<Fract>`](fixed::FixedI32) where Fract ≤ 29
/// - [`FixedI64<Fract>`](fixed::FixedI64) where Fract ≤ 61
/// - [`FixedI128<Fract>`](fixed::FixedI128) where Fract ≤ 125
///
/// Common type aliases like `I16F16`, `I32F32`, `I8F24`, `I24F8` all work.
pub trait CordicNumber:
    Copy
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
    /// Fractional bits. Determines CORDIC iteration count.
    fn frac_bits() -> u32;
    /// Total bits.
    fn total_bits() -> u32;
    /// Converts from a raw I1F63 representation (1 sign bit, 63 fractional bits).
    /// For constants in (-1, 1).
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
    /// Wrapping multiplication, for use only where overflow is provably
    /// impossible: then it matches [`saturating_mul`](Self::saturating_mul)
    /// bit for bit without the overflow check.
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
    /// Division by a positive integer, truncated toward zero. Matches
    /// [`div`](Self::div) bit for bit for non-negative `self` but divides the
    /// raw representation directly, which is several times cheaper on
    /// 128-bit types. A zero divisor saturates like [`div`](Self::div).
    #[must_use]
    fn div_int(self, divisor: u32) -> Self;
    /// Division.
    #[must_use]
    fn div(self, rhs: Self) -> Self;
    /// Integer part of the base-2 logarithm, `⌊log₂(self)⌋`, or `None` if
    /// `self ≤ 0`.
    fn checked_int_log2(self) -> Option<i32>;
    /// Square root, rounded to the nearest representable value. `self` must
    /// be non-negative; the `fixed` types return the root of the magnitude.
    #[must_use]
    fn sqrt_round(self) -> Self;
    /// Convert from numeric type.
    fn from_num<N: fixed::traits::ToFixed>(n: N) -> Self;
    /// Maximum value.
    fn max_value() -> Self;
    /// Minimum value.
    fn min_value() -> Self;
    /// Round to nearest integer (half away from zero), saturating at the
    /// type's bounds.
    #[must_use]
    fn round(self) -> Self;
    /// Convert to i32, rounding toward −∞ and saturating if the value does
    /// not fit.
    #[must_use]
    fn to_i32(self) -> i32;
}

// =============================================================================
// Square root helpers
// =============================================================================

/// Round-to-nearest square root when `X << (FRAC_NBITS + 2)` fits in
/// `$wide`: the raw result is `⌊√N + ½⌋ = (⌊2√N⌋ + 1) >> 1` with
/// `N = X·2^f`, and `⌊2√N⌋ = isqrt(4N)`.
macro_rules! sqrt_round_widened {
    ($self:expr, $bits:ty, $wide:ty) => {{
        let n4 = <$wide>::from($self.to_bits().unsigned_abs()) << (Self::FRAC_NBITS + 2);
        #[allow(
            clippy::cast_possible_truncation,
            clippy::cast_possible_wrap,
            reason = "the root has at most (total_bits + frac_bits) / 2 + 1 bits"
        )]
        Self::from_bits(((n4.isqrt() + 1) >> 1) as $bits)
    }};
}

/// Round-to-nearest square root for 128-bit raw representations, where
/// `X·2^f` does not fit in a machine integer. See [`sqrt_round_u128`].
macro_rules! sqrt_round_i128 {
    ($self:expr, $bits:ty, $wide:ty) => {{
        let root = sqrt_round_u128($self.to_bits().unsigned_abs(), Self::FRAC_NBITS, |a, b| {
            // ⌊a·2^f / b⌋ for b ≥ √N is below 2^127, so this cannot
            // overflow; the fallback only keeps the function total.
            #[allow(clippy::cast_possible_wrap, reason = "operands below 2^127")]
            let quotient =
                Fixed::checked_div(Self::from_bits(a as i128), Self::from_bits(b as i128));
            #[allow(clippy::cast_sign_loss, reason = "quotient of positive operands")]
            {
                quotient.map_or(i128::MAX, Self::to_bits) as u128
            }
        });
        #[allow(clippy::cast_possible_wrap, reason = "root below 2^127")]
        Self::from_bits(root as $bits)
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

/// Round-to-nearest square root of `x / 2^frac` (`x < 2^127`), in raw units:
/// `⌊√N + ½⌋` with `N = x·2^frac`, up to 254 bits.
///
/// `x` is shifted left by the largest `s ≤ leading_zeros(x)` with
/// `s ≡ frac (mod 2)`, giving `q = isqrt(x·2^s)` with 64 significant bits.
/// If `s ≥ frac`, `q >> ((s − frac)/2)` is exactly `⌊√N⌋`. Otherwise
/// `seed = (q + 1) << m` lies above `√N` by less than `2^(m+1)`, one Newton
/// step from above lands on `⌊√N⌋` or `⌊√N⌋ + 1` (its error is below
/// `2^(m − 64)`, `m ≤ 63`), and a 256-bit remainder settles floor and
/// rounding. `div_scaled(a, b)` must return `⌊a·2^frac / b⌋`; it is only
/// called with `b ≥ √N`.
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
            // ⌊2√N⌋ = root_hi >> (excess − 1); round half up.
            return ((root_hi >> (excess - 1)) + 1) >> 1;
        }
        // root_hi = ⌊√N⌋, and N fits in 128 bits since shift = frac ≤ lz.
        let root_sq = root_hi * root_hi;
        let rem = (x << frac) - root_sq;
        return if rem > root_hi { root_hi + 1 } else { root_hi };
    }

    let missing = (frac - shift) / 2;
    #[allow(clippy::cast_sign_loss, reason = "i128::MAX is positive")]
    let seed = ((root_hi + 1) << missing).min(i128::MAX as u128);
    let newton = u128::midpoint(seed, div_scaled(x, seed));

    let n = wide_shl(x, frac);
    let newton_sq = wide_mul(newton, newton);
    let (root, rem) = if wide_gt(newton_sq, n) {
        // newton = ⌊√N⌋ + 1: N − (newton − 1)² = (2·newton − 1) − (newton² − N).
        (newton - 1, (2 * newton - 1) - wide_sub_lo(newton_sq, n))
    } else {
        (newton, wide_sub_lo(n, newton_sq))
    };
    // rem = N − root² ∈ [0, 2·root]; round up exactly when rem > root.
    if rem > root { root + 1 } else { root }
}

// =============================================================================
// Generic implementations using macros
// =============================================================================

/// Macro to implement `CordicNumber` for `FixedI*` types generically.
///
/// The bounds ensure:
/// - `Fract` fits within the type (e.g., ≤ 8 for `FixedI8`)
/// - `Fract` allows π to be represented (needs ~2 integer bits)
/// - `Fract` allows π/2 to be represented (needs ~1 integer bit)
/// - `Fract` allows π/4 and ln(2) to be represented (needs ~1 integer bit)
macro_rules! impl_cordic_generic {
    (
        $fixed_type:ident,
        $bits_type:ty,
        $wide_type:ty,     // Unsigned type holding bits << (frac + 2) for sqrt
        $sqrt_impl:ident,  // sqrt_round_widened or sqrt_round_i128
        $total_bits:expr,
        $max_frac:ty,      // Maximum fractional bits for the type
        $pi_frac:ty,       // Max frac bits where PI fits (total - 2)
        $frac_pi_2:ty,     // Max frac bits where FRAC_PI_2 fits (total - 1)
        $frac_pi_4:ty      // Max frac bits where FRAC_PI_4 and LN_2 fit
    ) => {
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
            // Casts are safe: frac_bits ≤ 128, shift amounts bounded by type size
            #[allow(
                clippy::cast_possible_wrap,
                clippy::cast_lossless,
                reason = "frac_bits bounded by type size"
            )]
            fn from_i1f63(bits: i64) -> Self {
                // Convert from I1F63 representation to our type.
                // I1F63 has 63 fractional bits.
                // FRAC_NBITS is at most 128, which fits in i32.
                let our_frac = Self::FRAC_NBITS as i32;
                let shift = 63 - our_frac;

                if shift >= 0 {
                    // We have fewer frac bits than I1F63, shift right
                    #[allow(
                        clippy::cast_possible_truncation,
                        reason = "intentional truncation to target type"
                    )]
                    Self::from_bits((bits >> shift) as $bits_type)
                } else {
                    // We have more frac bits than I1F63, shift left
                    // Must cast first to avoid losing sign bit
                    #[allow(
                        clippy::cast_possible_truncation,
                        reason = "intentional truncation to target type"
                    )]
                    let wide = bits as $bits_type;
                    Self::from_bits(wide << (-shift))
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
                    // |rhs| exceeds the raw type's range (8- and 16-bit types
                    // only), so the exact product saturates unless self is 0.
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
                        // Division by zero: saturate based on sign.
                        None => {
                            if self.is_negative() {
                                Self::MIN
                            } else {
                                Self::MAX
                            }
                        }
                    },
                    // The divisor exceeds the raw type's range, so it exceeds
                    // |self| in raw units and the truncated quotient is zero.
                    Err(_) => Self::ZERO,
                }
            }

            #[inline]
            fn div(self, rhs: Self) -> Self {
                match Fixed::checked_div(self, rhs) {
                    Some(v) => v,
                    // Division by zero or overflow: saturate based on sign agreement.
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

// Import the specific U* types we need for bounds
use fixed::types::extra::{
    U5, U6, U7, U8, U13, U14, U15, U16, U29, U30, U31, U32, U61, U62, U63, U64, U125, U126, U127,
    U128,
};

// FixedI8<Fract>: 8 total bits
// - Max Fract: U8 (8 fractional bits = I0F8)
// - For PI (~3.14), need 2 integer bits, so Fract ≤ 6 (I2F6)
// - For FRAC_PI_2, FRAC_PI_4, LN_2, need 1 integer bit, so Fract ≤ 7 (I1F7)
// Being conservative: require Fract ≤ 5 so we have headroom
impl_cordic_generic!(FixedI8, i8, u16, sqrt_round_widened, 8, U8, U5, U6, U7);

// FixedI16<Fract>: 16 total bits
// - For PI, need Fract ≤ 14 (I2F14)
// - For FRAC_PI_2, FRAC_PI_4, LN_2, need Fract ≤ 15 (I1F15)
// - Conservative: Fract ≤ 13
impl_cordic_generic!(
    FixedI16,
    i16,
    u32,
    sqrt_round_widened,
    16,
    U16,
    U13,
    U14,
    U15
);

// FixedI32<Fract>: 32 total bits
// - For PI, need Fract ≤ 30
// - For FRAC_PI_2, FRAC_PI_4, LN_2, need Fract ≤ 31
// - Conservative: Fract ≤ 29
impl_cordic_generic!(
    FixedI32,
    i32,
    u64,
    sqrt_round_widened,
    32,
    U32,
    U29,
    U30,
    U31
);

// FixedI64<Fract>: 64 total bits
// - For PI, need Fract ≤ 62
// - For FRAC_PI_2, FRAC_PI_4, LN_2, need Fract ≤ 63
// - Conservative: Fract ≤ 61
impl_cordic_generic!(
    FixedI64,
    i64,
    u128,
    sqrt_round_widened,
    64,
    U64,
    U61,
    U62,
    U63
);

// FixedI128<Fract>: 128 total bits
// - For PI, need Fract ≤ 126
// - For FRAC_PI_2, FRAC_PI_4, LN_2, need Fract ≤ 127
// - Conservative: Fract ≤ 125
// No wider machine integer exists, so sqrt takes the Newton path.
impl_cordic_generic!(
    FixedI128,
    i128,
    u128,
    sqrt_round_i128,
    128,
    U128,
    U125,
    U126,
    U127
);
