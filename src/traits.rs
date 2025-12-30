//! Trait definitions for types compatible with CORDIC algorithms.

use core::ops::{Add, AddAssign, Mul, Neg, Shl, Shr, Sub, SubAssign};
use fixed::traits::Fixed;
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
    /// The zero value.
    fn zero() -> Self;

    /// The value one.
    fn one() -> Self;

    /// The value two.
    fn two() -> Self;

    /// The value one-half (0.5).
    fn half() -> Self;

    /// The mathematical constant π.
    fn pi() -> Self;

    /// The mathematical constant π/2.
    fn frac_pi_2() -> Self;

    /// The mathematical constant π/4.
    fn frac_pi_4() -> Self;

    /// Euler's number e.
    fn e() -> Self;

    /// The natural logarithm of 2.
    fn ln_2() -> Self;

    /// The natural logarithm of 10.
    fn ln_10() -> Self;

    /// Returns the absolute value of `self`.
    #[must_use]
    fn abs(self) -> Self;

    /// Returns the number of fractional bits in this type.
    fn frac_bits() -> u32;

    /// Returns the total number of bits in this type.
    fn total_bits() -> u32;

    /// Converts from a 64-bit signed fractional representation (I1F63).
    fn from_i64_frac(bits: i64) -> Self;

    /// Converts from a 64-bit signed fractional representation (I2F62).
    ///
    /// This is used for constants that are >= 1 but < 4, such as `1/K_h ≈ 1.2075`.
    fn from_i2f62_frac(bits: i64) -> Self;

    /// Checks if this value is negative.
    fn is_negative(self) -> bool;

    /// Checks if this value is positive.
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

    /// Division.
    #[must_use]
    fn div(self, rhs: Self) -> Self;

    /// Convert from a floating-point value.
    fn from_num<N: fixed::traits::ToFixed>(n: N) -> Self;

    /// The maximum representable value.
    fn max_value() -> Self;

    /// The minimum representable value.
    fn min_value() -> Self;
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
            fn two() -> Self {
                Self::from_num(2)
            }

            #[inline]
            fn half() -> Self {
                Self::from_num(0.5)
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
            fn frac_pi_4() -> Self {
                Self::FRAC_PI_4
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
                if self.is_negative() { -self } else { self }
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
            #[allow(clippy::cast_possible_wrap, clippy::cast_lossless)]
            fn from_i64_frac(bits: i64) -> Self {
                // Convert from I1F63 representation to our type.
                // I1F63 has 63 fractional bits.
                // FRAC_NBITS is at most 128, which fits in i32.
                let our_frac = Self::FRAC_NBITS as i32;
                let shift = 63 - our_frac;

                if shift >= 0 {
                    // We have fewer frac bits than I1F63, shift right
                    #[allow(clippy::cast_possible_truncation)]
                    Self::from_bits((bits >> shift) as $bits_type)
                } else {
                    // We have more frac bits than I1F63, shift left
                    // Must cast first to avoid losing sign bit
                    #[allow(clippy::cast_possible_truncation)]
                    let wide = bits as $bits_type;
                    Self::from_bits(wide << (-shift))
                }
            }

            #[inline]
            #[allow(clippy::cast_possible_wrap, clippy::cast_lossless)]
            fn from_i2f62_frac(bits: i64) -> Self {
                // Convert from I2F62 representation to our type.
                // I2F62 has 62 fractional bits.
                let our_frac = Self::FRAC_NBITS as i32;
                let shift = 62 - our_frac;

                if shift >= 0 {
                    // We have fewer frac bits than I2F62, shift right
                    #[allow(clippy::cast_possible_truncation)]
                    Self::from_bits((bits >> shift) as $bits_type)
                } else {
                    // We have more frac bits than I2F62, shift left
                    #[allow(clippy::cast_possible_truncation)]
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
            fn div(self, rhs: Self) -> Self {
                self / rhs
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
impl_cordic_generic!(FixedI8, i8, 8, U8, U5, U6, U7);

// FixedI16<Fract>: 16 total bits
// - For PI, need Fract ≤ 14 (I2F14)
// - For FRAC_PI_2, FRAC_PI_4, LN_2, need Fract ≤ 15 (I1F15)
// - Conservative: Fract ≤ 13
impl_cordic_generic!(FixedI16, i16, 16, U16, U13, U14, U15);

// FixedI32<Fract>: 32 total bits
// - For PI, need Fract ≤ 30
// - For FRAC_PI_2, FRAC_PI_4, LN_2, need Fract ≤ 31
// - Conservative: Fract ≤ 29
impl_cordic_generic!(FixedI32, i32, 32, U32, U29, U30, U31);

// FixedI64<Fract>: 64 total bits
// - For PI, need Fract ≤ 62
// - For FRAC_PI_2, FRAC_PI_4, LN_2, need Fract ≤ 63
// - Conservative: Fract ≤ 61
impl_cordic_generic!(FixedI64, i64, 64, U64, U61, U62, U63);

// FixedI128<Fract>: 128 total bits
// - For PI, need Fract ≤ 126
// - For FRAC_PI_2, FRAC_PI_4, LN_2, need Fract ≤ 127
// - Conservative: Fract ≤ 125
impl_cordic_generic!(FixedI128, i128, 128, U128, U125, U126, U127);
