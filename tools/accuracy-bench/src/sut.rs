//! The fixed-point types the worst-case sweep runs on, with exact conversion
//! of their values to MPFR.

use fixed::traits::Fixed;
use fixed_analytics::CordicNumber;
use rug::{Float, Integer};

/// Working precision of every reference value, enough for 128-bit inputs.
pub const PREC: u32 = 320;

/// A type under test: any signed `fixed` type the crate implements.
pub trait Sut: CordicNumber + Fixed + Send + Sync + core::fmt::Debug + 'static {
    fn raw(self) -> i128;
    fn from_raw(raw: i128) -> Self;

    fn bits() -> u32 {
        <Self as CordicNumber>::total_bits()
    }

    fn frac_n() -> u32 {
        <Self as CordicNumber>::frac_bits()
    }

    fn min_raw() -> i128 {
        <Self as CordicNumber>::min_value().raw()
    }

    fn max_raw() -> i128 {
        <Self as CordicNumber>::max_value().raw()
    }

    fn name() -> String {
        format!("I{}F{}", Self::bits() - Self::frac_n(), Self::frac_n())
    }

    fn to_float(self) -> Float {
        Float::with_val(PREC, Integer::from(self.raw())) >> Self::frac_n()
    }

    /// The nearest value to `x`, saturating at the ends of the range.
    fn from_float(x: &Float) -> Self {
        let scaled = Float::with_val(PREC, x << Self::frac_n());
        let raw = match scaled.to_integer() {
            Some(i) => i.to_i128().unwrap_or(if x.is_sign_negative() {
                i128::MIN
            } else {
                i128::MAX
            }),
            None => 0,
        };
        Self::from_raw(raw.clamp(Self::min_raw(), Self::max_raw()))
    }

    fn neg_sat(self) -> Self {
        CordicNumber::saturating_sub(<Self as CordicNumber>::zero(), self)
    }
}

impl<F> Sut for F
where
    F: CordicNumber + Fixed + Send + Sync + core::fmt::Debug + 'static,
    F::Bits: TryInto<i128> + TryFrom<i128>,
{
    fn raw(self) -> i128 {
        self.to_bits().try_into().ok().unwrap_or(0)
    }

    fn from_raw(raw: i128) -> Self {
        let clamped = raw.clamp(F::MIN.raw(), F::MAX.raw());
        F::from_bits(clamped.try_into().ok().unwrap_or_else(|| F::ZERO.to_bits()))
    }
}

/// One tested type: its name, width, and the sweep instantiated for it.
pub struct TypeEntry {
    pub name: String,
    pub bits: u32,
    pub run: fn(crate::sweep::Fun) -> crate::sweep::FunStats,
}

macro_rules! entries {
    ($($ty:ident < $($frac:ident),+ >);+ $(;)?) => {{
        use fixed::types::extra::*;
        use fixed::{FixedI8, FixedI16, FixedI32, FixedI64, FixedI128};
        let mut v: Vec<TypeEntry> = Vec::new();
        $($(
            v.push(TypeEntry {
                name: <$ty<$frac> as Sut>::name(),
                bits: <$ty<$frac> as Sut>::bits(),
                run: crate::sweep::run::<$ty<$frac>>,
            });
        )+)+
        v
    }};
}

/// Every 8-, 16- and 32-bit layout, and a spread of 64- and 128-bit ones
/// that covers both ends and the 63-bit edge of the I1F63 tables.
pub fn types() -> Vec<TypeEntry> {
    entries! {
        FixedI8<U0, U1, U2, U3, U4, U5>;
        FixedI16<U0, U1, U2, U3, U4, U5, U6, U7, U8, U9, U10, U11, U12, U13>;
        FixedI32<U0, U1, U2, U3, U4, U5, U6, U7, U8, U9, U10, U11, U12, U13, U14,
            U15, U16, U17, U18, U19, U20, U21, U22, U23, U24, U25, U26, U27, U28, U29>;
        FixedI64<U0, U1, U2, U3, U4, U8, U12, U16, U20, U24, U28, U31, U32, U33,
            U36, U40, U44, U48, U52, U56, U58, U59, U60, U61>;
        FixedI128<U0, U1, U2, U4, U8, U16, U24, U32, U40, U48, U56, U60, U61, U62,
            U63, U64, U65, U66, U72, U80, U88, U96, U104, U112, U118, U120, U122,
            U123, U124, U125>;
    }
}
