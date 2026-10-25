#![doc = include_str!("../../docs/chebyshev.md")]

use crate::traits::CordicNumber;

/// Evaluates I1F63 coefficients, highest degree first, by Horner's method.
///
/// The arithmetic wraps. For these tables and `0 ≤ x ≤ (π/4)²` every
/// intermediate stays below 1, where that matches saturating arithmetic.
#[inline]
pub fn horner<T: CordicNumber, const N: usize>(coeffs: &[i64; N], x: T) -> T {
    let mut iter = coeffs.iter();
    let mut result = T::from_i1f63(*iter.next().unwrap_or(&0));
    for &coeff in iter {
        result = T::from_i1f63(coeff).wrapping_add(x.wrapping_mul(result));
    }
    result
}

/// (sin(x) − x)/x³ in u = x², low precision.
#[rustfmt::skip]
pub const SIN_P_LO: [i64; 4] = [
     0x0000_16DB_E083_8A07,
    -0x0006_804E_9E4E_E633,
     0x0111_110D_EF2E_1C96,
    -0x1555_5555_45E0_ABF9,
];

/// (sin(x) − x)/x³ in u = x², high precision.
#[rustfmt::skip]
pub const SIN_P_HI: [i64; 7] = [
    -0x0000_0000_006A_C5F5,
     0x0000_0000_5848_5FC3,
    -0x0000_0035_CC8A_8259,
     0x0000_171D_E3A5_4609,
    -0x0006_8068_0680_664E,
     0x0111_1111_1111_1100,
    -0x1555_5555_5555_5555,
];

/// (cos(x) − 1)/x² in u = x², low precision.
#[rustfmt::skip]
pub const COS_Q_LO: [i64; 4] = [
     0x0000_CD37_95D8_0CFA,
    -0x002D_81C1_0FEA_FDD5,
     0x0555_5532_ED10_3C5F,
    -0x3FFF_FFFF_563B_4B19,
];

/// (cos(x) − 1)/x² in u = x², high precision.
#[rustfmt::skip]
pub const COS_Q_HI: [i64; 7] = [
    -0x0000_0000_063F_E751,
     0x0000_0004_7BA9_FC96,
    -0x0000_024F_C9F1_C946,
     0x0000_D00D_00CE_F099,
    -0x002D_82D8_2D82_BAF2,
     0x0555_5555_5555_5435,
    -0x3FFF_FFFF_FFFF_FFFE,
];
