//! Hyperbolic arctangent lookup table for hyperbolic CORDIC mode.
//!
//! Contains `atanh(2^-i)` for i = 1, 2, 3, ..., 64 in I1F63 format.
//!
//! # Important Notes
//!
//! - Hyperbolic CORDIC starts at i=1 (not i=0) because `atanh(1)` is undefined
//! - Certain iterations must be repeated for convergence (indices 4, 13, 40, ...)
//!
//! # Mathematical Background
//!
//! In hyperbolic CORDIC, we perform hyperbolic rotations using the identity:
//! ```text
//! [cosh(θ)  sinh(θ)] [x]   [x*cosh(θ) + y*sinh(θ)]
//! [sinh(θ)  cosh(θ)] [y] = [x*sinh(θ) + y*cosh(θ)]
//! ```
//!
//! By choosing θ = atanh(2^-i), the multiplication becomes a shift.
//!
//! # Table Generation
//!
//! Values are `round(atanh(2^-i) * 2^63)`, validated externally with 50-digit
//! precision. For indices 21+, `atanh(x) = x` to within I1F63 precision.

/// Precomputed hyperbolic arctangent values: `atanh(2^-i)` as I1F63 bit patterns.
///
/// Index 0 contains `atanh(2^-1)`, index 1 contains `atanh(2^-2)`, etc.
///
/// Note: `atanh(x) = 0.5 * ln((1+x)/(1-x))` for |x| < 1.
#[rustfmt::skip]
pub const ATANH_TABLE: [i64; 64] = [
    // atanh(2^-1) = atanh(0.5) ≈ 0.5493
    0x464F_A9EA_B40C_2A5E,  // 0.5493061443340549
    // atanh(2^-2) = atanh(0.25) ≈ 0.2554
    0x20B1_5DF5_0228_A34E,  // 0.25541281188299536
    // atanh(2^-3) = atanh(0.125) ≈ 0.1257
    0x1015_891C_9EAE_F76A,  // 0.12565721414045303
    // atanh(2^-4) ≈ 0.0626
    0x0802_AC45_69BA_D66E,  // 0.06258157147700301
    // atanh(2^-5)
    0x0400_5562_246B_B893,  // 0.031260178490666993
    // atanh(2^-6)
    0x0200_0AAB_1115_A393,  // 0.015626271752052209
    // atanh(2^-7)
    0x0100_0155_5888_91AD,  // 0.0078126589515404
    // atanh(2^-8)
    0x0080_002A_AAC4_4457,  // 0.003906269868396826
    // atanh(2^-9)
    0x0040_0005_5556_2222,  // 0.001953127483532550
    // atanh(2^-10)
    0x0020_0000_AAAA_B111,  // 0.0009765628104410357
    // atanh(2^-11)
    0x0010_0000_1555_5589,  // 0.0004882812888051129
    // atanh(2^-12)
    0x0008_0000_02AA_AAAC,  // 0.00024414062985063858
    // atanh(2^-13)
    0x0004_0000_0055_5555,  // 0.00012207031310632980
    // atanh(2^-14)
    0x0002_0000_000A_AAAB,  // 6.103515632579122e-05
    // atanh(2^-15)
    0x0001_0000_0001_5555,  // 3.051757813447391e-05
    // atanh(2^-16)
    0x0000_8000_0000_2AAB,  // 1.525878906368424e-05
    // atanh(2^-17)
    0x0000_4000_0000_0555,  // 7.629394531398029e-06
    // atanh(2^-18)
    0x0000_2000_0000_00AB,  // 3.814697265643503e-06
    // atanh(2^-19)
    0x0000_1000_0000_0015,  // 1.907348632814813e-06
    // atanh(2^-20)
    0x0000_0800_0000_0003,  // 9.536743164065390e-07
    // atanh(2^-21)
    0x0000_0400_0000_0000,  // 4.7683715820308847e-07
    // atanh(2^-22)
    0x0000_0200_0000_0000,  // 2.3841857910156701e-07
    // atanh(2^-23)
    0x0000_0100_0000_0000,  // 1.1920928955078181e-07
    // atanh(2^-24)
    0x0000_0080_0000_0000,  // 5.9604644775390696e-08
    // atanh(2^-25)
    0x0000_0040_0000_0000,  // 2.9802322387695321e-08
    // atanh(2^-26)
    0x0000_0020_0000_0000,  // 1.4901161193847657e-08
    // atanh(2^-27)
    0x0000_0010_0000_0000,  // 7.4505805969238283e-09
    // atanh(2^-28)
    0x0000_0008_0000_0000,  // 3.7252902984619141e-09
    // atanh(2^-29)
    0x0000_0004_0000_0000,  // 1.8626451492309570e-09
    // atanh(2^-30)
    0x0000_0002_0000_0000,  // 9.3132257461547852e-10
    // atanh(2^-31)
    0x0000_0001_0000_0000,  // 4.6566128730773926e-10
    // Remaining entries for high precision
    0x0000_0000_8000_0000,  // atanh(2^-32)
    0x0000_0000_4000_0000,  // atanh(2^-33)
    0x0000_0000_2000_0000,  // atanh(2^-34)
    0x0000_0000_1000_0000,  // atanh(2^-35)
    0x0000_0000_0800_0000,  // atanh(2^-36)
    0x0000_0000_0400_0000,  // atanh(2^-37)
    0x0000_0000_0200_0000,  // atanh(2^-38)
    0x0000_0000_0100_0000,  // atanh(2^-39)
    0x0000_0000_0080_0000,  // atanh(2^-40)
    0x0000_0000_0040_0000,  // atanh(2^-41)
    0x0000_0000_0020_0000,  // atanh(2^-42)
    0x0000_0000_0010_0000,  // atanh(2^-43)
    0x0000_0000_0008_0000,  // atanh(2^-44)
    0x0000_0000_0004_0000,  // atanh(2^-45)
    0x0000_0000_0002_0000,  // atanh(2^-46)
    0x0000_0000_0001_0000,  // atanh(2^-47)
    0x0000_0000_0000_8000,  // atanh(2^-48)
    0x0000_0000_0000_4000,  // atanh(2^-49)
    0x0000_0000_0000_2000,  // atanh(2^-50)
    0x0000_0000_0000_1000,  // atanh(2^-51)
    0x0000_0000_0000_0800,  // atanh(2^-52)
    0x0000_0000_0000_0400,  // atanh(2^-53)
    0x0000_0000_0000_0200,  // atanh(2^-54)
    0x0000_0000_0000_0100,  // atanh(2^-55)
    0x0000_0000_0000_0080,  // atanh(2^-56)
    0x0000_0000_0000_0040,  // atanh(2^-57)
    0x0000_0000_0000_0020,  // atanh(2^-58)
    0x0000_0000_0000_0010,  // atanh(2^-59)
    0x0000_0000_0000_0008,  // atanh(2^-60)
    0x0000_0000_0000_0004,  // atanh(2^-61)
    0x0000_0000_0000_0002,  // atanh(2^-62)
    0x0000_0000_0000_0001,  // atanh(2^-63)
    0x0000_0000_0000_0001,  // atanh(2^-64) - effectively rounds to 1 LSB
];

/// The CORDIC gain factor for hyperbolic mode (`K_h`).
///
/// `K_h` = ∏(i=1 to ∞) sqrt(1 - 2^(-2i)) * ∏(repeat factors)
///
/// With repeat iterations at indices 4, 13, 40, 121, 364:
/// `K_h` ≈ 0.82815936096021562708
/// `1/K_h` ≈ 1.2074970677630722
///
/// After hyperbolic CORDIC iterations, results are scaled by `1/K_h`.
/// To compensate, we divide by `K_h` (equivalent to multiplying by `1/K_h`).
///
/// Stored as I1F63 bit pattern.
pub const HYPERBOLIC_GAIN: i64 = 0x6A01_203D_99A6_3986; // ≈ 0.82815936096021563

/// Indices that must be repeated for hyperbolic CORDIC convergence.
///
/// In hyperbolic mode, the algorithm only converges if certain iterations
/// are repeated. These follow the pattern: 4, 13, 40, 121, 364, ...
/// where each subsequent index is 3*prev + 1.
///
/// For practical purposes with 64-bit precision, we only need the first few.
pub const REPEAT_INDICES: [u32; 5] = [4, 13, 40, 121, 364];

/// Check if an iteration index needs to be repeated in hyperbolic mode.
///
/// # Arguments
///
/// * `index` - The current iteration index (1-based for hyperbolic mode)
///
/// # Returns
///
/// `true` if this iteration should be performed twice for convergence.
#[inline]
#[must_use]
pub const fn needs_repeat(index: u32) -> bool {
    matches!(index, 4 | 13 | 40 | 121 | 364)
}

/// The inverse CORDIC gain factor for hyperbolic mode (`1/K_h`) in I2F62 format.
///
/// After hyperbolic CORDIC iterations, results are scaled by `1/K_h` ≈ 1.2075.
/// Precomputing this avoids a runtime division in `sinh_cosh`.
///
/// Value: `round((1/K_h) * 2^62)` where `K_h` includes repeat factors, validated externally.
/// Stored as I2F62 bit pattern (since `1/K_h` > 1, I1F63 would overflow).
pub const HYPERBOLIC_GAIN_INV: i64 = 0x4D47_A1C8_03BB_08CA; // ≈ 1.2075 (I2F62)

/// Precomputed value of `atanh(0.5)` for argument reduction.
///
/// Used in atanh argument reduction: `atanh(x) = atanh(0.5) + atanh((x-0.5)/(1-0.5*x))`
///
/// Value: `round(atanh(0.5) * 2^63)`, validated externally.
/// This equals `ATANH_TABLE[0]` but is provided as a named constant for clarity.
pub const ATANH_HALF: i64 = 0x464F_A9EA_B40C_2A5E; // ≈ 0.5493061443340549
