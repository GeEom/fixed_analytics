//! Arctangent lookup table for circular CORDIC mode.
//!
//! Contains `atan(2^-i)` for i = 0, 1, 2, ..., 63 in I1F63 format.
//!
//! # Mathematical Background
//!
//! In circular CORDIC, we rotate a vector by successively smaller angles.
//! The angle at iteration i is `atan(2^-i)`, which allows the rotation
//! to be performed using only shifts and additions.
//!
//! # Gain Factor
//!
//! After n iterations, the vector magnitude is scaled by the gain factor:
//! ```text
//! K = ∏(i=0 to n-1) sqrt(1 + 2^(-2i)) ≈ 1.6468
//! ```
//!
//! We store the inverse (1/K ≈ 0.6073) to multiply into the initial vector.
//!
//! # Table Generation
//!
//! Values are `round(atan(2^-i) * 2^63)`, validated externally with 50-digit
//! precision. For indices 21+, `atan(x) = x` to within I1F63 precision.

/// Precomputed arctangent values: `atan(2^-i)` as I1F63 bit patterns.
///
/// These values are used in the circular CORDIC kernel for computing
/// trigonometric functions. The values represent angles in radians.
///
/// Index 0: atan(1) = π/4 ≈ 0.7854
/// Index 1: atan(0.5) ≈ 0.4636
/// Index 2: atan(0.25) ≈ 0.2450
/// ... and so on
#[rustfmt::skip]
pub const ATAN_TABLE: [i64; 64] = [
    // atan(2^0) = atan(1) = π/4
    0x6487_ED51_10B4_611A,  // 0.7853981633974483
    // atan(2^-1) = atan(0.5)
    0x3B58_CE0A_C376_9ED1,  // 0.4636476090008061
    // atan(2^-2) = atan(0.25)
    0x1F5B_75F9_2C80_DD63,  // 0.24497866312686414
    // atan(2^-3) = atan(0.125)
    0x0FEA_DD4D_5617_B6E3,  // 0.12435499454676144
    // atan(2^-4)
    0x07FD_56ED_CB3F_7A72,  // 0.06241880999595735
    // atan(2^-5)
    0x03FF_AAB7_752E_C495,  // 0.031239833430268277
    // atan(2^-6)
    0x01FF_F555_BBB7_29AB,  // 0.015623728620476831
    // atan(2^-7)
    0x00FF_FEAA_ADDD_D4B9,  // 0.007812341060101111
    // atan(2^-8)
    0x007F_FFD5_556E_EEDD,  // 0.0039062301319669718
    // atan(2^-9)
    0x003F_FFFA_AAAB_7777,  // 0.0019531225164788188
    // atan(2^-10)
    0x001F_FFFF_5555_5BBC,  // 0.0009765621895593195
    // atan(2^-11)
    0x000F_FFFF_EAAA_AADE,  // 0.0004882812111948983
    // atan(2^-12)
    0x0007_FFFF_FD55_5557,  // 0.00024414062014936177
    // atan(2^-13)
    0x0003_FFFF_FFAA_AAAB,  // 0.00012207031189367021
    // atan(2^-14)
    0x0001_FFFF_FFF5_5555,  // 6.103515617420877e-05
    // atan(2^-15)
    0x0000_FFFF_FFFE_AAAB,  // 3.0517578115526096e-05
    // atan(2^-16)
    0x0000_7FFF_FFFF_D555,  // 1.5258789061315762e-05
    // atan(2^-17)
    0x0000_3FFF_FFFF_FAAB,  // 7.62939453110197e-06
    // atan(2^-18)
    0x0000_1FFF_FFFF_FF55,  // 3.814697265606496e-06
    // For indices 19+, atan(x) = x to within I1F63 precision (validated externally).
    // atan(2^-19)
    0x0000_0FFF_FFFF_FFEB,  // 1.907348632810187e-06
    // atan(2^-20)
    0x0000_07FF_FFFF_FFFD,  // 9.536743164059608e-07
    // atan(2^-21)
    0x0000_0400_0000_0000,  // 4.768371582030888e-07
    // atan(2^-22)
    0x0000_0200_0000_0000,  // 2.384185791015579e-07
    // atan(2^-23)
    0x0000_0100_0000_0000,  // 1.1920928955078068e-07
    // atan(2^-24)
    0x0000_0080_0000_0000,  // 5.960464477539055e-08
    // atan(2^-25)
    0x0000_0040_0000_0000,  // 2.9802322387695303e-08
    // atan(2^-26)
    0x0000_0020_0000_0000,  // 1.4901161193847655e-08
    // atan(2^-27)
    0x0000_0010_0000_0000,  // 7.450580596923828e-09
    // atan(2^-28)
    0x0000_0008_0000_0000,  // 3.725290298461914e-09
    // atan(2^-29)
    0x0000_0004_0000_0000,  // 1.862645149230957e-09
    // atan(2^-30)
    0x0000_0002_0000_0000,  // 9.313225746154785e-10
    // atan(2^-31)
    0x0000_0001_0000_0000,  // 4.656612873077393e-10
    // Remaining entries approach zero - for very high precision types
    0x0000_0000_8000_0000,  // atan(2^-32)
    0x0000_0000_4000_0000,  // atan(2^-33)
    0x0000_0000_2000_0000,  // atan(2^-34)
    0x0000_0000_1000_0000,  // atan(2^-35)
    0x0000_0000_0800_0000,  // atan(2^-36)
    0x0000_0000_0400_0000,  // atan(2^-37)
    0x0000_0000_0200_0000,  // atan(2^-38)
    0x0000_0000_0100_0000,  // atan(2^-39)
    0x0000_0000_0080_0000,  // atan(2^-40)
    0x0000_0000_0040_0000,  // atan(2^-41)
    0x0000_0000_0020_0000,  // atan(2^-42)
    0x0000_0000_0010_0000,  // atan(2^-43)
    0x0000_0000_0008_0000,  // atan(2^-44)
    0x0000_0000_0004_0000,  // atan(2^-45)
    0x0000_0000_0002_0000,  // atan(2^-46)
    0x0000_0000_0001_0000,  // atan(2^-47)
    0x0000_0000_0000_8000,  // atan(2^-48)
    0x0000_0000_0000_4000,  // atan(2^-49)
    0x0000_0000_0000_2000,  // atan(2^-50)
    0x0000_0000_0000_1000,  // atan(2^-51)
    0x0000_0000_0000_0800,  // atan(2^-52)
    0x0000_0000_0000_0400,  // atan(2^-53)
    0x0000_0000_0000_0200,  // atan(2^-54)
    0x0000_0000_0000_0100,  // atan(2^-55)
    0x0000_0000_0000_0080,  // atan(2^-56)
    0x0000_0000_0000_0040,  // atan(2^-57)
    0x0000_0000_0000_0020,  // atan(2^-58)
    0x0000_0000_0000_0010,  // atan(2^-59)
    0x0000_0000_0000_0008,  // atan(2^-60)
    0x0000_0000_0000_0004,  // atan(2^-61)
    0x0000_0000_0000_0002,  // atan(2^-62)
    0x0000_0000_0000_0001,  // atan(2^-63)
];

/// The inverse CORDIC gain factor for circular mode (1/K).
///
/// After CORDIC iterations, results are scaled by K ≈ 1.6468.
/// To compensate, we pre-multiply the initial vector by 1/K ≈ 0.6073.
///
/// K = ∏(i=0 to ∞) sqrt(1 + 2^(-2i)) ≈ 1.6467602581210656
/// 1/K ≈ 0.6072529350088812561694
///
/// Stored as I1F63 bit pattern.
pub const CIRCULAR_GAIN_INV: i64 = 0x4DBA_76D4_21AF_2D34; // ≈ 0.6072529350088813
