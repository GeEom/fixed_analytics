//! Hyperbolic arctangent table. Values: `atanh(2^-i)` as I1F63.
//!
//! Starts at i=1 (atanh(1) undefined). Requires repeats at 4, 13, 40, ...
//! For i ≥ 21, atanh(x) ≈ x to I1F63 precision.

/// `atanh(2^-i)` as I1F63. Index 0 = atanh(2^-1) = atanh(0.5).
#[rustfmt::skip]
pub const ATANH_TABLE: [i64; 64] = [
    0x464F_A9EA_B40C_2A5E,  // atanh(2^-1) = 0.549
    0x20B1_5DF5_0228_A34E,  // atanh(2^-2)
    0x1015_891C_9EAE_F76A,  // atanh(2^-3)
    0x0802_AC45_69BA_D66E,  // atanh(2^-4) - repeat this
    0x0400_5562_246B_B893,  // atanh(2^-5)
    0x0200_0AAB_1115_A393,  // atanh(2^-6)
    0x0100_0155_5888_91AD,  // atanh(2^-7)
    0x0080_002A_AAC4_4457,  // atanh(2^-8)
    0x0040_0005_5556_2222,  // atanh(2^-9)
    0x0020_0000_AAAA_B111,  // atanh(2^-10)
    0x0010_0000_1555_5589,  // atanh(2^-11)
    0x0008_0000_02AA_AAAC,  // atanh(2^-12)
    0x0004_0000_0055_5555,  // atanh(2^-13) - repeat this
    0x0002_0000_000A_AAAB,  // atanh(2^-14)
    0x0001_0000_0001_5555,  // atanh(2^-15)
    0x0000_8000_0000_2AAB,  // atanh(2^-16)
    0x0000_4000_0000_0555,  // atanh(2^-17)
    0x0000_2000_0000_00AB,  // atanh(2^-18)
    0x0000_1000_0000_0015,  // atanh(2^-19)
    0x0000_0800_0000_0003,  // atanh(2^-20)
    // For i >= 21: atanh(2^-i) ≈ 2^-i to I1F63 precision
    0x0000_0400_0000_0000,
    0x0000_0200_0000_0000,
    0x0000_0100_0000_0000,
    0x0000_0080_0000_0000,
    0x0000_0040_0000_0000,
    0x0000_0020_0000_0000,
    0x0000_0010_0000_0000,
    0x0000_0008_0000_0000,
    0x0000_0004_0000_0000,
    0x0000_0002_0000_0000,
    0x0000_0001_0000_0000,
    0x0000_0000_8000_0000,
    0x0000_0000_4000_0000,
    0x0000_0000_2000_0000,
    0x0000_0000_1000_0000,
    0x0000_0000_0800_0000,
    0x0000_0000_0400_0000,
    0x0000_0000_0200_0000,
    0x0000_0000_0100_0000,
    0x0000_0000_0080_0000,  // atanh(2^-40) - repeat this
    0x0000_0000_0040_0000,
    0x0000_0000_0020_0000,
    0x0000_0000_0010_0000,
    0x0000_0000_0008_0000,
    0x0000_0000_0004_0000,
    0x0000_0000_0002_0000,
    0x0000_0000_0001_0000,
    0x0000_0000_0000_8000,
    0x0000_0000_0000_4000,
    0x0000_0000_0000_2000,
    0x0000_0000_0000_1000,
    0x0000_0000_0000_0800,
    0x0000_0000_0000_0400,
    0x0000_0000_0000_0200,
    0x0000_0000_0000_0100,
    0x0000_0000_0000_0080,
    0x0000_0000_0000_0040,
    0x0000_0000_0000_0020,
    0x0000_0000_0000_0010,
    0x0000_0000_0000_0008,
    0x0000_0000_0000_0004,
    0x0000_0000_0000_0002,
    0x0000_0000_0000_0001,
    0x0000_0000_0000_0001,  // rounds to 1 LSB
];

/// `K_h` ≈ 0.828 (I1F63). Hyperbolic gain factor.
pub const HYPERBOLIC_GAIN: i64 = 0x6A01_203D_99A6_3986;

/// Returns true if iteration `i` must be repeated for hyperbolic CORDIC convergence.
///
/// The repeat sequence is 4, 13, 40, 121, 364, ... (each term is 3×previous + 1).
#[inline]
#[must_use]
pub const fn needs_repeat(index: u32) -> bool {
    matches!(index, 4 | 13 | 40 | 121 | 364)
}

/// `1/K_h` ≈ 1.2075 (I2F62). Pre-multiply to compensate for hyperbolic gain.
pub const HYPERBOLIC_GAIN_INV: i64 = 0x4D47_A1C8_03BB_08CA;

/// atanh(0.5) ≈ 0.549 (I1F63). Used for argument reduction.
pub const ATANH_HALF: i64 = 0x464F_A9EA_B40C_2A5E;
