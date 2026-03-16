//! Minimax (Chebyshev) polynomial coefficients for sin/cos evaluation.
//!
//! Generated via the Remez algorithm (`mpmath.chebyfit`). Each coefficient
//! approximates a reduced-form remainder,
//! factored so all values fit in I1F63 (magnitude < 1).
//!
//! Stored highest-degree first: `[cₙ, cₙ₋₁, …, c₀]`.
//! Evaluated via Horner's method:
//!
//! ```text
//! P(x) = cₙ·xⁿ + ··· + c₁·x + c₀ = c₀ + x·(c₁ + x·(c₂ + ··· + x·cₙ))
//! ```
//!
//! # Reconstruction
//!
//! | Function | Factored form | Variable | Domain |
//! |----------|---------------|----------|--------|
//! | sin | `x + x³·P(x²)` | u = x² | \[0, (π/4)²\] |
//! | cos | `1 + x²·Q(x²)` | u = x² | \[0, (π/4)²\] |

use crate::traits::CordicNumber;

/// Evaluate polynomial via Horner's method using precomputed I1F63 coefficients.
///
/// Stored highest-degree first: `[cₙ, cₙ₋₁, …, c₀]`.
/// Evaluates P(x) = cₙ·xⁿ + cₙ₋₁·xⁿ⁻¹ + ··· + c₀ via Horner's method.
#[inline]
pub fn horner<T: CordicNumber, const N: usize>(coeffs: &[i64; N], x: T) -> T {
    let mut iter = coeffs.iter();
    // First element is the highest-degree coefficient (N ≥ 3 for all tables).
    let mut result = T::from_i1f63(*iter.next().unwrap_or(&0));
    for &coeff in iter {
        result = T::from_i1f63(coeff).saturating_add(x.saturating_mul(result));
    }
    result
}

// =============================================================================
// Sin: P(u) = (sin(x) - x) / x³,  u = x²,  domain [0, (π/4)²]
// =============================================================================

/// Minimax coeffs for (sin(x)-x)/x³ in u=x² on [0,(π/4)²]. Low precision (I16F16-class).
#[rustfmt::skip]
pub const SIN_P_LO: [i64; 4] = [
     0x0000_16DB_E083_8A07,  // +2.725e-06
    -0x0006_804E_9E4E_E633,  // -1.984e-04
     0x0111_110D_EF2E_1C96,  // +8.333e-03
    -0x1555_5555_45E0_ABF9,  // -1.667e-01
];

/// Minimax coeffs for (sin(x)-x)/x³ in u=x² on [0,(π/4)²]. High precision (I32F32-class).
#[rustfmt::skip]
pub const SIN_P_HI: [i64; 7] = [
    -0x0000_0000_006A_C5F5,  // -7.587e-13
     0x0000_0000_5848_5FC3,  // +1.606e-10
    -0x0000_0035_CC8A_8259,  // -2.505e-08
     0x0000_171D_E3A5_4609,  // +2.756e-06
    -0x0006_8068_0680_664E,  // -1.984e-04
     0x0111_1111_1111_1100,  // +8.333e-03
    -0x1555_5555_5555_5555,  // -1.667e-01
];

// =============================================================================
// Cos: Q(u) = (cos(x) - 1) / x²,  u = x²,  domain [0, (π/4)²]
// =============================================================================

/// Minimax coeffs for (cos(x)-1)/x² in u=x² on [0,(π/4)²]. Low precision.
#[rustfmt::skip]
pub const COS_Q_LO: [i64; 4] = [
     0x0000_CD37_95D8_0CFA,  // +2.446e-05
    -0x002D_81C1_0FEA_FDD5,  // -1.389e-03
     0x0555_5532_ED10_3C5F,  // +4.167e-02
    -0x3FFF_FFFF_563B_4B19,  // -5.000e-01
];

/// Minimax coeffs for (cos(x)-1)/x² in u=x² on [0,(π/4)²]. High precision.
#[rustfmt::skip]
pub const COS_Q_HI: [i64; 7] = [
    -0x0000_0000_063F_E751,  // -1.137e-11
     0x0000_0004_7BA9_FC96,  // +2.088e-09
    -0x0000_024F_C9F1_C946,  // -2.756e-07
     0x0000_D00D_00CE_F099,  // +2.480e-05
    -0x002D_82D8_2D82_BAF2,  // -1.389e-03
     0x0555_5555_5555_5435,  // +4.167e-02
    -0x3FFF_FFFF_FFFF_FFFE,  // -5.000e-01
];
