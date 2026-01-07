//! Exponential and logarithmic functions.

use crate::error::{Error, Result};
use crate::ops::hyperbolic::sinh_cosh;
use crate::traits::CordicNumber;

/// Exponential (e^x). May overflow for large positive x.
#[must_use]
pub fn exp<T: CordicNumber>(x: T) -> T {
    let zero = T::zero();
    let one = T::one();
    let ln2 = T::ln_2();

    // Handle special case
    if x == zero {
        return one;
    }

    // For large |x|, use argument reduction: exp(x) = exp(x/2)²
    // Or better: exp(x) = 2^k * exp(r) where x = k*ln(2) + r
    let abs_x = x.abs();
    let threshold = ln2 + ln2; // About 1.386

    if abs_x > threshold {
        // Argument reduction using exp(x) = exp(x - ln2) * 2
        // Find k such that |x - k*ln2| < ln2
        let mut reduced = x;
        let mut scale_factor = one;

        let mut i = 0;
        if x.is_positive() {
            while reduced > ln2 && i < 128 {
                reduced -= ln2;
                scale_factor = scale_factor + scale_factor; // *= 2
                i += 1;
            }
        } else {
            while reduced < -ln2 && i < 128 {
                reduced += ln2;
                scale_factor = scale_factor >> 1; // /= 2
                i += 1;
            }
        }

        // Now compute exp(reduced) where |reduced| <= ln2
        let (sinh_r, cosh_r) = sinh_cosh(reduced);
        let exp_r = cosh_r.saturating_add(sinh_r);

        return scale_factor.saturating_mul(exp_r);
    }

    // For small x, use exp(x) = cosh(x) + sinh(x) directly
    let (sinh_x, cosh_x) = sinh_cosh(x);
    cosh_x.saturating_add(sinh_x)
}

/// Natural logarithm. Domain: `x > 0`.
///
/// # Errors
/// Returns `DomainError` if `x ≤ 0`.
///
/// # Panics
/// Panics if the internal atanh computation fails, which should never happen
/// as the normalized argument is always in the valid range (-1/3, 1/3).
#[must_use = "returns the natural logarithm result which should be handled"]
#[allow(clippy::missing_panics_doc)] // Panic only on internal invariant violation
pub fn ln<T: CordicNumber>(x: T) -> Result<T> {
    let zero = T::zero();
    let one = T::one();
    let two = T::two();

    if x <= zero {
        return Err(Error::domain("ln", "positive value"));
    }

    if x == one {
        return Ok(zero);
    }

    // For x very close to 1, the direct formula works well
    // ln(x) = 2 * atanh((x-1)/(x+1))

    // For x far from 1, use argument reduction:
    // ln(x) = ln(x * 2^(-k)) + k * ln(2)
    // where k is chosen so that x * 2^(-k) is close to 1

    let ln2 = T::ln_2();
    let mut normalized = x;
    let mut k_ln2 = zero;

    // Reduce to range [0.5, 2] for better convergence.
    let half = T::half();

    // For large x, divide by 2 repeatedly
    let mut i = 0;
    while normalized > two && i < 128 {
        normalized = normalized >> 1;
        k_ln2 += ln2;
        i += 1;
    }

    // For small x (< 0.5), multiply by 2 repeatedly
    i = 0;
    while normalized < half && i < 128 {
        normalized = normalized + normalized;
        k_ln2 -= ln2;
        i += 1;
    }

    // Now compute ln(normalized) where 0.5 <= normalized <= 2
    // Using ln(x) = 2 * atanh((x-1)/(x+1))
    let x_minus_1 = normalized - one;
    let x_plus_1 = normalized + one;
    let arg = x_minus_1.div(x_plus_1);

    // SAFETY: normalized ∈ [0.5, 2], so arg = (x-1)/(x+1) ∈ (-1/3, 1/3) ⊂ (-1, 1).
    // atanh cannot fail for this range.
    #[allow(clippy::expect_used)] // Invariant: normalized in [0.5, 2] guarantees valid atanh input
    let atanh_val = crate::ops::hyperbolic::atanh(arg)
        .expect("normalized in [0.5, 2] guarantees arg in (-1/3, 1/3)");
    let ln_normalized = atanh_val + atanh_val; // 2 * atanh

    Ok(ln_normalized + k_ln2)
}

/// Base-2 logarithm. Domain: `x > 0`.
///
/// # Errors
/// Returns `DomainError` if `x ≤ 0`.
#[must_use = "returns the base-2 logarithm result which should be handled"]
pub fn log2<T: CordicNumber>(x: T) -> Result<T> {
    let ln_x = ln(x)?;
    let ln_2 = T::ln_2();
    Ok(ln_x.div(ln_2))
}

/// Base-10 logarithm. Domain: `x > 0`.
///
/// # Errors
/// Returns `DomainError` if `x ≤ 0`.
#[must_use = "returns the base-10 logarithm result which should be handled"]
pub fn log10<T: CordicNumber>(x: T) -> Result<T> {
    let ln_x = ln(x)?;
    let ln_10 = T::ln_10();
    Ok(ln_x.div(ln_10))
}

/// Power of 2 (2^x).
#[must_use]
pub fn pow2<T: CordicNumber>(x: T) -> T {
    let ln_2 = T::ln_2();
    exp(x.saturating_mul(ln_2))
}
