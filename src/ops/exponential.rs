//! Exponential and logarithmic functions.

use crate::bounded::{NormalizedLnArg, OpenUnitInterval};
use crate::error::{Error, Result};
use crate::ops::hyperbolic::{atanh_open, sinh_cosh};
use crate::traits::CordicNumber;

/// Exponential function (e^x).
///
/// # Saturation Behavior
///
/// This function saturates for extreme inputs rather than returning an error:
///
/// | Condition | Result | Example (I16F16) |
/// |-----------|--------|------------------|
/// | x > `ln(T::MAX)` | `T::MAX` | x > ~10.4 → 32767.99 |
/// | x < `ln(T::MIN_POSITIVE)` | `T::ZERO` | x < ~-10.4 → 0 |
///
/// The exact thresholds depend on the type's range:
/// - **I16F16:** Saturates for x > ~10.4 or x < ~-20
/// - **I32F32:** Saturates for x > ~21.5 or x < ~-45
///
/// Saturation is silent and deterministic. If you need to detect overflow,
/// check the input range before calling:
///
/// ```
/// use fixed::types::I16F16;
/// use fixed_analytics::exp;
///
/// let x = I16F16::from_num(5.0);
/// let max_safe = I16F16::from_num(10.0);
///
/// if x < max_safe {
///     let result = exp(x);  // Safe
/// } else {
///     // Handle potential saturation
/// }
/// ```
#[must_use]
pub fn exp<T: CordicNumber>(x: T) -> T {
    let zero = T::zero();
    let one = T::one();
    let ln2 = T::ln_2();

    // Handle special case
    if x == zero {
        return one;
    }

    // For large |x|, use argument reduction: exp(x) = 2^k * exp(r)
    // where r is reduced to (-ln2, ln2) range
    let mut reduced = x;
    let mut scale: i32 = 0;

    // Reduce positive values
    let mut i = 0;
    while reduced > ln2 && i < 64 {
        reduced -= ln2;
        scale += 1;
        i += 1;
    }

    // Reduce negative values
    i = 0;
    while reduced < -ln2 && i < 64 {
        reduced += ln2;
        scale -= 1;
        i += 1;
    }

    // Compute exp(reduced) = cosh(reduced) + sinh(reduced)
    let (sinh_r, cosh_r) = sinh_cosh(reduced);
    let exp_r = cosh_r.saturating_add(sinh_r);

    // Scale by 2^scale using bit shifts
    #[allow(clippy::cast_possible_wrap, reason = "total_bits bounded by type size")]
    let max_shift = (T::total_bits() - 1) as i32;

    #[allow(clippy::cast_sign_loss, reason = "scale >= 0 checked before cast")]
    if scale >= 0 {
        if scale > max_shift {
            T::max_value()
        } else {
            let shift = scale as u32;
            exp_r << shift
        }
    } else {
        let neg_scale = -scale;
        if neg_scale > max_shift {
            zero
        } else {
            let shift = neg_scale as u32;
            exp_r >> shift
        }
    }
}

/// Natural logarithm. Domain: `x > 0`.
///
/// # Errors
/// Returns `DomainError` if `x ≤ 0`.
#[must_use = "returns the natural logarithm result which should be handled"]
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
    // NormalizedLnArg encodes that normalized ∈ [0.5, 2]
    let norm = NormalizedLnArg::from_normalized(normalized);

    // OpenUnitInterval::from_normalized_ln_arg computes (x-1)/(x+1),
    // which is in (-1/3, 1/3) ⊂ (-1, 1) for x ∈ [0.5, 2]
    let arg = OpenUnitInterval::from_normalized_ln_arg(norm);

    let atanh_val = atanh_open(arg);
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

/// Power of 2 (2^x). Computed as exp(x × ln(2)).
///
/// # Saturation Behavior
///
/// Saturates for extreme inputs:
///
/// | Condition | Result | Example (I16F16) |
/// |-----------|--------|------------------|
/// | x > `log2(T::MAX)` | `T::MAX` | x > ~15 → 32767.99 |
/// | x < `log2(T::MIN_POSITIVE)` | `T::ZERO` | x < ~-16 → 0 |
///
/// The exact thresholds:
/// - **I16F16:** Saturates for x > ~15 or x < ~-16
/// - **I32F32:** Saturates for x > ~31 or x < ~-32
#[must_use]
pub fn pow2<T: CordicNumber>(x: T) -> T {
    let ln_2 = T::ln_2();
    exp(x.saturating_mul(ln_2))
}
