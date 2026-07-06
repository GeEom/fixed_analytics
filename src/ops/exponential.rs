//! Exponential and logarithmic functions.

use crate::bounded::{NormalizedLnArg, OpenUnitInterval};
use crate::error::{Error, Result};
use crate::ops::hyperbolic::atanh_open;
use crate::traits::CordicNumber;

/// Right-shifts a non-negative value by `n`, rounding to nearest.
///
/// A plain `>>` truncates, which biases results low by up to a full ulp.
/// That matters for the final scaling of `exp`, where the true result may
/// be only a few ulps: rounding halves the worst case and removes the bias.
fn shr_rounded<T: CordicNumber>(x: T, n: u32) -> T {
    if n == 0 {
        return x;
    }
    let ulp = T::one() >> T::frac_bits();
    let shifted = x >> n;
    // Reconstruct the discarded low bits: rem ∈ [0, 2^n) ulps.
    let rem = x.saturating_sub(shifted << n);
    let half = ulp << (n - 1);
    if rem >= half {
        shifted.saturating_add(ulp)
    } else {
        shifted
    }
}

/// Exponential function (e^x).
///
/// # Saturation Behavior
///
/// This function saturates for extreme inputs rather than returning an error:
///
/// | Condition | Result | Example (I16F16) |
/// |-----------|--------|------------------|
/// | x > `ln(T::MAX)` | `T::MAX` | x > ~10.4 → 32767.99 |
/// | x < `ln(T::MIN_POSITIVE / 2)` | `T::ZERO` | x < ~-11.8 → 0 |
///
/// The zero threshold is where e^x rounds to nearest below half an ulp.
/// The exact thresholds depend on the type's range:
/// - **I16F16:** Saturates to MAX for x > ~10.4, to zero for x < ~-11.8
/// - **I32F32:** Saturates to MAX for x > ~21.5, to zero for x < ~-22.9
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
#[cfg_attr(feature = "verify-no-panic", no_panic::no_panic)]
pub fn exp<T: CordicNumber>(x: T) -> T {
    let zero = T::zero();
    let one = T::one();
    let ln2 = T::ln_2();

    // Handle special case
    if x == zero {
        return one;
    }

    // Argument reduction: exp(x) = 2^k * exp(r), where r ∈ [0, ln2).
    // Compute k = floor(x / ln2): start from the truncated quotient and
    // adjust once if the remainder is negative (the quotient error is
    // below one ulp, so a single adjustment suffices). Keeping r
    // non-negative gives exp(r) ∈ [1, 2) — a full significand ahead of
    // the final scaling shift.
    #[allow(clippy::cast_possible_wrap, reason = "total_bits bounded by type size")]
    let max_shift = (T::total_bits() - 1) as i32;
    let mut scale = x.div(ln2).to_i32();

    // Early exit for values that will saturate after scaling. This also
    // keeps `scale` small enough that scale * ln2 below cannot overflow.
    if scale > max_shift {
        return T::max_value();
    }
    if scale < -max_shift {
        return zero;
    }

    let mut r = x.saturating_sub(T::from_num(scale).saturating_mul(ln2));
    if r < zero {
        scale -= 1;
        r = r.saturating_add(ln2);
        if scale < -max_shift {
            return zero;
        }
    }

    // Factored Taylor: exp(r) = 1 + r*(1 + r/2*(1 + r/3*(1 + ... r/n)))
    let mut p = one;
    if T::frac_bits() >= 24 {
        // High precision: degree 12 Taylor
        // Truncation error: |r^13/13!| ≤ (ln2)^13/13! ≈ 3.4e-15
        p = one.saturating_add(r.div(T::from_num(12)).saturating_mul(p));
        p = one.saturating_add(r.div(T::from_num(11)).saturating_mul(p));
        p = one.saturating_add(r.div(T::from_num(10)).saturating_mul(p));
        p = one.saturating_add(r.div(T::from_num(9)).saturating_mul(p));
        p = one.saturating_add(r.div(T::from_num(8)).saturating_mul(p));
    }
    // Common terms (degree 7 base)
    // Low-precision truncation error: |r^8/8!| ≤ (ln2)^8/8! ≈ 8.9e-7
    p = one.saturating_add(r.div(T::from_num(7)).saturating_mul(p));
    p = one.saturating_add(r.div(T::from_num(6)).saturating_mul(p));
    p = one.saturating_add(r.div(T::from_num(5)).saturating_mul(p));
    p = one.saturating_add(r.div(T::from_num(4)).saturating_mul(p));
    p = one.saturating_add(r.div(T::from_num(3)).saturating_mul(p));
    p = one.saturating_add(r.div(T::from_num(2)).saturating_mul(p));
    let exp_r = one.saturating_add(r.saturating_mul(p));

    // Scale by 2^scale using bit shifts.
    // scale is already bounded to [-max_shift, max_shift] by the early exits above.
    #[allow(
        clippy::cast_sign_loss,
        reason = "sign of scale checked before each cast"
    )]
    if scale > 0 {
        let shift = scale as u32;
        // Detect overflow before shifting: if exp_r > MAX >> shift,
        // the left shift would wrap, so saturate to MAX instead.
        let headroom = T::max_value() >> shift;
        if exp_r > headroom {
            T::max_value()
        } else {
            exp_r << shift
        }
    } else {
        // scale ≤ 0: rounded right shift (shr_rounded(x, 0) is the identity)
        shr_rounded(exp_r, (-scale) as u32)
    }
}

/// Natural logarithm. Domain: `x > 0`.
///
/// # Errors
/// Returns `DomainError` if `x ≤ 0`.
#[must_use = "returns the natural logarithm result which should be handled"]
#[cfg_attr(feature = "verify-no-panic", no_panic::no_panic)]
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
        k_ln2 = k_ln2.saturating_add(ln2);
        i += 1;
    }

    // For small x (< 0.5), multiply by 2 repeatedly
    i = 0;
    while normalized < half && i < 128 {
        normalized = normalized.saturating_add(normalized);
        k_ln2 = k_ln2.saturating_sub(ln2);
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
    let ln_normalized = atanh_val.saturating_add(atanh_val); // 2 * atanh

    Ok(ln_normalized.saturating_add(k_ln2))
}

/// Base-2 logarithm. Domain: `x > 0`.
///
/// # Errors
/// Returns `DomainError` if `x ≤ 0`.
#[must_use = "returns the base-2 logarithm result which should be handled"]
#[cfg_attr(feature = "verify-no-panic", no_panic::no_panic)]
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
#[cfg_attr(feature = "verify-no-panic", no_panic::no_panic)]
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
/// | x < `log2(T::MIN_POSITIVE / 2)` | `T::ZERO` | x < ~-17 → 0 |
///
/// The zero threshold is where 2^x rounds to nearest below half an ulp.
/// The exact thresholds:
/// - **I16F16:** Saturates for x > ~15 or x < ~-17
/// - **I32F32:** Saturates for x > ~31 or x < ~-33
#[must_use]
#[cfg_attr(feature = "verify-no-panic", no_panic::no_panic)]
pub fn pow2<T: CordicNumber>(x: T) -> T {
    let ln_2 = T::ln_2();
    exp(x.saturating_mul(ln_2))
}
