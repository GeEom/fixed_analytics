//! Exponential and logarithmic functions.
//!
//! Provides exp, ln, log2, and log10 using CORDIC and related algorithms.
//!
//! # Implementation Notes
//!
//! - `exp(x)` uses the identity `exp(x) = cosh(x) + sinh(x)` with argument reduction
//! - `ln(x)` uses the identity `ln(x) = 2 * atanh((x-1)/(x+1))`
//! - `log2` and `log10` are derived from `ln` using change of base

use crate::error::{Error, Result};
use crate::ops::hyperbolic::sinh_cosh;
use crate::traits::CordicNumber;

/// Computes e^x (the exponential function).
///
/// Uses the identity `exp(x) = cosh(x) + sinh(x)` with argument reduction
/// for large values.
///
/// # Arguments
///
/// * `x` - The exponent
///
/// # Returns
///
/// e raised to the power x.
///
/// # Examples
///
/// ```
/// use fixed::types::I16F16;
/// use fixed_analytics::exp;
///
/// let x = I16F16::ZERO;
/// let result = exp(x);
/// // result ≈ 1.0
///
/// let x = I16F16::from_num(1.0);
/// let result = exp(x);
/// // result ≈ 2.718
/// ```
///
/// # Note
///
/// May overflow for large positive values of x. The exact overflow threshold
/// depends on the fixed-point format used.
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
        // Guard against infinite loops (limit to 128 iterations, enough for any practical value).
        let mut reduced = x;
        let mut scale_factor = one;
        let max_iterations = 128_u32;
        let mut iterations = 0_u32;

        if x.is_positive() {
            while reduced > ln2 && iterations < max_iterations {
                reduced -= ln2;
                scale_factor = scale_factor + scale_factor; // *= 2
                iterations += 1;
            }
        } else {
            while reduced < -ln2 && iterations < max_iterations {
                reduced += ln2;
                scale_factor = scale_factor >> 1; // /= 2
                iterations += 1;
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

/// Computes the natural logarithm (base e).
///
/// Uses the identity `ln(x) = 2 * atanh((x-1)/(x+1))` which is suitable
/// for CORDIC computation.
///
/// # Arguments
///
/// * `x` - A positive value
///
/// # Returns
///
/// The natural logarithm of x.
///
/// # Errors
///
/// Returns [`Error::DomainError`] if `x <= 0`.
///
/// # Examples
///
/// ```
/// use fixed::types::I16F16;
/// use fixed_analytics::ln;
///
/// let x = I16F16::from_num(1.0);
/// let result = ln(x).unwrap();
/// // result ≈ 0.0
///
/// let x = I16F16::E;
/// let result = ln(x).unwrap();
/// // result ≈ 1.0
/// ```
#[must_use = "returns the natural logarithm result which should be handled"]
pub fn ln<T: CordicNumber>(x: T) -> Result<T> {
    let zero = T::zero();
    let one = T::one();
    let two = T::two();

    if x <= zero {
        return Err(Error::DomainError {
            function: "ln",
            expected: "positive value",
        });
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
    // Guard against infinite loops (limit to 128 iterations, enough for any practical value).
    let half = T::half();
    let max_iterations = 128_u32;
    let mut iterations = 0_u32;

    // For large x, divide by 2 repeatedly
    while normalized > two && iterations < max_iterations {
        normalized = normalized >> 1;
        k_ln2 += ln2;
        iterations += 1;
    }

    // For small x (< 0.5), multiply by 2 repeatedly
    while normalized < half && iterations < max_iterations {
        normalized = normalized + normalized;
        k_ln2 -= ln2;
        iterations += 1;
    }

    // Now compute ln(normalized) where 0.5 <= normalized <= 2
    // Using ln(x) = 2 * atanh((x-1)/(x+1))
    let x_minus_1 = normalized - one;
    let x_plus_1 = normalized + one;
    let arg = x_minus_1.div(x_plus_1);

    // atanh is computed via CORDIC.
    // Since normalized is in [0.5, 2], arg = (x-1)/(x+1) is in (-1/3, 1/3) ⊂ (-1, 1),
    // so atanh will always succeed. The unwrap_or(zero) is defensive only.
    let atanh_val = crate::ops::hyperbolic::atanh(arg).unwrap_or(zero);
    let ln_normalized = atanh_val + atanh_val; // 2 * atanh

    Ok(ln_normalized + k_ln2)
}

/// Computes the base-2 logarithm.
///
/// `log2(x) = ln(x) / ln(2)`
///
/// # Arguments
///
/// * `x` - A positive value
///
/// # Returns
///
/// The base-2 logarithm of x.
///
/// # Errors
///
/// Returns [`Error::DomainError`] if `x <= 0`.
///
/// # Examples
///
/// ```
/// use fixed::types::I16F16;
/// use fixed_analytics::log2;
///
/// let x = I16F16::from_num(8.0);
/// let result = log2(x).unwrap();
/// // result ≈ 3.0
/// ```
#[must_use = "returns the base-2 logarithm result which should be handled"]
pub fn log2<T: CordicNumber>(x: T) -> Result<T> {
    let ln_x = ln(x)?;
    let ln_2 = T::ln_2();
    Ok(ln_x.div(ln_2))
}

/// Computes the base-10 logarithm.
///
/// `log10(x) = ln(x) / ln(10)`
///
/// # Arguments
///
/// * `x` - A positive value
///
/// # Returns
///
/// The base-10 logarithm of x.
///
/// # Errors
///
/// Returns [`Error::DomainError`] if `x <= 0`.
///
/// # Examples
///
/// ```
/// use fixed::types::I16F16;
/// use fixed_analytics::log10;
///
/// let x = I16F16::from_num(100.0);
/// let result = log10(x).unwrap();
/// // result ≈ 2.0
/// ```
#[must_use = "returns the base-10 logarithm result which should be handled"]
pub fn log10<T: CordicNumber>(x: T) -> Result<T> {
    let ln_x = ln(x)?;
    let ln_10 = T::ln_10();
    Ok(ln_x.div(ln_10))
}

/// Computes 2^x (power of 2).
///
/// `pow2(x) = exp(x * ln(2))`
///
/// # Arguments
///
/// * `x` - The exponent
///
/// # Returns
///
/// 2 raised to the power x.
///
/// # Examples
///
/// ```
/// use fixed::types::I16F16;
/// use fixed_analytics::ops::exponential::pow2;
///
/// let x = I16F16::from_num(3.0);
/// let result = pow2(x);
/// // result ≈ 8.0
/// ```
#[must_use]
pub fn pow2<T: CordicNumber>(x: T) -> T {
    let ln_2 = T::ln_2();
    exp(x.saturating_mul(ln_2))
}
