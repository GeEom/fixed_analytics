//! Hyperbolic functions.
//!
//! Provides sinh, cosh, tanh, coth, asinh, acosh, atanh, acoth using hyperbolic CORDIC.
//!
//! # Mathematical Definitions
//!
//! - `sinh(x) = (e^x - e^(-x)) / 2`
//! - `cosh(x) = (e^x + e^(-x)) / 2`
//! - `tanh(x) = sinh(x) / cosh(x)`
//!
//! # Precision
//!
//! Hyperbolic CORDIC has slightly lower precision than circular CORDIC due to
//! the need for iteration repetition. Expect ~3-4 decimal digits for `I16F16`.

use crate::error::{Error, Result};
use crate::kernel::{hyperbolic_gain_inv, hyperbolic_rotation, hyperbolic_vectoring};
use crate::traits::CordicNumber;

/// Fractional part of hyperbolic convergence limit (~0.1182).
/// The full limit is approximately 1.1182, computed as 1 + this fraction.
/// We store only the fractional part since I1F63 cannot represent values >= 1.
const HYPERBOLIC_LIMIT_FRAC: i64 = 0x0F22_3D70_A3D7_0A3D; // ~0.1182

/// Computes the hyperbolic sine and cosine simultaneously.
///
/// This is more efficient than calling [`sinh`] and [`cosh`] separately
/// since both values are produced by a single CORDIC iteration.
///
/// # Arguments
///
/// * `x` - The input value
///
/// # Returns
///
/// A tuple `(sinh(x), cosh(x))`.
///
/// # Examples
///
/// ```
/// use fixed::types::I16F16;
/// use fixed_analytics::sinh_cosh;
///
/// let x = I16F16::from_num(1.0);
/// let (s, c) = sinh_cosh(x);
/// // s ≈ 1.175, c ≈ 1.543
/// ```
#[must_use]
pub fn sinh_cosh<T: CordicNumber>(x: T) -> (T, T) {
    let zero = T::zero();
    let one = T::one();
    // Compute limit as 1 + fractional_part (~1.1182)
    let limit = one.saturating_add(T::from_i64_frac(HYPERBOLIC_LIMIT_FRAC));

    // Handle argument reduction for large values
    if x.abs() > limit {
        // Use the identities:
        // sinh(2x) = 2 * sinh(x) * cosh(x)
        // cosh(2x) = cosh²(x) + sinh²(x) = 2*cosh²(x) - 1
        let half_x = x >> 1;
        let (sh, ch) = sinh_cosh(half_x);

        let sinh_result = sh.saturating_mul(ch).saturating_mul(T::two());
        let cosh_result = ch.saturating_mul(ch).saturating_add(sh.saturating_mul(sh));

        return (sinh_result, cosh_result);
    }

    // For very small x, use Taylor series approximation to avoid CORDIC
    // overshoot on the first iteration (where atanh(0.5) ≈ 0.549 is larger than x).
    // sinh(x) ≈ x + x³/6 ≈ x for small x
    // cosh(x) ≈ 1 + x²/2 for small x
    let small_threshold = T::from_num(0.1); // Below first atanh table entry / 5
    if x.abs() < small_threshold {
        let x_sq = x.saturating_mul(x);
        // sinh(x) ≈ x (higher order terms negligible for |x| < 0.1)
        let sinh_approx = x;
        // cosh(x) ≈ 1 + x²/2
        let cosh_approx = one.saturating_add(x_sq >> 1);
        return (sinh_approx, cosh_approx);
    }

    // For moderate x, use CORDIC directly.
    // Hyperbolic CORDIC scales results by 1/K_h ≈ 1.2075.
    // To compensate, we pre-multiply by 1/K_h (using precomputed constant).
    let inv_gain = hyperbolic_gain_inv(); // 1/K_h ≈ 1.2075

    let (cosh_val, sinh_val, _) = hyperbolic_rotation(inv_gain, zero, x);

    (sinh_val, cosh_val)
}

/// Computes the hyperbolic sine.
///
/// `sinh(x) = (e^x - e^(-x)) / 2`
///
/// # Arguments
///
/// * `x` - The input value
///
/// # Returns
///
/// The hyperbolic sine of `x`.
///
/// # Examples
///
/// ```
/// use fixed::types::I16F16;
/// use fixed_analytics::sinh;
///
/// let x = I16F16::ZERO;
/// let result = sinh(x);
/// // result ≈ 0.0
/// ```
#[inline]
#[must_use]
pub fn sinh<T: CordicNumber>(x: T) -> T {
    sinh_cosh(x).0
}

/// Computes the hyperbolic cosine.
///
/// `cosh(x) = (e^x + e^(-x)) / 2`
///
/// # Arguments
///
/// * `x` - The input value
///
/// # Returns
///
/// The hyperbolic cosine of `x` (always >= 1).
///
/// # Examples
///
/// ```
/// use fixed::types::I16F16;
/// use fixed_analytics::cosh;
///
/// let x = I16F16::ZERO;
/// let result = cosh(x);
/// // result ≈ 1.0
/// ```
#[inline]
#[must_use]
pub fn cosh<T: CordicNumber>(x: T) -> T {
    sinh_cosh(x).1
}

/// Computes the hyperbolic tangent.
///
/// `tanh(x) = sinh(x) / cosh(x)`
///
/// # Arguments
///
/// * `x` - The input value
///
/// # Returns
///
/// The hyperbolic tangent of `x`, in the range (-1, 1).
///
/// # Examples
///
/// ```
/// use fixed::types::I16F16;
/// use fixed_analytics::tanh;
///
/// let x = I16F16::ZERO;
/// let result = tanh(x);
/// // result ≈ 0.0
/// ```
#[must_use]
pub fn tanh<T: CordicNumber>(x: T) -> T {
    let (s, c) = sinh_cosh(x);
    s.div(c)
}

/// Computes the hyperbolic cotangent.
///
/// `coth(x) = cosh(x) / sinh(x)`
///
/// # Arguments
///
/// * `x` - The input value (must be non-zero)
///
/// # Returns
///
/// The hyperbolic cotangent of `x`. Returns `T::MAX` for `x = 0` (mathematically undefined).
///
/// # Note
///
/// The function is mathematically undefined at `x = 0` where it has a pole.
/// For very small `|x|`, the result will be large and may saturate.
#[must_use]
pub fn coth<T: CordicNumber>(x: T) -> T {
    let zero = T::zero();
    if x == zero {
        // coth(0) is undefined (pole), return maximum value
        return T::max_value();
    }
    let (s, c) = sinh_cosh(x);
    c.div(s)
}

/// Computes the inverse hyperbolic sine.
///
/// `asinh(x) = ln(x + sqrt(x² + 1))`
///
/// # Arguments
///
/// * `x` - Any real value
///
/// # Returns
///
/// The inverse hyperbolic sine of `x`.
///
/// # Examples
///
/// ```
/// use fixed::types::I16F16;
/// use fixed_analytics::asinh;
///
/// let x = I16F16::ZERO;
/// let result = asinh(x);
/// // result ≈ 0.0
/// ```
#[must_use]
pub fn asinh<T: CordicNumber>(x: T) -> T {
    let zero = T::zero();
    let one = T::one();

    if x == zero {
        return zero;
    }

    // asinh(x) = sign(x) * ln(|x| + sqrt(x² + 1))
    // For CORDIC, we use: asinh(x) = atanh(x / sqrt(1 + x²))
    let x_sq = x.saturating_mul(x);
    let one_plus_x_sq = one.saturating_add(x_sq);
    let sqrt_term = crate::ops::algebraic::sqrt(one_plus_x_sq);

    // Compute x / sqrt(1 + x²), which is in (-1, 1)
    let arg = x.div(sqrt_term);

    atanh_inner(arg)
}

/// Computes the inverse hyperbolic cosine.
///
/// `acosh(x) = ln(x + sqrt(x² - 1))`
///
/// # Arguments
///
/// * `x` - A value >= 1
///
/// # Returns
///
/// The inverse hyperbolic cosine of `x` (always >= 0).
///
/// # Errors
///
/// Returns [`Error::DomainError`] if `x < 1`.
///
/// # Examples
///
/// ```
/// use fixed::types::I16F16;
/// use fixed_analytics::acosh;
///
/// let x = I16F16::from_num(1.0);
/// let result = acosh(x).unwrap();
/// // result ≈ 0.0
/// ```
#[must_use = "returns the inverse hyperbolic cosine result which should be handled"]
pub fn acosh<T: CordicNumber>(x: T) -> Result<T> {
    let one = T::one();

    if x < one {
        return Err(Error::DomainError {
            function: "acosh",
            expected: "value >= 1",
        });
    }

    if x == one {
        return Ok(T::zero());
    }

    // acosh(x) = ln(x + sqrt(x² - 1))
    // Using CORDIC: acosh(x) = atanh(sqrt(x² - 1) / x) for x > 0
    // But this requires |sqrt(x²-1)/x| < 1, which is true for x > 1
    let x_sq = x.saturating_mul(x);
    let x_sq_minus_one = x_sq.saturating_sub(one);
    let sqrt_term = crate::ops::algebraic::sqrt(x_sq_minus_one);

    let arg = sqrt_term.div(x);
    Ok(atanh_inner(arg))
}

/// Computes the inverse hyperbolic tangent.
///
/// `atanh(x) = 0.5 * ln((1 + x) / (1 - x))`
///
/// # Arguments
///
/// * `x` - A value in the range (-1, 1)
///
/// # Returns
///
/// The inverse hyperbolic tangent of `x`.
///
/// # Errors
///
/// Returns [`Error::DomainError`] if `|x| >= 1`.
///
/// # Examples
///
/// ```
/// use fixed::types::I16F16;
/// use fixed_analytics::atanh;
///
/// let x = I16F16::from_num(0.5);
/// let result = atanh(x).unwrap();
/// // result ≈ 0.549
/// ```
#[must_use = "returns the inverse hyperbolic tangent result which should be handled"]
pub fn atanh<T: CordicNumber>(x: T) -> Result<T> {
    let one = T::one();

    if x >= one || x <= -one {
        return Err(Error::DomainError {
            function: "atanh",
            expected: "value in range (-1, 1)",
        });
    }

    Ok(atanh_inner(x))
}

/// Inner atanh implementation without bounds checking.
///
/// Uses hyperbolic CORDIC vectoring mode with argument reduction for large inputs.
/// The CORDIC convergence limit is about 1.1182, meaning we can only compute
/// atanh(x) directly when |x| < tanh(1.1182) ≈ 0.807.
fn atanh_inner<T: CordicNumber>(x: T) -> T {
    let zero = T::zero();
    let one = T::one();

    if x == zero {
        return zero;
    }

    // Threshold for argument reduction: tanh(1.0) ≈ 0.762
    // Using 0.75 to stay safely within convergence
    let threshold = T::from_num(0.75);

    if x.abs() <= threshold {
        // Direct CORDIC computation
        let (_, _, z) = hyperbolic_vectoring(one, x, zero);
        return z;
    }

    // Argument reduction using the identity:
    // atanh(x) = atanh(a) + atanh((x - a) / (1 - a*x))
    // We use a = 0.5, for which atanh(0.5) is a precomputed constant.
    let half = T::half();
    let atanh_half = T::from_i64_frac(crate::tables::hyperbolic::ATANH_HALF);

    let sign = if x.is_negative() { -one } else { one };
    let abs_x = x.abs();

    // Compute reduced argument: (|x| - 0.5) / (1 - 0.5*|x|)
    let numerator = abs_x - half;
    let denominator = one - half.saturating_mul(abs_x);
    let reduced = numerator.div(denominator);

    // Recursively compute atanh of reduced argument
    let atanh_reduced = atanh_inner(reduced);

    // atanh(x) = sign * (atanh(0.5) + atanh(reduced))
    sign.saturating_mul(atanh_half.saturating_add(atanh_reduced))
}

/// Computes the inverse hyperbolic cotangent.
///
/// `acoth(x) = 0.5 * ln((x + 1) / (x - 1))`
///
/// # Arguments
///
/// * `x` - A value with `|x| > 1`
///
/// # Returns
///
/// The inverse hyperbolic cotangent of `x`.
///
/// # Errors
///
/// Returns [`Error::DomainError`] if `|x| <= 1`.
#[must_use = "returns the inverse hyperbolic cotangent result which should be handled"]
pub fn acoth<T: CordicNumber>(x: T) -> Result<T> {
    let one = T::one();

    if x.abs() <= one {
        return Err(Error::DomainError {
            function: "acoth",
            expected: "|value| > 1",
        });
    }

    // acoth(x) = atanh(1/x)
    let recip = one.div(x);
    Ok(atanh_inner(recip))
}
