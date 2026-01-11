//! Hyperbolic functions via hyperbolic CORDIC.

use crate::bounded::{AtLeastOne, NonNegative, OpenUnitInterval};
use crate::error::{Error, Result};
use crate::kernel::{hyperbolic_gain_inv, hyperbolic_rotation, hyperbolic_vectoring};
use crate::ops::algebraic::sqrt_nonneg;
use crate::traits::CordicNumber;

/// Hyperbolic CORDIC converges for |x| < sum of atanh table ≈ 1.1182.
/// Stored as fractional part (0.1182) since I1F63 can't hold 1.x.
/// Full limit = 1 + this value.
const HYPERBOLIC_CONVERGENCE_LIMIT_FRAC: i64 = 0x0F22_3D70_A3D7_0A3D;

/// Taylor series threshold for high precision (≥24 frac bits): 0.05 in I1F63 format.
/// Lower threshold means less Taylor usage, better for high precision where CORDIC excels.
const TAYLOR_THRESHOLD_HIGH_PREC_I1F63: i64 = 0x0666_6666_6666_6666; // 0.05 * 2^63

/// Taylor series threshold for low precision (<24 frac bits): 0.1 in I1F63 format.
/// Higher threshold uses Taylor more, which is sufficient for lower precision types.
const TAYLOR_THRESHOLD_LOW_PREC_I1F63: i64 = 0x0CCC_CCCC_CCCC_CCCD; // 0.1 * 2^63

/// atanh argument reduction threshold: 0.75 in I1F63 format.
/// tanh(1.0) ≈ 0.762; use 0.75 to stay within CORDIC convergence with margin.
const ATANH_REDUCTION_THRESHOLD_I1F63: i64 = 0x6000_0000_0000_0000; // 0.75 * 2^63

/// Hyperbolic sine and cosine. More efficient than separate calls.
#[must_use]
pub fn sinh_cosh<T: CordicNumber>(x: T) -> (T, T) {
    let zero = T::zero();
    let one = T::one();
    // Compute limit as 1 + fractional_part (~1.1182)
    let limit = one.saturating_add(T::from_i1f63(HYPERBOLIC_CONVERGENCE_LIMIT_FRAC));

    // Handle argument reduction for large values
    if x.abs() > limit {
        // Use the identities:
        // sinh(2x) = 2 * sinh(x) * cosh(x)
        // cosh(2x) = cosh²(x) + sinh²(x)
        let half_x = x >> 1;
        let (sh, ch) = sinh_cosh(half_x);

        let sinh_result = sh.saturating_mul(ch).saturating_mul(T::two());
        let cosh_result = ch.saturating_mul(ch).saturating_add(sh.saturating_mul(sh));

        return (sinh_result, cosh_result);
    }

    // For very small x, use Taylor series approximation to avoid CORDIC
    // overshoot on the first iteration (where atanh(0.5) ≈ 0.549 is larger than x).
    // Use precision-dependent threshold and order.
    let threshold_bits = if T::frac_bits() >= 24 {
        TAYLOR_THRESHOLD_HIGH_PREC_I1F63
    } else {
        TAYLOR_THRESHOLD_LOW_PREC_I1F63
    };
    let small_threshold = T::from_i1f63(threshold_bits);
    if x.abs() < small_threshold {
        let x_sq = x.saturating_mul(x);
        let x_cu = x_sq.saturating_mul(x);
        let x_qu = x_sq.saturating_mul(x_sq);

        // Higher precision benefits from higher-order Taylor terms
        if T::frac_bits() >= 24 {
            // sinh(x) ≈ x + x³/6 + x⁵/120, cosh(x) ≈ 1 + x²/2 + x⁴/24 + x⁶/720
            let x_5 = x_qu.saturating_mul(x);
            let x_6 = x_qu.saturating_mul(x_sq);
            // cosh base: 1 + x²/2 + x⁴/24
            let cosh_base = one
                .saturating_add(x_sq >> 1)
                .saturating_add(x_qu.div(T::from_num(24)));
            let cosh_approx = cosh_base.saturating_add(x_6.div(T::from_num(720)));
            // sinh base: x + x³/6
            let sinh_base = x.saturating_add(x_cu.div(T::from_num(6)));
            let sinh_approx = sinh_base.saturating_add(x_5.div(T::from_num(120)));
            return (sinh_approx, cosh_approx);
        }
        // sinh(x) ≈ x + x³/6, cosh(x) ≈ 1 + x²/2 + x⁴/24
        let c = one.saturating_add(x_sq >> 1);
        let cosh_approx = c.saturating_add(x_qu.div(T::from_num(24)));
        let sinh_approx = x.saturating_add(x_cu.div(T::from_num(6)));
        return (sinh_approx, cosh_approx);
    }

    // For moderate x, use CORDIC directly.
    // Hyperbolic CORDIC scales results by 1/K_h ≈ 1.2075.
    // To compensate, we pre-multiply by 1/K_h (using precomputed constant).
    let inv_gain = hyperbolic_gain_inv(); // 1/K_h ≈ 1.2075

    let (cosh_val, sinh_val, _) = hyperbolic_rotation(inv_gain, zero, x);

    (sinh_val, cosh_val)
}

/// Hyperbolic sine.
#[inline]
#[must_use]
pub fn sinh<T: CordicNumber>(x: T) -> T {
    sinh_cosh(x).0
}

/// Hyperbolic cosine. Always ≥1.
#[inline]
#[must_use]
pub fn cosh<T: CordicNumber>(x: T) -> T {
    sinh_cosh(x).1
}

/// Hyperbolic tangent. Result in `(-1, 1)`.
#[must_use]
pub fn tanh<T: CordicNumber>(x: T) -> T {
    let (s, c) = sinh_cosh(x);
    s.div(c)
}

/// Hyperbolic cotangent. Domain: `x ≠ 0`.
///
/// # Errors
/// Returns `DomainError` if `x = 0`.
#[must_use = "returns the hyperbolic cotangent result which should be handled"]
pub fn coth<T: CordicNumber>(x: T) -> Result<T> {
    if x == T::zero() {
        return Err(Error::domain("coth", "non-zero value"));
    }
    let (s, c) = sinh_cosh(x);
    Ok(c.div(s))
}

/// Inverse hyperbolic sine. Accepts any value.
#[must_use]
pub fn asinh<T: CordicNumber>(x: T) -> T {
    if x == T::zero() {
        return T::zero();
    }

    // asinh(x) = atanh(x / sqrt(1 + x²))
    // NonNegative::one_plus_square(x) returns 1 + x², which is always ≥ 1
    let sqrt_term = sqrt_nonneg(NonNegative::one_plus_square(x));

    // x / sqrt(1 + x²) is always in (-1, 1) since sqrt(1 + x²) > |x|
    let arg = OpenUnitInterval::from_div_by_sqrt_one_plus_square(x, sqrt_term);

    atanh_open(arg)
}

/// Inverse hyperbolic cosine. Domain: `x ≥ 1`.
///
/// # Errors
/// Returns `DomainError` if `x < 1`.
#[must_use = "returns the inverse hyperbolic cosine result which should be handled"]
pub fn acosh<T: CordicNumber>(x: T) -> Result<T> {
    let at_least_one = AtLeastOne::new(x).ok_or_else(|| Error::domain("acosh", "value >= 1"))?;

    if x == T::one() {
        return Ok(T::zero());
    }

    // acosh(x) = atanh(sqrt(x² - 1) / x) for x > 1
    // NonNegative::square_minus_one gives x² - 1, which is ≥ 0 since x ≥ 1
    let sqrt_term = sqrt_nonneg(NonNegative::square_minus_one(at_least_one));

    // sqrt(x² - 1) / x is in (-1, 1) for x > 1 since sqrt(x² - 1) < x
    let arg = OpenUnitInterval::from_sqrt_square_minus_one_div(sqrt_term, at_least_one);

    Ok(atanh_open(arg))
}

/// Inverse hyperbolic tangent. Domain: `(-1, 1)`.
///
/// # Errors
/// Returns `DomainError` if `|x| ≥ 1`.
#[must_use = "returns the inverse hyperbolic tangent result which should be handled"]
pub fn atanh<T: CordicNumber>(x: T) -> Result<T> {
    OpenUnitInterval::new(x)
        .map(atanh_open)
        .ok_or_else(|| Error::domain("atanh", "value in range (-1, 1)"))
}

/// Infallible inverse hyperbolic tangent for values in (-1, 1).
///
/// This function takes an [`OpenUnitInterval<T>`] wrapper, guaranteeing at the
/// type level that the input is valid. No domain check is performed at runtime.
///
/// Use this when the input is known to be in (-1, 1) through mathematical
/// invariants (e.g., `x / sqrt(1 + x²)`).
#[must_use]
pub fn atanh_open<T: CordicNumber>(x: OpenUnitInterval<T>) -> T {
    atanh_core(x.get())
}

/// Core atanh implementation. Caller must ensure |x| < 1.
fn atanh_core<T: CordicNumber>(x: T) -> T {
    let zero = T::zero();
    let one = T::one();

    if x == zero {
        return zero;
    }

    let threshold = T::from_i1f63(ATANH_REDUCTION_THRESHOLD_I1F63);

    if x.abs() <= threshold {
        // Direct CORDIC computation
        let (_, _, z) = hyperbolic_vectoring(one, x, zero);
        return z;
    }

    // Argument reduction using the identity:
    // atanh(x) = atanh(a) + atanh((x - a) / (1 - a*x))
    // We use a = 0.5, for which atanh(0.5) is a precomputed constant.
    let half = T::half();
    let atanh_half = T::from_i1f63(crate::tables::hyperbolic::ATANH_HALF);

    let sign = if x.is_negative() { -one } else { one };
    let abs_x = x.abs();

    // Compute reduced argument: (|x| - 0.5) / (1 - 0.5*|x|)
    let numerator = abs_x - half;
    let denominator = one - half.saturating_mul(abs_x);
    let reduced = numerator.div(denominator);

    // Recursively compute atanh of reduced argument
    let atanh_reduced = atanh_core(reduced);

    // atanh(x) = sign * (atanh(0.5) + atanh(reduced))
    sign.saturating_mul(atanh_half.saturating_add(atanh_reduced))
}

/// Inverse hyperbolic cotangent. Domain: `|x| > 1`.
///
/// # Errors
/// Returns `DomainError` if `|x| ≤ 1`.
#[must_use = "returns the inverse hyperbolic cotangent result which should be handled"]
pub fn acoth<T: CordicNumber>(x: T) -> Result<T> {
    let one = T::one();

    if x.abs() <= one {
        return Err(Error::domain("acoth", "|value| > 1"));
    }

    // acoth(x) = atanh(1/x)
    let recip = one.div(x);
    Ok(atanh_core(recip))
}
