//! Hyperbolic functions via hyperbolic CORDIC.

use crate::bounded::{AtLeastOne, NonNegative, OpenUnitInterval};
use crate::error::{Error, Result};
use crate::kernel::hyperbolic_vectoring;
use crate::ops::algebraic::sqrt_nonneg;
use crate::traits::CordicNumber;

/// Hyperbolic CORDIC converges for |x| < sum of atanh table ≈ 1.1182.
/// Stored as fractional part (0.1182) since I1F63 can't hold 1.x.
/// Full limit = 1 + this value.
const HYPERBOLIC_CONVERGENCE_LIMIT_FRAC_I1F63: i64 = 0x0F22_3D70_A3D7_0A3D;

/// Argument reduction threshold for atanh.
///
/// For |x| > this threshold, atanh uses the identity:
/// ```text
/// atanh(x) = atanh(0.5) + atanh((x - 0.5) / (1 - 0.5x))
/// ```
/// to reduce the argument into CORDIC's optimal convergence range.
///
/// Value 0.75 keeps reduced arguments within the convergent region.
const ATANH_REDUCTION_THRESHOLD_I1F63: i64 = 0x6000_0000_0000_0000;

/// Hyperbolic sine and cosine. More efficient than separate calls.
///
/// # Saturation Behavior
///
/// Both outputs can saturate for large |x|:
/// - sinh saturates to `T::MAX` or `T::MIN`
/// - cosh saturates to `T::MAX`
///
/// See [`sinh`] and [`cosh`] for threshold details.
///
/// When saturation occurs, both values saturate together (they grow
/// at the same rate), so the relationship cosh²(x) - sinh²(x) = 1
/// will not hold for saturated outputs.
#[must_use]
#[cfg_attr(feature = "verify-no-panic", no_panic::no_panic)]
pub fn sinh_cosh<T: CordicNumber>(x: T) -> (T, T) {
    let one = T::one();
    // Compute limit as 1 + fractional_part (~1.1182)
    let limit = one.saturating_add(T::from_i1f63(HYPERBOLIC_CONVERGENCE_LIMIT_FRAC_I1F63));

    // Iterative argument reduction for large values.
    // Count how many halvings are needed, then rebuild on the way back.
    let mut reduced = x;
    let mut depth: u32 = 0;
    // Max depth bounded by bit width: each halving shifts right by 1,
    // so after total_bits iterations the value is zero and below the limit.
    while reduced.abs() > limit && depth < T::total_bits() {
        reduced = reduced >> 1;
        depth += 1;
    }

    // Factored Taylor evaluation. After argument reduction, |reduced| ≤ 1.1182.
    //
    // sinh(x) = x * (1 + u/6 * (1 + u/20 * (1 + u/42 * ...)))
    // cosh(x) = 1 + u/2 * (1 + u/12 * (1 + u/30 * ...))
    // where u = x².
    //
    // Integer division (u/K) is used because it computes the quotient in one
    // step without pre-rounding a reciprocal, and the factored form keeps each
    // step's multiplicand u/K below 1 for K ≥ 2.
    let u = reduced.saturating_mul(reduced);

    // High-precision path requires enough integer bits for divisors up to 182.
    let mut sp = one;
    let (mut sh, mut ch) = if T::frac_bits() >= 24 && T::total_bits() >= T::frac_bits() + 9 {
        // High precision: degree 13 sinh, degree 14 cosh
        sp = one.saturating_add(u.div(T::from_num(156)).saturating_mul(sp));
        sp = one.saturating_add(u.div(T::from_num(110)).saturating_mul(sp));
        sp = one.saturating_add(u.div(T::from_num(72)).saturating_mul(sp));
        sp = one.saturating_add(u.div(T::from_num(42)).saturating_mul(sp));
        sp = one.saturating_add(u.div(T::from_num(20)).saturating_mul(sp));
        sp = one.saturating_add(u.div(T::from_num(6)).saturating_mul(sp));
        let sinh_approx = reduced.saturating_mul(sp);

        let mut cp = one;
        cp = one.saturating_add(u.div(T::from_num(182)).saturating_mul(cp));
        cp = one.saturating_add(u.div(T::from_num(132)).saturating_mul(cp));
        cp = one.saturating_add(u.div(T::from_num(90)).saturating_mul(cp));
        cp = one.saturating_add(u.div(T::from_num(56)).saturating_mul(cp));
        cp = one.saturating_add(u.div(T::from_num(30)).saturating_mul(cp));
        cp = one.saturating_add(u.div(T::from_num(12)).saturating_mul(cp));
        cp = one.saturating_add(u.div(T::from_num(2)).saturating_mul(cp));
        (sinh_approx, cp)
    } else {
        // Low precision: degree 9 sinh, degree 10 cosh
        sp = one;
        sp = one.saturating_add(u.div(T::from_num(72)).saturating_mul(sp));
        sp = one.saturating_add(u.div(T::from_num(42)).saturating_mul(sp));
        sp = one.saturating_add(u.div(T::from_num(20)).saturating_mul(sp));
        sp = one.saturating_add(u.div(T::from_num(6)).saturating_mul(sp));
        let sinh_approx = reduced.saturating_mul(sp);

        let mut cp = one;
        cp = one.saturating_add(u.div(T::from_num(90)).saturating_mul(cp));
        cp = one.saturating_add(u.div(T::from_num(56)).saturating_mul(cp));
        cp = one.saturating_add(u.div(T::from_num(30)).saturating_mul(cp));
        cp = one.saturating_add(u.div(T::from_num(12)).saturating_mul(cp));
        cp = one.saturating_add(u.div(T::from_num(2)).saturating_mul(cp));
        (sinh_approx, cp)
    };

    // Reconstruct via doubling: sinh(2x) = 2·sinh(x)·cosh(x),
    //                           cosh(2x) = cosh²(x) + sinh²(x)
    for _ in 0..depth {
        let new_sh = sh.saturating_mul(ch).saturating_mul(T::two());
        let new_ch = ch.saturating_mul(ch).saturating_add(sh.saturating_mul(sh));
        sh = new_sh;
        ch = new_ch;
    }

    (sh, ch)
}

/// Hyperbolic sine.
///
/// # Saturation Behavior
///
/// sinh(x) grows exponentially: sinh(x) ≈ e^x/2 for large |x|.
/// This function saturates for extreme inputs:
///
/// | Condition | Result | Example (I16F16) |
/// |-----------|--------|------------------|
/// | x > `asinh(T::MAX)` | `T::MAX` | x > ~11.1 → 32767.99 |
/// | x < `-asinh(T::MAX)` | `T::MIN` | x < ~-11.1 → -32768.0 |
///
/// The exact thresholds:
/// - **I16F16:** Saturates for |x| > ~11.1
/// - **I32F32:** Saturates for |x| > ~22.2
///
/// Within the non-saturating range, sinh is computed via polynomial
/// evaluation with argument reduction (halving/doubling) for |x| > 1.118.
#[inline]
#[must_use]
#[cfg_attr(feature = "verify-no-panic", no_panic::no_panic)]
pub fn sinh<T: CordicNumber>(x: T) -> T {
    sinh_cosh(x).0
}

/// Hyperbolic cosine. Always ≥ 1.
///
/// # Saturation Behavior
///
/// cosh(x) grows exponentially: cosh(x) ≈ e^|x|/2 for large |x|.
/// This function saturates for extreme inputs:
///
/// | Condition | Result | Example (I16F16) |
/// |-----------|--------|------------------|
/// | \|x\| > `acosh(T::MAX)` | `T::MAX` | \|x\| > ~11.1 → 32767.99 |
///
/// Unlike sinh, cosh is always positive, so it only saturates to `T::MAX`.
///
/// The exact thresholds:
/// - **I16F16:** Saturates for |x| > ~11.1
/// - **I32F32:** Saturates for |x| > ~22.2
#[inline]
#[must_use]
#[cfg_attr(feature = "verify-no-panic", no_panic::no_panic)]
pub fn cosh<T: CordicNumber>(x: T) -> T {
    sinh_cosh(x).1
}

/// Hyperbolic tangent. Result in `(-1, 1)`.
#[must_use]
#[cfg_attr(feature = "verify-no-panic", no_panic::no_panic)]
pub fn tanh<T: CordicNumber>(x: T) -> T {
    let (s, c) = sinh_cosh(x);
    s.div(c)
}

/// Hyperbolic cotangent. Domain: `x ≠ 0`.
///
/// # Errors
/// Returns `DomainError` if `x = 0`.
#[must_use = "returns the hyperbolic cotangent result which should be handled"]
#[cfg_attr(feature = "verify-no-panic", no_panic::no_panic)]
pub fn coth<T: CordicNumber>(x: T) -> Result<T> {
    if x == T::zero() {
        return Err(Error::domain("coth", "non-zero value"));
    }
    let (s, c) = sinh_cosh(x);
    Ok(c.div(s))
}

/// Inverse hyperbolic sine. Accepts any value.
#[must_use]
#[cfg_attr(feature = "verify-no-panic", no_panic::no_panic)]
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
#[cfg_attr(feature = "verify-no-panic", no_panic::no_panic)]
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
#[cfg_attr(feature = "verify-no-panic", no_panic::no_panic)]
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
#[cfg_attr(feature = "verify-no-panic", no_panic::no_panic)]
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

    // Fast path: no argument reduction needed, use CORDIC directly.
    if x.abs() <= threshold {
        let (_, _, z) = hyperbolic_vectoring(one, x, zero);
        return z;
    }

    // Argument reduction needed. Work with |x| and track sign.
    let half = T::half();
    let atanh_half = T::from_i1f63(crate::tables::hyperbolic::ATANH_HALF);
    let sign = if x.is_negative() { -one } else { one };
    let mut abs_x = x.abs();

    // Iterative argument reduction: each step reduces |x| and accumulates
    // atanh(0.5) into the result. Max iterations bounded by frac_bits since
    // each reduction roughly halves the distance from the threshold.
    let mut accumulated = zero;
    let mut i: u32 = 0;
    while abs_x > threshold && i < T::frac_bits() {
        // atanh(x) = atanh(0.5) + atanh((x - 0.5) / (1 - 0.5*x))
        let numerator = abs_x.saturating_sub(half);
        let denominator = one.saturating_sub(half.saturating_mul(abs_x));
        abs_x = numerator.div(denominator);
        accumulated = accumulated.saturating_add(atanh_half);
        i += 1;
    }

    // Direct CORDIC computation on the reduced argument
    let (_, _, z) = hyperbolic_vectoring(one, abs_x, zero);

    sign.saturating_mul(accumulated.saturating_add(z))
}

/// Inverse hyperbolic cotangent. Domain: `|x| > 1`.
///
/// # Errors
/// Returns `DomainError` if `|x| ≤ 1`.
#[must_use = "returns the inverse hyperbolic cotangent result which should be handled"]
#[cfg_attr(feature = "verify-no-panic", no_panic::no_panic)]
pub fn acoth<T: CordicNumber>(x: T) -> Result<T> {
    let one = T::one();

    if x.abs() <= one {
        return Err(Error::domain("acoth", "|value| > 1"));
    }

    // acoth(x) = atanh(1/x)
    let recip = one.div(x);
    Ok(atanh_core(recip))
}
