//! Exponential and logarithmic functions.

use crate::bounded::{NormalizedLnArg, OpenUnitInterval};
use crate::error::{Error, Result};
use crate::ops::hyperbolic::atanh_open;
use crate::traits::CordicNumber;

/// `1/ln 2 − 1` as I1F63 (`1/ln 2` itself exceeds 1).
const FRAC_1_LN2_MINUS_1_I1F63: i64 = 0x38AA_3B29_5C17_F0BC;

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

    // Argument reduction: exp(x) = 2^k * exp(r), where r ∈ [0, ln2), so
    // exp(r) ∈ [1, 2) carries a full significand into the scaling shift.
    // k = ⌊x / ln2⌋ is estimated to within one (`to_i32` floors; the
    // reciprocal multiply is within |x|·1.5·2^-frac + 2^-frac < 1 when the
    // integer bits do not outnumber the fractional bits, and other types
    // divide) and then corrected once in either direction.
    #[allow(clippy::cast_possible_wrap, reason = "bit counts bounded by type size")]
    let (int_bits, frac_bits) = (
        (T::total_bits() - T::frac_bits()) as i32,
        T::frac_bits() as i32,
    );
    let quotient = if T::total_bits() <= 2 * T::frac_bits() {
        x.saturating_add(x.saturating_mul(T::from_i1f63(FRAC_1_LN2_MINUS_1_I1F63)))
    } else {
        x.div(ln2)
    };
    let mut scale = quotient.to_i32();

    // 2^k·exp(r) exceeds MAX once k ≥ int_bits − 1 and rounds to zero once
    // k ≤ −(frac_bits + 2); allow one of slack until k is corrected.
    if scale >= int_bits {
        return T::max_value();
    }
    if scale < -(frac_bits + 2) {
        return zero;
    }

    let mut r = x.saturating_sub(ln2.mul_int(scale));
    if r < zero {
        scale -= 1;
        r = r.saturating_add(ln2);
    } else if r >= ln2 {
        scale += 1;
        r = r.saturating_sub(ln2);
    }
    if scale >= int_bits - 1 {
        return T::max_value();
    }
    if scale <= -(frac_bits + 2) {
        return zero;
    }

    // Factored Taylor: exp(r) = 1 + r*(1 + r/2*(1 + r/3*(1 + ... r/n))).
    // The omitted term (ln2)^(n+1)/(n+1)! is 1.3e-6 at degree 7 and 1.4e-12
    // at degree 12: below half an ulp up to 18 and 38 fractional bits.
    let mut p = one;
    if T::frac_bits() >= 24 {
        p = one.saturating_add(r.div_int(12).saturating_mul(p));
        p = one.saturating_add(r.div_int(11).saturating_mul(p));
        p = one.saturating_add(r.div_int(10).saturating_mul(p));
        p = one.saturating_add(r.div_int(9).saturating_mul(p));
        p = one.saturating_add(r.div_int(8).saturating_mul(p));
    }
    p = one.saturating_add(r.div_int(7).saturating_mul(p));
    p = one.saturating_add(r.div_int(6).saturating_mul(p));
    p = one.saturating_add(r.div_int(5).saturating_mul(p));
    p = one.saturating_add(r.div_int(4).saturating_mul(p));
    p = one.saturating_add(r.div_int(3).saturating_mul(p));
    p = one.saturating_add(r.div_int(2).saturating_mul(p));
    let exp_r = one.saturating_add(r.saturating_mul(p));

    // Scale by 2^scale using bit shifts.
    // scale is already bounded to (-(frac_bits + 2), int_bits - 1) by the exits above.
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

    // Reduce to [0.5, 2] from the leading bit e = ⌊log₂ x⌋: x > 2 shifts
    // right until ≤ 2 (stopping early when a shift lands exactly on 2),
    // x < 0.5 doubles until ≥ 0.5.
    let e = x.checked_int_log2().unwrap_or(0);
    #[allow(
        clippy::cast_sign_loss,
        reason = "shift counts are non-negative by construction"
    )]
    let (normalized, k) = if e >= 1 {
        let candidate = x >> (e - 1) as u32;
        if candidate == two {
            (candidate, e - 1)
        } else {
            (x >> e as u32, e)
        }
    } else if e <= -2 {
        (x << (-1 - e) as u32, e + 1)
    } else {
        (x, 0)
    };
    let k_ln2 = ln2.mul_int(k);

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
