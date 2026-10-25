#![doc = include_str!("../../docs/exponential.md")]

use crate::bounded::{NormalizedLnArg, OpenUnitInterval};
use crate::error::{Error, Result};
use crate::ops::hyperbolic::atanh_open;
use crate::traits::CordicNumber;

/// `1/ln 2 − 1` as I1F63 (`1/ln 2` itself exceeds 1).
const FRAC_1_LN2_MINUS_1_I1F63: i64 = 0x38AA_3B29_5C17_F0BC;

/// Right-shifts a non-negative value by `n`, rounding to nearest. A plain `>>`
/// would bias small `exp` results, only a few ulps, low by up to an ulp.
fn shr_rounded<T: CordicNumber>(x: T, n: u32) -> T {
    if n == 0 {
        return x;
    }
    let ulp = T::one() >> T::frac_bits();
    let shifted = x >> n;
    let rem = x.saturating_sub(shifted << n);
    let half = ulp << (n - 1);
    if rem >= half {
        shifted.saturating_add(ulp)
    } else {
        shifted
    }
}

/// Exponential function e^x. Saturates to `T::MAX` above `ln(T::MAX)` and to
/// zero far below 0; the [module docs](self) list the thresholds.
#[must_use]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn exp<T: CordicNumber>(x: T) -> T {
    let zero = T::zero();
    let one = T::one();
    let ln2 = T::ln_2();

    if x == zero {
        return one;
    }

    // exp(x) = 2^k·exp(r) with r in [0, ln 2). The estimate of k = ⌊x / ln 2⌋
    // is within one, so one correction below fixes it (docs/exp.md).
    let (int_bits, frac_bits) = (
        (T::total_bits() - T::frac_bits()).cast_signed(),
        T::frac_bits().cast_signed(),
    );
    let quotient = if T::total_bits() <= 2 * T::frac_bits() {
        x.saturating_add(x.saturating_mul(T::from_i1f63(FRAC_1_LN2_MINUS_1_I1F63)))
    } else {
        x.div(ln2)
    };
    let mut scale = quotient.to_i32();

    // The result exceeds MAX once k ≥ int_bits − 1 and rounds to zero once
    // k ≤ −(frac_bits + 2); these exits allow one of slack before correction.
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

    // Factored Taylor series, of degree 7, 12 or 19 by precision (docs/exp.md).
    let mut p = one;
    if T::frac_bits() >= 40 {
        p = one.saturating_add(r.div_int(19).saturating_mul(p));
        p = one.saturating_add(r.div_int(18).saturating_mul(p));
        p = one.saturating_add(r.div_int(17).saturating_mul(p));
        p = one.saturating_add(r.div_int(16).saturating_mul(p));
        p = one.saturating_add(r.div_int(15).saturating_mul(p));
        p = one.saturating_add(r.div_int(14).saturating_mul(p));
        p = one.saturating_add(r.div_int(13).saturating_mul(p));
    }
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

    if scale > 0 {
        let shift = scale.unsigned_abs();
        // A left shift past the headroom would wrap.
        let headroom = T::max_value() >> shift;
        if exp_r > headroom {
            T::max_value()
        } else {
            exp_r << shift
        }
    } else {
        shr_rounded(exp_r, scale.unsigned_abs())
    }
}

/// Power `base^exponent`, as `exp(exponent · ln|base|)`, saturating like [`exp`].
///
/// A negative base takes its sign from the exponent's parity. The relative
/// error grows with `|exponent|`, which scales the logarithm's rounding error.
/// # Errors
/// Returns `DomainError` for a zero base with a negative exponent, or a
/// negative base with a non-integer exponent, which has no real result.
#[must_use = "returns the power result which should be handled"]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn pow<T: CordicNumber>(base: T, exponent: T) -> Result<T> {
    let zero = T::zero();
    if exponent == zero {
        return Ok(T::one());
    }
    if base == zero {
        return if exponent > zero {
            Ok(zero)
        } else {
            Err(Error::domain(
                "pow",
                "non-negative exponent for a zero base",
            ))
        };
    }
    let magnitude = exp(exponent.saturating_mul(ln_positive(base.abs())));
    if !base.is_negative() {
        return Ok(magnitude);
    }
    // The shift isolates the integer part exactly at any magnitude.
    let integer_part = exponent >> T::frac_bits();
    if integer_part << T::frac_bits() != exponent {
        return Err(Error::domain("pow", "integer exponent for a negative base"));
    }
    let odd = (integer_part >> 1) << 1 != integer_part;
    Ok(if odd { -magnitude } else { magnitude })
}

/// Natural logarithm. Domain: `x > 0`.
/// # Errors
/// Returns `DomainError` if `x ≤ 0`.
#[must_use = "returns the natural logarithm result which should be handled"]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn ln<T: CordicNumber>(x: T) -> Result<T> {
    if x <= T::zero() {
        return Err(Error::domain("ln", "positive value"));
    }
    Ok(ln_positive(x))
}

/// Natural logarithm of a positive value. The caller must ensure `x > 0`.
pub(crate) fn ln_positive<T: CordicNumber>(x: T) -> T {
    let zero = T::zero();
    let one = T::one();
    let two = T::two();

    if x == one {
        return zero;
    }

    let ln2 = T::ln_2();

    // ln(x) = ln(x·2^-k) + k·ln 2, with x·2^-k in [0.5, 2] chosen from the
    // leading bit e = ⌊log₂ x⌋; a value landing exactly on 2 stops there.
    let e = x.checked_int_log2().unwrap_or(0);
    let (normalized, k) = if e >= 1 {
        let candidate = x >> (e - 1).unsigned_abs();
        if candidate == two {
            (candidate, e - 1)
        } else {
            (x >> e.unsigned_abs(), e)
        }
    } else if e <= -2 {
        (x << (-1 - e).unsigned_abs(), e + 1)
    } else {
        (x, 0)
    };
    let k_ln2 = ln2.mul_int(k);

    // ln(x) = 2·atanh((x − 1)/(x + 1)), whose argument is within 1/3 of 0.
    let norm = NormalizedLnArg::from_normalized(normalized);
    let arg = OpenUnitInterval::from_normalized_ln_arg(norm);

    let atanh_val = atanh_open(arg);
    let ln_normalized = atanh_val.saturating_add(atanh_val);

    ln_normalized.saturating_add(k_ln2)
}

/// Base-2 logarithm. Domain: `x > 0`.
/// # Errors
/// Returns `DomainError` if `x ≤ 0`.
#[must_use = "returns the base-2 logarithm result which should be handled"]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn log2<T: CordicNumber>(x: T) -> Result<T> {
    let ln_x = ln(x)?;
    let ln_2 = T::ln_2();
    Ok(ln_x.div(ln_2))
}

/// Base-10 logarithm. Domain: `x > 0`.
/// # Errors
/// Returns `DomainError` if `x ≤ 0`.
#[must_use = "returns the base-10 logarithm result which should be handled"]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn log10<T: CordicNumber>(x: T) -> Result<T> {
    let ln_x = ln(x)?;
    let ln_10 = T::ln_10();
    Ok(ln_x.div(ln_10))
}

/// Power of two 2^x, as `exp(x · ln 2)`. Saturates to `T::MAX` above
/// `log2(T::MAX)` and to zero far below 0; the [module docs](self) list the
/// thresholds.
#[must_use]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn pow2<T: CordicNumber>(x: T) -> T {
    let ln_2 = T::ln_2();
    exp(x.saturating_mul(ln_2))
}
