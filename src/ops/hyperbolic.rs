#![doc = include_str!("../../docs/hyperbolic.md")]

use crate::bounded::{AtLeastOne, NonNegative, OpenUnitInterval, UnitInterval};
use crate::error::{Error, Result};
use crate::kernel::hyperbolic_vectoring;
use crate::ops::algebraic::sqrt_nonneg;
use crate::ops::exponential::ln_positive;
use crate::traits::CordicNumber;

/// The fractional part of 1.1182, the sum of the atanh table, which bounds the
/// polynomial's argument; I1F63 cannot hold 1.1182 itself.
const HYPERBOLIC_CONVERGENCE_LIMIT_FRAC_I1F63: i64 = 0x0F22_3D70_A3D7_0A3D;

/// 0.75, above which atanh moves its argument into the CORDIC's convergent
/// range with `atanh(x) = atanh(0.5) + atanh((x − 0.5) / (1 − 0.5x))`.
const ATANH_REDUCTION_THRESHOLD_I1F63: i64 = 0x6000_0000_0000_0000;

/// Hyperbolic sine and cosine, cheaper than separate calls. Both saturate
/// together for large `|x|`, as the [module docs](self) describe.
#[must_use]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn sinh_cosh<T: CordicNumber>(x: T) -> (T, T) {
    let one = T::one();
    let limit = one.saturating_add(T::from_i1f63(HYPERBOLIC_CONVERGENCE_LIMIT_FRAC_I1F63));

    // Halve until within the limit, then double back; total_bits halvings
    // reach zero, so the bound never cuts the loop short.
    let mut reduced = x;
    let mut depth: u32 = 0;
    while reduced.abs() > limit && depth < T::total_bits() {
        reduced = reduced >> 1;
        depth += 1;
    }

    // Factored Taylor series in u = x², with divisors (2k)(2k+1) for sinh and
    // (2k−1)(2k) for cosh; docs/sinh_cosh.md gives the degrees and errors.
    let u = reduced.saturating_mul(reduced);

    let mut sp = one;
    let mut cp = one;
    if T::frac_bits() >= 40 {
        sp = one.saturating_add(u.div_int(420).saturating_mul(sp));
        sp = one.saturating_add(u.div_int(342).saturating_mul(sp));
        sp = one.saturating_add(u.div_int(272).saturating_mul(sp));
        sp = one.saturating_add(u.div_int(210).saturating_mul(sp));
        cp = one.saturating_add(u.div_int(462).saturating_mul(cp));
        cp = one.saturating_add(u.div_int(380).saturating_mul(cp));
        cp = one.saturating_add(u.div_int(306).saturating_mul(cp));
        cp = one.saturating_add(u.div_int(240).saturating_mul(cp));
    }
    if T::frac_bits() >= 24 {
        sp = one.saturating_add(u.div_int(156).saturating_mul(sp));
        sp = one.saturating_add(u.div_int(110).saturating_mul(sp));
        cp = one.saturating_add(u.div_int(182).saturating_mul(cp));
        cp = one.saturating_add(u.div_int(132).saturating_mul(cp));
    }
    sp = one.saturating_add(u.div_int(72).saturating_mul(sp));
    sp = one.saturating_add(u.div_int(42).saturating_mul(sp));
    sp = one.saturating_add(u.div_int(20).saturating_mul(sp));
    sp = one.saturating_add(u.div_int(6).saturating_mul(sp));
    cp = one.saturating_add(u.div_int(90).saturating_mul(cp));
    cp = one.saturating_add(u.div_int(56).saturating_mul(cp));
    cp = one.saturating_add(u.div_int(30).saturating_mul(cp));
    cp = one.saturating_add(u.div_int(12).saturating_mul(cp));
    cp = one.saturating_add(u.div_int(2).saturating_mul(cp));
    let (mut sh, mut ch) = (reduced.saturating_mul(sp), cp);

    // sinh(2x) = 2·sinh(x)·cosh(x) and cosh(2x) = cosh²(x) + sinh²(x).
    for _ in 0..depth {
        let new_sh = sh.saturating_mul(ch).saturating_mul(T::two());
        let new_ch = ch.saturating_mul(ch).saturating_add(sh.saturating_mul(sh));
        sh = new_sh;
        ch = new_ch;
    }

    (sh, ch)
}

/// Hyperbolic sine. Saturates to `T::MAX` or `T::MIN` past `±asinh(T::MAX)`;
/// the [module docs](self) list the thresholds.
#[inline]
#[must_use]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn sinh<T: CordicNumber>(x: T) -> T {
    sinh_cosh(x).0
}

/// Hyperbolic cosine, always at least 1. Saturates to `T::MAX` past
/// `acosh(T::MAX)`; the [module docs](self) list the thresholds.
#[inline]
#[must_use]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn cosh<T: CordicNumber>(x: T) -> T {
    sinh_cosh(x).1
}

/// Hyperbolic tangent. Result in `(-1, 1)`.
#[must_use]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn tanh<T: CordicNumber>(x: T) -> T {
    let (s, c) = sinh_cosh(x);
    s.div(c)
}

/// Hyperbolic cotangent. Domain: `x ≠ 0`.
/// # Errors
/// Returns `DomainError` if `x = 0`.
#[must_use = "returns the hyperbolic cotangent result which should be handled"]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn coth<T: CordicNumber>(x: T) -> Result<T> {
    if x == T::zero() {
        return Err(Error::domain("coth", "non-zero value"));
    }
    let (s, c) = sinh_cosh(x);
    Ok(c.div(s))
}

/// Inverse hyperbolic sine. Accepts any value.
#[must_use]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn asinh<T: CordicNumber>(x: T) -> T {
    let zero = T::zero();
    let one = T::one();

    if x == zero {
        return zero;
    }

    let abs_x = x.abs();
    if abs_x <= one {
        // asinh(x) = atanh(x / √(1 + x²)), whose argument is at most 1/√2.
        let sqrt_term = sqrt_nonneg(NonNegative::one_plus_square(x));
        let arg = OpenUnitInterval::from_div_by_sqrt_one_plus_square(x, sqrt_term);
        return atanh_open(arg);
    }

    // Past 1 that argument nears 1, where each atanh reduction step triples
    // the error, so use sign(x)·ln(|x|·(1 + √(1 + (1/|x|)²))).
    let recip = one.div(abs_x);
    let factor = one.saturating_add(sqrt_nonneg(NonNegative::one_plus_square(recip)));
    let magnitude = if abs_x > T::max_value() >> 2 {
        // |x|·factor would overflow, as factor < 2.42.
        ln_positive(abs_x).saturating_add(ln_positive(factor))
    } else {
        ln_positive(abs_x.saturating_mul(factor))
    };
    if x.is_negative() {
        -magnitude
    } else {
        magnitude
    }
}

/// Inverse hyperbolic cosine. Domain: `x ≥ 1`.
/// # Errors
/// Returns `DomainError` if `x < 1`.
#[must_use = "returns the inverse hyperbolic cosine result which should be handled"]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn acosh<T: CordicNumber>(x: T) -> Result<T> {
    let at_least_one = AtLeastOne::new(x).ok_or_else(|| Error::domain("acosh", "value >= 1"))?;
    let zero = T::zero();
    let one = T::one();

    if x == one {
        return Ok(zero);
    }

    if x < one.saturating_add(T::half()) {
        // acosh(x) = atanh(√(x² − 1) / x), whose argument is under 0.75 here.
        let sqrt_term = sqrt_nonneg(NonNegative::square_minus_one(at_least_one));
        let arg = OpenUnitInterval::from_sqrt_square_minus_one_div(sqrt_term, at_least_one);
        return Ok(atanh_open(arg));
    }

    // Past 1.5 that argument nears 1, where each atanh reduction step triples
    // the error, so use ln(x·(1 + √(1 − (1/x)²))).
    let recip = UnitInterval::from_reciprocal(at_least_one);
    let factor = one.saturating_add(sqrt_nonneg(NonNegative::one_minus_square(recip)));
    let result = if x > T::max_value() >> 1 {
        // x·factor would overflow, as factor < 2.
        ln_positive(x).saturating_add(ln_positive(factor))
    } else {
        ln_positive(x.saturating_mul(factor))
    };
    Ok(result)
}

/// Inverse hyperbolic tangent. Domain: `(-1, 1)`.
/// # Errors
/// Returns `DomainError` if `|x| ≥ 1`.
#[must_use = "returns the inverse hyperbolic tangent result which should be handled"]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn atanh<T: CordicNumber>(x: T) -> Result<T> {
    OpenUnitInterval::new(x)
        .map(atanh_open)
        .ok_or_else(|| Error::domain("atanh", "value in range (-1, 1)"))
}

/// Inverse hyperbolic tangent of an [`OpenUnitInterval`] value, with no
/// runtime domain check.
#[must_use]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn atanh_open<T: CordicNumber>(x: OpenUnitInterval<T>) -> T {
    atanh_core(x.get())
}

/// atanh for `|x| < 1`, which the caller must ensure.
fn atanh_core<T: CordicNumber>(x: T) -> T {
    let zero = T::zero();
    let one = T::one();

    if x == zero {
        return zero;
    }

    let threshold = T::from_i1f63(ATANH_REDUCTION_THRESHOLD_I1F63);

    if x.abs() <= threshold {
        let (_, _, z) = hyperbolic_vectoring(one, x, zero);
        return z;
    }

    let half = T::half();
    let atanh_half = T::from_i1f63(crate::tables::hyperbolic::ATANH_HALF);
    let sign = if x.is_negative() { -one } else { one };
    let mut abs_x = x.abs();

    // Each step takes atanh(0.5) off and roughly halves the distance to the
    // threshold, so frac_bits steps suffice.
    let mut accumulated = zero;
    let mut i: u32 = 0;
    while abs_x > threshold && i < T::frac_bits() {
        let numerator = abs_x.saturating_sub(half);
        let denominator = one.saturating_sub(half.saturating_mul(abs_x));
        abs_x = numerator.div(denominator);
        accumulated = accumulated.saturating_add(atanh_half);
        i += 1;
    }

    let (_, _, z) = hyperbolic_vectoring(one, abs_x, zero);

    sign.saturating_mul(accumulated.saturating_add(z))
}

/// Inverse hyperbolic cotangent. Domain: `|x| > 1`.
/// # Errors
/// Returns `DomainError` if `|x| ≤ 1`.
#[must_use = "returns the inverse hyperbolic cotangent result which should be handled"]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn acoth<T: CordicNumber>(x: T) -> Result<T> {
    let one = T::one();

    if x.abs() <= one {
        return Err(Error::domain("acoth", "|value| > 1"));
    }

    let recip = one.div(x);
    Ok(atanh_core(recip))
}
