#![doc = include_str!("../../docs/exponential.md")]

use crate::consts::{FRAC_1_LN_10, LN_2};
use crate::error::{Error, Result};
use crate::eval::{HUGE, W, exp as exp_work, fixed, ln as ln_work, log2 as log2_work, work};
use crate::real::Work;
use crate::traits::CordicNumber;

/// `e^v`, or for `|v|` past [`HUGE`] a value every type saturates or rounds
/// to zero from.
fn exp_wide<T: CordicNumber>(v: W<T>) -> W<T> {
    if v.abs() > W::<T>::from_int(HUGE) {
        let edge = W::<T>::one().ldexp(256);
        return if v.is_neg() { W::<T>::zero() } else { edge };
    }
    exp_work(v)
}

/// Exponential function e^x. Saturates to `T::MAX` above `ln(T::MAX)` and to
/// zero where e^x rounds below half an ulp; the [module docs](self) list the
/// thresholds.
#[must_use]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn exp<T: CordicNumber>(x: T) -> T {
    fixed(exp_wide::<T>(work(x)))
}

/// Power `base^exponent`, as `exp(exponent · ln|base|)`, saturating like [`exp`].
///
/// A negative base takes its sign from the exponent's parity.
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
    let magnitude = exp_wide::<T>(work(exponent) * ln_work(work(base).abs()));
    if !base.is_negative() {
        return Ok(fixed(magnitude));
    }
    // The shift isolates the integer part exactly at any magnitude.
    let integer_part = exponent >> T::frac_bits();
    if integer_part << T::frac_bits() != exponent {
        return Err(Error::domain("pow", "integer exponent for a negative base"));
    }
    let odd = (integer_part >> 1) << 1 != integer_part;
    Ok(fixed(if odd { -magnitude } else { magnitude }))
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
    Ok(fixed(ln_work(work(x))))
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
    if x <= T::zero() {
        return Err(Error::domain("log2", "positive value"));
    }
    Ok(fixed(log2_work(work(x))))
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
    if x <= T::zero() {
        return Err(Error::domain("log10", "positive value"));
    }
    Ok(fixed(ln_work(work(x)) * W::<T>::constant(&FRAC_1_LN_10)))
}

/// Power of two 2^x. Saturates to `T::MAX` above `log2(T::MAX)` and to zero
/// where 2^x rounds below half an ulp; the [module docs](self) list the
/// thresholds.
#[must_use]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn pow2<T: CordicNumber>(x: T) -> T {
    // 2^x = 2^k · e^(r·ln 2) with k = ⌊x⌋ and r its exact fractional part.
    let k = x.raw_i128() >> T::frac_bits();
    if k > 256 {
        return T::max_value();
    }
    if k < -256 {
        return T::zero();
    }
    let r = W::<T>::from_raw(x.raw_i128() - (k << T::frac_bits()), T::frac_bits());
    let k = i32::try_from(k).unwrap_or(0);
    fixed(exp_work(r * W::<T>::constant(&LN_2)).ldexp(k))
}
