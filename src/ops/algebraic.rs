//! Algebraic functions (sqrt).

use crate::bounded::NonNegative;
use crate::error::{Error, Result};
use crate::traits::CordicNumber;

/// Square root, rounded to the nearest representable value. Domain: `x ≥ 0`.
/// # Errors
/// Returns `DomainError` if `x < 0`.
#[must_use = "returns the square root result which should be handled"]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn sqrt<T: CordicNumber>(x: T) -> Result<T> {
    NonNegative::new(x)
        .map(sqrt_nonneg)
        .ok_or_else(|| Error::domain("sqrt", "non-negative value"))
}

/// Square root of a [`NonNegative`] value, with no runtime domain check.
///
/// The result is correctly rounded at every magnitude: raw bits
/// `⌊√(X·2^f) + ½⌋`, as [`CordicNumber::sqrt_round`] computes them.
#[must_use]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn sqrt_nonneg<T: CordicNumber>(x: NonNegative<T>) -> T {
    x.get().sqrt_round()
}
