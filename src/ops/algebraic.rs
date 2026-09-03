//! Algebraic functions (sqrt).

use crate::bounded::NonNegative;
use crate::error::{Error, Result};
use crate::traits::CordicNumber;

/// Square root, rounded to the nearest representable value. Domain: `x ≥ 0`.
///
/// # Errors
/// Returns `DomainError` if `x < 0`.
#[must_use = "returns the square root result which should be handled"]
#[cfg_attr(feature = "verify-no-panic", no_panic::no_panic)]
pub fn sqrt<T: CordicNumber>(x: T) -> Result<T> {
    NonNegative::new(x)
        .map(sqrt_nonneg)
        .ok_or_else(|| Error::domain("sqrt", "non-negative value"))
}

/// Infallible square root for non-negative values, rounded to the nearest
/// representable value.
///
/// This function takes a [`NonNegative<T>`] wrapper, guaranteeing at the type
/// level that the input is valid. No domain check is performed at runtime.
///
/// Use this when the non-negativity of the input is already established
/// through mathematical invariants (e.g., `1 + x²`, `1 - x²` for `|x| ≤ 1`).
///
/// The result is correctly rounded at every magnitude (raw bits
/// `⌊√(X·2^f) + ½⌋`), computed from an integer square root of the raw
/// representation; see [`CordicNumber::sqrt_round`].
#[must_use]
#[cfg_attr(feature = "verify-no-panic", no_panic::no_panic)]
pub fn sqrt_nonneg<T: CordicNumber>(x: NonNegative<T>) -> T {
    x.get().sqrt_round()
}
