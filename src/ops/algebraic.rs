//! Algebraic functions (sqrt).

use crate::bounded::NonNegative;
use crate::error::{Error, Result};
use crate::traits::CordicNumber;

/// Square root. Domain: `x ≥ 0`. Uses Newton-Raphson iteration.
///
/// # Errors
/// Returns `DomainError` if `x < 0`.
#[must_use = "returns the square root result which should be handled"]
pub fn sqrt<T: CordicNumber>(x: T) -> Result<T> {
    NonNegative::new(x)
        .map(sqrt_nonneg)
        .ok_or_else(|| Error::domain("sqrt", "non-negative value"))
}

/// Infallible square root for non-negative values.
///
/// This function takes a [`NonNegative<T>`] wrapper, guaranteeing at the type
/// level that the input is valid. No domain check is performed at runtime.
///
/// Use this when the non-negativity of the input is already established
/// through mathematical invariants (e.g., `1 + x²`, `1 - x²` for `|x| ≤ 1`).
#[must_use]
pub fn sqrt_nonneg<T: CordicNumber>(x: NonNegative<T>) -> T {
    let x = x.get();
    let zero = T::zero();
    let one = T::one();
    let half = T::half();

    if x == zero {
        return zero;
    }

    if x == one {
        return one;
    }

    // Initial guess: use bit-level estimation for faster convergence
    // For sqrt(x), a good initial estimate is 2^(floor(log2(x))/2)
    let mut guess = if x > one {
        // For large numbers, use bit estimation: sqrt(x) ≈ 2^(log2(x)/2)
        // We find the approximate position by successive squaring comparison
        let mut g = one;
        let four = T::two().saturating_mul(T::two());
        let mut test = x;
        // Guard: max 64 iterations sufficient for any representable value
        let mut iter_guard = 0u32;
        while test >= four && iter_guard < 64 {
            test = test >> 2;
            g = g << 1;
            iter_guard += 1;
        }
        // g is now approximately sqrt(x) rounded down to a power of 2
        // Refine: average with x/g for a better starting point
        let quotient = x.div(g);
        g.saturating_add(quotient) >> 1
    } else {
        // For numbers < 1, x is a reasonable starting point
        // since sqrt(x) > x when 0 < x < 1
        x
    };

    // Newton-Raphson iteration: x_new = (x_old + n/x_old) / 2
    // Number of iterations depends on precision needed.
    // Newton-Raphson for sqrt converges quadratically, so 8-20 iterations
    // is sufficient for any fixed-point precision up to 128 bits.
    let iterations = (T::frac_bits() / 2).clamp(8, 20);

    // Pre-compute epsilon: approximately 2^(-frac_bits/2) for convergence check.
    // frac_bits ≤ 128 for all supported types, so shift is in range [0, 63].
    #[allow(
        clippy::cast_possible_wrap,
        clippy::cast_sign_loss,
        reason = "frac_bits bounded by type size"
    )]
    let epsilon_shift = (63i32 - (T::frac_bits() / 2) as i32).max(0) as u32;
    let epsilon = T::from_i1f63(1i64 << epsilon_shift);

    // Run iterations - 1 times with early exit on convergence
    for _ in 0..iterations.saturating_sub(1) {
        let quotient = x.div(guess);
        let sum = guess.saturating_add(quotient);
        let new_guess = sum.saturating_mul(half);

        let diff = if new_guess > guess {
            new_guess - guess
        } else {
            guess - new_guess
        };

        if diff <= epsilon {
            return new_guess;
        }

        guess = new_guess;
    }

    // Final iteration - always performed, result always returned
    let quotient = x.div(guess);
    let sum = guess.saturating_add(quotient);
    sum.saturating_mul(half)
}
