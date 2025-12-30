//! Algebraic functions.
//!
//! Provides sqrt using an iterative algorithm optimized for fixed-point.
//!
//! # Algorithm
//!
//! Square root is computed using a digit-by-digit method similar to
//! long division, which is well-suited for fixed-point arithmetic.

use crate::traits::CordicNumber;

/// Computes the square root of a value.
///
/// Uses a Newton-Raphson iteration method that is efficient
/// for fixed-point numbers.
///
/// # Arguments
///
/// * `x` - A non-negative value
///
/// # Returns
///
/// The square root of `x`. Returns 0 for negative inputs.
///
/// # Examples
///
/// ```
/// use fixed::types::I16F16;
/// use fixed_analytics::sqrt;
///
/// let x = I16F16::from_num(4.0);
/// let result = sqrt(x);
/// // result ≈ 2.0
///
/// let x = I16F16::from_num(2.0);
/// let result = sqrt(x);
/// // result ≈ 1.414
/// ```
///
/// # Algorithm Details
///
/// This implementation uses Newton-Raphson iteration:
/// ```text
/// x_{n+1} = (x_n + S/x_n) / 2
/// ```
/// Starting from an initial guess, this converges quadratically to sqrt(S).
#[must_use]
pub fn sqrt<T: CordicNumber>(x: T) -> T {
    let zero = T::zero();
    let one = T::one();
    let half = T::half();

    // Handle edge cases
    if x <= zero {
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
    // Number of iterations depends on precision needed
    let iterations = (T::frac_bits() / 2).clamp(8, 20);

    for _ in 0..iterations {
        if guess == zero {
            return zero;
        }

        let quotient = x.div(guess);
        let sum = guess.saturating_add(quotient);
        let new_guess = sum.saturating_mul(half);

        // Check for convergence
        let diff = if new_guess > guess {
            new_guess - guess
        } else {
            guess - new_guess
        };

        // Epsilon is approximately 2^(-frac_bits/2) for reasonable convergence.
        // This scales with the type's precision: smaller for higher precision types.
        // frac_bits() is at most 64, so this cast is safe.
        #[allow(clippy::cast_possible_wrap)]
        let epsilon_shift = 63i32 - (T::frac_bits() / 2) as i32;
        let epsilon_bits = if epsilon_shift >= 0 {
            1i64 << epsilon_shift
        } else {
            1i64 >> (-epsilon_shift)
        };
        let epsilon = T::from_i64_frac(epsilon_bits);
        if diff <= epsilon {
            return new_guess;
        }

        guess = new_guess;
    }

    guess
}
