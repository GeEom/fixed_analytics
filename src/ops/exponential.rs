#![doc = include_str!("../../docs/exponential.md")]

use crate::consts::{FRAC_1_LN_10, LN_2};
use crate::error::{Error, Result};
use crate::eval::{
    HUGE, W, beyond, exp as exp_work, expm1 as expm1_work, fixed, ln as ln_work,
    log1p as log1p_work, log1p_minus, log1p_minus_ratio, log1p_ratio as log1p_ratio_work,
    log2 as log2_work, phi, work,
};
use crate::real::Work;
use crate::traits::CordicNumber;

/// `e^v`, or for `|v|` past [`HUGE`] a value every type saturates or rounds
/// to zero from.
fn exp_wide<T: CordicNumber>(v: W<T>) -> W<T> {
    if v.abs() > W::<T>::from_int(HUGE) {
        return if v.is_neg() { W::<T>::zero() } else { beyond() };
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

/// Power `base^exponent`, as `exp(exponent · ln|base|)`, saturating like
/// [`exp`]; a negative base takes its sign from the exponent's parity.
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

/// `e^x − 1`, to full relative precision near 0, saturating like [`exp`].
#[must_use]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn expm1<T: CordicNumber>(x: T) -> T {
    let v = work(x);
    if v > W::<T>::from_int(HUGE) {
        return T::max_value();
    }
    if v < W::<T>::from_int(-HUGE) {
        return -T::one();
    }
    fixed(expm1_work(v))
}

/// `(e^x − 1)/x`, which is 1 at 0, saturating for large `x`.
#[must_use]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn exprel<T: CordicNumber>(x: T) -> T {
    fixed(phi(1, work(x)))
}

/// `(e^x − 1 − x)/x²`, which is ½ at 0, saturating for large `x`.
#[must_use]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn exprel2<T: CordicNumber>(x: T) -> T {
    fixed(phi(2, work(x)))
}

/// `[E₁, E₂, E₃, e^(−x)]`, with `Eₖ(x) = ∫₀¹ s^(k−1)·e^(−xs) ds`.
#[must_use]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn exp_moments<T: CordicNumber>(x: T) -> [T; 4] {
    let v = work(x);
    if v > W::<T>::from_int(HUGE) {
        let one = W::<T>::one();
        let e1 = one / v;
        return [
            fixed(e1),
            fixed(e1 / v),
            fixed((e1 / v.square()).ldexp(1)),
            T::zero(),
        ];
    }
    if v < W::<T>::from_int(-HUGE) {
        return [
            T::max_value(),
            T::max_value(),
            T::max_value(),
            T::max_value(),
        ];
    }
    let decay = exp_work(-v);
    [
        fixed(decay * phi(1, v)),
        fixed(decay * phi(2, v)),
        fixed((decay * phi(3, v)).ldexp(1)),
        fixed(decay),
    ]
}

/// `ln(1 + v)`, to full relative precision near 0. Domain: `v > −1`.
/// # Errors
/// Returns `DomainError` if `v ≤ −1`.
#[must_use = "returns the logarithm which should be handled"]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn log1p<T: CordicNumber>(v: T) -> Result<T> {
    above_minus_one(v, "log1p").map(|w| fixed(log1p_work(w)))
}

/// `ln(1 + z)/z`, which is 1 at 0. Domain: `z > −1`.
/// # Errors
/// Returns `DomainError` if `z ≤ −1`.
#[must_use = "returns the ratio which should be handled"]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn log1p_ratio<T: CordicNumber>(z: T) -> Result<T> {
    above_minus_one(z, "log1p_ratio").map(|w| fixed(log1p_ratio_work(w)))
}

/// `ln(1 + v) − v`, without cancelling near 0. Domain: `v > −1`.
/// # Errors
/// Returns `DomainError` if `v ≤ −1`.
#[must_use = "returns the residual which should be handled"]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn log1p_residual<T: CordicNumber>(v: T) -> Result<T> {
    above_minus_one(v, "log1p_residual").map(|w| fixed(log1p_minus(w)))
}

/// `(ln(1 + u) − u)/u²`, which is −½ at 0. Domain: `u > −1`.
/// # Errors
/// Returns `DomainError` if `u ≤ −1`.
#[must_use = "returns the ratio which should be handled"]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn log1p_residual_ratio<T: CordicNumber>(u: T) -> Result<T> {
    above_minus_one(u, "log1p_residual_ratio").map(|w| fixed(log1p_minus_ratio(w)))
}

fn above_minus_one<T: CordicNumber>(v: T, name: &'static str) -> Result<W<T>> {
    if v <= -T::one() {
        return Err(Error::domain(name, "value > -1"));
    }
    Ok(work(v))
}
