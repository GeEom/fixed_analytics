//! The error functions and the Fresnel integrals of a real argument.

use crate::consts::FRAC_2_SQRT_PI;
use crate::cx::{Cx, erf as erf_work, faddeeva};
use crate::error::{Error, Result};
use crate::eval::{W, fixed, work};
use crate::real::Work;
use crate::traits::CordicNumber;

/// `erfc x` for `x ≥ 0`: `e^(−x²)·w(ix)`, which is real, past ½, else `1 − erf x`.
fn erfc_positive<V: Work>(x: V) -> V {
    if x <= V::one().ldexp(-1) {
        return V::one() - erf_work(Cx::real(x)).re;
    }
    crate::eval::exp_wide(-(x.square())) * faddeeva(Cx::new(V::zero(), x)).re
}

/// The error function `erf x`, in `[−1, 1]`.
#[must_use]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn erf<T: CordicNumber>(x: T) -> T {
    fixed(erf_work(Cx::real(work(x))).re)
}

/// The complementary error function `erfc x = 1 − erf x`, in `[0, 2]`.
#[must_use]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn erfc<T: CordicNumber>(x: T) -> T {
    let v = work(x);
    let c = erfc_positive(v.abs());
    fixed(if v.is_neg() {
        W::<T>::from_int(2) - c
    } else {
        c
    })
}

/// The Fresnel integrals `C(x) = ∫₀ˣ cos(πt²/2) dt` and `S(x) = ∫₀ˣ sin(πt²/2) dt`.
#[must_use]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn fresnel<T: CordicNumber>(x: T) -> (T, T) {
    let v = work(x);
    let scale = v / W::<T>::constant(&FRAC_2_SQRT_PI);
    let e = erf_work(Cx::new(scale, -scale));
    let half = Cx::new(W::<T>::one(), W::<T>::one()).ldexp(-1);
    let cs = half * e;
    (fixed(cs.re), fixed(cs.im))
}

/// The gamma function `Γ(x)` for `x > 0`, saturating as it grows.
/// # Errors
/// Returns `DomainError` if `x ≤ 0`.
#[must_use = "returns the gamma function which should be handled"]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn gamma<T: CordicNumber>(x: T) -> Result<T> {
    if x <= T::zero() {
        return Err(Error::domain("gamma", "positive value"));
    }
    Ok(fixed(crate::cx::gamma(work(x))))
}
