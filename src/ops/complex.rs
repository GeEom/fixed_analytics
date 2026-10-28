//! Complex numbers over a fixed-point type, each part rounded once.

use crate::consts::INV_FACT;
use crate::cx::{Cx, atan2};
use crate::error::{Error, Result};
use crate::eval::{W, by_width, fixed, work};
use crate::real::Work;
use crate::traits::CordicNumber;
use core::ops::{Add, Div, Mul, Neg, Sub};

/// A complex number `re + i·im`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Complex<T> {
    /// The real part.
    pub re: T,
    /// The imaginary part.
    pub im: T,
}

impl<T: CordicNumber> Complex<T> {
    /// `re + i·im`.
    #[must_use]
    pub const fn new(re: T, im: T) -> Self {
        Self { re, im }
    }

    fn work(self) -> Cx<W<T>> {
        Cx::new(work(self.re), work(self.im))
    }

    fn round(z: Cx<W<T>>) -> Self {
        Self::new(fixed(z.re), fixed(z.im))
    }

    /// The conjugate `re − i·im`.
    #[must_use]
    pub fn conj(self) -> Self {
        Self::new(self.re, T::zero().saturating_sub(self.im))
    }

    /// The modulus `|z|`, rounded once.
    #[must_use]
    pub fn norm(self) -> T {
        fixed(self.work().abs())
    }

    /// The argument in `(−π, π]`, 0 at the origin.
    #[must_use]
    pub fn arg(self) -> T {
        fixed(atan2(work(self.im), work(self.re)))
    }

    /// `e^z`.
    #[must_use]
    #[cfg_attr(
        all(feature = "verify-no-panic", not(debug_assertions)),
        no_panic::no_panic
    )]
    pub fn exp(self) -> Self {
        Self::round(self.work().exp())
    }

    /// `(e^z − 1)/z`, which is 1 at 0.
    #[must_use]
    #[cfg_attr(
        all(feature = "verify-no-panic", not(debug_assertions)),
        no_panic::no_panic
    )]
    pub fn exprel(self) -> Self {
        Self::round(phi(1, self.work()))
    }

    /// `(e^z − 1 − z)/z²`, which is ½ at 0.
    #[must_use]
    #[cfg_attr(
        all(feature = "verify-no-panic", not(debug_assertions)),
        no_panic::no_panic
    )]
    pub fn exprel2(self) -> Self {
        Self::round(phi(2, self.work()))
    }

    /// The Faddeeva function `w(z) = e^(−z²)·erfc(−iz)`.
    #[must_use]
    #[cfg_attr(
        all(feature = "verify-no-panic", not(debug_assertions)),
        no_panic::no_panic
    )]
    pub fn faddeeva(self) -> Self {
        Self::round(crate::cx::faddeeva(self.work()))
    }

    /// `E₁(z) = ∫_z^∞ e^(−t)/t dt` on the principal branch, from above on the negative axis.
    /// # Errors
    /// Returns `DomainError` at `z = 0`.
    #[must_use = "returns the integral which should be handled"]
    #[cfg_attr(
        all(feature = "verify-no-panic", not(debug_assertions)),
        no_panic::no_panic
    )]
    pub fn e1(self) -> Result<Self> {
        let z = self.work();
        if z.is_zero() {
            return Err(Error::domain("e1", "non-zero argument"));
        }
        Ok(Self::round(crate::cx::e1(z)))
    }

    /// `Γ(a, z) = ∫_z^∞ t^(a−1)·e^(−t) dt` for real `a` in `(0, 128]`, principal branch.
    /// # Errors
    /// Returns `DomainError` for `a` outside `(0, 128]`.
    #[must_use = "returns the integral which should be handled"]
    #[cfg_attr(
        all(feature = "verify-no-panic", not(debug_assertions)),
        no_panic::no_panic
    )]
    pub fn gamma_upper(self, a: T) -> Result<Self> {
        let wa = work(a);
        if wa <= W::<T>::zero() || wa > W::<T>::from_int(128) {
            return Err(Error::domain("gamma_upper", "a in (0, 128]"));
        }
        Ok(Self::round(crate::cx::gamma_upper(wa, self.work())))
    }

    /// The error function `erf z`.
    #[must_use]
    #[cfg_attr(
        all(feature = "verify-no-panic", not(debug_assertions)),
        no_panic::no_panic
    )]
    pub fn erf(self) -> Self {
        Self::round(crate::cx::erf(self.work()))
    }
}

/// `φₖ(z) = (e^z − Σ_{n<k} zⁿ/n!)/z^k` for `k` of 1 or 2: the series near 0.
fn phi<V: Work>(k: u32, z: Cx<V>) -> Cx<V> {
    if z.abs() <= V::one().ldexp(-1) {
        let n = by_width::<V>([18, 30, 39]);
        let mut p = Cx::zero();
        for j in (0..=n).rev() {
            let coef = INV_FACT
                .get((k + j) as usize)
                .map_or_else(V::zero, V::constant);
            p = Cx::real(coef) + z * p;
        }
        return p;
    }
    let mut num = z.exp_m1();
    let mut den = z;
    if k >= 2 {
        num = num - z;
        den = den * z;
    }
    num / den
}

impl<T: CordicNumber> Add for Complex<T> {
    type Output = Self;
    fn add(self, o: Self) -> Self {
        Self::new(self.re.saturating_add(o.re), self.im.saturating_add(o.im))
    }
}

impl<T: CordicNumber> Sub for Complex<T> {
    type Output = Self;
    fn sub(self, o: Self) -> Self {
        Self::new(self.re.saturating_sub(o.re), self.im.saturating_sub(o.im))
    }
}

impl<T: CordicNumber> Neg for Complex<T> {
    type Output = Self;
    fn neg(self) -> Self {
        Self::new(T::zero(), T::zero()) - self
    }
}

impl<T: CordicNumber> Mul for Complex<T> {
    type Output = Self;
    /// Each part formed exactly and rounded once.
    fn mul(self, o: Self) -> Self {
        let (xr, xi, yr, yi) = (work(self.re), work(self.im), work(o.re), work(o.im));
        Self::new(
            fixed(W::<T>::dot2(xr, yr, -xi, yi)),
            fixed(W::<T>::dot2(xr, yi, xi, yr)),
        )
    }
}

impl<T: CordicNumber> Div for Complex<T> {
    type Output = Self;
    /// From exact sums of products, so no `w ≠ 0` overflows; a zero divisor
    /// saturates each part by its sign.
    fn div(self, o: Self) -> Self {
        let (xr, xi, yr, yi) = (work(self.re), work(self.im), work(o.re), work(o.im));
        if yr.is_zero() && yi.is_zero() {
            let edge = |x: T| {
                if x == T::zero() {
                    x
                } else if x.is_negative() {
                    T::min_value()
                } else {
                    T::max_value()
                }
            };
            return Self::new(edge(self.re), edge(self.im));
        }
        let norm = W::<T>::dot2(yr, yr, yi, yi);
        Self::new(
            fixed(W::<T>::dot2(xr, yr, xi, yi) / norm),
            fixed(W::<T>::dot2(xi, yr, -xr, yi) / norm),
        )
    }
}
