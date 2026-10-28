//! Complex numbers over a fixed-point type, each part rounded once.

use crate::cx::{Cx, atan2, phi};
use crate::error::{Error, Result};
use crate::eval::{W, atan_ratio, exp_wide, fixed, sin_cos, sin_cos_octant, tiered, work};
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

    fn work_in<V: Work>(self, as_work: fn(T) -> V) -> Cx<V> {
        Cx::new(as_work(self.re), as_work(self.im))
    }

    fn round(z: Cx<W<T>>) -> Self {
        Self::new(fixed(z.re), fixed(z.im))
    }

    /// The conjugate `re − i·im`.
    #[must_use]
    #[cfg_attr(
        all(feature = "verify-no-panic", not(debug_assertions)),
        no_panic::no_panic
    )]
    pub fn conj(self) -> Self {
        Self::new(self.re, T::zero().saturating_sub(self.im))
    }

    /// The modulus `|z|`, rounded once.
    #[must_use]
    #[cfg_attr(
        all(feature = "verify-no-panic", not(debug_assertions)),
        no_panic::no_panic
    )]
    pub fn norm(self) -> T {
        fixed(self.work().abs())
    }

    /// The argument in `(−π, π]`, 0 at the origin.
    #[must_use]
    #[cfg_attr(
        all(feature = "verify-no-panic", not(debug_assertions)),
        no_panic::no_panic
    )]
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
        let (re, im) = tiered!(T, |v| self.work_in(v).exp().parts());
        Self::new(re, im)
    }

    /// `(e^z − 1)/z`, which is 1 at 0.
    #[must_use]
    #[cfg_attr(
        all(feature = "verify-no-panic", not(debug_assertions)),
        no_panic::no_panic
    )]
    pub fn exprel(self) -> Self {
        let (re, im) = tiered!(T, |v| phi(1, self.work_in(v)).parts());
        Self::new(re, im)
    }

    /// `(e^z − 1 − z)/z²`, which is ½ at 0.
    #[must_use]
    #[cfg_attr(
        all(feature = "verify-no-panic", not(debug_assertions)),
        no_panic::no_panic
    )]
    pub fn exprel2(self) -> Self {
        let (re, im) = tiered!(T, |v| phi(2, self.work_in(v)).parts());
        Self::new(re, im)
    }

    /// The Faddeeva function `w(z) = e^(−z²)·erfc(−iz)`.
    #[must_use]
    #[cfg_attr(
        all(feature = "verify-no-panic", not(debug_assertions)),
        no_panic::no_panic
    )]
    pub fn faddeeva(self) -> Self {
        let (re, im) = tiered!(T, 48, |v| crate::cx::faddeeva(self.work_in(v)).parts());
        Self::new(re, im)
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
        let (re, im) = tiered!(T, 24, |v| crate::cx::e1(self.work_in(v)).parts());
        Ok(Self::new(re, im))
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
        let (re, im) = tiered!(T, 24, |v| crate::cx::gamma_upper(v(a), self.work_in(v))
            .parts());
        Ok(Self::new(re, im))
    }

    /// The principal logarithm `ln|z| + i·arg z`, `arg z` in `(−π, π]`.
    /// # Errors
    /// Returns `DomainError` at `z = 0`.
    #[must_use = "returns the logarithm which should be handled"]
    #[cfg_attr(
        all(feature = "verify-no-panic", not(debug_assertions)),
        no_panic::no_panic
    )]
    pub fn ln(self) -> Result<Self> {
        let z = self.work();
        if z.is_zero() {
            return Err(Error::domain("ln", "non-zero argument"));
        }
        let (re, im) = tiered!(T, |v| self.work_in(v).ln().parts());
        Ok(Self::new(re, im))
    }

    /// `z^a = e^(a·ln z)` for real `a`, principal branch; the phase's multiple of `π/4` reduced
    /// exactly, so parts that vanish on an axis or diagonal stay 0, and near one keep their bits.
    /// # Errors
    /// Returns `DomainError` at `z = 0` for `a < 0`.
    #[must_use = "returns the power which should be handled"]
    #[cfg_attr(
        all(feature = "verify-no-panic", not(debug_assertions)),
        no_panic::no_panic
    )]
    pub fn pow(self, power: T) -> Result<Self> {
        let z = self.work();
        if z.is_zero() {
            return match power.partial_cmp(&T::zero()) {
                Some(core::cmp::Ordering::Greater) => Ok(Self::new(T::zero(), T::zero())),
                Some(core::cmp::Ordering::Equal) => Ok(Self::new(T::one(), T::zero())),
                _ => Err(Error::domain("pow", "non-negative exponent at zero")),
            };
        }
        let (x, y) = (self.re.raw_i128(), self.im.raw_i128());
        let (ax, ay) = (x.unsigned_abs(), y.unsigned_abs());
        let (wx, wy) = (z.re.abs(), z.im.abs());
        // atan2(|y|, |x|) = j·π/4 + ψ about the nearest axis or diagonal, ψ to
        // full relative precision; past π/2, π − that.
        let (line, psi) = if ax == ay {
            (1u8, W::<T>::zero())
        } else if wy.mul_int(5) <= wx.mul_int(2) {
            (0, atan_ratio(wy, wx))
        } else if wx.mul_int(5) <= wy.mul_int(2) {
            (2, -atan_ratio(wx, wy))
        } else {
            let gap = wy - wx;
            (1, atan_ratio(gap.abs(), wx + wy).with_sign(gap.is_neg()))
        };
        let (eighths, eta) = if x < 0 { (4 - line, -psi) } else { (line, psi) };
        // a·k/8 turns, from |a|·k reduced modulo 2^(frac + 3).
        let point = T::frac_bits() + 3;
        let mask = u128::MAX >> (128 - point);
        let turns = power.raw_i128().unsigned_abs().wrapping_mul(eighths.into()) & mask;
        let (oct, frac) = W::<T>::octant(turns, point, false);
        let (s1, c1) = sin_cos_octant(oct, frac);
        let s1 = if power.is_negative() { -s1 } else { s1 };
        let (s2, c2) = sin_cos(work(power) * eta);
        let sine = s1 * c2 + c1 * s2;
        let cosine = c1 * c2 - s1 * s2;
        let sine = if y < 0 { -sine } else { sine };
        let modulus = exp_wide(work(power) * z.ln().re);
        Ok(Self::round(Cx::new(cosine, sine).scale(modulus)))
    }

    /// [`exprel`](Self::exprel) and [`exprel2`](Self::exprel2) weighted by `wₜ(s) = ((1 − εt)/(1 − εs))^p/(1 − εs)`:
    /// `K₁ = ∫₀¹ e^(zs)·w₁(s) ds` and `K₂ = ∫₀¹ ∫₀ᵗ e^(zs)·wₜ(s) ds dt`, precise as `ε → 0` for any `p`.
    /// # Errors
    /// Returns `DomainError` unless `ε` is in `[0, 1)` and `p ≥ 0`.
    #[must_use = "returns the integrals which should be handled"]
    #[cfg_attr(
        all(feature = "verify-no-panic", not(debug_assertions)),
        no_panic::no_panic
    )]
    pub fn exprel_weighted(self, eps: T, p: T) -> Result<(Self, Self)> {
        let zero = T::zero();
        if eps < zero || eps >= T::one() || p < zero {
            return Err(Error::domain("exprel_weighted", "eps in [0, 1) and p >= 0"));
        }
        let ((r1, i1), (r2, i2)) = tiered!(T, 16, |v| {
            let z = self.work_in(v);
            let (k1, k2) = if eps == zero {
                (phi(1, z), phi(2, z))
            } else {
                crate::cx::exprel_weighted(z, v(eps), v(p))
            };
            (k1.parts(), k2.parts())
        });
        Ok((Self::new(r1, i1), Self::new(r2, i2)))
    }

    /// The error function `erf z`.
    #[must_use]
    #[cfg_attr(
        all(feature = "verify-no-panic", not(debug_assertions)),
        no_panic::no_panic
    )]
    pub fn erf(self) -> Self {
        let (re, im) = tiered!(T, 48, |v| crate::cx::erf(self.work_in(v)).parts());
        Self::new(re, im)
    }
}

impl<T: CordicNumber> Complex<T> {
    /// `self + o`, each part saturating; the `+` operator.
    #[must_use]
    #[cfg_attr(
        all(feature = "verify-no-panic", not(debug_assertions)),
        no_panic::no_panic
    )]
    pub fn saturating_add(self, o: Self) -> Self {
        Self::new(self.re.saturating_add(o.re), self.im.saturating_add(o.im))
    }

    /// `self − o`, each part saturating; the `−` operator.
    #[must_use]
    #[cfg_attr(
        all(feature = "verify-no-panic", not(debug_assertions)),
        no_panic::no_panic
    )]
    pub fn saturating_sub(self, o: Self) -> Self {
        Self::new(self.re.saturating_sub(o.re), self.im.saturating_sub(o.im))
    }

    /// `−self`, each part saturating; unary `−`.
    #[must_use]
    #[cfg_attr(
        all(feature = "verify-no-panic", not(debug_assertions)),
        no_panic::no_panic
    )]
    pub fn saturating_neg(self) -> Self {
        Self::new(T::zero(), T::zero()).saturating_sub(self)
    }

    /// `self·o`, each part formed exactly and rounded once; the `*` operator.
    #[must_use]
    #[cfg_attr(
        all(feature = "verify-no-panic", not(debug_assertions)),
        no_panic::no_panic
    )]
    pub fn saturating_mul(self, o: Self) -> Self {
        let (xr, xi, yr, yi) = (work(self.re), work(self.im), work(o.re), work(o.im));
        Self::new(
            fixed(W::<T>::dot2(xr, yr, -xi, yi)),
            fixed(W::<T>::dot2(xr, yi, xi, yr)),
        )
    }

    /// `self/o` from exact sums of products, so no `o ≠ 0` overflows, else saturating; `/`.
    #[must_use]
    #[cfg_attr(
        all(feature = "verify-no-panic", not(debug_assertions)),
        no_panic::no_panic
    )]
    pub fn saturating_div(self, o: Self) -> Self {
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

impl<T: CordicNumber> Add for Complex<T> {
    type Output = Self;
    fn add(self, o: Self) -> Self {
        self.saturating_add(o)
    }
}

impl<T: CordicNumber> Sub for Complex<T> {
    type Output = Self;
    fn sub(self, o: Self) -> Self {
        self.saturating_sub(o)
    }
}

impl<T: CordicNumber> Neg for Complex<T> {
    type Output = Self;
    fn neg(self) -> Self {
        self.saturating_neg()
    }
}

impl<T: CordicNumber> Mul for Complex<T> {
    type Output = Self;
    fn mul(self, o: Self) -> Self {
        self.saturating_mul(o)
    }
}

impl<T: CordicNumber> Div for Complex<T> {
    type Output = Self;
    fn div(self, o: Self) -> Self {
        self.saturating_div(o)
    }
}
