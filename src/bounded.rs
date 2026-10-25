//! Value types whose range is known from how they were built.
//!
//! Each holds a value checked by `new`, or the result of an identity that
//! cannot leave the range, such as [`NonNegative::one_plus_square`] for
//! `1 + x²`. Internal code then needs no domain checks, `unsafe` or `expect`.

use crate::traits::CordicNumber;

/// A value `≥ 0`, the input to an infallible square root.
#[derive(Clone, Copy, Debug)]
pub struct NonNegative<T>(T);

impl<T: CordicNumber> NonNegative<T> {
    /// Returns `Some` if `value ≥ 0`.
    #[inline]
    #[must_use]
    pub fn new(value: T) -> Option<Self> {
        (value >= T::zero()).then_some(Self(value))
    }

    /// Builds `1 + x²`, which is at least 1.
    #[inline]
    #[must_use]
    pub fn one_plus_square(x: T) -> Self {
        let x_sq = x.saturating_mul(x);
        Self(T::one().saturating_add(x_sq))
    }

    /// Builds `1 − x²`, which is non-negative for `|x| ≤ 1`.
    #[inline]
    #[must_use]
    pub fn one_minus_square(x: UnitInterval<T>) -> Self {
        let x_sq = x.0.saturating_mul(x.0);
        Self(T::one().saturating_sub(x_sq))
    }

    /// Builds `x² − 1`, which is non-negative for `x ≥ 1`.
    #[inline]
    #[must_use]
    pub fn square_minus_one(x: AtLeastOne<T>) -> Self {
        let x_sq = x.0.saturating_mul(x.0);
        Self(x_sq.saturating_sub(T::one()))
    }

    /// Returns the inner value.
    #[inline]
    #[must_use]
    pub const fn get(self) -> T {
        self.0
    }
}

/// A value in the closed interval `[-1, 1]`, as asin and acos take.
#[derive(Clone, Copy, Debug)]
pub struct UnitInterval<T>(T);

impl<T: CordicNumber> UnitInterval<T> {
    /// Returns `Some` if `value` is in `[-1, 1]`.
    #[inline]
    #[must_use]
    pub fn new(value: T) -> Option<Self> {
        let one = T::one();
        (value >= -one && value <= one).then_some(Self(value))
    }

    /// Builds `1 / x`, which is in `(0, 1]` for `x ≥ 1`.
    #[inline]
    #[must_use]
    pub fn from_reciprocal(x: AtLeastOne<T>) -> Self {
        Self(T::one().div(x.0))
    }

    /// Returns the inner value.
    #[inline]
    #[must_use]
    pub const fn get(self) -> T {
        self.0
    }
}

/// A value in the open interval `(-1, 1)`, as atanh takes.
#[derive(Clone, Copy, Debug)]
pub struct OpenUnitInterval<T>(T);

impl<T: CordicNumber> OpenUnitInterval<T> {
    /// Returns `Some` if `value` is in `(-1, 1)`.
    #[inline]
    #[must_use]
    pub fn new(value: T) -> Option<Self> {
        let one = T::one();
        (value > -one && value < one).then_some(Self(value))
    }

    /// Builds `x / √(1 + x²)` from the computed root; it is below 1 in magnitude as
    /// `√(1 + x²) > |x|`.
    #[inline]
    #[must_use]
    pub fn from_div_by_sqrt_one_plus_square(x: T, sqrt_one_plus_x_sq: T) -> Self {
        Self(x.div(sqrt_one_plus_x_sq))
    }

    /// Builds `√(x² − 1) / x` from the computed root; it is below 1 in magnitude as
    /// `√(x² − 1) < x`.
    #[inline]
    #[must_use]
    pub fn from_sqrt_square_minus_one_div(sqrt_x_sq_minus_one: T, x: AtLeastOne<T>) -> Self {
        Self(sqrt_x_sq_minus_one.div(x.0))
    }

    /// Builds `(x − 1) / (x + 1)`, which lies in `[-1/3, 1/3]` for `x` in `[0.5, 2]`.
    #[inline]
    #[must_use]
    pub(crate) fn from_normalized_ln_arg(x: NormalizedLnArg<T>) -> Self {
        let x_minus_1 = x.0 - T::one();
        let x_plus_1 = x.0 + T::one();
        Self(x_minus_1.div(x_plus_1))
    }

    /// Returns the inner value.
    #[inline]
    #[must_use]
    pub const fn get(self) -> T {
        self.0
    }
}

/// A value `≥ 1`, as acosh takes.
#[derive(Clone, Copy, Debug)]
pub struct AtLeastOne<T>(T);

impl<T: CordicNumber> AtLeastOne<T> {
    /// Returns `Some` if `value ≥ 1`.
    #[inline]
    #[must_use]
    pub fn new(value: T) -> Option<Self> {
        (value >= T::one()).then_some(Self(value))
    }

    /// Returns the inner value.
    #[inline]
    #[must_use]
    pub const fn get(self) -> T {
        self.0
    }
}

/// The ln argument after normalization, in `[0.5, 2]`.
#[derive(Clone, Copy, Debug)]
pub(crate) struct NormalizedLnArg<T>(T);

impl<T: CordicNumber> NormalizedLnArg<T> {
    /// Wraps the normalized value unchecked; `ln_positive` guarantees the range.
    #[inline]
    #[must_use]
    pub(crate) const fn from_normalized(value: T) -> Self {
        Self(value)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fixed::types::I16F16;

    #[test]
    fn non_negative_new() {
        assert!(NonNegative::new(I16F16::from_num(0)).is_some());
        assert!(NonNegative::new(I16F16::from_num(1)).is_some());
        assert!(NonNegative::new(I16F16::from_num(-1)).is_none());
    }

    #[test]
    fn non_negative_one_plus_square() {
        let x = I16F16::from_num(2);
        let nn = NonNegative::one_plus_square(x);
        assert_eq!(nn.get(), I16F16::from_num(5));
    }

    #[test]
    fn non_negative_one_minus_square() -> Result<(), &'static str> {
        let unit = UnitInterval::new(I16F16::from_num(0.5)).ok_or("out of range")?;
        let nn = NonNegative::one_minus_square(unit);
        let val: f32 = nn.get().to_num();
        assert!((val - 0.75).abs() < 0.01);
        Ok(())
    }

    #[test]
    fn non_negative_square_minus_one() -> Result<(), &'static str> {
        let at_least = AtLeastOne::new(I16F16::from_num(2)).ok_or("out of range")?;
        let nn = NonNegative::square_minus_one(at_least);
        assert_eq!(nn.get(), I16F16::from_num(3));
        Ok(())
    }

    #[test]
    fn unit_interval_new() {
        assert!(UnitInterval::new(I16F16::from_num(0)).is_some());
        assert!(UnitInterval::new(I16F16::from_num(1)).is_some());
        assert!(UnitInterval::new(I16F16::from_num(-1)).is_some());
        assert!(UnitInterval::new(I16F16::from_num(1.1)).is_none());
        assert!(UnitInterval::new(I16F16::from_num(-1.1)).is_none());
    }

    #[test]
    fn unit_interval_from_reciprocal() -> Result<(), &'static str> {
        let at_least = AtLeastOne::new(I16F16::from_num(4)).ok_or("out of range")?;
        assert_eq!(
            UnitInterval::from_reciprocal(at_least).get(),
            I16F16::from_num(0.25)
        );
        Ok(())
    }

    #[test]
    fn unit_interval_get() -> Result<(), &'static str> {
        let unit = UnitInterval::new(I16F16::from_num(0.5)).ok_or("out of range")?;
        assert_eq!(unit.get(), I16F16::from_num(0.5));
        Ok(())
    }

    #[test]
    fn open_unit_interval_new() {
        assert!(OpenUnitInterval::new(I16F16::from_num(0)).is_some());
        assert!(OpenUnitInterval::new(I16F16::from_num(0.5)).is_some());
        assert!(OpenUnitInterval::new(I16F16::from_num(1)).is_none());
        assert!(OpenUnitInterval::new(I16F16::from_num(-1)).is_none());
    }

    #[test]
    fn open_unit_interval_get() -> Result<(), &'static str> {
        let open = OpenUnitInterval::new(I16F16::from_num(0.5)).ok_or("out of range")?;
        assert_eq!(open.get(), I16F16::from_num(0.5));
        Ok(())
    }

    #[test]
    fn open_unit_interval_from_div() {
        let x = I16F16::from_num(1);
        let sqrt_1_plus_x_sq = I16F16::SQRT_2;
        let open = OpenUnitInterval::from_div_by_sqrt_one_plus_square(x, sqrt_1_plus_x_sq);
        let val: f32 = open.get().to_num();
        assert!((val - core::f32::consts::FRAC_1_SQRT_2).abs() < 0.01);
    }

    #[test]
    fn open_unit_interval_from_sqrt_div() -> Result<(), &'static str> {
        let at_least = AtLeastOne::new(I16F16::from_num(2)).ok_or("out of range")?;
        let sqrt_x_sq_minus_one = I16F16::from_num(1.732);
        let open = OpenUnitInterval::from_sqrt_square_minus_one_div(sqrt_x_sq_minus_one, at_least);
        let val: f32 = open.get().to_num();
        let expected = 0.866_f32;
        assert!((val - expected).abs() < 0.01);
        Ok(())
    }

    #[test]
    fn at_least_one_new() {
        assert!(AtLeastOne::new(I16F16::from_num(1)).is_some());
        assert!(AtLeastOne::new(I16F16::from_num(2)).is_some());
        assert!(AtLeastOne::new(I16F16::from_num(0.9)).is_none());
    }

    #[test]
    fn at_least_one_get() -> Result<(), &'static str> {
        let at_least = AtLeastOne::new(I16F16::from_num(2)).ok_or("out of range")?;
        assert_eq!(at_least.get(), I16F16::from_num(2));
        Ok(())
    }

    #[test]
    fn open_unit_interval_from_normalized_ln_arg() {
        let norm = NormalizedLnArg::from_normalized(I16F16::from_num(1.5));
        let open = OpenUnitInterval::from_normalized_ln_arg(norm);
        let val: f32 = open.get().to_num();
        assert!((val - 0.2).abs() < 0.01);
    }
}
