//! Bounded value types that encode mathematical invariants at the type level.
//!
//! These types provide compile-time guarantees about value ranges, eliminating
//! the need for runtime checks in internal computations where the domain is
//! already validated.
//!
//! # Types
//!
//! - [`NonNegative<T>`]: Values >= 0 (for sqrt inputs)
//! - [`UnitInterval<T>`]: Values in [-1, 1] (for asin/acos inputs)
//! - [`OpenUnitInterval<T>`]: Values in (-1, 1) (for atanh inputs)
//!
//! # Design Philosophy
//!
//! Rather than using `unsafe` or `expect`, these types encode the mathematical
//! relationships between operations. For example:
//!
//! - `1 + x^2` is always >= 1, so `NonNegative::one_plus_square(x)` is infallible
//! - If `|x| <= 1`, then `1 - x^2` is in [0, 1], so `NonNegative::one_minus_square(unit_x)` is infallible
//! - `x / sqrt(1 + x^2)` is always in (-1, 1), so `OpenUnitInterval::x_div_sqrt_one_plus_x_sq(x)` is infallible

use crate::traits::CordicNumber;

/// A value guaranteed to be non-negative (>= 0).
///
/// This type enables infallible sqrt operations by encoding the non-negativity
/// constraint at the type level.
///
/// # Construction
///
/// - [`NonNegative::new`]: Checked construction, returns `Option`
/// - [`NonNegative::one_plus_square`]: From `1 + x^2`, always valid
/// - [`NonNegative::one_minus_square`]: From `1 - x^2` where `|x| <= 1`, always valid
/// - [`NonNegative::square_minus_one`]: From `x^2 - 1` where `|x| >= 1`, always valid
#[derive(Clone, Copy, Debug)]
pub struct NonNegative<T>(T);

impl<T: CordicNumber> NonNegative<T> {
    /// Creates a new `NonNegative` value if the input is >= 0.
    ///
    /// Returns `None` if the value is negative.
    #[inline]
    #[must_use]
    pub fn new(value: T) -> Option<Self> {
        (value >= T::zero()).then_some(Self(value))
    }

    /// Constructs from `1 + x^2`, which is always >= 1.
    ///
    /// This is mathematically infallible: for any real `x`, `1 + x^2 >= 1`.
    #[inline]
    #[must_use]
    pub fn one_plus_square(x: T) -> Self {
        let x_sq = x.saturating_mul(x);
        Self(T::one().saturating_add(x_sq))
    }

    /// Constructs from `1 - x^2` where `|x| <= 1`.
    ///
    /// This is mathematically infallible: if `|x| <= 1`, then `x^2 <= 1`,
    /// so `1 - x^2 >= 0`.
    #[inline]
    #[must_use]
    pub fn one_minus_square(x: UnitInterval<T>) -> Self {
        let x_sq = x.0.saturating_mul(x.0);
        Self(T::one().saturating_sub(x_sq))
    }

    /// Constructs from `x^2 - 1` where `|x| >= 1`.
    ///
    /// This is mathematically infallible: if `|x| >= 1`, then `x^2 >= 1`,
    /// so `x^2 - 1 >= 0`.
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

/// A value guaranteed to be in the closed interval [-1, 1].
///
/// This type is used for inputs to functions like asin and acos that
/// require their argument to be in this range.
#[derive(Clone, Copy, Debug)]
pub struct UnitInterval<T>(T);

impl<T: CordicNumber> UnitInterval<T> {
    /// Creates a new `UnitInterval` value if the input is in [-1, 1].
    ///
    /// Returns `None` if the value is outside the interval.
    #[inline]
    #[must_use]
    pub fn new(value: T) -> Option<Self> {
        let one = T::one();
        (value >= -one && value <= one).then_some(Self(value))
    }

    /// Returns the inner value.
    #[inline]
    #[must_use]
    pub const fn get(self) -> T {
        self.0
    }
}

/// A value guaranteed to be in the open interval (-1, 1).
///
/// This type is used for inputs to atanh, which requires strict inequality.
#[derive(Clone, Copy, Debug)]
pub struct OpenUnitInterval<T>(T);

impl<T: CordicNumber> OpenUnitInterval<T> {
    /// Creates a new `OpenUnitInterval` value if the input is in (-1, 1).
    ///
    /// Returns `None` if the value is outside the interval or on the boundary.
    #[inline]
    #[must_use]
    pub fn new(value: T) -> Option<Self> {
        let one = T::one();
        (value > -one && value < one).then_some(Self(value))
    }

    /// Constructs from `x / sqrt(1 + x^2)`, which is always in (-1, 1).
    ///
    /// This is mathematically infallible: for any real `x`,
    /// `|x / sqrt(1 + x^2)| < 1` because `sqrt(1 + x^2) > |x|`.
    ///
    /// Note: Requires the sqrt to be computed first.
    #[inline]
    #[must_use]
    pub fn from_div_by_sqrt_one_plus_square(x: T, sqrt_one_plus_x_sq: T) -> Self {
        Self(x.div(sqrt_one_plus_x_sq))
    }

    /// Constructs from `sqrt(x^2 - 1) / x` where `|x| >= 1`.
    ///
    /// This is mathematically infallible: for `|x| >= 1`,
    /// `|sqrt(x^2 - 1) / x| < 1` because `sqrt(x^2 - 1) < |x|` for `|x| > 1`.
    ///
    /// Note: Requires the sqrt to be computed first.
    #[inline]
    #[must_use]
    pub fn from_sqrt_square_minus_one_div(sqrt_x_sq_minus_one: T, x: AtLeastOne<T>) -> Self {
        Self(sqrt_x_sq_minus_one.div(x.0))
    }

    /// Constructs from `(x - 1) / (x + 1)` where `x` is in [0.5, 2].
    ///
    /// This is mathematically infallible: for `x` in [0.5, 2],
    /// `(x - 1) / (x + 1)` is in `(-1/3, 1/3)` which is a subset of `(-1, 1)`.
    #[inline]
    #[must_use]
    pub fn from_normalized_ln_arg(x: NormalizedLnArg<T>) -> Self {
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

/// A value guaranteed to be >= 1 (or <= -1 for the absolute value).
///
/// This type is used for inputs to acosh which requires x >= 1.
#[derive(Clone, Copy, Debug)]
pub struct AtLeastOne<T>(T);

impl<T: CordicNumber> AtLeastOne<T> {
    /// Creates a new `AtLeastOne` value if the input is >= 1.
    ///
    /// Returns `None` if the value is less than 1.
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

/// A value guaranteed to be in [0.5, 2], used for ln argument normalization.
///
/// After normalizing the input for ln computation, the value is always
/// in this range, which guarantees that `(x-1)/(x+1)` is in `(-1/3, 1/3)`.
#[derive(Clone, Copy, Debug)]
pub struct NormalizedLnArg<T>(T);

impl<T: CordicNumber> NormalizedLnArg<T> {
    /// Creates a new `NormalizedLnArg` from the normalization loop result.
    ///
    /// The ln function's normalization loop guarantees the result is in [0.5, 2].
    /// This constructor trusts that invariant (used only in ln implementation).
    #[inline]
    #[must_use]
    pub(crate) const fn from_normalized(value: T) -> Self {
        Self(value)
    }

    /// Returns the inner value.
    #[inline]
    #[must_use]
    pub const fn get(self) -> T {
        self.0
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used, reason = "test code uses unwrap for conciseness")]
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
    fn non_negative_one_minus_square() {
        let unit = UnitInterval::new(I16F16::from_num(0.5)).unwrap();
        let nn = NonNegative::one_minus_square(unit);
        let val: f32 = nn.get().to_num();
        assert!((val - 0.75).abs() < 0.01);
    }

    #[test]
    fn non_negative_square_minus_one() {
        let at_least = AtLeastOne::new(I16F16::from_num(2)).unwrap();
        let nn = NonNegative::square_minus_one(at_least);
        assert_eq!(nn.get(), I16F16::from_num(3));
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
    fn unit_interval_get() {
        let unit = UnitInterval::new(I16F16::from_num(0.5)).unwrap();
        assert_eq!(unit.get(), I16F16::from_num(0.5));
    }

    #[test]
    fn open_unit_interval_new() {
        assert!(OpenUnitInterval::new(I16F16::from_num(0)).is_some());
        assert!(OpenUnitInterval::new(I16F16::from_num(0.5)).is_some());
        assert!(OpenUnitInterval::new(I16F16::from_num(1)).is_none());
        assert!(OpenUnitInterval::new(I16F16::from_num(-1)).is_none());
    }

    #[test]
    fn open_unit_interval_get() {
        let open = OpenUnitInterval::new(I16F16::from_num(0.5)).unwrap();
        assert_eq!(open.get(), I16F16::from_num(0.5));
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
    fn open_unit_interval_from_sqrt_div() {
        let at_least = AtLeastOne::new(I16F16::from_num(2)).unwrap();
        #[allow(
            clippy::approx_constant,
            reason = "testing with known approximation of sqrt(3)"
        )]
        let sqrt_x_sq_minus_one = I16F16::from_num(1.732);
        let open = OpenUnitInterval::from_sqrt_square_minus_one_div(sqrt_x_sq_minus_one, at_least);
        let val: f32 = open.get().to_num();
        #[allow(
            clippy::approx_constant,
            reason = "testing with known approximation of sqrt(3)/2"
        )]
        let expected = 0.866_f32;
        assert!((val - expected).abs() < 0.01);
    }

    #[test]
    fn at_least_one_new() {
        assert!(AtLeastOne::new(I16F16::from_num(1)).is_some());
        assert!(AtLeastOne::new(I16F16::from_num(2)).is_some());
        assert!(AtLeastOne::new(I16F16::from_num(0.9)).is_none());
    }

    #[test]
    fn at_least_one_get() {
        let at_least = AtLeastOne::new(I16F16::from_num(2)).unwrap();
        assert_eq!(at_least.get(), I16F16::from_num(2));
    }

    #[test]
    fn normalized_ln_arg_get() {
        let norm = NormalizedLnArg::from_normalized(I16F16::from_num(1.5));
        assert_eq!(norm.get(), I16F16::from_num(1.5));
    }

    #[test]
    fn open_unit_interval_from_normalized_ln_arg() {
        let norm = NormalizedLnArg::from_normalized(I16F16::from_num(1.5));
        let open = OpenUnitInterval::from_normalized_ln_arg(norm);
        let val: f32 = open.get().to_num();
        assert!((val - 0.2).abs() < 0.01);
    }
}
