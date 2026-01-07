//! MPFR reference implementations.

use crate::REFERENCE_PRECISION;
use rug::Float;

pub mod circular {
    use super::*;
    pub fn sin(x: &Float) -> Float { x.clone().sin() }
    pub fn cos(x: &Float) -> Float { x.clone().cos() }
    pub fn tan(x: &Float) -> Float { x.clone().tan() }
    pub fn asin(x: &Float) -> Float { x.clone().asin() }
    pub fn acos(x: &Float) -> Float { x.clone().acos() }
    pub fn atan(x: &Float) -> Float { x.clone().atan() }
}

pub mod hyperbolic {
    use super::*;
    pub fn sinh(x: &Float) -> Float { x.clone().sinh() }
    pub fn cosh(x: &Float) -> Float { x.clone().cosh() }
    pub fn tanh(x: &Float) -> Float { x.clone().tanh() }
    pub fn asinh(x: &Float) -> Float { x.clone().asinh() }
    pub fn acosh(x: &Float) -> Float { x.clone().acosh() }
    pub fn atanh(x: &Float) -> Float { x.clone().atanh() }
    pub fn coth(x: &Float) -> Float {
        Float::with_val(REFERENCE_PRECISION, 1.0) / x.clone().tanh()
    }
    pub fn acoth(x: &Float) -> Float {
        (Float::with_val(REFERENCE_PRECISION, 1.0) / x).atanh()
    }
}

pub mod exponential {
    use super::*;
    pub fn exp(x: &Float) -> Float { x.clone().exp() }
    pub fn ln(x: &Float) -> Float { x.clone().ln() }
    pub fn log2(x: &Float) -> Float { x.clone().log2() }
    pub fn log10(x: &Float) -> Float { x.clone().log10() }
    pub fn pow2(x: &Float) -> Float { x.clone().exp2() }
}

pub mod algebraic {
    use super::*;
    pub fn sqrt(x: &Float) -> Float { x.clone().sqrt() }
}
