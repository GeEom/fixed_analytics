use crate::{Domain, TestedFunction, reference};
use fixed::types::{I16F16, I32F32};
use rug::Float;

pub fn register() -> Vec<Box<dyn TestedFunction>> {
    vec![
        Box::new(Sinh),
        Box::new(Cosh),
        Box::new(Tanh),
        Box::new(Coth),
        Box::new(Asinh),
        Box::new(Acosh),
        Box::new(Atanh),
        Box::new(Acoth),
    ]
}

struct Sinh;
impl TestedFunction for Sinh {
    fn name(&self) -> &'static str {
        "sinh"
    }
    fn domain(&self) -> Domain {
        Domain::Closed(-8.0, 8.0)
    }
    fn reference(&self, x: &Float) -> Float {
        reference::hyperbolic::sinh(x)
    }
    fn compute_i16f16(&self, x: I16F16) -> I16F16 {
        fixed_analytics::sinh(x)
    }
    fn compute_i32f32(&self, x: I32F32) -> I32F32 {
        fixed_analytics::sinh(x)
    }
}

struct Cosh;
impl TestedFunction for Cosh {
    fn name(&self) -> &'static str {
        "cosh"
    }
    fn domain(&self) -> Domain {
        Domain::Closed(-8.0, 8.0)
    }
    fn reference(&self, x: &Float) -> Float {
        reference::hyperbolic::cosh(x)
    }
    fn compute_i16f16(&self, x: I16F16) -> I16F16 {
        fixed_analytics::cosh(x)
    }
    fn compute_i32f32(&self, x: I32F32) -> I32F32 {
        fixed_analytics::cosh(x)
    }
}

struct Tanh;
impl TestedFunction for Tanh {
    fn name(&self) -> &'static str {
        "tanh"
    }
    fn domain(&self) -> Domain {
        Domain::Closed(-10.0, 10.0)
    }
    fn reference(&self, x: &Float) -> Float {
        reference::hyperbolic::tanh(x)
    }
    fn compute_i16f16(&self, x: I16F16) -> I16F16 {
        fixed_analytics::tanh(x)
    }
    fn compute_i32f32(&self, x: I32F32) -> I32F32 {
        fixed_analytics::tanh(x)
    }
}

struct Coth;
impl TestedFunction for Coth {
    fn name(&self) -> &'static str {
        "coth"
    }
    fn domain(&self) -> Domain {
        Domain::Closed(0.1, 10.0)
    }
    fn reference(&self, x: &Float) -> Float {
        reference::hyperbolic::coth(x)
    }
    fn compute_i16f16(&self, x: I16F16) -> I16F16 {
        fixed_analytics::coth(x).unwrap_or(I16F16::ZERO)
    }
    fn compute_i32f32(&self, x: I32F32) -> I32F32 {
        fixed_analytics::coth(x).unwrap_or(I32F32::ZERO)
    }
}

struct Asinh;
impl TestedFunction for Asinh {
    fn name(&self) -> &'static str {
        "asinh"
    }
    fn domain(&self) -> Domain {
        Domain::Closed(-20.0, 20.0)
    }
    fn reference(&self, x: &Float) -> Float {
        reference::hyperbolic::asinh(x)
    }
    fn compute_i16f16(&self, x: I16F16) -> I16F16 {
        fixed_analytics::asinh(x)
    }
    fn compute_i32f32(&self, x: I32F32) -> I32F32 {
        fixed_analytics::asinh(x)
    }
}

struct Acosh;
impl TestedFunction for Acosh {
    fn name(&self) -> &'static str {
        "acosh"
    }
    fn domain(&self) -> Domain {
        Domain::Closed(1.01, 20.0)
    }
    fn reference(&self, x: &Float) -> Float {
        reference::hyperbolic::acosh(x)
    }
    fn compute_i16f16(&self, x: I16F16) -> I16F16 {
        fixed_analytics::acosh(x).unwrap_or(I16F16::ZERO)
    }
    fn compute_i32f32(&self, x: I32F32) -> I32F32 {
        fixed_analytics::acosh(x).unwrap_or(I32F32::ZERO)
    }
}

struct Atanh;
impl TestedFunction for Atanh {
    fn name(&self) -> &'static str {
        "atanh"
    }
    fn domain(&self) -> Domain {
        Domain::Open(-0.99, 0.99)
    }
    fn reference(&self, x: &Float) -> Float {
        reference::hyperbolic::atanh(x)
    }
    fn compute_i16f16(&self, x: I16F16) -> I16F16 {
        fixed_analytics::atanh(x).unwrap_or(I16F16::ZERO)
    }
    fn compute_i32f32(&self, x: I32F32) -> I32F32 {
        fixed_analytics::atanh(x).unwrap_or(I32F32::ZERO)
    }
}

struct Acoth;
impl TestedFunction for Acoth {
    fn name(&self) -> &'static str {
        "acoth"
    }
    fn domain(&self) -> Domain {
        Domain::OutsideUnit(1.01)
    }
    fn reference(&self, x: &Float) -> Float {
        reference::hyperbolic::acoth(x)
    }
    fn compute_i16f16(&self, x: I16F16) -> I16F16 {
        fixed_analytics::acoth(x).unwrap_or(I16F16::ZERO)
    }
    fn compute_i32f32(&self, x: I32F32) -> I32F32 {
        fixed_analytics::acoth(x).unwrap_or(I32F32::ZERO)
    }
}
