use crate::{Domain, TestedFunction, reference};
use fixed::types::{I16F16, I32F32};
use rug::Float;

pub fn register() -> Vec<Box<dyn TestedFunction>> {
    vec![
        Box::new(Sin),
        Box::new(Cos),
        Box::new(Tan),
        Box::new(Asin),
        Box::new(Acos),
        Box::new(Atan),
    ]
}

struct Sin;
impl TestedFunction for Sin {
    fn name(&self) -> &'static str {
        "sin"
    }
    fn domain(&self) -> Domain {
        Domain::Full
    }
    fn reference(&self, x: &Float) -> Float {
        reference::circular::sin(x)
    }
    fn compute_i16f16(&self, x: I16F16) -> I16F16 {
        fixed_analytics::sin(x)
    }
    fn compute_i32f32(&self, x: I32F32) -> I32F32 {
        fixed_analytics::sin(x)
    }
}

struct Cos;
impl TestedFunction for Cos {
    fn name(&self) -> &'static str {
        "cos"
    }
    fn domain(&self) -> Domain {
        Domain::Full
    }
    fn reference(&self, x: &Float) -> Float {
        reference::circular::cos(x)
    }
    fn compute_i16f16(&self, x: I16F16) -> I16F16 {
        fixed_analytics::cos(x)
    }
    fn compute_i32f32(&self, x: I32F32) -> I32F32 {
        fixed_analytics::cos(x)
    }
}

struct Tan;
impl TestedFunction for Tan {
    fn name(&self) -> &'static str {
        "tan"
    }
    fn domain(&self) -> Domain {
        Domain::Open(-1.5, 1.5)
    }
    fn reference(&self, x: &Float) -> Float {
        reference::circular::tan(x)
    }
    fn compute_i16f16(&self, x: I16F16) -> I16F16 {
        fixed_analytics::tan(x)
    }
    fn compute_i32f32(&self, x: I32F32) -> I32F32 {
        fixed_analytics::tan(x)
    }
}

struct Asin;
impl TestedFunction for Asin {
    fn name(&self) -> &'static str {
        "asin"
    }
    fn domain(&self) -> Domain {
        Domain::Closed(-0.99, 0.99)
    }
    fn reference(&self, x: &Float) -> Float {
        reference::circular::asin(x)
    }
    fn compute_i16f16(&self, x: I16F16) -> I16F16 {
        fixed_analytics::asin(x).unwrap_or(I16F16::ZERO)
    }
    fn compute_i32f32(&self, x: I32F32) -> I32F32 {
        fixed_analytics::asin(x).unwrap_or(I32F32::ZERO)
    }
}

struct Acos;
impl TestedFunction for Acos {
    fn name(&self) -> &'static str {
        "acos"
    }
    fn domain(&self) -> Domain {
        Domain::Closed(-0.99, 0.99)
    }
    fn reference(&self, x: &Float) -> Float {
        reference::circular::acos(x)
    }
    fn compute_i16f16(&self, x: I16F16) -> I16F16 {
        fixed_analytics::acos(x).unwrap_or(I16F16::ZERO)
    }
    fn compute_i32f32(&self, x: I32F32) -> I32F32 {
        fixed_analytics::acos(x).unwrap_or(I32F32::ZERO)
    }
}

struct Atan;
impl TestedFunction for Atan {
    fn name(&self) -> &'static str {
        "atan"
    }
    fn domain(&self) -> Domain {
        Domain::Closed(-100.0, 100.0)
    }
    fn reference(&self, x: &Float) -> Float {
        reference::circular::atan(x)
    }
    fn compute_i16f16(&self, x: I16F16) -> I16F16 {
        fixed_analytics::atan(x)
    }
    fn compute_i32f32(&self, x: I32F32) -> I32F32 {
        fixed_analytics::atan(x)
    }
}
