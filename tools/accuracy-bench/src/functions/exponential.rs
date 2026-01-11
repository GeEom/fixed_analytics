use crate::{Domain, TestedFunction, reference};
use fixed::types::{I16F16, I32F32};
use rug::Float;

pub fn register() -> Vec<Box<dyn TestedFunction>> {
    vec![
        Box::new(Exp),
        Box::new(Ln),
        Box::new(Log2),
        Box::new(Log10),
        Box::new(Pow2),
    ]
}

struct Exp;
impl TestedFunction for Exp {
    fn name(&self) -> &'static str {
        "exp"
    }
    fn domain(&self) -> Domain {
        Domain::Closed(-10.0, 8.0)
    }
    fn reference(&self, x: &Float) -> Float {
        reference::exponential::exp(x)
    }
    fn compute_i16f16(&self, x: I16F16) -> I16F16 {
        fixed_analytics::exp(x)
    }
    fn compute_i32f32(&self, x: I32F32) -> I32F32 {
        fixed_analytics::exp(x)
    }
}

struct Ln;
impl TestedFunction for Ln {
    fn name(&self) -> &'static str {
        "ln"
    }
    fn domain(&self) -> Domain {
        Domain::Closed(0.001, 1000.0)
    }
    fn reference(&self, x: &Float) -> Float {
        reference::exponential::ln(x)
    }
    fn compute_i16f16(&self, x: I16F16) -> I16F16 {
        fixed_analytics::ln(x).unwrap_or(I16F16::ZERO)
    }
    fn compute_i32f32(&self, x: I32F32) -> I32F32 {
        fixed_analytics::ln(x).unwrap_or(I32F32::ZERO)
    }
}

struct Log2;
impl TestedFunction for Log2 {
    fn name(&self) -> &'static str {
        "log2"
    }
    fn domain(&self) -> Domain {
        Domain::Closed(0.01, 1000.0)
    }
    fn reference(&self, x: &Float) -> Float {
        reference::exponential::log2(x)
    }
    fn compute_i16f16(&self, x: I16F16) -> I16F16 {
        fixed_analytics::log2(x).unwrap_or(I16F16::ZERO)
    }
    fn compute_i32f32(&self, x: I32F32) -> I32F32 {
        fixed_analytics::log2(x).unwrap_or(I32F32::ZERO)
    }
}

struct Log10;
impl TestedFunction for Log10 {
    fn name(&self) -> &'static str {
        "log10"
    }
    fn domain(&self) -> Domain {
        Domain::Closed(0.01, 1000.0)
    }
    fn reference(&self, x: &Float) -> Float {
        reference::exponential::log10(x)
    }
    fn compute_i16f16(&self, x: I16F16) -> I16F16 {
        fixed_analytics::log10(x).unwrap_or(I16F16::ZERO)
    }
    fn compute_i32f32(&self, x: I32F32) -> I32F32 {
        fixed_analytics::log10(x).unwrap_or(I32F32::ZERO)
    }
}

struct Pow2;
impl TestedFunction for Pow2 {
    fn name(&self) -> &'static str {
        "pow2"
    }
    fn domain(&self) -> Domain {
        Domain::Closed(-10.0, 10.0)
    }
    fn reference(&self, x: &Float) -> Float {
        reference::exponential::pow2(x)
    }
    fn compute_i16f16(&self, x: I16F16) -> I16F16 {
        fixed_analytics::pow2(x)
    }
    fn compute_i32f32(&self, x: I32F32) -> I32F32 {
        fixed_analytics::pow2(x)
    }
}
