use crate::{reference, Domain, TestedFunction};
use fixed::types::{I16F16, I32F32};
use rug::Float;

pub fn register() -> Vec<Box<dyn TestedFunction>> {
    vec![Box::new(Sqrt)]
}

struct Sqrt;
impl TestedFunction for Sqrt {
    fn name(&self) -> &'static str { "sqrt" }
    fn domain(&self) -> Domain { Domain::Closed(0.0, 10000.0) }
    fn reference(&self, x: &Float) -> Float { reference::algebraic::sqrt(x) }
    fn compute_i16f16(&self, x: I16F16) -> I16F16 { fixed_analytics::sqrt(x).unwrap_or(I16F16::ZERO) }
    fn compute_i32f32(&self, x: I32F32) -> I32F32 { fixed_analytics::sqrt(x).unwrap_or(I32F32::ZERO) }
}
