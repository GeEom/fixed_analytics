//! Accuracy benchmarking framework for fixed_analytics.

pub mod functions;
pub mod metrics;
pub mod reference;
pub mod report;
pub mod sampling;

use fixed::traits::Fixed;
use metrics::ErrorStats;
use rug::Float;
use sampling::SampleStrategy;

pub const REFERENCE_PRECISION: u32 = 256;

#[derive(Debug, Clone)]
pub enum Domain {
    Full,
    Open(f64, f64),
    Closed(f64, f64),
    Positive,
    OutsideUnit(f64),
}

impl Domain {
    pub fn contains(&self, x: f64) -> bool {
        match self {
            Domain::Full => true,
            Domain::Open(a, b) => x > *a && x < *b,
            Domain::Closed(a, b) => x >= *a && x <= *b,
            Domain::Positive => x > 0.0,
            Domain::OutsideUnit(bound) => x.abs() > *bound,
        }
    }

    pub fn sampling_bounds(&self) -> (f64, f64) {
        match self {
            Domain::Full => (-100.0, 100.0),
            Domain::Open(a, b) | Domain::Closed(a, b) => (*a, *b),
            Domain::Positive => (1e-6, 1000.0),
            Domain::OutsideUnit(bound) => (*bound + 0.01, 100.0),
        }
    }
}

pub trait TestedFunction: Send + Sync {
    fn name(&self) -> &'static str;
    fn domain(&self) -> Domain;
    fn reference(&self, x: &Float) -> Float;
    fn compute_i16f16(&self, x: fixed::types::I16F16) -> fixed::types::I16F16;
    fn compute_i32f32(&self, x: fixed::types::I32F32) -> fixed::types::I32F32;
}

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct FunctionResult {
    pub name: String,
    pub i16f16: ErrorStats,
    pub i32f32: ErrorStats,
    pub samples_tested: usize,
}

pub fn test_function(func: &dyn TestedFunction, strategy: &SampleStrategy) -> FunctionResult {
    let domain = func.domain();
    let (lo, hi) = domain.sampling_bounds();
    let points = strategy.generate(lo, hi);

    let mut i16f16_errors = Vec::new();
    let mut i32f32_errors = Vec::new();
    let mut tested = 0;

    for &x_f64 in &points {
        if !domain.contains(x_f64) {
            continue;
        }

        let x_mpfr = Float::with_val(REFERENCE_PRECISION, x_f64);
        let ref_f64 = func.reference(&x_mpfr).to_f64();

        if let Some(x) = try_from_f64::<fixed::types::I16F16>(x_f64) {
            let result: f64 = func.compute_i16f16(x).to_num();
            if let Some(err) = metrics::compute_error(result, ref_f64) {
                i16f16_errors.push(err);
            }
        }

        if let Some(x) = try_from_f64::<fixed::types::I32F32>(x_f64) {
            let result: f64 = func.compute_i32f32(x).to_num();
            if let Some(err) = metrics::compute_error(result, ref_f64) {
                i32f32_errors.push(err);
            }
        }

        tested += 1;
    }

    FunctionResult {
        name: func.name().to_string(),
        i16f16: ErrorStats::from_errors(&i16f16_errors),
        i32f32: ErrorStats::from_errors(&i32f32_errors),
        samples_tested: tested,
    }
}

fn try_from_f64<T: Fixed>(x: f64) -> Option<T> {
    let max: f64 = T::MAX.to_num();
    let min: f64 = T::MIN.to_num();
    if x > max || x < min || !x.is_finite() {
        return None;
    }
    Some(T::from_num(x))
}

pub type FunctionRegistry = Vec<Box<dyn TestedFunction>>;

pub fn build_registry() -> FunctionRegistry {
    let mut reg: FunctionRegistry = Vec::new();
    reg.extend(functions::circular::register());
    reg.extend(functions::hyperbolic::register());
    reg.extend(functions::exponential::register());
    reg.extend(functions::algebraic::register());
    reg
}
