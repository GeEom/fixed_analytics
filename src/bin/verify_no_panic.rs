//! Binary that instantiates every public function with concrete types.
//!
//! This exists solely to trigger monomorphization so that `no_panic`'s
//! linker-level check can verify that no panic paths survive optimization.
//! It is only compiled under the `verify-no-panic` feature, for three
//! layouts that compile to different code: `I16F16`, `I24F8` (more integer
//! than fractional bits) and `I64F64` (128-bit).

#[cfg(not(feature = "verify-no-panic"))]
compile_error!("this binary should only be built with --features verify-no-panic");

use fixed::types::{I16F16, I24F8, I64F64};
use fixed_analytics::CordicNumber;
use fixed_analytics::bounded::{NonNegative, OpenUnitInterval};
use fixed_analytics::ops::algebraic::sqrt_nonneg;
use fixed_analytics::ops::hyperbolic::atanh_open;
use fixed_analytics::{
    acos, acosh, acoth, asin, asinh, atan, atan2, atanh, cos, cosh, coth, exp, ln, log2, log10,
    pow, pow2, sin, sin_cos, sinh, sinh_cosh, sqrt, tan, tanh,
};

fn exercise<T: CordicNumber>(x: T, y: T, two: T) {
    // Total functions (return T)
    let _ = std::hint::black_box(sin(x));
    let _ = std::hint::black_box(cos(x));
    let _ = std::hint::black_box(tan(x));
    let _ = std::hint::black_box(sin_cos(x));
    let _ = std::hint::black_box(atan(x));
    let _ = std::hint::black_box(atan2(y, x));
    let _ = std::hint::black_box(exp(x));
    let _ = std::hint::black_box(pow2(x));
    let _ = std::hint::black_box(sinh(x));
    let _ = std::hint::black_box(cosh(x));
    let _ = std::hint::black_box(tanh(x));
    let _ = std::hint::black_box(sinh_cosh(x));
    let _ = std::hint::black_box(asinh(x));
    let _ = std::hint::black_box(asinh(two));

    // Fallible functions (return Result<T>)
    let _ = std::hint::black_box(asin(x));
    let _ = std::hint::black_box(acos(x));
    let _ = std::hint::black_box(sqrt(x));
    let _ = std::hint::black_box(ln(x));
    let _ = std::hint::black_box(log2(x));
    let _ = std::hint::black_box(log10(x));
    let _ = std::hint::black_box(acosh(two));
    let _ = std::hint::black_box(atanh(x));
    let _ = std::hint::black_box(coth(x));
    let _ = std::hint::black_box(acoth(two));
    let _ = std::hint::black_box(pow(two, x));

    // Type-safe wrapper functions
    if let Some(nn) = NonNegative::new(x) {
        let _ = std::hint::black_box(sqrt_nonneg(nn));
    }
    if let Some(ou) = OpenUnitInterval::new(x) {
        let _ = std::hint::black_box(atanh_open(ou));
    }
}

fn main() {
    // Use black_box to prevent the optimizer from eliminating calls entirely.
    exercise(
        std::hint::black_box(I16F16::from_num(0.5)),
        std::hint::black_box(I16F16::from_num(0.25)),
        std::hint::black_box(I16F16::from_num(2)),
    );
    exercise(
        std::hint::black_box(I24F8::from_num(0.5)),
        std::hint::black_box(I24F8::from_num(0.25)),
        std::hint::black_box(I24F8::from_num(2)),
    );
    exercise(
        std::hint::black_box(I64F64::from_num(0.5)),
        std::hint::black_box(I64F64::from_num(0.25)),
        std::hint::black_box(I64F64::from_num(2)),
    );
}
