//! Binary that instantiates every public function with a concrete type.
//!
//! This exists solely to trigger monomorphization so that `no_panic`'s
//! linker-level check can verify that no panic paths survive optimization.
//! It is only compiled under the `verify-no-panic` feature.

#[cfg(not(feature = "verify-no-panic"))]
compile_error!("this binary should only be built with --features verify-no-panic");

use fixed::types::I16F16;
use fixed_analytics::bounded::{NonNegative, OpenUnitInterval};
use fixed_analytics::ops::algebraic::sqrt_nonneg;
use fixed_analytics::ops::hyperbolic::atanh_open;
use fixed_analytics::{
    acos, acosh, acoth, asin, asinh, atan, atan2, atanh, cos, cosh, coth, exp, ln, log2, log10,
    pow2, sin, sin_cos, sinh, sinh_cosh, sqrt, tan, tanh,
};

fn main() {
    // Use black_box to prevent the optimizer from eliminating calls entirely.
    let x = std::hint::black_box(I16F16::from_num(0.5));
    let y = std::hint::black_box(I16F16::from_num(0.25));

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

    // Fallible functions (return Result<T>)
    let _ = std::hint::black_box(asin(x));
    let _ = std::hint::black_box(acos(x));
    let _ = std::hint::black_box(sqrt(x));
    let _ = std::hint::black_box(ln(x));
    let _ = std::hint::black_box(log2(x));
    let _ = std::hint::black_box(log10(x));
    let _ = std::hint::black_box(acosh(I16F16::from_num(2)));
    let _ = std::hint::black_box(atanh(x));
    let _ = std::hint::black_box(coth(x));
    let _ = std::hint::black_box(acoth(I16F16::from_num(2)));

    // Type-safe wrapper functions
    let nn = NonNegative::new(x).unwrap();
    let _ = std::hint::black_box(sqrt_nonneg(nn));
    let ou = OpenUnitInterval::new(x).unwrap();
    let _ = std::hint::black_box(atanh_open(ou));
}
