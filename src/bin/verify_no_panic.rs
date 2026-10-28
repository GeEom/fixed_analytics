//! Instantiates every public function so `no_panic`'s link check can prove
//! it panic-free, for three layouts that compile differently: `I16F16`,
//! `I24F8` (more integer than fractional bits) and `I64F64` (128-bit).

#[cfg(not(feature = "verify-no-panic"))]
compile_error!("this binary should only be built with --features verify-no-panic");

use fixed::types::{I16F16, I24F8, I64F64};
use fixed_analytics::CordicNumber;
use fixed_analytics::bounded::{NonNegative, OpenUnitInterval};
use fixed_analytics::ops::algebraic::sqrt_nonneg;
use fixed_analytics::ops::hyperbolic::atanh_open;
use fixed_analytics::{
    acos, acosh, acoth, asin, asin_from_gap, asinh, atan, atan_range, atan2, atanh, cos, cosh,
    coth, div_or_saturate, dot, enclose, exp, exp_range, hypot, ln, log2, log10, mul_add, mul_div,
    norm, normalize, pow, pow2, rsqrt, scale_by, scale_exponent, sin, sin_cos, sin_cos_range, sinh,
    sinh_cosh, sqrt, sqrt_or_zero, sqrt_product, tan, tanh, tanh_range,
};

/// Hands `v` to `black_box`, so the optimizer must compute it in full.
fn keep<V>(v: V) {
    std::hint::black_box(v);
}

fn exercise<T: CordicNumber>(x: T, y: T, two: T) {
    keep(sin(x));
    keep(cos(x));
    keep(tan(x));
    keep(sin_cos(x));
    keep(atan(x));
    keep(atan2(y, x));
    keep(exp(x));
    keep(pow2(x));
    keep(sinh(x));
    keep(cosh(x));
    keep(tanh(x));
    keep(sinh_cosh(x));
    keep(asinh(x));
    keep(asinh(two));

    keep(asin(x));
    keep(acos(x));
    keep(sqrt(x));
    keep(ln(x));
    keep(log2(x));
    keep(log10(x));
    keep(acosh(two));
    keep(atanh(x));
    keep(coth(x));
    keep(acoth(two));
    keep(pow(two, x));
    keep(sqrt_or_zero(x));
    keep(mul_add(x, y, two));
    keep(mul_div(x, y, two));
    keep(div_or_saturate(x, y));
    keep(dot(&[x, y], &[two, x]));
    keep(hypot(x, y));
    keep(norm(&[x, y, two]));
    keep(normalize([x, y, two]));
    keep(scale_exponent(x));
    keep(scale_by(x, 3));
    keep(sqrt_product(x, y));
    keep(rsqrt(x));
    keep(asin_from_gap(x));
    keep(enclose(x));
    keep(sin_cos_range(x, two));
    keep(exp_range(x, two));
    keep(atan_range(x, two));
    keep(tanh_range(x, two));

    if let Some(nn) = NonNegative::new(x) {
        keep(sqrt_nonneg(nn));
    }
    if let Some(ou) = OpenUnitInterval::new(x) {
        keep(atanh_open(ou));
    }
}

fn main() {
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
