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
use fixed_analytics::{Complex, erf, erfc, fresnel, gamma, kepler_universal, root_in};
use fixed_analytics::{
    Modulus, carlson_rc, carlson_rd, carlson_rf, carlson_rf_roots, carlson_rj, carlson_rj_roots,
    ellip_e, ellip_k, ellip_pi, ellip_pi_complete, jacobi_am, jacobi_inverse, jacobi_sn_cn_dn,
};
use fixed_analytics::{
    acos, acosh, acoth, asin, asin_from_gap, asinh, atan, atan_range, atan2, atanh, cos, cosh,
    coth, div_or_saturate, dot, enclose, exp, exp_range, hypot, ln, log2, log10, mul_add, mul_div,
    norm, normalize, pow, pow2, rsqrt, scale_by, scale_exponent, sin, sin_cos, sin_cos_range, sinh,
    sinh_cosh, sqrt, sqrt_or_zero, sqrt_product, tan, tanh, tanh_range,
};
use fixed_analytics::{
    atan_ratio, atanh_ratio, exp_moments, expm1, exprel, exprel2, log1p, log1p_ratio,
    log1p_residual, log1p_residual_ratio, stumpff, tan_ratio, tanh_ratio,
};
use fixed_analytics::{
    atan2_turns, div_by_tau, mul_by_tau, sin_cos_range_turns, sin_cos_turns, wrap_principal,
    wrap_radians,
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
    keep(expm1(x));
    keep(exprel(x));
    keep(exprel2(x));
    keep(exp_moments(x));
    keep(log1p(x));
    keep(log1p_ratio(x));
    keep(log1p_residual(x));
    keep(log1p_residual_ratio(x));
    keep(stumpff(x));
    keep(tan_ratio(x));
    keep(atan_ratio(x));
    keep(atanh_ratio(x));
    keep(tanh_ratio(x));
    if let Some(m) = Modulus::from_complement(y) {
        keep(ellip_k(m));
        keep(ellip_e(m));
        keep(ellip_pi_complete(x, m));
        keep(ellip_pi(x, y, x, y, m));
        keep(jacobi_sn_cn_dn(x, m));
        keep(jacobi_am(x, m));
        keep(jacobi_inverse(x, y, x, m));
        keep(m.parameter());
        keep(m.complement());
    }
    keep(Modulus::from_parameter(y));
    keep(carlson_rf(x, y, two));
    keep(carlson_rf_roots(x, y, two));
    keep(carlson_rj(x, y, two, x));
    keep(carlson_rj_roots(x, y, two, x));
    keep(carlson_rc(x, y));
    keep(carlson_rd(x, y, two));
    exercise_more(x, y, two);
}

/// The turns and complex functions, kept apart so neither list grows too long.
fn exercise_more<T: CordicNumber>(x: T, y: T, two: T) {
    keep(sin_cos_turns(x));
    keep(atan2_turns(y, x));
    keep(div_by_tau(x));
    keep(mul_by_tau(x));
    keep(sin_cos_range_turns(x, y));
    keep(wrap_radians(x));
    keep(wrap_principal(x));
    let z = Complex::new(x, y);
    let w = Complex::new(two, x);
    keep(z + w);
    keep(z - w);
    keep(z * w);
    keep(z / w);
    keep(-z);
    keep(z.saturating_add(w));
    keep(z.saturating_sub(w));
    keep(z.saturating_mul(w));
    keep(z.saturating_div(w));
    keep(z.saturating_neg());
    keep(z.ln());
    keep(z.pow(two));
    keep(z.exprel_weighted(y, two));
    keep(z.conj());
    keep(z.norm());
    keep(z.arg());
    keep(z.exp());
    keep(z.exprel());
    keep(z.exprel2());
    keep(z.faddeeva());
    keep(z.erf());
    keep(z.e1());
    keep(z.gamma_upper(two));
    keep(erf(x));
    keep(erfc(x));
    keep(fresnel(x));
    keep(gamma(x));
    keep(root_in(|t: T| t.saturating_sub(y), x, two, 40));
    keep(kepler_universal(two, y, x, two));

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
