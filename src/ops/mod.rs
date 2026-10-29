//! The mathematical functions, each also re-exported at the crate root.

pub mod algebraic;
pub mod arith;
pub mod circular;
pub mod complex;
pub mod elliptic;
pub mod exponential;
pub mod hyperbolic;
pub mod interval;
pub mod ratio;
pub mod solve;
pub mod special;
pub mod turns;

pub use algebraic::{hypot, norm, normalize, rsqrt, scale_by, scale_exponent, sqrt, sqrt_product};
pub use arith::{
    Rounding, div_or_saturate, div_rounded, dot, dot_rounded, mul_add, mul_add_rounded, mul_div,
    mul_div_rounded, mul_rounded, sqrt_or_zero,
};
pub use circular::{acos, asin, asin_from_gap, atan, atan2, cos, sin, sin_cos, tan};
pub use complex::Complex;
pub use elliptic::{
    Modulus, carlson_rc, carlson_rd, carlson_rf, carlson_rf_roots, carlson_rj, carlson_rj_roots,
    ellip_e, ellip_k, ellip_pi, ellip_pi_complete, ellip_pi_increment, jacobi_am, jacobi_inverse,
    jacobi_sn_cn_dn,
};
pub use exponential::{
    exp, exp_moments, expm1, exprel, exprel2, ln, log1p, log1p_ratio, log1p_residual,
    log1p_residual_ratio, log2, log10, pow, pow2,
};
pub use hyperbolic::{acosh, acoth, asinh, atanh, cosh, coth, sinh, sinh_cosh, tanh};
pub use interval::{
    atan_range, atan_ratio_range, atan2_range, atanh_ratio_range, enclose, exp_moments_range,
    exp_range, exprel_range, exprel2_range, hypot_range, log1p_range, log1p_ratio_range,
    log1p_residual_ratio_range, sin_cos_range, sin_cos_range_turns, sqrt_range, stumpff_range,
    tan_ratio_range, tanh_range, tanh_ratio_range,
};
pub use ratio::{atan_ratio, atanh_ratio, stumpff, tan_ratio, tanh_ratio};
pub use solve::{kepler_universal, root_in};
pub use special::{erf, erfc, fresnel, gamma};
pub use turns::{atan2_turns, div_by_tau, mul_by_tau, sin_cos_turns, wrap_principal, wrap_radians};
