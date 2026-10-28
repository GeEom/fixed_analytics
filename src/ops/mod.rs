//! The mathematical functions, each also re-exported at the crate root.

pub mod algebraic;
pub mod arith;
pub mod circular;
pub mod elliptic;
pub mod exponential;
pub mod hyperbolic;
pub mod interval;
pub mod ratio;

pub use algebraic::{hypot, norm, normalize, rsqrt, scale_by, scale_exponent, sqrt, sqrt_product};
pub use arith::{div_or_saturate, dot, mul_add, mul_div, sqrt_or_zero};
pub use circular::{acos, asin, asin_from_gap, atan, atan2, cos, sin, sin_cos, tan};
pub use elliptic::{
    Modulus, carlson_rc, carlson_rd, carlson_rf, carlson_rf_roots, carlson_rj, carlson_rj_roots,
    ellip_e, ellip_k, ellip_pi, ellip_pi_complete, jacobi_am, jacobi_inverse, jacobi_sn_cn_dn,
};
pub use exponential::{
    exp, exp_moments, expm1, exprel, exprel2, ln, log1p, log1p_ratio, log1p_residual,
    log1p_residual_ratio, log2, log10, pow, pow2,
};
pub use hyperbolic::{acosh, acoth, asinh, atanh, cosh, coth, sinh, sinh_cosh, tanh};
pub use interval::{atan_range, enclose, exp_range, sin_cos_range, tanh_range};
pub use ratio::{atan_ratio, atanh_ratio, stumpff, tan_ratio, tanh_ratio};
