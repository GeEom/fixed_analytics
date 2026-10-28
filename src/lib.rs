#![doc = include_str!("../docs/crate.md")]
#![no_std]
#![cfg_attr(docsrs, feature(doc_cfg))]

pub mod bounded;
mod consts;
mod cx;
pub mod error;
mod eval;
pub mod kernel;
pub mod ops;
mod real;
mod series;
mod stirling;
pub mod tables;
pub mod traits;

pub use fixed;

/// Every function's worst-case error in ulps, strict, wherever the exact
/// value is in range; out of range a result saturates.
pub const ERROR_ULPS: u32 = 1;

pub use error::{Error, Result};
pub use traits::CordicNumber;

pub use ops::algebraic::{
    hypot, norm, normalize, rsqrt, scale_by, scale_exponent, sqrt, sqrt_product,
};
pub use ops::arith::{div_or_saturate, dot, mul_add, mul_div, sqrt_or_zero};
pub use ops::circular::{acos, asin, asin_from_gap, atan, atan2, cos, sin, sin_cos, tan};
pub use ops::complex::Complex;
pub use ops::elliptic::{
    Modulus, carlson_rc, carlson_rd, carlson_rf, carlson_rf_roots, carlson_rj, carlson_rj_roots,
    ellip_e, ellip_k, ellip_pi, ellip_pi_complete, jacobi_am, jacobi_inverse, jacobi_sn_cn_dn,
};
pub use ops::exponential::{
    exp, exp_moments, expm1, exprel, exprel2, ln, log1p, log1p_ratio, log1p_residual,
    log1p_residual_ratio, log2, log10, pow, pow2,
};
pub use ops::hyperbolic::{acosh, acoth, asinh, atanh, cosh, coth, sinh, sinh_cosh, tanh};
pub use ops::interval::{
    atan_range, enclose, exp_range, sin_cos_range, sin_cos_range_turns, tanh_range,
};
pub use ops::ratio::{atan_ratio, atanh_ratio, stumpff, tan_ratio, tanh_ratio};
pub use ops::solve::{kepler_universal, root_in};
pub use ops::special::{erf, erfc, fresnel, gamma};
pub use ops::turns::{
    atan2_turns, div_by_tau, mul_by_tau, sin_cos_turns, wrap_principal, wrap_radians,
};
