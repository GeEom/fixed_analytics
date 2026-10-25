//! The mathematical functions, each also re-exported at the crate root.

pub mod algebraic;
pub mod circular;
pub mod exponential;
pub mod hyperbolic;

pub use algebraic::sqrt;
pub use circular::{acos, asin, atan, atan2, cos, sin, sin_cos, tan};
pub use exponential::{exp, ln, log2, log10, pow, pow2};
pub use hyperbolic::{acosh, acoth, asinh, atanh, cosh, coth, sinh, sinh_cosh, tanh};
