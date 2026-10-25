//! Lookup tables for CORDIC and polynomial evaluation, as `i64` holding signed
//! I1F63 values (1 sign bit, 63 fractional bits) that convert to the target type.

pub mod chebyshev;
pub mod circular;
pub mod hyperbolic;

pub use circular::ATAN_TABLE;
pub use hyperbolic::ATANH_TABLE;
