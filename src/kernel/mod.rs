#![doc = include_str!("../../docs/kernel.md")]

mod cordic;

pub use crate::kernel::cordic::{circular_vectoring, hyperbolic_vectoring};
