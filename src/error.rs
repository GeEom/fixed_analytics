//! Error types for CORDIC operations.

use core::fmt;

/// Errors that can occur during CORDIC computations.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Error {
    /// The input is outside the function's domain, as for `asin(2.0)`.
    DomainError {
        /// Name of the function that encountered the error.
        function: &'static str,
        /// Human-readable description of the valid domain.
        expected: &'static str,
    },
}

impl Error {
    /// A domain error for the given function.
    #[must_use]
    pub const fn domain(function: &'static str, expected: &'static str) -> Self {
        Self::DomainError { function, expected }
    }
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::DomainError { function, expected } => {
                write!(
                    f,
                    "{function}: input outside valid domain, expected {expected}"
                )
            }
        }
    }
}

#[cfg(feature = "std")]
extern crate std;

#[cfg(feature = "std")]
impl std::error::Error for Error {}

/// A specialized `Result` type for CORDIC operations.
pub type Result<T> = core::result::Result<T, Error>;
