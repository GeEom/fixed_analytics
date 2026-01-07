//! Error types for CORDIC operations.

use core::fmt;

/// Errors that can occur during CORDIC computations.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Error {
    /// Input value is outside the valid domain for the function.
    ///
    /// For example, `asin(2.0)` would produce this error since
    /// arcsine is only defined for inputs in [-1, 1].
    DomainError {
        /// Name of the function that encountered the error.
        function: &'static str,
        /// Human-readable description of the valid domain.
        expected: &'static str,
    },

    /// The computation would overflow the fixed-point representation.
    Overflow {
        /// Name of the function that encountered the error.
        function: &'static str,
    },
}

impl Error {
    /// Create a domain error for the given function.
    #[must_use]
    pub const fn domain(function: &'static str, expected: &'static str) -> Self {
        Self::DomainError { function, expected }
    }

    /// Create an overflow error for the given function.
    #[must_use]
    pub const fn overflow(function: &'static str) -> Self {
        Self::Overflow { function }
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
            Self::Overflow { function } => {
                write!(f, "{function}: result would overflow")
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
