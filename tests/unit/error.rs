//! Tests for error types

#[cfg(test)]
mod tests {
    use fixed_analytics::Error;

    #[test]
    fn error_display_domain_error() {
        let err = Error::DomainError {
            function: "test_fn",
            expected: "positive value",
        };
        let msg = format!("{err}");
        assert!(msg.contains("test_fn"));
        assert!(msg.contains("positive value"));
    }

    #[test]
    fn error_constructors() {
        // Test domain error constructor
        let err = Error::domain("asin", "value in range [-1, 1]");
        let msg = format!("{err}");
        assert!(msg.contains("asin"));
        assert!(msg.contains("[-1, 1]"));

        // Test overflow error constructor
        let err = Error::overflow("exp");
        let msg = format!("{err}");
        assert!(msg.contains("exp"));
        assert!(msg.contains("overflow"));
    }

    #[test]
    fn error_equality() {
        let err1 = Error::domain("test", "expected");
        let err2 = Error::domain("test", "expected");
        let err3 = Error::domain("other", "expected");

        assert_eq!(err1, err2);
        assert_ne!(err1, err3);
    }
}
