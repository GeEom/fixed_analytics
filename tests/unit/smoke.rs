//! Smoke tests and multi-type tests for the library API

#[cfg(test)]
#[allow(clippy::unwrap_used, reason = "test code uses unwrap for conciseness")]
mod tests {
    use fixed::types::I16F16;
    use fixed_analytics::{
        acos, acosh, acoth, asin, asinh, atan, atan2, atanh, cos, cosh, coth, exp, ln, log2, log10,
        sin, sin_cos, sinh, sinh_cosh, sqrt, tan, tanh,
    };

    #[test]
    fn smoke_test_trig() {
        let angle = I16F16::from_num(0.5);
        let _ = sin(angle);
        let _ = cos(angle);
        let _ = tan(angle);
        let _ = sin_cos(angle);
    }

    #[test]
    fn smoke_test_inverse_trig() {
        let x = I16F16::from_num(0.5);
        let _ = asin(x);
        let _ = acos(x);
        let _ = atan(x);
        let _ = atan2(x, I16F16::ONE);
    }

    #[test]
    fn smoke_test_hyperbolic() {
        let x = I16F16::from_num(0.5);
        let _ = sinh(x);
        let _ = cosh(x);
        let _ = tanh(x);
        let _ = coth(x);
        let _ = sinh_cosh(x);
    }

    #[test]
    fn smoke_test_inverse_hyperbolic() {
        let x = I16F16::from_num(0.5);
        let _ = asinh(x);
        let _ = atanh(x);

        let x_large = I16F16::from_num(1.5);
        let _ = acosh(x_large);
        let _ = acoth(x_large);
    }

    #[test]
    fn smoke_test_exponential() {
        let x = I16F16::from_num(0.5);
        let _ = exp(x);
        let _ = ln(x);
        let _ = log2(x);
        let _ = log10(x);
    }

    #[test]
    fn smoke_test_algebraic() {
        let x = I16F16::from_num(2.0);
        let _ = sqrt(x).unwrap();
    }
}

// ==========================================================================
// Multi-type tests: verify functions work with different fixed-point types
// ==========================================================================

#[cfg(test)]
#[allow(clippy::unwrap_used, reason = "test code uses unwrap for conciseness")]
mod multi_type {
    use fixed::types::{I8F24, I32F32};
    use fixed_analytics::{acos, asin, atan, exp, ln, sin_cos, sinh_cosh, sqrt};

    // I32F32 tests - higher precision (32 fractional bits)
    #[test]
    fn trig_i32f32() {
        let angle = I32F32::from_num(0.5);
        let (s, c) = sin_cos(angle);
        let sum_sq: f64 = (s * s + c * c).to_num();
        assert!(
            (sum_sq - 1.0).abs() < 0.001,
            "sin²+cos² = {sum_sq} (I32F32), expected ~1.0"
        );
    }

    #[test]
    fn trig_identity_i32f32() {
        // sin²(x) + cos²(x) = 1 with higher precision
        for i in -10..=10 {
            let angle = I32F32::from_num(i) * I32F32::from_num(0.3);
            let (s, c) = sin_cos(angle);
            let sum_sq: f64 = (s * s + c * c).to_num();
            assert!(
                (sum_sq - 1.0).abs() < 0.01,
                "sin²({}) + cos²({}) = {} (I32F32)",
                angle.to_num::<f64>(),
                angle.to_num::<f64>(),
                sum_sq
            );
        }
    }

    #[test]
    fn inverse_trig_i32f32() {
        let x = I32F32::from_num(0.5);
        let asin_val = asin(x).unwrap();
        let acos_val = acos(x).unwrap();
        let atan_val = atan(x);

        // Verify results are in expected range
        let asin_f64: f64 = asin_val.to_num();
        assert!(asin_f64 > 0.5 && asin_f64 < 0.6, "asin(0.5) ≈ 0.524");

        let acos_f64: f64 = acos_val.to_num();
        assert!(acos_f64 > 1.0 && acos_f64 < 1.1, "acos(0.5) ≈ 1.047");

        let atan_f64: f64 = atan_val.to_num();
        assert!(atan_f64 > 0.45 && atan_f64 < 0.5, "atan(0.5) ≈ 0.464");
    }

    #[test]
    fn hyperbolic_i32f32() {
        let x = I32F32::from_num(0.5);
        let (s, c) = sinh_cosh(x);
        // cosh²(x) - sinh²(x) = 1
        let diff: f64 = (c * c - s * s).to_num();
        assert!(
            (diff - 1.0).abs() < 0.01,
            "cosh²-sinh² = {diff} (I32F32), expected ~1.0"
        );
    }

    #[test]
    fn exp_ln_i32f32() {
        let x = I32F32::from_num(2.0);
        let result = exp(ln(x).unwrap());
        let result_f64: f64 = result.to_num();
        assert!(
            (result_f64 - 2.0).abs() < 0.3,
            "exp(ln(2)) = {result_f64} (I32F32), expected ~2.0"
        );
    }

    #[test]
    fn sqrt_i32f32() {
        let x = I32F32::from_num(4.0);
        let result: f64 = sqrt(x).unwrap().to_num();
        assert!(
            (result - 2.0).abs() < 0.001,
            "sqrt(4) = {result} (I32F32), expected 2.0"
        );
    }

    // I8F24 tests - small integer part, high fractional precision
    #[test]
    fn trig_i8f24() {
        // I8F24 can represent values from about -128 to 127 with 24 fractional bits
        let angle = I8F24::from_num(0.5);
        let (s, c) = sin_cos(angle);
        let sum_sq: f32 = (s * s + c * c).to_num();
        assert!(
            (sum_sq - 1.0).abs() < 0.02,
            "sin²+cos² = {sum_sq} (I8F24), expected ~1.0"
        );
    }

    #[test]
    fn hyperbolic_i8f24() {
        let x = I8F24::from_num(0.5);
        let (s, c) = sinh_cosh(x);
        let diff: f32 = (c * c - s * s).to_num();
        assert!(
            (diff - 1.0).abs() < 0.1,
            "cosh²-sinh² = {diff} (I8F24), expected ~1.0"
        );
    }

    #[test]
    fn sqrt_i8f24() {
        let x = I8F24::from_num(2.0);
        let result: f32 = sqrt(x).unwrap().to_num();
        assert!(
            (result - 1.414).abs() < 0.01,
            "sqrt(2) = {result} (I8F24), expected ~1.414"
        );
    }

    #[test]
    fn ln_i8f24() {
        let x = I8F24::from_num(2.0);
        let result: f32 = ln(x).unwrap().to_num();
        assert!(
            (result - 0.693).abs() < 0.2,
            "ln(2) = {result} (I8F24), expected ~0.693"
        );
    }
}

// I8F8 tests - lower precision, 16-bit total
#[cfg(test)]
#[allow(clippy::unwrap_used, reason = "test code uses unwrap for conciseness")]
mod i8f8 {
    use fixed::types::I8F8;
    use fixed_analytics::{sin_cos, sqrt};

    #[test]
    fn basic_trig() {
        // I8F8 has limited precision but should still work
        let angle = I8F8::from_num(0.5);
        let (s, c) = sin_cos(angle);
        let sum_sq: f32 = (s * s + c * c).to_num();
        // Lower precision type, allow more tolerance
        assert!(
            (sum_sq - 1.0).abs() < 0.2,
            "sin²+cos² = {sum_sq} (I8F8), expected ~1.0"
        );
    }

    #[test]
    fn basic_sqrt() {
        let x = I8F8::from_num(4.0);
        let result: f32 = sqrt(x).unwrap().to_num();
        assert!(
            (result - 2.0).abs() < 0.1,
            "sqrt(4) = {result} (I8F8), expected ~2.0"
        );
    }
}
