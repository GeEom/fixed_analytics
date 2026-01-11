//! Tests for hyperbolic functions

#[cfg(test)]
#[allow(clippy::unwrap_used, reason = "test code uses unwrap for conciseness")]
mod tests {
    use fixed::types::{I16F16, I32F32};
    use fixed_analytics::{acosh, acoth, asinh, atanh, cosh, coth, sinh, sinh_cosh, tanh};

    const TOLERANCE: f32 = 0.05;

    fn approx_eq(a: I16F16, b: f32) -> bool {
        (a.to_num::<f32>() - b).abs() < TOLERANCE
    }

    #[test]
    fn sinh_special_values() {
        assert!(approx_eq(sinh(I16F16::ZERO), 0.0));
    }

    #[test]
    fn cosh_special_values() {
        assert!(approx_eq(cosh(I16F16::ZERO), 1.0));
    }

    #[test]
    fn tanh_special_values() {
        assert!(approx_eq(tanh(I16F16::ZERO), 0.0));
    }

    #[test]
    fn hyperbolic_identity() {
        // cosh²(x) - sinh²(x) = 1
        for i in -5..=5 {
            let x = I16F16::from_num(i) * I16F16::from_num(0.2);
            let (sh, ch) = sinh_cosh(x);
            let diff: f32 = (ch * ch - sh * sh).to_num();
            assert!(
                (diff - 1.0).abs() < 0.1,
                "cosh²({}) - sinh²({}) = {}, expected ~1.0",
                x.to_num::<f32>(),
                x.to_num::<f32>(),
                diff
            );
        }
    }

    #[test]
    fn atanh_domain_check() {
        assert!(atanh(I16F16::from_num(1.5)).is_err());
        assert!(atanh(I16F16::from_num(-1.5)).is_err());
        assert!(atanh(I16F16::ONE).is_err());
        assert!(atanh(I16F16::from_num(0.5)).is_ok());
    }

    #[test]
    fn acosh_domain_check() {
        assert!(acosh(I16F16::from_num(0.5)).is_err());
        assert!(acosh(I16F16::ONE).is_ok());
        assert!(acosh(I16F16::from_num(2.0)).is_ok());
    }

    #[test]
    fn acoth_domain_check() {
        // acoth requires |x| > 1
        assert!(acoth(I16F16::from_num(0.5)).is_err());
        assert!(acoth(I16F16::ONE).is_err());
        assert!(acoth(-I16F16::ONE).is_err());
        assert!(acoth(I16F16::from_num(1.5)).is_ok());
        assert!(acoth(I16F16::from_num(-1.5)).is_ok());
    }

    #[test]
    fn acoth_values() {
        // acoth(x) = atanh(1/x)
        // acoth(2) = atanh(0.5) ≈ 0.5493
        let result = acoth(I16F16::from_num(2.0));
        assert!(result.is_ok());
        let val: f32 = result.unwrap().to_num();
        assert!(
            (val - 0.5493).abs() < TOLERANCE,
            "acoth(2) expected ~0.5493, got {val}"
        );

        // acoth(-2) = atanh(-0.5) ≈ -0.5493
        let result_neg = acoth(I16F16::from_num(-2.0));
        assert!(result_neg.is_ok());
        let val_neg: f32 = result_neg.unwrap().to_num();
        assert!(
            (val_neg + 0.5493).abs() < TOLERANCE,
            "acoth(-2) expected ~-0.5493, got {val_neg}"
        );
    }

    #[test]
    fn sinh_asinh_roundtrip() {
        // sinh(asinh(x)) ≈ x for various x
        for i in -10..=10 {
            let x = I16F16::from_num(i) * I16F16::from_num(0.3);
            let result = sinh(asinh(x));
            let x_f32: f32 = x.to_num();
            let result_f32: f32 = result.to_num();
            assert!(
                (result_f32 - x_f32).abs() < 0.15,
                "sinh(asinh({x_f32})) = {result_f32}, expected {x_f32}"
            );
        }
    }

    #[test]
    fn cosh_acosh_roundtrip() {
        // cosh(acosh(x)) ≈ x for x >= 1
        for i in 1..=10 {
            let x = I16F16::from_num(i) * I16F16::from_num(0.5);
            if x >= I16F16::ONE {
                let result = cosh(acosh(x).unwrap());
                let x_f32: f32 = x.to_num();
                let result_f32: f32 = result.to_num();
                assert!(
                    (result_f32 - x_f32).abs() < 0.2,
                    "cosh(acosh({x_f32})) = {result_f32}, expected {x_f32}"
                );
            }
        }
    }

    #[test]
    fn tanh_atanh_roundtrip() {
        // tanh(atanh(x)) ≈ x for x in (-1, 1)
        for i in -9..=9 {
            let x = I16F16::from_num(i) * I16F16::from_num(0.1);
            let result = tanh(atanh(x).unwrap());
            let x_f32: f32 = x.to_num();
            let result_f32: f32 = result.to_num();
            assert!(
                (result_f32 - x_f32).abs() < 0.1,
                "tanh(atanh({x_f32})) = {result_f32}, expected {x_f32}"
            );
        }
    }

    #[test]
    fn atanh_near_boundary() {
        // atanh approaches infinity as |x| approaches 1
        // Test values close to but not at the boundary
        let near_one = I16F16::from_num(0.99);
        let result = atanh(near_one);
        assert!(result.is_ok());
        let val: f32 = result.unwrap().to_num();
        // atanh(0.99) ≈ 2.647
        assert!(val > 2.0, "atanh(0.99) = {val}, expected > 2.0");

        let neg_near_one = I16F16::from_num(-0.99);
        let result_neg = atanh(neg_near_one);
        assert!(result_neg.is_ok());
        let val_neg: f32 = result_neg.unwrap().to_num();
        assert!(val_neg < -2.0, "atanh(-0.99) = {val_neg}, expected < -2.0");
    }

    #[test]
    fn acosh_at_boundary() {
        // acosh(1) should be exactly 0
        let result: f32 = acosh(I16F16::ONE).unwrap().to_num();
        assert!(result.abs() < 0.01, "acosh(1) = {result}, expected 0");

        // acosh near boundary (1.01)
        let near_one = I16F16::from_num(1.01);
        let result_near: f32 = acosh(near_one).unwrap().to_num();
        // acosh(1.01) ≈ 0.141
        assert!(
            result_near > 0.0 && result_near < 0.3,
            "acosh(1.01) = {result_near}, expected ~0.14"
        );
    }

    #[test]
    fn tanh_large_values() {
        // tanh should approach ±1 for large values
        let large = I16F16::from_num(10.0);
        let result: f32 = tanh(large).to_num();
        assert!(
            (result - 1.0).abs() < 0.01,
            "tanh(10) = {result}, expected ~1.0"
        );

        let neg_large = I16F16::from_num(-10.0);
        let result_neg: f32 = tanh(neg_large).to_num();
        assert!(
            (result_neg + 1.0).abs() < 0.01,
            "tanh(-10) = {result_neg}, expected ~-1.0"
        );
    }

    #[test]
    fn sinh_cosh_large_values() {
        // Test argument reduction for large values
        let large = I16F16::from_num(5.0);
        let (s, c) = sinh_cosh(large);

        // Verify cosh²-sinh² = 1 identity
        let diff: f32 = (c * c - s * s).to_num();
        assert!(
            (diff - 1.0).abs() < 0.2,
            "cosh²(5) - sinh²(5) = {diff}, expected ~1.0"
        );
    }

    #[test]
    fn coth_at_zero() {
        // coth(0) is undefined (pole), should return DomainError
        let result = coth(I16F16::ZERO);
        assert!(result.is_err(), "coth(0) should return Err");
    }

    #[test]
    fn coth_nonzero_values() {
        // coth(x) = cosh(x)/sinh(x)
        // coth(1) ≈ 1.3130
        let result: f32 = coth(I16F16::ONE).unwrap().to_num();
        assert!(
            (result - 1.3130).abs() < TOLERANCE,
            "coth(1) = {result}, expected ~1.3130"
        );

        // coth(-1) ≈ -1.3130
        let result_neg: f32 = coth(-I16F16::ONE).unwrap().to_num();
        assert!(
            (result_neg + 1.3130).abs() < TOLERANCE,
            "coth(-1) = {result_neg}, expected ~-1.3130"
        );
    }

    #[test]
    fn sinh_cosh_small_values_high_precision() {
        // Test Taylor series approximation for high-precision types (≥24 frac bits)
        // Uses fifth/sixth-order Taylor series
        let small = I32F32::from_num(0.03); // Below 0.05 threshold
        let (s, c) = sinh_cosh(small);

        // sinh(0.03) ≈ 0.03 (Taylor: sinh(x) ≈ x)
        let s_f32: f32 = s.to_num();
        assert!(
            (s_f32 - 0.03).abs() < 0.01,
            "sinh(0.03) = {s_f32}, expected ~0.03"
        );

        // cosh(0.03) ≈ 1.00045 (Taylor: 1 + x²/2)
        let c_f32: f32 = c.to_num();
        assert!(
            (c_f32 - 1.00045).abs() < 0.01,
            "cosh(0.03) = {c_f32}, expected ~1.00045"
        );

        // Test negative small value for high precision
        let small_neg = I32F32::from_num(-0.03);
        let (s_neg, c_neg) = sinh_cosh(small_neg);
        let s_neg_f32: f32 = s_neg.to_num();
        let c_neg_f32: f32 = c_neg.to_num();
        assert!(
            (s_neg_f32 + 0.03).abs() < 0.01,
            "sinh(-0.03) = {s_neg_f32}, expected ~-0.03"
        );
        assert!(
            (c_neg_f32 - 1.00045).abs() < 0.01,
            "cosh(-0.03) = {c_neg_f32}, expected ~1.00045"
        );

        // Additional test with even smaller value to ensure full Taylor path coverage
        let tiny = core::hint::black_box(I32F32::from_num(0.01));
        let (s_tiny, c_tiny) = sinh_cosh(tiny);
        // Use black_box to prevent optimization
        let s_tiny_f32: f32 = core::hint::black_box(s_tiny).to_num();
        let c_tiny_f32: f32 = core::hint::black_box(c_tiny).to_num();
        assert!(
            (s_tiny_f32 - 0.01).abs() < 0.001,
            "sinh(0.01) = {s_tiny_f32}, expected ~0.01"
        );
        assert!(
            (c_tiny_f32 - 1.0).abs() < 0.001,
            "cosh(0.01) = {c_tiny_f32}, expected ~1.0"
        );
    }

    #[test]
    fn sinh_cosh_small_values() {
        // Test Taylor series approximation for very small values (|x| < 0.1)
        // sinh(x) ≈ x for small x
        // cosh(x) ≈ 1 + x²/2 for small x
        let small = I16F16::from_num(0.05);
        let (s, c) = sinh_cosh(small);

        // sinh(0.05) ≈ 0.05 (Taylor: sinh(x) ≈ x)
        let s_f32: f32 = s.to_num();
        assert!(
            (s_f32 - 0.05).abs() < 0.01,
            "sinh(0.05) = {s_f32}, expected ~0.05"
        );

        // cosh(0.05) ≈ 1.00125 (Taylor: 1 + x²/2 = 1 + 0.00125)
        let c_f32: f32 = c.to_num();
        assert!(
            (c_f32 - 1.00125).abs() < 0.01,
            "cosh(0.05) = {c_f32}, expected ~1.00125"
        );

        // Also test negative small value
        let small_neg = I16F16::from_num(-0.05);
        let (s_neg, c_neg) = sinh_cosh(small_neg);
        let s_neg_f32: f32 = s_neg.to_num();
        let c_neg_f32: f32 = c_neg.to_num();
        assert!(
            (s_neg_f32 + 0.05).abs() < 0.01,
            "sinh(-0.05) = {s_neg_f32}, expected ~-0.05"
        );
        assert!(
            (c_neg_f32 - 1.00125).abs() < 0.01,
            "cosh(-0.05) = {c_neg_f32}, expected ~1.00125"
        );
    }
}
