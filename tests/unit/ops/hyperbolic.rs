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

    mod saturation {
        use super::*;
        use fixed::types::I32F32;

        /// Check if I16F16 value is saturated to MAX (within 0.01%)
        fn is_max_16(val: I16F16) -> bool {
            val.to_num::<f32>() >= I16F16::MAX.to_num::<f32>() * 0.9999
        }

        /// Check if I16F16 value is saturated to MIN (within 0.01%)
        fn is_min_16(val: I16F16) -> bool {
            val.to_num::<f32>() <= I16F16::MIN.to_num::<f32>() * 0.9999
        }

        /// Check if I32F32 value is saturated to MAX (within 0.01%)
        fn is_max_32(val: I32F32) -> bool {
            val.to_num::<f64>() >= I32F32::MAX.to_num::<f64>() * 0.9999
        }

        /// Check if I32F32 value is saturated to MIN (within 0.01%)
        fn is_min_32(val: I32F32) -> bool {
            val.to_num::<f64>() <= I32F32::MIN.to_num::<f64>() * 0.9999
        }

        // ===== sinh saturation thresholds =====
        // I16F16: saturates to MAX at x >= 11.1, to MIN at x <= -11.1
        // I32F32: saturates to MAX at x >= 22.2, to MIN at x <= -22.2

        #[test]
        fn sinh_i16f16_positive_threshold() {
            // Below threshold: should NOT saturate
            assert!(
                !is_max_16(sinh(I16F16::from_num(11.0))),
                "sinh(11.0) should not saturate"
            );
            // At threshold: should saturate
            assert!(
                is_max_16(sinh(I16F16::from_num(11.1))),
                "sinh(11.1) should saturate to MAX"
            );
        }

        #[test]
        fn sinh_i16f16_negative_threshold() {
            // Above threshold: should NOT saturate
            assert!(
                !is_min_16(sinh(I16F16::from_num(-11.0))),
                "sinh(-11.0) should not saturate to MIN"
            );
            // At threshold: should saturate
            assert!(
                is_min_16(sinh(I16F16::from_num(-11.1))),
                "sinh(-11.1) should saturate to MIN"
            );
        }

        #[test]
        fn sinh_i32f32_positive_threshold() {
            // Below threshold: should NOT saturate
            assert!(
                !is_max_32(sinh(I32F32::from_num(22.1))),
                "sinh(22.1) should not saturate"
            );
            // At threshold: should saturate
            assert!(
                is_max_32(sinh(I32F32::from_num(22.2))),
                "sinh(22.2) should saturate to MAX"
            );
        }

        #[test]
        fn sinh_i32f32_negative_threshold() {
            // Above threshold: should NOT saturate
            assert!(
                !is_min_32(sinh(I32F32::from_num(-22.1))),
                "sinh(-22.1) should not saturate to MIN"
            );
            // At threshold: should saturate
            assert!(
                is_min_32(sinh(I32F32::from_num(-22.2))),
                "sinh(-22.2) should saturate to MIN"
            );
        }

        // ===== cosh saturation thresholds =====
        // I16F16: saturates to MAX at |x| >= 11.1
        // I32F32: saturates to MAX at |x| >= 22.2
        // (cosh is even, so positive and negative thresholds are symmetric)

        #[test]
        fn cosh_i16f16_positive_threshold() {
            // Below threshold: should NOT saturate
            assert!(
                !is_max_16(cosh(I16F16::from_num(11.0))),
                "cosh(11.0) should not saturate"
            );
            // At threshold: should saturate
            assert!(
                is_max_16(cosh(I16F16::from_num(11.1))),
                "cosh(11.1) should saturate to MAX"
            );
        }

        #[test]
        fn cosh_i16f16_negative_threshold() {
            // Above threshold (less negative): should NOT saturate
            assert!(
                !is_max_16(cosh(I16F16::from_num(-11.0))),
                "cosh(-11.0) should not saturate"
            );
            // At threshold: should saturate
            assert!(
                is_max_16(cosh(I16F16::from_num(-11.1))),
                "cosh(-11.1) should saturate to MAX"
            );
        }

        #[test]
        fn cosh_i32f32_positive_threshold() {
            // Below threshold: should NOT saturate
            assert!(
                !is_max_32(cosh(I32F32::from_num(22.1))),
                "cosh(22.1) should not saturate"
            );
            // At threshold: should saturate
            assert!(
                is_max_32(cosh(I32F32::from_num(22.2))),
                "cosh(22.2) should saturate to MAX"
            );
        }

        #[test]
        fn cosh_i32f32_negative_threshold() {
            // Above threshold (less negative): should NOT saturate
            assert!(
                !is_max_32(cosh(I32F32::from_num(-22.1))),
                "cosh(-22.1) should not saturate"
            );
            // At threshold: should saturate
            assert!(
                is_max_32(cosh(I32F32::from_num(-22.2))),
                "cosh(-22.2) should saturate to MAX"
            );
        }

        // ===== sinh/cosh threshold consistency =====

        #[test]
        fn sinh_cosh_thresholds_match() {
            // sinh and cosh should have the same saturation threshold
            // (both grow as e^x/2 for large |x|)
            assert!(
                is_max_16(sinh(I16F16::from_num(11.1))),
                "sinh(11.1) should saturate"
            );
            assert!(
                is_max_16(cosh(I16F16::from_num(11.1))),
                "cosh(11.1) should saturate"
            );
            assert!(
                is_max_32(sinh(I32F32::from_num(22.2))),
                "sinh(22.2) should saturate"
            );
            assert!(
                is_max_32(cosh(I32F32::from_num(22.2))),
                "cosh(22.2) should saturate"
            );
        }
    }
}

/// Integer-bit-poor types, and the logarithmic inverse forms.
#[cfg(test)]
#[allow(clippy::unwrap_used, reason = "test code uses unwrap for conciseness")]
mod wide_and_narrow {
    use crate::unit::support::Lcg;
    use fixed::types::{I4F12, I4F60, I16F16, I32F32, I64F64};
    use fixed_analytics::{acosh, asinh, sinh_cosh};

    #[test]
    fn sinh_cosh_works_on_types_with_few_integer_bits() {
        // I4F12 and I4F60 (range ±8) cannot hold the Taylor divisors in T.
        for i in 0..=100 {
            let v = -1.5 + 3.0 * f64::from(i) / 100.0;
            let x12 = I4F12::from_num(v);
            let (s12, c12) = sinh_cosh(x12);
            let (s12, c12, v12): (f64, f64, f64) = (s12.to_num(), c12.to_num(), x12.to_num());
            assert!((s12 - v12.sinh()).abs() < 5e-3, "I4F12 sinh({v12}) = {s12}");
            assert!((c12 - v12.cosh()).abs() < 5e-3, "I4F12 cosh({v12}) = {c12}");
            let x60 = I4F60::from_num(v);
            let (s60, c60) = sinh_cosh(x60);
            let (s60, c60, v60): (f64, f64, f64) = (s60.to_num(), c60.to_num(), x60.to_num());
            assert!(
                (s60 - v60.sinh()).abs() < 4e-15,
                "I4F60 sinh({v60}) = {s60}"
            );
            assert!(
                (c60 - v60.cosh()).abs() < 4e-15,
                "I4F60 cosh({v60}) = {c60}"
            );
        }
    }

    #[test]
    fn asinh_acosh_track_f64_at_i64f64() {
        let mut rng = Lcg(0xA5);
        for _ in 0..2000 {
            let sign = if rng.unit() < 0.5 { -1.0 } else { 1.0 };
            let x = I64F64::from_num(sign * (rng.range(-20.0, 60.0)).exp2());
            let got: f64 = asinh(x).to_num();
            let want = x.to_num::<f64>().asinh();
            assert!(
                (got - want).abs() < 1e-15 * want.abs().max(1.0),
                "asinh({x}) = {got}, want {want}"
            );
            // acosh is ill-conditioned near 1; start where f64 is trustworthy.
            let w = I64F64::from_num(1.0 + (rng.range(-10.0, 61.0)).exp2());
            let got_w: f64 = acosh(w).unwrap().to_num();
            let want_w = w.to_num::<f64>().acosh();
            assert!(
                (got_w - want_w).abs() < 1e-13 * want_w.abs().max(1.0),
                "acosh({w}) = {got_w}, want {want_w}"
            );
        }
        // Both branches around the 1.5 threshold of acosh.
        for w in [1.01f64, 1.25, 1.49, 1.5, 1.51, 2.0] {
            let got: f64 = acosh(I64F64::from_num(w)).unwrap().to_num();
            assert!((got - w.acosh()).abs() < 1e-9, "acosh({w}) = {got}");
        }
    }

    #[test]
    fn asinh_acosh_near_the_type_bounds() {
        // The product would overflow, so the logarithm is split.
        for (x, want) in [
            (I16F16::MAX, (2.0 * I16F16::MAX.to_num::<f64>()).ln()),
            (I16F16::from_num(10_000), 10_000f64.asinh()),
            (I16F16::MIN, -(2.0 * I16F16::MAX.to_num::<f64>()).ln()),
        ] {
            let got: f64 = asinh(x).to_num();
            assert!((got - want).abs() < 2e-4, "asinh({x}) = {got}, want {want}");
        }
        for x in [I16F16::MAX, I16F16::from_num(20_000)] {
            let want = (2.0f64 * x.to_num::<f64>()).ln();
            let got: f64 = acosh(x).unwrap().to_num();
            assert!((got - want).abs() < 2e-4, "acosh({x}) = {got}, want {want}");
        }
        let want_32 = (2.0f64 * I32F32::MAX.to_num::<f64>()).ln();
        let got_32: f64 = acosh(I32F32::MAX).unwrap().to_num();
        assert!(
            (got_32 - want_32).abs() < 1e-8,
            "acosh(I32F32::MAX) = {got_32}"
        );
        let want_64 = (2.0f64 * I64F64::MAX.to_num::<f64>()).ln();
        let got_64: f64 = acosh(I64F64::MAX).unwrap().to_num();
        assert!(
            (got_64 - want_64).abs() < 1e-12,
            "acosh(I64F64::MAX) = {got_64}"
        );
        let got_asinh: f64 = asinh(I64F64::MAX).to_num();
        assert!(
            (got_asinh - want_64).abs() < 1e-12,
            "asinh(MAX) = {got_asinh}, want {want_64}"
        );
        assert_eq!(asinh(-I64F64::from_num(3)), -asinh(I64F64::from_num(3)));
    }

    #[test]
    fn asinh_is_accurate_for_moderate_arguments_at_i16f16() {
        for i in 1..=200 {
            let v = f64::from(i).mul_add(0.1, 1.0);
            let got: f64 = asinh(I16F16::from_num(v)).to_num();
            let want = v.asinh();
            assert!(
                (got - want).abs() < 16.0 / 65_536.0,
                "asinh({v}) = {got}, want {want}"
            );
        }
    }
}
