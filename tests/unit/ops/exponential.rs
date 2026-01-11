//! Tests for exponential and logarithmic functions

#[cfg(test)]
#[allow(clippy::unwrap_used, reason = "test code uses unwrap for conciseness")]
mod tests {
    use fixed::types::I16F16;
    use fixed_analytics::{exp, ln, log2, log10, pow2};

    const TOLERANCE: f32 = 0.15;

    fn approx_eq(a: I16F16, b: f32, tolerance: f32) -> bool {
        (a.to_num::<f32>() - b).abs() < tolerance
    }

    #[test]
    fn exp_special_values() {
        assert!(approx_eq(exp(I16F16::ZERO), 1.0, TOLERANCE));
        assert!(approx_eq(exp(I16F16::ONE), core::f32::consts::E, TOLERANCE));
    }

    #[test]
    fn exp_negative() {
        let result = exp(-I16F16::ONE);
        let expected = 1.0 / core::f32::consts::E;
        assert!(approx_eq(result, expected, TOLERANCE));
    }

    #[test]
    fn ln_special_values() {
        assert!(approx_eq(ln(I16F16::ONE).unwrap(), 0.0, TOLERANCE));
        let ln_e: f32 = ln(I16F16::E).unwrap().to_num();
        assert!((ln_e - 1.0).abs() < 0.25, "ln(e) = {ln_e}, expected ~1.0");
    }

    #[test]
    fn ln_domain_check() {
        assert!(ln(I16F16::ZERO).is_err());
        assert!(ln(I16F16::from_num(-1.0)).is_err());
        assert!(ln(I16F16::from_num(0.5)).is_ok());
    }

    #[test]
    fn log2_powers_of_two() {
        assert!(approx_eq(
            log2(I16F16::from_num(1.0)).unwrap(),
            0.0,
            TOLERANCE
        ));
        assert!(approx_eq(log2(I16F16::from_num(2.0)).unwrap(), 1.0, 0.25));
        assert!(approx_eq(log2(I16F16::from_num(4.0)).unwrap(), 2.0, 0.3));
        assert!(approx_eq(log2(I16F16::from_num(8.0)).unwrap(), 3.0, 0.4));
    }

    #[test]
    fn log10_powers_of_ten() {
        assert!(approx_eq(
            log10(I16F16::from_num(1.0)).unwrap(),
            0.0,
            TOLERANCE
        ));
        assert!(approx_eq(log10(I16F16::from_num(10.0)).unwrap(), 1.0, 0.25));
        assert!(approx_eq(
            log10(I16F16::from_num(100.0)).unwrap(),
            2.0,
            0.35
        ));
    }

    #[test]
    fn exp_ln_inverse() {
        // exp(ln(x)) ≈ x (with limited precision due to CORDIC)
        for i in 1..5 {
            let x = I16F16::from_num(i);
            let result = exp(ln(x).unwrap());
            let x_f32: f32 = x.to_num();
            let result_f32: f32 = result.to_num();
            assert!(
                (result_f32 - x_f32).abs() < 0.5,
                "exp(ln({x_f32})) = {result_f32}, expected {x_f32}"
            );
        }
    }

    #[test]
    fn pow2_log2_inverse() {
        // pow2(log2(x)) ≈ x for positive x
        for i in 1..8 {
            let x = I16F16::from_num(i);
            let result = pow2(log2(x).unwrap());
            let x_f32: f32 = x.to_num();
            let result_f32: f32 = result.to_num();
            assert!(
                (result_f32 - x_f32).abs() < 0.6,
                "pow2(log2({x_f32})) = {result_f32}, expected {x_f32}"
            );
        }
    }

    #[test]
    fn ln_exp_inverse() {
        // ln(exp(x)) ≈ x for small x (where exp doesn't overflow)
        for i in -3..=2 {
            let x = I16F16::from_num(i);
            let result = ln(exp(x));
            assert!(result.is_ok());
            let x_f32: f32 = x.to_num();
            let result_f32: f32 = result.unwrap().to_num();
            assert!(
                (result_f32 - x_f32).abs() < 0.3,
                "ln(exp({x_f32})) = {result_f32}, expected {x_f32}"
            );
        }
    }

    #[test]
    fn ln_near_zero() {
        // ln approaches -infinity as x approaches 0
        // Test small positive values
        let small = I16F16::from_num(0.01);
        let result = ln(small);
        assert!(result.is_ok());
        let val: f32 = result.unwrap().to_num();
        // ln(0.01) ≈ -4.605
        assert!(val < -3.0, "ln(0.01) = {val}, expected < -3.0");

        // Very small value
        let very_small = I16F16::from_num(0.001);
        let result2 = ln(very_small);
        assert!(result2.is_ok());
        let val2: f32 = result2.unwrap().to_num();
        // ln(0.001) ≈ -6.908
        assert!(val2 < -5.0, "ln(0.001) = {val2}, expected < -5.0");
    }

    #[test]
    fn ln_at_one() {
        // ln(1) should be exactly 0
        let result = ln(I16F16::ONE);
        assert!(result.is_ok());
        let val: f32 = result.unwrap().to_num();
        assert!(val.abs() < 0.01, "ln(1) = {val}, expected 0");
    }

    #[test]
    fn exp_large_negative() {
        // exp of large negative values should approach 0
        let neg_large = I16F16::from_num(-10.0);
        let result: f32 = exp(neg_large).to_num();
        // exp(-10) ≈ 0.0000454
        assert!(
            (0.0..0.01).contains(&result),
            "exp(-10) = {result}, expected ~0"
        );
    }

    #[test]
    fn exp_zero() {
        // exp(0) should be exactly 1
        let result: f32 = exp(I16F16::ZERO).to_num();
        assert!(
            (result - 1.0).abs() < 0.001,
            "exp(0) = {result}, expected 1"
        );
    }

    #[test]
    fn log2_at_one() {
        // log2(1) should be exactly 0
        let result = log2(I16F16::ONE);
        assert!(result.is_ok());
        let val: f32 = result.unwrap().to_num();
        assert!(val.abs() < 0.01, "log2(1) = {val}, expected 0");
    }

    #[test]
    fn log10_at_one() {
        // log10(1) should be exactly 0
        let result = log10(I16F16::ONE);
        assert!(result.is_ok());
        let val: f32 = result.unwrap().to_num();
        assert!(val.abs() < 0.01, "log10(1) = {val}, expected 0");
    }

    #[test]
    fn exp_large_positive() {
        // exp of large positive values should exercise the argument reduction loop
        let large = I16F16::from_num(5.0);
        let result: f32 = exp(large).to_num();
        // exp(5) ≈ 148.41
        assert!(
            result > 100.0 && result < 200.0,
            "exp(5) = {result}, expected ~148"
        );

        // Test even larger value to ensure multiple reduction iterations
        let larger = I16F16::from_num(8.0);
        let result2: f32 = exp(larger).to_num();
        // exp(8) ≈ 2981
        assert!(result2 > 2000.0, "exp(8) = {result2}, expected > 2000");
    }

    #[test]
    fn ln_large_values() {
        // ln of large values should exercise the argument reduction loop
        let large = I16F16::from_num(1000.0);
        let result = ln(large);
        assert!(result.is_ok());
        let val: f32 = result.unwrap().to_num();
        // ln(1000) ≈ 6.9
        assert!((val - 6.9).abs() < 0.5, "ln(1000) = {val}, expected ~6.9");

        // Test even larger value
        let larger = I16F16::from_num(10000.0);
        let result2 = ln(larger);
        assert!(result2.is_ok());
        let val2: f32 = result2.unwrap().to_num();
        // ln(10000) ≈ 9.2
        assert!(
            (val2 - 9.2).abs() < 0.5,
            "ln(10000) = {val2}, expected ~9.2"
        );
    }

    #[test]
    fn ln_very_small_values() {
        // Test very small values to ensure the multiply-by-2 loop iterates
        let tiny = I16F16::from_num(0.0001);
        let result = ln(tiny);
        assert!(result.is_ok());
        let val: f32 = result.unwrap().to_num();
        // ln(0.0001) ≈ -9.2
        assert!(val < -8.0, "ln(0.0001) = {val}, expected < -8.0");
    }

    #[test]
    fn exp_overflow_to_max() {
        // exp of very large positive values should return max when scale > max_shift
        // For I16F16: total_bits = 32, max_shift = 31
        // Need scale > 31, i.e., x > 31 * ln(2) ≈ 21.5
        // exp(25) triggers scale=36 which exceeds max_shift=31
        let very_large = I16F16::from_num(25.0);
        let result: f32 = exp(very_large).to_num();
        // Should return max value
        let max: f32 = I16F16::MAX.to_num();
        assert!(
            (result - max).abs() < 1.0,
            "exp(25) = {result}, expected max {max}"
        );
    }

    #[test]
    fn exp_underflow_to_zero() {
        // exp of very large negative values should return zero when -scale > max_shift
        // For I16F16: total_bits = 32, max_shift = 31
        // Need -scale > 31, i.e., x < -31 * ln(2) ≈ -21.5
        // exp(-25) triggers scale=-36 which exceeds -max_shift=-31
        let very_negative = I16F16::from_num(-25.0);
        let result: f32 = exp(very_negative).to_num();
        // Should return zero
        assert!(result == 0.0, "exp(-25) = {result}, expected 0");
    }

    mod saturation {
        use super::*;
        use fixed::types::I32F32;

        /// Check if I16F16 value is saturated to MAX (within 0.01%)
        fn is_max_16(val: I16F16) -> bool {
            val.to_num::<f32>() >= I16F16::MAX.to_num::<f32>() * 0.9999
        }

        /// Check if I16F16 value is saturated to zero
        fn is_zero_16(val: I16F16) -> bool {
            val == I16F16::ZERO || val.to_num::<f32>().abs() < 0.0001
        }

        /// Check if I32F32 value is saturated to MAX (within 0.01%)
        fn is_max_32(val: I32F32) -> bool {
            val.to_num::<f64>() >= I32F32::MAX.to_num::<f64>() * 0.9999
        }

        /// Check if I32F32 value is saturated to zero
        fn is_zero_32(val: I32F32) -> bool {
            val == I32F32::ZERO || val.to_num::<f64>().abs() < 0.000_000_1
        }

        // ===== exp saturation thresholds =====
        // I16F16: saturates to MAX at x >= 22.2, to zero at x <= -9.2
        // I32F32: saturates to MAX at x >= 44.4, to zero at x <= -16.2

        #[test]
        fn exp_i16f16_upper_threshold() {
            // Below threshold: should NOT saturate
            assert!(
                !is_max_16(exp(I16F16::from_num(22.1))),
                "exp(22.1) should not saturate"
            );
            // At threshold: should saturate
            assert!(
                is_max_16(exp(I16F16::from_num(22.2))),
                "exp(22.2) should saturate to MAX"
            );
        }

        #[test]
        fn exp_i16f16_lower_threshold() {
            // Above threshold: should NOT be zero
            assert!(
                !is_zero_16(exp(I16F16::from_num(-9.1))),
                "exp(-9.1) should not be zero"
            );
            // At threshold: should be zero
            assert!(
                is_zero_16(exp(I16F16::from_num(-9.2))),
                "exp(-9.2) should be zero"
            );
        }

        #[test]
        fn exp_i32f32_upper_threshold() {
            // Below threshold: should NOT saturate
            assert!(
                !is_max_32(exp(I32F32::from_num(44.3))),
                "exp(44.3) should not saturate"
            );
            // At threshold: should saturate
            assert!(
                is_max_32(exp(I32F32::from_num(44.4))),
                "exp(44.4) should saturate to MAX"
            );
        }

        #[test]
        fn exp_i32f32_lower_threshold() {
            // Above threshold: should NOT be zero
            assert!(
                !is_zero_32(exp(I32F32::from_num(-16.1))),
                "exp(-16.1) should not be zero"
            );
            // At threshold: should be zero
            assert!(
                is_zero_32(exp(I32F32::from_num(-16.2))),
                "exp(-16.2) should be zero"
            );
        }

        // ===== pow2 saturation thresholds =====
        // I16F16: saturates to MAX at x >= 15.0, to zero at x <= -13.2
        // I32F32: saturates to MAX at x >= 31.0, to zero at x <= -23.3

        #[test]
        fn pow2_i16f16_upper_threshold() {
            // Below threshold: should NOT saturate
            assert!(
                !is_max_16(pow2(I16F16::from_num(14.9))),
                "pow2(14.9) should not saturate"
            );
            // At threshold: should saturate
            assert!(
                is_max_16(pow2(I16F16::from_num(15.0))),
                "pow2(15.0) should saturate to MAX"
            );
        }

        #[test]
        fn pow2_i16f16_lower_threshold() {
            // Above threshold: should NOT be zero
            assert!(
                !is_zero_16(pow2(I16F16::from_num(-13.1))),
                "pow2(-13.1) should not be zero"
            );
            // At threshold: should be zero
            assert!(
                is_zero_16(pow2(I16F16::from_num(-13.2))),
                "pow2(-13.2) should be zero"
            );
        }

        #[test]
        fn pow2_i32f32_upper_threshold() {
            // Below threshold: should NOT saturate
            assert!(
                !is_max_32(pow2(I32F32::from_num(30.9))),
                "pow2(30.9) should not saturate"
            );
            // At threshold: should saturate
            assert!(
                is_max_32(pow2(I32F32::from_num(31.0))),
                "pow2(31.0) should saturate to MAX"
            );
        }

        #[test]
        fn pow2_i32f32_lower_threshold() {
            // Above threshold: should NOT be zero
            assert!(
                !is_zero_32(pow2(I32F32::from_num(-23.2))),
                "pow2(-23.2) should not be zero"
            );
            // At threshold: should be zero
            assert!(
                is_zero_32(pow2(I32F32::from_num(-23.3))),
                "pow2(-23.3) should be zero"
            );
        }
    }
}
