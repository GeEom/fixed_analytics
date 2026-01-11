//! Tests for circular trigonometric functions

#[cfg(test)]
#[allow(
    clippy::unwrap_used,
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    reason = "test code uses unwrap and f32/f64 casts for conciseness"
)]
mod tests {
    use fixed::types::{I16F16, I32F32};
    use fixed_analytics::{acos, asin, atan, atan2, cos, sin, sin_cos, tan};

    const TOLERANCE: f32 = 0.002;

    fn approx_eq(a: I16F16, b: f32) -> bool {
        (a.to_num::<f32>() - b).abs() < TOLERANCE
    }

    // Tests sin(0) = 0, sin(π/2) = 1, sin(-π/2) = -1, sin(π) = 0
    #[test]
    fn sin_special_values() {
        assert!(approx_eq(sin(I16F16::ZERO), 0.0));
        assert!(approx_eq(sin(I16F16::FRAC_PI_2), 1.0));
        assert!(approx_eq(sin(-I16F16::FRAC_PI_2), -1.0));
        assert!(approx_eq(sin(I16F16::PI), 0.0));
    }

    #[test]
    fn cos_special_values() {
        assert!(approx_eq(cos(I16F16::ZERO), 1.0));
        assert!(approx_eq(cos(I16F16::FRAC_PI_2), 0.0));
        assert!(approx_eq(cos(I16F16::PI), -1.0));
    }

    #[test]
    fn tan_special_values() {
        assert!(approx_eq(tan(I16F16::ZERO), 0.0));
        assert!(approx_eq(tan(I16F16::FRAC_PI_4), 1.0));
    }

    // Tests Pythagorean identity: sin²(x) + cos²(x) = 1
    #[test]
    fn sin_cos_pythagorean_identity() {
        for i in -20..=20 {
            let angle = I16F16::from_num(i) * I16F16::from_num(0.1);
            let (s, c) = sin_cos(angle);
            let sum_sq: f32 = (s * s + c * c).to_num();
            let angle_val = i as f32 * 0.1;
            assert!(
                (sum_sq - 1.0).abs() < 0.02,
                "sin²({angle_val}) + cos²({angle_val}) = {sum_sq}, expected ~1.0"
            );
        }
    }

    #[test]
    fn atan_special_values() {
        assert!(approx_eq(atan(I16F16::ZERO), 0.0));
        assert!(approx_eq(atan(I16F16::ONE), core::f32::consts::FRAC_PI_4));
        assert!(approx_eq(atan(-I16F16::ONE), -core::f32::consts::FRAC_PI_4));
    }

    // Tests atan2 quadrant correctness
    #[test]
    fn atan2_quadrants() {
        let one = I16F16::ONE;
        let neg_one = -I16F16::ONE;

        // Q1: both positive
        let q1 = atan2(one, one);
        assert!(q1 > I16F16::ZERO && q1 < I16F16::FRAC_PI_2);

        // Q2: y positive, x negative
        let q2 = atan2(one, neg_one);
        assert!(q2 > I16F16::FRAC_PI_2);

        // Q3: both negative
        let q3 = atan2(neg_one, neg_one);
        assert!(q3 < -I16F16::FRAC_PI_2);

        // Q4: y negative, x positive
        let q4 = atan2(neg_one, one);
        assert!(q4 < I16F16::ZERO && q4 > -I16F16::FRAC_PI_2);
    }

    #[test]
    fn atan2_precise_quadrant_values() {
        // Test precise values for each quadrant
        // atan2(1, 1) should be π/4 ≈ 0.785
        let q1: f32 = atan2(I16F16::ONE, I16F16::ONE).to_num();
        let expected_q1 = core::f32::consts::FRAC_PI_4;
        assert!(
            (q1 - expected_q1).abs() < TOLERANCE,
            "Q1: expected {expected_q1}, got {q1}"
        );

        // atan2(1, -1) should be 3π/4 ≈ 2.356
        let q2: f32 = atan2(I16F16::ONE, -I16F16::ONE).to_num();
        let expected_q2 = core::f32::consts::PI - core::f32::consts::FRAC_PI_4;
        assert!(
            (q2 - expected_q2).abs() < TOLERANCE,
            "Q2: expected {expected_q2}, got {q2}"
        );

        // atan2(-1, -1) should be -3π/4 ≈ -2.356
        let q3: f32 = atan2(-I16F16::ONE, -I16F16::ONE).to_num();
        let expected_q3 = -core::f32::consts::PI + core::f32::consts::FRAC_PI_4;
        assert!(
            (q3 - expected_q3).abs() < TOLERANCE,
            "Q3: expected {expected_q3}, got {q3}"
        );

        // atan2(-1, 1) should be -π/4 ≈ -0.785
        let q4: f32 = atan2(-I16F16::ONE, I16F16::ONE).to_num();
        let expected_q4 = -core::f32::consts::FRAC_PI_4;
        assert!(
            (q4 - expected_q4).abs() < TOLERANCE,
            "Q4: expected {expected_q4}, got {q4}"
        );
    }

    #[test]
    fn atan2_axis_values() {
        // Test values along the axes
        // atan2(0, 1) = 0
        assert!(approx_eq(atan2(I16F16::ZERO, I16F16::ONE), 0.0));

        // atan2(0, -1) = π
        assert!(approx_eq(
            atan2(I16F16::ZERO, -I16F16::ONE),
            core::f32::consts::PI
        ));

        // atan2(1, 0) = π/2
        assert!(approx_eq(
            atan2(I16F16::ONE, I16F16::ZERO),
            core::f32::consts::FRAC_PI_2
        ));

        // atan2(-1, 0) = -π/2
        assert!(approx_eq(
            atan2(-I16F16::ONE, I16F16::ZERO),
            -core::f32::consts::FRAC_PI_2
        ));

        // atan2(0, 0) = 0 (undefined but we return 0)
        assert!(approx_eq(atan2(I16F16::ZERO, I16F16::ZERO), 0.0));
    }

    // Tests asin domain validation rejects |x| > 1
    #[test]
    fn asin_domain_check() {
        assert!(asin(I16F16::from_num(1.5)).is_err());
        assert!(asin(I16F16::from_num(-1.5)).is_err());
        assert!(asin(I16F16::from_num(0.5)).is_ok());
    }

    #[test]
    fn acos_special_values() {
        let result: f32 = acos(I16F16::ONE).unwrap().to_num();
        assert!(result.abs() < 0.01);
    }

    #[test]
    fn sin_asin_roundtrip() {
        // sin(asin(x)) ≈ x for x in [-1, 1]
        for i in -9..=9 {
            let x = I16F16::from_num(i) * I16F16::from_num(0.1);
            let result = sin(asin(x).unwrap());
            let x_f32: f32 = x.to_num();
            let result_f32: f32 = result.to_num();
            assert!(
                (result_f32 - x_f32).abs() < 0.02,
                "sin(asin({x_f32})) = {result_f32}, expected {x_f32}"
            );
        }
    }

    #[test]
    fn cos_acos_roundtrip() {
        // cos(acos(x)) ≈ x for x in [-1, 1]
        for i in -9..=9 {
            let x = I16F16::from_num(i) * I16F16::from_num(0.1);
            let result = cos(acos(x).unwrap());
            let x_f32: f32 = x.to_num();
            let result_f32: f32 = result.to_num();
            assert!(
                (result_f32 - x_f32).abs() < 0.02,
                "cos(acos({x_f32})) = {result_f32}, expected {x_f32}"
            );
        }
    }

    #[test]
    fn tan_atan_roundtrip() {
        // tan(atan(x)) ≈ x for various x
        for i in -10..=10 {
            let x = I16F16::from_num(i) * I16F16::from_num(0.5);
            let result = tan(atan(x));
            let x_f32: f32 = x.to_num();
            let result_f32: f32 = result.to_num();
            assert!(
                (result_f32 - x_f32).abs() < 0.05,
                "tan(atan({x_f32})) = {result_f32}, expected {x_f32}"
            );
        }
    }

    #[test]
    fn asin_boundary_values() {
        // asin at exactly ±1 should return ±π/2
        let result_pos: f32 = asin(I16F16::ONE).unwrap().to_num();
        let expected_pos = core::f32::consts::FRAC_PI_2;
        assert!(
            (result_pos - expected_pos).abs() < 0.01,
            "asin(1) = {result_pos}, expected {expected_pos}"
        );

        let result_neg: f32 = asin(-I16F16::ONE).unwrap().to_num();
        let expected_neg = -core::f32::consts::FRAC_PI_2;
        assert!(
            (result_neg - expected_neg).abs() < 0.01,
            "asin(-1) = {result_neg}, expected {expected_neg}"
        );

        // Near boundary (0.999)
        let near_one = I16F16::from_num(0.999);
        assert!(asin(near_one).is_ok());
        let result: f32 = asin(near_one).unwrap().to_num();
        assert!(
            result > 1.4 && result < 1.6,
            "asin(0.999) should be near π/2"
        );
    }

    #[test]
    fn acos_boundary_values() {
        // acos at exactly ±1 should return 0 or π
        let result_pos: f32 = acos(I16F16::ONE).unwrap().to_num();
        assert!(
            result_pos.abs() < 0.01,
            "acos(1) = {result_pos}, expected 0"
        );

        let result_neg: f32 = acos(-I16F16::ONE).unwrap().to_num();
        let expected_neg = core::f32::consts::PI;
        assert!(
            (result_neg - expected_neg).abs() < 0.01,
            "acos(-1) = {result_neg}, expected {expected_neg}"
        );
    }

    #[test]
    fn sin_cos_large_angles() {
        // Test angle reduction with large angles
        let large_angle = I16F16::from_num(100.0); // Well beyond 2π
        let (s, c) = sin_cos(large_angle);

        // Should still satisfy sin²+cos² = 1
        let sum_sq: f32 = (s * s + c * c).to_num();
        assert!(
            (sum_sq - 1.0).abs() < 0.05,
            "sin²(100) + cos²(100) = {sum_sq}, expected ~1.0"
        );

        // Test negative large angle
        let neg_large_angle = I16F16::from_num(-100.0);
        let (s2, c2) = sin_cos(neg_large_angle);
        let sum_sq2: f32 = (s2 * s2 + c2 * c2).to_num();
        assert!(
            (sum_sq2 - 1.0).abs() < 0.05,
            "sin²(-100) + cos²(-100) = {sum_sq2}, expected ~1.0"
        );
    }

    #[test]
    fn atan_large_values() {
        // atan of large values should approach ±π/2
        let large = I16F16::from_num(1000.0);
        let result: f32 = atan(large).to_num();
        let expected = core::f32::consts::FRAC_PI_2;
        assert!(
            (result - expected).abs() < 0.01,
            "atan(1000) = {result}, expected ~{expected}"
        );

        let neg_large = I16F16::from_num(-1000.0);
        let result_neg: f32 = atan(neg_large).to_num();
        let expected_neg = -core::f32::consts::FRAC_PI_2;
        assert!(
            (result_neg - expected_neg).abs() < 0.01,
            "atan(-1000) = {result_neg}, expected ~{expected_neg}"
        );
    }

    #[test]
    fn asin_near_negative_one() {
        // Test asin for values very close to -1, exercising the boundary path
        // where sqrt_term is very small and x is negative
        let near_neg_one = I16F16::from_num(-0.9999);
        let result = asin(near_neg_one);
        assert!(result.is_ok());
        let val: f32 = result.unwrap().to_num();
        // asin(-0.9999) should be very close to -π/2
        let expected = -core::f32::consts::FRAC_PI_2;
        assert!(
            (val - expected).abs() < 0.1,
            "asin(-0.9999) = {val}, expected ~{expected}"
        );

        // Test exactly at -1
        let neg_one = -I16F16::ONE;
        let result_exact = asin(neg_one);
        assert!(result_exact.is_ok());
        let val_exact: f32 = result_exact.unwrap().to_num();
        assert!(
            (val_exact - expected).abs() < 0.01,
            "asin(-1) = {val_exact}, expected ~{expected}"
        );
    }

    #[test]
    fn asin_boundary_very_close() {
        // Test with higher precision type to hit the small sqrt_term boundary path
        // The threshold is 2^(-15) ≈ 0.0000305, so we need sqrt(1-x²) < 0.0000305
        // This requires 1-x² < 2^(-30), meaning x must be within 2^(-31) of ±1
        // Using I32F32 for higher precision

        // Create a value extremely close to -1 but not exactly -1
        // -1 + 2^(-32) (smallest representable distance)
        let near_neg_one = -I32F32::ONE + I32F32::from_bits(1);
        let result = asin(near_neg_one);
        assert!(result.is_ok());
        let val: f64 = result.unwrap().to_num();
        let expected = -core::f64::consts::FRAC_PI_2;
        assert!(
            (val - expected).abs() < 0.01,
            "asin(~-1) = {val}, expected ~{expected}"
        );

        // Also test positive boundary
        let near_pos_one = I32F32::ONE - I32F32::from_bits(1);
        let result_pos = asin(near_pos_one);
        assert!(result_pos.is_ok());
        let val_pos: f64 = result_pos.unwrap().to_num();
        let expected_pos = core::f64::consts::FRAC_PI_2;
        assert!(
            (val_pos - expected_pos).abs() < 0.01,
            "asin(~1) = {val_pos}, expected ~{expected_pos}"
        );
    }

    #[test]
    fn asin_domain_error_message() {
        // Test that domain errors are properly returned for out-of-range values
        let too_large = I16F16::from_num(1.5);
        let err = asin(too_large).unwrap_err();
        assert!(matches!(err, fixed_analytics::Error::DomainError { .. }));

        let too_small = I16F16::from_num(-1.5);
        let err2 = asin(too_small).unwrap_err();
        assert!(matches!(err2, fixed_analytics::Error::DomainError { .. }));
    }

    #[test]
    fn sin_cos_very_large_positive_angle() {
        // Test with very large positive angle to exercise the reduction loop
        // The loop iterates while reduced > pi, subtracting 2π each time
        // 200.0 radians is about 32 times 2π, so the loop will iterate many times
        let very_large = I16F16::from_num(200.0);
        let (s, c) = sin_cos(very_large);

        // Should still satisfy sin²+cos² = 1
        let sum_sq: f32 = (s * s + c * c).to_num();
        assert!(
            (sum_sq - 1.0).abs() < 0.1,
            "sin²(200) + cos²(200) = {sum_sq}, expected ~1.0"
        );
    }

    #[test]
    fn sin_cos_very_large_negative_angle() {
        // Test with very large negative angle to exercise the second reduction loop
        // The loop iterates while reduced < -pi, adding 2π each time
        let very_large_neg = I16F16::from_num(-200.0);
        let (s, c) = sin_cos(very_large_neg);

        // Should still satisfy sin²+cos² = 1
        let sum_sq: f32 = (s * s + c * c).to_num();
        assert!(
            (sum_sq - 1.0).abs() < 0.1,
            "sin²(-200) + cos²(-200) = {sum_sq}, expected ~1.0"
        );
    }

    #[test]
    fn sin_cos_extreme_angles() {
        // Test near I16F16::MAX - previously would produce garbage with iterative reduction
        let extreme = I16F16::MAX - I16F16::from_num(1.0);
        let (s, c) = sin_cos(extreme);

        // Must satisfy Pythagorean identity regardless of input magnitude
        let sum_sq: f32 = (s * s + c * c).to_num();
        assert!(
            (sum_sq - 1.0).abs() < 0.05,
            "sin²({}) + cos²({}) = {}, expected ~1.0",
            extreme.to_num::<f32>(),
            extreme.to_num::<f32>(),
            sum_sq
        );

        // Test negative extreme
        let neg_extreme = I16F16::MIN + I16F16::from_num(1.0);
        let (s2, c2) = sin_cos(neg_extreme);
        let sum_sq2: f32 = (s2 * s2 + c2 * c2).to_num();
        assert!(
            (sum_sq2 - 1.0).abs() < 0.05,
            "sin²({}) + cos²({}) = {}, expected ~1.0",
            neg_extreme.to_num::<f32>(),
            neg_extreme.to_num::<f32>(),
            sum_sq2
        );
    }

    #[test]
    fn sin_cos_known_large_values() {
        // 1000 radians ≈ 159.15 * 2π, remainder ≈ 0.96 rad
        let large = I16F16::from_num(1000.0);
        let (s, c) = sin_cos(large);

        // Compare against f64 reference
        let expected_sin = 1000.0_f64.sin() as f32;
        let expected_cos = 1000.0_f64.cos() as f32;

        assert!(
            (s.to_num::<f32>() - expected_sin).abs() < 0.01,
            "sin(1000) = {}, expected {}",
            s.to_num::<f32>(),
            expected_sin
        );
        assert!(
            (c.to_num::<f32>() - expected_cos).abs() < 0.01,
            "cos(1000) = {}, expected {}",
            c.to_num::<f32>(),
            expected_cos
        );
    }

    mod saturation {
        use super::*;
        use core::f64::consts::FRAC_PI_2;

        fn is_max_16(val: I16F16) -> bool {
            val.to_num::<f32>() >= I16F16::MAX.to_num::<f32>() * 0.999
        }

        #[test]
        fn tan_i16f16_near_positive_pole() {
            // tan saturates to MAX when approaching π/2 from below
            // Threshold is approximately 0.00008 (8e-5)
            let far_from_pole = I16F16::from_num(FRAC_PI_2 - 0.0001);
            let near_pole = I16F16::from_num(FRAC_PI_2 - 0.00005);

            assert!(
                !is_max_16(tan(far_from_pole)),
                "tan(π/2 - 0.0001) should not saturate"
            );
            assert!(
                is_max_16(tan(near_pole)),
                "tan(π/2 - 0.00005) should saturate to MAX"
            );
        }

        #[test]
        fn tan_i16f16_near_negative_pole() {
            // tan saturates to MIN when approaching -π/2 from below (more negative)
            // The saturation is asymmetric due to fixed-point representation of -π/2
            let far_from_pole = I16F16::from_num(-FRAC_PI_2 - 0.0001);
            let near_pole = I16F16::from_num(-FRAC_PI_2 - 0.00005);

            // Approaching from below (more negative), tan goes to +∞ then wraps to MIN
            // Test the basic property: extreme values near the pole
            let at_neg_pole_offset = tan(far_from_pole);
            let closer_to_neg_pole = tan(near_pole);

            // Values should have large absolute magnitude near the pole
            assert!(
                at_neg_pole_offset.to_num::<f32>().abs() > 1000.0,
                "tan near -π/2 should have large magnitude"
            );
            assert!(
                closer_to_neg_pole.to_num::<f32>().abs() > at_neg_pole_offset.to_num::<f32>().abs(),
                "tan closer to pole should have larger magnitude"
            );
        }

        #[test]
        fn tan_i16f16_at_3pi_over_2() {
            // tan also has a pole at 3π/2
            use core::f64::consts::PI;
            let pole_3pi2 = 3.0 * PI / 2.0;

            let far_from_pole = I16F16::from_num(pole_3pi2 - 0.0001);
            let near_pole = I16F16::from_num(pole_3pi2 - 0.00005);

            // Near 3π/2, tan approaches +∞ from below
            assert!(
                !is_max_16(tan(far_from_pole)),
                "tan(3π/2 - 0.0001) should not saturate"
            );
            assert!(
                is_max_16(tan(near_pole)),
                "tan(3π/2 - 0.00005) should saturate to MAX"
            );
        }
    }
}
