#[cfg(test)]
mod tests {
    use crate::unit::support::TestResult;
    use fixed::types::{I16F16, I32F32};
    use fixed_analytics::{acos, asin, atan, atan2, cos, sin, sin_cos, tan};

    const TOLERANCE: f32 = 0.002;

    fn approx_eq(a: I16F16, b: f32) -> bool {
        (a.to_num::<f32>() - b).abs() < TOLERANCE
    }

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

    #[test]
    fn sin_cos_pythagorean_identity() {
        for i in -20i16..=20 {
            let angle = I16F16::from_num(i) * I16F16::from_num(0.1);
            let (s, c) = sin_cos(angle);
            let sum_sq: f32 = (s * s + c * c).to_num();
            let angle_val = f32::from(i) * 0.1;
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

    #[test]
    fn atan2_quadrants() {
        let one = I16F16::ONE;
        let neg_one = -I16F16::ONE;

        let q1 = atan2(one, one);
        assert!(q1 > I16F16::ZERO && q1 < I16F16::FRAC_PI_2);

        let q2 = atan2(one, neg_one);
        assert!(q2 > I16F16::FRAC_PI_2);

        let q3 = atan2(neg_one, neg_one);
        assert!(q3 < -I16F16::FRAC_PI_2);

        let q4 = atan2(neg_one, one);
        assert!(q4 < I16F16::ZERO && q4 > -I16F16::FRAC_PI_2);
    }

    #[test]
    fn atan2_precise_quadrant_values() {
        let q1: f32 = atan2(I16F16::ONE, I16F16::ONE).to_num();
        let expected_q1 = core::f32::consts::FRAC_PI_4;
        assert!(
            (q1 - expected_q1).abs() < TOLERANCE,
            "Q1: expected {expected_q1}, got {q1}"
        );

        let q2: f32 = atan2(I16F16::ONE, -I16F16::ONE).to_num();
        let expected_q2 = core::f32::consts::PI - core::f32::consts::FRAC_PI_4;
        assert!(
            (q2 - expected_q2).abs() < TOLERANCE,
            "Q2: expected {expected_q2}, got {q2}"
        );

        let q3: f32 = atan2(-I16F16::ONE, -I16F16::ONE).to_num();
        let expected_q3 = -core::f32::consts::PI + core::f32::consts::FRAC_PI_4;
        assert!(
            (q3 - expected_q3).abs() < TOLERANCE,
            "Q3: expected {expected_q3}, got {q3}"
        );

        let q4: f32 = atan2(-I16F16::ONE, I16F16::ONE).to_num();
        let expected_q4 = -core::f32::consts::FRAC_PI_4;
        assert!(
            (q4 - expected_q4).abs() < TOLERANCE,
            "Q4: expected {expected_q4}, got {q4}"
        );
    }

    #[test]
    fn atan2_axis_values() {
        assert!(approx_eq(atan2(I16F16::ZERO, I16F16::ONE), 0.0));

        assert!(approx_eq(
            atan2(I16F16::ZERO, -I16F16::ONE),
            core::f32::consts::PI
        ));

        assert!(approx_eq(
            atan2(I16F16::ONE, I16F16::ZERO),
            core::f32::consts::FRAC_PI_2
        ));

        assert!(approx_eq(
            atan2(-I16F16::ONE, I16F16::ZERO),
            -core::f32::consts::FRAC_PI_2
        ));

        // atan2(0, 0) is undefined; the crate returns 0.
        assert!(approx_eq(atan2(I16F16::ZERO, I16F16::ZERO), 0.0));
    }

    #[test]
    fn asin_domain_check() {
        assert!(asin(I16F16::from_num(1.5)).is_err());
        assert!(asin(I16F16::from_num(-1.5)).is_err());
        assert!(asin(I16F16::from_num(0.5)).is_ok());
    }

    #[test]
    fn acos_special_values() -> TestResult {
        let result: f32 = acos(I16F16::ONE)?.to_num();
        assert!(result.abs() < 0.01);
        Ok(())
    }

    #[test]
    fn sin_asin_roundtrip() -> TestResult {
        for i in -9..=9 {
            let x = I16F16::from_num(i) * I16F16::from_num(0.1);
            let result = sin(asin(x)?);
            let x_f32: f32 = x.to_num();
            let result_f32: f32 = result.to_num();
            assert!(
                (result_f32 - x_f32).abs() < 0.02,
                "sin(asin({x_f32})) = {result_f32}, expected {x_f32}"
            );
        }
        Ok(())
    }

    #[test]
    fn cos_acos_roundtrip() -> TestResult {
        for i in -9..=9 {
            let x = I16F16::from_num(i) * I16F16::from_num(0.1);
            let result = cos(acos(x)?);
            let x_f32: f32 = x.to_num();
            let result_f32: f32 = result.to_num();
            assert!(
                (result_f32 - x_f32).abs() < 0.02,
                "cos(acos({x_f32})) = {result_f32}, expected {x_f32}"
            );
        }
        Ok(())
    }

    #[test]
    fn tan_atan_roundtrip() {
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
    fn asin_boundary_values() -> TestResult {
        let result_pos: f32 = asin(I16F16::ONE)?.to_num();
        let expected_pos = core::f32::consts::FRAC_PI_2;
        assert!(
            (result_pos - expected_pos).abs() < 0.01,
            "asin(1) = {result_pos}, expected {expected_pos}"
        );

        let result_neg: f32 = asin(-I16F16::ONE)?.to_num();
        let expected_neg = -core::f32::consts::FRAC_PI_2;
        assert!(
            (result_neg - expected_neg).abs() < 0.01,
            "asin(-1) = {result_neg}, expected {expected_neg}"
        );

        let near_one = I16F16::from_num(0.999);
        assert!(asin(near_one).is_ok());
        let result: f32 = asin(near_one)?.to_num();
        assert!(
            result > 1.4 && result < 1.6,
            "asin(0.999) should be near π/2"
        );
        Ok(())
    }

    #[test]
    fn acos_boundary_values() -> TestResult {
        let result_pos: f32 = acos(I16F16::ONE)?.to_num();
        assert!(
            result_pos.abs() < 0.01,
            "acos(1) = {result_pos}, expected 0"
        );

        let result_neg: f32 = acos(-I16F16::ONE)?.to_num();
        let expected_neg = core::f32::consts::PI;
        assert!(
            (result_neg - expected_neg).abs() < 0.01,
            "acos(-1) = {result_neg}, expected {expected_neg}"
        );
        Ok(())
    }

    #[test]
    fn sin_cos_large_angles() {
        let large_angle = I16F16::from_num(100.0);
        let (s, c) = sin_cos(large_angle);

        let sum_sq: f32 = (s * s + c * c).to_num();
        assert!(
            (sum_sq - 1.0).abs() < 0.05,
            "sin²(100) + cos²(100) = {sum_sq}, expected ~1.0"
        );

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
    fn asin_near_negative_one() -> TestResult {
        let near_neg_one = I16F16::from_num(-0.9999);
        let result = asin(near_neg_one);
        assert!(result.is_ok());
        let val: f32 = result?.to_num();
        let expected = -core::f32::consts::FRAC_PI_2;
        assert!(
            (val - expected).abs() < 0.1,
            "asin(-0.9999) = {val}, expected ~{expected}"
        );

        let neg_one = -I16F16::ONE;
        let result_exact = asin(neg_one);
        assert!(result_exact.is_ok());
        let val_exact: f32 = result_exact?.to_num();
        assert!(
            (val_exact - expected).abs() < 0.01,
            "asin(-1) = {val_exact}, expected ~{expected}"
        );
        Ok(())
    }

    #[test]
    fn asin_boundary_very_close() -> TestResult {
        // The boundary path needs √(1 - x²) < 2^-15, so x within 2^-31 of -1.
        let near_neg_one = -I32F32::ONE + I32F32::from_bits(1);
        let result = asin(near_neg_one);
        assert!(result.is_ok());
        let val: f64 = result?.to_num();
        let expected = -core::f64::consts::FRAC_PI_2;
        assert!(
            (val - expected).abs() < 0.01,
            "asin(~-1) = {val}, expected ~{expected}"
        );

        let near_pos_one = I32F32::ONE - I32F32::from_bits(1);
        let result_pos = asin(near_pos_one);
        assert!(result_pos.is_ok());
        let val_pos: f64 = result_pos?.to_num();
        let expected_pos = core::f64::consts::FRAC_PI_2;
        assert!(
            (val_pos - expected_pos).abs() < 0.01,
            "asin(~1) = {val_pos}, expected ~{expected_pos}"
        );
        Ok(())
    }

    #[test]
    fn asin_domain_error_message() {
        let too_large = I16F16::from_num(1.5);
        assert!(matches!(
            asin(too_large),
            Err(fixed_analytics::Error::DomainError { .. })
        ));

        let too_small = I16F16::from_num(-1.5);
        assert!(matches!(
            asin(too_small),
            Err(fixed_analytics::Error::DomainError { .. })
        ));
    }

    #[test]
    fn sin_cos_very_large_positive_angle() {
        let very_large = I16F16::from_num(200.0);
        let (s, c) = sin_cos(very_large);

        let sum_sq: f32 = (s * s + c * c).to_num();
        assert!(
            (sum_sq - 1.0).abs() < 0.1,
            "sin²(200) + cos²(200) = {sum_sq}, expected ~1.0"
        );
    }

    #[test]
    fn sin_cos_very_large_negative_angle() {
        let very_large_neg = I16F16::from_num(-200.0);
        let (s, c) = sin_cos(very_large_neg);

        let sum_sq: f32 = (s * s + c * c).to_num();
        assert!(
            (sum_sq - 1.0).abs() < 0.1,
            "sin²(-200) + cos²(-200) = {sum_sq}, expected ~1.0"
        );
    }

    #[test]
    fn sin_cos_extreme_angles() {
        // Iterative reduction once returned garbage near MAX.
        let extreme = I16F16::MAX - I16F16::from_num(1.0);
        let (s, c) = sin_cos(extreme);

        let sum_sq: f32 = (s * s + c * c).to_num();
        assert!(
            (sum_sq - 1.0).abs() < 0.05,
            "sin²({}) + cos²({}) = {}, expected ~1.0",
            extreme.to_num::<f32>(),
            extreme.to_num::<f32>(),
            sum_sq
        );

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
        let large = I16F16::from_num(1000.0);
        let (s, c) = sin_cos(large);

        let expected_sin = 1000.0_f64.sin();
        let expected_cos = 1000.0_f64.cos();

        assert!(
            (s.to_num::<f64>() - expected_sin).abs() < 0.01,
            "sin(1000) = {}, expected {}",
            s.to_num::<f64>(),
            expected_sin
        );
        assert!(
            (c.to_num::<f64>() - expected_cos).abs() < 0.01,
            "cos(1000) = {}, expected {}",
            c.to_num::<f64>(),
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
            // tan saturates within about 3 ulps of π/2.
            let far_from_pole = I16F16::from_num(FRAC_PI_2 - 0.0001);
            let near_pole = I16F16::from_num(FRAC_PI_2 - 0.00003);

            assert!(
                !is_max_16(tan(far_from_pole)),
                "tan(π/2 - 0.0001) should not saturate"
            );
            assert!(
                is_max_16(tan(near_pole)),
                "tan(π/2 - 0.00003) should saturate to MAX"
            );
        }

        #[test]
        fn tan_i16f16_near_negative_pole() {
            // Saturation is asymmetric, since -π/2 rounds differently from π/2.
            let far_from_pole = I16F16::from_num(-FRAC_PI_2 - 0.0001);
            let near_pole = I16F16::from_num(-FRAC_PI_2 - 0.00005);

            let at_neg_pole_offset = tan(far_from_pole);
            let closer_to_neg_pole = tan(near_pole);

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
            use core::f64::consts::PI;
            let pole_3pi2 = 3.0 * PI / 2.0;

            let far_from_pole = I16F16::from_num(pole_3pi2 - 0.0001);
            let near_pole = I16F16::from_num(pole_3pi2 - 0.00005);

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

#[cfg(test)]
mod reduction_and_accuracy {
    use crate::unit::support::Lcg;
    use crate::unit::support::TestResult;
    use fixed::types::{I24F8, I32F32, I64F64};
    use fixed_analytics::{CordicNumber, acos, asin, atan, atan2, sin_cos};

    fn reduced_f64<T: CordicNumber + fixed::traits::Fixed>(angle: T) -> f64 {
        let two_pi: f64 = (T::pi() + T::pi()).to_num();
        let a: f64 = angle.to_num();
        (a / two_pi).round().mul_add(-two_pi, a)
    }

    #[test]
    fn sin_cos_is_correct_within_pi_of_max() {
        // n·2π exceeds MAX here; saturation used to return (0, 1).
        for angle in [I32F32::MAX, I32F32::MIN, I32F32::MAX - I32F32::ONE] {
            let (s, c) = sin_cos(angle);
            let r = reduced_f64(angle);
            let (s, c): (f64, f64) = (s.to_num(), c.to_num());
            assert!(
                (s - r.sin()).abs() < 1e-5,
                "sin({angle}) = {s}, want {}",
                r.sin()
            );
            assert!(
                (c - r.cos()).abs() < 1e-5,
                "cos({angle}) = {c}, want {}",
                r.cos()
            );
        }
        for angle in [I24F8::MAX, I24F8::MIN] {
            let (s, c) = sin_cos(angle);
            let r = reduced_f64(angle);
            let (s, c): (f64, f64) = (s.to_num(), c.to_num());
            assert!(
                (s - r.sin()).abs() < 0.02,
                "sin({angle}) = {s}, want {}",
                r.sin()
            );
            assert!(
                (c - r.cos()).abs() < 0.02,
                "cos({angle}) = {c}, want {}",
                r.cos()
            );
        }
    }

    #[test]
    fn sin_cos_large_angles_i64f64() {
        let mut rng = Lcg(0x51);
        for _ in 0..500 {
            let v = rng.range(-1e6, 1e6);
            let angle = I64F64::from_num(v);
            let (s, c) = sin_cos(angle);
            let r = reduced_f64(angle);
            let (s, c): (f64, f64) = (s.to_num(), c.to_num());
            assert!(
                (s - r.sin()).abs() < 1e-9,
                "sin({v}) = {s}, want {}",
                r.sin()
            );
            assert!(
                (c - r.cos()).abs() < 1e-9,
                "cos({v}) = {c}, want {}",
                r.cos()
            );
        }
    }

    #[test]
    fn inverse_functions_track_f64_at_i64f64() -> TestResult {
        let mut rng = Lcg(0xA7A);
        for _ in 0..2000 {
            let v = rng.range(-1.0, 1.0);
            let x = I64F64::from_num(v);
            let asin_v: f64 = asin(x)?.to_num();
            assert!((asin_v - v.asin()).abs() < 2e-15, "asin({v}) = {asin_v}");
            let acos_v: f64 = acos(x)?.to_num();
            assert!((acos_v - v.acos()).abs() < 2e-15, "acos({v}) = {acos_v}");
            let big = (rng.range(-30.0, 30.0)).exp2() * v.signum();
            let atan_v: f64 = atan(I64F64::from_num(big)).to_num();
            assert!(
                (atan_v - big.atan()).abs() < 1e-15,
                "atan({big}) = {atan_v}"
            );
            let y = rng.range(-5.0, 5.0);
            let x2 = rng.range(-5.0, 5.0);
            let atan2_v: f64 = atan2(I64F64::from_num(y), I64F64::from_num(x2)).to_num();
            assert!(
                (atan2_v - y.atan2(x2)).abs() < 1e-15,
                "atan2({y}, {x2}) = {atan2_v}"
            );
        }
        Ok(())
    }
}
