#[cfg(test)]
mod tests {
    use crate::unit::support::TestResult;
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
    fn ln_special_values() -> TestResult {
        assert!(approx_eq(ln(I16F16::ONE)?, 0.0, TOLERANCE));
        let ln_e: f32 = ln(I16F16::E)?.to_num();
        assert!((ln_e - 1.0).abs() < 0.25, "ln(e) = {ln_e}, expected ~1.0");
        Ok(())
    }

    #[test]
    fn ln_domain_check() {
        assert!(ln(I16F16::ZERO).is_err());
        assert!(ln(I16F16::from_num(-1.0)).is_err());
        assert!(ln(I16F16::from_num(0.5)).is_ok());
    }

    #[test]
    fn log2_powers_of_two() -> TestResult {
        assert!(approx_eq(log2(I16F16::from_num(1.0))?, 0.0, TOLERANCE));
        assert!(approx_eq(log2(I16F16::from_num(2.0))?, 1.0, 0.25));
        assert!(approx_eq(log2(I16F16::from_num(4.0))?, 2.0, 0.3));
        assert!(approx_eq(log2(I16F16::from_num(8.0))?, 3.0, 0.4));
        Ok(())
    }

    #[test]
    fn log10_powers_of_ten() -> TestResult {
        assert!(approx_eq(log10(I16F16::from_num(1.0))?, 0.0, TOLERANCE));
        assert!(approx_eq(log10(I16F16::from_num(10.0))?, 1.0, 0.25));
        assert!(approx_eq(log10(I16F16::from_num(100.0))?, 2.0, 0.35));
        Ok(())
    }

    #[test]
    fn exp_ln_inverse() -> TestResult {
        for i in 1..5 {
            let x = I16F16::from_num(i);
            let result = exp(ln(x)?);
            let x_f32: f32 = x.to_num();
            let result_f32: f32 = result.to_num();
            assert!(
                (result_f32 - x_f32).abs() < 0.5,
                "exp(ln({x_f32})) = {result_f32}, expected {x_f32}"
            );
        }
        Ok(())
    }

    #[test]
    fn pow2_log2_inverse() -> TestResult {
        for i in 1..8 {
            let x = I16F16::from_num(i);
            let result = pow2(log2(x)?);
            let x_f32: f32 = x.to_num();
            let result_f32: f32 = result.to_num();
            assert!(
                (result_f32 - x_f32).abs() < 0.6,
                "pow2(log2({x_f32})) = {result_f32}, expected {x_f32}"
            );
        }
        Ok(())
    }

    #[test]
    fn ln_exp_inverse() -> TestResult {
        for i in -3..=2 {
            let x = I16F16::from_num(i);
            let result = ln(exp(x));
            assert!(result.is_ok());
            let x_f32: f32 = x.to_num();
            let result_f32: f32 = result?.to_num();
            assert!(
                (result_f32 - x_f32).abs() < 0.3,
                "ln(exp({x_f32})) = {result_f32}, expected {x_f32}"
            );
        }
        Ok(())
    }

    #[test]
    fn ln_near_zero() -> TestResult {
        let small = I16F16::from_num(0.01);
        let result = ln(small);
        assert!(result.is_ok());
        let val: f32 = result?.to_num();
        assert!(val < -3.0, "ln(0.01) = {val}, expected < -3.0");

        let very_small = I16F16::from_num(0.001);
        let result2 = ln(very_small);
        assert!(result2.is_ok());
        let val2: f32 = result2?.to_num();
        assert!(val2 < -5.0, "ln(0.001) = {val2}, expected < -5.0");
        Ok(())
    }

    #[test]
    fn ln_at_one() -> TestResult {
        let result = ln(I16F16::ONE);
        assert!(result.is_ok());
        let val: f32 = result?.to_num();
        assert!(val.abs() < 0.01, "ln(1) = {val}, expected 0");
        Ok(())
    }

    #[test]
    fn exp_large_negative() {
        let neg_large = I16F16::from_num(-10.0);
        let result: f32 = exp(neg_large).to_num();
        assert!(
            (0.0..0.01).contains(&result),
            "exp(-10) = {result}, expected ~0"
        );
    }

    #[test]
    fn exp_zero() {
        let result: f32 = exp(I16F16::ZERO).to_num();
        assert!(
            (result - 1.0).abs() < 0.001,
            "exp(0) = {result}, expected 1"
        );
    }

    #[test]
    fn log2_at_one() -> TestResult {
        let result = log2(I16F16::ONE);
        assert!(result.is_ok());
        let val: f32 = result?.to_num();
        assert!(val.abs() < 0.01, "log2(1) = {val}, expected 0");
        Ok(())
    }

    #[test]
    fn log10_at_one() -> TestResult {
        let result = log10(I16F16::ONE);
        assert!(result.is_ok());
        let val: f32 = result?.to_num();
        assert!(val.abs() < 0.01, "log10(1) = {val}, expected 0");
        Ok(())
    }

    #[test]
    fn exp_large_positive() {
        let large = I16F16::from_num(5.0);
        let result: f32 = exp(large).to_num();
        assert!(
            result > 100.0 && result < 200.0,
            "exp(5) = {result}, expected ~148"
        );

        let larger = I16F16::from_num(8.0);
        let result2: f32 = exp(larger).to_num();
        assert!(result2 > 2000.0, "exp(8) = {result2}, expected > 2000");
    }

    #[test]
    fn ln_large_values() -> TestResult {
        let large = I16F16::from_num(1000.0);
        let result = ln(large);
        assert!(result.is_ok());
        let val: f32 = result?.to_num();
        assert!((val - 6.9).abs() < 0.5, "ln(1000) = {val}, expected ~6.9");

        let larger = I16F16::from_num(10000.0);
        let result2 = ln(larger);
        assert!(result2.is_ok());
        let val2: f32 = result2?.to_num();
        assert!(
            (val2 - 9.2).abs() < 0.5,
            "ln(10000) = {val2}, expected ~9.2"
        );
        Ok(())
    }

    #[test]
    fn ln_very_small_values() -> TestResult {
        let tiny = I16F16::from_num(0.0001);
        let result = ln(tiny);
        assert!(result.is_ok());
        let val: f32 = result?.to_num();
        assert!(val < -8.0, "ln(0.0001) = {val}, expected < -8.0");
        Ok(())
    }

    #[test]
    fn exp_overflow_to_max() {
        let very_large = I16F16::from_num(25.0);
        let result: f32 = exp(very_large).to_num();
        let max: f32 = I16F16::MAX.to_num();
        assert!(
            (result - max).abs() < 1.0,
            "exp(25) = {result}, expected max {max}"
        );
    }

    #[test]
    fn exp_underflow_to_zero() {
        let very_negative = I16F16::from_num(-25.0);
        let result: f32 = exp(very_negative).to_num();
        assert!(result == 0.0, "exp(-25) = {result}, expected 0");
    }

    mod saturation {
        use super::*;
        use fixed::types::I32F32;

        fn is_max_16(val: I16F16) -> bool {
            val.to_num::<f32>() >= I16F16::MAX.to_num::<f32>() * 0.9999
        }

        fn is_zero_16(val: I16F16) -> bool {
            val == I16F16::ZERO
        }

        fn is_max_32(val: I32F32) -> bool {
            val.to_num::<f64>() >= I32F32::MAX.to_num::<f64>() * 0.9999
        }

        fn is_zero_32(val: I32F32) -> bool {
            val == I32F32::ZERO
        }

        // Zero once the result rounds below half an ulp.
        #[test]
        fn exp_i16f16_upper_threshold() {
            assert!(
                !is_max_16(exp(I16F16::from_num(10.39))),
                "exp(10.39) should not saturate"
            );
            assert!(
                is_max_16(exp(I16F16::from_num(10.4))),
                "exp(10.4) should saturate to MAX"
            );
        }

        #[test]
        fn exp_i16f16_lower_threshold() {
            assert!(
                !is_zero_16(exp(I16F16::from_num(-11.7))),
                "exp(-11.7) should not be zero"
            );
            assert!(
                is_zero_16(exp(I16F16::from_num(-11.8))),
                "exp(-11.8) should be zero"
            );
        }

        #[test]
        fn exp_deep_negative_band_is_zero() {
            // trunc(-21.6/ln2) = -31 for I16F16, and r < 0 then steps scale to -32.
            assert!(
                is_zero_16(exp(I16F16::from_num(-21.6))),
                "exp(-21.6) should be zero"
            );
            assert!(
                is_zero_32(exp(I32F32::from_num(-43.7))),
                "exp(-43.7) should be zero"
            );
        }

        #[test]
        fn exp_i32f32_upper_threshold() {
            assert!(
                !is_max_32(exp(I32F32::from_num(21.48))),
                "exp(21.48) should not saturate"
            );
            assert!(
                is_max_32(exp(I32F32::from_num(21.49))),
                "exp(21.49) should saturate to MAX"
            );
        }

        #[test]
        fn exp_i32f32_lower_threshold() {
            assert!(
                !is_zero_32(exp(I32F32::from_num(-22.8))),
                "exp(-22.8) should not be zero"
            );
            assert!(
                is_zero_32(exp(I32F32::from_num(-22.9))),
                "exp(-22.9) should be zero"
            );
        }

        #[test]
        fn pow2_i16f16_upper_threshold() {
            assert!(
                !is_max_16(pow2(I16F16::from_num(14.9))),
                "pow2(14.9) should not saturate"
            );
            assert!(
                is_max_16(pow2(I16F16::from_num(15.0))),
                "pow2(15.0) should saturate to MAX"
            );
        }

        #[test]
        fn pow2_i16f16_lower_threshold() {
            assert!(
                !is_zero_16(pow2(I16F16::from_num(-16.9))),
                "pow2(-16.9) should not be zero"
            );
            assert!(
                is_zero_16(pow2(I16F16::from_num(-17.1))),
                "pow2(-17.1) should be zero"
            );
        }

        #[test]
        fn pow2_i32f32_upper_threshold() {
            assert!(
                !is_max_32(pow2(I32F32::from_num(30.9))),
                "pow2(30.9) should not saturate"
            );
            assert!(
                is_max_32(pow2(I32F32::from_num(31.0))),
                "pow2(31.0) should saturate to MAX"
            );
        }

        #[test]
        fn pow2_i32f32_lower_threshold() {
            assert!(
                !is_zero_32(pow2(I32F32::from_num(-32.9))),
                "pow2(-32.9) should not be zero"
            );
            assert!(
                is_zero_32(pow2(I32F32::from_num(-33.1))),
                "pow2(-33.1) should be zero"
            );
        }
    }
}

#[cfg(test)]
mod reduction {
    use crate::unit::support::Lcg;
    use crate::unit::support::TestResult;
    use fixed::types::{I4F60, I16F16, I24F8, I32F32, I64F64};
    use fixed_analytics::{CordicNumber, exp, ln, log2, pow, pow2, sqrt};

    #[test]
    fn exp_saturates_at_the_extremes_of_wide_types() {
        // x / ln2 does not fit in i32 here; the conversion saturates.
        assert_eq!(exp(I32F32::MAX), I32F32::MAX);
        assert_eq!(exp(I32F32::MIN), I32F32::ZERO);
        assert_eq!(exp(I32F32::from_num(2e9)), I32F32::MAX);
        assert_eq!(exp(I64F64::MAX), I64F64::MAX);
        assert_eq!(exp(I64F64::MIN), I64F64::ZERO);
        assert_eq!(exp(I64F64::from_num(1.5e9)), I64F64::MAX);
        assert_eq!(exp(I64F64::from_num(-1.5e9)), I64F64::ZERO);
        assert_eq!(exp(I16F16::MAX), I16F16::MAX);
        assert_eq!(exp(I16F16::MIN), I16F16::ZERO);
        assert_eq!(pow2(I64F64::MAX), I64F64::MAX);
        assert_eq!(pow2(I64F64::MIN), I64F64::ZERO);
    }

    #[test]
    fn exp_thresholds_after_correction_i16f16() {
        // k reaches 15 (overflow) or -18 (zero) only after the correction.
        assert_eq!(exp(I16F16::from_num(10.5)), I16F16::MAX);
        assert!(exp(I16F16::from_num(10.39)) < I16F16::MAX);
        assert_eq!(exp(I16F16::from_num(-12.1)), I16F16::ZERO);
        assert!(exp(I16F16::from_num(-11.4)) > I16F16::ZERO);
    }

    #[test]
    fn exp_estimate_corrected_downward() {
        // Just below a negative multiple of ln2 the estimate rounds up; r < 0 steps k down.
        let x16 = I16F16::LN_2.mul_int(-3) - I16F16::from_bits(1);
        let got16: f64 = exp(x16).to_num();
        assert!((got16 - 0.125).abs() < 3e-5, "exp({x16}) = {got16}");
        let x64 = I64F64::LN_2.mul_int(-3) - I64F64::from_bits(1);
        let got64: f64 = exp(x64).to_num();
        assert!((got64 - 0.125).abs() < 1e-12, "exp({x64}) = {got64}");
    }

    #[test]
    fn exp_tracks_f64_at_i64f64() {
        let mut rng = Lcg(0xE4);
        for i in 0..2000 {
            let v = if i < 1000 {
                -40.0 + 80.0 * f64::from(i) / 1000.0
            } else {
                rng.range(-8.0, 8.0)
            };
            let x = I64F64::from_num(v);
            let got: f64 = exp(x).to_num();
            let want = x.to_num::<f64>().exp();
            // The f64 reference loses 2^-53·|x| of relative precision in x.
            let tol = 8.0f64.mul_add(2f64.powi(-64), 5e-16 * (1.0 + v.abs()) * want);
            assert!((got - want).abs() < tol, "exp({v}) = {got}, want {want}");
        }
    }

    #[test]
    fn exp_works_on_types_with_few_integer_bits() {
        // I4F60 (range ±8) cannot hold the Taylor divisors or k·ln2 in T.
        for i in 0..=200 {
            let v = -2.0 + 4.0 * f64::from(i) / 200.0;
            let x = I4F60::from_num(v);
            let got: f64 = exp(x).to_num();
            let want = x.to_num::<f64>().exp();
            assert!(
                ((got - want) / want).abs() < 5e-15,
                "I4F60 exp({v}) = {got}, want {want}"
            );
        }
        assert_eq!(exp(I4F60::from_num(3)), I4F60::MAX);
    }

    #[test]
    fn exp_on_a_type_with_more_integer_than_fractional_bits() {
        // I24F8 divides for the estimate; k·ln2 carries k·0.0017 of quantisation.
        for i in 0..=100 {
            let v = -5.0 + 8.0 * f64::from(i) / 100.0;
            let x = I24F8::from_num(v);
            let got: f64 = exp(x).to_num();
            let want = x.to_num::<f64>().exp();
            assert!(
                (got - want).abs() < 0.01 * want.max(1.0),
                "I24F8 exp({v}) = {got}, want {want}"
            );
        }
        assert_eq!(exp(I24F8::MAX), I24F8::MAX);
        assert_eq!(exp(I24F8::MIN), I24F8::ZERO);
    }

    #[test]
    fn ln_normalisation_branches() -> TestResult {
        let ulp = I64F64::from_bits(1);
        let cases: [(I64F64, &str); 8] = [
            (I64F64::from_num(2), "exactly two"),
            (I64F64::from_num(4) + ulp, "lands on two after one shift"),
            (
                I64F64::from_num(4) + ulp + ulp,
                "lands on 1 + ulp after two shifts",
            ),
            (I64F64::from_num(5), "in (2, 4)"),
            (I64F64::from_num(1000), "large"),
            (I64F64::from_num(0.3), "below one half"),
            (I64F64::from_num(0.25), "exactly one quarter"),
            (I64F64::from_bits(1), "one ulp"),
        ];
        for (x, label) in cases {
            let got: f64 = ln(x)?.to_num();
            let want: f64 = x.to_num::<f64>().ln();
            assert!(
                (got - want).abs() < 1e-15 * want.abs().max(1.0),
                "ln({x}) [{label}] = {got}, want {want}"
            );
        }
        for k in 1..=60 {
            let want = f64::from(k) * core::f64::consts::LN_2;
            let got_log2: f64 = log2(I64F64::from_num(1u64 << k))?.to_num();
            assert!(
                (got_log2 - f64::from(k)).abs() < 1e-14,
                "log2(2^{k}) = {got_log2}"
            );
            let got_ln: f64 = ln(I64F64::from_num(1u64 << k))?.to_num();
            assert!(
                (got_ln - want).abs() < 1e-14,
                "ln(2^{k}) = {got_ln}, want {want}"
            );
        }
        assert!(ln(I16F16::from_num(4) + I16F16::from_bits(1)).is_ok());
        Ok(())
    }

    #[test]
    fn pow_special_cases_and_domain() -> TestResult {
        assert_eq!(pow(I16F16::ZERO, I16F16::ZERO)?, I16F16::ONE);
        assert_eq!(pow(I16F16::from_num(7), I16F16::ZERO)?, I16F16::ONE);
        assert_eq!(pow(I16F16::ZERO, I16F16::from_num(2))?, I16F16::ZERO);
        assert_eq!(pow(I16F16::ONE, I16F16::from_num(-9))?, I16F16::ONE);
        assert!(pow(I16F16::ZERO, -I16F16::ONE).is_err());
        assert!(pow(-I16F16::ONE, I16F16::from_num(0.5)).is_err());
        assert!(pow(-I16F16::from_num(2), I16F16::from_num(-1.5)).is_err());
        assert_eq!(
            pow(I16F16::from_num(200), I16F16::from_num(3))?,
            I16F16::MAX
        );
        assert_eq!(
            pow(-I16F16::from_num(200), I16F16::from_num(3))?,
            -I16F16::MAX
        );
        assert_eq!(pow(-I16F16::from_num(3), I16F16::ZERO)?, I16F16::ONE);
        assert_eq!(pow(-I16F16::ONE, I16F16::from_num(7))?, -I16F16::ONE);
        assert_eq!(pow(-I16F16::ONE, I16F16::from_num(-8))?, I16F16::ONE);
        for (b, e, want) in [
            (-2.0, 3.0, -8.0),
            (-2.0, 2.0, 4.0),
            (-2.0, -1.0, -0.5),
            (-0.5, 3.0, -0.125),
        ] {
            let got16: f64 = pow(I16F16::from_num(b), I16F16::from_num(e))?.to_num();
            assert!(
                (got16 - want).abs() < 0.05,
                "I16F16 {b}^{e} = {got16}, want {want}"
            );
            let got64: f64 = pow(I64F64::from_num(b), I64F64::from_num(e))?.to_num();
            assert!(
                (got64 - want).abs() < 1e-11,
                "I64F64 {b}^{e} = {got64}, want {want}"
            );
        }
        // Parity is read from the raw bits, so it survives exponents beyond i32.
        let huge = I64F64::from_num(1u64 << 40);
        assert_eq!(pow(-I64F64::ONE, huge)?, I64F64::ONE);
        assert_eq!(pow(-I64F64::ONE, huge + I64F64::ONE)?, -I64F64::ONE);
        let got: f64 = pow(I16F16::from_num(2), I16F16::from_num(10))?.to_num();
        assert!((got - 1024.0).abs() < 0.5, "I16F16 2^10 = {got}");
        Ok(())
    }

    #[test]
    fn pow_tracks_f64_at_i64f64() -> TestResult {
        let mut rng = Lcg(0x90);
        for _ in 0..2000 {
            let base = I64F64::from_num((rng.range(-20.0, 20.0)).exp2());
            let exponent = I64F64::from_num(rng.range(-3.0, 3.0));
            let (b, e): (f64, f64) = (base.to_num(), exponent.to_num());
            let got: f64 = pow(base, exponent)?.to_num();
            let want = b.powf(e);
            // Small results are quantised to a few raw ulps.
            let tol = 8.0f64.mul_add(2f64.powi(-64), 1e-14 * want);
            assert!(
                (got - want).abs() < tol,
                "pow({b}, {e}) = {got}, want {want}"
            );
        }
        let x = I64F64::from_num(2.5);
        let root: f64 = pow(x, I64F64::from_num(0.5))?.to_num();
        let want: f64 = sqrt(x)?.to_num();
        assert!((root - want).abs() < 1e-15, "2.5^0.5 = {root}, want {want}");
        let inv: f64 = pow(x, -I64F64::ONE)?.to_num();
        assert!((inv - 0.4).abs() < 1e-15, "2.5^-1 = {inv}");
        Ok(())
    }

    #[test]
    fn ln_tracks_f64_at_i64f64() -> TestResult {
        let mut rng = Lcg(0x1A);
        for _ in 0..3000 {
            let x = I64F64::from_num((rng.range(-60.0, 62.0)).exp2());
            let got: f64 = ln(x)?.to_num();
            let want = x.to_num::<f64>().ln();
            assert!(
                (got - want).abs() < 1e-15 * want.abs().max(1.0),
                "ln({x}) = {got}, want {want}"
            );
        }
        Ok(())
    }
}
