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
}
