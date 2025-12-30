//! Tests for exponential and logarithmic functions

#[cfg(test)]
#[allow(clippy::unwrap_used)]
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
}
