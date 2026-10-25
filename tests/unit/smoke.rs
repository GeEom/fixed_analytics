#[cfg(test)]
mod tests {
    use crate::unit::support::TestResult;
    use core::hint::black_box;
    use fixed::types::I16F16;
    use fixed_analytics::{
        acos, acosh, acoth, asin, asinh, atan, atan2, atanh, cos, cosh, coth, exp, ln, log2, log10,
        sin, sin_cos, sinh, sinh_cosh, sqrt, tan, tanh,
    };

    #[test]
    fn smoke_test_trig() {
        let angle = I16F16::from_num(0.5);
        black_box(sin(angle));
        black_box(cos(angle));
        black_box(tan(angle));
        black_box(sin_cos(angle));
    }

    #[test]
    fn smoke_test_inverse_trig() -> TestResult {
        let x = I16F16::from_num(0.5);
        black_box(asin(x)?);
        black_box(acos(x)?);
        black_box(atan(x));
        black_box(atan2(x, I16F16::ONE));
        Ok(())
    }

    #[test]
    fn smoke_test_hyperbolic() -> TestResult {
        let x = I16F16::from_num(0.5);
        black_box(sinh(x));
        black_box(cosh(x));
        black_box(tanh(x));
        black_box(coth(x)?);
        black_box(sinh_cosh(x));
        Ok(())
    }

    #[test]
    fn smoke_test_inverse_hyperbolic() -> TestResult {
        let x = I16F16::from_num(0.5);
        black_box(asinh(x));
        black_box(atanh(x)?);

        let x_large = I16F16::from_num(1.5);
        black_box(acosh(x_large)?);
        black_box(acoth(x_large)?);
        Ok(())
    }

    #[test]
    fn smoke_test_exponential() -> TestResult {
        let x = I16F16::from_num(0.5);
        black_box(exp(x));
        black_box(ln(x)?);
        black_box(log2(x)?);
        black_box(log10(x)?);
        Ok(())
    }

    #[test]
    fn smoke_test_algebraic() -> TestResult {
        let x = I16F16::from_num(2.0);
        black_box(sqrt(x)?);
        Ok(())
    }
}

#[cfg(test)]
mod multi_type {
    use crate::unit::support::TestResult;
    use fixed::types::{I8F24, I32F32};
    use fixed_analytics::{acos, asin, atan, exp, ln, sin_cos, sinh_cosh, sqrt};

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
    fn inverse_trig_i32f32() -> TestResult {
        let x = I32F32::from_num(0.5);
        let asin_val = asin(x)?;
        let acos_val = acos(x)?;
        let atan_val = atan(x);

        let asin_f64: f64 = asin_val.to_num();
        assert!(asin_f64 > 0.5 && asin_f64 < 0.6, "asin(0.5) ≈ 0.524");

        let acos_f64: f64 = acos_val.to_num();
        assert!(acos_f64 > 1.0 && acos_f64 < 1.1, "acos(0.5) ≈ 1.047");

        let atan_f64: f64 = atan_val.to_num();
        assert!(atan_f64 > 0.45 && atan_f64 < 0.5, "atan(0.5) ≈ 0.464");
        Ok(())
    }

    #[test]
    fn hyperbolic_i32f32() {
        let x = I32F32::from_num(0.5);
        let (s, c) = sinh_cosh(x);
        let diff: f64 = (c * c - s * s).to_num();
        assert!(
            (diff - 1.0).abs() < 0.01,
            "cosh²-sinh² = {diff} (I32F32), expected ~1.0"
        );
    }

    #[test]
    fn exp_ln_i32f32() -> TestResult {
        let x = I32F32::from_num(2.0);
        let result = exp(ln(x)?);
        let result_f64: f64 = result.to_num();
        assert!(
            (result_f64 - 2.0).abs() < 0.3,
            "exp(ln(2)) = {result_f64} (I32F32), expected ~2.0"
        );
        Ok(())
    }

    #[test]
    fn sqrt_i32f32() -> TestResult {
        let x = I32F32::from_num(4.0);
        let result: f64 = sqrt(x)?.to_num();
        assert!(
            (result - 2.0).abs() < 0.001,
            "sqrt(4) = {result} (I32F32), expected 2.0"
        );
        Ok(())
    }

    #[test]
    fn trig_i8f24() {
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
    fn sqrt_i8f24() -> TestResult {
        let x = I8F24::from_num(2.0);
        let result: f32 = sqrt(x)?.to_num();
        assert!(
            (result - 1.414).abs() < 0.01,
            "sqrt(2) = {result} (I8F24), expected ~1.414"
        );
        Ok(())
    }

    #[test]
    fn ln_i8f24() -> TestResult {
        let x = I8F24::from_num(2.0);
        let result: f32 = ln(x)?.to_num();
        assert!(
            (result - 0.693).abs() < 0.2,
            "ln(2) = {result} (I8F24), expected ~0.693"
        );
        Ok(())
    }
}

#[cfg(test)]
mod i8f8 {
    use crate::unit::support::TestResult;
    use fixed::types::I8F8;
    use fixed_analytics::{sin_cos, sqrt};

    #[test]
    fn basic_trig() {
        let angle = I8F8::from_num(0.5);
        let (s, c) = sin_cos(angle);
        let sum_sq: f32 = (s * s + c * c).to_num();
        assert!(
            (sum_sq - 1.0).abs() < 0.2,
            "sin²+cos² = {sum_sq} (I8F8), expected ~1.0"
        );
    }

    #[test]
    fn basic_sqrt() -> TestResult {
        let x = I8F8::from_num(4.0);
        let result: f32 = sqrt(x)?.to_num();
        assert!(
            (result - 2.0).abs() < 0.1,
            "sqrt(4) = {result} (I8F8), expected ~2.0"
        );
        Ok(())
    }
}
