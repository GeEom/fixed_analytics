#[cfg(test)]
mod tests {
    use crate::unit::support::{Lcg, TestResult};
    use fixed::types::{I16F16, I32F32, I64F64};
    use fixed_analytics::{
        atan_ratio, atanh_ratio, exp_moments, expm1, exprel, exprel2, log1p, log1p_ratio,
        log1p_residual, log1p_residual_ratio, mul_add, stumpff, tan_ratio, tanh_ratio,
    };

    #[test]
    fn values_at_zero() -> TestResult {
        let z = I32F32::ZERO;
        assert_eq!(exprel(z), I32F32::ONE);
        assert_eq!(exprel2(z), I32F32::from_num(0.5));
        assert_eq!(exp_moments(z)[..2], [I32F32::ONE, I32F32::from_num(0.5)]);
        assert_eq!(log1p_ratio(z)?, I32F32::ONE);
        assert_eq!(log1p_residual_ratio(z)?, I32F32::from_num(-0.5));
        assert_eq!(
            stumpff(z)[..3],
            [I32F32::ONE, I32F32::ONE, I32F32::from_num(0.5)]
        );
        assert_eq!(tan_ratio(z), (I32F32::ONE, I32F32::ONE));
        assert_eq!(atan_ratio(z), I32F32::ONE);
        assert_eq!(atanh_ratio(z)?, I32F32::ONE);
        assert_eq!(tanh_ratio(z), I32F32::ONE);
        Ok(())
    }

    #[test]
    fn tiny_arguments_keep_their_bits() -> TestResult {
        let tiny = I64F64::from_bits(1 << 20);
        assert_eq!(expm1(tiny), tiny);
        assert_eq!(log1p(tiny)?, tiny);
        // ln(1 + v) − v ≈ −v²/2, far below an ulp here: it rounds to 0 or −0.
        assert!(log1p_residual(tiny)?.abs() <= I64F64::from_bits(1));
        let v = I64F64::from_num(1e-6);
        let r: f64 = log1p_residual(v)?.to_num();
        assert!((r - (-5e-13 + 1e-18 / 3.0)).abs() < 1e-19);
        Ok(())
    }

    #[test]
    fn stumpff_identities_hold() {
        let mut rng = Lcg(0x5717);
        for _ in 0..500 {
            let z = I32F32::from_num(rng.range(-30.0, 30.0));
            let [c0, c1, c2, c3] = stumpff(z);
            // c₀ + z·c₂ = 1 and c₁ + z·c₃ = 1, each to a few ulps.
            let e0 = mul_add(z, c2, c0) - I32F32::ONE;
            let e1 = mul_add(z, c3, c1) - I32F32::ONE;
            let tol = I32F32::from_bits(40);
            assert!(e0.abs() <= tol && e1.abs() <= tol, "z = {z}: {e0}, {e1}");
        }
    }

    #[test]
    fn moments_and_relatives_agree() {
        let x = I32F32::from_num(2);
        let [e1, e2, e3, decay] = exp_moments(x);
        let ex = core::f64::consts::E.powi(-2);
        let close = |a: I32F32, b: f64| (a.to_num::<f64>() - b).abs() < 1e-9;
        assert!(close(e1, (1.0 - ex) / 2.0));
        assert!(close(e2, ex.mul_add(-3.0, 1.0) / 4.0));
        assert!(close(e3, ex.mul_add(-10.0, 2.0) / 8.0));
        assert!(close(decay, ex));
        assert!(close(
            exprel(x),
            core::f64::consts::E.mul_add(core::f64::consts::E, -1.0) / 2.0
        ));
    }

    #[test]
    fn domains_and_poles() {
        assert!(log1p(-I16F16::ONE).is_err());
        assert!(log1p_ratio(I16F16::from_num(-2)).is_err());
        assert!(atanh_ratio(I16F16::ONE).is_err());
        // fixed's FRAC_PI_2 truncates, so it sits just below the pole.
        assert_eq!(tan_ratio(I16F16::FRAC_PI_2).0, I16F16::MAX);
        assert_eq!(tan_ratio(I16F16::FRAC_PI_2 + I16F16::DELTA).0, I16F16::MIN);
        assert_eq!(exprel(I16F16::from_num(30)), I16F16::MAX);
        assert_eq!(exp_moments(I16F16::from_num(-40))[3], I16F16::MAX);
    }
}
