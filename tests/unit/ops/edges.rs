#[cfg(test)]
mod tests {
    use crate::unit::support::TestResult;
    use fixed::types::{I16F16, I32F32, I64F64};
    use fixed_analytics::{
        Complex, Error, Modulus, carlson_rc, carlson_rd, carlson_rj, carlson_rj_roots, ellip_e,
        ellip_k, ellip_pi, ellip_pi_complete, jacobi_am, jacobi_inverse, jacobi_sn_cn_dn,
        kepler_universal, root_in, wrap_principal, wrap_radians,
    };

    #[test]
    fn wraps_that_round_onto_an_end_step_back_inside() {
        // Just below 2π, where the rounding first lands on 2π itself.
        let tau = core::f64::consts::TAU;
        let narrow = wrap_radians(-I16F16::DELTA).to_num::<f64>();
        assert!(narrow < tau && narrow > tau - 1e-3);
        let wide = wrap_radians(-I64F64::DELTA).to_num::<f64>();
        assert!(wide <= tau && wide > tau - 1e-9);
        // Exact wraps that round onto 2π, nearest 0 modulo 2π, and below −π.
        assert_eq!(wrap_radians(I16F16::from_bits(1_647_099)), I16F16::ZERO);
        let below = wrap_principal(I32F32::from_bits(1_524_713_260_611)).to_num::<f64>();
        let pi = core::f64::consts::PI;
        assert!(below > -pi && below < -pi + 1e-9);
        for x in [
            I32F32::from_num(pi),
            I32F32::from_num(-pi),
            I32F32::from_num(3.0 * pi),
        ] {
            let w = wrap_principal(x).to_num::<f64>();
            assert!(w > -pi && w <= pi + 1e-9, "{x}: {w}");
        }
    }

    #[test]
    fn solvers_take_their_edge_cases() -> TestResult {
        let line = |x: I32F32| x - I32F32::from_num(2);
        let two = I32F32::from_num(2);
        assert_eq!(root_in(line, I32F32::ZERO, two, 8)?, (two, two));
        assert_eq!(root_in(line, two, I32F32::from_num(5), 8)?, (two, two));
        let (a, b) = root_in(line, I32F32::from_num(9), I32F32::ONE, 100)?;
        assert!(a <= two && two <= b);
        let t = |x: f64| I64F64::from_num(x);
        assert!(kepler_universal(t(-1.0), t(0.0), t(1.0), t(1.0)).is_err());
        let far = kepler_universal(t(1.0), t(0.0), t(-1.0), t(1e6))?;
        assert!(far > t(10.0) && far < t(20.0));
        let back = kepler_universal(t(1.0), t(0.5), t(0.5), t(-3.0))?;
        assert!(back < t(0.0));
        Ok(())
    }

    #[test]
    fn elliptic_functions_take_their_edge_cases() -> TestResult {
        let m = |x: f64| {
            Modulus::from_parameter(I32F32::from_num(x)).ok_or(Error::domain("test", "modulus"))
        };
        let one = m(1.0)?;
        assert!(ellip_k(one).is_err());
        assert_eq!(ellip_e(one), I32F32::ONE);
        let (sn, cn, dn) = jacobi_sn_cn_dn(I32F32::from_num(0.5), one);
        assert!((sn.to_num::<f64>() - 0.5f64.tanh()).abs() < 1e-9 && cn == dn);
        assert!(jacobi_am(I32F32::from_num(5000), one) > I32F32::ONE);
        assert!(jacobi_inverse(I32F32::ONE, I32F32::ZERO, I32F32::ZERO, m(0.5)?).is_err());
        let half = m(0.5)?;
        assert_eq!(half.parameter(), I32F32::from_num(0.5));
        assert!((half.complement().to_num::<f64>() - 0.5f64.sqrt()).abs() < 1e-9);
        assert!(Modulus::from_parameter(I32F32::from_num(1.5)).is_none());
        assert!(Modulus::<I32F32>::from_complement(-I32F32::ONE).is_none());
        let kc = Modulus::from_complement(I32F32::from_num(0.6))
            .ok_or(Error::domain("test", "modulus"))?;
        assert!((kc.parameter().to_num::<f64>() - 0.64).abs() < 1e-9);
        assert!(carlson_rc(I32F32::ZERO, I32F32::ONE)? > I32F32::ONE);
        assert!(carlson_rc(-I32F32::ONE, I32F32::ONE).is_err());
        assert!(carlson_rd(I32F32::ZERO, I32F32::ZERO, I32F32::ONE).is_err());
        assert!(carlson_rj(I32F32::ONE, I32F32::ONE, I32F32::ONE, I32F32::ZERO).is_err());
        assert!(carlson_rj_roots(I32F32::ONE, I32F32::ONE, I32F32::ONE, I32F32::ZERO).is_err());
        assert!(carlson_rj_roots(I32F32::ZERO, I32F32::ZERO, I32F32::ONE, I32F32::ONE).is_err());
        let far = carlson_rj(
            I32F32::from_num(0.01),
            I32F32::from_num(0.5),
            I32F32::ONE,
            I32F32::from_num(400),
        )?;
        assert!(far > I32F32::ZERO && far < I32F32::from_num(0.05));
        let rj0 = carlson_rj(
            I32F32::ZERO,
            I32F32::from_num(0.5),
            I32F32::ONE,
            I32F32::from_num(400),
        )?;
        assert!(rj0 > I32F32::ZERO);
        assert!(ellip_pi_complete(I32F32::ONE, half).is_err());
        assert!(
            ellip_pi(
                I32F32::from_num(2),
                I32F32::ONE,
                I32F32::ZERO,
                I32F32::ONE,
                half
            )
            .is_err()
        );
        let wide = ellip_pi(
            I32F32::from_num(0.5),
            I32F32::from_num(0.6),
            I32F32::from_num(-0.8),
            I32F32::ONE,
            one,
        )?;
        assert_eq!(wide, I32F32::MAX);
        Ok(())
    }

    #[test]
    fn complex_edges() -> TestResult {
        let zero = Complex::new(I16F16::ZERO, I16F16::ZERO);
        assert_eq!(zero.arg(), I16F16::ZERO);
        assert_eq!(zero.pow(I16F16::ONE)?, zero);
        let up = Complex::new(I16F16::ZERO, I16F16::ONE);
        assert_eq!(up.arg(), fixed_analytics::atan2(I16F16::ONE, I16F16::ZERO));
        Ok(())
    }
}
