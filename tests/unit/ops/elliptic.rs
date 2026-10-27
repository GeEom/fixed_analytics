#[cfg(test)]
mod tests {
    use crate::unit::support::{Lcg, TestResult};
    use fixed::types::{I32F32, I64F64};
    use fixed_analytics::{
        Error, Modulus, carlson_rc, carlson_rd, carlson_rf, carlson_rf_roots, carlson_rj, ellip_e,
        ellip_k, ellip_pi, ellip_pi_complete, jacobi_am, jacobi_inverse, jacobi_sn_cn_dn,
    };

    fn half() -> Result<Modulus<I64F64>, Error> {
        Modulus::from_parameter(I64F64::from_num(0.5)).ok_or(Error::domain("test", "m"))
    }

    fn near(got: I64F64, want: f64, tol: f64) -> bool {
        (got.to_num::<f64>() - want).abs() <= tol
    }

    #[test]
    fn complete_integrals_match_mpmath() -> TestResult {
        let m = half()?;
        assert!(near(ellip_k(m)?, 1.854_074_677_301_372, 1e-15));
        assert!(near(ellip_e(m), 1.350_643_881_047_675_5, 1e-15));
        assert!(near(
            ellip_pi_complete(I64F64::from_num(-0.5), m)?,
            1.487_846_992_668_798_4,
            1e-15
        ));
        assert!(near(
            ellip_pi_complete(I64F64::from_num(-1000), m)?,
            0.050_138_952_552_005_41,
            1e-15
        ));
        let one = Modulus::from_complement(I64F64::ZERO).ok_or(Error::domain("test", "kc"))?;
        assert!(ellip_k(one).is_err());
        assert_eq!(ellip_e(one), I64F64::ONE);
        Ok(())
    }

    #[test]
    fn jacobi_functions_match_mpmath() -> TestResult {
        let (sn, cn, dn) = jacobi_sn_cn_dn(I64F64::from_num(0.5), half()?);
        assert!(near(sn, 0.470_750_473_655_657_3, 1e-15));
        assert!(near(cn, 0.882_266_394_890_440_3, 1e-15));
        assert!(near(dn, 0.942_972_425_777_385_7, 1e-15));
        let m9 =
            Modulus::from_parameter(I64F64::from_num(0.9)).ok_or(Error::domain("test", "m"))?;
        let far = jacobi_sn_cn_dn(I64F64::from_num(3), m9);
        assert!(near(far.0, 0.990_630_599_937_832_6, 1e-15));
        assert!(near(far.1, -0.136_568_717_013_853_34, 1e-15));
        assert!(near(far.2, 0.341_739_539_737_691_07, 1e-15));
        Ok(())
    }

    #[test]
    fn jacobi_identities_and_inverse() -> TestResult {
        let mut rng = Lcg(0xE11);
        for _ in 0..300 {
            let kc = I32F32::from_num(rng.range(0.0, 1.0));
            let m = Modulus::from_complement(kc).ok_or(Error::domain("test", "kc"))?;
            let k: f64 = ellip_k(m)?.to_num();
            let u = I32F32::from_num(rng.range(-1.9, 1.9) * k);
            let (sn, cn, dn) = jacobi_sn_cn_dn(u, m);
            let (sf, cf, df): (f64, f64, f64) = (sn.to_num(), cn.to_num(), dn.to_num());
            let kf = kc.to_num::<f64>();
            let mm = kf.mul_add(-kf, 1.0);
            assert!(sf.mul_add(sf, cf * cf) - 1.0 < 1e-8);
            assert!((mm * sf).mul_add(sf, df * df) - 1.0 < 1e-8);
            let back: f64 = jacobi_inverse(sn, cn, dn, m)?.to_num();
            assert!(
                (back - u.to_num::<f64>()).abs() < 1e-6 * k,
                "u = {u}, back = {back}"
            );
            let am = jacobi_am(u, m).to_num::<f64>();
            assert!((am.sin() - sf).abs() < 1e-8 && (am.cos() - cf).abs() < 1e-8);
        }
        Ok(())
    }

    #[test]
    fn separatrix_closed_forms() -> TestResult {
        let one = Modulus::from_complement(I64F64::ZERO).ok_or(Error::domain("test", "kc"))?;
        let (sn, cn, dn) = jacobi_sn_cn_dn(I64F64::from_num(2), one);
        assert!(near(sn, 2f64.tanh(), 1e-15));
        assert!(near(cn, 1.0 / 2f64.cosh(), 1e-15) && cn == dn);
        assert_eq!(
            jacobi_sn_cn_dn(I64F64::MAX, one),
            (I64F64::ONE, I64F64::ZERO, I64F64::ZERO)
        );
        Ok(())
    }

    #[test]
    fn carlson_integrals_match_mpmath() -> TestResult {
        let f = I64F64::from_num::<f64>;
        assert!(near(
            carlson_rf(f(1.0), f(2.0), f(3.0))?,
            0.726_945_935_468_908_2,
            1e-15
        ));
        assert!(near(
            carlson_rj(f(2.0), f(3.0), f(4.0), f(5.0))?,
            0.142_975_796_671_567_54,
            1e-15
        ));
        let big_p = carlson_rj(f(0.001), f(0.002), f(0.004), f(1000.0))?;
        assert!(near(big_p, 0.064_844_213_471_528_12, 1e-15));
        assert!(near(
            carlson_rc(f(2.25), f(2.0))?,
            core::f64::consts::LN_2,
            1e-15
        ));
        assert!(near(
            carlson_rd(f(0.0), f(2.0), f(1.0))?,
            1.797_210_352_103_388_3,
            1e-15
        ));
        // R_F(t², t², t²) = 1/t, though t² = 2^-108 underflows I64F64.
        let tiny = I64F64::from_bits(1 << 10);
        assert_eq!(
            carlson_rf_roots(tiny, tiny, tiny)?,
            I64F64::from_num(1u64 << 54)
        );
        assert!(carlson_rf(f(0.0), f(0.0), f(1.0)).is_err());
        assert!(carlson_rj(f(1.0), f(1.0), f(1.0), f(0.0)).is_err());
        Ok(())
    }

    #[test]
    fn incomplete_third_kind() -> TestResult {
        let (sn, cn, dn) = (
            I64F64::from_num(0.644_217_687_237_691),
            I64F64::from_num(0.764_842_187_284_488_4),
            I64F64::from_num(0.890_220_077_129_841_1),
        );
        let p = ellip_pi(I64F64::from_num(-0.5), sn, cn, dn, half()?)?;
        assert!(near(p, 0.679_470_265_670_900_9, 1e-14));
        assert!(
            ellip_pi(
                I64F64::from_num(2),
                I64F64::ONE,
                I64F64::ZERO,
                I64F64::from_num(0.7),
                half()?
            )
            .is_err()
        );
        Ok(())
    }
}
