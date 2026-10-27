#[cfg(test)]
mod tests {
    use crate::unit::support::TestResult;
    use fixed::types::{I8F0, I8F24, I10F118, I16F0, I16F16, I28F100, I32F0, I64F64, I128F0};
    use fixed_analytics::{Complex, erf, erfc, fresnel, gamma};

    fn z(re: f64, im: f64) -> Complex<I64F64> {
        Complex::new(I64F64::from_num(re), I64F64::from_num(im))
    }

    fn near(got: Complex<I64F64>, re: f64, im: f64) -> bool {
        let tol = 4e-16 * re.hypot(im).max(1.0);
        (got.re.to_num::<f64>() - re).abs() <= tol && (got.im.to_num::<f64>() - im).abs() <= tol
    }

    #[test]
    fn arithmetic_rounds_each_part_once() {
        let big = Complex::new(I16F16::from_num(200), I16F16::from_num(150));
        let q = big / big;
        assert_eq!(q, Complex::new(I16F16::ONE, I16F16::ZERO));
        let p = big * Complex::new(I16F16::from_num(0.25), I16F16::from_num(0.5));
        assert_eq!(
            p,
            Complex::new(I16F16::from_num(-25), I16F16::from_num(137.5))
        );
        let zero = Complex::new(I16F16::ZERO, I16F16::ZERO);
        assert_eq!(big / zero, Complex::new(I16F16::MAX, I16F16::MAX));
        assert_eq!(big.conj().im, I16F16::from_num(-150));
        assert_eq!(
            Complex::new(I16F16::from_num(3), I16F16::from_num(4)).norm(),
            I16F16::from_num(5)
        );
    }

    #[test]
    fn exponentials_match_mpmath() {
        let one = z(1.0, 0.0);
        assert!(near(one.exp(), core::f64::consts::E, 0.0));
        assert!(near(z(0.0, 0.0).exprel(), 1.0, 0.0));
        assert!(near(z(0.0, 0.0).exprel2(), 0.5, 0.0));
        let w = z(0.0, core::f64::consts::PI);
        assert!(near(w.exp(), -1.0, 0.0));
    }

    #[test]
    fn faddeeva_and_erf_match_mpmath() {
        assert!(near(
            z(0.5, 0.3).faddeeva(),
            0.614_851_539_146_991,
            0.303_124_349_647_351_1
        ));
        assert!(near(
            z(5.0, 0.2).faddeeva(),
            0.004_807_037_335_964_243,
            0.115_040_120_131_140_57
        ));
        assert!(near(
            z(-2.0, 1.0).faddeeva(),
            0.140_239_581_366_277_94,
            -0.222_213_440_179_899_1
        ));
        assert!(near(
            z(0.5, 0.3).erf(),
            0.561_565_188_524_213_2,
            0.267_605_864_957_603_6
        ));
        assert!(near(
            z(-2.0, 1.0).erf(),
            -1.003_606_342_725_651_8,
            -0.011_259_006_028_815_025
        ));
    }

    #[test]
    fn exponential_integral_matches_mpmath() -> TestResult {
        assert!(near(
            z(0.5, 0.3).e1()?,
            0.422_211_324_225_017_94,
            -0.305_371_136_174_266_65
        ));
        assert!(near(
            z(5.0, 0.2).e1()?,
            0.001_116_122_642_478_627,
            -0.000_266_867_640_371_248_9
        ));
        assert!(near(
            z(0.0, 6.0).e1()?,
            0.068_057_243_893_247_13,
            -0.146_108_775_514_390_08
        ));
        // On the negative axis, the value from above it.
        assert!(near(
            z(-3.0, 0.0).e1()?,
            -9.933_832_570_625_416,
            -core::f64::consts::PI
        ));
        assert!(z(0.0, 0.0).e1().is_err());
        Ok(())
    }

    #[test]
    fn incomplete_gamma_matches_mpmath() -> TestResult {
        let a = I64F64::from_num(2.5);
        assert!(near(
            z(0.5, 0.3).gamma_upper(a)?,
            1.299_139_215_932_174_5,
            -0.066_377_001_160_412_6
        ));
        assert!(near(
            z(-2.0, 1.0).gamma_upper(a)?,
            -11.043_544_259_962_623,
            4.945_042_286_802_596
        ));
        assert!(near(
            z(0.0, 6.0).gamma_upper(a)?,
            -10.838_494_867_793_583,
            10.170_960_665_022_385
        ));
        assert!(near(
            z(0.0, 0.0).gamma_upper(a)?,
            1.329_340_388_179_137,
            0.0
        ));
        assert!(z(1.0, 1.0).gamma_upper(I64F64::ZERO).is_err());
        Ok(())
    }

    #[test]
    fn wide_values_keep_exact_parts_and_phases() -> TestResult {
        let base: i128 = (1 << 100) + 12_345;
        let wide = |re: i128, im: i128| Complex::new(I128F0::from_num(re), I128F0::from_num(im));
        assert_eq!(
            (wide(base, base + 1) * wide(base + 2, base + 3)).re,
            I128F0::from_num(-2 * base - 3)
        );
        let far = wide(
            16_669_777_780_427_424 << 60,
            -(167_377_028_183_139_440 << 60),
        );
        assert_eq!(far.faddeeva(), wide(i128::MAX, i128::MAX));
        let axis = Complex::new(I8F0::ZERO, I8F0::from_num(-7)).erf();
        assert_eq!(axis, Complex::new(I8F0::ZERO, I8F0::MIN));
        let on_axis = Complex::new(I8F0::from_num(-50), I8F0::ZERO).gamma_upper(I8F0::ONE)?;
        assert_eq!(on_axis, Complex::new(I8F0::MAX, I8F0::ZERO));
        let far_left = Complex::new(I16F0::from_num(-12_570), I16F0::from_num(-2)).e1()?;
        assert_eq!(far_left, Complex::new(I16F0::MAX, I16F0::MIN));
        assert_eq!(erfc(I32F0::MIN), I32F0::from_num(2));
        let part = Complex::new(I8F24::from_bits(16), I8F24::from_num(4))
            .erf()
            .re;
        assert!((part.to_num::<f64>() - 9.562_398_898_839_027).abs() < 6e-8);
        let tiny = I28F100::from_bits(1 << 10);
        let near_axis =
            Complex::new(I28F100::from_num(-40), tiny).gamma_upper(I28F100::from_num(3))?;
        assert!((near_axis.im.to_num::<f64>() + 3.042_283_268_878_188_5e-7).abs() < 1e-22);
        let one = I10F118::ONE;
        let small = Complex::new(one + one, one).gamma_upper(I10F118::from_bits(1 << 8))?;
        assert!((small.re.to_bits() - 3_119_751_710_710_189_250_821_825_353_125_381).abs() <= 1);
        assert!((small.im.to_bits() + 14_775_364_147_287_538_117_270_599_311_667_794).abs() <= 1);
        Ok(())
    }

    #[test]
    fn real_special_functions() -> TestResult {
        let f = |x: f64| I64F64::from_num(x);
        assert!((erf(f(1.0)).to_num::<f64>() - 0.842_700_792_949_714_9).abs() < 3e-16);
        assert!((erfc(f(2.0)).to_num::<f64>() - 0.004_677_734_981_047_266).abs() < 1e-17);
        assert_eq!(erfc(f(-30.0)), f(2.0));
        let (c, s) = fresnel(f(1.0));
        assert!((c.to_num::<f64>() - 0.779_893_400_376_822_8).abs() < 3e-16);
        assert!((s.to_num::<f64>() - 0.438_259_147_390_354_77).abs() < 3e-16);
        assert!((gamma(f(4.5))?.to_num::<f64>() - 11.631_728_396_567_45).abs() < 1e-14);
        assert!((gamma(f(0.5))?.to_num::<f64>() - 1.772_453_850_905_516).abs() < 5e-16);
        assert!(gamma(f(0.0)).is_err());
        Ok(())
    }
}
