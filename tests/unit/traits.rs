#[cfg(test)]
mod tests {
    use fixed::types::{I4F12, I4F60, I8F8, I8F24, I16F16, I20F12, I24F8, I32F32, I48F16, I64F64};
    use fixed_analytics::CordicNumber;

    #[test]
    fn basic_operations_i16f16() {
        let x = I16F16::from_num(2.5);
        assert_eq!(I16F16::zero(), I16F16::ZERO);
        assert_eq!(I16F16::one(), I16F16::ONE);
        assert!(I16F16::pi() > I16F16::from_num(314) / 100);
        assert!(I16F16::pi() < I16F16::from_num(3.15));
        assert!(!x.is_negative());
        assert!((-x).is_negative());
    }

    #[test]
    fn generic_impl_works_for_various_types() {
        let _: I8F24 = I8F24::pi();
        let _: I8F24 = I8F24::frac_pi_2();

        let _: I24F8 = I24F8::pi();
        let _: I24F8 = I24F8::frac_pi_2();

        let _: I4F12 = I4F12::pi();

        let _: I20F12 = I20F12::pi();

        let _: I48F16 = I48F16::pi();

        let _: I4F60 = I4F60::pi();
    }

    #[test]
    fn from_i1f63_across_types() {
        let half_bits: i64 = 0x4000_0000_0000_0000;

        let i8f8_half: f32 = I8F8::from_i1f63(half_bits).to_num();
        assert!((i8f8_half - 0.5).abs() < 0.01);

        let i16f16_half: f32 = I16F16::from_i1f63(half_bits).to_num();
        assert!((i16f16_half - 0.5).abs() < 0.0001);

        let i32f32_half: f64 = I32F32::from_i1f63(half_bits).to_num();
        assert!((i32f32_half - 0.5).abs() < 1e-9);

        let i64f64_half: f64 = I64F64::from_i1f63(half_bits).to_num();
        assert!((i64f64_half - 0.5).abs() < 1e-15);

        let i24f8_half: f32 = I24F8::from_i1f63(half_bits).to_num();
        assert!((i24f8_half - 0.5).abs() < 0.01);

        let i4f60_half: f64 = I4F60::from_i1f63(half_bits).to_num();
        assert!((i4f60_half - 0.5).abs() < 1e-15);
    }

    #[test]
    fn frac_bits_correct() {
        assert_eq!(I8F8::frac_bits(), 8);
        assert_eq!(I16F16::frac_bits(), 16);
        assert_eq!(I32F32::frac_bits(), 32);
        assert_eq!(I64F64::frac_bits(), 64);

        assert_eq!(I8F24::frac_bits(), 24);
        assert_eq!(I24F8::frac_bits(), 8);
        assert_eq!(I4F12::frac_bits(), 12);
        assert_eq!(I48F16::frac_bits(), 16);
    }

    #[test]
    fn div_overflow_saturates_to_min() {
        // Signs disagree, since zero counts as positive.
        let neg = I16F16::from_num(-1);
        let zero = I16F16::ZERO;
        assert_eq!(neg.div(zero), I16F16::MIN);
    }

    #[test]
    fn div_by_zero_positive_saturates_to_max() {
        let pos = I16F16::from_num(1);
        let zero = I16F16::ZERO;
        assert_eq!(pos.div(zero), I16F16::MAX);
    }

    #[test]
    fn frac_pi_4_values() {
        let pi_4_16: f32 = I16F16::frac_pi_4().to_num();
        assert!(
            (pi_4_16 - core::f32::consts::FRAC_PI_4).abs() < 0.001,
            "I16F16::frac_pi_4() = {pi_4_16}, expected ~0.7854"
        );

        let pi_4_32: f64 = I32F32::frac_pi_4().to_num();
        assert!(
            (pi_4_32 - core::f64::consts::FRAC_PI_4).abs() < 1e-9,
            "I32F32::frac_pi_4() = {pi_4_32}, expected ~0.7854"
        );

        let pi_4_64: f64 = I64F64::frac_pi_4().to_num();
        assert!(
            (pi_4_64 - core::f64::consts::FRAC_PI_4).abs() < 1e-15,
            "I64F64::frac_pi_4() = {pi_4_64}, expected ~0.7854"
        );
    }
}

#[cfg(test)]
mod fast_path_helpers {
    use fixed::types::{I4F4, I8F8, I16F16, I32F32, I64F64};
    use fixed_analytics::CordicNumber;

    #[test]
    fn div_int_matches_fixed_point_division_for_non_negative() {
        for raw in [0i64, 1, 7, 12, 65_536, 100_000, i64::from(i32::MAX)] {
            let x = I32F32::from_bits(raw << 20);
            for k in 1..=200u32 {
                assert_eq!(
                    x.div_int(k),
                    CordicNumber::div(x, I32F32::from_num(k)),
                    "{x} / {k}"
                );
            }
        }
    }

    #[test]
    fn div_int_truncates_toward_zero_and_saturates_on_zero() {
        let x = I16F16::from_num(-1.5);
        assert_eq!(x.div_int(4), I16F16::from_num(-0.375));
        assert_eq!(I16F16::from_bits(-7).div_int(2), I16F16::from_bits(-3));
        assert_eq!(I16F16::ONE.div_int(0), I16F16::MAX);
        assert_eq!((-I16F16::ONE).div_int(0), I16F16::MIN);
        // A divisor beyond the raw type's range exceeds |bits|.
        assert_eq!(I4F4::MAX.div_int(200), I4F4::ZERO);
        assert_eq!(I4F4::MIN.div_int(200), I4F4::ZERO);
    }

    #[test]
    fn mul_int_is_exact_and_saturates() {
        assert_eq!(
            I16F16::LN_2.mul_int(3),
            I16F16::LN_2 + I16F16::LN_2 + I16F16::LN_2
        );
        assert_eq!(I16F16::LN_2.mul_int(-2), -(I16F16::LN_2 + I16F16::LN_2));
        assert_eq!(I16F16::LN_2.mul_int(0), I16F16::ZERO);
        assert_eq!(I16F16::from_num(2).mul_int(100_000), I16F16::MAX);
        assert_eq!(I16F16::from_num(-2).mul_int(100_000), I16F16::MIN);
        assert_eq!(I8F8::ONE.mul_int(40_000), I8F8::MAX);
        assert_eq!(I8F8::ONE.mul_int(-40_000), I8F8::MIN);
        assert_eq!((-I8F8::ONE).mul_int(40_000), I8F8::MIN);
        assert_eq!((-I8F8::ONE).mul_int(-40_000), I8F8::MAX);
        assert_eq!(I8F8::ZERO.mul_int(40_000), I8F8::ZERO);
        assert_eq!(I4F4::ONE.mul_int(200), I4F4::MAX);
    }

    #[test]
    fn wrapping_ops_match_saturating_ops_in_range() {
        let a = I64F64::from_num(0.617);
        let b = I64F64::from_num(-0.1667);
        assert_eq!(a.wrapping_mul(b), a.saturating_mul(b));
        assert_eq!(a.wrapping_add(b), a.saturating_add(b));
        assert_eq!(a.wrapping_sub(b), a.saturating_sub(b));
        assert_eq!(I16F16::MAX.wrapping_add(I16F16::from_bits(1)), I16F16::MIN);
        assert_eq!(I16F16::MIN.wrapping_sub(I16F16::from_bits(1)), I16F16::MAX);
    }

    #[test]
    fn checked_int_log2_is_floor_of_log2() {
        assert_eq!(I16F16::ONE.checked_int_log2(), Some(0));
        assert_eq!(I16F16::from_num(0.5).checked_int_log2(), Some(-1));
        assert_eq!(I16F16::from_num(0.49).checked_int_log2(), Some(-2));
        assert_eq!(I16F16::from_num(6).checked_int_log2(), Some(2));
        assert_eq!(I16F16::from_num(8).checked_int_log2(), Some(3));
        assert_eq!(I16F16::from_bits(1).checked_int_log2(), Some(-16));
        assert_eq!(I16F16::ZERO.checked_int_log2(), None);
        assert_eq!(I16F16::from_num(-1).checked_int_log2(), None);
        assert_eq!(I64F64::MAX.checked_int_log2(), Some(62));
    }

    #[test]
    fn to_i32_and_round_saturate() {
        assert_eq!(I32F32::MAX.to_i32(), i32::MAX);
        assert_eq!(I32F32::MIN.to_i32(), i32::MIN);
        assert_eq!(I64F64::from_num(1e12).to_i32(), i32::MAX);
        assert_eq!(I64F64::from_num(-1e12).to_i32(), i32::MIN);
        assert_eq!(I16F16::from_num(-2.7).to_i32(), -3);
        assert_eq!(I16F16::from_num(2.7).to_i32(), 2);
        assert_eq!(CordicNumber::round(I16F16::MAX), I16F16::MAX);
        assert_eq!(CordicNumber::round(I16F16::MIN), I16F16::MIN);
        assert_eq!(
            CordicNumber::round(I16F16::from_num(2.5)),
            I16F16::from_num(3)
        );
        assert_eq!(
            CordicNumber::round(I16F16::from_num(-2.5)),
            I16F16::from_num(-3)
        );
    }

    #[test]
    fn sqrt_round_of_negative_is_root_of_magnitude() {
        assert_eq!(I16F16::from_num(-4).sqrt_round(), I16F16::from_num(2));
        assert_eq!(I64F64::from_num(-4).sqrt_round(), I64F64::from_num(2));
        assert_eq!(I16F16::ZERO.sqrt_round(), I16F16::ZERO);
    }
}
