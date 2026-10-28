#[cfg(test)]
mod tests {
    use fixed::FixedU128;
    use fixed::types::extra::U64;
    use fixed::types::{I16F16, I32F32, I64F64};
    use fixed_analytics::{
        atan2_turns, div_by_tau, mul_by_tau, sin_cos, sin_cos_turns, wrap_principal, wrap_radians,
    };

    #[test]
    fn mul_by_tau_rounds_once() {
        // fixed's PI and FRAC_PI_2 truncate; these round to nearest.
        let pi = I32F32::from_num(core::f64::consts::PI);
        assert_eq!(mul_by_tau(I32F32::from_num(0.5)), pi);
        let quarter = fixed_analytics::atan2(I64F64::ONE, I64F64::ZERO);
        assert_eq!(mul_by_tau(I64F64::from_num(-0.25)), -quarter);
        assert_eq!(mul_by_tau(I16F16::from_num(10_000)), I16F16::MAX);
        let t = I64F64::from_num(0.123_456_789);
        assert_eq!(div_by_tau(mul_by_tau(t)), t);
    }

    #[test]
    fn quarter_turns_are_exact() {
        let quarter = I32F32::from_num(0.25);
        assert_eq!(sin_cos_turns(quarter), (I32F32::ONE, I32F32::ZERO));
        assert_eq!(sin_cos_turns(-quarter), (-I32F32::ONE, I32F32::ZERO));
        assert_eq!(
            sin_cos_turns(I32F32::from_num(1_000_000.5)),
            (I32F32::ZERO, -I32F32::ONE)
        );
    }

    #[test]
    fn turns_agree_with_radians() {
        let t = I64F64::from_num(0.1);
        let (s, c) = sin_cos_turns(t);
        let (sr, cr) = sin_cos(t * I64F64::TAU);
        assert!(
            (s - sr).abs() < I64F64::from_num(1e-17) && (c - cr).abs() < I64F64::from_num(1e-17)
        );
        assert_eq!(
            atan2_turns(I32F32::ONE, I32F32::ZERO),
            I32F32::from_num(0.25)
        );
        assert_eq!(
            atan2_turns(I32F32::ZERO, -I32F32::ONE),
            I32F32::from_num(0.5)
        );
        assert_eq!(
            atan2_turns(-I32F32::ONE, -I32F32::ONE),
            I32F32::from_num(-0.375)
        );
    }

    #[test]
    fn turns_from_an_unsigned_fraction() {
        // A U64F64 fraction of a turn is the same bits as an I64F64 one below ½.
        let frac = FixedU128::<U64>::from_num(0.125);
        let t = I64F64::from_bits(frac.to_bits().cast_signed());
        let (s, c) = sin_cos_turns(t);
        assert!((s - c).abs() <= I64F64::from_bits(1));
    }

    #[test]
    fn wraps_stay_in_range() {
        for k in -2000..2000 {
            let x = I16F16::from_num(k) / 7;
            let w = wrap_radians(x);
            assert!(
                w >= I16F16::ZERO && w.to_num::<f64>() < core::f64::consts::TAU,
                "{x} -> {w}"
            );
            let p = wrap_principal(x);
            let pf = p.to_num::<f64>();
            assert!(
                pf > -core::f64::consts::PI && pf <= core::f64::consts::PI,
                "{x} -> {p}"
            );
        }
        assert_eq!(div_by_tau(I32F32::TAU), I32F32::ONE);
    }
}
