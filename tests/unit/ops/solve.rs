#[cfg(test)]
mod tests {
    use crate::unit::support::TestResult;
    use fixed::types::extra::{U5, U13, U125};
    use fixed::types::{I16F16, I32F32, I64F64};
    use fixed::{FixedI8, FixedI16, FixedI128};
    use fixed_analytics::{
        CordicNumber, div_or_saturate, exp, kepler_universal, root_in, root_in_newton, sin, sin_cos,
    };

    #[test]
    fn root_in_closes_to_adjacent_values() -> TestResult {
        let (a, b) = root_in(sin, I32F32::from_num(3), I32F32::from_num(3.5), 64)?;
        assert_eq!((a, b), (I32F32::from_num(core::f64::consts::PI), a));
        assert_eq!(sin(a), I32F32::ZERO);
        let cube = |x: I64F64| x * x * x - I64F64::from_num(2);
        let (low, high) = root_in(cube, I64F64::from_num(2), I64F64::ZERO, 200)?;
        assert_eq!(high - low, I64F64::DELTA);
        let root = 2f64.cbrt();
        assert!((low.to_num::<f64>() - root).abs() < 1e-15);
        Ok(())
    }

    #[test]
    fn root_in_keeps_to_its_budget() -> TestResult {
        let mut calls = 0;
        let f = |x: I16F16| {
            calls += 1;
            x * I16F16::from_num(3) - I16F16::ONE
        };
        let (a, b) = root_in(f, I16F16::from_num(-100), I16F16::from_num(100), 5)?;
        assert_eq!(calls, 7);
        let third = 1.0 / 3.0;
        assert!(a.to_num::<f64>() < third && third < b.to_num::<f64>());
        let line = |x: I16F16| x - I16F16::from_num(0.5);
        assert_eq!(
            root_in(line, I16F16::ZERO, I16F16::ONE, 40)?,
            (I16F16::from_num(0.5), I16F16::from_num(0.5))
        );
        assert!(root_in(line, I16F16::ONE, I16F16::from_num(2), 40).is_err());
        Ok(())
    }

    #[test]
    fn root_in_probes_inside_its_bracket() -> TestResult {
        // Regula falsi lands within an ulp of an end here, where f steps by single ulps.
        let (e, m) = (I16F16::from_num(0.9), I16F16::from_num(1.05));
        let kepler = |x: I16F16| x - e * sin(x) - m;
        let (a, b) = root_in(kepler, I16F16::ZERO, I16F16::from_num(3), 40)?;
        assert_eq!((a.to_bits(), b.to_bits()), (124_602, 124_603));
        Ok(())
    }

    /// Whether both solvers close on one bracket, or on exact zeros.
    fn agrees<T: CordicNumber>(
        f: impl Fn(T) -> T,
        step: impl Fn(T, T) -> T,
        guess: T,
        (lo, hi): (T, T),
        budget: u32,
    ) -> Result<bool, fixed_analytics::Error> {
        let itp = root_in(&f, lo, hi, 400)?;
        let newton = root_in_newton(&f, step, guess, lo, hi, budget)?;
        let zero = |r: (T, T)| r == (r.0, r.0) && f(r.0) == T::zero();
        Ok(newton == itp || (zero(newton) && zero(itp)))
    }

    fn monotone_roots<T: CordicNumber>(name: &str) -> TestResult {
        let n = |v: f64| T::from_num(v);
        for c in [0.2, 0.6, 1.0, 1.25].map(n) {
            let cube = |x: T| x * x * x - c;
            let newton = |x: T, y: T| div_or_saturate(y, (x + x + x) * x);
            assert!(
                agrees(cube, newton, n(1.0), (n(0.0), n(1.1)), 200)?,
                "{name}: cube"
            );
        }
        for c in [1.5, 2.0, 3.5].map(n) {
            let rise = |x: T| exp(x) - c;
            // f/f′ = 1 − c·e^(−x).
            let newton = |x: T, _| T::one() - c * exp(-x);
            assert!(agrees(rise, newton, n(0.5), (n(0.0), n(1.3)), 200)?);
        }
        let e = n(0.9);
        for m in [0.05, 0.3, 1.0, 2.5].map(n) {
            let kepler = |x: T| x - e * sin(x) - m;
            let newton = |x: T, y: T| div_or_saturate(y, T::one() - e * sin_cos(x).1);
            let halley = |x: T, y: T| {
                let (s, c) = sin_cos(x);
                let slope = T::one() - e * c;
                let lean = div_or_saturate(y * e * s, slope + slope);
                div_or_saturate(y, slope - lean)
            };
            let span = (n(0.0), n(3.0));
            assert!(agrees(kepler, newton, m, span, 200)?, "{name}: kepler");
            assert!(agrees(kepler, halley, m, span, 200)?, "{name}: halley");
        }
        Ok(())
    }

    #[test]
    fn root_in_newton_closes_where_root_in_does() -> TestResult {
        monotone_roots::<FixedI8<U5>>("I3F5")?;
        monotone_roots::<FixedI16<U13>>("I3F13")?;
        monotone_roots::<I16F16>("I16F16")?;
        monotone_roots::<I32F32>("I32F32")?;
        monotone_roots::<I64F64>("I64F64")?;
        monotone_roots::<FixedI128<U125>>("I3F125")
    }

    #[test]
    fn root_in_newton_bisects_past_bad_steps() -> TestResult {
        let e = I64F64::from_num(0.9);
        let m = I64F64::from_num(0.3);
        let kepler = |x: I64F64| x - e * sin(x) - m;
        let span = (I64F64::ZERO, I64F64::from_num(3));
        let tiny = |_, y: I64F64| y >> 40;
        let huge = |_, _| I64F64::MAX;
        let wrong = |x: I64F64, _| x;
        let none = |_, _| I64F64::ZERO;
        for budget in [140, 200] {
            assert!(agrees(kepler, tiny, m, span, budget)?);
            assert!(agrees(kepler, huge, m, span, budget)?);
            assert!(agrees(kepler, wrong, m, span, budget)?);
            assert!(agrees(kepler, none, m, span, budget)?);
            assert!(agrees(kepler, none, I64F64::MIN, span, budget)?);
        }
        Ok(())
    }

    #[test]
    fn root_in_newton_keeps_to_its_budget() -> TestResult {
        let mut calls = 0;
        let f = |x: I16F16| {
            calls += 1;
            x * I16F16::from_num(3) - I16F16::ONE
        };
        let lo = I16F16::from_num(-100);
        let hi = I16F16::from_num(100);
        let (a, b) = root_in_newton(f, |_, _| I16F16::ZERO, lo, lo, hi, 5)?;
        assert_eq!(calls, 7);
        let third = 1.0 / 3.0;
        assert!(a.to_num::<f64>() < third && third < b.to_num::<f64>());
        let line = |x: I16F16| x - I16F16::from_num(0.5);
        let same = |_, y| y;
        assert_eq!(root_in_newton(line, same, lo, lo, hi, 0)?, (lo, hi));
        assert_eq!(
            root_in_newton(line, same, lo, I16F16::ZERO, I16F16::ONE, 40)?,
            (I16F16::from_num(0.5), I16F16::from_num(0.5))
        );
        Ok(())
    }

    #[test]
    fn root_in_newton_takes_its_edge_cases() -> TestResult {
        let two = I32F32::from_num(2);
        let line = |x: I32F32| x - two;
        let same = |_, y| y;
        assert_eq!(
            root_in_newton(line, same, two, I32F32::ZERO, two, 8)?,
            (two, two)
        );
        assert_eq!(
            root_in_newton(line, same, two, two, I32F32::from_num(5), 8)?,
            (two, two)
        );
        let swapped = root_in_newton(line, same, two, I32F32::from_num(9), I32F32::ONE, 8)?;
        assert_eq!(swapped, (two, two));
        let falling = |x: I32F32| two - x;
        let wide = I32F32::from_num(1000);
        let outside = root_in_newton(falling, same, wide, I32F32::ZERO, wide, 8)?;
        assert_eq!(outside, (two, two));
        let step = |x: I32F32, _| x;
        assert!(root_in_newton(line, step, two, I32F32::ONE, I32F32::ZERO, 40).is_err());
        assert!(root_in_newton(line, step, two, I32F32::from_num(3), wide, 40).is_err());
        Ok(())
    }

    #[test]
    fn kepler_matches_mpmath() -> TestResult {
        let t = |x: f64| I64F64::from_num(x);
        let close = |got: I64F64, want: f64| {
            (got.to_num::<f64>() - want).abs() < 4e-16 * want.abs().max(1.0)
        };
        // A circle, whose anomaly is its time, a parabola, a hyperbola and an ellipse.
        assert_eq!(kepler_universal(t(1.0), t(0.0), t(1.0), t(2.5))?, t(2.5));
        assert!(close(
            kepler_universal(t(1.0), t(0.5), t(0.0), t(2.0))?,
            1.263_997_823_528_725_5
        ));
        assert!(close(
            kepler_universal(t(1.0), t(0.25), t(-0.5), t(3.0))?,
            1.586_900_711_280_019
        ));
        assert!(close(
            kepler_universal(t(1.0), t(-0.2), t(0.75), t(10.0))?,
            7.591_239_334_586_656
        ));
        assert_eq!(kepler_universal(t(1.0), t(0.5), t(0.0), t(0.0))?, t(0.0));
        Ok(())
    }

    #[test]
    fn kepler_needs_a_real_orbit() {
        let t = |x: f64| I32F32::from_num(x);
        assert!(kepler_universal(t(0.0), t(0.0), t(1.0), t(1.0)).is_err());
        assert!(kepler_universal(t(1.0), t(2.0), t(1.0), t(1.0)).is_err());
        assert!(kepler_universal(t(1.0), t(1.0), t(1.0), t(1.0)).is_ok());
    }
}
