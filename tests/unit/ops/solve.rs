#[cfg(test)]
mod tests {
    use crate::unit::support::TestResult;
    use fixed::types::{I16F16, I32F32, I64F64};
    use fixed_analytics::{kepler_universal, root_in, sin};

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
