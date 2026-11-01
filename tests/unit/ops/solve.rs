#[cfg(test)]
mod tests {
    use crate::unit::support::{Lcg, TestResult};
    use fixed::types::extra::{U5, U13, U125};
    use fixed::types::{I16F16, I32F32, I64F64};
    use fixed::{FixedI8, FixedI16, FixedI128};
    use fixed_analytics::{
        CordicNumber, div_or_saturate, exp, kepler_universal, kepler_universal_near,
        kepler_universal_stumpff, root_in, root_in_newton, sin, sin_cos, stumpff,
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

    #[test]
    fn kepler_stumpff_matches_mpmath() -> TestResult {
        let t = |x: f64| I64F64::from_num(x);
        let close = |got: I64F64, want: f64| {
            (got.to_num::<f64>() - want).abs() < 4e-16 * want.abs().max(1.0)
        };
        // Hyperbolic, elliptic and the brief's case.
        let cases = [
            (
                0.243_580_071_287_150_6,
                -0.971_229_468_923_029_6,
                16.413_642_390_256_246,
            ),
            (
                -0.177_366_616_630_719_7,
                1.481_145_779_631_915_2,
                0.558_012_059_195_549_4,
            ),
            (1.9171, -2.4429, 0.58454),
        ];
        let want = [
            (
                2.877_069_856_311_604,
                [
                    8.547_788_020_915_303,
                    2.993_986_921_750_068_5,
                    0.938_851_770_480_668_6,
                    0.248_027_388_502_801_38,
                ],
            ),
            (
                0.606_608_030_904_304_9,
                [
                    0.739_643_309_593_761_4,
                    0.911_606_510_903_198_3,
                    0.477_699_322_693_433_6,
                    0.162_183_310_158_715_15,
                ],
            ),
            (
                0.394_580_713_472_777_6,
                [
                    1.196_276_890_168_411_3,
                    1.064_607_279_852_589_3,
                    0.516_049_986_861_274_4,
                    0.169_865_060_988_403_6,
                ],
            ),
        ];
        for ((s, a, tau), (chi, c)) in cases.into_iter().zip(want) {
            let (got, stumpffs) = kepler_universal_stumpff(t(1.0), t(s), t(a), t(tau))?;
            assert_eq!(got, kepler_universal(t(1.0), t(s), t(a), t(tau))?);
            assert!(close(got, chi));
            for (g, w) in stumpffs.into_iter().zip(c) {
                assert!(close(g, w));
            }
        }
        let (zero, at_zero) = kepler_universal_stumpff(t(1.0), t(0.5), t(1.0), t(0.0))?;
        assert_eq!((zero, at_zero), (t(0.0), stumpff(t(0.0))));
        Ok(())
    }

    #[test]
    fn kepler_stumpff_is_stumpff_at_its_anomaly() -> TestResult {
        let mut rng = Lcg(17);
        let ulp = I16F16::DELTA.to_num::<I64F64>();
        for _ in 0..400 {
            let r0 = rng.range(0.25, 4.0);
            let w = rng.range(0.0, 3.0);
            let alpha = I16F16::from_num((2.0 - w) / r0);
            let sigma = I16F16::from_num(rng.range(-0.9, 0.9) * (r0 * w).sqrt());
            let tau = I16F16::from_num(rng.range(-20.0, 20.0));
            let r0 = I16F16::from_num(r0);
            let Ok((chi, got)) = kepler_universal_stumpff(r0, sigma, alpha, tau) else {
                continue;
            };
            assert_eq!(chi, kepler_universal(r0, sigma, alpha, tau)?);
            // α·χ² is exact in I64F64.
            let chi = chi.to_num::<I64F64>();
            let want = stumpff(alpha.to_num::<I64F64>() * chi * chi);
            let (lo, hi) = (
                I16F16::MIN.to_num::<I64F64>(),
                I16F16::MAX.to_num::<I64F64>(),
            );
            for (c, exact) in got.into_iter().zip(want) {
                assert!((c.to_num::<I64F64>() - exact.clamp(lo, hi)).abs() <= ulp);
            }
        }
        Ok(())
    }

    #[test]
    fn kepler_near_ignores_its_start() -> TestResult {
        fn starts<T: CordicNumber + core::fmt::Debug>(
            r0: T,
            sigma: T,
            alpha: T,
            tau: T,
        ) -> TestResult {
            let want = kepler_universal_stumpff(r0, sigma, alpha, tau)?;
            let chi = want.0;
            let neg = |v: T| T::zero().saturating_sub(v);
            let mut ulp = T::one();
            for _ in 0..T::frac_bits() {
                ulp = ulp.div_int(2);
            }
            let seeds = [
                T::zero(),
                chi,
                chi.saturating_mul(T::two()),
                chi.div_int(2),
                chi.saturating_mul(T::from_num(1000)),
                chi.div_int(1000),
                neg(chi),
                ulp,
                T::max_value(),
                T::min_value(),
                T::one(),
                neg(T::one()),
            ];
            for near in seeds {
                assert_eq!(kepler_universal_near(r0, sigma, alpha, tau, near)?, want);
            }
            Ok(())
        }
        let orbits = [
            (1.0, 0.25, -0.5, 3.0),
            (1.0, -0.2, 0.75, 10.0),
            (1.0, 0.5, 0.0, -2.0),
            (0.5, 0.1, 3.5, 40.0),
            (1.0, 1.9171, -2.4429, 0.58454),
            (2.0, -0.3, -1.5, -25.0),
        ];
        for (r0, s, a, t) in orbits {
            starts(
                I64F64::from_num(r0),
                I64F64::from_num(s),
                I64F64::from_num(a),
                I64F64::from_num(t),
            )?;
            starts(
                I32F32::from_num(r0),
                I32F32::from_num(s),
                I32F32::from_num(a),
                I32F32::from_num(t),
            )?;
            starts(
                I16F16::from_num(r0),
                I16F16::from_num(s),
                I16F16::from_num(a),
                I16F16::from_num(t),
            )?;
        }
        Ok(())
    }
}
