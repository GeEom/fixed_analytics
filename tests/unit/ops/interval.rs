#[cfg(test)]
mod tests {
    use crate::unit::support::Lcg;
    use fixed::types::{I16F16, I32F32};
    use fixed_analytics::{
        atan_range, atan_ratio, atan_ratio_range, atan2, atan2_range, atanh_ratio,
        atanh_ratio_range, enclose, exp, exp_moments, exp_moments_range, exp_range, exprel,
        exprel_range, exprel2, exprel2_range, hypot, hypot_range, log1p, log1p_range, log1p_ratio,
        log1p_ratio_range, log1p_residual_ratio, log1p_residual_ratio_range, sin, sin_cos_range,
        sin_cos_range_turns, sqrt_or_zero, sqrt_range, stumpff, stumpff_range, tan_ratio,
        tan_ratio_range, tanh_range, tanh_ratio, tanh_ratio_range,
    };

    fn across(lo: I32F32, hi: I32F32, n: i64) -> impl Iterator<Item = I32F32> {
        let span = i128::from(hi.to_bits()) - i128::from(lo.to_bits());
        (0..=n).map(move |k| {
            let step = span * i128::from(k) / i128::from(n);
            I32F32::from_bits(i64::try_from(i128::from(lo.to_bits()) + step).unwrap_or(0))
        })
    }

    fn draw(rng: &mut Lcg, options: &[f64]) -> f64 {
        let at = usize::try_from(rng.next_u64() % options.len() as u64).unwrap_or(0);
        options.get(at).copied().unwrap_or(0.0)
    }

    fn holds(range: (I32F32, I32F32), v: I32F32) -> bool {
        range.0 <= v && v <= range.1
    }

    #[test]
    fn enclose_spans_an_ulp_each_way() {
        let y = exp(I32F32::ONE);
        let (lo, hi) = enclose(y);
        let e = core::f64::consts::E;
        assert!(lo.to_num::<f64>() < e && e < hi.to_num::<f64>());
        assert_eq!(enclose(I16F16::MAX).1, I16F16::MAX);
    }

    #[test]
    fn sin_cos_range_finds_interior_extremes() {
        let ((smin, smax), (cmin, cmax)) = sin_cos_range(I32F32::ZERO, I32F32::PI);
        assert_eq!(smax, I32F32::ONE);
        assert!(smin <= I32F32::ZERO);
        assert_eq!((cmin, cmax), (-I32F32::ONE, I32F32::ONE));
        let wide = sin_cos_range(I32F32::from_num(-10), I32F32::from_num(10));
        assert_eq!(
            wide,
            ((-I32F32::ONE, I32F32::ONE), (-I32F32::ONE, I32F32::ONE))
        );
        let (s, _) = sin_cos_range(I32F32::from_num(0.25), I32F32::from_num(0.5));
        assert!(s.0 <= sin(I32F32::from_num(0.25)) && sin(I32F32::from_num(0.5)) <= s.1);
        assert!(s.1 < I32F32::ONE);
    }

    #[test]
    fn turn_ranges_hold_every_sampled_value() {
        let quarter = I32F32::from_num(0.25);
        let below = -I32F32::DELTA;
        assert_eq!(
            sin_cos_range_turns(I32F32::ZERO, quarter),
            ((below, I32F32::ONE), (below, I32F32::ONE))
        );
        assert_eq!(
            sin_cos_range_turns(I32F32::from_num(-3), I32F32::from_num(-2)),
            ((-I32F32::ONE, I32F32::ONE), (-I32F32::ONE, I32F32::ONE))
        );
        let mut rng = Lcg(0x7E57);
        for _ in 0..2000 {
            let a = rng.range(-8.0, 8.0);
            let b = a + rng.range(0.0, 0.6);
            let (lo, hi) = (I32F32::from_num(a), I32F32::from_num(b));
            let ((smin, smax), (cmin, cmax)) = sin_cos_range_turns(hi, lo);
            for k in 0..=16 {
                let t = lo.to_num::<f64>().mul_add(
                    1.0 - f64::from(k) / 16.0,
                    hi.to_num::<f64>() * f64::from(k) / 16.0,
                );
                let x = t * core::f64::consts::TAU;
                let inside = |v: f64, l: I32F32, h: I32F32| {
                    l.to_num::<f64>() - 1e-12 <= v && v <= h.to_num::<f64>() + 1e-12
                };
                assert!(
                    inside(x.sin(), smin, smax),
                    "sin at {t} turns outside [{smin}, {smax}]"
                );
                assert!(
                    inside(x.cos(), cmin, cmax),
                    "cos at {t} turns outside [{cmin}, {cmax}]"
                );
            }
        }
    }

    #[test]
    fn ranges_hold_every_sampled_value() {
        let mut rng = Lcg(0x5EED);
        for _ in 0..2000 {
            let a = rng.range(-50.0, 50.0);
            let b = a + rng.range(0.0, 3.0);
            let (lo, hi) = (I32F32::from_num(a), I32F32::from_num(b));
            let ((smin, smax), (cmin, cmax)) = sin_cos_range(lo, hi);
            let (lo4, hi4) = (lo / 4, hi / 4);
            let (emin, emax) = exp_range(lo4, hi4);
            let (amin, amax) = atan_range(lo, hi);
            let (tmin, tmax) = tanh_range(lo, hi);
            for k in 0..=16 {
                let x = (lo.to_num::<f64>()).mul_add(
                    1.0 - f64::from(k) / 16.0,
                    hi.to_num::<f64>() * f64::from(k) / 16.0,
                );
                let inside = |v: f64, l: I32F32, h: I32F32| {
                    l.to_num::<f64>() - 1e-12 <= v && v <= h.to_num::<f64>() + 1e-12
                };
                assert!(
                    inside(x.sin(), smin, smax),
                    "sin({x}) outside [{smin}, {smax}]"
                );
                assert!(
                    inside(x.cos(), cmin, cmax),
                    "cos({x}) outside [{cmin}, {cmax}]"
                );
                let x4 = (lo4.to_num::<f64>()).mul_add(
                    1.0 - f64::from(k) / 16.0,
                    hi4.to_num::<f64>() * f64::from(k) / 16.0,
                );
                assert!(inside(x4.exp(), emin, emax), "exp outside at {x4}");
                assert!(inside(x.atan(), amin, amax), "atan outside at {x}");
                assert!(inside(x.tanh(), tmin, tmax), "tanh outside at {x}");
            }
        }
    }

    #[test]
    fn kernel_ranges_hold_every_value_between() -> crate::unit::support::TestResult {
        let mut rng = Lcg(0x4A7E);
        for _ in 0..1500 {
            let a = rng.range(-30.0, 30.0);
            let width = draw(&mut rng, &[1e-6, 1e-3, 0.1, 2.0, 20.0]);
            let (lo, hi) = (
                I32F32::from_num(a),
                I32F32::from_num(a + rng.range(0.0, width)),
            );
            let (ex, ex2, th, at, tn, sq) = (
                exprel_range(hi, lo),
                exprel2_range(lo, hi),
                tanh_ratio_range(lo, hi),
                atan_ratio_range(lo, hi),
                tan_ratio_range(lo, hi),
                sqrt_range(lo, hi),
            );
            let moments = exp_moments_range(lo, hi);
            for x in across(lo, hi, 64) {
                assert!(
                    holds(ex, exprel(x)) && holds(ex2, exprel2(x)),
                    "exprel at {x}"
                );
                assert!(
                    holds(th, tanh_ratio(x)) && holds(at, atan_ratio(x)),
                    "ratio at {x}"
                );
                assert!(
                    holds(tn, tan_ratio(x).0),
                    "tan_ratio at {x} in [{lo}, {hi}]"
                );
                assert!(holds(sq, sqrt_or_zero(x)), "sqrt at {x}");
                for (range, v) in moments.iter().zip(exp_moments(x)) {
                    assert!(holds(*range, v), "moment at {x}");
                }
            }
            let (ul, uh) = (lo / 32, hi / 32);
            if ul > -I32F32::ONE {
                let (l1, lr, lq) = (
                    log1p_range(ul, uh)?,
                    log1p_ratio_range(ul, uh)?,
                    log1p_residual_ratio_range(ul, uh)?,
                );
                for x in across(ul, uh, 64) {
                    assert!(
                        holds(l1, log1p(x)?) && holds(lr, log1p_ratio(x)?),
                        "log1p at {x}"
                    );
                    assert!(holds(lq, log1p_residual_ratio(x)?), "residual at {x}");
                }
            }
            if ul > -I32F32::ONE && uh < I32F32::ONE {
                let range = atanh_ratio_range(ul, uh)?;
                assert!(across(ul, uh, 64).all(|x| atanh_ratio(x).is_ok_and(|v| holds(range, v))));
            }
        }
        assert!(log1p_range(-I32F32::ONE, I32F32::ZERO).is_err());
        assert!(atanh_ratio_range(I32F32::ZERO, I32F32::ONE).is_err());
        Ok(())
    }

    #[test]
    fn stumpff_ranges_hold_through_their_extremes() {
        let mut rng = Lcg(0x57F);
        // Extremes of c₀ at (jπ)², of c₁ near 20.19 and 59.68, of c₂ at 4π² and near 80.8.
        let hot = [9.8696, 20.19, 39.478, 59.68, 80.76, 88.83, 157.9, 4.0, 0.0];
        for i in 0..3000 {
            let centre = if i % 3 == 0 {
                rng.range(-40.0, 400.0)
            } else {
                draw(&mut rng, &hot)
            };
            let width = draw(&mut rng, &[1e-7, 1e-4, 0.01, 1.0, 30.0]);
            let lo = I32F32::from_num(centre - rng.range(0.0, width));
            let hi = I32F32::from_num(centre + rng.range(0.0, width));
            let ranges = stumpff_range(lo, hi);
            for x in across(lo, hi, 200) {
                for (k, (range, v)) in ranges.iter().zip(stumpff(x)).enumerate() {
                    assert!(
                        holds(*range, v),
                        "c{k}({x}) = {v} outside {range:?} over [{lo}, {hi}]"
                    );
                }
            }
        }
        let (c1lo, c1hi) = stumpff_range(I32F32::from_num(20.18), I32F32::from_num(20.2))[1];
        assert!(c1lo.to_num::<f64>() > -0.21724 && c1hi.to_num::<f64>() < -0.21723);
        let (lo, hi) = (I32F32::from_num(-50), I32F32::MAX);
        let whole = stumpff_range(lo, hi);
        for x in across(lo, hi, 4000).chain(across(I32F32::ZERO, I32F32::from_num(500), 4000)) {
            for (range, v) in whole.iter().zip(stumpff(x)) {
                assert!(holds(*range, v), "{v} at {x} outside {range:?}");
            }
        }
        let wide = stumpff_range(I32F32::from_num(4), I32F32::from_num(300))[1];
        assert!(wide.0.to_num::<f64>() > -0.21731 && wide.1.to_num::<f64>() < 0.455);
        let (c0lo, c0hi) = stumpff_range(I32F32::from_num(-3), I32F32::from_num(2))[0];
        assert!(
            c0lo == stumpff(I32F32::from_num(2))[0] - I32F32::DELTA
                && c0hi == stumpff(I32F32::from_num(-3))[0] + I32F32::DELTA
        );
    }

    #[test]
    fn tan_ratio_range_spans_its_poles() {
        let (lo, hi) = (I32F32::from_num(1.5), I32F32::from_num(1.6));
        assert_eq!(tan_ratio_range(lo, hi), (I32F32::MIN, I32F32::MAX));
        assert_eq!(
            tan_ratio_range(I32F32::from_num(-2), I32F32::from_num(2)),
            (I32F32::MIN, I32F32::MAX)
        );
        let (one, _) = tan_ratio_range(I32F32::from_num(-1), I32F32::from_num(0.5));
        assert_eq!(one, I32F32::ONE);
        let far = tan_ratio_range(I32F32::from_num(3.5), I32F32::from_num(4.5));
        assert!(far.0 > I32F32::ZERO && far.1 < I32F32::from_num(1.1));
    }

    #[test]
    fn box_ranges_hold_every_point() {
        let mut rng = Lcg(0xB0C5);
        for _ in 0..1500 {
            let corner = |g: &mut Lcg| {
                let a = g.range(-4.0, 4.0);
                (I32F32::from_num(a), I32F32::from_num(a + g.range(0.0, 2.0)))
            };
            let (x, y) = (corner(&mut rng), corner(&mut rng));
            let (h, angle) = (hypot_range(x, y), atan2_range(y, x));
            for px in across(x.0, x.1, 16) {
                for py in across(y.0, y.1, 16) {
                    assert!(holds(h, hypot(px, py)), "hypot({px}, {py})");
                    assert!(
                        holds(angle, atan2(py, px)),
                        "atan2({py}, {px}) outside {angle:?}"
                    );
                }
            }
        }
        let pi = I32F32::PI + I32F32::DELTA;
        let across_cut = atan2_range(
            (I32F32::from_num(-1), I32F32::ONE),
            (I32F32::from_num(-2), I32F32::from_num(-1)),
        );
        assert!(across_cut.0 <= -pi + I32F32::DELTA && across_cut.1 >= pi - I32F32::DELTA);
        let quadrant = atan2_range(
            (I32F32::ONE, I32F32::from_num(2)),
            (I32F32::ONE, I32F32::from_num(2)),
        );
        assert!(quadrant.0 > I32F32::from_num(0.46) && quadrant.1 < I32F32::from_num(1.11));
    }
}
