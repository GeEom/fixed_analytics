#[cfg(test)]
mod tests {
    use crate::unit::support::Lcg;
    use fixed::types::{I16F16, I32F32};
    use fixed_analytics::{atan_range, enclose, exp, exp_range, sin, sin_cos_range, tanh_range};

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
}
