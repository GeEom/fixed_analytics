//! Tests for CORDIC core algorithms

#[cfg(test)]
mod tests {
    use fixed::types::{I16F16, I64F64};
    use fixed_analytics::kernel::circular_vectoring;

    #[test]
    fn circular_vectoring_atan_one() {
        // vectoring mode with x=1, y=1 should give z ≈ π/4
        let (_, _, z) = circular_vectoring(I16F16::ONE, I16F16::ONE, I16F16::ZERO);
        let z_f32: f32 = z.to_num();
        let expected = core::f32::consts::FRAC_PI_4;
        assert!((z_f32 - expected).abs() < 0.01);
    }

    #[test]
    fn circular_vectoring_high_precision_type() {
        // I64F64 has 64 fractional bits, more than the I1F63 tables:
        // exercises the exact (no rounding needed) conversion path.
        let (_, _, z) = circular_vectoring(I64F64::ONE, I64F64::ONE, I64F64::ZERO);
        let z_f64: f64 = z.to_num();
        let expected = core::f64::consts::FRAC_PI_4;
        assert!(
            (z_f64 - expected).abs() < 1e-15,
            "atan(1) at I64F64 = {z_f64}, expected {expected}"
        );
    }
}

/// Sweep both early-terminating kernels against `f64` at three widths.
#[cfg(test)]
mod early_termination {
    use crate::unit::support::Lcg;
    use fixed::types::{I16F16, I32F32, I64F64};
    use fixed_analytics::CordicNumber;
    use fixed_analytics::kernel::{circular_vectoring, hyperbolic_vectoring};

    fn sweep<T: CordicNumber + fixed::traits::Fixed>(tol_circ: f64, tol_hyp: f64) {
        let mut rng = Lcg(0xC0DE);
        for i in 0..=2000 {
            let v = if i <= 1000 {
                -1.0 + 2.0 * f64::from(i) / 1000.0
            } else {
                rng.range(-1.0, 1.0)
            };
            let y = <T as fixed::traits::Fixed>::from_num(v);
            let z: f64 = circular_vectoring(T::one(), y, T::zero()).2.to_num();
            let want = v.atan();
            assert!(
                (z - want).abs() < tol_circ,
                "atan({v}) via kernel = {z}, want {want} (tol {tol_circ})"
            );
            if v.abs() <= 0.75 {
                let zh: f64 = hyperbolic_vectoring(T::one(), y, T::zero()).2.to_num();
                let want_h = v.atanh();
                assert!(
                    (zh - want_h).abs() < tol_hyp,
                    "atanh({v}) via kernel = {zh}, want {want_h} (tol {tol_hyp})"
                );
            }
        }
    }

    #[test]
    fn kernels_track_f64_at_i16f16() {
        sweep::<I16F16>(8e-5, 8e-5);
    }

    #[test]
    fn kernels_track_f64_at_i32f32() {
        sweep::<I32F32>(3e-9, 3e-9);
    }

    #[test]
    fn kernels_track_f64_at_i64f64() {
        // Above an f64 ulp: the libm reference varies by platform.
        sweep::<I64F64>(6e-16, 6e-16);
    }
}
