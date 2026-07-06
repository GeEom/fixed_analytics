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
