//! Tests for CORDIC core algorithms

#[cfg(test)]
mod tests {
    use fixed::types::I16F16;
    use fixed_analytics::kernel::circular_vectoring;

    #[test]
    fn circular_vectoring_atan_one() {
        // vectoring mode with x=1, y=1 should give z ≈ π/4
        let (_, _, z) = circular_vectoring(I16F16::ONE, I16F16::ONE, I16F16::ZERO);
        let z_f32: f32 = z.to_num();
        let expected = core::f32::consts::FRAC_PI_4;
        assert!((z_f32 - expected).abs() < 0.01);
    }
}
