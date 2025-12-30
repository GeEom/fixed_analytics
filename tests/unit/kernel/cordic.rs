//! Tests for CORDIC core algorithms

#[cfg(test)]
mod tests {
    use fixed::types::I16F16;
    use fixed_analytics::kernel::{circular_gain_inv, circular_rotation, circular_vectoring};

    #[test]
    fn circular_rotation_zero_angle() {
        let inv_gain = circular_gain_inv::<I16F16>();
        let (x, y, z) = circular_rotation(inv_gain, I16F16::ZERO, I16F16::ZERO);
        // After rotation by 0, x should be close to 1 (after gain compensation), y should be ~0
        let x_f32: f32 = x.to_num();
        let y_f32: f32 = y.to_num();
        let z_f32: f32 = z.to_num();
        // The result depends on gain compensation; x should be close to 1
        assert!((x_f32 - 1.0).abs() < 0.02, "x = {x_f32}, expected ~1.0");
        assert!(y_f32.abs() < 0.01, "y = {y_f32}, expected ~0");
        assert!(z_f32.abs() < 0.01, "z = {z_f32}, expected ~0");
    }

    #[test]
    fn circular_vectoring_atan_one() {
        // vectoring mode with x=1, y=1 should give z ≈ π/4
        let (_, _, z) = circular_vectoring(I16F16::ONE, I16F16::ONE, I16F16::ZERO);
        let z_f32: f32 = z.to_num();
        let expected = core::f32::consts::FRAC_PI_4;
        assert!((z_f32 - expected).abs() < 0.01);
    }
}
