//! Tests for `CordicNumber` trait implementations

#[cfg(test)]
mod tests {
    use fixed::types::{I4F12, I4F60, I8F8, I8F24, I16F16, I20F12, I24F8, I32F32, I48F16, I64F64};
    use fixed_analytics::CordicNumber;
    use fixed_analytics::kernel::hyperbolic_gain_inv;

    #[test]
    #[allow(clippy::approx_constant)]
    fn basic_operations_i16f16() {
        let x = I16F16::from_num(2.5);
        assert_eq!(I16F16::zero(), I16F16::ZERO);
        assert_eq!(I16F16::one(), I16F16::ONE);
        assert!(I16F16::pi() > I16F16::from_num(3.14));
        assert!(I16F16::pi() < I16F16::from_num(3.15));
        assert!(!x.is_negative());
        assert!((-x).is_negative());
    }

    #[test]
    fn generic_impl_works_for_various_types() {
        // Test that the generic impl works for various type configurations

        // I8F24: 8 integer bits, 24 fractional bits (32-bit total)
        let _: I8F24 = I8F24::pi();
        let _: I8F24 = I8F24::frac_pi_2();

        // I24F8: 24 integer bits, 8 fractional bits
        let _: I24F8 = I24F8::pi();
        let _: I24F8 = I24F8::frac_pi_2();

        // I4F12: 4 integer bits, 12 fractional bits (16-bit total)
        let _: I4F12 = I4F12::pi();

        // I20F12: 20 integer bits, 12 fractional bits (32-bit total)
        let _: I20F12 = I20F12::pi();

        // I48F16: 48 integer bits, 16 fractional bits (64-bit total)
        let _: I48F16 = I48F16::pi();

        // I4F60: 4 integer bits, 60 fractional bits (64-bit total) - high precision
        let _: I4F60 = I4F60::pi();
    }

    #[test]
    fn from_i1f63_across_types() {
        // Test that from_i1f63 works correctly across different types
        // 0.5 in I1F63 format
        let half_bits: i64 = 0x4000_0000_0000_0000;

        let i8f8_half: f32 = I8F8::from_i1f63(half_bits).to_num();
        assert!((i8f8_half - 0.5).abs() < 0.01);

        let i16f16_half: f32 = I16F16::from_i1f63(half_bits).to_num();
        assert!((i16f16_half - 0.5).abs() < 0.0001);

        let i32f32_half: f64 = I32F32::from_i1f63(half_bits).to_num();
        assert!((i32f32_half - 0.5).abs() < 1e-9);

        let i64f64_half: f64 = I64F64::from_i1f63(half_bits).to_num();
        assert!((i64f64_half - 0.5).abs() < 1e-15);

        // Also test a non-standard type
        let i24f8_half: f32 = I24F8::from_i1f63(half_bits).to_num();
        assert!((i24f8_half - 0.5).abs() < 0.01);

        let i4f60_half: f64 = I4F60::from_i1f63(half_bits).to_num();
        assert!((i4f60_half - 0.5).abs() < 1e-15);
    }

    #[test]
    fn frac_bits_correct() {
        assert_eq!(I8F8::frac_bits(), 8);
        assert_eq!(I16F16::frac_bits(), 16);
        assert_eq!(I32F32::frac_bits(), 32);
        assert_eq!(I64F64::frac_bits(), 64);

        // Non-standard types
        assert_eq!(I8F24::frac_bits(), 24);
        assert_eq!(I24F8::frac_bits(), 8);
        assert_eq!(I4F12::frac_bits(), 12);
        assert_eq!(I48F16::frac_bits(), 16);
    }

    #[test]
    fn frac_pi_4_values() {
        // Test the frac_pi_4() default implementation
        // π/4 ≈ 0.7854
        let pi_4_16: f32 = I16F16::frac_pi_4().to_num();
        assert!(
            (pi_4_16 - core::f32::consts::FRAC_PI_4).abs() < 0.001,
            "I16F16::frac_pi_4() = {pi_4_16}, expected ~0.7854"
        );

        let pi_4_32: f64 = I32F32::frac_pi_4().to_num();
        assert!(
            (pi_4_32 - core::f64::consts::FRAC_PI_4).abs() < 1e-9,
            "I32F32::frac_pi_4() = {pi_4_32}, expected ~0.7854"
        );

        let pi_4_64: f64 = I64F64::frac_pi_4().to_num();
        assert!(
            (pi_4_64 - core::f64::consts::FRAC_PI_4).abs() < 1e-15,
            "I64F64::frac_pi_4() = {pi_4_64}, expected ~0.7854"
        );
    }

    #[test]
    fn from_i2f62_high_precision() {
        // Test from_i2f62 with high-precision types (frac_bits > 62)
        // This exercises the shift-left branch in from_i2f62
        // 1/K_h ≈ 1.2075 is stored in I2F62 format
        // We call hyperbolic_gain_inv which uses from_i2f62

        // I64F64 has 64 fractional bits > 62, so it should use the shift-left branch
        let gain_inv_64: f64 = hyperbolic_gain_inv::<I64F64>().to_num();
        assert!(
            (gain_inv_64 - 1.2075).abs() < 0.01,
            "hyperbolic_gain_inv::<I64F64>() = {gain_inv_64}, expected ~1.2075"
        );

        // I4F60 has 60 fractional bits < 62, uses shift-right branch (control test)
        let gain_inv_60: f64 = hyperbolic_gain_inv::<I4F60>().to_num();
        assert!(
            (gain_inv_60 - 1.2075).abs() < 0.01,
            "hyperbolic_gain_inv::<I4F60>() = {gain_inv_60}, expected ~1.2075"
        );
    }
}
