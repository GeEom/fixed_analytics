//! Tests for circular CORDIC lookup tables

#[cfg(test)]
#[allow(clippy::indexing_slicing, clippy::cast_precision_loss)]
mod tests {
    use fixed_analytics::tables::circular::{ATAN_TABLE, CIRCULAR_GAIN_INV};

    #[test]
    fn atan_table_has_64_entries() {
        assert_eq!(ATAN_TABLE.len(), 64);
    }

    #[test]
    fn atan_table_spot_check() {
        // ATAN_TABLE[0] = atan(2^0) = atan(1) = π/4 ≈ 0.7853981633974483
        // In I1F63 format: π/4 * 2^63
        // Use u64 for 2^63 to avoid i64 overflow
        let scale = (1_u64 << 63) as f64;
        let atan_1: f64 = (ATAN_TABLE[0] as f64) / scale;
        let expected = core::f64::consts::FRAC_PI_4;
        assert!(
            (atan_1 - expected).abs() < 1e-15,
            "atan(1) = {atan_1}, expected {expected}"
        );

        // ATAN_TABLE[1] = atan(2^(-1)) = atan(0.5) ≈ 0.4636476090008061
        let atan_half: f64 = (ATAN_TABLE[1] as f64) / scale;
        let expected_half = 0.5_f64.atan();
        assert!(
            (atan_half - expected_half).abs() < 1e-15,
            "atan(0.5) = {atan_half}, expected {expected_half}"
        );
    }

    #[test]
    fn atan_table_decreasing_until_convergence() {
        // Each entry should be smaller than or equal to the previous
        // (last few entries may be equal due to precision limits)
        for i in 1..ATAN_TABLE.len() {
            assert!(
                ATAN_TABLE[i] <= ATAN_TABLE[i - 1],
                "ATAN_TABLE[{}] = {} should be <= ATAN_TABLE[{}] = {}",
                i,
                ATAN_TABLE[i],
                i - 1,
                ATAN_TABLE[i - 1]
            );
        }
    }

    #[test]
    fn atan_table_strictly_decreasing_early() {
        // First 60 entries should be strictly decreasing
        for i in 1..60 {
            assert!(
                ATAN_TABLE[i] < ATAN_TABLE[i - 1],
                "ATAN_TABLE[{}] = {} should be < ATAN_TABLE[{}] = {}",
                i,
                ATAN_TABLE[i],
                i - 1,
                ATAN_TABLE[i - 1]
            );
        }
    }

    #[test]
    fn cordic_scale_factor_value() {
        // 1/K ≈ 0.6073 in I1F63 format
        // Expected value: 0x4DBA_76D4_21AF_2D34
        assert_eq!(CIRCULAR_GAIN_INV, 0x4DBA_76D4_21AF_2D34);
    }
}
