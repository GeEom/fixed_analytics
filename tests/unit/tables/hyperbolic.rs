//! Tests for hyperbolic CORDIC lookup tables

#[cfg(test)]
#[allow(
    clippy::indexing_slicing,
    clippy::cast_precision_loss,
    reason = "test code uses direct indexing and f64 casts"
)]
mod tests {
    use fixed_analytics::tables::hyperbolic::{ATANH_HALF, ATANH_TABLE, needs_repeat};

    /// Repeat indices for hyperbolic CORDIC convergence (used only in tests).
    const REPEAT_INDICES: [u32; 5] = [4, 13, 40, 121, 364];

    #[test]
    fn atanh_table_has_64_entries() {
        assert_eq!(ATANH_TABLE.len(), 64);
    }

    #[test]
    fn atanh_table_decreasing_until_convergence() {
        // Each entry should be smaller than or equal to the previous
        // (last few entries may be equal due to precision limits)
        for i in 1..ATANH_TABLE.len() {
            assert!(
                ATANH_TABLE[i] <= ATANH_TABLE[i - 1],
                "ATANH_TABLE[{}] = {} should be <= ATANH_TABLE[{}] = {}",
                i,
                ATANH_TABLE[i],
                i - 1,
                ATANH_TABLE[i - 1]
            );
        }
    }

    #[test]
    fn atanh_table_strictly_decreasing_early() {
        // First 60 entries should be strictly decreasing
        for i in 1..60 {
            assert!(
                ATANH_TABLE[i] < ATANH_TABLE[i - 1],
                "ATANH_TABLE[{}] = {} should be < ATANH_TABLE[{}] = {}",
                i,
                ATANH_TABLE[i],
                i - 1,
                ATANH_TABLE[i - 1]
            );
        }
    }

    #[test]
    fn needs_repeat_correct_indices() {
        // Verify needs_repeat returns true for expected indices
        for &idx in &REPEAT_INDICES {
            assert!(needs_repeat(idx), "needs_repeat({idx}) should be true");
        }
        // Verify needs_repeat returns false for non-repeat indices
        assert!(!needs_repeat(5));
        assert!(!needs_repeat(100));
        assert!(!needs_repeat(0));
        assert!(!needs_repeat(3));
    }

    #[test]
    fn atanh_half_matches_table() {
        // ATANH_HALF should equal ATANH_TABLE[0]
        assert_eq!(ATANH_HALF, ATANH_TABLE[0]);
    }

    #[test]
    fn atanh_table_spot_check() {
        // ATANH_TABLE[0] = atanh(2^(-1)) = atanh(0.5) ≈ 0.5493061443340548
        // In I1F63 format: atanh(0.5) * 2^63
        // Use u64 for 2^63 to avoid i64 overflow
        let scale = (1_u64 << 63) as f64;
        let atanh_half: f64 = (ATANH_TABLE[0] as f64) / scale;
        let expected = 0.5_f64.atanh();
        assert!(
            (atanh_half - expected).abs() < 1e-15,
            "atanh(0.5) = {atanh_half}, expected {expected}"
        );

        // ATANH_TABLE[1] = atanh(2^(-2)) = atanh(0.25) ≈ 0.2554128118829953
        let atanh_quarter: f64 = (ATANH_TABLE[1] as f64) / scale;
        let expected_quarter = 0.25_f64.atanh();
        assert!(
            (atanh_quarter - expected_quarter).abs() < 1e-15,
            "atanh(0.25) = {atanh_quarter}, expected {expected_quarter}"
        );
    }
}
