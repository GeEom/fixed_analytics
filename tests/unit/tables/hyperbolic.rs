#[cfg(test)]
mod tests {
    use fixed_analytics::tables::hyperbolic::{ATANH_HALF, ATANH_TABLE, needs_repeat};

    const REPEAT_INDICES: [u32; 5] = [4, 13, 40, 121, 364];

    #[test]
    fn atanh_table_has_64_entries() {
        assert_eq!(ATANH_TABLE.len(), 64);
    }

    #[test]
    fn atanh_table_decreasing_until_convergence() {
        for (i, (prev, cur)) in ATANH_TABLE
            .iter()
            .zip(ATANH_TABLE.iter().skip(1))
            .enumerate()
        {
            assert!(
                cur <= prev,
                "ATANH_TABLE[{}] = {cur} should be <= ATANH_TABLE[{i}] = {prev}",
                i + 1
            );
        }
    }

    #[test]
    fn atanh_table_strictly_decreasing_early() {
        for (i, (prev, cur)) in ATANH_TABLE
            .iter()
            .zip(ATANH_TABLE.iter().skip(1))
            .enumerate()
            .take(59)
        {
            assert!(
                cur < prev,
                "ATANH_TABLE[{}] = {cur} should be < ATANH_TABLE[{i}] = {prev}",
                i + 1
            );
        }
    }

    #[test]
    fn needs_repeat_correct_indices() {
        for &idx in &REPEAT_INDICES {
            assert!(needs_repeat(idx), "needs_repeat({idx}) should be true");
        }
        assert!(!needs_repeat(5));
        assert!(!needs_repeat(100));
        assert!(!needs_repeat(0));
        assert!(!needs_repeat(3));
    }

    #[test]
    fn atanh_half_matches_table() {
        assert_eq!(Some(&ATANH_HALF), ATANH_TABLE.first());
    }

    #[test]
    fn atanh_table_spot_check() {
        let [t0, t1, ..] = ATANH_TABLE;
        let atanh_half: f64 = fixed::types::I1F63::from_bits(t0).to_num();
        let expected = 0.5_f64.atanh();
        assert!(
            (atanh_half - expected).abs() < 1e-15,
            "atanh(0.5) = {atanh_half}, expected {expected}"
        );

        let atanh_quarter: f64 = fixed::types::I1F63::from_bits(t1).to_num();
        let expected_quarter = 0.25_f64.atanh();
        assert!(
            (atanh_quarter - expected_quarter).abs() < 1e-15,
            "atanh(0.25) = {atanh_quarter}, expected {expected_quarter}"
        );
    }
}
