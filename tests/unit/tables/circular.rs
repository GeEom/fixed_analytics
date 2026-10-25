#[cfg(test)]
mod tests {
    use fixed_analytics::tables::circular::ATAN_TABLE;

    #[test]
    fn atan_table_has_64_entries() {
        assert_eq!(ATAN_TABLE.len(), 64);
    }

    #[test]
    fn atan_table_spot_check() {
        let [t0, t1, ..] = ATAN_TABLE;
        let atan_1: f64 = fixed::types::I1F63::from_bits(t0).to_num();
        let expected = core::f64::consts::FRAC_PI_4;
        assert!(
            (atan_1 - expected).abs() < 1e-15,
            "atan(1) = {atan_1}, expected {expected}"
        );

        let atan_half: f64 = fixed::types::I1F63::from_bits(t1).to_num();
        let expected_half = 0.5_f64.atan();
        assert!(
            (atan_half - expected_half).abs() < 1e-15,
            "atan(0.5) = {atan_half}, expected {expected_half}"
        );
    }

    #[test]
    fn atan_table_decreasing_until_convergence() {
        for (i, (prev, cur)) in ATAN_TABLE.iter().zip(ATAN_TABLE.iter().skip(1)).enumerate() {
            assert!(
                cur <= prev,
                "ATAN_TABLE[{}] = {cur} should be <= ATAN_TABLE[{i}] = {prev}",
                i + 1
            );
        }
    }

    #[test]
    fn atan_table_strictly_decreasing_early() {
        for (i, (prev, cur)) in ATAN_TABLE
            .iter()
            .zip(ATAN_TABLE.iter().skip(1))
            .enumerate()
            .take(59)
        {
            assert!(
                cur < prev,
                "ATAN_TABLE[{}] = {cur} should be < ATAN_TABLE[{i}] = {prev}",
                i + 1
            );
        }
    }
}
