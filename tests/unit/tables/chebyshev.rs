//! Tests for Chebyshev polynomial coefficient tables

#[cfg(test)]
#[allow(
    clippy::cast_precision_loss,
    reason = "test code uses f64 casts for verification"
)]
mod tests {
    use fixed_analytics::tables::chebyshev::{COS_Q_HI, COS_Q_LO, SIN_P_HI, SIN_P_LO};

    const SCALE: f64 = (1_u64 << 63) as f64;

    fn i1f63_to_f64(bits: i64) -> f64 {
        (bits as f64) / SCALE
    }

    // The constant term (last element, c₀) of each polynomial should match
    // the Taylor series leading coefficient to high precision. This guards
    // against coefficient ordering mistakes and regeneration drift.

    #[test]
    fn sin_constant_term_is_neg_one_sixth() {
        // (sin(x)-x)/x³ → -1/6 at x=0
        let expected = -1.0 / 6.0;
        let lo = i1f63_to_f64(*SIN_P_LO.last().unwrap_or(&0));
        let hi = i1f63_to_f64(*SIN_P_HI.last().unwrap_or(&0));
        assert!(
            (lo - expected).abs() < 1e-6,
            "SIN_P_LO constant = {lo}, expected {expected}"
        );
        assert!(
            (hi - expected).abs() < 1e-15,
            "SIN_P_HI constant = {hi}, expected {expected}"
        );
    }

    #[test]
    fn cos_constant_term_is_neg_one_half() {
        // (cos(x)-1)/x² → -1/2 at x=0
        let expected = -0.5;
        let lo = i1f63_to_f64(*COS_Q_LO.last().unwrap_or(&0));
        let hi = i1f63_to_f64(*COS_Q_HI.last().unwrap_or(&0));
        assert!(
            (lo - expected).abs() < 1e-6,
            "COS_Q_LO constant = {lo}, expected {expected}"
        );
        assert!(
            (hi - expected).abs() < 1e-15,
            "COS_Q_HI constant = {hi}, expected {expected}"
        );
    }

    #[test]
    fn all_coefficients_magnitude_below_one() {
        for (name, table) in [
            ("SIN_P_LO", SIN_P_LO.as_slice()),
            ("SIN_P_HI", SIN_P_HI.as_slice()),
            ("COS_Q_LO", COS_Q_LO.as_slice()),
            ("COS_Q_HI", COS_Q_HI.as_slice()),
        ] {
            for (i, &bits) in table.iter().enumerate() {
                let val = i1f63_to_f64(bits).abs();
                assert!(val < 1.0, "{name}[{i}] = {val}, exceeds I1F63 range");
            }
        }
    }

    #[test]
    fn expected_array_lengths() {
        assert_eq!(SIN_P_LO.len(), 4);
        assert_eq!(SIN_P_HI.len(), 7);
        assert_eq!(COS_Q_LO.len(), 4);
        assert_eq!(COS_Q_HI.len(), 7);
    }
}
