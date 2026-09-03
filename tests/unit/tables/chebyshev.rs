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

/// `horner` wraps on the strength of every intermediate staying below 1 in
/// magnitude: check the bound and the bit-identity with saturating Horner.
#[cfg(test)]
#[allow(
    clippy::cast_precision_loss,
    clippy::unwrap_used,
    reason = "test code uses f64 for verification"
)]
mod wrapping_invariant {
    use fixed::types::{I16F16, I32F32, I64F64};
    use fixed_analytics::CordicNumber;
    use fixed_analytics::tables::chebyshev::{COS_Q_HI, COS_Q_LO, SIN_P_HI, SIN_P_LO, horner};

    const SCALE: f64 = (1_u64 << 63) as f64;
    const U_MAX: f64 = core::f64::consts::FRAC_PI_4 * core::f64::consts::FRAC_PI_4;

    fn max_intermediate(coeffs: &[i64]) -> f64 {
        let mut worst: f64 = 0.0;
        for i in 0..=1000 {
            let u = U_MAX * f64::from(i) / 1000.0;
            let mut iter = coeffs.iter();
            let mut acc = *iter.next().unwrap_or(&0) as f64 / SCALE;
            worst = worst.max(acc.abs());
            for &c in iter {
                acc = u.mul_add(acc, c as f64 / SCALE);
                worst = worst.max(acc.abs());
                worst = worst.max((u * acc).abs());
            }
        }
        worst
    }

    #[test]
    fn intermediates_stay_well_below_one() {
        for table in [&SIN_P_LO[..], &SIN_P_HI[..], &COS_Q_LO[..], &COS_Q_HI[..]] {
            let worst = max_intermediate(table);
            assert!(worst < 0.6, "Horner intermediate reached {worst}");
        }
    }

    fn saturating_horner<T: CordicNumber, const N: usize>(coeffs: &[i64; N], x: T) -> T {
        let mut iter = coeffs.iter();
        let mut result = T::from_i1f63(*iter.next().unwrap());
        for &c in iter {
            result = T::from_i1f63(c).saturating_add(x.saturating_mul(result));
        }
        result
    }

    fn check<T: CordicNumber + core::fmt::Debug>() {
        for i in 0..=2000 {
            let u = T::from_num(U_MAX * f64::from(i) / 2000.0);
            assert_eq!(horner(&SIN_P_LO, u), saturating_horner(&SIN_P_LO, u));
            assert_eq!(horner(&SIN_P_HI, u), saturating_horner(&SIN_P_HI, u));
            assert_eq!(horner(&COS_Q_LO, u), saturating_horner(&COS_Q_LO, u));
            assert_eq!(horner(&COS_Q_HI, u), saturating_horner(&COS_Q_HI, u));
        }
    }

    #[test]
    fn wrapping_matches_saturating_bit_for_bit() {
        check::<I16F16>();
        check::<I32F32>();
        check::<I64F64>();
    }
}
