//! Mathematical correctness verification tests.
//!
//! These tests verify the CORDIC implementations against known mathematical
//! properties: reference values, identities, inverse roundtrips, monotonicity,
//! and output bounds.

// Test-specific lints - these are acceptable in test code
#![allow(
    clippy::unwrap_used,
    clippy::cast_possible_truncation,
    clippy::suboptimal_flops,
    clippy::cast_lossless,
    clippy::manual_range_contains,
    reason = "test code uses these patterns for conciseness and clarity"
)]

#[cfg(test)]
mod reference_comparison {
    //! Compare against f64 reference implementations across sampled inputs.

    use fixed::types::I16F16;
    use fixed_analytics::{
        acos, acosh, asin, asinh, atan, atan2, atanh, exp, ln, log2, log10, sin_cos, sinh_cosh,
        sqrt, tan, tanh,
    };

    /// Deterministic pseudo-random bit generator for reproducible sampling.
    fn sample_bits(seed: u64, index: u64) -> i32 {
        // Simple LCG-style mixing
        let mut x = seed.wrapping_add(index.wrapping_mul(0x9E37_79B9_7F4A_7C15));
        x = x.wrapping_mul(0x5851_F42D_4C95_7F2D);
        x = x ^ (x >> 32);
        x as i32
    }

    const SAMPLES: u64 = 2000;
    const SEED: u64 = 0xDEAD_BEEF_CAFE_BABE;

    // Error tolerance: I16F16 has 16 fractional bits, so ~4-5 decimal digits
    // Allow 2^-12 ≈ 0.000244 absolute error for well-conditioned functions
    const TRIG_TOL: f64 = 0.001;
    const SQRT_TOL: f64 = 0.0005;
    const EXP_TOL: f64 = 0.01; // exp amplifies errors
    const LOG_TOL: f64 = 0.01;
    const HYPER_TOL: f64 = 0.02; // hyperbolic has lower precision due to iteration repeats

    #[test]
    fn sin_cos_vs_f64() {
        let mut max_sin_err: f64 = 0.0;
        let mut max_cos_err: f64 = 0.0;

        for i in 0..SAMPLES {
            let bits = sample_bits(SEED, i);
            let x = I16F16::from_bits(bits);
            let x_f64: f64 = x.to_num();

            // Skip extreme values where f64 comparison becomes meaningless
            if x_f64.abs() > 100.0 {
                continue;
            }

            let (s, c) = sin_cos(x);
            let expected_sin = x_f64.sin();
            let expected_cos = x_f64.cos();

            let sin_err = (s.to_num::<f64>() - expected_sin).abs();
            let cos_err = (c.to_num::<f64>() - expected_cos).abs();

            max_sin_err = max_sin_err.max(sin_err);
            max_cos_err = max_cos_err.max(cos_err);

            assert!(
                sin_err < TRIG_TOL,
                "sin({x_f64}): got {}, expected {expected_sin}, err {sin_err}",
                s.to_num::<f64>()
            );
            assert!(
                cos_err < TRIG_TOL,
                "cos({x_f64}): got {}, expected {expected_cos}, err {cos_err}",
                c.to_num::<f64>()
            );
        }

        // Uncomment to see actual max errors during development:
        // println!("Max sin error: {max_sin_err}, max cos error: {max_cos_err}");
    }

    #[test]
    fn tan_vs_f64() {
        for i in 0..SAMPLES {
            let bits = sample_bits(SEED, i);
            let x = I16F16::from_bits(bits);
            let x_f64: f64 = x.to_num();

            // Skip near poles and extreme values
            let frac_pi_2 = core::f64::consts::FRAC_PI_2;
            if x_f64.abs() > 50.0 {
                continue;
            }
            // Skip near ±π/2 where tan explodes
            let reduced = x_f64 % core::f64::consts::PI;
            if (reduced.abs() - frac_pi_2).abs() < 0.1 {
                continue;
            }

            let result = tan(x);
            let expected = x_f64.tan();

            // Use relative error for large values
            let err = if expected.abs() > 1.0 {
                (result.to_num::<f64>() - expected).abs() / expected.abs()
            } else {
                (result.to_num::<f64>() - expected).abs()
            };

            assert!(
                err < TRIG_TOL * 2.0, // tan accumulates more error
                "tan({x_f64}): got {}, expected {expected}",
                result.to_num::<f64>()
            );
        }
    }

    #[test]
    fn atan_vs_f64() {
        for i in 0..SAMPLES {
            let bits = sample_bits(SEED, i);
            let x = I16F16::from_bits(bits);
            let x_f64: f64 = x.to_num();

            if x_f64.abs() > 1000.0 {
                continue;
            }

            let result = atan(x);
            let expected = x_f64.atan();
            let err = (result.to_num::<f64>() - expected).abs();

            assert!(
                err < TRIG_TOL,
                "atan({x_f64}): got {}, expected {expected}",
                result.to_num::<f64>()
            );
        }
    }

    #[test]
    fn atan2_vs_f64() {
        for i in 0..SAMPLES {
            let y_bits = sample_bits(SEED, i);
            let x_bits = sample_bits(SEED ^ 0xFFFF, i);
            let y = I16F16::from_bits(y_bits);
            let x = I16F16::from_bits(x_bits);
            let y_f64: f64 = y.to_num();
            let x_f64: f64 = x.to_num();

            // Skip extreme values
            if y_f64.abs() > 1000.0 || x_f64.abs() > 1000.0 {
                continue;
            }

            let result = atan2(y, x);
            let expected = y_f64.atan2(x_f64);
            let err = (result.to_num::<f64>() - expected).abs();

            assert!(
                err < TRIG_TOL,
                "atan2({y_f64}, {x_f64}): got {}, expected {expected}",
                result.to_num::<f64>()
            );
        }
    }

    #[test]
    fn asin_acos_vs_f64() {
        // Test over domain [-1, 1]
        for i in 0..500 {
            // Map to [-1, 1]
            let t = (i as f64) / 499.0 * 2.0 - 1.0;
            // Avoid exact ±1 where precision is worst
            let t = t * 0.999;
            let x = I16F16::from_num(t);
            let x_f64: f64 = x.to_num();

            if let Ok(result) = asin(x) {
                let expected = x_f64.asin();
                let err = (result.to_num::<f64>() - expected).abs();
                assert!(
                    err < TRIG_TOL,
                    "asin({x_f64}): got {}, expected {expected}",
                    result.to_num::<f64>()
                );
            }

            if let Ok(result) = acos(x) {
                let expected = x_f64.acos();
                let err = (result.to_num::<f64>() - expected).abs();
                assert!(
                    err < TRIG_TOL,
                    "acos({x_f64}): got {}, expected {expected}",
                    result.to_num::<f64>()
                );
            }
        }
    }

    #[test]
    fn sqrt_vs_f64() {
        for i in 0..SAMPLES {
            let bits = sample_bits(SEED, i);
            // Only test non-negative values
            let bits = bits & 0x7FFF_FFFF;
            let x = I16F16::from_bits(bits);
            let x_f64: f64 = x.to_num();

            let result = sqrt(x).unwrap();
            let expected = x_f64.sqrt();
            let err = (result.to_num::<f64>() - expected).abs();

            // Use relative error for larger values
            let tol = if expected > 1.0 {
                SQRT_TOL * expected
            } else {
                SQRT_TOL
            };

            assert!(
                err < tol,
                "sqrt({x_f64}): got {}, expected {expected}",
                result.to_num::<f64>()
            );
        }
    }

    #[test]
    fn exp_vs_f64() {
        for i in 0..SAMPLES {
            let bits = sample_bits(SEED, i);
            let x = I16F16::from_bits(bits);
            let x_f64: f64 = x.to_num();

            // Skip values where exp would overflow I16F16 (max ~32768)
            // ln(32768) ≈ 10.4
            if x_f64 > 10.0 || x_f64 < -20.0 {
                continue;
            }

            let result = exp(x);
            let expected = x_f64.exp();

            // Use relative error since exp spans many orders of magnitude
            let err = if expected > 0.01 {
                (result.to_num::<f64>() - expected).abs() / expected
            } else {
                (result.to_num::<f64>() - expected).abs()
            };

            assert!(
                err < EXP_TOL,
                "exp({x_f64}): got {}, expected {expected}",
                result.to_num::<f64>()
            );
        }
    }

    #[test]
    fn ln_vs_f64() {
        for i in 0..SAMPLES {
            let bits = sample_bits(SEED, i);
            // Only positive values
            let bits = (bits & 0x7FFF_FFFF).max(1);
            let x = I16F16::from_bits(bits);
            let x_f64: f64 = x.to_num();

            if x_f64 < 0.001 {
                continue; // ln near 0 is extreme
            }

            if let Ok(result) = ln(x) {
                let expected = x_f64.ln();
                let err = (result.to_num::<f64>() - expected).abs();

                assert!(
                    err < LOG_TOL,
                    "ln({x_f64}): got {}, expected {expected}",
                    result.to_num::<f64>()
                );
            }
        }
    }

    #[test]
    fn log2_log10_vs_f64() {
        for i in 0..500 {
            let bits = sample_bits(SEED, i);
            let bits = (bits & 0x7FFF_FFFF).max(1);
            let x = I16F16::from_bits(bits);
            let x_f64: f64 = x.to_num();

            if x_f64 < 0.01 {
                continue;
            }

            if let Ok(result) = log2(x) {
                let expected = x_f64.log2();
                let err = (result.to_num::<f64>() - expected).abs();
                assert!(
                    err < LOG_TOL,
                    "log2({x_f64}): got {}, expected {expected}",
                    result.to_num::<f64>()
                );
            }

            if let Ok(result) = log10(x) {
                let expected = x_f64.log10();
                let err = (result.to_num::<f64>() - expected).abs();
                assert!(
                    err < LOG_TOL,
                    "log10({x_f64}): got {}, expected {expected}",
                    result.to_num::<f64>()
                );
            }
        }
    }

    #[test]
    fn sinh_cosh_vs_f64() {
        for i in 0..SAMPLES {
            let bits = sample_bits(SEED, i);
            let x = I16F16::from_bits(bits);
            let x_f64: f64 = x.to_num();

            // Skip large values where overflow occurs
            if x_f64.abs() > 8.0 {
                continue;
            }

            let (sh, ch) = sinh_cosh(x);
            let expected_sinh = x_f64.sinh();
            let expected_cosh = x_f64.cosh();

            // Use relative error for large results
            let sinh_err = if expected_sinh.abs() > 1.0 {
                (sh.to_num::<f64>() - expected_sinh).abs() / expected_sinh.abs()
            } else {
                (sh.to_num::<f64>() - expected_sinh).abs()
            };
            let cosh_err = (ch.to_num::<f64>() - expected_cosh).abs() / expected_cosh;

            assert!(
                sinh_err < HYPER_TOL,
                "sinh({x_f64}): got {}, expected {expected_sinh}",
                sh.to_num::<f64>()
            );
            assert!(
                cosh_err < HYPER_TOL,
                "cosh({x_f64}): got {}, expected {expected_cosh}",
                ch.to_num::<f64>()
            );
        }
    }

    #[test]
    fn tanh_vs_f64() {
        for i in 0..SAMPLES {
            let bits = sample_bits(SEED, i);
            let x = I16F16::from_bits(bits);
            let x_f64: f64 = x.to_num();

            if x_f64.abs() > 20.0 {
                continue;
            }

            let result = tanh(x);
            let expected = x_f64.tanh();
            let err = (result.to_num::<f64>() - expected).abs();

            assert!(
                err < HYPER_TOL,
                "tanh({x_f64}): got {}, expected {expected}",
                result.to_num::<f64>()
            );
        }
    }

    #[test]
    fn asinh_vs_f64() {
        for i in 0..SAMPLES {
            let bits = sample_bits(SEED, i);
            let x = I16F16::from_bits(bits);
            let x_f64: f64 = x.to_num();

            // asinh involves sqrt which loses precision for large values
            if x_f64.abs() > 20.0 {
                continue;
            }

            let result = asinh(x);
            let expected = x_f64.asinh();
            let err = (result.to_num::<f64>() - expected).abs();

            assert!(
                err < HYPER_TOL,
                "asinh({x_f64}): got {}, expected {expected}",
                result.to_num::<f64>()
            );
        }
    }

    #[test]
    fn acosh_vs_f64() {
        // Test over domain [1, moderate] - precision degrades for large values
        for i in 0..500 {
            // Map to [1, 20]
            let t = 1.0 + (i as f64) / 499.0 * 19.0;
            let x = I16F16::from_num(t);
            let x_f64: f64 = x.to_num();

            if let Ok(result) = acosh(x) {
                let expected = x_f64.acosh();
                let err = (result.to_num::<f64>() - expected).abs();
                assert!(
                    err < HYPER_TOL,
                    "acosh({x_f64}): got {}, expected {expected}",
                    result.to_num::<f64>()
                );
            }
        }
    }

    #[test]
    fn atanh_vs_f64() {
        // Test over domain (-1, 1)
        for i in 0..500 {
            // Map to (-0.99, 0.99)
            let t = (i as f64) / 499.0 * 1.98 - 0.99;
            let x = I16F16::from_num(t);
            let x_f64: f64 = x.to_num();

            if let Ok(result) = atanh(x) {
                let expected = x_f64.atanh();
                let err = (result.to_num::<f64>() - expected).abs();
                assert!(
                    err < HYPER_TOL,
                    "atanh({x_f64}): got {}, expected {expected}",
                    result.to_num::<f64>()
                );
            }
        }
    }
}

#[cfg(test)]
mod identities {
    //! Verify fundamental mathematical identities.

    use fixed::types::I16F16;
    use fixed_analytics::{exp, ln, sin_cos, sinh_cosh};

    fn sample_bits(seed: u64, index: u64) -> i32 {
        let mut x = seed.wrapping_add(index.wrapping_mul(0x9E37_79B9_7F4A_7C15));
        x = x.wrapping_mul(0x5851_F42D_4C95_7F2D);
        x = x ^ (x >> 32);
        x as i32
    }

    const SAMPLES: u64 = 1000;
    const SEED: u64 = 0x1234_5678_9ABC_DEF0;

    #[test]
    fn pythagorean_identity() {
        // sin²(x) + cos²(x) = 1
        for i in 0..SAMPLES {
            let bits = sample_bits(SEED, i);
            let x = I16F16::from_bits(bits);
            let x_f64: f64 = x.to_num();

            if x_f64.abs() > 100.0 {
                continue;
            }

            let (s, c) = sin_cos(x);
            let s_f64: f64 = s.to_num();
            let c_f64: f64 = c.to_num();
            let sum_sq = s_f64 * s_f64 + c_f64 * c_f64;

            assert!(
                (sum_sq - 1.0).abs() < 0.01,
                "sin²({x_f64}) + cos²({x_f64}) = {sum_sq}, expected 1.0"
            );
        }
    }

    #[test]
    fn hyperbolic_identity() {
        // cosh²(x) - sinh²(x) = 1
        for i in 0..SAMPLES {
            let bits = sample_bits(SEED, i);
            let x = I16F16::from_bits(bits);
            let x_f64: f64 = x.to_num();

            // Skip large values where overflow occurs
            if x_f64.abs() > 8.0 {
                continue;
            }

            let (sh, ch) = sinh_cosh(x);
            let sh_f64: f64 = sh.to_num();
            let ch_f64: f64 = ch.to_num();
            let diff = ch_f64 * ch_f64 - sh_f64 * sh_f64;

            assert!(
                (diff - 1.0).abs() < 0.02,
                "cosh²({x_f64}) - sinh²({x_f64}) = {diff}, expected 1.0"
            );
        }
    }

    #[test]
    fn angle_addition() {
        // sin(a + b) = sin(a)cos(b) + cos(a)sin(b)
        for i in 0..200 {
            for j in 0..200 {
                let a_bits = sample_bits(SEED, i);
                let b_bits = sample_bits(SEED ^ 0xFFFF, j);

                // Use smaller values to stay in good precision range
                let a = I16F16::from_bits(a_bits >> 8);
                let b = I16F16::from_bits(b_bits >> 8);

                let (sin_a, cos_a) = sin_cos(a);
                let (sin_b, cos_b) = sin_cos(b);
                let (sin_ab, _) = sin_cos(a + b);

                // sin(a+b) should equal sin(a)cos(b) + cos(a)sin(b)
                let expected = sin_a.to_num::<f64>() * cos_b.to_num::<f64>()
                    + cos_a.to_num::<f64>() * sin_b.to_num::<f64>();
                let actual: f64 = sin_ab.to_num();

                assert!(
                    (actual - expected).abs() < 0.02,
                    "sin({}+{}) = {actual}, expected {expected}",
                    a.to_num::<f64>(),
                    b.to_num::<f64>()
                );
            }
        }
    }

    #[test]
    fn double_angle() {
        // sin(2x) = 2*sin(x)*cos(x)
        // cos(2x) = cos²(x) - sin²(x)
        for i in 0..SAMPLES {
            let bits = sample_bits(SEED, i);
            // Keep x small to avoid angle reduction errors compounding
            let x = I16F16::from_bits(bits >> 8); // Much smaller range
            let x_f64: f64 = x.to_num();

            // Skip large angles where errors compound
            if x_f64.abs() > 50.0 {
                continue;
            }

            let two_x = x + x;

            let (sin_x, cos_x) = sin_cos(x);
            let (sin_2x, cos_2x) = sin_cos(two_x);

            let expected_sin_2x = 2.0 * sin_x.to_num::<f64>() * cos_x.to_num::<f64>();
            let expected_cos_2x = cos_x.to_num::<f64>().powi(2) - sin_x.to_num::<f64>().powi(2);

            let sin_err = (sin_2x.to_num::<f64>() - expected_sin_2x).abs();
            let cos_err = (cos_2x.to_num::<f64>() - expected_cos_2x).abs();

            assert!(
                sin_err < 0.02,
                "sin(2*{}) = {}, expected {expected_sin_2x}",
                x.to_num::<f64>(),
                sin_2x.to_num::<f64>()
            );
            assert!(
                cos_err < 0.02,
                "cos(2*{}) = {}, expected {expected_cos_2x}",
                x.to_num::<f64>(),
                cos_2x.to_num::<f64>()
            );
        }
    }

    #[test]
    fn exp_addition() {
        // exp(a + b) = exp(a) * exp(b)
        for i in 0..200 {
            for j in 0..200 {
                let a_bits = sample_bits(SEED, i);
                let b_bits = sample_bits(SEED ^ 0xFFFF, j);

                // Keep values small to avoid overflow
                // Map to roughly [-2, 2]
                let a = I16F16::from_bits(a_bits >> 14);
                let b = I16F16::from_bits(b_bits >> 14);

                let a_f64: f64 = a.to_num();
                let b_f64: f64 = b.to_num();

                // Skip if sum would overflow
                if (a_f64 + b_f64).abs() > 8.0 {
                    continue;
                }

                let exp_a = exp(a);
                let exp_b = exp(b);
                let exp_ab = exp(a + b);

                let expected = exp_a.to_num::<f64>() * exp_b.to_num::<f64>();
                let actual: f64 = exp_ab.to_num();

                let rel_err = (actual - expected).abs() / expected.max(0.001);

                assert!(
                    rel_err < 0.05,
                    "exp({a_f64}+{b_f64}) = {actual}, expected {expected}"
                );
            }
        }
    }

    #[test]
    fn ln_multiplication() {
        // ln(a * b) = ln(a) + ln(b)
        for i in 0..200 {
            for j in 0..200 {
                let a_bits = sample_bits(SEED, i);
                let b_bits = sample_bits(SEED ^ 0xFFFF, j);

                // Keep values in range [0.5, 128] for reasonable precision
                let a_bits = (a_bits & 0x007F_FFFF).max(0x0000_8000);
                let b_bits = (b_bits & 0x007F_FFFF).max(0x0000_8000);
                let a = I16F16::from_bits(a_bits);
                let b = I16F16::from_bits(b_bits);

                let a_f64: f64 = a.to_num();
                let b_f64: f64 = b.to_num();

                // Skip if product would overflow or underflow significantly
                let product = a_f64 * b_f64;
                if product > 10000.0 || product < 0.1 {
                    continue;
                }

                let ab = a.saturating_mul(b);

                if let (Ok(ln_a), Ok(ln_b), Ok(ln_ab)) = (ln(a), ln(b), ln(ab)) {
                    let expected = ln_a.to_num::<f64>() + ln_b.to_num::<f64>();
                    let actual: f64 = ln_ab.to_num();

                    assert!(
                        (actual - expected).abs() < 0.05,
                        "ln({a_f64}*{b_f64}) = {actual}, expected {expected}"
                    );
                }
            }
        }
    }
}

#[cfg(test)]
mod roundtrips {
    //! Verify that inverse functions are actually inverses.

    use fixed::types::I16F16;
    use fixed_analytics::{
        acos, acosh, asin, asinh, atan, atanh, cos, cosh, exp, ln, sin, sinh, tan, tanh,
    };

    const TOL: f64 = 0.01;

    #[test]
    fn sin_asin_roundtrip() {
        // For x in [-1, 1]: sin(asin(x)) ≈ x
        for i in 0..200 {
            let t = (i as f64) / 199.0 * 1.98 - 0.99; // (-0.99, 0.99)
            let x = I16F16::from_num(t);

            if let Ok(asin_x) = asin(x) {
                let roundtrip = sin(asin_x);
                let err = (roundtrip.to_num::<f64>() - t).abs();
                assert!(err < TOL, "sin(asin({t})) = {}", roundtrip.to_num::<f64>());
            }
        }
    }

    #[test]
    fn asin_sin_roundtrip() {
        // For x in [-π/2, π/2]: asin(sin(x)) ≈ x
        let half_pi = core::f64::consts::FRAC_PI_2 * 0.95;
        for i in 0..200 {
            let t = (i as f64) / 199.0 * 2.0 * half_pi - half_pi;
            let x = I16F16::from_num(t);
            let sin_x = sin(x);

            if let Ok(roundtrip) = asin(sin_x) {
                let err = (roundtrip.to_num::<f64>() - t).abs();
                assert!(err < TOL, "asin(sin({t})) = {}", roundtrip.to_num::<f64>());
            }
        }
    }

    #[test]
    fn cos_acos_roundtrip() {
        // For x in [-1, 1]: cos(acos(x)) ≈ x
        for i in 0..200 {
            let t = (i as f64) / 199.0 * 1.98 - 0.99;
            let x = I16F16::from_num(t);

            if let Ok(acos_x) = acos(x) {
                let roundtrip = cos(acos_x);
                let err = (roundtrip.to_num::<f64>() - t).abs();
                assert!(err < TOL, "cos(acos({t})) = {}", roundtrip.to_num::<f64>());
            }
        }
    }

    #[test]
    fn tan_atan_roundtrip() {
        // atan(tan(x)) ≈ x for x in (-π/2, π/2)
        let half_pi = core::f64::consts::FRAC_PI_2 * 0.9;
        for i in 0..200 {
            let t = (i as f64) / 199.0 * 2.0 * half_pi - half_pi;
            let x = I16F16::from_num(t);
            let tan_x = tan(x);
            let roundtrip = atan(tan_x);

            let err = (roundtrip.to_num::<f64>() - t).abs();
            assert!(err < TOL, "atan(tan({t})) = {}", roundtrip.to_num::<f64>());
        }
    }

    #[test]
    fn exp_ln_roundtrip() {
        // For x > 0: exp(ln(x)) ≈ x
        for i in 1..200 {
            let t = (i as f64) / 20.0; // (0.05, 10)
            let x = I16F16::from_num(t);

            if let Ok(ln_x) = ln(x) {
                let roundtrip = exp(ln_x);
                let rel_err = (roundtrip.to_num::<f64>() - t).abs() / t;
                assert!(
                    rel_err < TOL,
                    "exp(ln({t})) = {}",
                    roundtrip.to_num::<f64>()
                );
            }
        }
    }

    #[test]
    fn ln_exp_roundtrip() {
        // ln(exp(x)) ≈ x for small x
        for i in 0..200 {
            let t = (i as f64) / 20.0 - 5.0; // (-5, 5)
            let x = I16F16::from_num(t);
            let exp_x = exp(x);

            if let Ok(roundtrip) = ln(exp_x) {
                let err = (roundtrip.to_num::<f64>() - t).abs();
                assert!(err < TOL, "ln(exp({t})) = {}", roundtrip.to_num::<f64>());
            }
        }
    }

    #[test]
    fn sinh_asinh_roundtrip() {
        // Hyperbolic functions have lower precision, use smaller range
        for i in 0..200 {
            let t = (i as f64) / 50.0 - 2.0; // (-2, 2)
            let x = I16F16::from_num(t);
            let sinh_x = sinh(x);
            let roundtrip = asinh(sinh_x);

            let err = (roundtrip.to_num::<f64>() - t).abs();
            assert!(
                err < TOL * 2.0,
                "asinh(sinh({t})) = {}",
                roundtrip.to_num::<f64>()
            );
        }
    }

    #[test]
    fn cosh_acosh_roundtrip() {
        // For x >= 0: acosh(cosh(x)) ≈ x (cosh is not injective, so test non-negative)
        // Use smaller range for better precision
        for i in 0..200 {
            let t = (i as f64) / 50.0; // (0, 4)
            let x = I16F16::from_num(t);
            let cosh_x = cosh(x);

            if let Ok(roundtrip) = acosh(cosh_x) {
                let err = (roundtrip.to_num::<f64>() - t).abs();
                assert!(
                    err < TOL * 2.0,
                    "acosh(cosh({t})) = {}",
                    roundtrip.to_num::<f64>()
                );
            }
        }
    }

    #[test]
    fn tanh_atanh_roundtrip() {
        // For x in (-1, 1): tanh(atanh(x)) ≈ x
        for i in 0..200 {
            let t = (i as f64) / 199.0 * 1.9 - 0.95; // (-0.95, 0.95)
            let x = I16F16::from_num(t);

            if let Ok(atanh_x) = atanh(x) {
                let roundtrip = tanh(atanh_x);
                let err = (roundtrip.to_num::<f64>() - t).abs();
                assert!(
                    err < TOL,
                    "tanh(atanh({t})) = {}",
                    roundtrip.to_num::<f64>()
                );
            }
        }
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used, reason = "test code uses unwrap for conciseness")]
mod monotonicity {
    //! Verify that monotonic functions are actually monotonic.

    use fixed::types::I16F16;
    use fixed_analytics::{asin, atan, exp, ln, sin, sqrt, tanh};

    #[test]
    fn sqrt_is_increasing() {
        let mut prev = I16F16::ZERO;
        for i in 0..1000 {
            let x = I16F16::from_bits(i * 1000);
            let y = sqrt(x).unwrap();
            assert!(
                y >= prev,
                "sqrt({}) = {} < sqrt({}) = {}",
                x.to_num::<f64>(),
                y.to_num::<f64>(),
                I16F16::from_bits((i.saturating_sub(1)) * 1000).to_num::<f64>(),
                prev.to_num::<f64>()
            );
            prev = y;
        }
    }

    #[test]
    fn exp_is_increasing() {
        let mut prev = I16F16::ZERO;
        for i in 0..500 {
            // Range from -5 to 5
            let t = (i as f64) / 50.0 - 5.0;
            let x = I16F16::from_num(t);
            let y = exp(x);
            assert!(
                y >= prev,
                "exp({}) = {} should be >= {}",
                t,
                y.to_num::<f64>(),
                prev.to_num::<f64>()
            );
            prev = y;
        }
    }

    #[test]
    fn ln_is_increasing() {
        let mut prev = I16F16::MIN;
        for i in 1..500 {
            let x = I16F16::from_bits(i * 500);
            if let Ok(y) = ln(x) {
                assert!(
                    y >= prev,
                    "ln({}) = {} should be >= {}",
                    x.to_num::<f64>(),
                    y.to_num::<f64>(),
                    prev.to_num::<f64>()
                );
                prev = y;
            }
        }
    }

    #[test]
    fn sin_increasing_on_neg_half_pi_to_half_pi() {
        // sin is increasing on [-π/2, π/2]
        let half_pi = core::f64::consts::FRAC_PI_2;
        let mut prev = I16F16::MIN;
        for i in 0..200 {
            let t = (i as f64) / 199.0 * 2.0 * half_pi - half_pi;
            let x = I16F16::from_num(t);
            let y = sin(x);
            assert!(
                y >= prev,
                "sin({t}) = {} should be >= {}",
                y.to_num::<f64>(),
                prev.to_num::<f64>()
            );
            prev = y;
        }
    }

    #[test]
    fn atan_is_increasing() {
        let mut prev = I16F16::MIN;
        for i in 0..500 {
            // Range from -100 to 100
            let t = (i as f64) / 2.5 - 100.0;
            let x = I16F16::from_num(t);
            let y = atan(x);
            assert!(
                y >= prev,
                "atan({t}) = {} should be >= {}",
                y.to_num::<f64>(),
                prev.to_num::<f64>()
            );
            prev = y;
        }
    }

    #[test]
    fn asin_is_increasing() {
        let mut prev = I16F16::MIN;
        for i in 0..200 {
            let t = (i as f64) / 199.0 * 1.98 - 0.99; // (-0.99, 0.99)
            let x = I16F16::from_num(t);
            if let Ok(y) = asin(x) {
                assert!(
                    y >= prev,
                    "asin({t}) = {} should be >= {}",
                    y.to_num::<f64>(),
                    prev.to_num::<f64>()
                );
                prev = y;
            }
        }
    }

    #[test]
    fn tanh_is_increasing() {
        let mut prev = I16F16::MIN;
        for i in 0..500 {
            // Use smaller range to avoid saturation at ±1
            let t = (i as f64) / 50.0 - 5.0; // (-5, 5)
            let x = I16F16::from_num(t);
            let y = tanh(x);
            assert!(
                y >= prev,
                "tanh({t}) = {} should be >= {}",
                y.to_num::<f64>(),
                prev.to_num::<f64>()
            );
            prev = y;
        }
    }
}

#[cfg(test)]
mod bounds {
    //! Verify output bounds for functions with known ranges.

    use fixed::types::I16F16;
    use fixed_analytics::{atan, atan2, cos, exp, sin, sqrt, tanh};

    fn sample_bits(seed: u64, index: u64) -> i32 {
        let mut x = seed.wrapping_add(index.wrapping_mul(0x9E37_79B9_7F4A_7C15));
        x = x.wrapping_mul(0x5851_F42D_4C95_7F2D);
        x = x ^ (x >> 32);
        x as i32
    }

    const SAMPLES: u64 = 2000;
    const SEED: u64 = 0xBEEF_CAFE_DEAD_F00D;

    #[test]
    fn sin_in_bounds() {
        // sin(x) ∈ [-1, 1]
        for i in 0..SAMPLES {
            let bits = sample_bits(SEED, i);
            let x = I16F16::from_bits(bits);
            let y = sin(x);
            let y_f64: f64 = y.to_num();

            assert!(
                y_f64 >= -1.01 && y_f64 <= 1.01,
                "sin({}) = {} out of bounds",
                x.to_num::<f64>(),
                y_f64
            );
        }
    }

    #[test]
    fn cos_in_bounds() {
        // cos(x) ∈ [-1, 1]
        for i in 0..SAMPLES {
            let bits = sample_bits(SEED, i);
            let x = I16F16::from_bits(bits);
            let y = cos(x);
            let y_f64: f64 = y.to_num();

            assert!(
                y_f64 >= -1.01 && y_f64 <= 1.01,
                "cos({}) = {} out of bounds",
                x.to_num::<f64>(),
                y_f64
            );
        }
    }

    #[test]
    fn tanh_in_bounds() {
        // tanh(x) ∈ (-1, 1)
        for i in 0..SAMPLES {
            let bits = sample_bits(SEED, i);
            let x = I16F16::from_bits(bits);
            let y = tanh(x);
            let y_f64: f64 = y.to_num();

            assert!(
                y_f64 > -1.01 && y_f64 < 1.01,
                "tanh({}) = {} out of bounds",
                x.to_num::<f64>(),
                y_f64
            );
        }
    }

    #[test]
    fn sqrt_non_negative() {
        for i in 0..SAMPLES {
            let bits = sample_bits(SEED, i);
            let x = I16F16::from_bits(bits);

            // sqrt returns Err for negative inputs, Ok for non-negative
            match sqrt(x) {
                Ok(y) => {
                    assert!(
                        y >= I16F16::ZERO,
                        "sqrt({}) = {} should be non-negative",
                        x.to_num::<f64>(),
                        y.to_num::<f64>()
                    );
                }
                Err(_) => {
                    assert!(
                        x < I16F16::ZERO,
                        "sqrt({}) returned Err but input is non-negative",
                        x.to_num::<f64>()
                    );
                }
            }
        }
    }

    #[test]
    fn exp_positive() {
        for i in 0..SAMPLES {
            let bits = sample_bits(SEED, i);
            let x = I16F16::from_bits(bits);
            let x_f64: f64 = x.to_num();

            // Skip extreme values
            if x_f64 > 10.0 || x_f64 < -20.0 {
                continue;
            }

            let y = exp(x);
            assert!(
                y > I16F16::ZERO,
                "exp({}) = {} should be positive",
                x_f64,
                y.to_num::<f64>()
            );
        }
    }

    #[test]
    fn atan_in_bounds() {
        // atan(x) ∈ (-π/2, π/2)
        let half_pi = core::f64::consts::FRAC_PI_2;
        for i in 0..SAMPLES {
            let bits = sample_bits(SEED, i);
            let x = I16F16::from_bits(bits);
            let y = atan(x);
            let y_f64: f64 = y.to_num();

            assert!(
                y_f64 > -half_pi - 0.01 && y_f64 < half_pi + 0.01,
                "atan({}) = {} out of bounds",
                x.to_num::<f64>(),
                y_f64
            );
        }
    }

    #[test]
    fn atan2_in_bounds() {
        // atan2(y, x) ∈ [-π, π]
        // Allow slightly wider tolerance for extreme input values
        let pi = core::f64::consts::PI;
        for i in 0..SAMPLES {
            let y_bits = sample_bits(SEED, i);
            let x_bits = sample_bits(SEED ^ 0xFFFF, i);
            let y = I16F16::from_bits(y_bits);
            let x = I16F16::from_bits(x_bits);

            // Skip extreme values near I16F16::MAX that cause precision issues
            if y.to_num::<f64>().abs() > 30000.0 || x.to_num::<f64>().abs() > 30000.0 {
                continue;
            }

            let result = atan2(y, x);
            let r_f64: f64 = result.to_num();

            assert!(
                r_f64 >= -pi - 0.15 && r_f64 <= pi + 0.15,
                "atan2({}, {}) = {} out of bounds",
                y.to_num::<f64>(),
                x.to_num::<f64>(),
                r_f64
            );
        }
    }
}
