#[cfg(test)]
mod tests {
    use crate::unit::support::TestResult;
    use fixed::types::I16F16;
    use fixed_analytics::sqrt;

    const TOLERANCE: f32 = 0.02;

    fn approx_eq(a: I16F16, b: f32) -> bool {
        (a.to_num::<f32>() - b).abs() < TOLERANCE
    }

    #[test]
    fn sqrt_perfect_squares() -> TestResult {
        assert!(approx_eq(sqrt(I16F16::from_num(0.0))?, 0.0));
        assert!(approx_eq(sqrt(I16F16::from_num(1.0))?, 1.0));
        assert!(approx_eq(sqrt(I16F16::from_num(4.0))?, 2.0));
        assert!(approx_eq(sqrt(I16F16::from_num(9.0))?, 3.0));
        assert!(approx_eq(sqrt(I16F16::from_num(16.0))?, 4.0));
        assert!(approx_eq(sqrt(I16F16::from_num(25.0))?, 5.0));
        Ok(())
    }

    #[test]
    fn sqrt_common_values() -> TestResult {
        assert!(approx_eq(
            sqrt(I16F16::from_num(2.0))?,
            core::f32::consts::SQRT_2
        ));
        assert!(approx_eq(sqrt(I16F16::from_num(3.0))?, 1.7321));
        assert!(approx_eq(
            sqrt(I16F16::from_num(0.5))?,
            core::f32::consts::FRAC_1_SQRT_2
        ));
        assert!(approx_eq(sqrt(I16F16::from_num(0.25))?, 0.5));
        Ok(())
    }

    #[test]
    fn sqrt_negative_returns_error() {
        assert!(sqrt(I16F16::from_num(-1.0)).is_err());
        assert!(sqrt(I16F16::from_num(-100.0)).is_err());
    }

    #[test]
    fn sqrt_squared_gives_original() -> TestResult {
        for i in 1..20 {
            let x = I16F16::from_num(i) * I16F16::from_num(0.5);
            let root = sqrt(x)?;
            let squared: f32 = (root * root).to_num();
            let original: f32 = x.to_num();
            assert!(
                (squared - original).abs() < 0.1,
                "sqrt({original})² = {squared}, expected {original}"
            );
        }
        Ok(())
    }
}

/// Each root must be the representable value nearest the true root.
#[cfg(test)]
mod exact {
    use crate::unit::support::Lcg;
    use crate::unit::support::TestResult;
    use fixed::traits::Fixed;
    use fixed::types::{
        I3F125, I4F4, I4F60, I8F8, I8F24, I16F16, I24F8, I32F32, I48F16, I64F64, I128F0, U128F0,
    };
    use fixed_analytics::sqrt;

    type U256 = (u128, u128);

    fn wide_mul(a: u128, b: u128) -> U256 {
        const MASK: u128 = u64::MAX as u128;
        let (a_hi, a_lo) = (a >> 64, a & MASK);
        let (b_hi, b_lo) = (b >> 64, b & MASK);
        let (mid, mid_carry) = (a_lo * b_hi).overflowing_add(a_hi * b_lo);
        let (lo, lo_carry) = (a_lo * b_lo).overflowing_add(mid << 64);
        let hi = a_hi * b_hi + (mid >> 64) + (u128::from(mid_carry) << 64) + u128::from(lo_carry);
        (hi, lo)
    }

    fn wide_shl(x: u128, shift: u32) -> U256 {
        if shift == 0 {
            (0, x)
        } else {
            (x >> (128 - shift), x << shift)
        }
    }

    fn wide_sub(a: U256, b: U256) -> U256 {
        let (lo, borrow) = a.1.overflowing_sub(b.1);
        (a.0 - b.0 - u128::from(borrow), lo)
    }

    fn expected_bits<T: Fixed>(x: T) -> i128
    where
        T::Bits: Into<i128>,
    {
        let floor: i128 = x.sqrt().to_bits().into();
        let xb: i128 = x.to_bits().into();
        let n = wide_shl(xb.cast_unsigned(), T::FRAC_NBITS);
        let f = floor.cast_unsigned();
        let rem = wide_sub(n, wide_mul(f, f));
        assert_eq!(rem.0, 0, "remainder exceeds 2·root + 1");
        if rem.1 > f { floor + 1 } else { floor }
    }

    fn samples<T: Fixed>(random: usize) -> Vec<T>
    where
        T::Bits: TryFrom<i128> + Into<i128>,
    {
        let bits = T::INT_NBITS + T::FRAC_NBITS;
        let max: i128 = T::MAX.to_bits().into();
        let mut raw: Vec<i128> = vec![0, 1, 2, 3, max, max - 1];
        for e in 0..(bits - 1) {
            for d in [-2, -1, 0, 1, 2] {
                raw.push((1i128 << e) + d);
            }
        }
        // X = k² ± j puts the root near an integer: both rounding directions and the Newton fix.
        for p in 1..64u32 {
            let k = 1i128 << p;
            if let Some(k2) = k.checked_mul(k) {
                for j in [-3, -2, -1, 0, 1, 2, 3] {
                    raw.push(k2 + j);
                }
            }
        }
        let mut rng = Lcg(0xABCD_EF01 ^ u64::from(bits));
        for _ in 0..random {
            let r = (u128::from(rng.next_u64()) | (u128::from(rng.next_u64()) << 64)).cast_signed();
            let full = (r >> (128 - bits)).abs().min(max);
            raw.extend([
                full,
                full >> (bits / 2),
                full >> (bits / 3),
                full >> (2 * bits / 3),
            ]);
            let k = full >> (bits / 2 + 1);
            for j in [-1, 0, 1] {
                raw.push(k * k + k + j);
            }
        }
        raw.into_iter()
            .filter(|&v| (0..=max).contains(&v))
            .filter_map(|v| T::Bits::try_from(v).ok())
            .map(T::from_bits)
            .collect()
    }

    fn check<T>(random: usize) -> TestResult
    where
        T: Fixed + fixed_analytics::CordicNumber + core::fmt::Display,
        T::Bits: TryFrom<i128> + Into<i128>,
    {
        for x in samples::<T>(random) {
            let got: i128 = sqrt(x)?.to_bits().into();
            let want = expected_bits(x);
            assert_eq!(got, want, "sqrt({x}) raw bits: got {got}, want {want}");
        }
        Ok(())
    }

    #[test]
    fn i4f4_and_i8f8_are_correctly_rounded() -> TestResult {
        check::<I4F4>(500)?;
        check::<I8F8>(2000)?;
        Ok(())
    }

    #[test]
    fn i16f16_i8f24_i24f8_are_correctly_rounded() -> TestResult {
        check::<I16F16>(5000)?;
        check::<I8F24>(2000)?;
        check::<I24F8>(2000)?;
        Ok(())
    }

    #[test]
    fn i32f32_i4f60_i48f16_are_correctly_rounded() -> TestResult {
        check::<I32F32>(5000)?;
        check::<I4F60>(2000)?;
        check::<I48F16>(2000)?;
        Ok(())
    }

    #[test]
    fn i64f64_i3f125_i128f0_are_correctly_rounded() -> TestResult {
        // 128-bit types take the Newton path above 2^(128 - 2·frac).
        check::<I64F64>(5000)?;
        check::<I3F125>(3000)?;
        check::<I128F0>(1000)?;
        Ok(())
    }

    /// Below 0.25 on I64F64 the exact root fits in `u128` arithmetic.
    #[test]
    fn i64f64_sub_unit_matches_integer_root() -> TestResult {
        let mut rng = Lcg(0x5EED);
        let mut raw: Vec<u128> = vec![1, 2, 3, u64::MAX.into()];
        for i in 0..4000 {
            let exponent = 62.0 * f64::from(i) / 4000.0;
            raw.push(U128F0::saturating_from_num(exponent.exp2().max(1.0).trunc()).to_bits());
            raw.push(rng.next_u64().into());
            raw.push(u128::from(rng.next_u64()) >> (rng.next_u64() % 62));
        }
        for x_raw in raw {
            let x_raw = x_raw.min((1u128 << 62) - 1);
            let x = I64F64::from_bits(x_raw.cast_signed());
            // round(√(X·2^64)) = (isqrt(X·2^66) + 1) >> 1
            let want = (((x_raw << 66).isqrt() + 1) >> 1).cast_signed();
            let got = sqrt(x)?.to_bits();
            assert_eq!(got, want, "sqrt({x}) raw bits: got {got}, want {want}");
        }
        Ok(())
    }

    /// Just below a perfect square, the Newton step lands one above the floor.
    #[test]
    fn i64f64_newton_correction_cases() -> TestResult {
        for p in 33..63u32 {
            for j in 1..=3i128 {
                let x = I64F64::from_bits((1i128 << (2 * p - 64)) - j);
                let got: i128 = sqrt(x)?.to_bits();
                let want = expected_bits(x);
                assert_eq!(got, want, "sqrt({x}) raw bits: got {got}, want {want}");
            }
        }
        Ok(())
    }

    #[test]
    fn special_values() -> TestResult {
        assert_eq!(sqrt(I64F64::ZERO)?, I64F64::ZERO);
        assert_eq!(sqrt(I64F64::ONE)?, I64F64::ONE);
        assert_eq!(sqrt(I64F64::from_num(4))?, I64F64::from_num(2));
        assert_eq!(sqrt(I64F64::from_num(0.25))?, I64F64::from_num(0.5));
        assert_eq!(sqrt(I64F64::from_bits(1))?, I64F64::from_bits(1 << 32));
        assert_eq!(sqrt(I64F64::MAX)?.to_bits(), expected_bits(I64F64::MAX));
        assert_eq!(sqrt(I3F125::MAX)?.to_bits(), expected_bits(I3F125::MAX));
        assert_eq!(
            i128::from(sqrt(I16F16::MAX)?.to_bits()),
            expected_bits(I16F16::MAX)
        );
        Ok(())
    }

    #[test]
    fn asin_acos_still_saturate_near_one() -> TestResult {
        let x = I64F64::ONE - I64F64::from_bits(1);
        assert!((fixed_analytics::asin(x)? - I64F64::FRAC_PI_2).abs() < I64F64::from_num(1e-9));
        Ok(())
    }
}
