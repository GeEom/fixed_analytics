#[cfg(test)]
mod tests {
    use crate::unit::support::Lcg;
    use fixed::types::{I16F16, I32F32, I64F64};
    use fixed_analytics::Rounding::{Down, Nearest, Up};
    use fixed_analytics::{
        Rounding, div_or_saturate, div_rounded, dot, dot_rounded, mul_add, mul_add_rounded,
        mul_div, mul_div_rounded, mul_rounded, sqrt_or_zero,
    };

    /// `num/den` rounded by `mode` to a raw `I32F32`, saturating, for `den > 0`.
    fn exact(num: i128, den: i128, mode: Rounding) -> I32F32 {
        let floor = num.div_euclid(den);
        let rem = num.rem_euclid(den);
        let raw = match mode {
            Down => floor,
            Up => floor + i128::from(rem != 0),
            Nearest => floor + i128::from(2 * rem > den || (2 * rem == den && floor % 2 != 0)),
        };
        I32F32::from_bits(i64::try_from(raw).unwrap_or(if raw < 0 { i64::MIN } else { i64::MAX }))
    }

    #[test]
    fn mul_div_keeps_the_product_whole() {
        let big = I16F16::from_num(30_000);
        // 30000² overflows I16F16 long before the division brings it back.
        assert_eq!(mul_div(big, big, big), big);
        let third = mul_div(I64F64::ONE, I64F64::ONE, I64F64::from_num(3));
        assert_eq!(third, I64F64::from_bits(0x5555_5555_5555_5555));
    }

    #[test]
    fn mul_div_by_zero_saturates_by_sign() {
        let (two, zero) = (I32F32::from_num(2), I32F32::ZERO);
        assert_eq!(mul_div(two, two, zero), I32F32::MAX);
        assert_eq!(mul_div(-two, two, zero), I32F32::MIN);
        assert_eq!(mul_div(zero, two, zero), I32F32::ZERO);
        assert_eq!(mul_div(I32F32::MAX, I32F32::MAX, I32F32::ONE), I32F32::MAX);
    }

    #[test]
    fn mul_add_rounds_once_to_nearest() {
        let ulp = I16F16::from_bits(1);
        let half = I16F16::from_bits(0x8000);
        // ulp·½ is half an ulp: ties to the even neighbour, 0 then 2 ulps.
        assert_eq!(mul_add(ulp, half, I16F16::ZERO), I16F16::ZERO);
        assert_eq!(mul_add(ulp, half, ulp), I16F16::from_bits(2));
        // 182² = 33124 is past MAX; the sum is not.
        let big = I16F16::from_num(182);
        assert_eq!(mul_add(big, big, I16F16::MIN), I16F16::from_num(356));
    }

    #[test]
    fn div_or_saturate_never_panics() {
        let one = I16F16::ONE;
        assert_eq!(div_or_saturate(one, I16F16::ZERO), I16F16::MAX);
        assert_eq!(div_or_saturate(-one, I16F16::ZERO), I16F16::MIN);
        assert_eq!(div_or_saturate(I16F16::ZERO, I16F16::ZERO), I16F16::ZERO);
        assert_eq!(
            div_or_saturate(one, I16F16::from_num(3)),
            I16F16::from_bits(0x5555)
        );
        assert_eq!(
            div_or_saturate(I16F16::from_num(2), I16F16::from_num(3)),
            I16F16::from_bits(0xAAAB)
        );
    }

    #[test]
    fn dot_sums_exactly_before_rounding() {
        let big = I16F16::from_num(200);
        let a = [big, big, big];
        let b = [big, -big, I16F16::from_num(0.25)];
        assert_eq!(dot(&a, &b), I16F16::from_num(50));
        assert_eq!(dot(&a, &b[..1]), I16F16::MAX);
        assert_eq!(dot::<I16F16>(&[], &[]), I16F16::ZERO);
        let m = [I64F64::MIN; 8];
        let n = [
            I64F64::MIN,
            I64F64::MIN,
            I64F64::MIN,
            I64F64::MIN,
            -I64F64::MAX,
            -I64F64::MAX,
            -I64F64::MAX,
            -I64F64::MAX,
        ];
        // Eight products of about 2^126 sum past 2^256 before cancelling.
        assert_eq!(dot(&m[..4], &n[..4]), I64F64::MAX);
    }

    #[test]
    fn sqrt_or_zero_is_total() {
        assert_eq!(sqrt_or_zero(I16F16::from_num(-4)), I16F16::ZERO);
        assert_eq!(sqrt_or_zero(I16F16::ZERO), I16F16::ZERO);
        assert_eq!(sqrt_or_zero(I16F16::from_num(9)), I16F16::from_num(3));
    }

    #[test]
    fn directed_forms_round_the_exact_value_their_way() {
        let mut rng = Lcg(0xD1);
        let one = 1i128 << 32;
        for _ in 0..20_000 {
            let scale = |g: &mut Lcg| g.next_u64().cast_signed() >> (g.next_u64() % 40);
            let (a, b, c) = (scale(&mut rng), scale(&mut rng), scale(&mut rng));
            let (fa, fb, fc) = (
                I32F32::from_bits(a),
                I32F32::from_bits(b),
                I32F32::from_bits(c),
            );
            let (ra, rb, rc) = (i128::from(a), i128::from(b), i128::from(c));
            for mode in [Down, Nearest, Up] {
                assert_eq!(mul_rounded(fa, fb, mode), exact(ra * rb, one, mode));
                assert_eq!(
                    mul_add_rounded(fa, fb, fc, mode),
                    exact(ra * rb + rc * one, one, mode)
                );
                if c != 0 {
                    let (num, den) = if c < 0 {
                        (-ra * rb, -rc)
                    } else {
                        (ra * rb, rc)
                    };
                    assert_eq!(mul_div_rounded(fa, fb, fc, mode), exact(num, den, mode));
                    let (scaled, by) = if c < 0 {
                        (-ra * one, -rc)
                    } else {
                        (ra * one, rc)
                    };
                    assert_eq!(div_rounded(fa, fc, mode), exact(scaled, by, mode));
                }
            }
            assert_eq!(mul_div_rounded(fa, fb, fc, Nearest), mul_div(fa, fb, fc));
            assert_eq!(div_rounded(fa, fc, Nearest), div_or_saturate(fa, fc));
            assert_eq!(mul_add_rounded(fa, fb, fc, Nearest), mul_add(fa, fb, fc));
        }
    }

    #[test]
    fn directed_dot_rounds_its_exact_sum() {
        let mut rng = Lcg(0xD2);
        for _ in 0..5000 {
            let v: [i32; 6] = core::array::from_fn(|_| {
                i32::try_from(rng.next_u64().cast_signed() >> 32).unwrap_or(0)
            });
            let [a0, a1, a2, b0, b1, b2] = v;
            let (fa, fb) = (
                [a0, a1, a2].map(I16F16::from_bits),
                [b0, b1, b2].map(I16F16::from_bits),
            );
            let wide = |x: i32, y: i32| i128::from(x) * i128::from(y);
            let sum = wide(a0, b0) + wide(a1, b1) + wide(a2, b2);
            let floor = sum.div_euclid(1 << 16);
            let clamp = |r: i128| {
                I16F16::from_bits(i32::try_from(r).unwrap_or(if r < 0 {
                    i32::MIN
                } else {
                    i32::MAX
                }))
            };
            let (down, up) = (dot_rounded(&fa, &fb, Down), dot_rounded(&fa, &fb, Up));
            assert_eq!(down, clamp(floor));
            assert_eq!(up, clamp(floor + i128::from(sum.rem_euclid(1 << 16) != 0)));
            assert_eq!(dot_rounded(&fa, &fb, Nearest), dot(&fa, &fb));
        }
    }

    #[test]
    fn directed_forms_bracket_on_128_bits() {
        let ulp = I64F64::DELTA;
        let half = I64F64::from_num(0.5);
        assert_eq!(
            (mul_rounded(ulp, half, Down), mul_rounded(ulp, half, Up)),
            (I64F64::ZERO, ulp)
        );
        assert_eq!(
            (mul_rounded(-ulp, half, Down), mul_rounded(-ulp, half, Up)),
            (-ulp, I64F64::ZERO)
        );
        let three = I64F64::from_num(3);
        assert_eq!(
            div_rounded(I64F64::ONE, three, Down),
            I64F64::from_bits(0x5555_5555_5555_5555)
        );
        assert_eq!(
            div_rounded(I64F64::ONE, three, Up),
            I64F64::from_bits(0x5555_5555_5555_5556)
        );
        assert_eq!(
            div_rounded(-I64F64::ONE, three, Down),
            -I64F64::from_bits(0x5555_5555_5555_5556)
        );
        assert_eq!(mul_rounded(I64F64::MAX, three, Down), I64F64::MAX);
        assert_eq!(div_rounded(-I64F64::ONE, I64F64::ZERO, Up), I64F64::MIN);
        assert_eq!(
            mul_div_rounded(I64F64::ZERO, three, I64F64::ZERO, Up),
            I64F64::ZERO
        );
        let mut rng = Lcg(0xD3);
        for _ in 0..5000 {
            let pick = |g: &mut Lcg| {
                I64F64::from_bits(
                    (i128::from(g.next_u64().cast_signed()) << 64 | i128::from(g.next_u64()))
                        >> (g.next_u64() % 100),
                )
            };
            let (a, b, c) = (pick(&mut rng), pick(&mut rng), pick(&mut rng));
            let forms: [&dyn Fn(Rounding) -> I64F64; 4] = [
                &|m| mul_rounded(a, b, m),
                &|m| mul_div_rounded(a, b, c, m),
                &|m| mul_add_rounded(a, b, c, m),
                &|m| dot_rounded(&[a, b], &[c, a], m),
            ];
            for f in forms {
                let (down, near, up) = (f(Down), f(Nearest), f(Up));
                assert!(
                    down <= near && near <= up && up.saturating_sub(down) <= ulp,
                    "{a} {b} {c}"
                );
            }
        }
    }
}
