#[cfg(test)]
mod tests {
    use fixed::types::{I16F16, I32F32, I64F64};
    use fixed_analytics::{div_or_saturate, dot, mul_add, mul_div, sqrt_or_zero};

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
}
