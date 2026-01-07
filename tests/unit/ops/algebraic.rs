//! Tests for algebraic functions (sqrt)

#![allow(clippy::unwrap_used)]

#[cfg(test)]
mod tests {
    use fixed::types::I16F16;
    use fixed_analytics::sqrt;

    const TOLERANCE: f32 = 0.02;

    fn approx_eq(a: I16F16, b: f32) -> bool {
        (a.to_num::<f32>() - b).abs() < TOLERANCE
    }

    #[test]
    fn sqrt_perfect_squares() {
        assert!(approx_eq(sqrt(I16F16::from_num(0.0)).unwrap(), 0.0));
        assert!(approx_eq(sqrt(I16F16::from_num(1.0)).unwrap(), 1.0));
        assert!(approx_eq(sqrt(I16F16::from_num(4.0)).unwrap(), 2.0));
        assert!(approx_eq(sqrt(I16F16::from_num(9.0)).unwrap(), 3.0));
        assert!(approx_eq(sqrt(I16F16::from_num(16.0)).unwrap(), 4.0));
        assert!(approx_eq(sqrt(I16F16::from_num(25.0)).unwrap(), 5.0));
    }

    #[test]
    fn sqrt_common_values() {
        assert!(approx_eq(
            sqrt(I16F16::from_num(2.0)).unwrap(),
            core::f32::consts::SQRT_2
        ));
        assert!(approx_eq(sqrt(I16F16::from_num(3.0)).unwrap(), 1.7321));
        assert!(approx_eq(
            sqrt(I16F16::from_num(0.5)).unwrap(),
            core::f32::consts::FRAC_1_SQRT_2
        ));
        assert!(approx_eq(sqrt(I16F16::from_num(0.25)).unwrap(), 0.5));
    }

    #[test]
    fn sqrt_negative_returns_error() {
        assert!(sqrt(I16F16::from_num(-1.0)).is_err());
        assert!(sqrt(I16F16::from_num(-100.0)).is_err());
    }

    #[test]
    fn sqrt_squared_gives_original() {
        for i in 1..20 {
            let x = I16F16::from_num(i) * I16F16::from_num(0.5);
            let root = sqrt(x).unwrap();
            let squared: f32 = (root * root).to_num();
            let original: f32 = x.to_num();
            assert!(
                (squared - original).abs() < 0.1,
                "sqrt({original})² = {squared}, expected {original}"
            );
        }
    }
}
