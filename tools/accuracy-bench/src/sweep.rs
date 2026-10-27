//! The worst-case sweep: each function's largest error in ulps on every
//! tested type, from exact inputs against MPFR, with symmetry and
//! monotonicity checks on the same inputs.
//!
//! Types of 16 bits or fewer are swept exhaustively (8 bits for functions of
//! two arguments); wider ones on a seeded sample that crowds round the points
//! where each function is hardest.

use crate::sut::{PREC, Sut};
use fixed_analytics::CordicNumber;
use rug::Float;
use rug::float::Constant;
use rug::ops::Pow;

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Fun {
    Sin,
    Cos,
    Tan,
    Asin,
    Acos,
    Atan,
    Atan2,
    Sinh,
    Cosh,
    Tanh,
    Coth,
    Asinh,
    Acosh,
    Atanh,
    Acoth,
    Exp,
    Pow2,
    Ln,
    Log2,
    Log10,
    Pow,
    Sqrt,
    SqrtOrZero,
    MulAdd,
    MulDiv,
    DivOrSaturate,
    Dot3,
    Hypot,
    Norm3,
    Normalize3,
    ScaleBy,
    SqrtProduct,
    AsinFromGap,
    Rsqrt,
    Expm1,
    Exprel,
    Exprel2,
    MomentE1,
    MomentE2,
    MomentE3,
    Log1p,
    Log1pRatio,
    Log1pResidual,
    Log1pResidualRatio,
    StumpffC0,
    StumpffC1,
    StumpffC2,
    StumpffC3,
    TanRatio,
    AtanRatio,
    AtanhRatio,
    TanhRatio,
}

use Fun::*;

pub const ALL: &[Fun] = &[
    Sin, Cos, Tan, Asin, Acos, Atan, Atan2, Sinh, Cosh, Tanh, Coth, Asinh, Acosh, Atanh, Acoth,
    Exp, Pow2, Ln, Log2, Log10, Pow, Sqrt, SqrtOrZero, MulAdd, MulDiv, DivOrSaturate, Dot3,
    Hypot, Norm3, Normalize3, ScaleBy, SqrtProduct, AsinFromGap, Rsqrt,
    Expm1, Exprel, Exprel2, MomentE1, MomentE2, MomentE3, Log1p, Log1pRatio, Log1pResidual, Log1pResidualRatio, StumpffC0, StumpffC1, StumpffC2, StumpffC3, TanRatio, AtanRatio, AtanhRatio, TanhRatio,
];

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Parity {
    Odd,
    Even,
    Neither,
}

/// Sampled points per function and wide type; pairs for two arguments.
const SAMPLES: usize = 40_000;
/// Pairs drawn per function of two arguments on 16-bit types.
const PAIRS_16: usize = 1 << 18;

pub const fn sample_count() -> usize {
    SAMPLES
}

/// References for the functions whose closed forms cancel near 0, worked at
/// three times the precision, which covers the cancellation for any input of
/// a 128-bit type.
fn special(fun: Fun, x: &Float) -> Float {
    let p = PREC * 3;
    let one = Float::with_val(p, 1);
    let zero = x.is_zero();
    let at_zero = |v: f64| Float::with_val(p, v);
    match fun {
        Expm1 => x.clone().exp_m1(),
        Exprel if zero => at_zero(1.0),
        Exprel => x.clone().exp_m1() / x,
        Exprel2 if zero => at_zero(0.5),
        Exprel2 => (x.clone().exp_m1() - x) / x.clone().square(),
        MomentE1 | MomentE2 | MomentE3 if zero => {
            one / match fun { MomentE1 => 1u32, MomentE2 => 2, _ => 3 }
        }
        MomentE1 => -(-x.clone()).exp_m1() / x,
        MomentE2 => (one.clone() - (-x.clone()).exp() * (one + x)) / x.clone().square(),
        MomentE3 => {
            let poly = Float::with_val(p, x * x) + Float::with_val(p, x * 2u32) + 2u32;
            (Float::with_val(p, 2) - (-x.clone()).exp() * poly) / Float::with_val(p, x * x) / x
        }
        Log1p => x.clone().ln_1p(),
        Log1pRatio if zero => at_zero(1.0),
        Log1pRatio => x.clone().ln_1p() / x,
        Log1pResidual => x.clone().ln_1p() - x,
        Log1pResidualRatio if zero => at_zero(-0.5),
        Log1pResidualRatio => (x.clone().ln_1p() - x) / x.clone().square(),
        StumpffC0 | StumpffC1 | StumpffC2 | StumpffC3 => {
            let k = match fun { StumpffC0 => 0, StumpffC1 => 1, StumpffC2 => 2, _ => 3 };
            if zero {
                return one / [1u32, 1, 2, 6][k];
            }
            let w = x.clone().abs().sqrt();
            let (c, s) = if x.is_sign_positive() {
                (w.clone().cos(), w.clone().sin())
            } else {
                (w.clone().cosh(), w.clone().sinh())
            };
            let c1 = s / &w;
            match k {
                0 => c,
                1 => c1,
                2 => (one - c) / x,
                _ => (one - c1) / x,
            }
        }
        TanRatio if zero => at_zero(1.0),
        TanRatio => x.clone().tan() / x,
        AtanRatio if zero => at_zero(1.0),
        AtanRatio => x.clone().atan() / x,
        AtanhRatio if zero => at_zero(1.0),
        AtanhRatio => x.clone().atanh() / x,
        TanhRatio if zero => at_zero(1.0),
        TanhRatio => x.clone().tanh() / x,
        _ => at_zero(f64::NAN),
    }
}

/// `num / den`, with 0 for 0/0 and an infinity of the numerator's sign for a
/// zero denominator, as the saturating forms define it.
fn ratio(num: Float, den: &Float) -> Float {
    if num.is_zero() {
        return f(0.0);
    }
    if den.is_zero() {
        return f(if num.is_sign_negative() { f64::NEG_INFINITY } else { f64::INFINITY });
    }
    num / den
}

fn f(x: f64) -> Float {
    Float::with_val(PREC, x)
}

fn pi() -> Float {
    Float::with_val(PREC, Constant::Pi)
}

impl Fun {
    pub fn name(self) -> &'static str {
        match self {
            Sin => "sin",
            Cos => "cos",
            Tan => "tan",
            Asin => "asin",
            Acos => "acos",
            Atan => "atan",
            Atan2 => "atan2",
            Sinh => "sinh",
            Cosh => "cosh",
            Tanh => "tanh",
            Coth => "coth",
            Asinh => "asinh",
            Acosh => "acosh",
            Atanh => "atanh",
            Acoth => "acoth",
            Exp => "exp",
            Pow2 => "pow2",
            Ln => "ln",
            Log2 => "log2",
            Log10 => "log10",
            Pow => "pow",
            Sqrt => "sqrt",
            SqrtOrZero => "sqrt_or_zero",
            MulAdd => "mul_add",
            MulDiv => "mul_div",
            DivOrSaturate => "div_or_saturate",
            Dot3 => "dot3",
            Hypot => "hypot",
            Norm3 => "norm3",
            Normalize3 => "normalize3",
            ScaleBy => "scale_by",
            SqrtProduct => "sqrt_product",
            AsinFromGap => "asin_from_gap",
            Rsqrt => "rsqrt",
            Expm1 => "expm1",
            Exprel => "exprel",
            Exprel2 => "exprel2",
            MomentE1 => "exp_moments_e1",
            MomentE2 => "exp_moments_e2",
            MomentE3 => "exp_moments_e3",
            Log1p => "log1p",
            Log1pRatio => "log1p_ratio",
            Log1pResidual => "log1p_residual",
            Log1pResidualRatio => "log1p_residual_ratio",
            StumpffC0 => "stumpff_c0",
            StumpffC1 => "stumpff_c1",
            StumpffC2 => "stumpff_c2",
            StumpffC3 => "stumpff_c3",
            TanRatio => "tan_ratio",
            AtanRatio => "atan_ratio",
            AtanhRatio => "atanh_ratio",
            TanhRatio => "tanh_ratio",
        }
    }

    pub fn arity(self) -> usize {
        match self {
            Atan2 | Pow | DivOrSaturate | Hypot | ScaleBy | SqrtProduct => 2,
            MulAdd | MulDiv | Norm3 | Normalize3 => 3,
            Dot3 => 6,
            _ => 1,
        }
    }

    pub fn eval<T: Sut>(self, a: &[T]) -> Option<T> {
        use fixed_analytics as fa;
        let x = a[0];
        match self {
            Sin => Some(fa::sin(x)),
            Cos => Some(fa::cos(x)),
            Tan => Some(fa::tan(x)),
            Asin => fa::asin(x).ok(),
            Acos => fa::acos(x).ok(),
            Atan => Some(fa::atan(x)),
            Atan2 => Some(fa::atan2(x, a[1])),
            Sinh => Some(fa::sinh(x)),
            Cosh => Some(fa::cosh(x)),
            Tanh => Some(fa::tanh(x)),
            Coth => fa::coth(x).ok(),
            Asinh => Some(fa::asinh(x)),
            Acosh => fa::acosh(x).ok(),
            Atanh => fa::atanh(x).ok(),
            Acoth => fa::acoth(x).ok(),
            Exp => Some(fa::exp(x)),
            Pow2 => Some(fa::pow2(x)),
            Ln => fa::ln(x).ok(),
            Log2 => fa::log2(x).ok(),
            Log10 => fa::log10(x).ok(),
            Pow => fa::pow(x, a[1]).ok(),
            Sqrt => fa::sqrt(x).ok(),
            SqrtOrZero => Some(fa::sqrt_or_zero(x)),
            MulAdd => Some(fa::mul_add(x, a[1], a[2])),
            MulDiv => Some(fa::mul_div(x, a[1], a[2])),
            DivOrSaturate => Some(fa::div_or_saturate(x, a[1])),
            Dot3 => Some(fa::dot(&a[0..3], &a[3..6])),
            Hypot => Some(fa::hypot(x, a[1])),
            Norm3 => Some(fa::norm(&a[0..3])),
            Normalize3 => fa::normalize([x, a[1], a[2]]).ok().map(|v| v[0]),
            ScaleBy => Some(fa::scale_by(x, exponent_of(a[1]))),
            SqrtProduct => fa::sqrt_product(x, a[1]).ok(),
            AsinFromGap => fa::asin_from_gap(x).ok(),
            Rsqrt => fa::rsqrt(x).ok(),
            Expm1 => Some(fa::expm1(x)),
            Exprel => Some(fa::exprel(x)),
            Exprel2 => Some(fa::exprel2(x)),
            MomentE1 => Some(fa::exp_moments(x)[0]),
            MomentE2 => Some(fa::exp_moments(x)[1]),
            MomentE3 => Some(fa::exp_moments(x)[2]),
            Log1p => fa::log1p(x).ok(),
            Log1pRatio => fa::log1p_ratio(x).ok(),
            Log1pResidual => fa::log1p_residual(x).ok(),
            Log1pResidualRatio => fa::log1p_residual_ratio(x).ok(),
            StumpffC0 => Some(fa::stumpff(x)[0]),
            StumpffC1 => Some(fa::stumpff(x)[1]),
            StumpffC2 => Some(fa::stumpff(x)[2]),
            StumpffC3 => Some(fa::stumpff(x)[3]),
            TanRatio => Some(fa::tan_ratio(x).0),
            AtanRatio => Some(fa::atan_ratio(x)),
            AtanhRatio => fa::atanh_ratio(x).ok(),
            TanhRatio => Some(fa::tanh_ratio(x)),
        }
    }

    /// Whether the fallible form should succeed at `a`.
    pub fn defined<T: Sut>(self, a: &[T]) -> bool {
        let x = a[0];
        let zero = T::zero();
        let one = T::one();
        match self {
            Asin | Acos => x >= -one && x <= one,
            Coth => x != zero,
            Acosh => x >= one,
            Atanh => x > -one && x < one,
            Acoth => x < -one || x > one,
            Ln | Log2 | Log10 | Rsqrt => x > zero,
            Log1p | Log1pRatio | Log1pResidual | Log1pResidualRatio => x > -one,
            AtanhRatio => x > -one && x < one,
            Sqrt => x >= zero,
            AsinFromGap => x >= zero && x <= one + one,
            SqrtProduct => !((x < zero && a[1] > zero) || (x > zero && a[1] < zero)),
            Normalize3 => a.iter().any(|&v| v != zero),
            ScaleBy => {
                let e = a[1];
                (e >> T::frac_n()) << T::frac_n() == e && exponent_of(e).abs() <= 2 * T::bits() as i32
            }
            Pow => {
                let e = a[1];
                if x == zero {
                    e >= zero
                } else if x < zero {
                    e == zero || (e >> T::frac_n()) << T::frac_n() == e
                } else {
                    true
                }
            }
            _ => true,
        }
    }

    /// The exact value at `PREC` bits; infinite past the range of `Float`.
    pub fn reference(self, a: &[Float]) -> Float {
        let x = a[0].clone();
        let one = f(1.0);
        match self {
            Sin => x.sin(),
            Cos => x.cos(),
            Tan => x.tan(),
            Asin => x.asin(),
            Acos => x.acos(),
            Atan => x.atan(),
            Atan2 => x.atan2(&a[1]),
            Sinh => x.sinh(),
            Cosh => x.cosh(),
            Tanh => x.tanh(),
            Coth => one / x.tanh(),
            Asinh => x.asinh(),
            Acosh => x.acosh(),
            Atanh => x.atanh(),
            Acoth => (one / x).atanh(),
            Exp => x.exp(),
            Pow2 => x.exp2(),
            Ln => x.ln(),
            Log2 => x.log2(),
            Log10 => x.log10(),
            Pow => {
                if a[1].is_zero() {
                    one
                } else {
                    x.pow(&a[1])
                }
            }
            Sqrt => x.sqrt(),
            SqrtOrZero => {
                if x.is_sign_positive() && !x.is_zero() {
                    x.sqrt()
                } else {
                    f(0.0)
                }
            }
            MulAdd => Float::with_val(PREC * 3, &x * &a[1]) + &a[2],
            MulDiv => ratio(Float::with_val(PREC * 3, &x * &a[1]), &a[2]),
            DivOrSaturate => ratio(x, &a[1]),
            Dot3 => {
                let mut acc = Float::with_val(PREC * 3, 0);
                for i in 0..3 {
                    acc += Float::with_val(PREC * 3, &a[i] * &a[i + 3]);
                }
                acc
            }
            Hypot => (Float::with_val(PREC * 3, x.clone().square()) + a[1].clone().square()).sqrt(),
            Norm3 | Normalize3 => {
                let mut acc = Float::with_val(PREC * 3, 0);
                for v in &a[0..3] {
                    acc += v.clone().square();
                }
                if self == Norm3 { acc.sqrt() } else { x / acc.sqrt() }
            }
            ScaleBy => {
                let e = a[1].to_i32_saturating().unwrap_or(0);
                Float::with_val(PREC * 3, x) * Float::with_val(PREC, Float::i_exp(1, e))
            }
            SqrtProduct => Float::with_val(PREC * 3, &x * &a[1]).abs().sqrt(),
            AsinFromGap => (one - x).asin(),
            Rsqrt => one / x.sqrt(),
            _ => special(self, &Float::with_val(PREC * 3, &x)),
        }
    }

    /// Whether the result is exact before its one rounding, so its error is at
    /// most half an ulp.
    pub fn rounds_once(self) -> bool {
        matches!(self, Sqrt | SqrtOrZero | MulAdd | MulDiv | DivOrSaturate | Dot3)
    }

    pub fn parity(self) -> Parity {
        match self {
            Sin | Tan | Asin | Atan | Sinh | Tanh | Coth | Asinh | Atanh | Acoth => Parity::Odd,
            Cos | Cosh | TanRatio | AtanRatio | AtanhRatio | TanhRatio => Parity::Even,
            _ => Parity::Neither,
        }
    }

    /// Open intervals on which the function is monotone, and whether it rises.
    pub fn monotone(self) -> Vec<(Float, Float, bool)> {
        let inf = f(f64::INFINITY);
        let ninf = f(f64::NEG_INFINITY);
        let half_pi = pi() / 2u32;
        let one = f(1.0);
        let zero = f(0.0);
        match self {
            Sin | Tan => vec![(-half_pi.clone(), half_pi, true)],
            Cos => vec![(-pi(), zero.clone(), true), (zero, pi(), false)],
            Asin => vec![(f(-1.0) - f(1e-300), one + f(1e-300), true)],
            Acos => vec![(f(-1.0) - f(1e-300), one + f(1e-300), false)],
            Atan | Sinh | Tanh | Asinh | Exp | Pow2 => vec![(ninf, inf, true)],
            Cosh => vec![(ninf, zero.clone() + f(1e-300), false), (zero - f(1e-300), inf, true)],
            Coth => vec![(ninf, zero.clone(), false), (zero, inf, false)],
            Acosh => vec![(one - f(1e-300), inf, true)],
            Atanh => vec![(f(-1.0), one, true)],
            Acoth => vec![(ninf, f(-1.0), false), (one, inf, false)],
            Ln | Log2 | Log10 => vec![(zero, inf, true)],
            Rsqrt => vec![(zero, inf, false)],
            Expm1 | Exprel | Exprel2 | Log1p => vec![(ninf, inf, true)],
            MomentE1 | MomentE2 | MomentE3 | Log1pRatio => vec![(ninf, inf, false)],
            StumpffC0 | StumpffC1 | StumpffC2 | StumpffC3 => vec![(ninf, zero.clone() + f(1e-300), false)],
            Log1pResidual => vec![(f(-1.0), zero.clone(), true), (zero, inf, false)],
            AsinFromGap => vec![(zero.clone() - f(1e-300), f(2.0) + f(1e-300), false)],
            Sqrt => vec![(zero - f(1e-300), inf, true)],
            _ => vec![],
        }
    }

    /// Points where the function is hardest, for the sampled types.
    fn hot<T: Sut>(self) -> Vec<Float> {
        let max = T::max_value().to_float();
        let int_bits = T::bits() - T::frac_n();
        let ln2 = Float::with_val(PREC, Constant::Log2);
        let min_exp = -(i64::from(T::frac_n()) + 1);
        let mut h = vec![f(1.0), f(-1.0), f(0.5), f(-0.5), f(2.0), f(-2.0)];
        match self {
            Sin | Cos | Tan => {
                for j in 0..int_bits + 1 {
                    let k = Float::with_val(PREC, 1u32) << j;
                    h.push(pi() * k.clone() / 2u32);
                    h.push(-(pi() * k.clone() / 2u32));
                    h.push(pi() * k.clone() * 3u32 / 2u32);
                    h.push(pi() * k);
                }
            }
            Exp => {
                h.push(max.clone().ln());
                h.push(ln2 * min_exp);
            }
            Pow2 => {
                h.push(max.clone().log2());
                h.push(f(min_exp as f64));
            }
            Sinh | Cosh => {
                h.push(max.clone().asinh());
                h.push(-(max.clone().asinh()));
            }
            _ => {}
        }
        h
    }
}

/// splitmix64.
pub struct Rng(u64);

impl Rng {
    pub fn new(seed: u64) -> Self {
        Self(seed)
    }

    pub fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    pub fn below(&mut self, n: usize) -> usize {
        (self.next() % n as u64) as usize
    }

    fn wide(&mut self) -> u128 {
        (u128::from(self.next()) << 64) | u128::from(self.next())
    }

    /// Uniform over the raw range of `T`.
    fn uniform<T: Sut>(&mut self) -> i128 {
        (self.wide() as i128) >> (128 - T::bits())
    }

    /// A magnitude whose bit length is uniform, with a random sign.
    fn log_uniform<T: Sut>(&mut self) -> i128 {
        let len = self.below(T::bits() as usize) as u32;
        let m = if len == 0 {
            0
        } else {
            ((self.wide() >> (128 - len)) | (1u128 << (len - 1))) as i128
        };
        if self.next() & 1 == 0 { m } else { -m }
    }
}

fn seed(fun: Fun, ty: &str) -> u64 {
    let mut h: u64 = 0xCBF2_9CE4_8422_2325;
    for b in fun.name().bytes().chain(ty.bytes()) {
        h = (h ^ u64::from(b)).wrapping_mul(0x0100_0000_01B3);
    }
    h
}

fn clamp_raw<T: Sut>(r: i128) -> i128 {
    r.clamp(T::min_raw(), T::max_raw())
}

/// Candidate arguments: the ends of the range, each hot point and its
/// neighbours at every power-of-two distance, then uniform, log-uniform and
/// small values, and values crowding the hot points.
fn candidates<T: Sut>(fun: Fun, rng: &mut Rng, n: usize) -> Vec<T> {
    let mut centers = vec![0, T::min_raw(), T::max_raw()];
    let pi_t = T::from_float(&pi()).rb();
    centers.extend([pi_t, -pi_t, pi_t / 2, -pi_t / 2]);
    centers.extend(fun.hot::<T>().iter().map(|h| T::from_float(h).rb()));
    let mut raws = Vec::with_capacity(n + centers.len() * 2 * 130);
    for &c in &centers {
        for d in -3..=3 {
            raws.push(c.saturating_add(d));
        }
        for k in 0..T::bits() - 1 {
            raws.push(c.saturating_add(1 << k));
            raws.push(c.saturating_sub(1 << k));
        }
    }
    let small = (8i128 << T::frac_n()).min(T::max_raw());
    for i in 0..n {
        let r = match i % 4 {
            0 => rng.uniform::<T>(),
            1 => rng.log_uniform::<T>(),
            2 => (rng.wide() % (2 * small as u128 + 1)) as i128 - small,
            _ => {
                let c = centers[rng.below(centers.len())];
                c.saturating_add(rng.log_uniform::<T>() >> rng.below(T::bits() as usize))
            }
        };
        raws.push(r);
    }
    raws.into_iter()
        .map(|r| T::from_rb(clamp_raw::<T>(r)))
        .collect()
}

/// The integer part of `e`, the exponent `scale_by` reads from its sample.
fn exponent_of<T: Sut>(e: T) -> i32 {
    e.to_float().to_i32_saturating().unwrap_or(0)
}

fn exponent_candidates<T: Sut>(rng: &mut Rng, n: usize) -> Vec<T> {
    let mut v: Vec<T> = (-16..=16).map(|k: i32| <T as CordicNumber>::from_num(k)).collect();
    v.extend(candidates::<T>(Pow, rng, n));
    v
}

fn all<T: Sut>() -> Vec<T> {
    (T::min_raw()..=T::max_raw()).map(T::from_rb).collect()
}

fn points<T: Sut>(fun: Fun, rng: &mut Rng) -> Vec<Vec<T>> {
    if fun.arity() == 1 {
        let xs = if T::bits() <= 16 {
            all::<T>()
        } else {
            candidates::<T>(fun, rng, SAMPLES)
        };
        return xs.into_iter().map(|x| vec![x]).collect();
    }
    if fun.arity() > 2 {
        let xs = if T::bits() <= 16 {
            all::<T>()
        } else {
            candidates::<T>(fun, rng, SAMPLES)
        };
        let n = if T::bits() <= 16 { PAIRS_16 } else { SAMPLES };
        return (0..n)
            .map(|_| (0..fun.arity()).map(|_| xs[rng.below(xs.len())]).collect())
            .collect();
    }
    if T::bits() <= 8 {
        let xs = all::<T>();
        return xs
            .iter()
            .flat_map(|&x| xs.iter().map(move |&y| vec![x, y]))
            .collect();
    }
    let n = if T::bits() <= 16 { PAIRS_16 } else { SAMPLES };
    let xs = if T::bits() <= 16 {
        all::<T>()
    } else {
        candidates::<T>(fun, rng, n)
    };
    let ys = if fun == Pow || fun == ScaleBy {
        exponent_candidates::<T>(rng, n)
    } else if T::bits() <= 16 {
        all::<T>()
    } else {
        candidates::<T>(fun, rng, n)
    };
    (0..n)
        .map(|_| vec![xs[rng.below(xs.len())], ys[rng.below(ys.len())]])
        .collect()
}

#[derive(Clone, Debug, Default)]
pub struct FunStats {
    pub fun: String,
    pub ty: String,
    pub samples: u64,
    pub max_ulps: f64,
    pub worst: String,
    pub max_rel_ulps: f64,
    pub worst_rel: String,
    pub parity_bad: u64,
    pub parity_n: u64,
    pub mono_bad: u64,
    pub mono_n: u64,
    pub errors: u64,
    pub example: String,
}

fn show<T: Sut>(a: &[T]) -> String {
    a.iter()
        .map(|x| format!("{:e}", x.to_float().to_f64()))
        .collect::<Vec<_>>()
        .join(",")
}

/// Every check on one function and type.
pub fn run<T: Sut>(fun: Fun) -> FunStats {
    let name = T::name();
    let mut rng = Rng::new(seed(fun, &name));
    let pts = points::<T>(fun, &mut rng);
    let mut s = FunStats {
        fun: fun.name().to_string(),
        ty: name,
        ..FunStats::default()
    };
    let lo = T::min_value().to_float();
    let hi = T::max_value().to_float();

    for p in &pts {
        if !fun.defined(p) {
            continue;
        }
        s.samples += 1;
        let Some(c) = fun.eval(p) else {
            s.errors += 1;
            if s.example.is_empty() {
                s.example = format!("error at {}", show(p));
            }
            continue;
        };
        let args: Vec<Float> = p.iter().map(|x| x.to_float()).collect();
        let r = fun.reference(&args);
        if r.is_nan() {
            continue;
        }
        let expected = if r < lo {
            lo.clone()
        } else if r > hi {
            hi.clone()
        } else {
            r
        };
        let err = Float::with_val(PREC, c.to_float() - &expected).abs() << T::frac_n();
        let ulps = err.to_f64();
        let scale = expected.abs().to_f64().max(1.0);
        if ulps > s.max_ulps {
            s.max_ulps = ulps;
            s.worst = show(p);
        }
        if ulps / scale > s.max_rel_ulps {
            s.max_rel_ulps = ulps / scale;
            s.worst_rel = show(p);
        }
    }

    if fun.arity() == 1 {
        parity::<T>(fun, &pts, &mut s);
        monotone::<T>(fun, &pts, &mut s);
    }
    s
}

fn parity<T: Sut>(fun: Fun, pts: &[Vec<T>], s: &mut FunStats) {
    let kind = fun.parity();
    if kind == Parity::Neither {
        return;
    }
    for p in pts {
        let x = p[0];
        let nx = x.neg_sat();
        if x == T::min_value() || !fun.defined(&[x]) || !fun.defined(&[nx]) {
            continue;
        }
        let (Some(y), Some(z)) = (fun.eval(&[x]), fun.eval(&[nx])) else {
            continue;
        };
        s.parity_n += 1;
        let want = if kind == Parity::Odd { y.neg_sat() } else { y };
        // −MAX is MIN + ulp, so an odd function saturating at MAX meets MIN.
        let saturated = kind == Parity::Odd && y == T::max_value() && z == T::min_value();
        if z != want && !saturated {
            s.parity_bad += 1;
            if s.example.is_empty() {
                s.example = format!(
                    "f({}) = {:e}, f(-x) = {:e}",
                    show(&[x]),
                    y.to_float().to_f64(),
                    z.to_float().to_f64()
                );
            }
        }
    }
}

fn monotone<T: Sut>(fun: Fun, pts: &[Vec<T>], s: &mut FunStats) {
    for (lo, hi, rising) in fun.monotone() {
        let mut raws: Vec<i128> = Vec::new();
        for p in pts {
            let r = p[0].rb();
            raws.push(r);
            if T::bits() > 16 {
                raws.push(clamp_raw::<T>(r.saturating_add(1)));
                raws.push(clamp_raw::<T>(r.saturating_add(2)));
            }
        }
        raws.sort_unstable();
        raws.dedup();
        let mut prev: Option<(T, T)> = None;
        for r in raws {
            let x = T::from_rb(r);
            let xf = x.to_float();
            if xf <= lo || xf >= hi || !fun.defined(&[x]) {
                continue;
            }
            let Some(y) = fun.eval(&[x]) else { continue };
            if let Some((px, py)) = prev {
                s.mono_n += 1;
                let bad = if rising { y < py } else { y > py };
                if bad {
                    s.mono_bad += 1;
                    if s.example.is_empty() {
                        s.example = format!(
                            "f({}) = {:e} then f({}) = {:e}",
                            show(&[px]),
                            py.to_float().to_f64(),
                            show(&[x]),
                            y.to_float().to_f64()
                        );
                    }
                }
            }
            prev = Some((x, y));
        }
    }
}
