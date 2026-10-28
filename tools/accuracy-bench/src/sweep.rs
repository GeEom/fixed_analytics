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
    MomentE4,
    MomentE5,
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
    EllipK,
    EllipKm,
    EllipE,
    JacobiSn,
    JacobiCn,
    JacobiDn,
    JacobiAm,
    JacobiInverse,
    CarlsonRf,
    CarlsonRfRoots,
    CarlsonRj,
    CarlsonRc,
    CarlsonRd,
    EllipPi,
    EllipPiComplete,
    SinTurns,
    CosTurns,
    Atan2Turns,
    DivByTau,
    WrapRadians,
    WrapPrincipal,
    MulByTau,
    CMulRe,
    CMulIm,
    CDivRe,
    CDivIm,
    CExpRe,
    CExpIm,
    CExprelRe,
    CExprelIm,
    CExprel2Re,
    CExprel2Im,
    FaddeevaRe,
    FaddeevaIm,
    CErfRe,
    CErfIm,
    E1Re,
    E1Im,
    GammaUpperRe,
    GammaUpperIm,
    CLnRe,
    CLnIm,
    CPowRe,
    CPowIm,
    WeightedK1Re,
    WeightedK1Im,
    WeightedK2Re,
    WeightedK2Im,
    KeplerUniversal,
    Erf,
    Erfc,
    FresnelC,
    FresnelS,
    Gamma,
    MulDown,
    MulUp,
    MulDivDown,
    MulDivUp,
    MulAddDown,
    MulAddUp,
    DivDown,
    DivUp,
    Dot3Down,
    Dot3Up,
    EllipPiIncrement,
}

use Fun::*;

pub const ALL: &[Fun] = &[
    Sin, Cos, Tan, Asin, Acos, Atan, Atan2, Sinh, Cosh, Tanh, Coth, Asinh, Acosh, Atanh, Acoth,
    Exp, Pow2, Ln, Log2, Log10, Pow, Sqrt, SqrtOrZero, MulAdd, MulDiv, DivOrSaturate, Dot3,
    Hypot, Norm3, Normalize3, ScaleBy, SqrtProduct, AsinFromGap, Rsqrt,
    Expm1, Exprel, Exprel2, MomentE1, MomentE2, MomentE3, MomentE4, MomentE5, Log1p, Log1pRatio, Log1pResidual, Log1pResidualRatio, StumpffC0, StumpffC1, StumpffC2, StumpffC3, TanRatio, AtanRatio, AtanhRatio, TanhRatio,
    EllipK, EllipKm, EllipE, JacobiSn, JacobiCn, JacobiDn, JacobiAm, JacobiInverse, CarlsonRf, CarlsonRfRoots, CarlsonRj, CarlsonRc, CarlsonRd, EllipPi, EllipPiComplete,
    SinTurns, CosTurns, Atan2Turns, DivByTau, WrapRadians, WrapPrincipal, MulByTau,
    CMulRe, CMulIm, CDivRe, CDivIm, CExpRe, CExpIm, CExprelRe, CExprelIm, CExprel2Re, CExprel2Im, FaddeevaRe, FaddeevaIm, CErfRe, CErfIm, E1Re, E1Im, GammaUpperRe, GammaUpperIm, CLnRe, CLnIm, CPowRe, CPowIm, WeightedK1Re, WeightedK1Im, WeightedK2Re, WeightedK2Im, KeplerUniversal, Erf, Erfc, FresnelC, FresnelS, Gamma,
    MulDown, MulUp, MulDivDown, MulDivUp, MulAddDown, MulAddUp, DivDown, DivUp, Dot3Down, Dot3Up, EllipPiIncrement,
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
/// Samples per type for the costly functions.
const HEAVY_SAMPLES: usize = 3_000;

/// CI's quick mode: an eighth of the samples, and only 8-bit types swept
/// whole, so the run fits a small runner. Local runs, which write worst.txt,
/// stay full.
pub static QUICK: std::sync::atomic::AtomicBool = std::sync::atomic::AtomicBool::new(false);

fn quick() -> bool {
    QUICK.load(std::sync::atomic::Ordering::Relaxed)
}

fn scaled(n: usize) -> usize {
    if quick() { n / 8 } else { n }
}

/// Types this narrow are swept whole for one argument.
fn exhaustive_bits() -> u32 {
    if quick() { 8 } else { 16 }
}

pub const fn sample_count() -> usize {
    SAMPLES
}

fn fun_roots(fun: Fun) -> bool {
    fun == CarlsonRfRoots
}

fn elliptic(fun: Fun, a: &[Float]) -> Float {
    use crate::elliptic_ref as er;
    let one = Float::with_val(er::P, 1);
    let up = |v: &Float| Float::with_val(er::P, v);
    match fun {
        EllipK => er::k(&a[0]),
        EllipKm => er::k(&(one - &a[0]).sqrt()),
        EllipE => er::e(&a[0]),
        JacobiSn | JacobiCn | JacobiDn | JacobiAm => {
            let (sn, cn, dn, am) = er::jacobi(&up(&a[0]), &up(&a[1]));
            match fun { JacobiSn => sn, JacobiCn => cn, JacobiDn => dn, _ => am }
        }
        JacobiInverse => er::inverse(&a[0], &a[1], &a[2], &a[3]),
        CarlsonRf => er::rf(&a[0], &a[1], &a[2]),
        CarlsonRfRoots => er::rf(&up(&a[0]).square(), &up(&a[1]).square(), &up(&a[2]).square()),
        CarlsonRj => er::rj(&a[0], &a[1], &a[2], &a[3]),
        CarlsonRc => er::rc(&up(&a[0]), &up(&a[1])),
        CarlsonRd => er::rd(&a[0], &a[1], &a[2]),
        EllipPi => er::pi(&up(&a[0]), &a[1], &a[2], &a[3], &a[4]),
        EllipPiIncrement => er::pi_increment(&up(&a[0]), &up(&a[1]), &up(&a[2]), &up(&a[3])),
        _ => er::pi_complete(&up(&a[0]), &a[1]),
    }
}

/// Points for the elliptic functions: moduli crowded toward both ends,
/// non-negative Carlson arguments, and `(sn, cn, dn)` rounded from an angle.
fn elliptic_points<T: Sut>(fun: Fun, rng: &mut Rng) -> Vec<Vec<T>> {
    let any = candidates::<T>(fun, rng, 4 * HEAVY_SAMPLES);
    let one = <T as CordicNumber>::one();
    let unit: Vec<T> = any.iter().map(|&v| { let a = CordicNumber::abs(v); if a > one { T::from_rb(a.rb() % one.rb().max(1)) } else { a } }).collect();
    let pos: Vec<T> = any.iter().map(|&v| CordicNumber::abs(v)).collect();
    let exps = exponent_candidates::<T>(rng, 4 * HEAVY_SAMPLES);
    let pick = |rng: &mut Rng, v: &[T]| v[rng.below(v.len())];
    let count = if is_weighted(fun) { HEAVY_SAMPLES / 8 } else { HEAVY_SAMPLES };
    (0..scaled(count))
        .map(|_| match fun {
            EllipK | EllipKm | EllipE => vec![pick(rng, &unit)],
            JacobiSn | JacobiCn | JacobiDn | JacobiAm => vec![pick(rng, &any), pick(rng, &unit)],
            CarlsonRf | CarlsonRd => (0..3).map(|_| pick(rng, &pos)).collect(),
            CarlsonRfRoots => (0..3).map(|_| pick(rng, &any)).collect(),
            CarlsonRj => (0..4).map(|_| pick(rng, &pos)).collect(),
            CarlsonRc => vec![pick(rng, &pos), pick(rng, &pos)],
            EllipPiComplete => vec![pick(rng, &any), pick(rng, &unit)],
            EllipPiIncrement => {
                // Characteristics mostly within (−1, 1), and a third of the steps tiny.
                let n = if rng.below(4) == 0 { pick(rng, &any) } else if rng.below(2) == 0 { pick(rng, &unit) } else { -pick(rng, &unit) };
                let n = if n >= one { n - one - one } else { n };
                let step = pick(rng, &any);
                let step = match rng.below(3) { 0 => T::from_rb(step.rb() >> rng.below(T::bits() as usize)), 1 => T::from_rb(step.rb() % (one.rb() << 3).max(1)), _ => step };
                vec![n, pick(rng, &any), step, pick(rng, &unit)]
            }
            FaddeevaRe | FaddeevaIm | CErfRe | CErfIm | E1Re | E1Im => vec![pick(rng, &any), pick(rng, &any)],
            GammaUpperRe | GammaUpperIm => vec![pick(rng, &pos), pick(rng, &any), pick(rng, &any)],
            FresnelC | FresnelS => vec![pick(rng, &any)],
            WeightedK1Re | WeightedK1Im | WeightedK2Re | WeightedK2Im => {
                vec![pick(rng, &any), pick(rng, &any), pick(rng, &unit), pick(rng, &pos)]
            }
            KeplerUniversal => {
                // A real orbit: α = (2 − w)/r₀ for w = r₀v₀²/μ ≥ 0, σ₀ = s·√(r₀w), |s| ≤ 1.
                let r0 = pick(rng, &pos).to_float();
                let w = Float::with_val(PREC, pick(rng, &pos).to_float()) * f([0.01, 0.5, 1.0, 1.9, 2.0, 3.0, 100.0][rng.below(7)]);
                let al = Float::with_val(PREC, Float::with_val(PREC, 2) - &w) / &r0;
                let s = (rng.next() as f64 / u64::MAX as f64) * 2.0 - 1.0;
                let sigma = Float::with_val(PREC, Float::with_val(PREC, &r0 * &w).sqrt() * s);
                vec![T::from_float(&r0), T::from_float(&sigma), T::from_float(&al), pick(rng, &any)]
            }
            CPowRe | CPowIm => {
                // A third of the points on an axis or diagonal, where parts vanish exactly.
                let (x, y) = (pick(rng, &any), pick(rng, &any));
                let y = match rng.below(9) { 0 => <T as CordicNumber>::zero(), 1 => x, 2 => -x, _ => y };
                let x = if rng.below(9) == 0 { <T as CordicNumber>::zero() } else { x };
                vec![x, y, pick(rng, &exps)]
            }
            _ => {
                let kc = pick(rng, &unit);
                let phi = Float::with_val(PREC, rng.next() as f64 / u64::MAX as f64 * 6.2 - 3.1);
                let m = Float::with_val(PREC, 1) - kc.to_float().square();
                let s = phi.clone().sin();
                let d = (Float::with_val(PREC, 1) - m * s.clone().square()).sqrt();
                let parts = vec![T::from_float(&s), T::from_float(&phi.cos()), T::from_float(&d), kc];
                if fun == JacobiInverse { parts } else { [vec![pick(rng, &any)], parts].concat() }
            }
        })
        .collect()
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
        MomentE1 | MomentE2 | MomentE3 | MomentE4 | MomentE5 => {
            let k = match fun { MomentE1 => 1u32, MomentE2 => 2, MomentE3 => 3, MomentE4 => 4, _ => 5 };
            // Σ (−x)ⁿ/(n!(n + k)) below 1, else (k − 1)!/x^k·(1 − e^(−x)·Σ_{j<k} x^j/j!).
            if x.clone().abs() < 1 {
                let mut term = Float::with_val(p, 1);
                let mut sum = Float::with_val(p, 1) / k;
                for n in 1..400u32 {
                    term = -term * x / n;
                    sum += Float::with_val(p, &term / (n + k));
                }
                return sum;
            }
            let mut power = Float::with_val(p, 1);
            let mut partial = Float::with_val(p, 0);
            for j in 0..k {
                partial += &power;
                power = power * x / (j + 1);
            }
            let fact = Float::with_val(p, Float::factorial(k - 1));
            fact * (one - (-x.clone()).exp() * partial) / x.clone().pow(k)
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

/// `χ` with `F(χ) = r₀χ + σ₀χ²c₂(αχ²) + (1 − αr₀)χ³c₃(αχ²) = τ`, by a doubled
/// bracket and Newton's method at three times the precision, the Stumpff
/// functions from MPFR's sines and cosines or their series below 1.
fn kepler(a: &[Float]) -> Float {
    let p = PREC * 3;
    let up = |v: &Float| Float::with_val(p, v);
    let (r0, s0, al, tau) = (up(&a[0]), up(&a[1]), up(&a[2]), up(&a[3]));
    if tau.is_zero() {
        return f(0.0);
    }
    let beta = Float::with_val(p, 1) - Float::with_val(p, &al * &r0);
    let eval = |chi: &Float| {
        let x2 = Float::with_val(p, chi.square_ref());
        let z = Float::with_val(p, &al * &x2);
        if z < -1e8 {
            // e^√|z| swamps the rest, its coefficient of χ's sign on a real orbit.
            let inf = Float::with_val(p, if chi.is_sign_negative() { f64::NEG_INFINITY } else { f64::INFINITY });
            return (inf, Float::with_val(p, f64::INFINITY));
        }
        let (c0, c1, c2, c3) = if z.clone().abs() < 1 {
            let mut c = [Float::with_val(p, 0), Float::with_val(p, 0), Float::with_val(p, 0), Float::with_val(p, 0)];
            for (k, ck) in c.iter_mut().enumerate() {
                let mut term = Float::with_val(p, 1) / Float::with_val(p, Float::factorial(k as u32));
                let mut n = k as u32;
                for _ in 0..400 {
                    *ck += &term;
                    term = -term * &z / ((n + 1) * (n + 2));
                    n += 2;
                }
            }
            let [a0, a1, a2, a3] = c;
            (a0, a1, a2, a3)
        } else {
            let w = Float::with_val(p, z.clone().abs().sqrt());
            let (c, s) = if z.is_sign_positive() { (w.clone().cos(), w.clone().sin()) } else { (w.clone().cosh(), w.clone().sinh()) };
            let c1 = Float::with_val(p, &s / &w);
            let c2 = Float::with_val(p, Float::with_val(p, 1) - &c) / &z;
            let c3 = Float::with_val(p, Float::with_val(p, 1) - &c1) / &z;
            (c, c1, c2, c3)
        };
        let fx = Float::with_val(p, &r0 * chi) + Float::with_val(p, &s0 * &x2) * &c2
            + Float::with_val(p, &beta * &x2) * chi * &c3 - &tau;
        let rx = Float::with_val(p, &r0 * &c0) + Float::with_val(p, &s0 * chi) * &c1 + Float::with_val(p, &x2 * &c2);
        (fx, rx)
    };
    let neg = tau.is_sign_negative();
    let past = |m: &Float| {
        let chi = if neg { Float::with_val(p, -m) } else { m.clone() };
        let (fx, _) = eval(&chi);
        fx.is_zero() || fx.is_sign_negative() == neg
    };
    let start = Float::with_val(p, &tau / &r0).abs();
    let (mut near, mut far) = (start.clone(), start.clone());
    if past(&start) {
        for _ in 0..1400 {
            near = Float::with_val(p, &far / 2u32);
            if !past(&near) { break; }
            far = near.clone();
        }
    } else {
        for _ in 0..1400 {
            far = Float::with_val(p, &near * 2u32);
            if past(&far) { break; }
            if far > Float::with_val(p, Float::i_exp(1, 300)) {
                return f(if neg { f64::NEG_INFINITY } else { f64::INFINITY });
            }
            near = far.clone();
        }
    }
    let (mut lo, mut hi) = if neg { (Float::with_val(p, -&far), Float::with_val(p, -&near)) } else { (near, far) };
    let mut chi = Float::with_val(p, &lo + &hi) / 2u32;
    let mut stride = Float::with_val(64, Float::with_val(p, &hi - &lo).abs());
    for _ in 0..4000 {
        let (fx, rx) = eval(&chi);
        if fx.is_zero() {
            break;
        }
        if fx.is_sign_negative() { lo = chi.clone() } else { hi = chi.clone() }
        let newton = Float::with_val(p, &chi - Float::with_val(p, &fx / &rx));
        let slow = Float::with_val(64, Float::with_val(p, &fx / &rx).abs()) * 2u32 > stride;
        let next = if slow || !rx.is_finite() || rx <= 0 || newton <= lo || newton >= hi { Float::with_val(p, &lo + &hi) / 2u32 } else { newton };
        stride = Float::with_val(64, Float::with_val(p, &next - &chi).abs());
        let step = Float::with_val(64, &next - &chi).abs();
        chi = next;
        let scale = Float::with_val(64, chi.abs_ref()).max(&Float::with_val(64, Float::i_exp(1, -400)));
        if step.is_zero() || step < scale * Float::with_val(64, Float::i_exp(1, -(p as i32) + 8)) {
            break;
        }
    }
    chi
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
            MomentE4 => "exp_moments_e4",
            MomentE5 => "exp_moments_e5",
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
            EllipK => "ellip_k",
            EllipKm => "ellip_k_from_m",
            EllipE => "ellip_e",
            JacobiSn => "jacobi_sn",
            JacobiCn => "jacobi_cn",
            JacobiDn => "jacobi_dn",
            JacobiAm => "jacobi_am",
            JacobiInverse => "jacobi_inverse",
            CarlsonRf => "carlson_rf",
            CarlsonRfRoots => "carlson_rf_roots",
            CarlsonRj => "carlson_rj",
            CarlsonRc => "carlson_rc",
            CarlsonRd => "carlson_rd",
            EllipPi => "ellip_pi",
            EllipPiComplete => "ellip_pi_complete",
            SinTurns => "sin_turns",
            CosTurns => "cos_turns",
            Atan2Turns => "atan2_turns",
            DivByTau => "div_by_tau",
            WrapRadians => "wrap_radians",
            WrapPrincipal => "wrap_principal",
            MulByTau => "mul_by_tau",
            CMulRe => "c_mul_re",
            CMulIm => "c_mul_im",
            CDivRe => "c_div_re",
            CDivIm => "c_div_im",
            CExpRe => "c_exp_re",
            CExpIm => "c_exp_im",
            CExprelRe => "c_exprel_re",
            CExprelIm => "c_exprel_im",
            CExprel2Re => "c_exprel2_re",
            CExprel2Im => "c_exprel2_im",
            FaddeevaRe => "faddeeva_re",
            FaddeevaIm => "faddeeva_im",
            CErfRe => "c_erf_re",
            CErfIm => "c_erf_im",
            E1Re => "e1_re",
            E1Im => "e1_im",
            GammaUpperRe => "gamma_upper_re",
            GammaUpperIm => "gamma_upper_im",
            CLnRe => "c_ln_re",
            CLnIm => "c_ln_im",
            CPowRe => "c_pow_re",
            CPowIm => "c_pow_im",
            WeightedK1Re => "exprel_weighted_k1_re",
            WeightedK1Im => "exprel_weighted_k1_im",
            WeightedK2Re => "exprel_weighted_k2_re",
            WeightedK2Im => "exprel_weighted_k2_im",
            KeplerUniversal => "kepler_universal",
            Erf => "erf",
            Erfc => "erfc",
            FresnelC => "fresnel_c",
            FresnelS => "fresnel_s",
            Gamma => "gamma",
            MulDown => "mul_down",
            MulUp => "mul_up",
            MulDivDown => "mul_div_down",
            MulDivUp => "mul_div_up",
            MulAddDown => "mul_add_down",
            MulAddUp => "mul_add_up",
            DivDown => "div_down",
            DivUp => "div_up",
            Dot3Down => "dot3_down",
            Dot3Up => "dot3_up",
            EllipPiIncrement => "ellip_pi_increment",
        }
    }

    pub fn arity(self) -> usize {
        match self {
            Atan2 | Pow | DivOrSaturate | Hypot | ScaleBy | SqrtProduct => 2,
            MulAdd | MulDiv | Norm3 | Normalize3 | CarlsonRf | CarlsonRfRoots | CarlsonRd => 3,
            JacobiSn | JacobiCn | JacobiDn | JacobiAm | CarlsonRc | EllipPiComplete | Atan2Turns => 2,
            JacobiInverse | CarlsonRj | CMulRe | CMulIm | CDivRe | CDivIm => 4,
            CExpRe | CExpIm | CExprelRe | CExprelIm | CExprel2Re | CExprel2Im | FaddeevaRe
            | FaddeevaIm | CErfRe | CErfIm | E1Re | E1Im => 2,
            GammaUpperRe | GammaUpperIm | CPowRe | CPowIm => 3,
            CLnRe | CLnIm => 2,
            WeightedK1Re | WeightedK1Im | WeightedK2Re | WeightedK2Im | KeplerUniversal => 4,
            EllipPi => 5,
            Dot3 | Dot3Down | Dot3Up => 6,
            MulDown | MulUp | DivDown | DivUp => 2,
            MulDivDown | MulDivUp | MulAddDown | MulAddUp => 3,
            EllipPiIncrement => 4,
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
            MomentE4 => Some(fa::exp_moments(x)[3]),
            MomentE5 => Some(fa::exp_moments(x)[4]),
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
            EllipK => fa::ellip_k(fa::Modulus::from_complement(x)?).ok(),
            EllipKm => fa::ellip_k(fa::Modulus::from_parameter(x)?).ok(),
            EllipE => Some(fa::ellip_e(fa::Modulus::from_complement(x)?)),
            JacobiSn | JacobiCn | JacobiDn => {
                let (sn, cn, dn) = fa::jacobi_sn_cn_dn(x, fa::Modulus::from_complement(a[1])?);
                Some(match self { JacobiSn => sn, JacobiCn => cn, _ => dn })
            }
            JacobiAm => Some(fa::jacobi_am(x, fa::Modulus::from_complement(a[1])?)),
            JacobiInverse => fa::jacobi_inverse(x, a[1], a[2], fa::Modulus::from_complement(a[3])?).ok(),
            CarlsonRf => fa::carlson_rf(x, a[1], a[2]).ok(),
            CarlsonRfRoots => fa::carlson_rf_roots(x, a[1], a[2]).ok(),
            CarlsonRj => fa::carlson_rj(x, a[1], a[2], a[3]).ok(),
            CarlsonRc => fa::carlson_rc(x, a[1]).ok(),
            CarlsonRd => fa::carlson_rd(x, a[1], a[2]).ok(),
            EllipPi => fa::ellip_pi(x, a[1], a[2], a[3], fa::Modulus::from_complement(a[4])?).ok(),
            EllipPiComplete => fa::ellip_pi_complete(x, fa::Modulus::from_complement(a[1])?).ok(),
            SinTurns => Some(fa::sin_cos_turns(x).0),
            CosTurns => Some(fa::sin_cos_turns(x).1),
            Atan2Turns => Some(fa::atan2_turns(x, a[1])),
            DivByTau => Some(fa::div_by_tau(x)),
            WrapRadians => Some(fa::wrap_radians(x)),
            WrapPrincipal => Some(fa::wrap_principal(x)),
            CMulRe | CMulIm | CDivRe | CDivIm => {
                let (u, v) = (fa::Complex::new(x, a[1]), fa::Complex::new(a[2], a[3]));
                let r = if matches!(self, CMulRe | CMulIm) { u * v } else { u / v };
                Some(if matches!(self, CMulRe | CDivRe) { r.re } else { r.im })
            }
            CExpRe | CExpIm | CExprelRe | CExprelIm | CExprel2Re | CExprel2Im | FaddeevaRe
            | FaddeevaIm | CErfRe | CErfIm | E1Re | E1Im => {
                let z = fa::Complex::new(x, a[1]);
                let r = match self {
                    CExpRe | CExpIm => z.exp(),
                    CExprelRe | CExprelIm => z.exprel(),
                    CExprel2Re | CExprel2Im => z.exprel2(),
                    FaddeevaRe | FaddeevaIm => z.faddeeva(),
                    CErfRe | CErfIm => z.erf(),
                    _ => z.e1().ok()?,
                };
                Some(if (self as u32) % 2 == (CExpRe as u32) % 2 { r.re } else { r.im })
            }
            GammaUpperRe | GammaUpperIm => {
                let r = fa::Complex::new(a[1], a[2]).gamma_upper(x).ok()?;
                Some(if self == GammaUpperRe { r.re } else { r.im })
            }
            CLnRe | CLnIm | CPowRe | CPowIm => {
                let z = fa::Complex::new(x, a[1]);
                let r = if matches!(self, CLnRe | CLnIm) { z.ln() } else { z.pow(a[2]) }.ok()?;
                Some(if matches!(self, CLnRe | CPowRe) { r.re } else { r.im })
            }
            WeightedK1Re | WeightedK1Im | WeightedK2Re | WeightedK2Im => {
                let (k1, k2) = fa::Complex::new(x, a[1]).exprel_weighted(a[2], a[3]).ok()?;
                Some(match self {
                    WeightedK1Re => k1.re,
                    WeightedK1Im => k1.im,
                    WeightedK2Re => k2.re,
                    _ => k2.im,
                })
            }
            KeplerUniversal => fa::kepler_universal(x, a[1], a[2], a[3]).ok(),
            MulByTau => Some(fa::mul_by_tau(x)),
            Erf => Some(fa::erf(x)),
            Erfc => Some(fa::erfc(x)),
            FresnelC => Some(fa::fresnel(x).0),
            FresnelS => Some(fa::fresnel(x).1),
            Gamma => fa::gamma(x).ok(),
            MulDown | MulUp | MulDivDown | MulDivUp | MulAddDown | MulAddUp | DivDown | DivUp
            | Dot3Down | Dot3Up => {
                let m = if self.direction() == Some(true) { fa::Rounding::Up } else { fa::Rounding::Down };
                Some(match self {
                    MulDown | MulUp => fa::mul_rounded(x, a[1], m),
                    MulDivDown | MulDivUp => fa::mul_div_rounded(x, a[1], a[2], m),
                    MulAddDown | MulAddUp => fa::mul_add_rounded(x, a[1], a[2], m),
                    DivDown | DivUp => fa::div_rounded(x, a[1], m),
                    _ => fa::dot_rounded(&a[0..3], &a[3..6], m),
                })
            }
            EllipPiIncrement => fa::ellip_pi_increment(x, a[1], a[2], fa::Modulus::from_complement(a[3])?).ok(),
        }
    }

    /// For the directed forms, whether they round up rather than down.
    pub fn direction(self) -> Option<bool> {
        match self {
            MulDown | MulDivDown | MulAddDown | DivDown | Dot3Down => Some(false),
            MulUp | MulDivUp | MulAddUp | DivUp | Dot3Up => Some(true),
            _ => None,
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
            EllipK => x > zero && x <= one,
            EllipKm => x >= zero && x < one,
            EllipE => x >= zero && x <= one,
            JacobiSn | JacobiCn | JacobiDn | JacobiAm => a[1] >= zero && a[1] <= one,
            JacobiInverse => a[3] >= zero && a[3] <= one && (a[1] != zero || a[2] != zero),
            CarlsonRf | CarlsonRfRoots => a.iter().filter(|&&v| v == zero).count() <= 1 && (fun_roots(self) || a.iter().all(|&v| v >= zero)),
            CarlsonRj => a[..3].iter().filter(|&&v| v == zero).count() <= 1 && a.iter().all(|&v| v >= zero) && a[3] > zero,
            CarlsonRc => x >= zero && a[1] > zero,
            CarlsonRd => x >= zero && a[1] >= zero && a[2] > zero && (x != zero || a[1] != zero),
            EllipPi => {
                let sn = a[1].to_float();
                let pole = Float::with_val(PREC, 1) - Float::with_val(PREC, a[0].to_float() * sn.square());
                let reachable = a[2] >= zero || a[0] < one;
                a[4] >= zero && a[4] <= one && pole.is_sign_positive() && !pole.is_zero() && reachable && (a[2] != zero || a[3] != zero)
            }
            EllipPiComplete => x < one && a[1] > zero && a[1] <= one,
            EllipPiIncrement => x < one && a[3] >= zero && a[3] <= one,
            E1Re | E1Im => x != zero || a[1] != zero,
            GammaUpperRe | GammaUpperIm => x > zero && x.to_float() <= 128,
            CLnRe | CLnIm => x != zero || a[1] != zero,
            CPowRe | CPowIm => x != zero || a[1] != zero || a[2] >= zero,
            WeightedK1Re | WeightedK1Im | WeightedK2Re | WeightedK2Im => a[2] >= zero && a[2] < one && a[3] >= zero,
            KeplerUniversal => {
                let (r, s, al) = (a[0].to_float(), a[1].to_float(), a[2].to_float());
                let h2 = Float::with_val(PREC * 3, &r * (Float::with_val(PREC * 3, 2) - Float::with_val(PREC * 3, &al * &r))) - Float::with_val(PREC * 3, s.square_ref());
                x > zero && h2 > 0
            }
            Gamma => x > zero,
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
            MulAdd | MulAddDown | MulAddUp => Float::with_val(PREC * 3, &x * &a[1]) + &a[2],
            MulDiv | MulDivDown | MulDivUp => ratio(Float::with_val(PREC * 3, &x * &a[1]), &a[2]),
            DivOrSaturate | DivDown | DivUp => ratio(x, &a[1]),
            MulDown | MulUp => Float::with_val(PREC * 3, &x * &a[1]),
            Dot3 | Dot3Down | Dot3Up => {
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
            EllipK | EllipKm | EllipE | JacobiSn | JacobiCn | JacobiDn | JacobiAm | JacobiInverse
            | CarlsonRf | CarlsonRfRoots | CarlsonRj | CarlsonRc | CarlsonRd | EllipPi
            | EllipPiComplete | EllipPiIncrement => elliptic(self, a),
            SinTurns | CosTurns | DivByTau | WrapRadians | WrapPrincipal => {
                let p = PREC * 3;
                let tau = Float::with_val(p, Constant::Pi) * 2u32;
                let t = Float::with_val(p, &x);
                match self {
                    SinTurns => (t * tau).sin(),
                    CosTurns => (t * tau).cos(),
                    DivByTau => t / tau,
                    _ => {
                        let r = Float::with_val(p, &t - (t.clone() / &tau).floor() * &tau);
                        let half = Float::with_val(p, &tau / 2u32);
                        if self == WrapPrincipal && r > half { r - tau } else { r }
                    }
                }
            }
            CMulRe | CMulIm | CDivRe | CDivIm => {
                let p = PREC * 3;
                let u = rug::Complex::with_val(p, (&a[0], &a[1]));
                let v = rug::Complex::with_val(p, (&a[2], &a[3]));
                let r = if matches!(self, CMulRe | CMulIm) {
                    u * v
                } else if v.is_zero() {
                    let part = |c: &Float| ratio(c.clone(), &Float::with_val(PREC, 0));
                    return if self == CDivRe { part(u.real()) } else { part(u.imag()) };
                } else {
                    u / v
                };
                if matches!(self, CMulRe | CDivRe) { r.real().clone() } else { r.imag().clone() }
            }
            CExpRe | CExpIm | CExprelRe | CExprelIm | CExprel2Re | CExprel2Im => {
                let p = PREC * 3;
                let z = rug::Complex::with_val(p, (&a[0], &a[1]));
                let r = if z.is_zero() {
                    let at = match self { CExpRe | CExpIm => 1.0, CExprelRe | CExprelIm => 1.0, _ => 0.5 };
                    rug::Complex::with_val(p, (at, 0))
                } else {
                    let e = crate::complex_ref::exp_capped(&z);
                    let em1 = rug::Complex::with_val(p, &e - 1u32);
                    match self {
                        CExpRe | CExpIm => e,
                        CExprelRe | CExprelIm => em1 / &z,
                        _ => (em1 - &z) / rug::Complex::with_val(p, &z * &z),
                    }
                };
                if (self as u32) % 2 == (CExpRe as u32) % 2 { r.real().clone() } else { r.imag().clone() }
            }
            FaddeevaRe | FaddeevaIm | CErfRe | CErfIm | E1Re | E1Im => {
                use crate::complex_ref as cr;
                let r = match self {
                    FaddeevaRe | FaddeevaIm => cr::faddeeva(&a[0], &a[1]),
                    CErfRe | CErfIm => cr::erf(&a[0], &a[1]),
                    _ => cr::e1(&a[0], &a[1]),
                };
                if (self as u32) % 2 == (CExpRe as u32) % 2 { r.real().clone() } else { r.imag().clone() }
            }
            GammaUpperRe | GammaUpperIm => {
                let r = crate::complex_ref::gamma_upper(&a[0], &a[1], &a[2]);
                if self == GammaUpperRe { r.real().clone() } else { r.imag().clone() }
            }
            CLnRe | CLnIm => {
                let r = rug::Complex::with_val(PREC * 3, (&a[0], &a[1])).ln();
                if self == CLnRe { r.real().clone() } else { r.imag().clone() }
            }
            CPowRe | CPowIm => {
                let r = crate::complex_ref::pow(&a[0], &a[1], &a[2]);
                if self == CPowRe { r.real().clone() } else { r.imag().clone() }
            }
            WeightedK1Re | WeightedK1Im | WeightedK2Re | WeightedK2Im => {
                let (k1, k2) = weighted_reference(a);
                let r = match self {
                    WeightedK1Re | WeightedK1Im => k1,
                    _ => k2,
                };
                if matches!(self, WeightedK1Re | WeightedK2Re) { r.real().clone() } else { r.imag().clone() }
            }
            KeplerUniversal => kepler(a),
            MulByTau => {
                let p = PREC * 3;
                Float::with_val(p, &x) * (Float::with_val(p, Constant::Pi) * 2u32)
            }
            Erf => Float::with_val(PREC * 3, &x).erf(),
            Erfc => Float::with_val(PREC * 3, &x).erfc(),
            FresnelC | FresnelS => {
                let p = PREC * 3;
                let c = Float::with_val(p, Constant::Pi).sqrt() / 2u32 * &x;
                let e = crate::complex_ref::erf(&c, &Float::with_val(p, -&c));
                let half = rug::Complex::with_val(p, (0.5, 0.5));
                let r = half * e;
                if self == FresnelC { r.real().clone() } else { r.imag().clone() }
            }
            Gamma => Float::with_val(PREC * 3, &x).gamma(),
            Atan2Turns => {
                let p = PREC * 3;
                Float::with_val(p, x.atan2(&a[1])) / (Float::with_val(p, Constant::Pi) * 2u32)
            }
            _ => special(self, &Float::with_val(PREC * 3, &x)),
        }
    }

    /// The period an error is taken modulo, for the wraps.
    pub fn period(self) -> Option<Float> {
        matches!(self, WrapRadians | WrapPrincipal).then(|| Float::with_val(PREC, Constant::Pi) * 2u32)
    }

    /// Too costly to sweep exhaustively; sampled on every type instead.
    pub fn heavy(self) -> bool {
        matches!(
            self,
            FaddeevaRe | FaddeevaIm | CErfRe | CErfIm | E1Re | E1Im | GammaUpperRe | GammaUpperIm
                | FresnelC | FresnelS | CPowRe | CPowIm | WeightedK1Re | WeightedK1Im | WeightedK2Re
                | WeightedK2Im |
            EllipK | EllipKm | EllipE | JacobiSn | JacobiCn | JacobiDn | JacobiAm | JacobiInverse
                | CarlsonRf | CarlsonRfRoots | CarlsonRj | CarlsonRc | CarlsonRd | EllipPi
                | EllipPiComplete | EllipPiIncrement
        )
    }

    /// Whether the result is exact before its one rounding, so its error is at
    /// most half an ulp.
    pub fn rounds_once(self) -> bool {
        matches!(self, Sqrt | SqrtOrZero | MulAdd | MulDiv | DivOrSaturate | Dot3)
    }

    pub fn parity(self) -> Parity {
        match self {
            Sin | Tan | Asin | Atan | Sinh | Tanh | Coth | Asinh | Atanh | Acoth | SinTurns
            | DivByTau | MulByTau | Erf | FresnelC | FresnelS => Parity::Odd,
            Cos | Cosh | TanRatio | AtanRatio | AtanhRatio | TanhRatio | CosTurns => Parity::Even,
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
            SinTurns => vec![(f(-0.25), f(0.25), true)],
            DivByTau | MulByTau | Erf => vec![(ninf, inf, true)],
            Erfc => vec![(ninf, inf, false)],
            Expm1 | Exprel | Exprel2 | Log1p => vec![(ninf, inf, true)],
            MomentE1 | MomentE2 | MomentE3 | MomentE4 | MomentE5 | Log1pRatio => vec![(ninf, inf, false)],
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

fn is_weighted(fun: Fun) -> bool {
    matches!(fun, WeightedK1Re | WeightedK1Im | WeightedK2Re | WeightedK2Im)
}

/// The four parts of `exprel_weighted` share their points and references.
fn weighted_reference(a: &[Float]) -> (rug::Complex, rug::Complex) {
    use std::collections::HashMap;
    use std::sync::Mutex;
    static CACHE: std::sync::OnceLock<Mutex<HashMap<String, (rug::Complex, rug::Complex)>>> = std::sync::OnceLock::new();
    let key = a.iter().map(|v| v.to_string_radix(16, None)).collect::<Vec<_>>().join(",");
    let cache = CACHE.get_or_init(|| Mutex::new(HashMap::new()));
    if let Some(hit) = cache.lock().unwrap().get(&key) {
        return hit.clone();
    }
    let value = crate::complex_ref::exprel_weighted(&a[0], &a[1], &a[2], &a[3]);
    cache.lock().unwrap().insert(key, value.clone());
    value
}

fn seed(fun: Fun, ty: &str) -> u64 {
    let mut h: u64 = 0xCBF2_9CE4_8422_2325;
    let name = if is_weighted(fun) { "exprel_weighted" } else { fun.name() };
    for b in name.bytes().chain(ty.bytes()) {
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
    if fun.heavy() {
        return elliptic_points(fun, rng);
    }
    if fun.arity() == 1 {
        let xs = if T::bits() <= exhaustive_bits() {
            all::<T>()
        } else {
            candidates::<T>(fun, rng, scaled(SAMPLES))
        };
        return xs.into_iter().map(|x| vec![x]).collect();
    }
    if fun.arity() > 2 {
        let xs = if T::bits() <= 16 {
            all::<T>()
        } else {
            candidates::<T>(fun, rng, scaled(SAMPLES))
        };
        let n = scaled(if T::bits() <= 16 { PAIRS_16 } else { SAMPLES });
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
    let n = scaled(if T::bits() <= 16 { PAIRS_16 } else { SAMPLES });
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
        let mut diff = Float::with_val(PREC, c.to_float() - &expected);
        if let Some(period) = fun.period() {
            let half = Float::with_val(PREC, &period / 2u32);
            if diff > half {
                diff -= &period;
            } else if diff < -half.clone() {
                diff += &period;
            }
        }
        if let Some(up) = fun.direction() {
            // Exactly: the wrong way, or an ulp or more off, is a failure.
            let exact = Float::with_val(PREC * 4, c.to_float() - &expected);
            let wrong = !exact.is_zero() && exact.is_sign_negative() == up;
            if wrong || (Float::with_val(PREC * 4, exact.abs()) << T::frac_n()) >= 1 {
                s.errors += 1;
                if s.example.is_empty() {
                    s.example = format!("not rounded {} at {}", if up { "up" } else { "down" }, show(p));
                }
            }
        }
        let err = diff.abs() << T::frac_n();
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
