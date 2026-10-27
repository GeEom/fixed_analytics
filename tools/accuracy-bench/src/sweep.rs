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
}

use Fun::*;

pub const ALL: &[Fun] = &[
    Sin, Cos, Tan, Asin, Acos, Atan, Atan2, Sinh, Cosh, Tanh, Coth, Asinh, Acosh, Atanh, Acoth,
    Exp, Pow2, Ln, Log2, Log10, Pow, Sqrt,
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
        }
    }

    pub fn arity(self) -> usize {
        match self {
            Atan2 | Pow => 2,
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
            Ln | Log2 | Log10 => x > zero,
            Sqrt => x >= zero,
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
        }
    }

    pub fn parity(self) -> Parity {
        match self {
            Sin | Tan | Asin | Atan | Sinh | Tanh | Coth | Asinh | Atanh | Acoth => Parity::Odd,
            Cos | Cosh => Parity::Even,
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
            Sqrt => vec![(zero - f(1e-300), inf, true)],
            Atan2 | Pow => vec![],
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
    let pi_t = T::from_float(&pi()).raw();
    centers.extend([pi_t, -pi_t, pi_t / 2, -pi_t / 2]);
    centers.extend(fun.hot::<T>().iter().map(|h| T::from_float(h).raw()));
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
        .map(|r| T::from_raw(clamp_raw::<T>(r)))
        .collect()
}

fn exponent_candidates<T: Sut>(rng: &mut Rng, n: usize) -> Vec<T> {
    let mut v: Vec<T> = (-16..=16).map(|k: i32| <T as CordicNumber>::from_num(k)).collect();
    v.extend(candidates::<T>(Pow, rng, n));
    v
}

fn all<T: Sut>() -> Vec<T> {
    (T::min_raw()..=T::max_raw()).map(T::from_raw).collect()
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
    let ys = if fun == Pow {
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
        if z != want {
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
            let r = p[0].raw();
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
            let x = T::from_raw(r);
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
