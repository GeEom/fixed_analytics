//! Reference complex special functions in MPC, at precisions that cover each
//! series' cancellation, by methods unlike the crate's where they can be.

use rug::float::Constant;
use rug::{Complex, Float};


fn cx(re: &Float, im: &Float, p: u32) -> Complex {
    Complex::with_val(p, (re, im))
}

/// `e^z` with its real part capped at 10⁴, so that MPFR, whose exponents stop
/// near 2^(2^30), keeps the phase of values every type saturates from.
pub fn exp_capped(z: &Complex) -> Complex {
    let p = z.prec().0;
    if *z.real() > 10_000 {
        return Complex::with_val(p, (10_000, z.imag())).exp();
    }
    z.clone().exp()
}

/// `z^s·e^(−z)`, capped like [`exp_capped`]; on the negative axis its phase
/// `πs` is exact, where a rounded one times the huge modulus would swamp a
/// zero part.
fn power_exp(s: &Float, z: &Complex, p: u32) -> Complex {
    if z.imag().is_zero() && z.real().is_sign_negative() {
        let x = Float::with_val(p, z.real().abs_ref());
        let size = (Float::with_val(p, x.ln_ref()) * s + &x).min(&Float::with_val(p, 10_000)).exp();
        return Complex::with_val(p, (s.cos_pi_ref(), s.sin_pi_ref())) * size;
    }
    exp_capped(&(Complex::with_val(p, z.ln_ref()) * s - z))
}

/// `Γ(a, z) ~ z^(a−1)·e^(−z)·Σ (a − 1)…(a − k)/z^k`, for `|z|` past 300 near
/// the negative axis, where the continued fraction converges too slowly.
fn gamma_asymptotic(a: &Float, z: &Complex, p: u32) -> Complex {
    let inv = Complex::with_val(p, z.recip_ref());
    let mut term = Complex::with_val(p, (1, 0));
    let mut sum = term.clone();
    for k in 1..5000u32 {
        let next = Complex::with_val(p, &term * &inv) * Float::with_val(p, a - k);
        if next.is_zero() || Float::with_val(64, next.abs_ref()) > Float::with_val(64, term.abs_ref()) {
            break;
        }
        sum += &next;
        term = next;
        if settled(&term, &sum, 300) {
            break;
        }
    }
    power_exp(&Float::with_val(p, a - 1u32), z, p) * sum
}

fn settled(term: &Complex, sum: &Complex, bits: i32) -> bool {
    let t = Float::with_val(64, term.abs_ref());
    let s = Float::with_val(64, sum.abs_ref());
    t.is_zero() || (s.get_exp().unwrap_or(0) - t.get_exp().unwrap_or(0) > bits)
}

/// `erf u` by its Maclaurin series `(2/√π)·Σ (−1)ⁿu^(2n+1)/(n!(2n + 1))`.
fn erf_series(u: &Complex, p: u32) -> Complex {
    let u2 = Complex::with_val(p, u * u);
    let mut power = u.clone();
    let mut sum = u.clone();
    let mut n: u32 = 0;
    loop {
        n += 1;
        power = -(power * &u2) / n;
        let term = Complex::with_val(p, &power / (2 * n + 1));
        sum += &term;
        if n as f64 > Float::with_val(64, u2.abs_ref()).to_f64() && settled(&term, &sum, 300) {
            break;
        }
    }
    sum * (Float::with_val(p, 2) / Float::with_val(p, Constant::Pi).sqrt())
}

/// `w(z)` from `(i/√π)/z·Σ (2n − 1)!!/(2z²)ⁿ`, for `|z| ≥ 20`.
fn w_asymptotic(z: &Complex, p: u32) -> Complex {
    let inv = Complex::with_val(p, Complex::with_val(p, z * z) * 2u32).recip();
    let mut term = Complex::with_val(p, (1, 0));
    let mut sum = term.clone();
    for n in 1..2000u32 {
        let next = Complex::with_val(p, &term * &inv) * (2 * n - 1);
        if Float::with_val(64, next.abs_ref()) > Float::with_val(64, term.abs_ref()) {
            break;
        }
        sum += &next;
        term = next;
        if settled(&term, &sum, 300) {
            break;
        }
    }
    let root_pi = Float::with_val(p, Constant::Pi).sqrt();
    (sum / z) * Complex::with_val(p, (0, 1)) / root_pi
}

/// Past this `|z|` the asymptotic series for `w`, whose smallest term is
/// about `e^(−|z|²)`, is good to 2^-174.
const FAR: u32 = 11;

/// `|z|`, roughly.
fn modulus(re: &Float, im: &Float) -> f64 {
    (Float::with_val(64, re.square_ref()) + Float::with_val(64, im.square_ref())).sqrt().to_f64()
}

/// Bits that cover `e^(|z|²)` of cancellation with 320 to spare.
fn bits_for(re: &Float, im: &Float) -> u32 {
    let r = modulus(re, im).min(f64::from(FAR));
    (320.0 + 1.5 * r * r) as u32
}

pub fn faddeeva(re: &Float, im: &Float) -> Complex {
    let p = bits_for(re, im);
    let z = cx(re, im, p);
    if Float::with_val(64, z.abs_ref()) >= FAR {
        if im.is_sign_negative() && !im.is_zero() {
            let e = exp_capped(&Complex::with_val(p, -Complex::with_val(p, &z * &z))) * 2u32;
            return e - w_asymptotic(&(-z), p);
        }
        return w_asymptotic(&z, p);
    }
    // w = e^(−z²)·erfc(−iz) = e^(−z²)·(1 − erf(−iz)).
    let u = Complex::with_val(p, &z * Complex::with_val(p, (0, -1)));
    let one = Complex::with_val(p, (1, 0));
    Complex::with_val(p, -Complex::with_val(p, &z * &z)).exp() * (one - erf_series(&u, p))
}

/// `erf z`, imaginary on the imaginary axis, where the asymptotic `w` leaves
/// out the real part `e^(−x²)` that `1 − e^(−z²)·w` would cancel.
pub fn erf(re: &Float, im: &Float) -> Complex {
    let p = bits_for(re, im);
    let z = cx(re, im, p);
    if Float::with_val(64, z.abs_ref()) < FAR {
        return erf_series(&z, p);
    }
    let flip = re.is_sign_negative();
    let z = if flip { -z } else { z };
    let iz = Complex::with_val(p, &z * Complex::with_val(p, (0, 1)));
    let w = faddeeva(iz.real(), iz.imag());
    let one = Complex::with_val(p, (1, 0));
    let v = one - exp_capped(&Complex::with_val(p, -Complex::with_val(p, &z * &z))) * w;
    let v = if flip { -v } else { v };
    if re.is_zero() { Complex::with_val(p, (0, v.imag())) } else { v }
}

/// A continued fraction `b₀ + a₁/(b₁ + …)` by Lentz's method at `p` bits.
fn lentz(b0: Complex, p: u32, term: impl Fn(u32) -> (Complex, Complex)) -> Complex {
    let tiny = Complex::with_val(p, (Float::with_val(p, Float::i_exp(1, -3000)), 0));
    let mut f = if b0.is_zero() { tiny.clone() } else { b0 };
    let mut c = f.clone();
    let mut d = Complex::with_val(p, (0, 0));
    for k in 1..100_000u32 {
        let (a, b) = term(k);
        d = Complex::with_val(p, &b + Complex::with_val(p, &a * &d));
        if d.is_zero() {
            d = tiny.clone();
        }
        c = Complex::with_val(p, &b + Complex::with_val(p, &a / &c));
        if c.is_zero() {
            c = tiny.clone();
        }
        d = d.recip();
        let delta = Complex::with_val(p, &c * &d);
        f *= &delta;
        let diff = Complex::with_val(p, &delta - 1u32);
        let off = Float::with_val(64, diff.abs_ref());
        if off.is_zero() || off.get_exp().unwrap_or(0) < -300 {
            break;
        }
    }
    f
}

pub fn e1(re: &Float, im: &Float) -> Complex {
    let z = cx(re, im, (400.0 + 2.0 * modulus(re, im).min(1000.0)) as u32);
    let r = Float::with_val(64, z.abs_ref()).to_f64();
    let near_neg = re.is_sign_negative() && Float::with_val(64, im.abs_ref()) <= Float::with_val(64, re.abs_ref());
    if near_neg && r > 300.0 {
        let v = gamma_asymptotic(&Float::with_val(512, 0), &cx(re, im, 512), 512);
        // On the axis, the value from above it: −Ei(−z) − iπ.
        let pi = Float::with_val(512, Constant::Pi);
        return if im.is_zero() { Complex::with_val(512, (v.real(), -pi)) } else { v };
    }
    if r <= 60.0 || (near_neg && r <= 300.0) {
        let p = z.prec().0;
        let mut term = Complex::with_val(p, (1, 0));
        let mut sum = Complex::with_val(p, (0, 0));
        let mut n: u32 = 0;
        loop {
            n += 1;
            term = -(term * &z) / n;
            let part = Complex::with_val(p, &term / n);
            sum += &part;
            if n as f64 > r && settled(&part, &sum, 300) {
                break;
            }
        }
        let gamma = Float::with_val(p, Constant::Euler);
        return -(Complex::with_val(p, (gamma, 0)) + z.ln() + sum);
    }
    let p = 512;
    let z = cx(re, im, p);
    let one = Complex::with_val(p, (1, 0));
    let f = lentz(Complex::with_val(p, &z + &one), p, |k| {
        let i = Float::with_val(p, k);
        (Complex::with_val(p, (-(i.clone() * &i), 0)), Complex::with_val(p, &z + &one) + Complex::with_val(p, (i * 2u32, 0)))
    });
    exp_capped(&Complex::with_val(p, -z)) / f
}

pub fn gamma_upper(a: &Float, re: &Float, im: &Float) -> Complex {
    let p = (512.0 + 2.0 * modulus(re, im).min(1000.0)) as u32;
    let z = cx(re, im, p);
    if z.is_zero() {
        return Complex::with_val(p, (Float::with_val(p, a).gamma(), 0));
    }
    if im.is_zero() && re.is_sign_positive() && *re <= 100 {
        let v = Float::with_val(p, a).gamma_inc(&Float::with_val(p, re));
        return Complex::with_val(p, (v, 0));
    }
    let r = Float::with_val(64, z.abs_ref()).to_f64();
    let near_neg = re.is_sign_negative() && Float::with_val(64, im.abs_ref()) <= Float::with_val(64, re.abs_ref());
    if near_neg && r > 300.0 {
        let v = gamma_asymptotic(a, &z, p);
        // On the axis Γ(a, z) = Γ(a) − e^(iπa)·G with G real; the series gives
        // the second part.
        if im.is_zero() && !a.is_integer() {
            return Complex::with_val(p, (Float::with_val(p, v.real() + Float::with_val(p, a).gamma()), v.imag()));
        }
        return v;
    }
    let lead = power_exp(&Float::with_val(p, a), &z, p);
    if r <= (a.to_f64() + 1.0).max(8.0) || (near_neg && r <= 300.0) {
        let mut denom = Float::with_val(p, a);
        let mut term = Complex::with_val(p, (Float::with_val(p, 1) / &denom, 0));
        let mut sum = term.clone();
        loop {
            denom += 1u32;
            term = Complex::with_val(p, &term * &z) / &denom;
            sum += &term;
            if denom.to_f64() > r && settled(&term, &sum, 300) {
                break;
            }
        }
        return Complex::with_val(p, (Float::with_val(p, a).gamma(), 0)) - lead * sum;
    }
    let one_minus_a = Float::with_val(p, 1) - a;
    let f = lentz(Complex::with_val(p, &z + &one_minus_a), p, |k| {
        let i = Float::with_val(p, k);
        let ai = -(i.clone() * (i.clone() - a));
        let b = Complex::with_val(p, &z + &one_minus_a) + Float::with_val(p, i * 2u32);
        (Complex::with_val(p, (ai, 0)), b)
    });
    lead / f
}

/// `z^a` for real `a`: on the axes and diagonals, `|z|^a` turned by an exact
/// `a·k·π/4`, else `e^(a·ln z)`, both capped like [`exp_capped`].
pub fn pow(re: &Float, im: &Float, a: &Float) -> Complex {
    let p = 960;
    let z = cx(re, im, p);
    if z.is_zero() {
        return Complex::with_val(p, (u32::from(a.is_zero()), 0));
    }
    let eighths = if im.is_zero() {
        Some(if re.is_sign_negative() { 4 } else { 0 })
    } else if re.is_zero() || re.clone().abs() == im.clone().abs() {
        let k = if re.is_zero() { 2 } else if re.is_sign_positive() { 1 } else { 3 };
        Some(if im.is_sign_negative() { -k } else { k })
    } else {
        None
    };
    let Some(k) = eighths else {
        return exp_capped(&(Complex::with_val(p, z.ln_ref()) * a));
    };
    let ln_abs = Float::with_val(p, z.abs_ref()).ln();
    let size = (ln_abs * a).min(&Float::with_val(p, 10_000)).exp();
    let phase = Float::with_val(p, a * k) / 4u32;
    Complex::with_val(p, (phase.cos_pi_ref(), phase.sin_pi_ref())) * size
}

/// Bits for [`exprel_weighted`]'s quadrature.
const WP: u32 = 320;

/// Gauss–Legendre nodes and weights on `[0, 1]`, 40 of them, by Newton's method.
fn legendre() -> &'static [(Float, Float)] {
    static NODES: std::sync::OnceLock<Vec<(Float, Float)>> = std::sync::OnceLock::new();
    NODES.get_or_init(|| {
        let n = 40u32;
        let p = WP + 64;
        let pi = Float::with_val(p, Constant::Pi);
        (1..=n)
            .map(|i| {
                let mut x = Float::with_val(p, (Float::with_val(p, i) - 0.25) * &pi / (f64::from(n) + 0.5)).cos();
                let mut deriv = Float::with_val(p, 0);
                for _ in 0..100 {
                    let (mut p0, mut p1) = (Float::with_val(p, 1), x.clone());
                    for k in 2..=n {
                        let p2 = (Float::with_val(p, &x * &p1) * (2 * k - 1) - Float::with_val(p, &p0 * (k - 1))) / k;
                        p0 = p1;
                        p1 = p2;
                    }
                    deriv = Float::with_val(p, Float::with_val(p, &x * &p1) - &p0) * n / (Float::with_val(p, &x * &x) - 1u32);
                    let dx = Float::with_val(p, &p1 / &deriv);
                    x -= &dx;
                    if dx.is_zero() || dx.get_exp().unwrap_or(0) < -(p as i32) {
                        break;
                    }
                }
                let w = Float::with_val(p, 2) / ((Float::with_val(p, 1) - Float::with_val(p, &x * &x)) * deriv.square());
                ((x + 1u32) / 2u32, w / 2u32)
            })
            .collect()
    })
}

/// `(K₁, K₂)` of `Complex::exprel_weighted`: past `|z| = 8(a + 170)·δ` by
/// parts, else Gauss–Legendre panels on `e^(z(1 − u))·{W, V}(u)/(1 − ε)` from
/// each end toward the least of the integrand; NaN where that would take more
/// than 800 panels.
pub fn exprel_weighted(re: &Float, im: &Float, eps: &Float, pw: &Float) -> (Complex, Complex) {
    let p = WP;
    let z = cx(re, im, p);
    let one = Float::with_val(p, 1);
    if eps.is_zero() {
        // exprel and exprel2, at the bits their cancellation near 0 needs.
        let p = 1024;
        let z = cx(re, im, p);
        if z.is_zero() {
            return (Complex::with_val(p, (1, 0)), Complex::with_val(p, (0.5, 0)));
        }
        let em1 = exp_capped(&z) - Complex::with_val(p, (1, 0));
        let k1 = Complex::with_val(p, &em1 / &z);
        let k2 = (em1 - &z) / Complex::with_val(p, &z * &z);
        return (k1, k2);
    }
    let a = Float::with_val(p, pw + 1u32);
    let keep = Float::with_val(p, &one - eps);
    let delta = Float::with_val(p, eps / &keep);
    let log_keep = Float::with_val(p, -eps).ln_1p();
    let r = Float::with_val(64, z.abs_ref()).to_f64();
    let nan = || (Complex::with_val(p, (f64::NAN, f64::NAN)), Complex::with_val(p, (f64::NAN, f64::NAN)));
    let a64 = a.to_f64();
    let d64 = delta.to_f64();
    if r >= (8.0 * (a64 + 170.0) * d64).max(1.0) {
        let ez = exp_capped(&z);
        let sum = |shift: &Float, q: &Float, first: Complex| {
            let step = Complex::with_val(p, z.recip_ref()) * q;
            let mut term = first;
            let mut total = term.clone();
            for k in 0..2000u32 {
                term = -(term * &step) * Float::with_val(p, shift + k);
                total += &term;
                if settled(&term, &total, p as i32 + 20) {
                    break;
                }
            }
            total
        };
        let inv = Complex::with_val(p, z.recip_ref());
        let inv2 = Complex::with_val(p, &inv * &inv);
        let pow_p = (Float::with_val(p, pw * &log_keep)).exp();
        let pow_a = (Float::with_val(p, &a * &log_keep)).exp();
        let k1 = Complex::with_val(p, &ez * sum(&a, &delta, inv.clone())) / &keep - sum(&a, eps, inv.clone()) * pow_p;
        let a1 = Float::with_val(p, &a + 1u32);
        let tail = -(Float::with_val(p, &a * &log_keep)).exp_m1() / Float::with_val(p, &a * eps);
        let k2 = -(inv * tail) + Complex::with_val(p, &ez * sum(&a1, &delta, inv2.clone())) / &keep
            - sum(&a1, eps, inv2) * pow_a;
        return (k1, k2);
    }
    // Integrand parts at u: e^(z(from))·W(u) and e^(z(from))·V(u), from = −u
    // forward and 1 − u backward, so neither grows; e^z multiplies the first.
    let at = |u: &Float, forward: bool| {
        let l = Float::with_val(p, Float::with_val(p, &delta * u).ln_1p() * &a);
        let from = if forward { Float::with_val(p, -u) } else { Float::with_val(p, &one - u) };
        let lead = Complex::with_val(p, (Float::with_val(p, &from * re), Float::with_val(p, &from * im)));
        let base = lead.exp();
        let w = Float::with_val(p, -&l).exp();
        let v = -Float::with_val(p, -&l).exp_m1() / Float::with_val(p, &a * &delta);
        (Complex::with_val(p, &base * &w), base * v, Float::with_val(p, &l + Float::with_val(p, re * u)))
    };
    let meet = if re.is_sign_negative() && !re.is_zero() {
        Float::with_val(p, -(Float::with_val(p, &a / re)) - Float::with_val(p, &one / &delta)).clamp(&0, &1)
    } else {
        one.clone()
    };
    let nodes = legendre();
    let (mut k1, mut k2) = (Complex::with_val(p, (0, 0)), Complex::with_val(p, (0, 0)));
    let (mut k1f, mut k2f) = (k1.clone(), k2.clone());
    let mut panels = 0;
    for forward in [true, false] {
        let mut u = if forward { Float::with_val(p, 0) } else { one.clone() };
        loop {
            let rest = if forward { Float::with_val(p, &meet - &u) } else { Float::with_val(p, &u - &meet) };
            if rest <= 0 {
                break;
            }
            // The integrand's size only falls toward the meeting point, so its
            // size here times what is left bounds the rest.
            let (f1, f2, _) = at(&u, forward);
            let sum = if forward { (&k1f, &k2f) } else { (&k1, &k2) };
            let small = |f: &Complex, s: &Complex| {
                let left = Float::with_val(64, f.abs_ref()) * rest.to_f64();
                let got = Float::with_val(64, s.abs_ref());
                !got.is_zero() && left < got * Float::with_val(64, Float::i_exp(1, -(p as i32) - 20))
            };
            if small(&f1, sum.0) && small(&f2, sum.1) {
                break;
            }
            let rr = d64 / (1.0 + d64 * u.to_f64());
            let slope = (-re.to_f64() - a64 * rr).abs();
            let h = (5.0 / (im.to_f64().abs() + slope + 1.0)).min(0.25 / rr).min(1.5 / (rr * a64.sqrt()));
            let h = Float::with_val(p, h).min(&rest);
            panels += 1;
            if panels > 800 {
                return nan();
            }
            let lo = if forward { u.clone() } else { Float::with_val(p, &u - &h) };
            for (x, w) in nodes {
                let (f1, f2, _) = at(&Float::with_val(p, &lo + Float::with_val(p, x * &h)), forward);
                let wh = Float::with_val(p, w * &h);
                let (s1, s2) = if forward { (&mut k1f, &mut k2f) } else { (&mut k1, &mut k2) };
                *s1 += f1 * &wh;
                *s2 += f2 * &wh;
            }
            if forward { u += &h } else { u -= &h }
        }
    }
    let ez = exp_capped(&z);
    let k1 = (Complex::with_val(p, &ez * &k1f) + k1) / &keep;
    // V's quadrature stops where W's does, so past a small aδ, K₂ from the
    // identity K₂ = (φ₁(z) − (1 − ε)K₁)/(aε) at this precision.
    let ad = Float::with_val(p, &a * &delta);
    let k2 = if ad > Float::with_val(p, Float::i_exp(1, -8)) {
        let phi1 = if z.is_zero() { Complex::with_val(p, (1, 0)) } else { (Complex::with_val(p, &ez - 1u32)) / &z };
        (phi1 - Complex::with_val(p, &k1 * &keep)) / Float::with_val(p, &a * eps)
    } else {
        (ez * k2f + k2) / &keep
    };
    (k1, k2)
}
