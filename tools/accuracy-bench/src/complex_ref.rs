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
