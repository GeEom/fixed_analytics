//! Reference elliptic functions at 512 bits, independent of the crate's: the
//! AGM from MPFR, Carlson's duplication run until its deviations are below
//! 2^-70, and the Jacobi functions by the AGM with plain arcsines.

use rug::Float;
use rug::float::Constant;

pub const P: u32 = 512;

fn f(x: f64) -> Float {
    Float::with_val(P, x)
}

fn big(x: &Float) -> Float {
    Float::with_val(P, x)
}

/// Duplication until the largest deviation `|1 − v/A|` is below 2^-70.
fn duplicate(mut v: Vec<Float>, weights: &[u32], mut each: impl FnMut(&[Float], &[Float], &Float)) -> (Vec<Float>, Float) {
    let total: u32 = weights.iter().sum();
    let mut a = v.iter().zip(weights).fold(f(0.0), |acc, (x, &w)| acc + Float::with_val(P, x * w)) / total;
    let tiny = Float::with_val(P, Float::i_exp(1, -70));
    for _ in 0..400 {
        let dev = v.iter().map(|x| (Float::with_val(P, &a - x) / &a).abs()).fold(f(0.0), |m, d| if d > m { d } else { m });
        if dev < tiny {
            break;
        }
        let roots: Vec<Float> = v.iter().map(|x| x.clone().sqrt()).collect();
        let lam = Float::with_val(P, &roots[0] * &roots[1]) + Float::with_val(P, &roots[0] * &roots[2]) + Float::with_val(P, &roots[1] * &roots[2]);
        each(&v, &roots, &lam);
        for x in v.iter_mut() {
            *x = (x.clone() + &lam) / 4u32;
        }
        a = (a + &lam) / 4u32;
    }
    (v, a)
}

pub fn rf(x: &Float, y: &Float, z: &Float) -> Float {
    let (v, a) = duplicate(vec![big(x), big(y), big(z)], &[1, 1, 1], |_, _, _| {});
    let dx = Float::with_val(P, &a - &v[0]) / &a;
    let dy = Float::with_val(P, &a - &v[1]) / &a;
    let dz = -(dx.clone() + &dy);
    let e2 = Float::with_val(P, &dx * &dy) - dz.clone().square();
    let e3 = dx * dy * dz;
    let s = f(1.0) - Float::with_val(P, &e2 / 10u32) + Float::with_val(P, &e3 / 14u32)
        + e2.clone().square() / 24u32
        - Float::with_val(P, &e2 * &e3) * 3u32 / 44u32;
    s / a.sqrt()
}

pub fn rc(x: &Float, y: &Float) -> Float {
    if x.is_zero() {
        return Float::with_val(P, Constant::Pi) / 2u32 / y.clone().sqrt();
    }
    if x == y {
        return f(1.0) / x.clone().sqrt();
    }
    let d = Float::with_val(P, y - x);
    if d.is_sign_positive() {
        let t = (d.clone() / x).sqrt();
        t.atan() / d.sqrt()
    } else {
        let d = -d;
        let t = (d.clone() / x).sqrt();
        t.atanh() / d.sqrt()
    }
}

pub fn rj(x: &Float, y: &Float, z: &Float, p: &Float) -> Float {
    let mut sum = f(0.0);
    let mut scale = f(1.0);
    let (v, a) = duplicate(vec![big(x), big(y), big(z), big(p)], &[1, 1, 1, 2], |v, r, lam| {
        let alpha = (Float::with_val(P, &v[3] * (r[0].clone() + &r[1] + &r[2])) + Float::with_val(P, &r[0] * &r[1]) * &r[2]).square();
        let beta = Float::with_val(P, &v[3] * (v[3].clone() + lam).square());
        sum += Float::with_val(P, &scale * rc(&alpha, &beta));
        scale /= 4u32;
    });
    let dx = Float::with_val(P, &a - &v[0]) / &a;
    let dy = Float::with_val(P, &a - &v[1]) / &a;
    let dz = Float::with_val(P, &a - &v[2]) / &a;
    let dp = -(dx.clone() + &dy + &dz) / 2u32;
    let xyz = Float::with_val(P, &dx * &dy) * &dz;
    let e2 = Float::with_val(P, &dx * &dy) + Float::with_val(P, &dx * &dz) + Float::with_val(P, &dy * &dz) - dp.clone().square() * 3u32;
    let p3 = Float::with_val(P, dp.clone().square() * &dp);
    let e3 = xyz.clone() + Float::with_val(P, &e2 * &dp) * 2u32 + p3.clone() * 4u32;
    let e4 = (xyz.clone() * 2u32 + Float::with_val(P, &e2 * &dp) + p3 * 3u32) * &dp;
    let e5 = xyz * dp.square();
    let series = f(1.0) - Float::with_val(P, &e2 * 3u32) / 14u32 + Float::with_val(P, &e3 / 6u32)
        + e2.clone().square() * 9u32 / 88u32
        - e4 * 3u32 / 22u32
        - Float::with_val(P, &e2 * &e3) * 9u32 / 52u32
        + e5 * 3u32 / 26u32;
    let a32 = Float::with_val(P, &a * a.clone().sqrt());
    scale * series / a32 + sum * 3u32
}

pub fn rd(x: &Float, y: &Float, z: &Float) -> Float {
    rj(x, y, z, z)
}

pub fn k(kc: &Float) -> Float {
    Float::with_val(P, Constant::Pi) / 2u32 / f(1.0).agm(kc)
}

pub fn e(kc: &Float) -> Float {
    if kc.is_zero() {
        return f(1.0);
    }
    let kc2 = kc.clone().square();
    let m = f(1.0) - &kc2;
    rf(&f(0.0), &kc2, &f(1.0)) - m * rd(&f(0.0), &kc2, &f(1.0)) / 3u32
}

/// `(sn, cn, dn, am)` of `u`.
pub fn jacobi(u: &Float, kc: &Float) -> (Float, Float, Float, Float) {
    let one = f(1.0);
    if kc.is_zero() {
        let s = u.clone().tanh();
        let c = one.clone() / u.clone().cosh();
        let am = u.clone().sinh().atan();
        return (s, c.clone(), c, am);
    }
    let m = one.clone() - kc.clone().square();
    let (mut a, mut b, mut c) = (vec![one.clone()], vec![big(kc)], vec![m.clone().sqrt()]);
    let tiny = Float::with_val(P, Float::i_exp(1, -500));
    while c.last().map_or(false, |x| x.clone().abs() > tiny) && a.len() < 64 {
        let (an, bn) = (a.last().unwrap().clone(), b.last().unwrap().clone());
        a.push((an.clone() + &bn) / 2u32);
        b.push((an.clone() * &bn).sqrt());
        c.push((an - bn) / 2u32);
    }
    let n = a.len() - 1;
    let mut phi = Float::with_val(P, &a[n] * u) * Float::with_val(P, Float::i_exp(1, n as i32));
    for i in (1..=n).rev() {
        let x = Float::with_val(P, &c[i] / &a[i]) * phi.clone().sin();
        phi = (phi + x.asin()) / 2u32;
    }
    let (s, cn) = (phi.clone().sin(), phi.clone().cos());
    let dn = (one - m * s.clone().square()).sqrt();
    (s, cn, dn, phi)
}

/// `Π(n; φ|m)` from `sn`, `cn`, `dn` of `φ`, over `(−π, π]`.
pub fn pi(n: &Float, sn: &Float, cn: &Float, dn: &Float, kc: &Float) -> Float {
    let one = f(1.0);
    let s2 = sn.clone().square();
    let pole = one.clone() - Float::with_val(P, n * &s2);
    let (c2, d2) = (cn.clone().square(), dn.clone().square());
    let part = sn.clone().abs() * (rf(&c2, &d2, &one) + Float::with_val(P, n * &s2) * rj(&c2, &d2, &one, &pole) / 3u32);
    let total = if cn.is_sign_negative() && !cn.is_zero() { pi_complete(n, kc) * 2u32 - part } else { part };
    if sn.is_sign_negative() { -total } else { total }
}

pub fn pi_complete(n: &Float, kc: &Float) -> Float {
    let one = f(1.0);
    let kc2 = kc.clone().square();
    let pole = one.clone() - n;
    rf(&f(0.0), &kc2, &one) + Float::with_val(P, n * rj(&f(0.0), &kc2, &one, &pole)) / 3u32
}

/// `u` with these `(sn, cn, dn)`, in `(−2K, 2K]`.
pub fn inverse(sn: &Float, cn: &Float, dn: &Float, kc: &Float) -> Float {
    let part = sn.clone().abs() * rf(&cn.clone().square(), &dn.clone().square(), &f(1.0));
    let u = if cn.is_sign_negative() && !cn.is_zero() { k(kc) * 2u32 - part } else { part };
    if sn.is_sign_negative() { -u } else { u }
}

/// `u` less the nearest multiple of `2K`, or `u` itself on the separatrix.
fn reduced(u: &Float, kc: &Float) -> Float {
    if kc.is_zero() {
        return big(u);
    }
    let period = k(kc) * 2u32;
    let turns = Float::with_val(P, u / &period).round();
    Float::with_val(P, u - turns * period)
}

/// `Π(n; am v|m)` for `v` within a period of 0, continued on the separatrix past 10⁴
/// by `1/(1 − n)`, below which `sech²` has vanished.
fn pi_of_u(n: &Float, v: &Float, kc: &Float) -> Float {
    let cap = f(1e4);
    if kc.is_zero() && Float::with_val(P, v.clone().abs()) > cap {
        let edge = if v.is_sign_negative() { -cap } else { cap };
        let rest = Float::with_val(P, v - &edge) / (f(1.0) - n);
        return pi_of_u(n, &edge, kc) + rest;
    }
    if kc.is_zero() {
        let (sn, cn, dn, _) = jacobi(v, kc);
        return pi(n, &sn, &cn, &dn, kc);
    }
    // Π(v) = 2Π(n|m)·round(v/2K) + Π of the rest, within (−K, K].
    let period = k(kc) * 2u32;
    let turns = Float::with_val(P, v / &period).round();
    let rest = Float::with_val(P, v - Float::with_val(P, &turns * &period));
    let (sr, cr, dr, _) = jacobi(&rest, kc);
    let local = pi(n, &sr, &cr, &dr, kc);
    if turns.is_zero() { local } else { local + turns * pi_complete(n, kc) * 2u32 }
}

/// `∫ dv/(1 − n·sn²v)` from `u` to `u + s`.
pub fn pi_increment(n: &Float, u: &Float, s: &Float, kc: &Float) -> Float {
    let from = reduced(u, kc);
    let to = Float::with_val(P, &from + s);
    pi_of_u(n, &to, kc) - pi_of_u(n, &from, kc)
}
