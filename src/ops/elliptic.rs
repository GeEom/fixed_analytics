//! Elliptic integrals and Jacobi elliptic functions of an exactly held [`Modulus`].

use crate::error::{Error, Result};
use crate::eval::{W, atan_ratio, by_width, fixed, half_pi, log1p, pi, sin_cos, sinh_cosh, work};
use crate::real::Work;
use crate::traits::CordicNumber;

/// An elliptic modulus, from the parameter `m` or the complement `k'`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Modulus<T> {
    value: T,
    complement: bool,
}

impl<T: CordicNumber> Modulus<T> {
    /// From `m = k²`, in `[0, 1]`.
    #[must_use]
    pub fn from_parameter(m: T) -> Option<Self> {
        (m >= T::zero() && m <= T::one()).then_some(Self {
            value: m,
            complement: false,
        })
    }

    /// From the complement `k' = √(1 − m)`, in `[0, 1]`, exact where `m` is near 1.
    #[must_use]
    pub fn from_complement(kc: T) -> Option<Self> {
        (kc >= T::zero() && kc <= T::one()).then_some(Self {
            value: kc,
            complement: true,
        })
    }

    /// The parameter `m`, rounded.
    #[must_use]
    pub fn parameter(self) -> T {
        fixed(parts(self).m)
    }

    /// The complement `k'`, rounded.
    #[must_use]
    pub fn complement(self) -> T {
        fixed(parts(self).kc)
    }
}

/// `m`, `k'` and `k'²` in the working arithmetic, each from the given value.
struct Parts<V> {
    m: V,
    kc: V,
    kc2: V,
}

fn parts<T: CordicNumber>(modulus: Modulus<T>) -> Parts<W<T>> {
    let one = W::<T>::one();
    let v = work(modulus.value);
    if modulus.complement {
        Parts {
            m: (one - v) * (one + v),
            kc: v,
            kc2: v.square(),
        }
    } else {
        let kc2 = one - v;
        Parts {
            m: v,
            kc: kc2.sqrt(),
            kc2,
        }
    }
}

/// `R_C(1, 1 + e)`, `atan √e/√e` or `atanh √−e/√−e`, from `e` and a maybe tiny `1 + e`.
fn rc_one<V: Work>(e: V, one_plus_e: V) -> V {
    let t = e.abs().sqrt();
    if t.is_zero() {
        return V::one();
    }
    if e.is_neg() {
        // 1 − t = (1 + e)/(1 + t), so atanh t = log1p(2t(1 + t)/(1 + e))/2.
        let one = V::one();
        log1p((t * (one + t) / one_plus_e).ldexp(1)).ldexp(-1) / t
    } else {
        atan_ratio(t, V::one()) / t
    }
}

/// `R_C(x, y)` for `x ≥ 0` and `y > 0`, by ratios so `x` near `y` does not cancel.
fn rc<V: Work>(x: V, y: V) -> V {
    if x.is_zero() {
        return half_pi::<V>() / y.sqrt();
    }
    rc_one((y - x) / x, y / x) / x.sqrt()
}

/// The deviations `(A₀ − v₀)/(4^M A_M)` after `M` duplication steps.
fn deviation<V: Work>(a0: V, v0: V, scale: V, a: V) -> V {
    (a0 - v0) * scale / a
}

/// `R_F(x, y, z)` by duplication, any `roots` standing in for the first square roots.
fn rf<V: Work>(x0: V, y0: V, z0: V, roots: Option<[V; 3]>) -> V {
    let (mut vx, mut vy, mut vz) = (x0, y0, z0);
    let a0 = (vx + vy + vz).div_int(3);
    let mut mean = a0;
    let mut scale = V::one();
    // Octaves apart, about log₂log₂ of the ratio steps to meet; each then quarters the rest.
    for step in 0..by_width::<V>([13, 18, 23]) {
        let [sx, sy, sz] = match roots {
            Some(r) if step == 0 => r.map(Work::abs),
            _ => [vx.sqrt(), vy.sqrt(), vz.sqrt()],
        };
        let lam = sx * sy + sx * sz + sy * sz;
        mean = (mean + lam).ldexp(-2);
        vx = (vx + lam).ldexp(-2);
        vy = (vy + lam).ldexp(-2);
        vz = (vz + lam).ldexp(-2);
        scale = scale.ldexp(-2);
    }
    let dx = deviation(a0, x0, scale, mean);
    let dy = deviation(a0, y0, scale, mean);
    let dz = -(dx + dy);
    let e2 = dx * dy - dz.square();
    let e3 = dx * dy * dz;
    // DLMF 19.36.1, through the seventh order.
    let one = V::one();
    let series = one - e2.div_int(10) + e3.div_int(14) + e2.square().div_int(24)
        - (e2 * e3).mul_int(3).div_int(44)
        - (e2.square() * e2).mul_int(5).div_int(208)
        + e3.square().mul_int(3).div_int(104)
        + (e2.square() * e3).div_int(16);
    series / mean.sqrt()
}

/// The closing series of `R_J` and `R_D`, Carlson (1995).
fn rj_series<V: Work>(e2: V, e3: V, e4: V, e5: V) -> V {
    V::one() - e2.mul_int(3).div_int(14) + e3.div_int(6) + e2.square().mul_int(9).div_int(88)
        - e4.mul_int(3).div_int(22)
        - (e2 * e3).mul_int(9).div_int(52)
        + e5.mul_int(3).div_int(26)
}

/// `R_J(x, y, z, p)` for `p > 0`. A `p` far above the rest, which duplication
/// brings down only fourfold a step, DLMF 19.21.12 trades for `q` near the least,
/// `x`: `(p − x)R_J(p) + (q − x)R_J(q) = 3R_F − 3R_C(yz/x, pq/x)`.
fn rj<V: Work>(x0: V, y0: V, z0: V, p0: V, roots: Option<[V; 4]>) -> V {
    let mut sorted = [(x0, 0), (y0, 1), (z0, 2)];
    sorted.sort_unstable_by_key(|entry| entry.0);
    let [(lo, _), (mid, _), (hi, _)] = sorted;
    if p0 <= hi.ldexp(4) {
        return rj_duplicate(x0, y0, z0, p0, roots);
    }
    let gap = (mid - lo) * (hi - lo) / (p0 - lo);
    let near = lo + gap;
    let picked = roots.map(|picked| {
        let at = |i: usize| picked.get(i).copied().unwrap_or_else(V::zero);
        [at(sorted[0].1), at(sorted[1].1), at(sorted[2].1)]
    });
    let near_roots = picked.map(|[a, b, c]| [a, b, c, near.sqrt()]);
    let mut rest =
        rf(lo, mid, hi, picked).mul_int(3) - gap * rj_duplicate(lo, mid, hi, near, near_roots);
    if !lo.is_zero() {
        let xi = mid * hi / lo;
        let excess = (p0 - mid) * (p0 - hi) / ((p0 - lo) * xi);
        rest = rest - (rc_one(excess, p0 * near / (lo * xi)) / xi.sqrt()).mul_int(3);
    }
    rest / (p0 - lo)
}

/// [`rj`] by duplication, its `R_C` terms `4^(−3m)·δ/d_m²` from one `δ = Π(p − v)`.
fn rj_duplicate<V: Work>(x0: V, y0: V, z0: V, p0: V, roots: Option<[V; 4]>) -> V {
    let (mut vx, mut vy, mut vz, mut vp) = (x0, y0, z0, p0);
    let a0 = (vx + vy + vz + vp.ldexp(1)).div_int(5);
    let delta = (vp - vx) * (vp - vy) * (vp - vz);
    let mut mean = a0;
    let mut scale = V::one();
    let mut sum = V::zero();
    for step in 0..by_width::<V>([15, 21, 28]) {
        let [sx, sy, sz, sp] = match roots {
            Some(r) if step == 0 => r.map(Work::abs),
            _ => [vx.sqrt(), vy.sqrt(), vz.sqrt(), vp.sqrt()],
        };
        let lam = sx * sy + sx * sz + sy * sz;
        let den = (sp + sx) * (sp + sy) * (sp + sz);
        let ex = delta * scale.square() * scale / den.square();
        // 1 + ex = 2√vp(vp + λ)/den, as Π(√vp + √v) + Π(√vp − √v) = 2√vp(vp + λ).
        let one_plus_e = (sp * (vp + lam)).ldexp(1) / den;
        sum = sum + scale * rc_one(ex, one_plus_e) / den;
        mean = (mean + lam).ldexp(-2);
        vx = (vx + lam).ldexp(-2);
        vy = (vy + lam).ldexp(-2);
        vz = (vz + lam).ldexp(-2);
        vp = (vp + lam).ldexp(-2);
        scale = scale.ldexp(-2);
    }
    let dx = deviation(a0, x0, scale, mean);
    let dy = deviation(a0, y0, scale, mean);
    let dz = deviation(a0, z0, scale, mean);
    let dp = -(dx + dy + dz).ldexp(-1);
    let e2 = dx * dy + dx * dz + dy * dz - dp.square().mul_int(3);
    let xyz = dx * dy * dz;
    let e3 = xyz + (e2 * dp).ldexp(1) + (dp.square() * dp).ldexp(2);
    let e4 = (xyz.ldexp(1) + e2 * dp + (dp.square() * dp).mul_int(3)) * dp;
    let e5 = xyz * dp.square();
    scale * rj_series(e2, e3, e4, e5) / (mean * mean.sqrt()) + sum.mul_int(6)
}

/// `R_D(x, y, z)` by duplication.
fn rd<V: Work>(x0: V, y0: V, z0: V) -> V {
    let (mut vx, mut vy, mut vz) = (x0, y0, z0);
    let a0 = (vx + vy + vz.mul_int(3)).div_int(5);
    let mut mean = a0;
    let mut scale = V::one();
    let mut sum = V::zero();
    for _ in 0..by_width::<V>([15, 21, 28]) {
        let (sx, sy, sz) = (vx.sqrt(), vy.sqrt(), vz.sqrt());
        let lam = sx * sy + sx * sz + sy * sz;
        sum = sum + scale / (sz * (vz + lam));
        mean = (mean + lam).ldexp(-2);
        vx = (vx + lam).ldexp(-2);
        vy = (vy + lam).ldexp(-2);
        vz = (vz + lam).ldexp(-2);
        scale = scale.ldexp(-2);
    }
    let dx = deviation(a0, x0, scale, mean);
    let dy = deviation(a0, y0, scale, mean);
    let dz = -(dx + dy).div_int(3);
    let xy = dx * dy;
    let z2 = dz.square();
    let e2 = xy - z2.mul_int(6);
    let e3 = (xy.mul_int(3) - z2.mul_int(8)) * dz;
    let e4 = (xy - z2) * z2.mul_int(3);
    let e5 = xy * z2 * dz;
    scale * rj_series(e2, e3, e4, e5) / (mean * mean.sqrt()) + sum.mul_int(3)
}

/// The arithmetic-geometric mean of 1 and `b`, for `b` in `(0, 1]`.
fn agm<V: Work>(b: V) -> V {
    let (mut a, mut b) = (V::one(), b);
    for _ in 0..32 {
        if V::negligible(&(a - b), &a) {
            break;
        }
        (a, b) = ((a + b).ldexp(-1), (a * b).sqrt());
    }
    a
}

/// `K(m) = π/(2·AGM(1, k'))`, for `k' > 0`.
fn complete_k<V: Work>(kc: V) -> V {
    half_pi::<V>() / agm(kc)
}

/// `(sn, cn)` of `v` in `[0, K/2]` by the AGM, `φ_N = 2^N a_N v` brought down by
/// `φ_{n−1} = (φ_n + asin((c_n/a_n)·sin φ_n))/2`, each from its exact gap to 1.
fn sn_cn_core<V: Work>(arg: V, k: V, kc: V) -> (V, V) {
    const DEPTH: usize = 32;
    let mut arith = [V::zero(); DEPTH];
    let mut geo = [V::zero(); DEPTH];
    let mut half_gap = [V::zero(); DEPTH];
    let (mut an, mut bn, mut cn) = (V::one(), kc, k);
    let mut depth_used = 0;
    for i in 0..DEPTH {
        if let (Some(x), Some(y), Some(z)) = (arith.get_mut(i), geo.get_mut(i), half_gap.get_mut(i))
        {
            (*x, *y, *z) = (an, bn, cn);
        }
        depth_used = i;
        if V::negligible(&cn, &an) || i + 1 == DEPTH {
            break;
        }
        let next = (an + bn).ldexp(-1);
        (an, bn, cn) = (next, (an * bn).sqrt(), cn.square().ldexp(-2) / next);
    }
    let mut phi = an.ldexp(i32::try_from(depth_used).unwrap_or(0)) * arg;
    for i in (1..=depth_used).rev() {
        let at = |table: &[V; DEPTH], idx: usize| table.get(idx).copied().unwrap_or_else(V::zero);
        let (sine, co) = sin_cos(phi);
        let ratio = at(&half_gap, i) / at(&arith, i);
        let one = V::one();
        let gap = at(&geo, i - 1) / at(&arith, i) + ratio * co.square() / (one + sine.abs());
        let rest = one - gap;
        let angle = atan_ratio(rest.abs(), (gap * (one + rest)).sqrt()).with_sign(sine.is_neg());
        phi = (phi + angle).ldexp(-1);
    }
    sin_cos(phi)
}

/// `(sn, cn, dn, am)` of `u` for the modulus parts.
fn jacobi_work<V: Work>(u: V, parts_: &Parts<V>) -> (V, V, V, V) {
    let one = V::one();
    if parts_.kc.is_zero() {
        // The separatrix: sn = tanh u, cn = dn = sech u, am = atan(sinh u).
        if u.abs() > V::from_int(crate::eval::HUGE) {
            let edge = one.with_sign(u.is_neg());
            return (
                edge,
                V::zero(),
                V::zero(),
                half_pi::<V>().with_sign(u.is_neg()),
            );
        }
        let (sh, ch) = sinh_cosh(u.abs());
        let sech = one / ch;
        let am = atan_ratio(sh, one).with_sign(u.is_neg());
        return ((sh / ch).with_sign(u.is_neg()), sech, sech, am);
    }
    let kk = parts_.m.sqrt();
    let quarter = complete_k(parts_.kc);
    let period = quarter.ldexp(1);
    // u = (2^124·top + 2^62·high + low)·2K + rest, in parts so a 128-bit u fits.
    let cycles = u / period;
    let top = cycles.ldexp(-124).round_int();
    let cycles = cycles - V::from_int(top).ldexp(124);
    let high = cycles.ldexp(-62).round_int();
    let low = (cycles - V::from_int(high).ldexp(62)).round_int();
    let whole =
        |unit: V| unit.mul_int(top).ldexp(124) + unit.mul_int(high).ldexp(62) + unit.mul_int(low);
    let rest = u - whole(period);
    let dist = rest.abs();
    let (sn, cn, dn) = if dist <= quarter.ldexp(-1) {
        let (sine, cosine) = sn_cn_core(dist, kk, parts_.kc);
        (
            sine,
            cosine,
            (cosine.square() + parts_.kc2 * sine.square()).sqrt(),
        )
    } else {
        // sn(K − w) = cn/dn, cn(K − w) = kk'·sn/dn and dn(K − w) = kk'/dn.
        let (sine, cosine) = sn_cn_core(quarter - dist, kk, parts_.kc);
        let delta = (cosine.square() + parts_.kc2 * sine.square()).sqrt();
        (cosine / delta, parts_.kc * sine / delta, parts_.kc / delta)
    };
    let am = atan_ratio(sn, cn.abs());
    let am = if cn.is_neg() { pi::<V>() - am } else { am };
    let flip = low % 2 != 0;
    let sign = |x: V, neg: bool| if neg { -x } else { x };
    let winding = whole(pi::<V>());
    (
        sign(sign(sn, rest.is_neg()), flip),
        sign(cn, flip),
        dn,
        sign(am, rest.is_neg()) + winding,
    )
}

/// The Jacobi elliptic functions `(sn, cn, dn)` of `u`, any `u`.
#[must_use]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn jacobi_sn_cn_dn<T: CordicNumber>(u: T, modulus: Modulus<T>) -> (T, T, T) {
    let (sn, cn, dn, _) = jacobi_work(work(u), &parts(modulus));
    (fixed(sn), fixed(cn), fixed(dn))
}

/// The Jacobi amplitude `am(u)`, continuous: `am(u + 2K) = am(u) + π`.
#[must_use]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn jacobi_am<T: CordicNumber>(u: T, modulus: Modulus<T>) -> T {
    fixed(jacobi_work(work(u), &parts(modulus)).3)
}

/// `R_F(cn², dn², 1)`, finite unless `cn` and `dn` are both 0.
fn rf_sq<V: Work>(cn: V, dn: V) -> Option<V> {
    if cn.is_zero() && dn.is_zero() {
        return None;
    }
    Some(rf(
        cn.square(),
        dn.square(),
        V::one(),
        Some([cn, dn, V::one()]),
    ))
}

/// The `u` in `(−2K, 2K]` with these [`jacobi_sn_cn_dn`] values, all given so none is derived.
/// # Errors
/// Returns `DomainError` when `cn` and `dn` are both 0, where `u` is infinite.
#[must_use = "returns the argument which should be handled"]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn jacobi_inverse<T: CordicNumber>(sn: T, cn: T, dn: T, modulus: Modulus<T>) -> Result<T> {
    let (sw, cw, dw) = (work(sn), work(cn), work(dn));
    let Some(rf_value) = rf_sq(cw, dw) else {
        return Err(Error::domain("jacobi_inverse", "cn or dn non-zero"));
    };
    let part = sw.abs() * rf_value;
    let arg = if cw.is_neg() {
        let parts_ = parts(modulus);
        if parts_.kc.is_zero() {
            crate::eval::beyond()
        } else {
            complete_k(parts_.kc).ldexp(1) - part
        }
    } else {
        part
    };
    Ok(fixed(arg.with_sign(sw.is_neg())))
}

/// The complete elliptic integral of the first kind `K(m)`.
/// # Errors
/// Returns `DomainError` at `m = 1`, where `K` is infinite.
#[must_use = "returns the integral which should be handled"]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn ellip_k<T: CordicNumber>(modulus: Modulus<T>) -> Result<T> {
    let p = parts(modulus);
    if p.kc.is_zero() {
        return Err(Error::domain("ellip_k", "m below 1"));
    }
    Ok(fixed(complete_k(p.kc)))
}

/// The complete elliptic integral of the second kind `E(m)`, 1 at `m = 1`.
#[must_use]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn ellip_e<T: CordicNumber>(modulus: Modulus<T>) -> T {
    let p = parts(modulus);
    if p.kc.is_zero() {
        return T::one();
    }
    let one = W::<T>::one();
    let zero = W::<T>::zero();
    let f = rf(zero, p.kc2, one, Some([zero, p.kc, one]));
    fixed(f - (p.m * rd(zero, p.kc2, one)).div_int(3))
}

/// `Π(n; φ|m) = ∫₀^φ dθ/((1 − n·sin²θ)√(1 − m·sin²θ))`, `φ` in `(−π, π]` by `(sn, cn, dn)`.
/// # Errors
/// Returns `DomainError` when `1 − n·sn² ≤ 0`, or `cn` and `dn` are both 0.
#[must_use = "returns the integral which should be handled"]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn ellip_pi<T: CordicNumber>(n: T, sn: T, cn: T, dn: T, modulus: Modulus<T>) -> Result<T> {
    let (nw, sw, cw, dw) = (work(n), work(sn), work(cn), work(dn));
    let one = W::<T>::one();
    let s2 = sw.square();
    let pole = one - nw * s2;
    // Past π/2 the integral crosses its pole unless n < 1.
    let reachable = !cw.is_neg() || nw < one;
    let Some(rf_value) = rf_sq(cw, dw).filter(|_| !pole.is_neg() && !pole.is_zero() && reachable)
    else {
        return Err(Error::domain(
            "ellip_pi",
            "1 - n sn^2 > 0, n < 1 past pi/2, cn or dn non-zero",
        ));
    };
    let rj_value = rj(
        cw.square(),
        dw.square(),
        one,
        pole,
        Some([cw, dw, one, pole.sqrt()]),
    );
    let part = sw.abs() * (rf_value + (nw * s2 * rj_value).div_int(3));
    let total = if cw.is_neg() {
        let parts_ = parts(modulus);
        if parts_.kc.is_zero() {
            crate::eval::beyond()
        } else {
            complete_pi(nw, &parts_).ldexp(1) - part
        }
    } else {
        part
    };
    Ok(fixed(total.with_sign(sw.is_neg())))
}

fn complete_pi<V: Work>(n: V, p: &Parts<V>) -> V {
    let (one, zero) = (V::one(), V::zero());
    let f = rf(zero, p.kc2, one, Some([zero, p.kc, one]));
    let pole = one - n;
    f + (n * rj(zero, p.kc2, one, pole, Some([zero, p.kc, one, pole.sqrt()]))).div_int(3)
}

/// The complete elliptic integral of the third kind `Π(n|m)`.
/// # Errors
/// Returns `DomainError` for `n ≥ 1`, or at `m = 1`.
#[must_use = "returns the integral which should be handled"]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn ellip_pi_complete<T: CordicNumber>(n: T, modulus: Modulus<T>) -> Result<T> {
    let p = parts(modulus);
    if n >= T::one() || p.kc.is_zero() {
        return Err(Error::domain(
            "ellip_pi_complete",
            "n below 1 and m below 1",
        ));
    }
    Ok(fixed(complete_pi(work(n), &p)))
}

/// Carlson's `R_F(x, y, z)`. Domain: `x, y, z ≥ 0`, at most one of them 0.
/// # Errors
/// Returns `DomainError` outside the domain.
#[must_use = "returns the integral which should be handled"]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn carlson_rf<T: CordicNumber>(x: T, y: T, z: T) -> Result<T> {
    let v = [work(x), work(y), work(z)];
    check_rf(&v, "carlson_rf")?;
    Ok(fixed(rf(v[0], v[1], v[2], None)))
}

/// [`carlson_rf`] of the squares of `rx`, `ry` and `rz`, for arguments too small to hold.
/// # Errors
/// Returns `DomainError` if two of the roots are 0.
#[must_use = "returns the integral which should be handled"]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn carlson_rf_roots<T: CordicNumber>(rx: T, ry: T, rz: T) -> Result<T> {
    let r = [work(rx), work(ry), work(rz)];
    let v = r.map(Work::square);
    check_rf(&v, "carlson_rf_roots")?;
    Ok(fixed(rf(v[0], v[1], v[2], Some(r))))
}

fn check_rf<V: Work>(v: &[V; 3], name: &'static str) -> Result<()> {
    let zeros = v.iter().filter(|x| x.is_zero()).count();
    if v.iter().any(Work::is_neg) || zeros > 1 {
        return Err(Error::domain(name, "x, y, z >= 0 with at most one zero"));
    }
    Ok(())
}

/// Carlson's `R_J(x, y, z, p)`. Domain: `x, y, z ≥ 0`, at most one 0, and `p > 0`.
/// # Errors
/// Returns `DomainError` outside the domain.
#[must_use = "returns the integral which should be handled"]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn carlson_rj<T: CordicNumber>(x: T, y: T, z: T, p: T) -> Result<T> {
    let args = [work(x), work(y), work(z)];
    check_rf(&args, "carlson_rj")?;
    if p <= T::zero() {
        return Err(Error::domain("carlson_rj", "p > 0"));
    }
    Ok(fixed(rj(args[0], args[1], args[2], work(p), None)))
}

/// [`carlson_rj`] of the squares of `rx`, `ry`, `rz` and `rp`.
/// # Errors
/// Returns `DomainError` if two of `rx`, `ry`, `rz` are 0, or `rp` is.
#[must_use = "returns the integral which should be handled"]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn carlson_rj_roots<T: CordicNumber>(rx: T, ry: T, rz: T, rp: T) -> Result<T> {
    let roots = [work(rx), work(ry), work(rz), work(rp)];
    let squares = [roots[0].square(), roots[1].square(), roots[2].square()];
    check_rf(&squares, "carlson_rj_roots")?;
    if roots[3].is_zero() {
        return Err(Error::domain("carlson_rj_roots", "rp non-zero"));
    }
    Ok(fixed(rj(
        squares[0],
        squares[1],
        squares[2],
        roots[3].square(),
        Some(roots),
    )))
}

/// Carlson's `R_C(x, y)`. Domain: `x ≥ 0`, `y > 0`.
/// # Errors
/// Returns `DomainError` outside the domain.
#[must_use = "returns the integral which should be handled"]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn carlson_rc<T: CordicNumber>(x: T, y: T) -> Result<T> {
    if x < T::zero() || y <= T::zero() {
        return Err(Error::domain("carlson_rc", "x >= 0, y > 0"));
    }
    Ok(fixed(rc(work(x), work(y))))
}

/// Carlson's `R_D(x, y, z)`. Domain: `x, y ≥ 0`, not both 0, and `z > 0`.
/// # Errors
/// Returns `DomainError` outside the domain.
#[must_use = "returns the integral which should be handled"]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn carlson_rd<T: CordicNumber>(x: T, y: T, z: T) -> Result<T> {
    let zero = T::zero();
    if x < zero || y < zero || z <= zero || (x == zero && y == zero) {
        return Err(Error::domain(
            "carlson_rd",
            "x, y >= 0 not both zero, z > 0",
        ));
    }
    Ok(fixed(rd(work(x), work(y), work(z))))
}
