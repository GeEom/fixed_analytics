//! Elliptic integrals and Jacobi elliptic functions of an exactly held [`Modulus`].

use crate::error::{Error, Result};
use crate::eval::{
    W, atan_ratio, by_width, fixed, half_pi, log1p, pi, sin_cos_octant, sinh_cosh, tiered, work,
};
use crate::real::{Work, mul_q126};
use crate::series::ODD;
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
    #[cfg_attr(
        all(feature = "verify-no-panic", not(debug_assertions)),
        no_panic::no_panic
    )]
    pub fn from_parameter(m: T) -> Option<Self> {
        (m >= T::zero() && m <= T::one()).then_some(Self {
            value: m,
            complement: false,
        })
    }

    /// From the complement `k' = √(1 − m)`, in `[0, 1]`, exact where `m` is near 1.
    #[must_use]
    #[cfg_attr(
        all(feature = "verify-no-panic", not(debug_assertions)),
        no_panic::no_panic
    )]
    pub fn from_complement(kc: T) -> Option<Self> {
        (kc >= T::zero() && kc <= T::one()).then_some(Self {
            value: kc,
            complement: true,
        })
    }

    fn at_separatrix(self) -> bool {
        if self.complement {
            self.value == T::zero()
        } else {
            self.value == T::one()
        }
    }

    /// The parameter `m`, rounded.
    #[must_use]
    #[cfg_attr(
        all(feature = "verify-no-panic", not(debug_assertions)),
        no_panic::no_panic
    )]
    pub fn parameter(self) -> T {
        fixed(parts(self).m)
    }

    /// The complement `k'`, rounded.
    #[must_use]
    #[cfg_attr(
        all(feature = "verify-no-panic", not(debug_assertions)),
        no_panic::no_panic
    )]
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
    parts_in(modulus, work)
}

fn parts_in<T: CordicNumber, V: Work>(modulus: Modulus<T>, as_work: fn(T) -> V) -> Parts<V> {
    let one = V::one();
    let v = as_work(modulus.value);
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
fn rc_one<V: Work>(e: V, one_plus_e: impl FnOnce() -> V) -> V {
    // Σ (−e)^k/(2k + 1), as many terms as e's size needs below 2^−4.
    let size = -e.exponent();
    if size >= 4 {
        let terms = V::PREC / size.cast_unsigned() + 1;
        return V::horner(e, &ODD, 0, 1, terms, true);
    }
    let t = e.abs().sqrt();
    if t.is_zero() {
        return V::one();
    }
    if e.is_neg() {
        // 1 − t = (1 + e)/(1 + t), so atanh t = log1p(2t(1 + t)/(1 + e))/2.
        let one = V::one();
        log1p((t * (one + t) / one_plus_e()).ldexp(1)).ldexp(-1) / t
    } else {
        atan_ratio(t, V::one()) / t
    }
}

/// `R_C(x, y)` for `x ≥ 0` and `y > 0`, by ratios so `x` near `y` does not cancel.
fn rc<V: Work>(x: V, y: V) -> V {
    if x.is_zero() {
        return half_pi::<V>() / y.sqrt();
    }
    rc_one((y - x) / x, || y / x) / x.sqrt()
}

/// Bits by which deviations up to `spread·4^−steps/mean` are surely below 1.
fn closeness<V: Work>(spread: V, steps: i32, mean: V) -> i32 {
    let gap = spread.ldexp(-2 * steps);
    if gap.is_zero() {
        return V::PREC.cast_signed() + 4;
    }
    let bits = mean.exponent() - gap.exponent() - 1;
    if gap.ldexp(bits + 1) < mean.abs() {
        bits + 1
    } else {
        bits
    }
}

/// `2^126/(2n(2n + c))`, the closing series' divisors in fixed point.
const fn quotients(c: u128) -> [i128; 40] {
    let mut out = [0; 40];
    let mut rest: &mut [i128] = &mut out;
    let mut n = 0;
    while let [slot, tail @ ..] = rest {
        if n > 0 {
            *slot = ((1 << 126) / (2 * n * (2 * n + c))).cast_signed();
        }
        rest = tail;
        n += 1;
    }
    out
}

const QUOTIENTS_F: [i128; 40] = quotients(1);
const QUOTIENTS_J: [i128; 40] = quotients(3);

/// The closing series `Σ c·T_N/(2N + c)`, `T_N` the `t^N` term of `Π(1 − Z·t)^(−1/2)` by its
/// recurrence in `e₂…e₅`, `c` 1 for `R_F` else 3, to degree `(PREC + 4)/bits`; scaled by
/// `2^(N(bits − 3))` the terms stay below 1, so most run in fixed point.
fn closing<V: Work>([e2, e3, e4, e5]: [V; 4], c: i64, bits: i32) -> V {
    let degree = i64::from((V::PREC.cast_signed() + 3 + bits) / bits.max(1) - 1);
    let h = [-e2, e3, -e4, e5];
    let shift = i64::from(bits - 3);
    let fixed = if shift > 0 && degree < 40 {
        ((i64::from(V::PREC) - 118 + shift - 1) / shift).max(2)
    } else {
        degree + 1
    };
    let mut last = [V::zero(), V::zero(), V::zero(), V::zero(), V::one()];
    let mut sum = V::one();
    let coefficient = |n: i64, k: i64| (2 * n - k) * (2 * (n - k) + c);
    for n in 1..fixed.min(degree + 1) {
        let mut acc = V::zero();
        for ((k, &hk), &term) in (2..).zip(&h).zip(last.iter().rev().skip(1)) {
            acc = acc + (hk * term).mul_int(coefficient(n, k));
        }
        let next = acc.div_int((2 * n * (2 * n + c)).unsigned_abs());
        let [_, l1, l2, l3, l4] = last;
        last = [l1, l2, l3, l4, next];
        sum = sum + next;
    }
    if fixed > degree {
        return sum;
    }
    let raw = |v: V, n: i64| {
        let frac = u32::try_from(126 + n * shift).unwrap_or(0);
        v.to_raw(frac, i128::MIN, i128::MAX, false)
    };
    let hq = [raw(h[0], 2), raw(h[1], 3), raw(h[2], 4), raw(h[3], 5)];
    let mut window = [0; 5];
    for ((slot, &v), n) in window.iter_mut().zip(&last).zip(fixed - 5..) {
        *slot = raw(v, n);
    }
    let quotients = if c == 1 { &QUOTIENTS_F } else { &QUOTIENTS_J };
    let mut tail = 0i128;
    for n in fixed..=degree {
        let mut acc = 0i128;
        for ((k, &hk), &term) in (2..).zip(&hq).zip(window.iter().rev().skip(1)) {
            if hk != 0 && term != 0 {
                acc += mul_q126(hk, term) * i128::from(coefficient(n, k));
            }
        }
        let quotient = usize::try_from(n).ok().and_then(|i| quotients.get(i));
        let next = mul_q126(acc, quotient.copied().unwrap_or(0));
        let drop = u32::try_from((n - fixed) * shift).unwrap_or(u32::MAX);
        tail += next.checked_shr(drop).unwrap_or(0);
        let [_, w1, w2, w3, w4] = window;
        window = [w1, w2, w3, w4, next];
    }
    sum + V::from_raw(tail, u32::try_from(126 + fixed * shift).unwrap_or(0))
}

/// An integral's mean `A_m`, from `A₀`, and its arguments' widest spread about `A₀`.
#[derive(Clone, Copy)]
struct Mean<V> {
    a0: V,
    now: V,
    spread: V,
}

impl<V: Work> Mean<V> {
    fn new(a0: V, args: &[V]) -> Self {
        let spread = args.iter().fold(V::zero(), |m, &v| m.max((a0 - v).abs()));
        Self {
            a0,
            now: a0,
            spread,
        }
    }

    fn bits(&self, steps: i32) -> i32 {
        closeness(self.spread, steps, self.now)
    }

    fn deviation(&self, v: V, scale: V) -> V {
        (self.a0 - v) * scale
    }

    /// `R_F` after `steps` steps, from `inv = 1/A_m`.
    fn rf(&self, inv: V, steps: i32, [x0, y0]: [V; 2]) -> V {
        let scale = inv.ldexp(-2 * steps);
        let (dx, dy) = (self.deviation(x0, scale), self.deviation(y0, scale));
        let dz = -(dx + dy);
        let xy = dx * dy;
        let e = [xy - dz.square(), xy * dz, V::zero(), V::zero()];
        closing(e, 1, self.bits(steps)) * inv.sqrt()
    }

    /// `R_D` past `steps` steps: `4^−m·A_m^(−3/2)` times its closing series.
    fn rd_tail(&self, inv: V, steps: i32, [x0, y0]: [V; 2]) -> V {
        let scale = inv.ldexp(-2 * steps);
        let (dx, dy) = (self.deviation(x0, scale), self.deviation(y0, scale));
        let dz = -(dx + dy).div_int(3);
        let xy = dx * dy;
        let z2 = dz.square();
        let e = [
            xy - z2.mul_int(6),
            (xy.mul_int(3) - z2.mul_int(8)) * dz,
            (xy - z2) * z2.mul_int(3),
            xy * z2 * dz,
        ];
        (closing(e, 3, self.bits(steps)) * inv * inv.sqrt()).ldexp(-2 * steps)
    }

    fn rj_tail(&self, inv: V, steps: i32, [x0, y0, z0]: [V; 3]) -> V {
        let scale = inv.ldexp(-2 * steps);
        let (dx, dy, dz) = (
            self.deviation(x0, scale),
            self.deviation(y0, scale),
            self.deviation(z0, scale),
        );
        let dp = -(dx + dy + dz).ldexp(-1);
        let e2 = dx * dy + dx * dz + dy * dz - dp.square().mul_int(3);
        let xyz = dx * dy * dz;
        let p2 = dp.square();
        let e = [
            e2,
            xyz + (e2 * dp).ldexp(1) + (p2 * dp).ldexp(2),
            (xyz.ldexp(1) + e2 * dp + (p2 * dp).mul_int(3)) * dp,
            xyz * p2,
        ];
        (closing(e, 3, self.bits(steps)) * inv * inv.sqrt()).ldexp(-2 * steps)
    }
}

/// Each value's reciprocal in its place, all from one division.
fn invert<V: Work>(values: &mut [V]) {
    if let [only] = values {
        *only = V::one() / *only;
        return;
    }
    let mut before = [V::one(); 32];
    let mut all = V::one();
    for (i, (b, &v)) in before.iter_mut().zip(values.iter()).enumerate() {
        *b = all;
        all = if i == 0 { v } else { all * v };
    }
    let mut inv = V::one() / all;
    for (i, (v, &b)) in values.iter_mut().zip(&before).enumerate().rev() {
        if i == 0 {
            *v = inv;
        } else {
            let r = inv * b;
            inv = inv * *v;
            *v = r;
        }
    }
}

/// `R_F` and `R_D` as `want` asks (else 0), and `R_J` given `p`, from one duplication.
fn duplicate<V: Work>(
    [x0, y0, z0]: [V; 3],
    p: Option<V>,
    roots: Option<[V; 4]>,
    [want_f, want_d]: [bool; 2],
) -> [V; 3] {
    let (zero, one) = (V::zero(), V::one());
    let xyz = x0 + y0 + z0;
    let mut for_f = want_f.then(|| Mean::new(xyz.div_int(3), &[x0, y0, z0]));
    let mut for_d = want_d.then(|| Mean::new((xyz + z0.ldexp(1)).div_int(5), &[x0, y0, z0]));
    let mut for_j = p.map(|p0| Mean::new((xyz + p0.ldexp(1)).div_int(5), &[x0, y0, z0, p0]));
    let delta = p.map(|p0| (p0 - x0) * (p0 - y0) * (p0 - z0));
    let (mut vx, mut vy, mut vz, mut vp) = (x0, y0, z0, p.unwrap_or(zero));
    // R_D's Σ 4^−m/(√z(z + λ)) as num/den; R_J's d_m, inverted at the end with the means.
    let (mut num, mut den, mut sum) = (zero, one, zero);
    let mut values = [one; 32];
    let mut held = 0u32;
    let mut steps = 0;
    let enough = by_width::<V>([5, 6, 9]).cast_signed();
    while steps < by_width::<V>([15, 21, 28]).cast_signed() {
        if [for_f, for_d, for_j]
            .iter()
            .flatten()
            .all(|m| m.bits(steps) >= enough)
        {
            break;
        }
        let [sx, sy, sz, sp] = match roots {
            Some(r) if steps == 0 => r.map(Work::abs),
            _ => [
                vx.sqrt(),
                vy.sqrt(),
                vz.sqrt(),
                p.map_or(zero, |_| vp.sqrt()),
            ],
        };
        let lam = sx * sy + sx * sz + sy * sz;
        if want_d {
            let q = sz * (vz + lam);
            num = num * q + den.ldexp(-2 * steps);
            den = den * q;
        }
        if let Some(delta) = delta {
            let dm = (sp + sx) * (sp + sy) * (sp + sz);
            let dm_slot = usize::try_from(steps).ok().and_then(|i| values.get_mut(i));
            let at = delta.ldexp(-6 * steps);
            match dm_slot {
                Some(slot) if at.abs() < dm.square().ldexp(-4) => {
                    *slot = dm;
                    held |= 1 << steps;
                }
                _ => {
                    let inv = one / dm;
                    let ex = at * inv.square();
                    // 1 + ex = 2√vp(vp + λ)/d_m, as Π(√vp + √v) + Π(√vp − √v) = 2√vp(vp + λ).
                    let one_plus_e = || (sp * (vp + lam) * inv).ldexp(1);
                    sum = sum + (rc_one(ex, one_plus_e) * inv).ldexp(-2 * steps);
                }
            }
            vp = (vp + lam).ldexp(-2);
        }
        for m in [&mut for_f, &mut for_d, &mut for_j].into_iter().flatten() {
            m.now = (m.now + lam).ldexp(-2);
        }
        vx = (vx + lam).ldexp(-2);
        vy = (vy + lam).ldexp(-2);
        vz = (vz + lam).ldexp(-2);
        steps += 1;
    }
    let mut used = if delta.is_some() {
        usize::try_from(steps).unwrap_or(0)
    } else {
        0
    };
    let now = |m: Option<Mean<V>>| m.map(|m| m.now);
    let means = [now(for_f), now(for_d).map(|a| a * den), now(for_j)];
    let slots = means.map(|mean| {
        let at = used;
        let slot = values.get_mut(at)?;
        *slot = mean?;
        used += 1;
        Some(at)
    });
    if let Some(all) = values.get_mut(..used) {
        invert(all);
    }
    if let Some(delta) = delta {
        for (m, &inv) in (0..steps).zip(&values).filter(|&(m, _)| held >> m & 1 == 1) {
            let ex = delta.ldexp(-6 * m) * inv.square();
            sum = sum + (rc_one(ex, || one + ex) * inv).ldexp(-2 * m);
        }
    }
    let [inv_f, inv_d, inv_j] = slots.map(|at| at.and_then(|i| values.get(i).copied()));
    let rf = for_f
        .zip(inv_f)
        .map_or(zero, |(m, r)| m.rf(r, steps, [x0, y0]));
    let rd = for_d.zip(inv_d).map_or(zero, |(m, r)| {
        m.rd_tail(r * den, steps, [x0, y0]) + (num * r * m.now).mul_int(3)
    });
    let rj = for_j.zip(inv_j).map_or(zero, |(m, r)| {
        m.rj_tail(r, steps, [x0, y0, z0]) + sum.mul_int(6)
    });
    [rf, rd, rj]
}

/// [`duplicate`], but DLMF 19.21.12 trades a `p` far above the rest for `q` near the least,
/// `x`: `(p − x)R_J(p) + (q − x)R_J(q) = 3R_F − 3R_C(yz/x, pq/x)`.
fn carlson<V: Work>(args: [V; 3], p: Option<V>, roots: Option<[V; 4]>, want: [bool; 2]) -> [V; 3] {
    let Some(p0) = p else {
        return duplicate(args, None, roots, want);
    };
    let mut sorted = args;
    sorted.sort_unstable();
    let [lo, mid, hi] = sorted;
    if p0 <= hi.ldexp(4) {
        return duplicate(args, p, roots, want);
    }
    let gap = (mid - lo) * (hi - lo) / (p0 - lo);
    let near = lo + gap;
    let near_roots = roots.map(|[a, b, c, _]| [a, b, c, near.sqrt()]);
    let [rf, rd, near_rj] = duplicate(args, Some(near), near_roots, [true, want[1]]);
    let mut rest = rf.mul_int(3) - gap * near_rj;
    if !lo.is_zero() {
        let xi = mid * hi / lo;
        let excess = (p0 - mid) * (p0 - hi) / ((p0 - lo) * xi);
        rest = rest - (rc_one(excess, || p0 * near / (lo * xi)) / xi.sqrt()).mul_int(3);
    }
    [rf, rd, rest / (p0 - lo)]
}

/// `R_F(x, y, z)`, any `roots` standing in for the first square roots.
fn rf<V: Work>(x0: V, y0: V, z0: V, roots: Option<[V; 3]>) -> V {
    let roots = roots.map(|[a, b, c]| [a, b, c, V::zero()]);
    duplicate([x0, y0, z0], None, roots, [true, false])[0]
}

/// `R_J(x, y, z, p)` for `p > 0`.
fn rj<V: Work>(x0: V, y0: V, z0: V, p0: V, roots: Option<[V; 4]>) -> V {
    carlson([x0, y0, z0], Some(p0), roots, [false, false])[2]
}

/// `R_D(x, y, z)`.
fn rd<V: Work>(x0: V, y0: V, z0: V) -> V {
    duplicate([x0, y0, z0], None, None, [false, true])[1]
}

/// The arithmetic-geometric mean of 1 and `b`, for `b` in `(0, 1]`.
fn agm<V: Work>(b: V) -> V {
    let (mut a, mut b) = (V::one(), b);
    for _ in 0..32 {
        let mean = (a + b).ldexp(-1);
        if (a - b).abs().ldexp(V::PREC.cast_signed() / 2 + 2) <= a {
            return mean;
        }
        (a, b) = (mean, (a * b).sqrt());
    }
    a
}

/// `K(m) = π/(2·AGM(1, k'))`, for `k' > 0`.
fn complete_k<V: Work>(kc: V) -> V {
    half_pi::<V>() / agm(kc)
}

const DEPTH: usize = 16;

/// The AGM of 1 and `k′`, each level's means kept, `mean` its limit.
struct Means<V> {
    arith: [V; DEPTH],
    geo: [V; DEPTH],
    next: [V; DEPTH],
    used: usize,
    mean: V,
}

fn means<V: Work>(kc: V) -> Means<V> {
    let one = V::one();
    let mut out = Means {
        arith: [V::zero(); DEPTH],
        geo: [V::zero(); DEPTH],
        next: [V::zero(); DEPTH],
        used: 0,
        mean: one,
    };
    let (mut a, mut b) = (one, kc);
    let levels = out
        .arith
        .iter_mut()
        .zip(out.geo.iter_mut().zip(out.next.iter_mut()));
    for (i, (arith, (geo, next))) in levels.enumerate() {
        out.used = i + 1;
        out.mean = (a + b).ldexp(-1);
        (*arith, *geo, *next) = (a, b, out.mean);
        if (a - b).abs().ldexp(V::PREC.cast_signed() / 2 + 2) <= a {
            break;
        }
        b = (a * b).sqrt();
        a = out.mean;
    }
    out
}

/// `(sn, cn, dn)·h` and `h` for `w`, `c·w = (π/4)·frac`, by descending Landen steps.
fn landen<V: Work>(frac: V, means: &Means<V>) -> (V, V, V, V) {
    let one = V::one();
    let (sine, cosine) = sin_cos_octant(0, frac);
    let at = |t: &[V; DEPTH], i: usize| t.get(i).copied().unwrap_or(one);
    // 1 ± k·sn² as a_n·sn² or b_n·sn² plus a_{n+1}·cn²: sums, so small cn and dn keep their bits.
    let low = means.used.saturating_sub(1);
    let (mut s, mut c) = (at(&means.arith, low) * sine, at(&means.next, low) * cosine);
    let cos_part = c * cosine;
    let (mut d, mut h) = (
        at(&means.geo, low) * sine.square() + cos_part,
        s * sine + cos_part,
    );
    for i in (0..low).rev() {
        let (lead, tail) = (at(&means.arith, i) * s, at(&means.next, i) * c);
        let rest = tail * c;
        (s, c, d, h) = (
            lead * h,
            tail * d,
            at(&means.geo, i) * s.square() + rest,
            lead * s + rest,
        );
    }
    (s, c, d, h)
}

/// `u = (2^124·c₀ + 2^62·c₁ + c₂)·period + rest`, the count in parts so a 128-bit `u` fits.
fn reduce<V: Work>(u: V, period: V) -> ([i64; 3], V) {
    let cycles = u / period;
    let top = cycles.ldexp(-124).round_int();
    let cycles = cycles - V::from_int(top).ldexp(124);
    let high = cycles.ldexp(-62).round_int();
    let count = [
        top,
        high,
        (cycles - V::from_int(high).ldexp(62)).round_int(),
    ];
    (count, u - times(count, period))
}

/// `count·unit` for a count in [`reduce`]'s parts.
fn times<V: Work>([top, high, low]: [i64; 3], unit: V) -> V {
    unit.mul_int(top).ldexp(124) + unit.mul_int(high).ldexp(62) + unit.mul_int(low)
}

/// `(sn, cn, dn, c·u)` of `u` for the modulus parts; reducing `c·u` loses bits in proportion.
fn jacobi_work<V: Work>(u: V, parts_: &Parts<V>) -> (V, V, V, V) {
    let one = V::one();
    if parts_.kc.is_zero() {
        if u.abs() > V::from_int(crate::eval::HUGE) {
            return (one.with_sign(u.is_neg()), V::zero(), V::zero(), V::zero());
        }
        let (sh, ch) = sinh_cosh(u.abs());
        let sech = one / ch;
        return ((sh / ch).with_sign(u.is_neg()), sech, sech, V::zero());
    }
    let levels = means(parts_.kc);
    let scaled = u * levels.mean;
    let (oct, frac) = scaled.octant_of();
    let (sine, cosine, delta, den) = landen(frac, &levels);
    let (sn, cn, dn) = if (oct + 1) & 2 == 0 {
        let inv = one / den;
        (sine * inv, cosine * inv, delta * inv)
    } else {
        // sn(K ∓ w) = cn/dn, cn(K ∓ w) = ±k'·sn/dn and dn(K ∓ w) = k'/dn.
        let inv = one / delta;
        (cosine * inv, parts_.kc * sine * inv, parts_.kc * den * inv)
    };
    (
        sn.with_sign((oct & 4 != 0) != scaled.is_neg()),
        cn.with_sign((oct + 2) & 4 != 0),
        dn,
        scaled,
    )
}

/// The Jacobi elliptic functions `(sn, cn, dn)` of `u`, any `u`.
#[must_use]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn jacobi_sn_cn_dn<T: CordicNumber>(u: T, modulus: Modulus<T>) -> (T, T, T) {
    let ((sn, cn), (dn, _)) = tiered!(T, |v| {
        let (sn, cn, dn, scaled) = jacobi_work(v(u), &parts_in(modulus, v));
        ((sn, cn), (dn, scaled))
    });
    (sn, cn, dn)
}

/// The Jacobi amplitude `am(u)`, continuous: `am(u + 2K) = am(u) + π`.
#[must_use]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn jacobi_am<T: CordicNumber>(u: T, modulus: Modulus<T>) -> T {
    tiered!(T, |v| {
        let (sn, cn, _, scaled) = jacobi_work(v(u), &parts_in(modulus, v));
        amplitude(sn, cn, scaled)
    })
}

/// `nπ` nearest `c·u` with `(−1)^n·cn ≥ 0`, plus `atan(sn/cn)`.
fn amplitude<V: Work>(sn: V, cn: V, scaled: V) -> V {
    let (mut count, rest) = reduce(scaled, pi::<V>());
    if cn.is_neg() == (count[2] % 2 == 0) && !cn.is_zero() {
        count[2] += if rest.is_neg() { -1 } else { 1 };
    }
    let lean = atan_ratio(sn.abs(), cn.abs()).with_sign(sn.is_neg() != (count[2] % 2 != 0));
    times(count, pi::<V>()) + lean
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
    if cn == T::zero() && dn == T::zero() {
        return Err(Error::domain("jacobi_inverse", "cn or dn non-zero"));
    }
    Ok(tiered!(T, |v| inverse_work(
        v(sn),
        v(cn),
        v(dn),
        modulus,
        v
    )))
}

fn inverse_work<T: CordicNumber, V: Work>(
    sw: V,
    cw: V,
    dw: V,
    modulus: Modulus<T>,
    as_work: fn(T) -> V,
) -> V {
    let part = sw.abs() * rf_sq(cw, dw).unwrap_or_else(V::zero);
    let arg = if !cw.is_neg() {
        part
    } else if modulus.at_separatrix() {
        crate::eval::beyond()
    } else {
        complete_k(parts_in(modulus, as_work).kc).ldexp(1) - part
    };
    arg.with_sign(sw.is_neg())
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
    if modulus.at_separatrix() {
        return Err(Error::domain("ellip_k", "m below 1"));
    }
    Ok(tiered!(T, |v| complete_k(parts_in(modulus, v).kc)))
}

/// The complete elliptic integral of the second kind `E(m)`, 1 at `m = 1`.
#[must_use]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn ellip_e<T: CordicNumber>(modulus: Modulus<T>) -> T {
    if modulus.at_separatrix() {
        return T::one();
    }
    tiered!(T, |v| complete_e(&parts_in(modulus, v)))
}

fn complete_e<V: Work>(p: &Parts<V>) -> V {
    let (one, zero) = (V::one(), V::zero());
    let roots = Some([zero, p.kc, one, zero]);
    let [f, d, _] = carlson([zero, p.kc2, one], None, roots, [true, true]);
    f - (p.m * d).div_int(3)
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
    let (nw, sw) = (work(n), work(sn));
    let one = W::<T>::one();
    let pole = one - nw * sw.square();
    // Past π/2 the integral crosses its pole unless n < 1.
    let reachable = !cn.is_negative() || nw < one;
    if (cn == T::zero() && dn == T::zero()) || pole.is_neg() || pole.is_zero() || !reachable {
        return Err(Error::domain(
            "ellip_pi",
            "1 - n sn^2 > 0, n < 1 past pi/2, cn or dn non-zero",
        ));
    }
    Ok(tiered!(T, |v| {
        let (s_v, n_v) = (v(sn), v(n));
        let gap = v(T::one()) - n_v * s_v.square();
        pi_work(n_v, s_v, v(cn), v(dn), gap, &parts_in(modulus, v))
    }))
}

/// `Π(n; φ|m)` from `φ`'s `(sn, cn, dn)` and `1 − n·sn²`.
fn pi_work<V: Work>(nw: V, sw: V, cw: V, dw: V, pole: V, parts_: &Parts<V>) -> V {
    let one = V::one();
    let s2 = sw.square();
    let (roots, xyz) = ([cw, dw, one, pole.sqrt()], [cw.square(), dw.square(), one]);
    let [f, _, rj_value] = carlson(xyz, Some(pole), Some(roots), [true, false]);
    let both_zero = cw.is_zero() && dw.is_zero();
    let rf_value = if both_zero { V::zero() } else { f };
    let part = sw.abs() * (rf_value + (nw * s2 * rj_value).div_int(3));
    let total = if cw.is_neg() {
        if parts_.kc.is_zero() {
            crate::eval::beyond()
        } else {
            complete_pi(nw, parts_).ldexp(1) - part
        }
    } else {
        part
    };
    total.with_sign(sw.is_neg())
}

fn complete_pi<V: Work>(n: V, p: &Parts<V>) -> V {
    let (one, zero) = (V::one(), V::zero());
    let pole = one - n;
    let roots = Some([zero, p.kc, one, pole.sqrt()]);
    let [f, _, j] = carlson([zero, p.kc2, one], Some(pole), roots, [true, false]);
    f + (n * j).div_int(3)
}

/// `∫ dv/(1 − n·sn²v)` from `u` to `u + step`, the change in `Π(n; am v|m)`, from the ends
/// reduced to `[−K, K]` with `2Π(n|m)` a period, so a small step keeps its bits.
/// # Errors
/// Returns `DomainError` unless `n < 1`.
#[must_use = "returns the integral which should be handled"]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn ellip_pi_increment<T: CordicNumber>(n: T, u: T, step: T, modulus: Modulus<T>) -> Result<T> {
    if n >= T::one() {
        return Err(Error::domain("ellip_pi_increment", "n below 1"));
    }
    let (rise, _) = tiered!(T, |v| increment_work(
        v(n),
        v(u),
        v(step),
        &parts_in(modulus, v)
    ));
    Ok(rise)
}

/// The increment, and the size its error scales with: the ends' integrals, `u`, and `K/(1 − n·sn²)`.
fn increment_work<V: Work>(n: V, u: V, s: V, p: &Parts<V>) -> (V, V) {
    let one = V::one();
    let local = |v: V| {
        let (sn, cn, dn, _) = jacobi_work(v, p);
        let gap = cn.square() + (one - n) * sn.square();
        (pi_work(n, sn, cn, dn, gap, p), one / gap)
    };
    let reach = u.abs().max((u + s).abs());
    let ((from, near), (to, far), whole, scale) = if p.kc.is_zero() {
        // Past the cap sech² v is below the working precision and the integrand 1/(1 − n).
        let cap = V::from_int(i64::from(by_width::<V>([29, 59, 88])));
        if let Some(rise) = gauss(n, u, s, p) {
            return (rise, rise.abs().max(reach));
        }
        let capped = |v: V| {
            let edge = cap.with_sign(v.is_neg());
            if v.abs() > cap {
                let (at, steep) = local(edge);
                (at + (v - edge) / (one - n), steep)
            } else {
                local(v)
            }
        };
        (capped(u), capped(u + s), V::zero(), V::zero())
    } else {
        let quarter = complete_k(p.kc);
        let period = quarter.ldexp(1);
        let (_, start) = reduce(u, period);
        if let Some(rise) = gauss(n, start, s, p) {
            return (rise, rise.abs().max(reach));
        }
        let (count, end) = reduce(start + s, period);
        let whole = if count == [0; 3] {
            V::zero()
        } else {
            times(count, complete_pi(n, p).ldexp(1))
        };
        (local(start), local(end), whole, quarter)
    };
    let size = from.abs().max(to.abs()).max(whole.abs()).max(reach);
    (to - from + whole, size.max(scale * near.max(far)))
}

/// The increment by five-point Gauss–Legendre where the step is too short for a difference not to
/// cancel, so far inside the nearest zero of `1 − n·sn²`, by two derivatives, that it is exact.
fn gauss<V: Work>(n: V, from: V, step: V, p: &Parts<V>) -> Option<V> {
    let one = V::one();
    let half = step.ldexp(-1);
    let mid = from + half;
    let (sn, cn, dn, _) = jacobi_work(mid, p);
    let gap = cn.square() + (one - n) * sn.square();
    let slope = (n * sn * cn * dn).ldexp(1);
    let bend = (n * ((cn * dn).square() - (sn * dn).square() - p.m * (sn * cn).square())).ldexp(1);
    let room = (gap / slope).abs().min((gap / bend).abs().sqrt()).min(one);
    if step
        .abs()
        .ldexp(i32::try_from(by_width::<V>([7, 14, 20])).unwrap_or(0))
        > room
    {
        return None;
    }
    let at = |v: V| {
        let (sine, cosine, _, _) = jacobi_work(v, p);
        one / (cosine.square() + (one - n) * sine.square())
    };
    let root = V::from_int(10).div_int(7).sqrt().ldexp(1);
    let (inner, outer) = (
        (V::from_int(5) - root).sqrt().div_int(3),
        (V::from_int(5) + root).sqrt().div_int(3),
    );
    let lean = V::from_int(70).sqrt().mul_int(13);
    let pair = |x: V| at(mid - half * x) + at(mid + half * x);
    let sum = V::from_int(128).div_int(225) / gap
        + (V::from_int(322) + lean).div_int(900) * pair(inner)
        + (V::from_int(322) - lean).div_int(900) * pair(outer);
    Some(half * sum)
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
    if n >= T::one() || modulus.at_separatrix() {
        return Err(Error::domain(
            "ellip_pi_complete",
            "n below 1 and m below 1",
        ));
    }
    Ok(tiered!(T, |v| complete_pi(v(n), &parts_in(modulus, v))))
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
    check_rf(&[x, y, z], "carlson_rf")?;
    Ok(tiered!(T, |v| rf(v(x), v(y), v(z), None)))
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
    let zeros = [rx, ry, rz].iter().filter(|&&r| r == T::zero()).count();
    if zeros > 1 {
        return Err(Error::domain(
            "carlson_rf_roots",
            "x, y, z >= 0 with at most one zero",
        ));
    }
    Ok(tiered!(T, |v| {
        let r = [v(rx), v(ry), v(rz)];
        let [a, b, c] = r.map(Work::square);
        rf(a, b, c, Some(r))
    }))
}

fn check_rf<T: CordicNumber>(v: &[T; 3], name: &'static str) -> Result<()> {
    let zeros = v.iter().filter(|&&x| x == T::zero()).count();
    if v.iter().any(|x| x.is_negative()) || zeros > 1 {
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
    check_rf(&[x, y, z], "carlson_rj")?;
    if p <= T::zero() {
        return Err(Error::domain("carlson_rj", "p > 0"));
    }
    Ok(tiered!(T, |v| rj(v(x), v(y), v(z), v(p), None)))
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
    let zeros = [rx, ry, rz].iter().filter(|&&r| r == T::zero()).count();
    if zeros > 1 {
        return Err(Error::domain(
            "carlson_rj_roots",
            "x, y, z >= 0 with at most one zero",
        ));
    }
    if rp == T::zero() {
        return Err(Error::domain("carlson_rj_roots", "rp non-zero"));
    }
    Ok(tiered!(T, |v| {
        let roots = [v(rx), v(ry), v(rz), v(rp)];
        let [a, b, c, d] = roots.map(Work::square);
        rj(a, b, c, d, Some(roots))
    }))
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
    Ok(tiered!(T, |v| rc(v(x), v(y))))
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
    Ok(tiered!(T, |v| rd(v(x), v(y), v(z))))
}

/// [`carlson`] of the squares of `r` and `rp`, its first square roots exact.
fn squares<V: Work>(r: [V; 3], rp: Option<V>, want: [bool; 2]) -> [V; 3] {
    let [a, b, c] = r;
    let roots = [a, b, c, rp.unwrap_or_else(V::zero)];
    carlson(r.map(Work::square), rp.map(Work::square), Some(roots), want)
}

fn check_rd_roots<T: CordicNumber>(rx: T, ry: T, rz: T, name: &'static str) -> Result<()> {
    let zero = T::zero();
    if rz == zero || (rx == zero && ry == zero) {
        return Err(Error::domain(name, "rx, ry not both zero, rz non-zero"));
    }
    Ok(())
}

/// [`carlson_rd`] of the squares of `rx`, `ry` and `rz`, for arguments too small to hold.
/// # Errors
/// Returns `DomainError` if `rx` and `ry` are both 0, or `rz` is.
#[must_use = "returns the integral which should be handled"]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn carlson_rd_roots<T: CordicNumber>(rx: T, ry: T, rz: T) -> Result<T> {
    check_rd_roots(rx, ry, rz, "carlson_rd_roots")?;
    Ok(tiered!(T, |v| squares(
        [v(rx), v(ry), v(rz)],
        None,
        [false, true]
    )[1]))
}

/// [`carlson_rf`] and [`carlson_rd`] of the squares of `rx`, `ry` and `rz`, from one duplication.
/// # Errors
/// Returns `DomainError` if `rx` and `ry` are both 0, or `rz` is.
#[must_use = "returns the integrals which should be handled"]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn carlson_rf_rd_roots<T: CordicNumber>(rx: T, ry: T, rz: T) -> Result<(T, T)> {
    check_rd_roots(rx, ry, rz, "carlson_rf_rd_roots")?;
    Ok(tiered!(T, |v| {
        let [f, d, _] = squares([v(rx), v(ry), v(rz)], None, [true, true]);
        (f, d)
    }))
}

/// [`carlson_rd`] and [`carlson_rj`] of the squares of the arguments, from one duplication.
/// # Errors
/// Returns `DomainError` if `rx` and `ry` are both 0, or `rz` or `rp` is.
#[must_use = "returns the integrals which should be handled"]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn carlson_rd_rj_roots<T: CordicNumber>(rx: T, ry: T, rz: T, rp: T) -> Result<(T, T)> {
    check_rd_roots(rx, ry, rz, "carlson_rd_rj_roots")?;
    if rp == T::zero() {
        return Err(Error::domain("carlson_rd_rj_roots", "rp non-zero"));
    }
    Ok(tiered!(T, |v| {
        let [_, d, j] = squares([v(rx), v(ry), v(rz)], Some(v(rp)), [false, true]);
        (d, j)
    }))
}

/// [`carlson_rf`], [`carlson_rd`] and [`carlson_rj`] of the squares of the arguments, from one
/// duplication.
/// # Errors
/// Returns `DomainError` if `rx` and `ry` are both 0, or `rz` or `rp` is.
#[must_use = "returns the integrals which should be handled"]
#[cfg_attr(
    all(feature = "verify-no-panic", not(debug_assertions)),
    no_panic::no_panic
)]
pub fn carlson_rf_rd_rj_roots<T: CordicNumber>(rx: T, ry: T, rz: T, rp: T) -> Result<(T, T, T)> {
    check_rd_roots(rx, ry, rz, "carlson_rf_rd_rj_roots")?;
    if rp == T::zero() {
        return Err(Error::domain("carlson_rf_rd_rj_roots", "rp non-zero"));
    }
    let all: [T; 3] = tiered!(T, |v| squares(
        [v(rx), v(ry), v(rz)],
        Some(v(rp)),
        [true, true]
    ));
    Ok(all.into())
}
