//! Complex arithmetic in the working type, and the functions built on it.

use crate::consts::{EULER_GAMMA, FRAC_2_SQRT_PI, INV_FACT};
use crate::eval::{
    HUGE, atan_ratio, beyond, by_width, exp, exp_wide, expm1, half_pi, ln, log1p, log1p_ratio,
    phi as phi_real, pi, sin_cos, sin_cos_octant,
};
use crate::real::Work;
use crate::stirling::{LN_GAMMA_1P, LN_SQRT_2PI, STIRLING};
use core::ops::{Add, Div, Mul, Neg, Sub};

#[derive(Clone, Copy, Debug)]
pub struct Cx<V> {
    pub re: V,
    pub im: V,
}

pub fn atan2<V: Work>(y: V, x: V) -> V {
    let angle = if x.is_zero() {
        if y.is_zero() { V::zero() } else { half_pi() }
    } else {
        let base = atan_ratio(y.abs(), x.abs());
        if x.is_neg() { pi::<V>() - base } else { base }
    };
    angle.with_sign(y.is_neg())
}

impl<V: Work> Cx<V> {
    pub const fn new(re: V, im: V) -> Self {
        Self { re, im }
    }

    pub fn real(re: V) -> Self {
        Self::new(re, V::zero())
    }

    pub fn zero() -> Self {
        Self::real(V::zero())
    }

    pub const fn parts(self) -> (V, V) {
        (self.re, self.im)
    }

    #[must_use]
    pub fn scale(self, k: V) -> Self {
        Self::new(self.re * k, self.im * k)
    }

    #[must_use]
    pub fn ldexp(self, k: i32) -> Self {
        Self::new(self.re.ldexp(k), self.im.ldexp(k))
    }

    #[must_use]
    pub fn div_int(self, k: u64) -> Self {
        Self::new(self.re.div_int(k), self.im.div_int(k))
    }

    #[must_use]
    pub fn times_i(self) -> Self {
        Self::new(-self.im, self.re)
    }

    #[must_use]
    pub fn conj(self) -> Self {
        Self::new(self.re, -self.im)
    }

    pub fn norm2(self) -> V {
        self.re.square() + self.im.square()
    }

    pub fn abs(self) -> V {
        self.norm2().sqrt()
    }

    pub fn is_zero(&self) -> bool {
        self.re.is_zero() && self.im.is_zero()
    }

    #[must_use]
    pub fn square(self) -> Self {
        Self::new(
            self.re.square() - self.im.square(),
            (self.re * self.im).ldexp(1),
        )
    }

    #[must_use]
    pub fn recip(self) -> Self {
        let n = self.norm2();
        Self::new(self.re / n, -self.im / n)
    }

    /// `e^self`, saturating or vanishing past [`HUGE`].
    #[must_use]
    pub fn exp(self) -> Self {
        let modulus = exp_wide(self.re);
        let (s, c) = sin_cos(self.im);
        Self::new(modulus * c, modulus * s)
    }

    /// `e^self − 1`, its real part `expm1(x)·cos y − 2sin²(y/2)` to hold near 0.
    #[must_use]
    pub fn exp_m1(self) -> Self {
        let (s, c) = sin_cos(self.im);
        if self.re > V::from_int(HUGE) {
            return Self::new(beyond::<V>() * c, beyond::<V>() * s);
        }
        let (half, _) = sin_cos(self.im.ldexp(-1));
        let em1 = if self.re < V::from_int(-HUGE) {
            -V::one()
        } else {
            expm1(self.re)
        };
        Self::new(em1 * c - half.square().ldexp(1), (em1 + V::one()) * s)
    }

    /// `e^(−self²)` from exact `y² − x²` and `2xy`, its phase right at any size.
    #[must_use]
    pub fn exp_neg_square(self) -> Self {
        let modulus = exp_wide(V::dot2(self.im, self.im, -self.re, self.re));
        let (oct, frac) = V::octant_of_product(self.re.ldexp(1), self.im);
        let (s, c) = sin_cos_octant(oct, frac);
        let s = if self.re.is_neg() == self.im.is_neg() {
            -s
        } else {
            s
        };
        Self::new(modulus * c, modulus * s)
    }

    /// The principal logarithm, `ln|z|` kept exact in form near `|z| = 1`.
    #[must_use]
    pub fn ln(self) -> Self {
        let one = V::one();
        let near = (self.re - one) * (self.re + one) + self.im.square();
        let modulus = if near.abs() <= one.ldexp(-2) {
            log1p(near).ldexp(-1)
        } else {
            ln(self.norm2()).ldexp(-1)
        };
        Self::new(modulus, atan2(self.im, self.re))
    }
}

impl<V: Work> Add for Cx<V> {
    type Output = Self;
    fn add(self, o: Self) -> Self {
        Self::new(self.re + o.re, self.im + o.im)
    }
}

impl<V: Work> Sub for Cx<V> {
    type Output = Self;
    fn sub(self, o: Self) -> Self {
        Self::new(self.re - o.re, self.im - o.im)
    }
}

impl<V: Work> Neg for Cx<V> {
    type Output = Self;
    fn neg(self) -> Self {
        Self::new(-self.re, -self.im)
    }
}

impl<V: Work> Mul for Cx<V> {
    type Output = Self;
    fn mul(self, o: Self) -> Self {
        Self::new(
            self.re * o.re - self.im * o.im,
            self.re * o.im + self.im * o.re,
        )
    }
}

impl<V: Work> Div for Cx<V> {
    type Output = Self;
    fn div(self, o: Self) -> Self {
        let n = o.norm2();
        let t = self * o.conj();
        Self::new(t.re / n, t.im / n)
    }
}

fn settled<V: Work>(t: Cx<V>, sum: Cx<V>) -> bool {
    V::negligible(&t.norm2().sqrt(), &sum.norm2().sqrt())
}

fn two_over_root_pi<V: Work>() -> V {
    V::constant(&FRAC_2_SQRT_PI)
}

/// Terms past which every series and fraction here stops.
const CAP: u32 = 600;

/// The radius, 4 or 5.5, within which a series cancelling like `e^|z|²` keeps
/// 9 bits past a 32-bit ulp with one limb, or 20 past a 128-bit one.
fn kummer_reach<V: Work>() -> V {
    V::from_int(i64::from(crate::eval::by_width::<V>([8, 11, 11]))).ldexp(-1)
}

/// The same for `e^(2|z|)`, 6 or 8. Twice it bounds `r + Re z`, the power of
/// `e` by which the terms for `E₁` and `γ(a, z)` outgrow their sum; past it
/// Lentz's fraction converges like `e^(−2√(2k(r + Re z)))`.
fn linear_reach<V: Work>() -> V {
    V::from_int(i64::from(crate::eval::by_width::<V>([6, 8, 8])))
}

/// `D(z) = z·Σ (−2z²)ⁿ/(2n + 1)!!`, where `w(z) = e^(−z²) + (2i/√π)·D(z)`.
pub fn dawson_series<V: Work>(z: Cx<V>) -> Cx<V> {
    let step = -(z.square().ldexp(1));
    let mut term = z;
    let mut sum = z;
    let past = V::from_int(2) * z.norm2();
    for n in 1..CAP {
        term = (term * step).div_int(u64::from(2 * n + 1));
        sum = sum + term;
        if V::from_int(i64::from(n)) > past && settled(term, sum) {
            break;
        }
    }
    sum
}

/// `w(z) ~ (i/√π)/z·Σ (2n − 1)!!/(2z²)ⁿ`, for `|z|` where the smallest term,
/// about `e^(−|z|²)`, is past the working precision.
fn asymptotic<V: Work>(z: Cx<V>) -> Cx<V> {
    let inv = (z.square().ldexp(1)).recip();
    let mut term = Cx::real(V::one());
    let mut sum = term;
    let mut last = term.norm2();
    for n in 1..CAP {
        let next = (term * inv).scale(V::from_int(i64::from(2 * n - 1)));
        let size = next.norm2();
        if size > last || settled(next, sum) {
            break;
        }
        sum = sum + next;
        term = next;
        last = size;
    }
    (sum / z).scale(two_over_root_pi::<V>().ldexp(-1)).times_i()
}

/// `w(z)` for `Im z ≥ 1`, `(i/√π)/(z − ½/(z − 1/(z − (3/2)/(z − …))))`.
fn continued<V: Work>(z: Cx<V>) -> Cx<V> {
    let g = lentz(z, |k| (Cx::real(V::from_int(-i64::from(k)).ldexp(-1)), z));
    g.recip().scale(two_over_root_pi::<V>().ldexp(-1)).times_i()
}

fn upper<V: Work>(z: Cx<V>) -> Cx<V> {
    let r = z.abs();
    if r <= kummer_reach() {
        return z.exp_neg_square() + dawson_series(z).scale(two_over_root_pi()).times_i();
    }
    let far = V::from_int(i64::from(crate::eval::by_width::<V>([8, 11, 13])));
    if r >= far {
        return asymptotic(z);
    }
    let one = V::one();
    if z.im >= one {
        return continued(z);
    }
    // Near the real axis, Taylor from z₀ = x + i by w' = −2zw + 2i/√π.
    let z0 = Cx::new(z.re, one);
    let w0 = if z0.abs() >= far {
        asymptotic(z0)
    } else {
        continued(z0)
    };
    let delta = Cx::new(V::zero(), z.im - one);
    let lead = (-(z0 * w0).ldexp(1) + Cx::real(two_over_root_pi()).times_i()) * delta;
    let (mut prev, mut cur) = (w0, lead);
    let mut sum = w0 + lead;
    let a = -(z0 * delta).ldexp(1);
    let b = -(delta.square()).ldexp(1);
    for k in 1..CAP {
        let next = (a * cur + b * prev).div_int(u64::from(k + 1));
        sum = sum + next;
        if settled(next, sum) && settled(cur, sum) {
            break;
        }
        (prev, cur) = (cur, next);
    }
    sum
}

/// `w(z) = e^(−z²)·erfc(−iz)`, real on the imaginary axis.
pub fn faddeeva<V: Work>(z: Cx<V>) -> Cx<V> {
    let w = if z.im.is_neg() {
        // w(z) = 2e^(−z²) − w(−z).
        z.exp_neg_square().ldexp(1) - upper(-z)
    } else {
        upper(z)
    };
    if z.re.is_zero() { Cx::real(w.re) } else { w }
}

/// `Re erf(x + iy)` for `|2xy| ≤ ½`, which `|erf z| ≈ e^(y²)/y` swamps below:
/// `(2/√π)·e^(y²)·∫₀ˣ e^(−s²)·cos 2ys ds` in `x`, by `gₖ₊₁ = (−2iy·gₖ − 2gₖ₋₁)/(k + 1)`
/// for `e^(−2iys − s²)`, `gₖ/i` held for odd `k`.
fn erf_re_near_axis<V: Work>(x: V, y: V) -> V {
    let two_y = y.ldexp(1);
    let (mut prev, mut cur) = (V::one(), -two_y);
    let x2 = x.square();
    let mut power = x;
    let mut sum = x;
    for k in 1..CAP {
        let turned = if k % 2 == 1 {
            two_y * cur
        } else {
            -(two_y * cur)
        };
        (prev, cur) = (cur, (turned - prev.ldexp(1)).div_int(u64::from(k + 1)));
        if k % 2 == 1 {
            power = power * x2;
            let term = cur.div_int(u64::from(k + 2)) * power;
            sum = sum + term;
            if V::negligible(&term, &sum) {
                break;
            }
        }
    }
    exp_wide(y.square()) * sum * two_over_root_pi()
}

/// `erf z`, odd and real on the real axis: `(2/√π)·e^(−z²)·Σ (2z²)ⁿz/(2n + 1)!!`
/// near 0, and past it `1 − e^(−z²)·w(iz)`.
pub fn erf<V: Work>(z: Cx<V>) -> Cx<V> {
    let flip = z.re.is_neg();
    let u = if flip { -z } else { z };
    let r = u.abs();
    let value = if r <= kummer_reach() {
        // −i·D(iu) = u·Σ (2u²)ⁿ/(2n + 1)!!.
        let s = -(dawson_series(u.times_i()).times_i());
        (u.exp_neg_square() * s).scale(two_over_root_pi())
    } else {
        Cx::real(V::one()) - u.exp_neg_square() * faddeeva(u.times_i())
    };
    let value = if u.im.abs() > V::one() && (u.re * u.im).abs().ldexp(2) <= V::one() {
        Cx::new(erf_re_near_axis(u.re, u.im), value.im)
    } else {
        value
    };
    let value = if flip { -value } else { value };
    match (z.re.is_zero(), z.im.is_zero()) {
        (true, _) => Cx::new(V::zero(), value.im),
        (_, true) => Cx::real(value.re),
        _ => value,
    }
}

/// `b₀ + a₁/(b₁ + a₂/(b₂ + …))` by Lentz's method, `(aₖ, bₖ)` from `term(k)`.
fn lentz<V: Work>(b0: Cx<V>, term: impl Fn(u32) -> (Cx<V>, Cx<V>)) -> Cx<V> {
    let tiny = Cx::real(V::one().ldexp(-1000));
    let mut value = if b0.is_zero() { tiny } else { b0 };
    let mut upper = value;
    let mut lower = Cx::zero();
    let limit = V::one().ldexp(20 - 2 * V::PREC.cast_signed());
    for k in 1..8 * CAP {
        let (num, den) = term(k);
        lower = den + num * lower;
        if lower.is_zero() {
            lower = tiny;
        }
        upper = den + num / upper;
        if upper.is_zero() {
            upper = tiny;
        }
        lower = lower.recip();
        let delta = upper * lower;
        value = value * delta;
        if (delta - Cx::real(V::one())).norm2() < limit {
            break;
        }
    }
    value
}

/// `Γ(a, z) ~ z^(a−1)·e^(−z)·Σ (a − 1)…(a − k)/z^k`, for `|z|` past 256.
fn gamma_asymptotic<V: Work>(a: V, z: Cx<V>) -> Cx<V> {
    let inv = z.recip();
    let mut term = Cx::real(V::one());
    let mut sum = term;
    let mut last = V::one();
    for k in 1..CAP {
        let next = (term * inv).scale(a - V::from_int(i64::from(k)));
        let size = next.norm2();
        if size > last || settled(next, sum) {
            break;
        }
        sum = sum + next;
        term = next;
        last = size;
    }
    power(z, a - V::one()) * (-z).exp() * sum
}

/// `E₁(z)`, `z ≠ 0`: `−γ − ln z − Σ (−z)ⁿ/(n·n!)` near 0 or the negative axis,
/// asymptotic far along it, else `e^(−z)/(z + 1 − 1/(z + 3 − 4/(z + 5 − …)))`.
pub fn e1<V: Work>(z: Cx<V>) -> Cx<V> {
    let r = z.abs();
    let series = r + z.re <= linear_reach::<V>().ldexp(1);
    if series && r > V::from_int(256) {
        let value = gamma_asymptotic(V::zero(), z);
        // On the axis, the value from above it: −Ei(−z) − iπ.
        return if z.im.is_zero() {
            Cx::new(value.re, -pi::<V>())
        } else {
            value
        };
    }
    if series {
        let mut term = Cx::real(V::one());
        let mut sum = Cx::zero();
        for n in 1..8 * CAP {
            term = (-(term * z)).div_int(u64::from(n));
            let part = term.div_int(u64::from(n));
            sum = sum + part;
            if V::from_int(i64::from(n)) > r && settled(part, sum) {
                break;
            }
        }
        return -(Cx::real(V::constant(&EULER_GAMMA)) + z.ln() + sum);
    }
    let one = Cx::real(V::one());
    let f = lentz(z + one, |k| {
        let i = V::from_int(i64::from(k));
        (Cx::real(-i.square()), z + one + Cx::real(i.ldexp(1)))
    });
    (-z).exp() / f
}

/// `ln Γ(x)` for real `x > 0`: shifted to 40 or past, then Stirling's series.
fn ln_gamma<V: Work>(x: V) -> V {
    let mut shift = V::one();
    let mut y = x;
    for _ in 0..40 {
        if y >= V::from_int(40) {
            break;
        }
        shift = shift * y;
        y = y + V::one();
    }
    let inv = V::one() / y;
    let inv2 = inv.square();
    let mut power = inv;
    let mut sum = V::zero();
    for c in STIRLING
        .iter()
        .take(crate::eval::by_width::<V>([12, 24, 40]) as usize)
    {
        sum = sum + V::constant(c) * power;
        power = power * inv2;
    }
    (y - V::one().ldexp(-1)) * ln(y) - y + V::constant(&LN_SQRT_2PI) + sum - ln(shift)
}

pub fn gamma<V: Work>(x: V) -> V {
    let lg = ln_gamma(x);
    if lg > V::from_int(HUGE) {
        beyond()
    } else {
        exp(lg)
    }
}

/// `(Γ(1 + a) − 1)/a`, below 2^−8 by the Taylor series of `ln Γ(1 + a)`.
fn gamma_1p_ratio<V: Work>(a: V) -> V {
    if a >= V::one().ldexp(-8) {
        return (gamma(a + V::one()) - V::one()) / a;
    }
    let mut sum = V::zero();
    for c in LN_GAMMA_1P.iter().rev() {
        sum = V::constant(c) + a * sum;
    }
    expm1(a * sum) / a
}

/// `z^a` for real `a`; left of the imaginary axis `e^(±iπa)`, reduced exactly,
/// turned by `∓a·atan|y/x|`, so parts near the negative axis keep precision.
pub fn power<V: Work>(z: Cx<V>, a: V) -> Cx<V> {
    let log = z.ln();
    if !z.re.is_neg() {
        return log.scale(a).exp();
    }
    let whole = a.round_int();
    let rest = a - V::from_int(whole);
    let (sin, cos) = if rest.abs() == V::one().ldexp(-1) {
        (rest.ldexp(1), V::zero())
    } else {
        sin_cos(rest * pi::<V>())
    };
    let (st, ct) = sin_cos(a * atan_ratio(z.im.abs(), z.re.abs()));
    let turn = Cx::new(cos, sin) * Cx::new(ct, -st);
    let unit = if whole % 2 == 0 { turn } else { -turn };
    let unit = if z.im.is_neg() { unit.conj() } else { unit };
    unit.scale(exp_wide(a * log.re))
}

/// `Γ(a, z)` for `a` in `(0, 128]`: `Γ(a) − z^a·Σ (−z)ⁿ/(n!(a + n))` near the
/// negative axis or for `a ≤ ½`, asymptotic far along it, `Γ(a) − z^a·e^(−z)·Σ zⁿ/(a)ₙ₊₁`
/// within `a + 1` and more of 0, else Legendre's continued fraction.
pub fn gamma_upper<V: Work>(a: V, z: Cx<V>) -> Cx<V> {
    if z.is_zero() {
        return Cx::real(gamma(a));
    }
    let r = z.abs();
    let one = V::one();
    let small = a <= one.ldexp(-1);
    let left = z.re.is_neg() && r > a + one;
    if (small || left) && r + z.re <= linear_reach::<V>().ldexp(1) {
        if r > V::from_int(256) {
            let value = gamma_asymptotic(a, z);
            // On the axis Γ(a, z) = Γ(a) − e^(iπa)·G, G real: the series is the latter.
            let whole = V::from_int(a.round_int()) == a;
            return if z.im.is_zero() && !whole {
                Cx::new(value.re + gamma(a), value.im)
            } else {
                value
            };
        }
        let mut term = Cx::real(one);
        let mut sum = Cx::zero();
        for n in 1..8 * CAP {
            term = (-(term * z)).div_int(u64::from(n));
            let part = term.scale(one / (a + V::from_int(i64::from(n))));
            sum = sum + part;
            if V::from_int(i64::from(n)) > r && settled(part, sum) {
                break;
            }
        }
        if small {
            // Γ(a) − z^a/a as (Γ(1 + a) − 1)/a − (z^a − 1)/a.
            let head = gamma_1p_ratio(a);
            let ratio = z.ln().scale(a).exp_m1().scale(one / a);
            return Cx::real(head) - ratio - power(z, a) * sum;
        }
        return Cx::real(gamma(a)) - power(z, a) * (sum + Cx::real(one / a));
    }
    let reach = (a + V::one()).max(linear_reach());
    let lead = power(z, a) * (-z).exp();
    if r <= reach {
        let mut term = Cx::real(V::one() / a);
        let mut sum = term;
        let mut denom = a;
        for _ in 1..8 * CAP {
            denom = denom + V::one();
            term = (term * z).scale(V::one() / denom);
            sum = sum + term;
            if denom > r && settled(term, sum) {
                break;
            }
        }
        return Cx::real(gamma(a)) - lead * sum;
    }
    let f = lentz(z + Cx::real(V::one() - a), |k| {
        let i = V::from_int(i64::from(k));
        (
            Cx::real(-(i * (i - a))),
            z + Cx::real(i.ldexp(1) + V::one() - a),
        )
    });
    lead / f
}

/// `φₖ(z) = (e^z − Σ_{n<k} zⁿ/n!)/z^k` for `k` of 1 or 2: the series near 0.
pub fn phi<V: Work>(k: u32, z: Cx<V>) -> Cx<V> {
    if z.abs() <= V::one().ldexp(-1) {
        let n = by_width::<V>([18, 30, 39]);
        let mut p = Cx::zero();
        for j in (0..=n).rev() {
            let coef = INV_FACT
                .get((k + j) as usize)
                .map_or_else(V::zero, V::constant);
            p = Cx::real(coef) + z * p;
        }
        return p;
    }
    let mut num = z.exp_m1();
    let mut den = z;
    if k >= 2 {
        num = num - z;
        den = den * z;
    }
    num / den
}

const TERMS: usize = 320;

/// `Σ cₙ·Jₙ(w)`, `Jₙ(w) = ∫₀¹ sⁿ·e^(ws) ds`, `Re w ≤ 0`: forward to `|w|` by `Jₙ = (e^w − nJₙ₋₁)/w`,
/// past it backward from `J_M = e^w·Σ (−w)^k·M!/(M + k + 1)!`, `M ≥ 2|w|`, each way damping.
fn moment_sum<V: Work>(w: Cx<V>, coefs: &[V]) -> Cx<V> {
    let n = coefs.len().saturating_sub(1);
    let reach = usize::try_from(w.abs().floor_int().max(0)).unwrap_or(usize::MAX);
    let ew = w.exp();
    let coef = |k: usize| coefs.get(k).copied().unwrap_or_else(V::zero);
    let mut cur = phi(1, w);
    let mut sum = cur.scale(coef(0));
    let forward = reach.min(n);
    if forward > 0 {
        let inv = w.recip();
        for k in 1..=forward {
            cur = (ew - cur.scale(V::from_int(i64::try_from(k).unwrap_or(0)))) * inv;
            sum = sum + cur.scale(coef(k));
        }
    }
    if forward == n {
        return sum;
    }
    let top = n.max(reach.saturating_mul(2).saturating_add(1));
    let top_v = V::from_int(i64::try_from(top).unwrap_or(0));
    let mut term = Cx::real(V::one() / (top_v + V::one()));
    let mut series = term;
    for k in 2..TERMS {
        let kv = V::from_int(i64::try_from(k).unwrap_or(0));
        term = -(term * w).scale(V::one() / (top_v + kv));
        series = series + term;
        if settled(term, series) {
            break;
        }
    }
    let mut down = ew * series;
    for k in (forward + 1..=top).rev() {
        if k <= n {
            sum = sum + down.scale(coef(k));
        }
        down = (ew - w * down).div_int(u64::try_from(k).unwrap_or(1));
    }
    sum
}

/// `W(u) = (1 + δu)^(−a)` and `V(u) = (1 − W(u))/(aδ)`, `a = p + 1`, `δ = ε/(1 − ε)`.
struct Weight<V> {
    power: V,
    delta: V,
}

impl<V: Work> Weight<V> {
    fn log_at(&self, u: V) -> V {
        -(self.power * log1p(self.delta * u))
    }

    fn complement_at(&self, u: V) -> V {
        u * log1p_ratio(self.delta * u) * phi_real(1, self.log_at(u))
    }
}

/// `∫₀¹ e^(z(1 − u))·W(u) du`, or of `V` with `second`, from each end toward the least of the
/// log-convex `|e^(−zu)·W(u)|`: each step's `e^(ct)` whole by [`moment_sum`], `W`'s curve as series.
fn weighted_passes<V: Work>(z: Cx<V>, wt: &Weight<V>, second: bool) -> Cx<V> {
    let one = V::one();
    let (power, delta) = (wt.power, wt.delta);
    let meet = if z.re.is_neg() {
        (-(power / z.re) - one / delta).max(V::zero()).min(one)
    } else {
        one
    };
    let tiny = V::one().ldexp(-(V::PREC.cast_signed() + 8));
    let most = wt.complement_at(one);
    let curve = one.ldexp(-1) / power.sqrt();
    let (mut ahead, mut behind) = (Cx::zero(), Cx::zero());
    for forward in [true, false] {
        let pass = if forward { &mut ahead } else { &mut behind };
        let sign = if forward { one } else { -one };
        let mut u = if forward { V::zero() } else { one };
        for _ in 0..1024 {
            let rest = if forward { meet - u } else { u - meet };
            if rest <= V::zero() {
                break;
            }
            let log_w = wt.log_at(u);
            // Forward the integrand is e^(−zu)·W, e^z applied at the end.
            let from = if forward { -u } else { one - u };
            let lead = Cx::new(z.re * from + log_w, z.im * from).exp();
            let base = Cx::new(z.re * from, z.im * from).exp();
            // Its size falls toward the meeting point, so this bounds what is left.
            let left = if second {
                base.abs() * most
            } else {
                lead.abs()
            } * rest;
            if left < tiny.ldexp(-8) * pass.abs() {
                break;
            }
            let ratio = delta / (one + delta * u);
            let h = rest.min(curve / ratio).min(one.ldexp(-1) / ratio);
            let (rh, mut coefs) = (ratio * h, [V::zero(); TERMS]);
            let used = if second {
                // 1 − (1 ± rt)^(∓a) over aδ: rt/(δuⱼ + 1)·…, ratio ∓(a + k)rt/(k + 1).
                let mut c = h * h / (one + delta * u);
                let mut used = 1;
                for (k, slot) in coefs.iter_mut().enumerate().skip(1) {
                    *slot = c;
                    used = k;
                    if c.abs() < tiny * h {
                        break;
                    }
                    let kv = V::from_int(i64::try_from(k).unwrap_or(0));
                    c = -(sign * c * (power + kv) * rh).div_int(u64::try_from(k + 1).unwrap_or(1));
                }
                used
            } else {
                // e^(aρ(±rt)): b₀ = 1, b₁ = 0, bₖ₊₁ = (a·r²·bₖ₋₁ ∓ r·k·bₖ)/(k + 1), times h^(k+1).
                let curvature = power * rh.square();
                let (mut prev, mut cur) = (h, V::zero());
                let mut used = 1;
                if let Some(slot) = coefs.first_mut() {
                    *slot = h;
                }
                for (k, slot) in coefs.iter_mut().enumerate().skip(1) {
                    *slot = cur;
                    used = k;
                    if cur.abs() < tiny * h && prev.abs() < tiny * h {
                        break;
                    }
                    let kv = V::from_int(i64::try_from(k).unwrap_or(0));
                    let next = (curvature * prev - sign * rh * kv * cur)
                        .div_int(u64::try_from(k + 1).unwrap_or(1));
                    (prev, cur) = (cur, next);
                }
                used
            };
            let series = coefs.get(..=used).unwrap_or(&coefs);
            if second {
                let rate = (-z).scale(sign).scale(h);
                let flat = phi(1, rate).scale(h * wt.complement_at(u));
                let bend = moment_sum(rate, series).scale(sign * exp_wide(log_w));
                *pass = *pass + base * (flat + bend);
            } else {
                let rate = (-(z + Cx::real(power * ratio))).scale(sign * h);
                *pass = *pass + lead * moment_sum(rate, series);
            }
            u = if forward { u + h } else { u - h };
        }
    }
    z.exp() * ahead + behind
}

fn parts_sum<V: Work>(z: Cx<V>, a: V, q: V, first: Cx<V>) -> Cx<V> {
    let step = z.recip().scale(q);
    let mut term = first;
    let mut sum = term;
    for k in 0..TERMS {
        term = -(term * step).scale(a + V::from_int(i64::try_from(k).unwrap_or(0)));
        sum = sum + term;
        if settled(term, sum) {
            break;
        }
    }
    sum
}

/// `(K₁, K₂)` of `Complex::exprel_weighted`, `ε > 0`: `e^z/(1 − ε)·∫₀¹ e^(−zu)·{W, V}(u) du`, by
/// parts past `|z| = 4(a + PREC)δ` and 1, else stepped, `K₂ = (φ₁ − (1 − ε)K₁)/(εa)` for `εa ≥ ¼`.
pub fn exprel_weighted<V: Work>(z: Cx<V>, eps: V, p: V) -> (Cx<V>, Cx<V>) {
    let one = V::one();
    let a = p + one;
    let keep = one - eps;
    let delta = eps / keep;
    let log_keep = log1p(-eps);
    let ez = z.exp();
    let far = V::from_int(4) * (a + V::from_int(i64::from(V::PREC))) * delta;
    if z.abs() >= far.max(one) {
        let inv2 = z.square().recip();
        let k1 = (ez * parts_sum(z, a, delta, z.recip())).scale(one / keep)
            - parts_sum(z, a, eps, z.recip()).scale(exp_wide(p * log_keep));
        // V(1)/(1 − ε) = (1 − (1 − ε)^a)/(aε) = (Λ/ε)·φ₁(−aΛ), Λ = −log1p(−ε).
        let tail = log1p_ratio(-eps) * phi_real(1, a * log_keep);
        let k2 = -(z.recip().scale(tail))
            + (ez * parts_sum(z, a + one, delta, inv2)).scale(one / keep)
            - parts_sum(z, a + one, eps, inv2).scale(exp_wide(a * log_keep));
        return (k1, k2);
    }
    let wt = Weight { power: a, delta };
    let k1 = weighted_passes(z, &wt, false).scale(one / keep);
    let k2 = if eps * a >= one.ldexp(-2) {
        (phi(1, z) - k1.scale(keep)).scale(one / (eps * a))
    } else {
        weighted_passes(z, &wt, true).scale(one / keep)
    };
    (k1, k2)
}
