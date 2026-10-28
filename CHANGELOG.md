# Changelog

A release that changes any output bit of an existing function names that
function under the release, even for a one-ulp change. `tests/golden.txt`
pins a hash of every function's outputs on six types, and the tests fail when
one moves; `GOLDEN_UPDATE=1 cargo test golden` re-pins it and stamps each
changed row with the crate version, whose section here must then name the
function in backticks.

Releases before 3.1.1 are described on the
[releases page](https://github.com/GeEom/fixed_analytics/releases).

## 4.0.0

Every function now evaluates in binary floating point at twice its type's
width or more and rounds once, so each is within one ulp of the exact result
on every type, and nearest in every case the accuracy tool has met. Odd and
even functions are exactly so, monotone functions stay monotone, and sin,
cos and tan reduce their argument exactly at any magnitude.

Output changes, all of them corrections: `sin`, `cos`, `tan`, `sin_cos`,
`asin`, `acos`, `atan`, `atan2`, `sinh`, `cosh`, `tanh`, `sinh_cosh`, `coth`,
`asinh`, `acosh`, `atanh`, `acoth`, `exp`, `pow2`, `ln`, `log2`, `log10` and
`pow`. `sqrt` is unchanged. `tan` saturates only where the tangent is out of
range, and an odd function's result past `MIN` saturates to `MIN` rather than
`-MAX`. `acos` reports domain errors under its own name.

Breaking: `CordicNumber` is sealed, as its implementations now carry the
working arithmetic.

New: `mul_add`, `mul_div`, `div_or_saturate` and `dot`, which round an exact
result once, `sqrt_or_zero`, and `ERROR_ULPS`, the bound every function meets.

New: `enclose`, the interval of one ulp each way that holds a result's exact
value, and outward-rounded ranges over an interval, `sin_cos_range`,
`exp_range`, `atan_range` and `tanh_range`; `hypot`, `norm` and `normalize`,
which form their squares without underflow; `scale_exponent` and `scale_by`
for exact power-of-two scaling; `sqrt_product`, `rsqrt`, and `asin_from_gap`,
which takes the gap to 1 in place of a value that would round.

New, each to full relative precision near 0: `expm1`; `exprel` and
`exprel2`, `(e^x − 1)/x` and `(e^x − 1 − x)/x²`; `exp_moments`, the
integrals `∫₀¹ s^(k−1)·e^(−xs) ds` for k up to 3 with `e^(−x)`; `log1p`,
`log1p_ratio`, `log1p_residual` and `log1p_residual_ratio`; `stumpff`, the
Stumpff functions `c₀` to `c₃` of either sign; and `tan_ratio`, `atan_ratio`,
`atanh_ratio` and `tanh_ratio`.

New: elliptic integrals and functions. `Modulus` holds `m` or its
complement `k'` exact. `ellip_k`, `ellip_e` and `ellip_pi_complete` are the
complete integrals, `ellip_pi` the incomplete third kind from `(sn, cn, dn)`,
`jacobi_sn_cn_dn` and `jacobi_am` the Jacobi functions of any argument, and
`jacobi_inverse` their inverse. `carlson_rf`, `carlson_rj`, `carlson_rc` and
`carlson_rd` are Carlson's integrals, with `carlson_rf_roots` and
`carlson_rj_roots` taking square roots for arguments too small to hold.

New: angles in turns. `sin_cos_turns` takes a fraction of a turn, whose
reduction is its own bits with no rounded π, and `atan2_turns` returns one.
`div_by_tau` is `x/τ` rounded once, and `wrap_radians` and `wrap_principal`
reduce to `[0, 2π)` and `(−π, π]` by exact multiples of 2π.

New: `Complex`, a complex number over any of the types, whose arithmetic
rounds each part once and whose division cannot overflow for a non-zero
divisor, with `exp`, `exprel` and `exprel2`. The Faddeeva function
`faddeeva`, the complex `erf`, the exponential integral `e1` and the upper
incomplete gamma function `gamma_upper`, for `a` in `(0, 128]`, work on any
complex argument, on the principal branch. `erf`, `erfc`, `fresnel` and
`gamma` take a real one.

`exp_moments` returns `[E₁, …, E₅, e^(−x)]`, the moments up to `k = 5` all
from one exponential.

New: `Complex::exprel_weighted`, `exprel` and `exprel2` under the weight
`((1 − εt)/(1 − εs))^p/(1 − εs)` for complex `z`, `ε` in `[0, 1)` and `p ≥ 0`,
which keep their precision as `ε → 0`; the complex logarithm `Complex::ln`
and real power `Complex::pow`, exact in phase on the axes and diagonals;
`Complex::saturating_add`, `saturating_sub`, `saturating_mul`,
`saturating_div` and `saturating_neg`, the operators by name; `mul_by_tau`,
turns to radians rounded once, and `sin_cos_range_turns`; `root_in`, a
bracketed root finder by the ITP method with a fixed budget; and
`kepler_universal`, the universal anomaly of a two-body orbit.

Speed against 3.1.1, per call: I16F16's elementary functions take 3.4 to 9
times as long and I32F32's 1.6 to 6 times, while I64F64's run 0.65 to 2.4
times as long, `ln` 1.5 and `atan2` 1.2 times faster, as `sin_cos`, `exp`,
`ln` and `atan2` work in 126-bit fixed point on 128-bit types with at most 96
fractional bits; `sqrt` is unchanged. A function works in a float as wide as
its type when its result leaves 24 bits past the ulp, Carlson's duplication
stops once its closing series is exact, and the Jacobi functions take
Bulirsch's descent: on I64F64 `carlson_rj_roots` takes 6.7 µs and
`jacobi_sn_cn_dn` 2.5 µs.

## 3.1.1

No output changes. The package leaves out the toolchain and CI files.
