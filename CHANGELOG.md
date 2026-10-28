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

Slower: on the benchmarks I16F16 functions take 2.4 to 11 times as long as
in 3.1.1, and I64F64 ones 1.3 to 9 times; `sqrt` is as fast or faster.

## 3.1.1

No output changes. The package leaves out the toolchain and CI files.
