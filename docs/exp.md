# How `exp` is computed

This page records the derivations behind the comments in
`src/ops/exponential.rs`.

## Argument reduction

`exp(x) = 2^k · exp(r)` with `r ∈ [0, ln 2)`, so `exp(r) ∈ [1, 2)` carries a
full significand into the final shift.

`k = ⌊x / ln 2⌋` is estimated to within one and then corrected once in either
direction. `to_i32` floors the quotient. When the integer bits do not
outnumber the fractional bits, the quotient comes from a multiply by
`1/ln 2 − 1` (as I1F63, since `1/ln 2` itself exceeds 1) plus `x`, which is
within `|x|·1.5·2^-frac + 2^-frac < 1` of the true value. Other types divide.

`2^k · exp(r)` exceeds `MAX` once `k ≥ int_bits − 1` and rounds to zero once
`k ≤ −(frac_bits + 2)`. The early exits allow one of slack, and repeat the test
exactly after the correction.

## Polynomial

`exp(r) = 1 + r·(1 + r/2·(1 + r/3·(1 + … r/n)))`, evaluated innermost first
with integer division for each `r/K`. The omitted term `(ln 2)^(n+1)/(n+1)!`
is:

| Degree | Omitted term | Below half an ulp up to |
|--------|--------------|-------------------------|
| 7 | 1.3e-6 | 18 fractional bits |
| 12 | 1.4e-12 | 38 fractional bits |
| 19 | 2.7e-22 | 70 fractional bits |

The code uses degree 7 below 24 fractional bits, 12 from 24 and 19 from 40.

## Final scaling

A positive `k` shifts left after checking `exp(r)` against `MAX >> k`, since
the shift would wrap. A non-positive `k` shifts right with rounding to
nearest: truncation would bias results that are only a few ulps low by up to
a full ulp.
