# How `sin_cos` reduces its argument

This page records the derivation behind the comments in
`src/ops/circular.rs`.

An angle outside `[-π, π]` is reduced by `n = round(angle / 2π)` to
`angle − n·2π`, at any magnitude and without iteration limits.

When the integer bits do not outnumber the fractional bits, the quotient
comes from a multiply by `1/(2π)` as I1F63. That is within
`|angle|·1.5·2^-frac + 2^-frac < 1` of the true quotient, so `n` is off by at
most one, and the clamp to `[-π, π]` that follows absorbs it bit for bit.
Other types divide by 2π.

`n·2π` is exact and the difference is small, so wrapping arithmetic recovers
it even when `n·2π` exceeds `MAX`, for angles within π of `MAX`.

The reduced angle then folds to `[-π/2, π/2]`, negating both results. Past
`π/4` the polynomials take `π/2 − |x|` and swap roles, sin from the cosine
polynomial and cos from the sine one, which avoids cancellation near π/2.
