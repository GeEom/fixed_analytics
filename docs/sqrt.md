# Correctly rounded square root for 128-bit types

`CordicNumber::sqrt_round` returns the raw bits `⌊√N + ½⌋` with `N = X·2^f`,
where `X` is the raw value and `f` the number of fractional bits.

## Types up to 64 bits

When `X << (f + 2)` fits in a wider machine integer, the root comes from one
integer square root: `⌊√N + ½⌋ = (⌊2√N⌋ + 1) >> 1`, and `⌊2√N⌋ = isqrt(4N)`.
The root has at most `(total_bits + f) / 2 + 1` bits, so it fits the raw type.

## 128-bit types

`N` can need up to 254 bits, so `sqrt_round_u128` in `src/traits.rs` works in
two steps.

1. Shift `x` left by the largest `s ≤ leading_zeros(x)` with `s ≡ f (mod 2)`,
   giving `q = isqrt(x·2^s)` with 64 significant bits.
2. If `s ≥ f`, `q >> ((s − f)/2)` is exactly `⌊√N⌋`, and the rounding bit
   comes from one more bit of the root, or from the remainder `N − q²` when
   `s = f`.
3. Otherwise `seed = (q + 1) << m`, with `m = (f − s)/2`, lies above `√N` by
   less than `2^(m+1)`. One Newton step from above,
   `(seed + ⌊N / seed⌋) / 2`, lands on `⌊√N⌋` or `⌊√N⌋ + 1`, since its error
   is below `2^(m − 64)` and `m ≤ 63`.
4. A 256-bit product `newton²` against `N` settles which, and the remainder
   `rem = N − root²`, in `[0, 2·root]`, rounds up exactly when `rem > root`.

The division `⌊x·2^f / seed⌋` is the fixed-point division of the raw values.
It is called only with `seed ≥ √N`, so the quotient is below `2^127` and
cannot overflow; the saturating fallback only keeps the function total.

## The working arithmetic

`Real::sqrt` in `src/real.rs` rounds the same way. With the exponent made
even, the mantissa is an integer `X` of 256 bits, or 384 for three limbs,
whose top 128 bits `x` are at least `2^126`.

1. `s = ⌊√x⌋` and `x − s²` come from a tabled `1/√x`, two Newton steps on
   it, the product `x·(1/√x)` corrected once, and a step or two to the exact
   root.
2. With `a1` and `a0` the next two limbs of `X`, a Karatsuba step,
   `q = ⌊((x − s²)·2^64 + a1)/2s⌋` with remainder `u`, gives `s·2^64 + q`,
   the root of `x·2^128 + a1·2^64 + a0` or one more. The remainder
   `u·2^64 + a0 − q²` says which, without a full square.
3. Two limbs are done: the root rounds up when its remainder exceeds it.
   Three limbs take one more such step, `s` now the 128-bit root and `q` a
   digit of long division.

`q` reaches `2^64` only when `x − s² = 2s`, and capping it gives the root
`(s + 1)·2^64 − 1`.
