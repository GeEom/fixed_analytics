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
