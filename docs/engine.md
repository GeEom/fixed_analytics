# How the functions evaluate

This page records the design behind `src/real.rs` and `src/eval.rs`.

## The working arithmetic

Every function converts its arguments exactly into a binary floating-point
value, evaluates there, and rounds once to the fixed-point type. The mantissa
has 64 bits for 8-, 16- and 32-bit types, 128 for 64-bit ones and 192 for
128-bit ones, so each type works with at least 64 bits to spare above its
own. Each operation rounds to nearest even. Division by more than one limb
uses a Newton reciprocal, good to a few units in the last place.

Sign and magnitude are kept apart, so `f(−x)` repeats the work of `f(x)` on
the same magnitude and odd and even functions are exactly so. The final
rounding is monotone, and the working error, below 2^-20 ulp of the type, is
far smaller than any gap between neighbouring results except where the
function is flat, so monotone functions give monotone results. The accuracy
tool checks both on every input of the 8- and 16-bit types.

## Reductions

`sin` and `cos` reduce exactly. The magnitude of the argument times
`1/(2π)` to 384 bits gives the fraction of a turn with at least 257 bits
past the point, enough for any 128-bit integer argument. Its top three bits
pick the octant; odd octants count back from their top, so the remainder
`w` lies in `[0, 1]` and `(π/4)·w` in `[0, π/4]`.

The remaining work is short series around tabled points, the tables holding
384 bits:

| Function | Table | Remainder |
|----------|-------|-----------|
| `sin`, `cos` | `sin(jπ/128)`, `cos(jπ/128)` | `|d| ≤ π/256` |
| `exp` | `2^(j/32)` | `|r| ≤ ln 2/64` |
| `ln` | `ln(1 + j/16)` on `[√½, √2)` | `|s| ≤ 0.023` in `2·atanh(s)` |
| `atan` | `atan(j/16)` on `[0, 1]` | `|t| ≤ 1/32` |

`expm1` and `log1p` sum their series directly near 0, so they keep relative
precision, and the hyperbolic functions and their inverses are written in
them so that none cancels: `sinh = (u + u/(1 + u))/2` with `u = expm1(x)`,
`tanh = u/(u + 2)` with `u = expm1(2x)`, `asinh(z) = log1p(z + z²/(1 +
√(1 + z²)))`, `acosh(1 + t) = log1p(t + √(t(t + 2)))`, and `atanh(a) =
log1p(2a/(1 − a))/2`. `asin` and `acos` take `√((1 − a)(1 + a))`, whose
factors are exact.

## Series

Each series is a polynomial in the square of the remainder, its coefficients
`1/n!` or `1/(2k + 1)`, summed by Horner's rule in fixed point with 62, 126
or 190 fractional bits for the three widths. Its length is fixed per width,
long enough for the largest remainder with two bits to spare, so every
function does a bounded amount of work.
