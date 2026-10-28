Hyperbolic functions and their inverses.

# Saturation

[`sinh`](self::sinh) and [`cosh`](self::cosh) grow as `e^|x|/2` and saturate
for extreme inputs: `sinh` to `T::MAX` or `T::MIN` past `±asinh(T::MAX)`, and
`cosh`, which is always positive, to `T::MAX` past `acosh(T::MAX)`.
[`sinh_cosh`](self::sinh_cosh) saturates both together, so
`cosh²(x) − sinh²(x) = 1` does not hold for saturated outputs.

| Type | Saturates for |
|------|---------------|
| I16F16 | \|x\| > ~11.1 |
| I32F32 | \|x\| > ~22.2 |

Below that, both come from `u = expm1(|x|)`, as `sinh = (u + u/(1 + u))/2`
and `cosh = 1 + u²/(2(1 + u))`, which keep their relative precision at any
size.
