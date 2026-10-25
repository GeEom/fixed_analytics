Exponential and logarithmic functions.

# Saturation

[`exp`](self::exp), [`pow2`](self::pow2) and [`pow`](self::pow) saturate for
extreme inputs rather than returning an error. The result is `T::MAX` above
`ln(T::MAX)`, or `log2(T::MAX)` for `pow2`, and zero where the true value
rounds to nearest below half an ulp:

| Function | Result | I16F16 | I32F32 |
|----------|--------|--------|--------|
| `exp` | `T::MAX` | x > ~10.4 | x > ~21.5 |
| `exp` | zero | x < ~-11.8 | x < ~-22.9 |
| `pow2` | `T::MAX` | x > ~15 | x > ~31 |
| `pow2` | zero | x < ~-17 | x < ~-33 |

Saturation is silent and deterministic. To detect it, check the input range
before the call:

```
use fixed::types::I16F16;
use fixed_analytics::exp;

let x = I16F16::from_num(5.0);
let max_safe = I16F16::from_num(10.0);

if x < max_safe {
    let result = exp(x);
} else {
    // Handle potential saturation
}
```
