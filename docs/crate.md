Fixed-point mathematical functions: accurate, deterministic, and guaranteed not to panic.

`no_std` compatible with no floating-point operations.

# Quick Start

```rust
use fixed::types::I16F16;
use fixed_analytics::{sin, cos, sqrt, ln};

let angle = I16F16::from_num(0.5);
let (s, c) = (sin(angle), cos(angle));

let root = sqrt(I16F16::from_num(2.0)).unwrap();
assert!((root.to_num::<f32>() - 1.414).abs() < 0.001);

let log = ln(I16F16::E).unwrap();
assert!((log.to_num::<f32>() - 1.0).abs() < 0.01);
```

# Available Functions

**Total functions** return `T` directly, saturating on overflow.
**Fallible functions** return [`Result<T, Error>`] on domain violations.

| Category | Total | Fallible |
|--------------|-------|----------|
| Trigonometric | [`sin`], [`cos`], [`tan`], [`sin_cos`], [`atan`], [`atan2`] | [`asin`], [`acos`], [`asin_from_gap`] |
| Hyperbolic | [`sinh`], [`cosh`], [`tanh`], [`sinh_cosh`], [`asinh`] | [`acosh`], [`atanh`], [`acoth`], [`coth`] |
| Exponential | [`exp`], [`pow2`] | [`ln`], [`log2`], [`log10`], [`pow`] |
| Algebraic | [`sqrt_or_zero`], [`hypot`], [`norm`], [`scale_exponent`], [`scale_by`] | [`sqrt`], [`sqrt_product`], [`rsqrt`], [`normalize`] |
| Arithmetic | [`mul_add`], [`mul_div`], [`div_or_saturate`], [`dot`] | — |
| Enclosures | [`enclose`], [`sin_cos_range`], [`exp_range`], [`atan_range`], [`tanh_range`] | — |

Each function works in binary floating point at twice its type's width or
more, from its arguments converted exactly, and rounds once. Complete absence
of panic is verified at the linker level via the
[`no-panic`](https://github.com/dtolnay/no-panic) crate.

# Accuracy

Every function is within [`ERROR_ULPS`], one ulp, of the exact result on
every type, from `I8F0` to `I3F125`, and the arithmetic functions within half
an ulp; results out of range saturate to the nearer end. Odd and even
functions are exactly odd and even, except that an odd one saturating at
`MAX` meets `MIN` on the other side, and functions that are monotone stay
monotone.

The accuracy tool checks all of this against MPFR on 104 layouts: every
value of the 8- and 16-bit ones, and seeded samples of the wider ones.
Accuracy regressions are not permitted across releases, and a release that
changes any output bit says so in the changelog. The enclosure functions
turn the bound into intervals certain to hold the exact value.

# Features

- **`std`** (default): Enables `std::error::Error` impl on [`Error`]

See the [`kernel`] module for algorithm details.
