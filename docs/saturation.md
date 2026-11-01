# Saturation and domain errors

A total function returns `T` for every input. Past the type's range it
returns the nearer end, `T::MAX` or `T::MIN`. A result smaller than half an
ulp rounds to zero. Nothing reports either case, and both are deterministic.

A fallible function returns `Err(Error::DomainError { function, expected })`
only where no finite real result exists. Inside its domain it saturates like
a total function; `pow`, `coth` and `gamma` can.

`T::MIN` lies one ulp further from zero than `T::MAX`, so an odd function
stays odd only until it saturates. Where `sinh(x)` is `T::MAX`, `sinh(-x)` is
`T::MIN`, not `-T::MAX`.

## Thresholds

The thresholds follow from the range alone, because every result is rounded
from the exact value. Take a type with `i` integer bits, sign included, and
`f` fractional bits. Its `MAX` is about `2^(i−1)` and its ulp `2^−f`. The
values below were measured on the crate and rounded as shown. Each function
saturates from the first representable input past its value.

| Function | Saturates | I16F16 | I32F32 | I64F64 |
|----------|-----------|--------|--------|--------|
| `exp` | to `MAX` from `x = ln MAX` | 10.397 | 21.488 | 43.668 |
| `exp` | to 0 from `x = ln(ulp/2)` down | −11.784 | −22.874 | −45.055 |
| `expm1` | to `MAX` from `x = ln(MAX + 1)` | 10.397 | 21.488 | 43.668 |
| `pow2` | to `MAX` from `x = log₂ MAX` | 15 | 31 | 63 |
| `pow2` | to 0 from `x = log₂(ulp/2)` down | −17 | −33 | −65 |
| `exprel` | to `MAX` where `e^x/x` passes `MAX` | 12.959 | 24.694 | 47.530 |
| `exprel2` | to `MAX` where `e^x/x²` passes `MAX` | 15.934 | 28.164 | 51.554 |
| `exp_moments`, `e^(−x)` | to `MAX` from `x = −ln MAX` down | −10.397 | −21.488 | −43.668 |
| `exp_moments`, `E₁` | to `MAX` where `e^\|x\|/\|x\|` passes `MAX`, `x < 0` | −12.959 | −24.694 | −47.530 |
| `exp_moments`, `E₅` | to `MAX`, a little further out | −13.260 | −24.855 | −47.614 |
| `sinh`, `cosh`, `sinh_cosh` | to `MAX`, or `MIN` for `sinh` of `x < 0`, from `\|x\| = ln(2·MAX)` | 11.090 | 22.181 | 44.361 |
| `tan` | to `MAX` or `MIN` within `1/MAX` of a pole `π/2 + kπ` | 3.1e-5 | 4.7e-10 | 1.1e-19 |
| `tan_ratio`, `tan x / x` | to `MAX` within `1/(MAX·\|x\|)` of a pole, here `π/2` | 1.9e-5 | 3.0e-10 | 6.9e-20 |
| `coth` | to `MAX` or `MIN` for `0 < \|x\| ≤ 1/MAX` | 3.1e-5 | 4.7e-10 | 1.1e-19 |
| `gamma` | to `MAX` where `Γ(x)` passes `MAX` | 8.903 | 13.589 | 21.440 |
| `gamma` | to `MAX` at the smallest input, where `Γ(x) ≈ 1/x` | 1.5e-5 | 2.3e-10 | 5.4e-20 |
| `stumpff`, `c₀` | to `MAX` where `cosh √−z` passes `MAX`, `z < 0` | −123.00 | −491.98 | −1967.94 |
| `stumpff`, `c₁` | to `MAX` where `sinh √−z / √−z` passes `MAX` | −187.92 | −645.98 | −2326.86 |
| `stumpff`, `c₂` | to `MAX`, further out | −279.69 | −835.73 | −2732.62 |
| `stumpff`, `c₃` | to `MAX`, further out | −403.66 | −1065.18 | −3187.97 |
| `mul_by_tau` | to `MAX` or `MIN` from `\|x\| = MAX/2π` | 5215.19 | 3.418e8 | 1.468e18 |
| `carlson_rj`, `carlson_rd` | to `MAX` with every argument at or below `a = MAX^(−2/3)` | 9.8e-4 | 6.0e-7 | 2.3e-13 |
| `carlson_rf_roots` | to `MAX` with every root at or below `1/MAX` | 3.1e-5 | 4.7e-10 | 1.1e-19 |

`E₂` to `E₄` fall between `E₁` and `E₅`.

Carlson's integrals are homogeneous. Scaling every argument by `λ` scales
`R_J` and `R_D` by `λ^(−3/2)`, but `R_F` and `R_C` only by `λ^(−1/2)`, so
small arguments saturate the first pair long before the second. `rsqrt` never
saturates on these three types, since `1/√ulp` is in range.

## Without a fixed threshold

These saturate wherever the exact result leaves the range.

| Function | Saturates when |
|----------|----------------|
| `pow` | `\|base\|^exponent` passes `MAX`, signed by the exponent's parity for a negative base; it is 0 below half an ulp |
| `mul_add`, `mul_div`, `dot`, `scale_by`, and the `_rounded` forms | The exact result is out of range |
| `div_or_saturate`, `div_rounded` | The quotient is out of range, or the divisor is 0, when the result takes the numerator's sign; 0/0 is 0 |
| `mul_div` | `c` is 0, when the result takes the sign of `a·b`, or is 0 if `a·b` is |
| `hypot`, `norm`, `Complex::norm` | The length passes `MAX`, as `hypot(MAX, MAX)` does |
| `Complex` arithmetic | A part leaves the range; each part saturates on its own |
| `Complex::exp`, `exprel`, `exprel2`, `exprel_weighted`, `pow` | A part passes `MAX` |
| `Complex::erf`, `faddeeva`, `e1`, `gamma_upper` | A part passes `MAX`: `erf` grows away from the real axis, `faddeeva` below it, `e1` along the negative axis, and `gamma_upper` for large `a` or along the negative axis |
| `enclose` and the `_range` functions | An end passes the range; `tan_ratio_range` returns the whole range when a pole may lie inside |

Every other total function stays in range on I16F16, I32F32 and I64F64.

## Narrow layouts

The fewer integer bits a layout has, the more functions reach its ends.
`I3F125` spans only `[−4, 4)`. There `exp` saturates from `ln 4`, `ln` below
`e^−4`, and `ellip_k` near `m = 1`. The crate accepts every signed layout
that can represent π, from `I8F0` to `I3F125`, so `acos` and `atan2` never
saturate. The formulas in the table give the thresholds on any layout.

## Domain errors

| Function | Domain |
|----------|--------|
| `asin`, `acos` | `\|x\| ≤ 1` |
| `asin_from_gap` | `gap` in `[0, 2]` |
| `acosh` | `x ≥ 1` |
| `atanh`, `atanh_ratio` | `\|x\| < 1` |
| `acoth` | `\|x\| > 1` |
| `coth` | `x ≠ 0` |
| `ln`, `log2`, `log10` | `x > 0` |
| `log1p`, `log1p_ratio`, `log1p_residual`, `log1p_residual_ratio` | `x > −1` |
| `pow` | Not a zero base with a negative exponent, nor a negative base with a non-integer exponent |
| `sqrt` | `x ≥ 0` |
| `rsqrt` | `x > 0` |
| `sqrt_product` | `a` and `b` not of opposite signs |
| `normalize` | Not the zero vector |
| `gamma` | `x > 0` |
| `ellip_k` | `m < 1` |
| `ellip_pi_complete` | `n < 1` and `m < 1` |
| `ellip_pi` | `1 − n·sn² > 0`, and `cn` and `dn` not both 0 |
| `ellip_pi_increment` | `n < 1` |
| `jacobi_inverse` | `cn` and `dn` not both 0 |
| `carlson_rf` | `x, y, z ≥ 0`, at most one of them 0 |
| `carlson_rf_roots` | At most one root 0 |
| `carlson_rj` | `x, y, z ≥ 0`, at most one of them 0, and `p > 0` |
| `carlson_rj_roots` | At most one of `rx`, `ry`, `rz` 0, and `rp ≠ 0` |
| `carlson_rc` | `x ≥ 0` and `y > 0` |
| `carlson_rd` | `x, y ≥ 0`, not both 0, and `z > 0` |
| `Complex::ln`, `Complex::e1` | `z ≠ 0` |
| `Complex::pow` | `z ≠ 0` when `a < 0` |
| `Complex::gamma_upper` | `a` in `(0, 128]` |
| `Complex::exprel_weighted` | `ε` in `[0, 1)` and `p ≥ 0` |
| `root_in` | `f(lo)` and `f(hi)` not both non-zero with one sign |
| `kepler_universal` | `r₀ > 0` and `r₀(2 − αr₀) ≥ σ₀²`, a real orbit |
| `log1p_range`, `log1p_ratio_range`, `log1p_residual_ratio_range` | Both ends above −1 |
| `atanh_ratio_range` | Both ends in `(−1, 1)` |

`Modulus::from_parameter` and `Modulus::from_complement` return `None`
outside `[0, 1]`. The wrappers in `bounded` check a domain once, when built,
so `sqrt_nonneg` and `atanh_open` return `T` with no check of their own.

## Detecting saturation

A result of `T::MAX` or `T::MIN` may be saturated, or may be the exact value
rounded there. Where the difference matters, bound the input by the table
above before the call.

```rust
use fixed::types::I16F16;
use fixed_analytics::exp;

let x = I16F16::from_num(5.0);
if x < I16F16::from_num(10.39) {
    let y = exp(x); // in range
}
```

The `_range` functions give an interval certain to hold the exact value, and
an end at `MAX` or `MIN` shows where it may have saturated.
