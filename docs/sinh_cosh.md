# How `sinh_cosh` is computed

This page records the derivation behind the comments in
`src/ops/hyperbolic.rs`.

The argument is halved until `|x| ≤ 1.1182`, the sum of the atanh table, and
both functions come from factored Taylor series in `u = x²`:

```text
sinh(x) = x·(1 + u/6·(1 + u/20·(1 + u/42·…)))
cosh(x) = 1 + u/2·(1 + u/12·(1 + u/30·…))
```

The divisors are `(2k)(2k+1)` for sinh and `(2k−1)(2k)` for cosh. Integer
division computes each `u/K` in one step without pre-rounding a reciprocal,
and the factored form keeps each multiplicand `u/K` below 1 for `K ≥ 2`.

At `|x| ≤ 1.1182` the omitted relative term is:

| Degrees (sinh / cosh) | Omitted term | Below half an ulp up to |
|-----------------------|--------------|-------------------------|
| 9 / 10 | 7.7e-8 / 8.0e-9 | 22 fractional bits |
| 13 / 14 | 3.7e-12 / 2.9e-13 | 37 fractional bits |
| 21 / 22 | 4.5e-22 / 2.4e-23 | 70 fractional bits |

The code uses the first row below 24 fractional bits, the second from 24 and
the third from 40. The results are rebuilt by doubling,
`sinh(2x) = 2·sinh(x)·cosh(x)` and `cosh(2x) = cosh²(x) + sinh²(x)`, once
per halving.
