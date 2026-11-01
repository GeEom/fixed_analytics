# Accuracy

Every function returns a result within one ulp of the exact value, on every
layout the crate accepts, wherever that value is in range. `ERROR_ULPS` states
this bound. It is strict. Out of range a result saturates instead, as
[saturation.md](saturation.md) describes.

## How the bound holds

A function converts its arguments into binary floating point without loss and
rounds once at the end. The working mantissa has 64 bits for 8-, 16- and
32-bit types, 128 for 64-bit types and 192 for 128-bit types. Its error is
below 2^-20 ulp of the type, so a result misses the nearest value only when
the exact value lies within that distance of a midpoint.
[engine.md](engine.md) gives the design.

| Functions | Bound |
|-----------|-------|
| `mul_add`, `mul_div`, `div_or_saturate`, `dot`, `sqrt`, `sqrt_or_zero` | Half an ulp, the nearest value: the exact result is formed before the one rounding |
| `mul_rounded`, `div_rounded`, `mul_div_rounded`, `mul_add_rounded`, `dot_rounded` | Under one ulp, on the side `Rounding::Down` or `Rounding::Up` asks for; `Rounding::Nearest` gives the half-ulp result |
| `wrap_principal` | Under one ulp, and in `(−π, π]` even where the nearest value lies outside it |
| Everything else | Under one ulp |

Sign and magnitude are kept apart, which makes odd functions odd and even
functions even. An odd function that saturates at `T::MAX` meets `T::MIN` on
the other side, the one exception. Monotone functions give monotone results:
the final rounding is monotone, and the working error is far below the gap
between neighbouring results wherever a function is not flat.

## How it is measured

The accuracy tool in `tools/accuracy-bench` checks the bound against MPFR,
through the `rug` crate. It has two parts.

The worst-case sweep runs every function on 104 layouts. These are every
layout of 8, 16 and 32 bits, with 24 of 64 bits and 30 of 128 bits spread
across their fractional widths. Its inputs are

- one argument: every input on the 8- and 16-bit layouts, and 40000 seeded
  inputs on the wider ones;
- two arguments: every pair on the 8-bit layouts, 2^18 pairs on the 16-bit
  ones and 40000 on the wider ones;
- three or more: 2^18 tuples up to 16 bits and 40000 past them;
- the elliptic functions, the Fresnel integrals and the complex special
  functions and powers, which cost too much to sweep whole: about 3000
  inputs on every layout, and 375 for `exprel_weighted`.

MPFR computes each reference at 320 bits or more, and the error is taken in
ulps of the type. The sweep also tests symmetry, and tests monotone functions
for any reversal between successive inputs. On the wider layouts each
sampled input brings its next two neighbours for this.
[worst-case.md](worst-case.md) gives the largest errors.

The second part, the relative error table below, covers 21 functions on
I16F16 and I32F32. Each gets 5000 grid points, 50000 seeded random points and
2000 points crowded toward each end of its domain.

## Regressions

CI runs the tool on every pull request against the base branch's
`baseline.json`, and on every push to main against the previous commit's.
To fit a small runner it takes an eighth of the inputs and sweeps only the
8-bit layouts whole. It fails when a mean relative error rises by more than
0.1%, when any result breaks its bound, a symmetry or a monotonicity, or when
the table below is out of date. A full local run with `--worst-baseline` also
fails on any worst case that rises.

`tests/golden.txt` pins a hash of every function's outputs on six layouts. A
release that changes any output bit names the function in the changelog, even
for a one-ulp change.

## Running the tool

```sh
cd tools/accuracy-bench
cargo run --release                                # rewrites worst.txt and both docs pages
cargo run --release -- --baseline old.json         # compares, and checks the pages instead
cargo run --release -- --quick --baseline old.json # as CI runs it
cargo run --release --bin sweep_some -- exp ln     # the sweep for named functions only
```

`--worst-baseline old.txt` adds a comparison with an earlier `worst.txt`.
Each run writes its report to `reports/`.

## Relative error

Mean, median and 95th percentile relative error against MPFR. A result
rounded to nearest has a relative error of at most half an ulp over its own
magnitude, so these figures measure the type's resolution more than the
function. They are largest where results come close to zero: `exp`, `pow2`
and `pow` on the lower part of their domain, and `acoth` far from ±1.

<!-- ACCURACY_START -->
| Function | I16F16 Mean | I16F16 Median | I16F16 P95 | I32F32 Mean | I32F32 Median | I32F32 P95 |
|----------|-------------|---------------|------------|-------------|---------------|------------|
| sin | 3.11e-5 | 6.12e-6 | 6.18e-5 | 4.65e-10 | 9.35e-11 | 9.59e-10 |
| cos | 3.14e-5 | 6.17e-6 | 6.29e-5 | 4.73e-10 | 9.33e-11 | 9.89e-10 |
| tan | 3.28e-5 | 1.26e-5 | 7.51e-5 | 5.37e-10 | 1.91e-10 | 1.13e-9 |
| asin | 3.89e-5 | 1.12e-5 | 9.52e-5 | 6.56e-10 | 1.72e-10 | 1.50e-9 |
| acos | 8.39e-6 | 3.69e-6 | 2.85e-5 | 1.29e-10 | 5.64e-11 | 4.25e-10 |
| atan | 2.72e-6 | 2.52e-6 | 4.78e-6 | 4.19e-11 | 3.88e-11 | 7.29e-11 |
| sinh | 7.41e-6 | 4.17e-6 | 1.26e-5 | 1.18e-10 | 6.33e-11 | 1.91e-10 |
| cosh | 3.79e-6 | 3.73e-6 | 7.35e-6 | 5.76e-11 | 5.65e-11 | 1.12e-10 |
| tanh | 5.38e-6 | 2.10e-6 | 9.51e-6 | 1.09e-10 | 6.18e-11 | 1.41e-10 |
| coth | 3.53e-6 | 2.24e-6 | 8.26e-6 | 7.05e-11 | 5.99e-11 | 1.28e-10 |
| asinh | 2.85e-6 | 1.37e-6 | 4.96e-6 | 4.43e-11 | 2.08e-11 | 7.56e-11 |
| acosh | 2.29e-6 | 1.36e-6 | 4.80e-6 | 3.59e-11 | 2.09e-11 | 7.45e-11 |
| atanh | 4.22e-5 | 1.40e-5 | 1.01e-4 | 7.11e-10 | 2.13e-10 | 1.56e-9 |
| acoth | 1.92e-4 | 1.40e-4 | 5.49e-4 | 2.94e-9 | 2.16e-9 | 8.40e-9 |
| exp | 5.98e-3 | 9.19e-6 | 4.13e-2 | 9.45e-8 | 1.39e-10 | 6.50e-7 |
| ln | 1.02e-6 | 6.57e-7 | 1.46e-6 | 1.39e-11 | 1.01e-11 | 2.23e-11 |
| log2 | 6.76e-7 | 4.55e-7 | 1.01e-6 | 9.69e-12 | 6.89e-12 | 1.53e-11 |
| log10 | 1.92e-6 | 1.50e-6 | 3.35e-6 | 3.00e-11 | 2.30e-11 | 5.07e-11 |
| pow2 | 3.54e-4 | 5.03e-6 | 2.37e-3 | 5.47e-9 | 7.69e-11 | 3.67e-8 |
| pow | 6.47e-4 | 1.33e-5 | 3.13e-3 | 1.01e-8 | 2.02e-10 | 4.85e-8 |
| sqrt | 8.88e-8 | 5.80e-8 | 2.42e-7 | 1.37e-12 | 8.85e-13 | 3.62e-12 |
<!-- ACCURACY_END -->
