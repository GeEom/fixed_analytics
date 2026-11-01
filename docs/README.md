# Documentation

For users of the crate:

- [accuracy.md](accuracy.md): the error bound, how the accuracy tool measures
  it, and relative error per function.
- [worst-case.md](worst-case.md): the largest error of every function, on
  every layout swept.
- [saturation.md](saturation.md): where each function saturates or fails,
  with thresholds for I16F16, I32F32 and I64F64.

The tool keeps the tables in the first two up to date. The API reference is
on [docs.rs](https://docs.rs/fixed_analytics), built in part from the pages
below.

| Page | Shown as |
|------|----------|
| [crate.md](crate.md) | The crate root |
| [exponential.md](exponential.md) | `ops::exponential` |
| [hyperbolic.md](hyperbolic.md) | `ops::hyperbolic` |
| [kernel.md](kernel.md) | `kernel`, the CORDIC kernels |
| [chebyshev.md](chebyshev.md) | `tables::chebyshev` |

For maintainers, [engine.md](engine.md) records how the functions evaluate,
and [sqrt.md](sqrt.md) derives the 128-bit square root. Neither ships in the
package.
