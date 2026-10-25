CORDIC (Coordinate Rotation Digital Computer) kernels.

Users should call the functions in [`crate::ops`] rather than these kernels.

# Algorithm

CORDIC rotates a vector through a fixed sequence of angles using only shifts
and adds. Each iteration performs one micro-rotation:

```text
x' = x - σ·d·y·2^(-i)
y' = y + σ·x·2^(-i)
z' = z - σ·angle[i]
```

where

- σ = ±1 is the direction of rotation,
- d = +1 for circular and −1 for hyperbolic mode,
- angle\[i\] = atan(2^-i) for circular and atanh(2^-i) for hyperbolic mode.

**Vectoring mode** drives y toward zero while z accumulates the angle:

| Mode | Vectoring (y → 0) |
|------|-------------------|
| Circular | atan |
| Hyperbolic | atanh, ln |

Hyperbolic mode uses the `atanh(2^-i)` table and repeats iterations 4, 13,
40, … so that the sequence converges.

# Termination

Both kernels stop after about a third of the fractional bits and close the
residual angle with one division: for `|w| = |y/x| < 2^-⌈frac/3⌉`, `atan(w)`
and `atanh(w)` differ from `w` by under `|w|³/3 < 2^-frac/3`.
