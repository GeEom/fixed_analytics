Minimax (Chebyshev) polynomial coefficients for sin/cos evaluation.

The coefficients come from the Remez algorithm (`mpmath.chebyfit`). Each
table approximates a reduced remainder, factored so every coefficient fits in
I1F63 (magnitude below 1).

Tables are stored highest degree first, `[cₙ, cₙ₋₁, …, c₀]`, and evaluated by
Horner's method:

```text
P(x) = cₙ·xⁿ + ··· + c₁·x + c₀ = c₀ + x·(c₁ + x·(c₂ + ··· + x·cₙ))
```

# Reconstruction

| Function | Factored form | Variable | Domain |
|----------|---------------|----------|--------|
| sin | `x + x³·P(x²)` | u = x² | \[0, (π/4)²\] |
| cos | `1 + x²·Q(x²)` | u = x² | \[0, (π/4)²\] |

The `_LO` tables suit I16F16-class precision (fewer than 24 fractional bits)
and the `_HI` tables I32F32-class and wider.
