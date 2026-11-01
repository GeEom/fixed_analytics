# Worst-case error

The largest error the accuracy tool has found for each function, in ulps of
the type, against MPFR on exact inputs. Each row gives three common layouts
and the largest over all 104. [accuracy.md](accuracy.md#how-it-is-measured)
describes the inputs, and `tools/accuracy-bench/worst.txt` holds the row for
every layout, with the input behind each maximum.

Rows ending `_re` and `_im` are the parts of a `Complex` result. `dot3`,
`norm3` and `normalize3` take three components. `stumpff_c0` to `stumpff_c3`,
`exp_moments_e1` to `exp_moments_e5`, `jacobi_sn` to `jacobi_dn` and
`fresnel_c` and `fresnel_s` are elements of an array or tuple result, and
`sin_turns` and `cos_turns` those of `sin_cos_turns`. `ellip_k` builds its
`Modulus` from `k'` and `ellip_k_from_m` from `m`. `exprel_weighted_k1` and
`_k2` are the two results of `Complex::exprel_weighted`.
`carlson_rf_rd_roots_rf` and `_rd`, `carlson_rd_rj_roots_rd` and `_rj`, and
`carlson_rf_rd_rj_roots_rf` to `_rj` are the results of the joint Carlson
functions, and `sin_ratio_versine` is the second of `sin_ratio`.
`kepler_universal_stumpff_c0` to `_c3` are the Stumpff values that function
returns, measured at its own `χ`, which is `kepler_universal`'s. The `_down`
and `_up` rows are the `_rounded` functions with `Rounding::Down` and
`Rounding::Up`. Their error approaches a whole ulp and the table rounds it to
1.00, but the sweep checks that each is below one and on its side.

<!-- WORST_START -->
488720248 inputs over 135 functions and 104 layouts; 0 of 102948745 symmetry checks and 0 of 253753023 monotonicity checks failed.

| Function | I16F16 | I32F32 | I64F64 | All layouts |
|----------|------|------|------|------|
| sin | 0.50 | 0.50 | 0.50 | 0.50 |
| cos | 0.50 | 0.50 | 0.50 | 0.50 |
| tan | 0.50 | 0.50 | 0.50 | 0.50 |
| asin | 0.50 | 0.50 | 0.50 | 0.50 |
| acos | 0.50 | 0.50 | 0.50 | 0.50 |
| atan | 0.50 | 0.50 | 0.50 | 0.50 |
| atan2 | 0.50 | 0.50 | 0.50 | 0.50 |
| sinh | 0.50 | 0.50 | 0.50 | 0.50 |
| cosh | 0.50 | 0.50 | 0.50 | 0.50 |
| tanh | 0.50 | 0.50 | 0.50 | 0.50 |
| coth | 0.50 | 0.50 | 0.50 | 0.50 |
| asinh | 0.50 | 0.50 | 0.50 | 0.50 |
| acosh | 0.50 | 0.50 | 0.50 | 0.50 |
| atanh | 0.50 | 0.50 | 0.50 | 0.50 |
| acoth | 0.50 | 0.50 | 0.50 | 0.50 |
| exp | 0.50 | 0.50 | 0.50 | 0.50 |
| pow2 | 0.50 | 0.50 | 0.50 | 0.50 |
| ln | 0.50 | 0.50 | 0.50 | 0.50 |
| log2 | 0.50 | 0.50 | 0.50 | 0.50 |
| log10 | 0.50 | 0.50 | 0.50 | 0.50 |
| pow | 0.50 | 0.50 | 0.50 | 0.50 |
| sqrt | 0.50 | 0.50 | 0.50 | 0.50 |
| sqrt_or_zero | 0.50 | 0.50 | 0.50 | 0.50 |
| mul_add | 0.50 | 0.50 | 0.50 | 0.50 |
| mul_div | 0.50 | 0.50 | 0.50 | 0.50 |
| div_or_saturate | 0.50 | 0.50 | 0.50 | 0.50 |
| dot3 | 0.50 | 0.50 | 0.50 | 0.50 |
| hypot | 0.50 | 0.50 | 0.50 | 0.50 |
| norm3 | 0.50 | 0.50 | 0.50 | 0.50 |
| normalize3 | 0.50 | 0.50 | 0.50 | 0.50 |
| scale_by | 0.50 | 0.50 | 0.50 | 0.50 |
| sqrt_product | 0.50 | 0.50 | 0.50 | 0.50 |
| asin_from_gap | 0.50 | 0.50 | 0.50 | 0.50 |
| rsqrt | 0.50 | 0.50 | 0.50 | 0.50 |
| expm1 | 0.50 | 0.50 | 0.50 | 0.50 |
| exprel | 0.50 | 0.50 | 0.50 | 0.50 |
| exprel2 | 0.50 | 0.50 | 0.50 | 0.50 |
| exp_moments_e1 | 0.50 | 0.50 | 0.50 | 0.50 |
| exp_moments_e2 | 0.50 | 0.50 | 0.50 | 0.50 |
| exp_moments_e3 | 0.50 | 0.50 | 0.50 | 0.50 |
| exp_moments_e4 | 0.50 | 0.50 | 0.50 | 0.50 |
| exp_moments_e5 | 0.50 | 0.50 | 0.50 | 0.50 |
| log1p | 0.50 | 0.50 | 0.50 | 0.50 |
| log1p_ratio | 0.50 | 0.50 | 0.50 | 0.50 |
| log1p_residual | 0.50 | 0.50 | 0.50 | 0.50 |
| log1p_residual_ratio | 0.50 | 0.50 | 0.50 | 0.50 |
| stumpff_c0 | 0.50 | 0.50 | 0.50 | 0.50 |
| stumpff_c1 | 0.50 | 0.50 | 0.50 | 0.50 |
| stumpff_c2 | 0.50 | 0.50 | 0.50 | 0.50 |
| stumpff_c3 | 0.50 | 0.50 | 0.50 | 0.50 |
| sin_ratio | 0.50 | 0.50 | 0.50 | 0.50 |
| sin_ratio_versine | 0.50 | 0.50 | 0.50 | 0.50 |
| tan_ratio | 0.50 | 0.50 | 0.50 | 0.50 |
| atan_ratio | 0.50 | 0.50 | 0.50 | 0.50 |
| atanh_ratio | 0.50 | 0.50 | 0.50 | 0.50 |
| tanh_ratio | 0.50 | 0.50 | 0.50 | 0.50 |
| ellip_k | 0.50 | 0.50 | 0.50 | 0.50 |
| ellip_k_from_m | 0.50 | 0.50 | 0.50 | 0.50 |
| ellip_e | 0.50 | 0.50 | 0.50 | 0.50 |
| jacobi_sn | 0.50 | 0.50 | 0.50 | 0.50 |
| jacobi_cn | 0.50 | 0.50 | 0.50 | 0.50 |
| jacobi_dn | 0.50 | 0.50 | 0.50 | 0.50 |
| jacobi_am | 0.50 | 0.50 | 0.50 | 0.50 |
| jacobi_inverse | 0.50 | 0.50 | 0.50 | 0.50 |
| carlson_rf | 0.50 | 0.50 | 0.50 | 0.50 |
| carlson_rf_roots | 0.50 | 0.50 | 0.50 | 0.50 |
| carlson_rj | 0.50 | 0.50 | 0.50 | 0.50 |
| carlson_rc | 0.50 | 0.50 | 0.50 | 0.50 |
| carlson_rd | 0.50 | 0.50 | 0.50 | 0.50 |
| carlson_rj_roots | 0.50 | 0.50 | 0.50 | 0.50 |
| carlson_rd_roots | 0.50 | 0.50 | 0.50 | 0.50 |
| carlson_rf_rd_roots_rf | 0.50 | 0.50 | 0.50 | 0.50 |
| carlson_rf_rd_roots_rd | 0.50 | 0.50 | 0.50 | 0.50 |
| carlson_rd_rj_roots_rd | 0.50 | 0.50 | 0.50 | 0.50 |
| carlson_rd_rj_roots_rj | 0.50 | 0.50 | 0.50 | 0.50 |
| carlson_rf_rd_rj_roots_rf | 0.50 | 0.50 | 0.50 | 0.50 |
| carlson_rf_rd_rj_roots_rd | 0.50 | 0.50 | 0.50 | 0.50 |
| carlson_rf_rd_rj_roots_rj | 0.50 | 0.50 | 0.50 | 0.50 |
| ellip_pi | 0.50 | 0.50 | 0.50 | 0.50 |
| ellip_pi_complete | 0.50 | 0.50 | 0.50 | 0.50 |
| sin_turns | 0.50 | 0.50 | 0.50 | 0.50 |
| cos_turns | 0.50 | 0.50 | 0.50 | 0.50 |
| atan2_turns | 0.50 | 0.50 | 0.50 | 0.50 |
| div_by_tau | 0.50 | 0.50 | 0.50 | 0.50 |
| wrap_radians | 0.50 | 0.50 | 0.50 | 0.50 |
| wrap_principal | 0.50 | 0.50 | 0.50 | 0.98 |
| mul_by_tau | 0.50 | 0.50 | 0.50 | 0.50 |
| c_mul_re | 0.50 | 0.50 | 0.50 | 0.50 |
| c_mul_im | 0.50 | 0.50 | 0.50 | 0.50 |
| c_div_re | 0.50 | 0.50 | 0.50 | 0.50 |
| c_div_im | 0.50 | 0.50 | 0.50 | 0.50 |
| c_exp_re | 0.50 | 0.50 | 0.50 | 0.50 |
| c_exp_im | 0.50 | 0.50 | 0.50 | 0.50 |
| c_exprel_re | 0.50 | 0.50 | 0.50 | 0.50 |
| c_exprel_im | 0.50 | 0.50 | 0.50 | 0.50 |
| c_exprel2_re | 0.50 | 0.50 | 0.50 | 0.50 |
| c_exprel2_im | 0.50 | 0.50 | 0.50 | 0.50 |
| faddeeva_re | 0.50 | 0.50 | 0.50 | 0.50 |
| faddeeva_im | 0.50 | 0.50 | 0.50 | 0.50 |
| c_erf_re | 0.50 | 0.50 | 0.50 | 0.80 |
| c_erf_im | 0.50 | 0.50 | 0.50 | 0.50 |
| e1_re | 0.50 | 0.50 | 0.50 | 0.50 |
| e1_im | 0.50 | 0.50 | 0.50 | 0.50 |
| gamma_upper_re | 0.50 | 0.50 | 0.50 | 0.50 |
| gamma_upper_im | 0.50 | 0.50 | 0.50 | 0.50 |
| c_ln_re | 0.50 | 0.50 | 0.50 | 0.50 |
| c_ln_im | 0.50 | 0.50 | 0.50 | 0.50 |
| c_pow_re | 0.50 | 0.50 | 0.50 | 0.50 |
| c_pow_im | 0.50 | 0.50 | 0.50 | 0.50 |
| exprel_weighted_k1_re | 0.50 | 0.50 | 0.50 | 0.50 |
| exprel_weighted_k1_im | 0.50 | 0.50 | 0.49 | 0.50 |
| exprel_weighted_k2_re | 0.50 | 0.50 | 0.49 | 0.50 |
| exprel_weighted_k2_im | 0.50 | 0.50 | 0.50 | 0.50 |
| kepler_universal | 0.50 | 0.50 | 0.50 | 0.50 |
| kepler_universal_stumpff_c0 | 0.50 | 0.50 | 0.50 | 0.50 |
| kepler_universal_stumpff_c1 | 0.50 | 0.50 | 0.50 | 0.50 |
| kepler_universal_stumpff_c2 | 0.50 | 0.50 | 0.50 | 0.50 |
| kepler_universal_stumpff_c3 | 0.50 | 0.50 | 0.50 | 0.50 |
| kepler_universal_near | 0.50 | 0.50 | 0.50 | 0.50 |
| erf | 0.50 | 0.50 | 0.50 | 0.50 |
| erfc | 0.50 | 0.50 | 0.50 | 0.50 |
| fresnel_c | 0.50 | 0.50 | 0.50 | 0.50 |
| fresnel_s | 0.50 | 0.50 | 0.50 | 0.50 |
| gamma | 0.50 | 0.50 | 0.50 | 0.50 |
| mul_down | 1.00 | 1.00 | 1.00 | 1.00 |
| mul_up | 1.00 | 1.00 | 1.00 | 1.00 |
| mul_div_down | 1.00 | 1.00 | 1.00 | 1.00 |
| mul_div_up | 1.00 | 1.00 | 1.00 | 1.00 |
| mul_add_down | 1.00 | 1.00 | 1.00 | 1.00 |
| mul_add_up | 1.00 | 1.00 | 1.00 | 1.00 |
| div_down | 1.00 | 1.00 | 1.00 | 1.00 |
| div_up | 1.00 | 1.00 | 1.00 | 1.00 |
| dot3_down | 1.00 | 1.00 | 1.00 | 1.00 |
| dot3_up | 1.00 | 1.00 | 1.00 | 1.00 |
| ellip_pi_increment | 0.50 | 0.50 | 0.50 | 0.50 |

Past the nearest value by more than the working error, 2^-20 ulp:

| Function | Layout | Ulps | Input |
|----------|--------|------|-------|
| wrap_principal | I6F2 | 0.53 | -2.2e1 |
| wrap_principal | I14F2 | 0.57 | -3.55e2 |
| wrap_principal | I11F5 | 0.53 | -3.55e2 |
| wrap_principal | I6F10 | 0.98 | -3.1416015625e0 |
| wrap_principal | I5F11 | 0.96 | -3.1416015625e0 |
| wrap_principal | I4F12 | 0.93 | -3.1416015625e0 |
| wrap_principal | I3F13 | 0.85 | -3.1416015625e0 |
| wrap_principal | I30F2 | 0.57 | -3.55e2 |
| wrap_principal | I27F5 | 0.53 | -3.48285538125e7 |
| wrap_principal | I22F10 | 0.98 | 3.1416015625e0 |
| wrap_principal | I21F11 | 0.96 | 3.1416015625e0 |
| wrap_principal | I20F12 | 0.93 | 3.1416015625e0 |
| wrap_principal | I19F13 | 0.85 | 3.1416015625e0 |
| wrap_principal | I18F14 | 0.71 | 3.1416015625e0 |
| wrap_principal | I15F17 | 0.66 | 3.1415939331054688e0 |
| wrap_principal | I62F2 | 0.57 | 3.55e2 |
| wrap_principal | I52F12 | 0.93 | 3.1416015625e0 |
| wrap_principal | I8F56 | 0.65 | 3.141592653589793e0 |
| wrap_principal | I126F2 | 0.57 | -2.485e3 |
| wrap_principal | I72F56 | 0.65 | 3.141592653589793e0 |
| wrap_principal | I66F62 | 0.54 | 3.141592653589793e0 |
| wrap_principal | I10F118 | 0.63 | 3.141592653589793e0 |
| c_erf_re | I80F48 | 0.80 | -1.170949760664719e-1,-5.469643837406263e0 |
<!-- WORST_END -->

## The exceptions

Two functions miss the nearest value by more than the working error, both
within `ERROR_ULPS`.

`wrap_principal` returns a value in `(−π, π]`. Next to ±π the nearest value
of the type can fall outside, and the function returns the neighbour inside
instead. On I6F10, `−3.1416015625` wraps to `3.14158374…`, whose nearest
value `3.1416015625` is past π, so the result is `3.140625`, 0.98 ulp away.

The real part of `Complex::erf` reaches 0.80 ulp at one input on I80F48.
