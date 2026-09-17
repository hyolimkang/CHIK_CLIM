# Pernambuco Curved-Ridge Reparameterisation Test (Model D)

**Status: computational-coordinate experiment only. No changes to the scientific model, priors, data, or serology.**

## Result: **QUADRATIC PILOT FAIL**

## 1. Geometry pre-check (before touching Stan)

Using non-divergent post-warmup draws from the existing affine-reparameterised
full run (Model C, adapt_delta=0.95, 26 divergences), re-centred at
`u0 = median(logit_q) = -4.608278` (NOT the old affine reference
`qlogis(0.10) = -2.197`, which the user correctly identified as merely a
computational coordinate constant with no scientific status):

| Quantity | Value |
|---|---|
| Free quadratic fit | `alpha_R = 0.4185 - 0.5459*d + 0.3284*d^2` |
| RMSE (linear only -> quadratic) | 0.03349 -> 0.02654 (**20.7% reduction**) |
| Residual correlation with d | -0.0000 |
| Residual correlation with d^2 | -0.0000 |
| Binned mean range (9 bins) vs overall residual SD | 0.0119 vs 0.0265 (ratio 0.45) |

**GATE: PASS.** `PE_alphaR_logitq_linear_vs_quadratic.png` shows the free
quadratic curve tracking the alpha_R-logit_q relationship almost exactly
across the full observed range, where the old fixed affine slope
(b1=-0.45) visibly diverges from the data at both tails.

## 2. Curved-ridge implementation

`renewal_pernambuco_v4_9_global_q_curved_reparam.stan`: `alpha_R = c0 +
gamma_curve + b1*d + b2*d^2`, `d = logit_q - u0`, with `c0=0.418533,
u0=-4.608278, b1=-0.545921, b2=0.328419` supplied as fixed DATA (not
estimated in Stan). Original priors on `alpha_R` and `logit_q` unchanged;
no independent prior on `gamma_curve`. Unit Jacobian confirmed
(`d(alpha_R)/d(gamma_curve) = 1`). Y-1 orthonormal year-effect basis
retained unchanged from Model C.

## 3. Pilot (4 chains, 300 warmup, 300 sampling, adapt_delta=0.95, max_treedepth=14, dense_e)

| Chain | Divergent | Max-treedepth hits | Mean stepsize | Mean leapfrog | Mean accept |
|---|---|---|---|---|---|
| 1 | 0 | 0 | 0.00571 | 954.7 | 0.985 |
| 2 | 0 | 0 | 0.00516 | 1046.9 | 0.980 |
| 3 | 2 | 0 | 0.00550 | 992.6 | 0.967 |
| **4** | **30** | 0 | 0.00563 | 821.1 | 0.899 |

Overall: **32 divergences**, max Rhat=1.04, min bulk ESS=90, BFMI
0.77-0.99 (all chains healthy). `cor(logit_q, gamma_curve) = -0.087` --
the near-linear correlation is excellently resolved (affine model:
-0.16; original unreparameterised: -0.95).

**Posterior location is scientifically unchanged from Model C** (confirms
the unit-Jacobian transform preserved the target):

| Quantity | Model D (curved) | Model C (affine) |
|---|---|---|
| q (median) | 0.0097 | 0.0099 |
| immune_2025 | 0.564 [0.356, 0.756] | 0.554 [0.342, 0.738] |
| eta_geo (median) | 1.008 | 1.025 |
| U14 predicted (geo-adjusted) | 0.348 | 0.350 (observed 0.372) |

## 4. Why this pilot FAILS the pre-registered PASS rule

`PE_gamma_curve_pilot_diagnostics.png` shows that chain 4's 30 divergences
are **not scattered across the typical posterior range**. They form a
**tight, isolated cluster** at `logit_q ~= -4.9, gamma_curve ~= 0.12,
alpha_R ~= 0.70` -- visibly offset from, not embedded within, the main
posterior cloud at that same logit_q value (where non-divergent draws sit
at gamma_curve ~= -0.02 to 0.03). This is a different signature from a
single under-adapted chain with a globally miscalibrated step size (chain
4's stepsize, 0.00563, is in fact comparable to the other three chains'):
it looks like the chain became trapped in, or repeatedly attempted and
failed to cross into/out of, a narrow, locally isolated region -- i.e.
either residual local multimodality or an additional narrow pathological
pocket that the quadratic transform does not remove.

This does not satisfy the pre-registered PASS conditions ("no obvious new
funnel", "divergences = 0"). Per instruction, this is reported as a FAIL
without automatically escalating to adapt_delta=0.99 or a higher-order
polynomial/spline reparameterisation.

## 5. Outputs

Figures (`03_Output/figures/pernambuco_v4_9_replication/modelD_curved_reparam/`):
`PE_alphaR_logitq_linear_vs_quadratic.png`,
`PE_gamma_curve_logitq_geometry.png`,
`PE_curved_ridge_divergence_map.png`,
`PE_gamma_curve_pilot_diagnostics.png`.

Tables (`03_Output/tables/pernambuco_v4_9_curved_reparam/`):
`PE_curved_ridge_regression_fit.csv` (superseded diagnostic -- retained
existing b1, later abandoned per user instruction),
`PE_free_quadratic_ridge_fit.csv`, `PE_free_quadratic_ridge_binned.csv`,
`PE_free_quadratic_ridge_constants.rds`.

Fit: `outputs/modelD/pe_curved_U14_pilot.rds`.

## 6. Stop rule (honoured)

No automatic escalation to adapt_delta=0.99, no higher-order (cubic/spline)
reparameterisation, no change to priors, model structure, serology, or
seasonality. Returned for review before any further computational or
scientific model changes.
