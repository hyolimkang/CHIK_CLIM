# Pernambuco v4.9 global-q ridge reparameterisation

## Purpose

`renewal_pernambuco_v4_9_global_q_ridge_reparam.stan` is a new Stan source
for a computational test of the frozen v4.9 Pernambuco global-q model. It
does **not** add dynamic transmission, age-stratified transmission, changing
reporting, new serology parameters, or new biological mechanisms.

The preceding global-q fits had strong posterior q--alpha_R correlation and
failed HMC diagnostics. The reparameterised model tests whether this is an
avoidable coordinate-geometry problem before concluding that absolute PE
susceptibility is not identifiable.

## Exact changes

1. The former `z_year - mean(z_year)` representation is replaced by an
   orthonormal Y-1 sum-to-zero contrast basis. This preserves the exact
   induced distribution of annual effects but removes the unused mean
   direction.
2. `alpha_R` is represented as

   `alpha_R = gamma_ridge + b * (logit_q - logit_q_reference)`.

   `b = -0.45` is a fixed numerical preconditioner obtained from the prior
   PE case-only pilot's local posterior ridge. The original independent
   priors are still evaluated on `alpha_R` and `logit_q`; the map has a unit
   Jacobian. Thus `b` cannot change the likelihood, prior, or posterior.
3. Chains start at q = 0.006, 0.010, 0.020 and 0.050 along the rotated ridge.
   This is an intentional multimodal/initialisation diagnostic, not a
   scientific prior on q.

## Interpretation rule

A clean HMC run establishes only that the existing global-q posterior can be
sampled. It does not establish that q or absolute susceptibility is identified.
The high-depletion PE interpretation requires a clean fit, q-prior
sensitivity, and an explicit statement that the Recife U14 survey is a
local/age-specific transport check rather than a direct state-level q anchor.

## Outputs

- Fits: `33_pernambuco_v4_9_replication/outputs/modelC_global_q_ridge_reparam/`
- Tables: `03_Output/tables/pernambuco_v4_9_ridge_reparam/`
