# Ceará renewal-model development log -- v4.1 modular q

See `../17_v4_0_minimal_no_vaccine/MODEL_DEVELOPMENT.md` for v4.0 (base run
unconverged at max_treedepth=12; td14 control run PASSED the strict HMC gate
and is the proven reference used here).

## Why the joint-q attempt was abandoned

`../18_v4_1_minimal_no_vaccine/` attempted to simultaneously (a) estimate a
single constant reporting fraction q as a Stan parameter with an informative
Beta(16.2, 108.7) prior, and (b) replace the fixed `year_effect_prior_sd`
with a learned `sigma_year ~ half-normal(0, 0.35)` via an explicit
orthonormal sum-to-zero annual-effect basis. That run became more than 5x
slower than the HMC-passing v4.0/td14 benchmark (which took ~1478s) and was
aborted after ~2.7 hours of wall-clock time without producing any usable
diagnostics, because the fit script did not write incremental per-iteration
output (`sample_file` was added only afterward, and a restart with it
enabled also showed early divergences and treedepth-14 hits during warmup
before being superseded by this modular redesign).

**This is NOT evidence that q is definitely non-identifiable jointly with
transmission parameters and sigma_year.** Two structural changes were
introduced at once, so the slowdown cannot be attributed to q alone, to
sigma_year alone, or to their interaction. What IS established is that the
joint attempt was operationally impractical before any diagnostic could be
obtained, and that this is unacceptable for a working backbone that must
support iterative development and, eventually, vaccine-impact inference.

## Why the modular design was chosen

- v4.0/td14 already demonstrated stable HMC (divergences=0, treedepth
  hits=0, max Rhat=1.00897-1.00566 depending on which parameter subset is
  summarised, min bulk ESS ~473-481, BFMI in [0.926, 1.083]) with q fixed as
  data at 0.10.
- q's uncertainty is genuinely external/epidemiological (a product of two
  independently-sourced Beta priors for symptomatic proportion and
  detection-given-symptomatic), not something the Ceará case series alone
  can be expected to identify well -- q trades off against the transmission
  scale (alpha_R), the annual effects, and the entire latent infection/
  susceptible-depletion trajectory in a way that plausibly creates a long,
  difficult ridge in the joint posterior.
- Modular propagation keeps q's genuine external uncertainty in the final
  answer while never asking NUTS to explore that ridge: each conditional fit
  only has to explore `p(theta | cases, q = q_k)`, which is exactly the
  already-proven-tractable v4.0/td14 geometry, for five fixed values of q_k.

## What "modular" means here, precisely

For each node k (q10, q30, q50, q70, q90), Stan estimates
`p(theta | cases, q = q_k)` where theta = (alpha_R, year_effect, beta_sin,
beta_cos, phi_obs) and all derived quantities (X, S, U, R0_t, R_eff_t, ...).
The five fits are combined post-hoc into

    p_modular(theta, q | cases) = p_external(q) * p(theta | cases, q)

by drawing equal numbers of posterior draws from each node's fit and
tagging each draw with its q value (`05_combine_modular_posterior.R`). This
is a **numerical quadrature approximation to a continuous external
distribution**, not a fully joint Bayesian posterior -- the case likelihood
is deliberately prevented from updating q. Report this explicitly as
"external reporting-fraction uncertainty propagated through conditional
transmission fits" / "modular uncertainty propagation," never as "the joint
posterior of q and transmission parameters."

## Fit records

(Filled in as each node completes -- see `outputs/<node>/hmc_diagnostics_<node>.csv`,
`outputs/<node>/runtime_<node>.csv`.)

- **q50** (q_external = 0.12708, the median node): **STRICT HMC GATE: FAIL.**
  Divergences = 0, max-treedepth(14) hits = 0, runtime normal (2472s, 1.67x
  the td14 benchmark of 1478s -- not flagged COMPUTATIONALLY_EXPENSIVE). But
  max Rhat = 1.535 and min bulk ESS = 7.2 (target <=1.01 / >=100) across
  essentially every parameter except beta_sin (Rhat=1.037, ESS=94, still
  failing). Root cause identified precisely: this is a 3-vs-1 CHAIN SPLIT,
  not a diffuse mixing problem. Chains 2-4 converged to
  alpha_R~0.108, beta_cos~0.255, phi_obs~1.81, z_year[1]~-0.12 -- closely
  matching the td14 reference (alpha_R~0.146-0.147, beta_cos~0.251,
  phi_obs~1.85-1.86 across all 4 chains at q_fixed=0.10). Chain 1 instead
  converged to a qualitatively different, implausible regime:
  alpha_R=2.627 (implied baseline R0=exp(2.627)=13.8), phi_obs=0.644 (much
  heavier NB overdispersion), z_year[1]=-6.35 (an enormous compensating
  negative annual deviation for 2015) -- each internally smooth (hence zero
  divergences/treedepth hits within that chain) but a completely different
  basin from chains 2-4.
  **This is NOT a prior-tail-q artifact**: q50=0.127 is the median node, and
  the td14 comparison shows ALL FOUR chains agreed tightly at q_fixed=0.10
  with the exact same initialisation strategy (offsets
  -0.06/-0.02/0.02/0.06 by chain id) -- chain 1's -0.06 offset did not reach
  this alternate mode at q=0.10, but did at q=0.127, a ~27% change in a
  single fixed scalar. The renewal recursion's likelihood surface therefore
  appears to have (at least) a second local mode whose basin of attraction
  is sensitive to the fixed q value, reachable from at least one of the four
  prespecified initial offsets. Output preserved at `outputs/q50/`, NOT
  deleted; figure saved labelled `..._UNCONVERGED.pdf` per the established
  v4.0 convention. Per Section 19, stopping here to report before running
  q10/q90 or attempting any redesign.
- **q10**: not yet run (blocked pending review of q50's failure).
- **q90**: not yet run (blocked pending review of q50's failure).
- **q30**: not yet run.
- **q70**: not yet run.

## Decision-rule status

q50 (the canary) FAILED the strict HMC gate for a reason unrelated to q
being a prior-tail value. Per Section 19's explicit instruction: do NOT
delete this node's output (preserved), do NOT loosen the gate (not done),
report back before redesigning the q approximation (this entry + the
conversation report is that report). v4.1 modular-q cannot yet be called a
STABLE MODULAR HISTORICAL BACKBONE.
