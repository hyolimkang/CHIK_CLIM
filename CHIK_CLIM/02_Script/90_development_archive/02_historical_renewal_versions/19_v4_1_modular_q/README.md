# Ceará renewal v4.1 -- modular q

All-age, no-vaccination weekly renewal model for Ceará, 2015-01-04 through
2025-12-21. This version propagates reporting-fraction (q) uncertainty
**modularly**, not jointly: q is fixed DATA within every conditional Stan
fit, and is varied ACROSS five separate fits (one per external-prior
quadrature node) to represent the external uncertainty in q. See
`DEVELOPMENT_LOG.md` for why this replaced an earlier joint-q attempt.

## Why this exists

An earlier v4.1 attempt (`../18_v4_1_minimal_no_vaccine/`) estimated q
jointly with transmission parameters and a learned annual-effect scale
`sigma_year`. That run became more than 5x slower than the proven,
HMC-passing v4.0/td14 benchmark and was aborted before diagnostics could be
obtained -- an operationally unacceptable outcome, regardless of whether it
would eventually have converged. This version deliberately keeps the proven
v4.0/td14 Stan model **byte-for-byte unchanged** (aside from renaming
`q_fixed` to `q_external`) and handles q uncertainty outside Stan instead.

## Pipeline

Run in order (see each script's header for details):

| script | purpose |
|---|---|
| `scripts/00_audit_v4_0_reference.R` | Verifies this model is byte-identical to the HMC-passing v4.0/td14 Stan source, aside from the q_fixed -> q_external rename. Run this FIRST and after any edit to the `.stan` file. |
| `scripts/01_build_q_nodes.R` | Monte Carlo (2e6 draws) of the induced q = p_symptomatic * p_detect distribution; builds the 5 stratified q nodes (`q_nodes.csv`). |
| `scripts/02_run_q_node.R` | Fits ONE q node (`Q_NODE_ID=q50 Rscript ...`). Same settings as v4.0/td14: 4 chains, 1000 warmup + 1000 sampling, adapt_delta=0.95, max_treedepth=14. Writes live per-chain CSVs to `outputs/<node>/chains/`. |
| `scripts/03_monitor_v4_1_fit.R` | Run in a SEPARATE terminal while a node is sampling. Reports per-chain progress/diagnostics every 60s from the live CSVs; writes `outputs/<node>/live_sampler_status.tsv`. Read-only, safe to interrupt. |
| `scripts/04_check_hmc_gate.R` | Post-hoc strict HMC gate check for one completed node. |
| `scripts/05_combine_modular_posterior.R` | After all 5 nodes pass, combines them into the modular ensemble (equal draws per node, weight 0.20 each). |
| `scripts/06_make_ppc.R` | Per-node and combined posterior predictive checks, with explicit 2017/2021/2022/2023 rows. |
| `scripts/07_q_sensitivity_plots.R` | Plots key outputs vs q to check whether 5 nodes are sufficient. |
| `scripts/08_audit_generation_interval.R` | Documentation-only audit of the generation interval `w`/`G` (unchanged from v4.0). |

## Fit order (do not deviate)

q50 first (canary) -> q10, q90 (test extremes) -> only if all three pass,
q30, q70. Do not launch all five blindly.

## What this is not

This is **not** a fully joint Bayesian posterior in q. The Ceará case
likelihood never updates q. See `DEVELOPMENT_LOG.md` and the header comment
in `stan/` (i.e. `../../stan/renewal_ceara_v4_1_modular_q.stan`) for the
exact interpretation.

## No new parameters

This version does not add: jointly sampled q, sigma_year, time-varying/
yearly/weekly q, AR(1), weekly latent R, climate coefficients, serology
likelihood, fitted importation, extra seasonal harmonics, age structure, or
vaccination.
