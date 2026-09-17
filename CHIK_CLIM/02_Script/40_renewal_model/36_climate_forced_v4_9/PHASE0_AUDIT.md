# Phase 0 Audit — Climate-Forced v4.9 Extension

Purpose: identify the exact frozen v4.9 architecture, the currently
accepted state reconstructions, and the national climate dataset, before
any new Stan code is written. **No existing model or file was modified.**

## 1. The frozen v4.9 transmission core (two data-interface variants, ONE structure)

Two Stan files share a **byte-identical transmission core** — same
generation-interval convolution, same `log_R0[t]` formula, same
demographic S/U recursion with lifelong immunity, same NB2 case
likelihood, same transmission priors. They differ ONLY in how `q` and
serology enter as data/parameters:

| file | q | serology |
|---|---|---|
| `02_Script/stan/renewal_ceara_v4_9_hierarchical_seasonality.stan` | **fixed** (`real<lower=0,upper=1> q` in `data`) | single anchor (Juazeiro), beta-binomial |
| `02_Script/stan/renewal_bahia_v4_9_global_q_multisite_serology.stan` | **estimated** (`logit_q` parameter, `logit_q ~ normal(logit_q_prior_mean, logit_q_prior_sd)`) | 0..J multi-site, with a fixed geographic offset `eta_geo` per site |

The exact shared `log_R0[t]` line (present verbatim in both files) that
Phase 2 will extend:

```stan
log_R0[t] = alpha_R + year_effect[year_id[t]] +
  A_year[year_id[t]] * (beta_sin1 * seasonal_sin[t] + beta_cos1 * seasonal_cos[t]) +
  beta_sin2 * seasonal_sin2[t] + beta_cos2 * seasonal_cos2[t];
R0_t[t] = exp(log_R0[t]);
```

`A_year[y] = exp(sigma_season_year * z_season_year[y])` is the
hierarchical year-specific first-harmonic amplitude multiplier
(non-centred). `alpha_R`, `year_effect`, seasonal harmonics, generation
interval (`G=8`, `w` = discretised Gamma(4,2)), demographic accounting
(`N_end - N_start - births + deaths` reconciliation), and the seed
mechanism (`is_seed`/`X_seed`) are unchanged from v4.0-v4.8 throughout.

## 2. Canonical state fits currently accepted

| UF | q treatment | serology | canonical fit | HMC | identification status |
|---|---|---|---|---|---|
| CE | **fixed**, grid 0.05/0.10/0.15/0.20 (v4.6-v4.9 convention; q=0.05 established as "the Ceará reference" in the Bahia assessment) | single anchor (Juazeiro) | `28_v4_9_hierarchical_seasonality/outputs/q0.05/renewal_ceara_v4_9_fit_q0.05.rds` | **PASS** (div=0, max Rhat=1.0008 — cleanest of the 4; q=0.20 FAILS at 1.0125) | conditional on fixed q, not uniquely identified |
| BA | **estimated** (global logit_q) | 6-site, geographic-offset | `30_bahia_v4_9_replication/outputs_global_q_multisite_serology/on_full2/renewal_bahia_global_q_fit_on_full2.rds` | **PASS** | serology_anchored |
| RJ | **estimated** (global logit_q) | none (`J_sero=0`) | `34_rio_de_janeiro_v4_9_replication/outputs/caseonly/rj_global_q_case_only.rds` (ad098 rescue) | **PASS** (0 div after adapt_delta=0.98 rescue) | partially_identified_case_only (strong q-alpha_R ridge, \|cor\|~0.84) |
| MT | **estimated** (global logit_q) | none (`J_sero=0`) | `35_mato_grosso_v4_9_replication/outputs/caseonly/mt_global_q_case_only.rds` (ad098 rescue) | **PASS** (0 div after adapt_delta=0.98 rescue) | case_based_partial_identification (fixed-q profile: PARTIAL IDENTIFICATION, q=0.03-0.05 plateau) |
| PE | (multiple spatial rescues attempted, all failed) | — | — | — | **EXCLUDED** — absolute S not robustly identified (REGIONAL-R FEASIBILITY FAIL) |

All four candidate fits (CE q=0.05, BA on_full2, RJ ad098, MT ad098)
expose `S_prop`, `R0_t`, and weekly `cases`/`expected_reported_cases` in
their stanfit objects, and all share the identical 573-week window
(2015-01-04 to 2025-12-21) and `weekly_data` schema (`week_start,
week_of_year, t, year, cases, Tmean, PRCP, ...` plus demography columns).

**Recommended Phase 2 canary: CE @ q=0.05.** It is the cleanest HMC of any
candidate (max Rhat 1.0008), uses the ORIGINAL frozen v4.9 file (fixed q
removes one confound axis before introducing climate), and q=0.05 is
already the state's established reference scenario elsewhere in this
project — no new privileging decision is being made.

## 3. National state-week climate file

`03_Output/tables/national_wave_analysis/brazil_chik_uf_weekly_climate.csv`,
built by `09_build_brazil_uf_weekly_climate.R` (loops the existing,
unmodified `build_state_weekly()` over all 27 UFs against the national
ERA5-Land-derived DLNM panel `chik_dlnm_panel_muni_week_2015_2025.rds`).

- Columns: `state` (UF abbreviation), `week_start` (Date), `year`, `Tmean`
  (population-weighted mean daily temperature, deg C), `PRCP`
  (population-weighted weekly total precipitation, mm).
- Coverage: **27/27 states**, **575 weeks each** (2014-12-28 to
  2025-12-28), no gaps, no duplicate state-week rows.
- Missingness: **0** NA in `Tmean`, **0** NA in `PRCP`.
- `week_of_year` is not itself a column here but is defined identically
  everywhere upstream as `lubridate::isoweek(week_start)`
  (`02_Script/30_climate_covariates_dlnm/04_build_dlnm_panel.R`) — reused
  verbatim for the epidemiological-week climatology in Phase 1.

No date-alignment or climate-data problems found. Safe to proceed to
Phase 1.

## 4. Independent verification addendum (fresh session, before Phase 2)

Re-derived rather than assumed, per this repo's own audit norms:

- **Provenance of the canonical CE q=0.05 RDS**: `outputs/fit_q0.05_run.log`
  (mtime 15:01) records an EARLIER attempt at 600 total iterations that
  finished with `HMC gate: FAIL` (low bulk/tail ESS warnings).
  `outputs/fit_q0.05_rerun.log` (mtime 15:22:51) reran at 800 warmup + 400
  sampling = 1200 total iterations and finished `HMC gate: PASS`; the on-disk
  `q0.05/renewal_ceara_v4_9_fit_q0.05.rds` mtime matches the RERUN exactly
  (15:22:51.38), confirming the saved object is the passing rerun, not the
  earlier failing attempt. **Implication for Phase 2**: budget >=1200 total
  iterations/chain from the start for the climate canary; 600 was
  insufficient even for the frozen (non-climate) v4.9 structure at this q.
- **HMC gate re-run directly from the stanfit** (not just trusting the
  stored `bundle$hmc`): `rstan::get_num_divergent` = 0,
  `rstan::get_num_max_treedepth` = 0, BFMI (4 chains) = 0.887 / 0.947 /
  0.908 / 0.722 (all >= 0.3) -- matches the stored gate exactly.
- **Rhat caveat (new finding, not in the original audit table)**: the
  repo's `compute_hmc_gate()` uses `rstan::summary()$Rhat` (classic,
  non-rank-normalised) -> max Rhat = **1.0008** over the gate parameters,
  comfortably under the 1.01 threshold. Recomputing the SAME parameters
  with `posterior::summarise_draws()` (rank-normalised split-Rhat, the
  current Vehtari et al. 2021 default and what recent CmdStan reports by
  default) gives max Rhat = **1.0148** -- marginally OVER the repo's own
  1.01 threshold. Bulk/tail ESS are comfortably fine under either method
  (bulk ~527-570, tail ~581, vs the gate's 100 floor). This does not reverse
  CE q=0.05's status as the cleanest of the 4 candidates (BA/RJ/MT all
  needed an adapt_delta=0.98 rescue where CE passed at 0.95), but "PASS"
  should be read as "passes the repo's classic-Rhat gate," not as
  unambiguously clean by the modern statistic. **Recommendation**: Phase 3
  canary diagnostics should report BOTH Rhat statistics and gate on the
  stricter (rank-normalised) one.
- Only the CE q=0.05 canonical fit was re-verified at this depth (full
  reload + independent rstan/posterior recomputation) in this session; the
  BA/RJ/MT "PASS" statuses are taken from their own assessment
  docs/logs, not independently re-loaded here.
