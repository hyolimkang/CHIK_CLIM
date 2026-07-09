# Reporting-Delay (Nowcasting) Pipeline — Analysis Plan

Plan for the task defined in [inst.md](inst.md).

## Context

`inst.md` specifies:

> Estimate the reporting-delay distribution for chikungunya cases using SINAN
> line-list data, within a Dirichlet-multinomial Bayesian framework.
> Input: weekly case data at the state level for 11 Brazilian states, from
> `readRDS("01_Data/chik_sinan_individual_2015_2024.rds")`.
> Include a validation method for nowcasting results.
> Constraint: account for the existing population conservation issue.

This is a new component alongside the existing SIR transmission-model pipeline
(`fit_chik_ceara_stan_weekly.R`). Its job is to characterize *how quickly chik
cases become visible* after they occur, and to validate that a nowcast built
from a partial (as-of-today) reporting picture tracks the eventual fully-reported
count. Two scope decisions were made up front:

- **Delay definition**: onset (`DT_SIN_PRI`) → notification (`DT_NOTIFIC`).
  Both fields already exist in the cleaned individual-level dataset
  (`clean_chik_sinan_brazil.R:274-275`), so no raw-data re-extraction is needed.
- **State selection**: top 11 UFs by cumulative confirmed cases 2015–2024,
  computed directly from data (reproducible), rather than a fixed hand list.

## Findings from the existing codebase

- **Individual-level data** (`01_Data/chik_sinan_individual_2015_2024.rds`,
  produced by `02_Script/clean_chik_sinan_brazil.R`) has: `date_onset`,
  `date_notification`, `uf_notif`, `uf_residence`, `muni_notif6`,
  `muni_residence6`, `is_confirmed_chik`, `is_lab_confirmed`. UF is a 2-digit
  IBGE code (`substr(muni6, 1, 2)`, same convention used in
  `build_uf_weekly_panel.R:74`).
- **Existing weekly aggregation**: `02_Script/build_uf_weekly_panel.R` builds
  one state's weekly panel at a time from the *municipality-week* panel
  (`chik_brazil_muni_week_2015_2024.rds`), which only stores final case counts
  — it has no delay stratification, so it cannot be reused as-is for the
  reporting triangle, but its state-selection/aggregation pattern
  (`substr(muni6,1,2) == uf_code`, ISO week indexing via `make_week_id()`) is
  the right template.
- **Existing Stan workflow convention**
  (`02_Script/fit_chik_ceara_stan_weekly.R`): numbered config → priors →
  data-list-per-model → `check_stan_panel()`-style dimension validation →
  `rstan`/`cmdstanr` fit → diagnostics. Stan files live in `02_Script/stan/`,
  matching fitted objects are cached as `.rds` next to the `.stan` file.
- **Known "population conservation" bug**
  (`02_Script/stan/ceara_weekly_v12f.stan:7-16`): the v12e SIR model added
  weekly births to S but never removed deaths, so S+I+R silently drifted from
  the observed `pop[t]`; v12f fixed it with explicit demographic in/out flow.
  This bug lives in the *transmission* model, not in case-count data — the
  delay/nowcast model estimated here doesn't consume population denominators,
  so "accounting for" this constraint means: don't reintroduce an analogous
  silent-drift bug, and if nowcasted (delay-corrected) counts are ever handed
  back to `fit_chik_ceara_stan_weekly.R` as fitting data, they must still be
  checked against `pop[t]` the same way `check_stan_panel()` does for the SIR
  models. This is a downstream-reconciliation note, flagged in the diagnostics
  script below, not a modeling requirement of the delay model itself.
- **Diagnostics pattern to reuse**: `02_Script/check_v12e_mcmc_diagnostics.R`
  (divergences, max-treedepth hits, E-BFMI, Rhat/ESS via `posterior`/`bayesplot`).
- **Functions/**: `fixNAs_adj.R` and `thinning_chik.R` are spatial-FOI-only
  helpers (buffer imputation, background-point thinning) — not applicable here.
- **Folder convention**: `01_Data/` (data), `02_Script/` (+ `stan/`, `Functions/`
  subdirs), `03_Output/` (`figures/`, `tables/`). No existing README/docs beyond
  task briefs like `inst.md`.

## Proposed pipeline

New files, following the existing naming/structure conventions:

1. **`02_Script/select_top11_states.R`**
   Rank UFs by cumulative `cases_confirmed` 2015–2024 (from
   `chik_brazil_muni_week_2015_2024.rds`, grouped by `substr(muni6,1,2)`), take
   the top 11, save the UF-code list to `01_Data/top11_uf_states.csv` so every
   downstream script reads a single source of truth instead of re-deriving it.

2. **`02_Script/build_reporting_triangle.R`**
   From the individual-level RDS: for each state × onset epi-week × delay
   (`notification week − onset week`, in whole weeks, capped at a max horizon
   `D_max` with an overflow bucket for delays beyond it), tabulate case counts
   into a reporting-triangle array `[state, week, delay]`. Save to
   `01_Data/chik_reporting_triangle_top11.rds`. Row-sums by state/week must
   equal the state's known weekly confirmed-case totals — assert this
   (mirrors the fail-fast style of `check_stan_panel()`).

3. **`02_Script/stan/chik_delay_dirmult.stan`**
   Dirichlet-multinomial likelihood: for each state (and optionally
   season/year stratum, to let the delay distribution drift over time), the
   vector of case counts across delay bins is Multinomial(N, p), with
   `p ~ Dirichlet(alpha)` — the standard reporting-triangle/nowcast structure
   (Höhle & an der Heiden; Bastos et al.). Partial pooling across the 11
   states via a shared hyperprior on `alpha`, consistent with the
   partial-pooling style already used for `s0_alpha/beta`, `rho_alpha/beta` in
   `UF_PRIORS`.

4. **`02_Script/fit_chik_delay_nowcast.R`**
   Same shape as `fit_chik_ceara_stan_weekly.R`: config block (states, year
   window, `D_max`), build the Stan data list with a `check_*`-style validator,
   fit via `rstan`/`cmdstanr`, save `chik_delay_dirmult.rds` next to the
   `.stan` file.

5. **`02_Script/check_delay_mcmc_diagnostics.R`**
   Same pattern as `check_v12e_mcmc_diagnostics.R` (divergences, treedepth,
   E-BFMI, Rhat/ESS) applied to the new fit — run this before trusting any
   nowcast output.

6. **`02_Script/validate_nowcast.R`** (the validation method inst.md asks for)
   Retrospective backtesting: for a set of historical "as-of" cutoff weeks,
   truncate the reporting triangle to only the delay bins observable as of
   that cutoff, use the fitted delay distribution to nowcast the eventual
   fully-reported count for the most recent weeks, and compare against the
   actual eventual count (now fully resolved with hindsight). Report, per
   state and per weeks-since-occurrence horizon: bias/MAPE of the nowcast vs.
   eventual total, and empirical coverage of the nowcast's 50%/95% intervals.
   Save results to `03_Output/tables/nowcast_backtest_top11.csv` and
   `03_Output/figures/nowcast_backtest_<uf>.png` (fan chart: nowcasted vs.
   eventual weekly cases, one panel per state).

## Verification

R/Stan analytical code has no UI, so verification is "run and inspect":

- Run `select_top11_states.R` and confirm the printed ranking looks sane
  (Ceará and Minas Gerais, the two states already fitted in
  `fit_chik_ceara_stan_weekly.R`, should appear near the top).
- Run `build_reporting_triangle.R` for 1–2 states first and confirm the
  row-sum-equals-known-totals assertion passes.
- Run `fit_chik_delay_nowcast.R` with a short MCMC run (small `iter`/`warmup`)
  on those same 1–2 states to confirm the Stan model compiles and samples
  without divergences before scaling to all 11.
- Run `check_delay_mcmc_diagnostics.R` and confirm Rhat ≈ 1, no divergences.
- Run `validate_nowcast.R` on the test states and sanity-check that backtest
  calibration numbers are non-degenerate (coverage roughly near nominal, not
  0% or 100%) before scaling to the full top-11 set.
