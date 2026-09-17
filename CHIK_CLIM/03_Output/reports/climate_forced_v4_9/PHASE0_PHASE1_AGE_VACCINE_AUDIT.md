# Age-Structured Vaccine-Counterfactual Extension — Phase 0 + Phase 1 Audit

## Phase 0 — Accepted baseline architecture

**Canonical fit**: `36_climate_forced_v4_9/outputs/ce_canary/ce_climate_forced_canary_q0.05.rds`
(CE, q=0.05 fixed, HMC PASS, accepted per `CLIMATE_V49_CANARY_ASSESSMENT.md`).
**Frozen (pre-climate) comparator**: `28_v4_9_hierarchical_seasonality/outputs/q0.05/renewal_ceara_v4_9_fit_q0.05.rds`.
Neither file is modified by this work.

Quantities available for forward simulation (all already present in the
stanfit's transformed parameters / generated quantities, or in
`weekly_data` / `stan_data`, per posterior draw):

| quantity | source | notes |
|---|---|---|
| `R0_t[t]` | transformed parameter | **already includes the fitted climate_multiplier(t)** — this IS `R0_historical[t]` for Phase 2 |
| `climate_multiplier[t]` | generated quantity | `exp(beta_T*z_T_anom[t] + beta_P*z_P_anom[t])`, for the Phase 7 factorial |
| `S[t]`, `U[t]`, `X[t]` | transformed parameters (absolute counts) | `S+U = N_start[t]` exactly, every t, every draw |
| `S_prop[t]`, `immune_prop[t]` | transformed parameters | `= S[t]/N`, `U[t]/N` |
| `N_start[t]`, `N_end[t]` | `stan_data` | population at week start/end |
| `births[t]`, `deaths[t]` | `stan_data` | weekly totals (SINASC-derived births; IBGE all-cause deaths, day-weighted) |
| `net_population_reconciliation` | `weekly_data` (also derivable: `N_end-N_start-births+deaths`) | closes the demographic identity exactly |
| `imports_per_week` | `stan_data`, fixed = 1 | constant, unchanged |
| `w[1..8]` (generation interval) | `stan_data` | Gamma(4,2) discretised, sums to 1 |
| `q` | `stan_data`, fixed = 0.05 | unchanged |
| `is_seed[t]`, `X_seed[t]`, `fit_index` | `stan_data` | 2022 conditioned seed window — `X[t]=X_seed[t]` when `is_seed[t]==1`, otherwise the renewal formula |
| `year_id[t]`, `year_effect[y]`, `A_year[y]`, `beta_sin1/cos1/sin2/cos2` | transformed parameters / `stan_data` | already baked into `R0_t[t]`; not needed separately for Phase 2 (R0_t is used directly) |
| `week_start` (573 weeks, 2015-01-04 to 2025-12-21) | `weekly_data` | the simulation clock |

**Key simplification for Phase 2**: because `R0_t[t]` already has the
fitted climate forcing folded in, the age-expanded simulator does **not**
need to reconstruct `log_R0` from its components — it reuses `R0_t[t]`
verbatim per draw, exactly as instructed ("R0_historical(t) is EXACTLY
the fitted climate-forced R0(t)").

## Phase 1 — Age-demographic data audit

**Found, already in the repository** (`01_Data/`, built by
`02_Script/00_data_prep/04_fetch_ibge_uf_age_population.R`):

| file | content |
|---|---|
| `ibge_pop_uf_single_age_expanded_2015_2024.rds` | single-year age 0-100 (100 = open-ended top group), by UF, by year |
| `ibge_pop_uf_age_group_2015_2024.rds` / `..._wide_...` | coarser standard IBGE age bands (same source, aggregated) |
| `ibge_pop_uf_age_group_split12_2015_2024.rds` | age bands with age 12 isolated (vaccine-eligibility-relevant cut, pre-existing) |

Answers to the six required questions:

1. **Available years**: 2015-2024 (single-year file). **Gap**: the v4.9
   fitting window extends to 2025-12-21; the age file does not cover 2025.
2. **Age resolution**: single-year, ages 0-100 (101 rows/state-year, age
   100 is an open-ended "100+" top group per standard IBGE convention).
3. **State coverage**: **27/27 UFs present**, including CE (`uf_code=23`,
   confirmed).
4. **Sex stratification**: **none** (`uf_code, uf_name, year, age,
   population` only). Not required by this task (no sex-specific
   mechanism requested).
5. **Age-specific mortality**: **not available anywhere in the repo.**
   Only aggregate UF-level all-cause deaths (`deaths_total`, annual, from
   `ibge_population_projection_uf_2024revision.rds`) exist — the same
   quantity already used, undifferentiated by age, in the fitted v4.9
   demographic accounting.
6. **Births available separately**: **yes** — SINASC individual birth
   records, already aggregated to weekly totals per state (e.g.
   `ce_sinasc_births_weekly.rds`), already used unmodified in the fitted
   v4.9 model. Not age-specific (not needed: births are a single
   age-0-entry event by definition).

**Also found** (context, not reused): an older, structurally unrelated
age-vaccine prototype (`02_Script/20_vaccine_impact/01-03_*.R`,
`02_Script/10_sir_transmission_model/05-07_*.R`) built against the
**superseded** `ceara_weekly_v12f` classic-SIR model (explicit `beta_t`,
`gamma` recovery rate, `S0_frac`, `rho` — a different parameterisation
entirely, predating the renewal-equation v4.x line). **Not used here** —
it targets a different, no-longer-frozen transmission model, and reusing
its mechanics would violate "do not refit an age-specific transmission
model" by importing an incompatible structure. Only its general
cohort-bookkeeping *idea* (age vector, ageing, routine-schedule coverage
lookup) is a useful precedent, not its code.

### Assessment: data ADEQUATE, with two disclosed minimal-choice gaps (not a stop condition)

Both gaps have a resolution that is **consistent with the task's own
explicit preference for parsimony** (Phase 3: "prioritize exact aggregate
consistency... rather than introduce complex age-specific mortality"):

- **No age-specific mortality rates** → deaths are allocated across ages
  **in proportion to each age's population share** (not a real
  age-specific rate). This is not a new epidemiological assumption, just
  a bookkeeping split of the SAME already-fitted `deaths[t]` total, and
  it makes `sum_a deaths[a,t] = deaths[t]` hold **exactly by
  construction** — satisfying Phase 3's stated requirement directly.
- **Age shares stop at 2024, fit window needs 2025** → the age
  *proportions* (not absolute counts) for 2024 are held constant and
  applied to 2025 weeks. Because these are proportions applied to the
  model's own already-fitted `N(t)` (not the IBGE file's own population
  totals, which differ slightly from the demographic-accounting
  population used in the v4.9 fit — see below), this introduces no
  inconsistency with the fitted aggregate.

**Reconciliation note**: the age-structure file's own population totals
(e.g. CE 2024 = 9,390,560) do not exactly match the population trajectory
actually used in the v4.9 fit (`ibge_population_projection_uf_2024revision.rds`,
CE 2024 = 9,233,656 — a ~1.7% difference, expected: different IBGE
release vintages/methodologies for the two products). **Resolution**: the
age file is used ONLY for its within-year age **shares**
(`population[a,year] / sum_a population[·,year]`, self-normalised to sum
to 1), never for its absolute counts. Shares are then multiplied by the
already-fitted `N(t)` from the v4.9 posterior draw. This guarantees
`sum_a pop_age[a,t] = N(t)` exactly, independent of the reconciliation
gap between the two IBGE products.

**Proceeding to Phase 2-4** on this basis.
