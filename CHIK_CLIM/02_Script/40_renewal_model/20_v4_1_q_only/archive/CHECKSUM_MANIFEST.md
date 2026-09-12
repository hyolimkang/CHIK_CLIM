# Archived v4.0/td14 source -- checksum manifest

Archived on 2026-09-11 as the frozen reference for the v4.1 q-only ablation.
These are byte-for-byte copies (SHA-256 verified) of the exact files used
for the HMC-passing td14 run -- not reconstructed from memory.

| file | original path | SHA-256 |
|---|---|---|
| `renewal_ceara_v4_0_minimal_no_vaccine_ARCHIVED.stan` | `02_Script/stan/renewal_ceara_v4_0_minimal_no_vaccine.stan` | `c21c3d9e4e745bd30a12f735caab597f7b5a0d3ed8776ff98db12b54f99b75b5` |
| `01_fit_v4_0_minimal_no_vaccine_ARCHIVED.R` | `02_Script/40_renewal_model/17_v4_0_minimal_no_vaccine/01_fit_v4_0_minimal_no_vaccine.R` | `10d4226e03fa3159a700e3ee0670d51b5b04b9e7d2aab3ea6e2c26b9fb074216` |

## Verified: this is the run that PASSED the strict HMC gate

Source: `02_Script/stan/renewal_ceara_v4_0_minimal_no_vaccine_fit_td14.rds`,
`config$run_tag == "td14"`, `config$max_treedepth == 14`, `hmc$hmc_pass == TRUE`.

Stored HMC gate values (re-confirmed 2026-09-11):

- divergences = 0
- max_treedepth_hits = 0
- maximum Rhat = 1.00897 (recomputed independently via posterior::rhat on
  core parameters: 1.00566-1.00897 depending on exact parameter subset
  summarised -- both comfortably <=1.01)
- minimum bulk ESS = 472.7-480.8 (>= 100)
- BFMI range = [0.926, 1.083] (all >= 0.30)

Fixed inputs at time of that run: `q_fixed = 0.10`, `imports_per_week = 1`,
`year_effect_prior_sd = 0.40`.

## Exact values to preserve unchanged in v4.1 q-only

- `alpha_R ~ normal(log(1.2), 0.5)`
- `z_year ~ std_normal()`; `year_effect = year_effect_prior_sd * (z_year - mean(z_year))`
- `year_effect_prior_sd = 0.40` (fixed data, NOT estimated)
- `beta_sin ~ normal(0, 0.25)`
- `beta_cos ~ normal(0, 0.25)`
- `phi_obs ~ gamma(2, 0.1)`
- generation interval `w`: `diff(pgamma(0:8, shape=4, rate=2))`, renormalised, `G=8`
- `imports_per_week = 1` (fixed data, NOT estimated)
- Negative-binomial observation model: `C ~ neg_binomial_2(expected_reported_cases, phi_obs)`
- S/U demographic bookkeeping recursion (composition-based death/reconciliation split)
- Chain initialisation strategy: `v4_initial_values()` in the archived fit
  script (deterministic small offsets per chain: -0.06/-0.02/0.02/0.06)
