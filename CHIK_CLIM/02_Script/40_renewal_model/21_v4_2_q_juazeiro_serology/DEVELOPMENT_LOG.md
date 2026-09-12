# Ceará renewal-model development log -- v4.2 (q + Juazeiro serology)

## Lineage

- **v4.0 td14**: q fixed at 0.10, no serology. **HMC PASS.** Proven
  computational reference (`../17_v4_0_minimal_no_vaccine/`).
- **v4.1 q-only** (`../20_v4_1_q_only/`): q estimated (Beta(16.2,108.7)
  prior), everything else identical to td14. **HMC FAIL** (Stage B: 106/600
  post-warmup treedepth==14 hits). Root cause: all 4 chains agree with each
  other but converge to q~0.028 (far below the prior's 2.5th percentile of
  0.077), with alpha_R correspondingly inflated to ~0.69 (vs td14's ~0.146)
  -- a genuine prior-likelihood tension/ridge in the joint (q, alpha_R, ...)
  posterior, not a chain-separation artefact.
- **v4.2 (this experiment)**: v4.0/td14 backbone + estimated q (identical
  prior/setup to v4.1 q-only) + ONE additional likelihood term: the
  Juazeiro do Norte 2018 seroprevalence survey, via a beta-binomial
  observation model with FIXED kappa_sero=50 (not estimated). Tests whether
  one independent cumulative-infection observation is enough to break the
  q<->alpha_R ridge identified in v4.1 q-only.

## Serology source verification (Section 4)

Verified against the primary source (not reconstructed from memory):
**Barreto FKA et al., "Seroprevalence, spatial dispersion and factors
associated with flavivirus and chikungunya infection in a risk area: a
population-based seroprevalence study in Brazil," BMC Infectious Diseases
2020;20:881, DOI 10.1186/s12879-020-05611-5** (open-access full text:
PMC7685300).

Confirmed exactly:
- Location: Juazeiro do Norte, Ceará, Brazil.
- Collection period: **June-December 2018** (matches the task spec exactly).
- Design: cross-sectional, random spatial point sampling proportional to
  neighbourhood population (ArcGIS-based), all residents of selected houses
  invited -- population-based, not a convenience sample.
- Assay: **ELISA, IgM and/or IgG** (Euroimmun kits); CHIKV seropositive =
  IgM and/or IgG detected.
- **103/404 CHIKV positive (25.5%)** -- exact match.
- Age range 5-91 years (median 45) -- general population, not age-restricted.
- No mention of Quixadá or any other Ceará municipality in this study --
  independently confirms the earlier finding that "Quixada" (409/289, used
  elsewhere in this project's diagnose scripts) is NOT from this study and
  has no traceable source anywhere in this repository.

This is also independently corroborated within the repo: `CountryModel.xlsx`
(sheet `inclusion`), study_no=165, author "FKA, Barreto," 2018, sums to
N=404/N.pos=103 across its 4 age strata.

## Survey timing (Section 5)

No individual sample dates available in this repository -- used the
collection-window midpoint per the task's explicit fallback rule: (2018-06-01
+ 2018-12-31)/2 = **2018-09-15**, matched to the nearest weekly model index
-> **week_start = 2018-09-16 (index 194 of 573)**. Confirmed in Stage A log.

## kappa_sero = 50 calibration (Section 8)

Prior-predictive check (`01_kappa_calibration.R`, `kappa_sero_calibration.csv`):
if the true Ceará state prevalence were exactly 25% (matching Juazeiro's own
observed rate), a repeat 404-person Juazeiro-like sample would plausibly
range from 13.4% to 38.6% (95% predictive interval) under kappa_sero=50 --
a substantial representativeness-weakening relative to a naive
binomial(404, p), as intended, without being so loose as to carry no
information. Judged NOT obviously inconsistent with the intended
interpretation; used as specified without modification.

## Stan model identity (Section 1)

`renewal_ceara_v4_2_q_juazeiro_serology.stan` verified
(`00_verify_ablation_identity.R`, PASS) to be the archived, checksum-verified
v4.0/td14 source (SHA-256 `c21c3d9e...` matching
`../20_v4_1_q_only/archive/CHECKSUM_MANIFEST.md`'s record) plus EXACTLY: (a)
the same q_fixed->q change as v4.1 q-only, (b) one additive serology block
(new data: t_sero/sero_pos/sero_n/kappa_sero; new transformed parameters:
p_state_sero_at_anchor/p_sero_safe/alpha_sero/beta_sero; one new model-block
sampling statement; two new generated quantities). No other line changed.
Confirmed absent: sigma_year, B_year, z_year_free, any Quixada reference.

## Fit records

- **Stage A** (2026-09-11): 5+5, all 4 chains. Compile OK, 573 weeks
  confirmed, serology anchor correctly resolved to week index 194
  (2018-09-16), beta-binomial likelihood evaluated without numerical
  errors, sample_file CSVs grew correctly. NOT assessed for convergence
  (expected divergences from a 5-iteration run are meaningless).
- **Stage B**: IN PROGRESS / TBD.

## Decision-rule status

Pending Stage B result.
