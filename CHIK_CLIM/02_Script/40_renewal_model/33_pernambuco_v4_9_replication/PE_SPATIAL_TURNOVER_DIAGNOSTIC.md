# Pernambuco Spatial Turnover Diagnostic

**Status: descriptive-only diagnostic. No Stan model fitting, no changes to
q, priors, serology, or any frozen renewal-model component.** Mirrors
`31_bahia_spatial_turnover_diagnostic/` (state BA), adapted for Pernambuco
(UF 26, 184 municipalities), plus two PE-specific additions requested
directly: an explicit 2016-vs-2021/22 hotspot comparison, and a
cumulative-municipality-attack-under-candidate-q analysis (adapted from
`32_bahia_global_q_local_consistency_audit/`).

## Result: **MIXED** (same category as Bahia, with somewhat weaker
long-range persistence)

## 1. Municipality-week panel audit (Section 1)

184 municipalities x 573 weeks (2015-01-04 to 2025-12-21) = 105,432 rows,
0 NA cases/population, 0 municipalities with incomplete week coverage.
Summing the panel reproduces the exact PE state-level weekly series used by
the renewal model (`build_state_weekly(panel, "26")`) with **max abs diff =
0** across all 573 weeks. Total PE cases in this window: 58,692.

## 2. Wave definitions (Section 2)

Reused the existing objective national wave census
(`brazil_chik_wave_census_v2.csv`, `state == "PE"`,
`major_epidemic_primary == TRUE`) -- 9 major waves, `PE_wave_03` (onset
2015-11-29) through `PE_wave_11` (onset 2024-12-29). Recomputed wave totals
from the muni-panel differ from the census by at most 38 cases (out of
19,884, wave_07) -- a boundary-inclusivity artifact, not a data error;
recomputed totals used throughout for internal consistency (same
convention as Bahia).

| wave_id | onset | peak | end | total_cases (recomputed) |
|---|---|---|---|---|
| PE_wave_03 | 2015-11-29 | 2016-02-28 | 2017-12-17 | 15,849 |
| PE_wave_04 | 2017-12-31 | 2018-05-06 | 2018-12-16 | 655 |
| PE_wave_05 | 2018-12-23 | 2019-09-08 | 2019-12-15 | 1,360 |
| PE_wave_06 | 2020-05-03 | 2020-07-05 | 2020-12-20 | 3,005 |
| PE_wave_07 | 2021-02-07 | 2021-05-30 | 2021-12-19 | 19,846 |
| PE_wave_08 | 2022-01-16 | 2022-04-24 | 2022-12-18 | 12,994 |
| PE_wave_09 | 2022-12-25 | 2023-03-19 | 2023-06-25 | 777 |
| PE_wave_10 | 2023-07-09 | 2024-03-10 | 2024-12-22 | 2,209 |
| PE_wave_11 | 2024-12-29 | 2025-03-23 | 2025-12-28 | 1,204 |

## 3. Adjacent-wave turnover metrics (Section 4)

| Transition | Spearman | Top-10% Jaccard | Top-20% Jaccard | 80%-contributor Jaccard | Share from bottom-50% (prior) | Share from zero-case (prior) |
|---|---|---|---|---|---|---|
| wave_03 -> wave_04 | 0.267 | 0.086 | 0.194 | 0.250 | 0.133 | 0.002 |
| wave_04 -> wave_05 | 0.407 | 0.407 | 0.298 | 0.333 | 0.084 | 0.084 |
| wave_05 -> wave_06 | 0.391 | 0.152 | 0.345 | 0.200 | 0.043 | 0.043 |
| wave_06 -> wave_07 | 0.419 | 0.226 | 0.298 | 0.429 | 0.056 | 0.056 |
| **wave_07 -> wave_08** | **0.123** | **0.027** | **0.121** | **0.087** | **0.566** | **0.441** |
| wave_08 -> wave_09 | 0.370 | 0.310 | 0.298 | 0.083 | 0.143 | 0.024 |
| wave_09 -> wave_10 | 0.309 | 0.188 | 0.254 | 0.278 | 0.150 | 0.150 |
| wave_10 -> wave_11 | 0.458 | 0.267 | 0.298 | 0.471 | 0.051 | 0.017 |

Median adjacent Spearman = **0.380** (range 0.123-0.458); median top-10%
Jaccard = 0.207, top-20% = 0.298, 80%-contributor Jaccard = 0.264; median
share of a wave's cases from municipalities that were bottom-50%
incidence in the prior wave = 0.108, from zero-case municipalities = 0.049.

**The 2021->2022 transition (wave_07 -> wave_08) is a clear outlier**: only
0.123 Spearman correlation, and a striking **56.6%** of 2022's cases came
from municipalities that were in the bottom half of 2021's incidence
ranking, **44.1%** from municipalities with literally zero cases in 2021.
This is the sharpest single spatial-relocation event in the PE record --
directly relevant to the 2016-vs-2021/22 question below, since 2021 and
2022 (commonly treated together as "the 2021/22 epidemic") are themselves
spatially quite different from each other.

## 4. Full pairwise turnover matrix (Section 7)

Across all 36 wave pairs: Spearman ranges from **-0.001** (wave_06 vs
wave_08) to **0.458** (wave_10 vs wave_11), median 0.226; only one pair is
non-positive. Top-20% Jaccard ranges 0.088-0.345 (median 0.194). No pair
reaches strong persistence (>0.5); nearly every pair remains weakly
positive rather than crossing into pure turnover (~0). This matches the
same "weak-but-consistently-positive" pattern found for Bahia (all-pairwise
range there: 0.222-0.386), except PE has one essentially-null pair
(wave_06/2020 vs wave_08/2022) that Bahia did not.

## 5. Prior burden vs future wave (Section 5)

Pooled (all analysis waves) Spearman correlation between cumulative prior
incidence and subsequent-wave incidence = **+0.281** (per-wave range
0.123-0.418, all positive). Same qualitative direction as Bahia
(+0.352): municipalities with more PAST cumulative burden do not show less
subsequent incidence -- if anything, mildly more, consistent with
persistent structural risk (population density, vector ecology, prior
under-ascertained transmission correlating with future transmission) rather
than short-term depletion-driven anti-correlation. The correlation is
somewhat weaker for PE than Bahia, consistent with PE's greater long-range
turnover (Section 6 below).

## 6. Wave maps (Section 6)

`pe_wave_incidence_maps.png` (9-panel, shared colour scale) and
`pe_wave_max_incidence_summary_map.png` show two clearly dominant "peak
wave" clusters rather than one: **63 municipalities** had their maximum
incidence in the 2016 index epidemic (wave_03), and **62 municipalities**
had it in the 2022 epidemic (wave_08) -- together accounting for 125 of
184 municipalities (68%). The remaining waves each account for only 2-8
municipalities' peak. This bimodal "either 2016 or 2022, rarely anything
else" pattern is itself a form of MIXED behaviour: it is not spatially
diffuse turnover (many different waves each dominating different places),
but it is also not simple persistence (the SAME two waves recur as each
municipality's peak, but which of the two differs strongly by
municipality -- see Section 7's low 2016-vs-2021/22 correlation).

## 7. 2016 vs 2021/22 hotspot comparison (explicit user request)

Comparing the 2016 index epidemic (PE_wave_03) against the pooled 2021+2022
epidemics (PE_wave_07 + PE_wave_08):

- **Spearman correlation (municipality incidence), 2016 vs pooled
  2021/22 = 0.137** -- markedly weaker than the median adjacent-wave
  correlation (0.380), confirming that spatial pattern decays with time
  separation, not just wave-to-wave noise.
- **Jaccard overlap of top-20% hotspots = 0.156.**
- Hotspot-status transition (top-20% incidence in each era, 184 municipalities):

| Status | n municipalities |
|---|---|
| Persistent hotspot (both eras) | 10 |
| 2016-only hotspot (faded) | 27 |
| 2021/22-only hotspot (new/emerging) | 27 |
| Not top-20% in either era | 120 |

Only **10 of 184** municipalities (5.4%) were top-20%-incidence hotspots in
BOTH the 2016 epidemic and the pooled 2021/22 epidemics; 27 faded out and
27 new ones emerged -- a near-exact 1:1 churn of the hotspot set between
the two eras, with the persistent core being small relative to either
era's full hotspot list (10 of 37 in each era, ~27%). This is the clearest
single piece of evidence in this diagnostic that **PE's hotspot geography
is not fixed** across the decade, even though nearer-term (adjacent-wave)
correlations remain weakly positive.

Figures: `pe_2016_vs_2021_2022_incidence_maps.png` (3-panel, shared scale),
`pe_2016_vs_2021_2022_hotspot_transition_map.png` (4-colour status map),
`pe_2016_vs_2021_2022_scatter.png`.

## 8. Cumulative municipality attack fraction under candidate q

CRUDE, descriptive only: `A_i(2025) = cumulative_reported_incidence_i / q`,
evaluated at PE's fixed-q sweep grid (0.05-0.30) and at PE's estimated
global-q posterior (Model B, U14 serology: median q=0.0098, 95% CrI
[0.0075, 0.0151]).

| q scenario | median A_i | pop-weighted mean A_i | max A_i | % munis with A_i>0.5 | % munis with A_i>1.0 |
|---|---|---|---|---|---|
| global-q 2.5% (q=0.0151) | 0.251 | 0.824 | 11.27 | 29.9% | 16.3% |
| **global-q median (q=0.0098)** | **0.190** | **0.625** | **8.56** | **24.5%** | **11.4%** |
| global-q 97.5% (q=0.0075) | 0.124 | 0.407 | 5.56 | 15.8% | 7.1% |
| fixed q=0.05 | 0.037 | 0.123 | 1.68 | 5.4% | 1.1% |
| fixed q=0.10 | 0.019 | 0.061 | 0.84 | 1.1% | 0.0% |
| fixed q=0.15-0.30 | <0.013 | <0.041 | <0.56 | <0.6% | 0.0% |

**At PE's estimated global-q (q~0.01), 11-16% of municipality x checkpoint
combinations imply a biologically impossible cumulative attack fraction
(A_i > 1)** -- e.g. Cedro (A_i=8.6), Verdejante (7.3), Salgueiro (4.5),
Mirandiba (4.0), all interior Sertão-region municipalities far from
Recife. This is a substantially more severe version of the tension already
flagged in the Bahia consistency audit (there, q~0.14 produced NO
municipality with A_i>1). The mechanism is the same as in Bahia -- PE's
global q is estimated almost entirely from ONE hyperlocal, high-prevalence
serosurvey (U14, Recife, 37.2% observed) via `eta_geo`, and the resulting
very low state-wide q (needed to reconcile Recife's high seroprevalence
with the state's low case counts) is then applied uniformly to municipalities
whose case series look nothing like Recife's -- but the effect is far more
extreme in PE because the fitted q (~0.01) is roughly an order of magnitude
lower than Bahia's (~0.14). Fixed q in the 0.10-0.30 range does not produce
this problem (0-1% of municipalities exceed A_i=1). This is a local-
consistency warning about extrapolating one geographically offset serosurvey's
implied ascertainment to the whole state, not a new finding about
transmission dynamics -- flagged, not resolved, consistent with this
project's stop-rule discipline.

Figures: `pe_municipality_implied_attack_by_q_maps.png` (5-panel: 4 fixed-q
+ global-q median), `pe_municipality_A_gt1_flag_map.png`.

## 9. Classification: MIXED

Applying the same rule used for Bahia (all pairwise correlations
consistently positive but below ~0.4-0.5 -> MIXED; near-zero everywhere ->
TURNOVER; high everywhere -> PERSISTENCE):

- Adjacent-wave and full-pairwise Spearman correlations are consistently
  weak-to-moderate and (with one exception) positive -- not TURNOVER.
- No pair reaches strong persistence (>0.5); the explicit 2016-vs-2021/22
  comparison (0.137) and the wave_06-vs-wave_08 pair (-0.001) show that
  correlation decays toward zero as the time gap between compared waves
  grows -- not PERSISTENCE.
- The prior-burden-vs-future-wave correlation is positive (structural risk
  persists, not simple depletion), but hotspot IDENTITY still churns
  substantially over the long term (only 27% of each era's top-20%
  hotspots recur in the other era).

**PE shows the same qualitative MIXED signature as Bahia -- short/medium-term
spatial persistence (adjacent waves share structure) coexisting with
substantial long-term hotspot turnover (2016 vs 2021/22 hotspots overlap
only 15.6% by Jaccard) -- with PE's long-range turnover somewhat more
pronounced than Bahia's.**

## 10. Implication for the state-aggregated renewal model

Like Bahia, this MIXED spatial signature means PE's single state-wide
transmission/susceptibility trajectory is a population-weighted average
over a shifting mosaic of local epidemics, not a literal description of a
single homogeneously-mixing population. The renewal model's `eta_geo`
(fixed at `sero_geographic_sd=1.0`) already allows Recife's serology point
to sit far from the state-average trajectory, and the Section 8 finding
above shows why this offset is doing a lot of work: the state series alone
gives essentially no local information to compute a plausible weight for
Recife, so `eta_geo` is effectively unconstrained locally.

## 11. Output files

Tables (`03_Output/tables/pernambuco_v4_9_replication/spatial_turnover_diagnostic/`):
`pe_municipality_week_panel.rds`, `pe_municipality_week_panel_audit.csv`,
`pe_municipality_incomplete_weeks.csv`, `pe_wave_definitions.csv`,
`pe_municipality_wave_burden.csv`, `pe_wave_turnover_metrics.csv`,
`pe_prior_burden_future_wave.csv`, `pe_prior_burden_per_wave_spearman.csv`,
`pe_2016_vs_2021_22_hotspot_transition.csv`, `pe_2016_vs_2021_22_summary.csv`,
`pe_municipality_implied_attack_by_q.csv`,
`pe_municipality_implied_attack_summary.csv`,
`pe_municipality_flagged_A_gt1.csv`.

Figures (`03_Output/figures/pernambuco_v4_9_replication/spatial_turnover_diagnostic/`):
`pe_prior_burden_vs_future_incidence.png`, `pe_prior_burden_quintile_summary.png`,
`pe_wave_incidence_maps.png`, `pe_wave_max_incidence_summary_map.png`,
`pe_2016_vs_2021_2022_incidence_maps.png`,
`pe_2016_vs_2021_2022_hotspot_transition_map.png`,
`pe_2016_vs_2021_2022_scatter.png`,
`pe_wave_pairwise_correlation_heatmap.png`, `pe_top_municipality_overlap_heatmap.png`,
`pe_wave_case_contributor_turnover.png`,
`pe_municipality_implied_attack_by_q_maps.png`, `pe_municipality_A_gt1_flag_map.png`.

Scripts: `scripts/20_build_pe_municipality_week_panel.R` through
`scripts/27_pe_cumulative_municipality_attack_under_q.R`.

## 12. Stop rule (honoured)

Descriptive diagnostic only -- no Stan model touched, no q/prior/serology
changed, no wave definitions altered. Returned for review before any
further model or scope changes.
