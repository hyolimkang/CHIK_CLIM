# Number Needed to Vaccinate (NNV): Historical Age-12 Counterfactual, Ceará 2015–2021

**Status: mechanism-stress-test NNV, historical window only. Not a
policy-ready NNV estimate.** This document extends
`HISTORICAL_AGE12_VACCINE_COUNTERFACTUAL_ASSESSMENT.md` with a
dose-denominated efficiency analysis. It reuses the same paired Arm A/B/C
posterior simulations (`results/three_arm_results.rds`) — no refitting, no
new simulation, no change to any existing result.

**Revision note**: this version supersedes the first draft of the NNV
module on visualisation and terminology only (no numbers changed): the
NNV-over-time line plot was replaced by a discrete year-end cumulative bar
chart, the effect-decomposition violin plot was replaced by a point +
interval estimate for the main figure (violin kept as a supplementary
diagnostic), the age-specific dose/effect ratio was removed from all main
figures/tables (kept only as an explicitly-labelled, non-main diagnostic
file) and replaced by an age-specific *benefit* figure/table, and every
figure axis/label now uses comma-formatted numbers with no scientific
notation.

## Terminology used consistently in this document

- **"Number needed to vaccinate per infection averted"** and **"...per
  reported case averted"** are the two PRIMARY, PROGRAMME-LEVEL quantities.
  Both use total *population* infections/reported cases averted in the
  denominator (Arm A − Arm B), not a per-age-group denominator.
- For age groups, only **"infections averted"** and **"percentage reduction
  in infections"** are used. The phrases "NNV in age 0-11", "NNV in age
  65+", or "age-specific NNV" are **not used anywhere** in this analysis,
  because no vaccine doses are administered to those groups — see the
  direct/indirect interpretation section below.

## Numerator and denominator conventions (as specified)

- **Numerator**: cumulative vaccine **doses administered** (`doses_administered`,
  i.e. `coverage × entrants_to_target_total`) — every age-12 entrant offered
  the vaccine, whether or not they were susceptible or ultimately protected.
  Never "effectively immunised."
- **Primary denominator**: infections averted, Arm A − Arm B (total effect,
  full transmission feedback).
- **Secondary denominator**: reported cases averted = q × infections averted,
  under the accepted fit's fixed ascertainment q = 0.05 (this fit does not
  carry a posterior over q, so this is a **constant-q observation mapping**,
  not a draw-varying one — labelled "reported cases averted," never
  "symptomatic" or "clinical" cases, since no symptomatic-fraction model
  exists yet).
- All ratios computed **draw-by-draw first** (`NNV[d] = doses[d]/effect[d]`),
  then summarised — never `median(doses)/median(effect)` (validated
  explicitly, QA-NNV item 7).

## 1. How many vaccine doses are administered?

**913,032 doses** over 2015–2021 (identical to 1 decimal place across all
200 posterior draws — doses depend only on demography/age-cohort size and
the fixed coverage=1.0 assumption, neither of which varies across posterior
draws of the transmission model; only the *effect* denominator carries
posterior uncertainty).

## 2. How many infections are prevented?

**538,773 [95% CrI 473,374 – 567,351]** infections averted over the full
2015–2021 window (median, 200 posterior draws) — identical to Table 1/3 of
the main assessment, recomputed independently here for the NNV module.
Corresponding **reported cases averted (q=0.05, fixed): 26,939 [23,669 –
28,368]**.

## 3. Posterior NNV per infection averted

**NNV = 1.69 [95% CrI 1.61 – 1.93] doses per infection averted**, over the
full 7-year window. Interquartile range 1.66–1.75. This is the **primary,
programme-level NNV** for this stress-test scenario.

## 4. Posterior NNV per reported case averted

**NNV = 33.9 [95% CrI 32.2 – 38.6] doses per reported case averted**
(reported cases averted = q × infections averted, q = 0.05 fixed) — 1/q
times larger than the infection-based NNV, as expected from a constant
under-ascertainment factor.

**Probability of non-positive benefit: 0% for both metrics** — in all 200
draws, infections averted and reported cases averted were both strictly
positive (Figure NNV-1/QA-NNV item 5). This is expected under the
100%-coverage/100%-VE stress test in a population that is not already near
zero susceptibility; it would not necessarily hold under low coverage or
weak VE, which have not been evaluated here.

## 5. How does NNV evolve over the follow-up period? (Figures NNV-1, NNV-2)

**Substantially, and non-monotonically** — this is one of the most
important findings of the NNV module. Two complementary views are provided:
Figure NNV-1 shows cumulative NNV at *programme-start-relative* horizons
(1/2/3/5 years after 2015-01-04); Figure NNV-2 gives the calendar-year
picture in four panels, deliberately mixing annual (year-only) and
cumulative (since-programme-start) quantities so the reader can see both
the raw yearly pattern and the accumulating NNV it produces:

- **Panel A** — annual infections, no-vaccine vs vaccine (counts, not a
  ratio), by calendar year.
- **Panel B** — annual infections averted, `I^NV_y - I^V_y` (no-vaccine
  minus vaccine, that single year only), with a zero baseline so that
  negative (rebound) years are immediately visible by colour.
- **Panel C** — cumulative infections averted (entire population), since
  programme start, through each year-end.
- **Panel D** — cumulative NNV per infection averted: the posterior median
  of the **draw-level ratio** `cumulative_doses[d]/cumulative_infections_averted[d]`
  computed through each year-end, never a ratio of medians and never an
  annual (year-only) doses/effect ratio (panels A-B use raw annual counts
  and their difference only, exactly to avoid that).

| Year-end | Cumulative doses | Cumulative infections averted | NNV per infection averted |
|---|---|---|---|
| 2015 | 134,539 | 30 [19–45] | **unstable — not shown as a bar (see below)** |
| 2016 | 270,334 | 151,447 [133,339–168,414] | 1.79 [1.61–2.03] |
| 2017 | 406,370 | 730,941 [684,333–793,851] | **0.56 [0.51–0.59]** |
| 2018 | 537,003 | 706,189 [665,374–757,533] | 0.76 [0.71–0.81] |
| 2019 | 664,829 | 643,200 [620,236–662,307] | 1.03 [1.00–1.07] |
| 2020 | 790,071 | 579,102 [535,112–597,628] | 1.36 [1.32–1.48] |
| 2021 | 913,032 | 538,773 [473,374–567,351] | 1.69 [1.61–1.93] |

**2015 is flagged as unstable and deliberately not drawn as a bar.**
By the end of 2015 only ~30 infections had been averted against ~135,000
cumulative doses — every draw technically gives a finite ratio (so it is
not literally `NA`), but the ratio is not a meaningfully precise estimate
of programme efficiency this early, before the epidemic that the
programme's benefit depends on has even happened. Rather than plot a
seemingly precise but practically meaningless bar (implied 2015 NNV in the
thousands), Figure NNV-2 annotates 2015 with "Insufficient accrued benefit
for stable NNV" (Section 7/8 requirement) and shows only 2016–2021 as bars.
**This is "Option B" from the task specification** (mark the unstable early
period separately, plot the interpretable years on a linear scale) rather
than forcing a log axis over the whole series — chosen because 2016–2021
NNV values span less than a 4-fold range and are fully legible on a linear
axis, whereas a log axis spanning 2015's value would have compressed all
the *informative* years into an unreadably narrow band.

By the end of 2017 — immediately after Ceará's largest epidemic in this
window — cumulative infections averted (730,941) actually *exceeds* the
eventual 2021 total (538,773), giving the **lowest (most efficient)
year-end NNV of the whole follow-up period, ≈0.56** (Figure NNV-2, panel D).
Panels A-B make the underlying mechanism explicit: 2016 and 2017 are the
only two years where the vaccine arm has meaningfully fewer infections than
the no-vaccine arm (annual infections averted +151,419 and +580,216
respectively — panel B, blue). In **every year from 2018 to 2021, annual
infections averted is *negative*** (panel B, orange: −25,211 in 2018,
−60,824 in 2019, −62,303 in 2020, −42,126 in 2021 — median values, with
100% of the 200 draws negative in each of these four years, i.e. this is
not a chance sampling artefact). This is the post-epidemic transmission
rebound already documented in the main report: blunting the 2017 wave
leaves more natural (non-vaccine) susceptibles un-infected, who then
sustain modestly larger secondary waves in the vaccine arm than in the
no-vaccine arm during 2018–2021. Because cumulative infections averted
(panel C) is the running total of panel B's year-by-year values, it rises
sharply through 2017 and then **declines** every subsequent year as the
negative annual contributions accumulate, while cumulative doses keep
climbing linearly (Table `table_nnv_year_end_cumulative.csv`) as new
cohorts turn 12 every year regardless. The two trends combine to push
cumulative NNV (panel D) from 0.56 (2017) back up to 1.69 (2021). **A
single "the NNV of this programme" number is therefore highly sensitive to
when it is measured — this is a real property of the transmission dynamics
in this window, not noise**, and is the same phenomenon shown from the
start-relative-horizon angle in Figure NNV-1 (Section 9 below). The
underlying annual figures are saved in full in `tables/nnv/table_nnv_annual.csv`.

## 6. How much does indirect protection contribute to programme efficiency? (Figure NNV-4, Figure S-violin)

Using the existing three-arm decomposition (Total = A−B; Direct-only = A−C;
Indirect = C−B), and treating `doses / component` purely as a **mechanism
diagnostic ratio** (doses cannot literally be split between direct and
indirect recipients):

| Component | Infections averted (median) | Doses / effect (diagnostic, median) |
|---|---|---|
| Total (A−B) | 538,773 | 1.69 |
| Direct-only (A−C) | 82,295 | **11.1** |
| Transmission-mediated indirect (C−B) | 456,069 | **1.36** |

The direct-only diagnostic ratio (≈11 doses per direct-only infection
averted) is roughly 8× higher than the indirect diagnostic ratio (≈1.4).
Read together with the 84.7% indirect fraction already reported in the main
assessment, this shows the programme's overall efficiency (NNV≈1.69) is
overwhelmingly carried by transmission-mediated protection of people who
were never vaccinated — direct protection of the vaccinated cohort alone
would be a far less efficient use of the same doses in isolation.
**Restated once more: this decomposition is diagnostic, not an operational
claim that some doses "cause" direct effects and others "cause" indirect
effects — every dose contributes to both simultaneously.**

Figure NNV-4's main-figure panels show these three quantities as posterior
point estimates (median) with 50% (thick) and 95% (thin) credible
intervals, rather than full violin distributions — chosen because there are
only three categories being compared and an interval estimate is easier to
read at a glance for that purpose. The full posterior shapes are preserved
as a supplementary diagnostic
(`figures/supplementary/fig_s_nnv_posterior_distributions_violin`) for
readers who want to see e.g. skewness or multi-modality that a median/CrI
summary would not show.

## Direct/indirect interpretation, stated explicitly

The age-12 routine programme administers vaccine **only** to individuals
entering age 12 (the weekly age-11-to-12 ageing inflow). Nobody in any other
age group ever receives a dose in this scenario. Yet the **primary,
programme-level NNV** (Sections 3–4) divides total programme doses by
**total population** infections/reported cases averted — deliberately
including benefit realised in age groups that were never vaccinated. The
causal chain that justifies this is:

> vaccinate entrants to age 12 → fewer infections among those vaccinees →
> lower population infectiousness → lower future force of infection for
> **everyone** → additional infections averted in unvaccinated age groups
> (0–11, 13–17 carryover beyond direct protection, 18–64, 65+)

This is exactly what Section 6's three-arm decomposition quantifies (84.7%
of the total effect is this transmission-mediated component). **This must
not be read as "unvaccinated individuals have their own NNV."** No doses
were given to the 0–11 or 65+ groups, so there is no meaningful numerator to
divide their benefit by — Figure NNV-5 and Table `table_nnv_age_specific_benefit.csv`
therefore report only infections averted and percentage reduction for each
age group, never a group-specific dose ratio. The single, diagnostic-only
exception — "programme doses per infection averted *in* a given age group"
— is preserved solely as an internal sanity-check file
(`diagnostics/nnv/DIAGNOSTIC_doses_per_infection_by_group.csv`), explicitly
labelled as not-NNV, and excluded from every main table and figure.

## 7. How does NNV vary with the susceptible population at programme introduction? (Figure NNV-3)

**No usable within-window relationship could be estimated.** Baseline
susceptible fraction (S/N) at the first week of the window (2015-01-04) is
**exactly 1.0 for every one of the 200 posterior draws**, for the overall
population and for the age-12 band specifically — because this date
predates Ceará's first major chikungunya epidemic in this reconstruction
(2016–2017), the accepted model places the entire population as susceptible
at that starting point with no posterior uncertainty. With zero cross-draw
variance in the predictor, a Spearman correlation with NNV is undefined
(returned `NA`), and Figure NNV-3 shows all 200 draws collapsed onto a
single vertical line. **This is a genuine, disclosed limitation of using
this particular 7-year window to answer this particular question** — the
model may still show a susceptibility–NNV relationship, but the
2015–2021 window as defined here does not exercise it, since the programme
"introduction" in this counterfactual coincides with the immediate
pre-epidemic period rather than a moment of accumulated, variable immunity.
Re-examining this question would require either starting the counterfactual
programme after some accumulated transmission history (out of scope: the
task fixed 2015 as the start), or comparing susceptibility across the
2018–2021 sub-period across draws (not attempted here, since the primary
2015-introduction design was specified).

## 8. How many posterior draws show zero or negative estimated benefit?

**Zero of 200 draws** show infections_averted ≤ 0 or reported_cases_averted
≤ 0 in the main stress-test run (Pr(benefit > 0) = 100% for both). The
`NNV = NA` guard for non-positive benefit was still exercised and confirmed
correct via two dedicated structural unit tests (not from the main 200
draws, since none triggered it there): a coverage = 0 run (doses = 0,
impact = 0, NNV correctly `NA`, never 0) and a VE = 0 run (doses > 0 but
effective_protected = 0 and impact ≈ 0 to numerical precision, NNV
correctly `NA`, never an artificially small or falsely "successful" value).
See `diagnostics/nnv/NNV_VALIDATION.md` items 3–4.

## 9. How sensitive is NNV to the evaluation horizon?

**Extremely sensitive — see Section 5.** Median NNV per infection averted
ranges from ≈4,100 (1 year) to 0.56 (3 years) to 1.69 (full 7-year window)
within the very same programme and the very same posterior draws. Any
single-horizon NNV reported in isolation, without stating the evaluation
window, would be materially misleading for this transmission pattern
(one dominant epidemic followed by a partial rebound).

## 10. Are these NNV values suitable for policy interpretation yet?

**No, not yet, and for the same reason as the main assessment.** Every
number above is generated under a **mechanism stress test** — 100%
coverage, 100% infection-blocking VE, lifelong protection — deliberately
chosen as an upper bound to test the simulator's mechanics, not a forecast
of a real programme. Concretely:
- Real coverage will be well below 100%, which will proportionally shrink
  the numerator (doses) roughly linearly, but shrink the (mostly
  transmission-mediated) denominator **non-linearly** — herd protection
  from a partially-vaccinated cohort is not a fixed fraction of the
  100%-coverage effect. Realistic-coverage NNV cannot be inferred by simple
  rescaling of the figures in this report.
- Real VE against infection is very unlikely to be 100% or lifelong;
  waning is not modelled at all currently.
- The strong horizon-sensitivity documented in Section 5/9 means any
  future realistic-parameter NNV estimate must also state its evaluation
  horizon explicitly, and ideally should not be judged from a window that
  happens to end shortly after a large epidemic (this window's rebound
  dynamic is itself dependent on the arbitrary 2021 cutoff).

**The code is already fully parameterised** (coverage, VE_infection are
plain arguments to `run_arm_B`/`run_arm_C`, not hard-coded) — producing a
realistic-value NNV requires only choosing evidence-based coverage/VE inputs
in `00_config.R$SCENARIOS` and re-running scripts 06 → 08 → 11 → 12; no
change to the NNV calculation logic itself (11–13) is needed.

## QA summary

All 10 required NNV QA checks **PASS** — see
`diagnostics/nnv/NNV_VALIDATION.md` and `diagnostics/nnv/NNV_QA_summary.csv`
for full detail and the two structural unit tests (coverage=0, VE=0) run
specifically to exercise the undefined-NNV code path, which the main 200
posterior draws never trigger on their own.

## Files produced by this module (current, post-revision)

- Scripts: `scripts/11_calculate_nnv.R`, `scripts/12_plot_nnv.R`, `scripts/13_validate_nnv.R`
- Results: `results/nnv/nnv_posterior_draws.{rds,csv}`, `nnv_by_horizon.rds`, `nnv_year_end_cumulative_by_draw.rds`, `nnv_annual_by_draw.rds`, `nnv_cumulative_by_week.rds`, `nnv_age_specific_by_draw.rds`, `nnv_susceptibility_link.rds`
- Tables (main): `tables/nnv/table_nnv_primary.{csv,rds}`, `table_nnv_probability_benefit.csv`, `table_nnv_by_horizon.csv`, `table_nnv_year_end_cumulative.csv`, `table_nnv_annual.csv`, `table_nnv_by_week_summary.csv`, `table_nnv_age_specific_benefit.csv`, `table_nnv_susceptibility_link.csv`, `table_nnv_susceptibility_spearman.csv`
- Figures (main, `03_Output/figures/ceara_historical_age12_vaccine_counterfactual/nnv/`): `fig_nnv1_primary_nnv`, `fig_nnv2_year_end_cumulative_nnv`, `fig_nnv3_susceptibility_vs_nnv`, `fig_nnv4_effect_decomposition_interval`, `fig_nnv5_age_specific_benefit` (each `.pdf`/`.svg`/`.png`, ≥600dpi)
- Figures (supplementary, `.../supplementary/`): `fig_s_nnv_posterior_distributions_violin`
- Diagnostics: `diagnostics/nnv/NNV_VALIDATION.md`, `NNV_QA_summary.csv`, `DIAGNOSTIC_doses_per_infection_by_group.csv` (explicitly non-NNV, not a main output)

**Removed/superseded by this revision** (no longer present under these
names): `fig_nnv2_nnv_over_time.*` (continuous line plot), `fig_nnv4_direct_indirect_nnv.*`
(violin-only main figure), `fig_nnv5_age_specific_efficiency.*` (dose/effect
ratio figure), `tables/nnv/table_nnv_age_specific_efficiency.csv` (dose/effect
ratio table, main-tables location).

No file outside this analysis directory, and no existing file within it
prior to this module, was modified. No transmission or vaccine simulation
was rerun for this revision — only summaries, figures, tables, and
terminology changed, from the already-saved posterior draw-level results.
