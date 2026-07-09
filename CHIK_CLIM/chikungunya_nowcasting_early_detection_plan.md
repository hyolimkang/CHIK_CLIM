# KAKENHI Wakate Grant Plan (Mini-Proposal Draft)
## Source-to-recipient chikungunya outbreak intelligence: nowcasting source-country epidemic intensity and estimating Japan's introduction opportunity

**This document supersedes `grant_topic_assessment.md` and
`idea2_climate_immunity_importation_plan.md`**, which recorded earlier stages
of this idea (an initial three-way comparison, then a climate/immunity/
phylodynamics design that proved too broad, then an "early-warning alarm
system" framing that overreached given Japan has no historical local
chikungunya transmission to validate against). Keep those as history only.
This is the current plan.

## Rationale (why this framing, not the earlier ones)

Two corrections got us here:

1. **An operational "alarm system, validated by detecting past local
   outbreaks" cannot be built for Japan** — Japan has zero historical local
   chikungunya transmission and only one autochthonous dengue event (2014,
   plus a smaller 2019 recurrence). There is no ground truth to validate a
   detection algorithm against locally, and imported case counts (~0–9/year)
   are too sparse to fit any statistical model on directly — there is no
   continuous case stream to "nowcast" in the classical sense.
2. **But the premise that Japan has negligible outbreak-risk relevance is
   also wrong.** *Aedes albopictus* is confirmed established across
   essentially all major Japanese urban regions — Kanto, Kansai, Tokai, the
   Seto Inland Sea, and western Kyushu — per a 2024 Japan-specific
   distribution model
   ([Yang, Higa, Kasai et al., *PLOS ONE* 2024](https://journals.plos.org/plosone/article?id=10.1371%2Fjournal.pone.0303137)),
   and local Aedes-borne transmission has already happened twice: Tokyo 2014
   (162 autochthonous dengue cases,
   [PMC4344289](https://pmc.ncbi.nlm.nih.gov/articles/PMC4344289/)) and a
   smaller 2019 recurrence
   ([NIID IASR](https://www.niid.go.jp/niid/en/basic-science/865-iasr/9916-484te.html)).
   So the vector and the proof-of-principle both exist; what's missing is
   not "can it happen" but "when do the necessary conditions actually line
   up."

The resolution: reframe away from *predicting an outbreak* (overclaiming)
and away from *modelling Japan's own sparse case counts* (statistically
unsupported), toward **quantifying when source-country epidemic activity,
travel volume into Japan, and Japan's already-established seasonal vector
window overlap** — an "introduction opportunity" index. All three inputs are
either your own data (Brazil) or existing published data/models (travel
statistics, the 2024 vector distribution model) — nothing here requires
fitting a model to Japan's own sparse chikungunya case series.

This also has a direct, already-published precedent for dengue in Japan —
[Sakamoto, Yamauchi & Kokaze, "Mathematical model estimation of dengue fever
transmission risk from Southeast and South Asia into Japan between 2016 and
2018"](https://pmc.ncbi.nlm.nih.gov/articles/PMC10495242/) — which combined
source-country case counts, travel volume, and temperature-dependent vector
survival into exactly two endpoints: probability of importation, and
probability of autochthonous transmission given importation. Aim 2 below
adapts this design to chikungunya, using Aim 1's nowcasted source-country
intensity as an improved input in place of raw reported case counts.

## Grant title

**"Source-to-recipient chikungunya outbreak intelligence: nowcasting
source-country epidemic intensity and estimating introduction opportunity
in Japan"**

## Aim 1 — Nowcast source-country (Brazil) chikungunya epidemic intensity from line-list surveillance data

**Purpose**: produce a delay-corrected, real-time-accurate estimate of how
much chikungunya activity is actually occurring in a source country at any
given time — the input Aim 2 needs to assess Japan's importation pressure.
Also: quantify the minimum case volume needed for this method to work at
all, since that number is what justifies *not* attempting the same kind of
model on Japan's own sparse case series.

**Data**: `01_Data/chik_sinan_individual_2015_2024.rds` — individual-level
Brazilian SINAN records with `date_onset` and `date_notification` per case
(already in this repo; no new acquisition).

**Method**:
- Build the reporting triangle (state × onset-week × delay) —
  `02_Script/build_reporting_triangle.R`.
- Fit a Dirichlet-multinomial nowcasting model —
  `02_Script/stan/chik_delay_dirmult.stan`, `02_Script/fit_chik_delay_nowcast.R`.
- Retrospective validation: truncate the triangle at past cutoff dates,
  nowcast, compare to the eventual true count —
  `02_Script/validate_nowcast.R`.
- **New**: subsample the Brazilian data down to progressively lower annual
  case volumes and rerun validation, to find the minimum volume at which
  nowcasting remains usably calibrated.

**Output (feeds directly into Aim 2)**: a nowcasted, delay-corrected weekly/
monthly chikungunya incidence series for the source country/region, with
uncertainty — plus an outbreak-phase/growth indicator (e.g. rising, peak,
declining), reconstructed from onset dates rather than notification dates so
it reflects true epidemic timing, not reporting artefacts.

## Aim 2 — Estimate Japan's chikungunya introduction opportunity by combining source epidemic intensity, travel connectivity, and established seasonal vector suitability

**Purpose**: not to predict a Japanese outbreak, but to quantify **when**
overseas epidemic activity, traveler volume, and Japan's already-confirmed
Aedes-albopictus season actually overlap — the "introduction opportunity"
window — following the Sakamoto et al. dengue-in-Japan precedent above,
substituting Aim 1's nowcasted chikungunya intensity for their raw reported
dengue case counts.

**Data** (all existing/published — no raw Japanese case-level microdata
required for the core model):
- Aim 1's nowcasted source-country epidemic-intensity series.
- Monthly traveler volume by origin country/region into Japan — Japan
  Immigration Services Agency / JNTO statistics.
- Japan-specific *Aedes albopictus* seasonal suitability — the 2024 1km-
  resolution distribution model
  ([Yang, Higa, Kasai et al.](https://journals.plos.org/plosone/article?id=10.1371%2Fjournal.pone.0303137)),
  combined with monthly temperature (existing Open-Meteo fetch pattern) to
  get a temperature-dependent vector-activity season by region, following
  the survival-rate approach in Sakamoto et al.
- Aggregate, published Japanese import-case context only (NIID/IASR
  bulletins — counts and months, not individual-level data) — used to
  contextualize/sanity-check the model, not to fit it.

**Method**: combine three layers multiplicatively/probabilistically,
directly following the Sakamoto et al. structure —
1. **P(importation)**: source-country nowcasted intensity (Aim 1) × travel
   volume from that source into Japan, by month.
2. **P(local transmission given introduction)**: temperature-dependent
   *Ae. albopictus* activity/survival in each Japanese region, by month,
   using the 2024 distribution model's regional suitability as the spatial
   layer.
3. **Introduction opportunity index** = overlap of (1) and (2) in space and
   time — where and when do non-trivial importation pressure and an active
   vector season coincide?

**Output**:
- A **source-country prioritization** ranking — which source countries
  matter most for Japan's risk isn't necessarily the ones with the largest
  outbreaks; it's the product of outbreak size and travel volume, which can
  favor a moderate Southeast Asian outbreak over a larger, more distant one.
- A **seasonal risk-window** result — quantifying how much importation
  timing does or doesn't overlap with Japan's vector-active months
  (historically June–September, per the 2014/2019 dengue precedents).
- A **regional risk map** — Kanto/Kansai/Tokai/Seto Inland Sea/western
  Kyushu stratified by introduction-opportunity level, using the existing
  Japan-specific vector model rather than a newly built one.
- A **surveillance-trigger framework** — an operational threshold (e.g.
  nowcasted source intensity × travel volume × in-season vector suitability
  above a defined level → recommend enhanced chikungunya/dengue differential
  testing) rather than a retrospective detection claim.

**Explicit scope limits (stated in the proposal, not hidden)**: this is not
outbreak prediction. Three conditions must all align for local transmission
— an infectious imported case, sufficient concurrent vector activity, and
permissive temperature — and even where conditions look permissive,
established Japanese surveillance work has found sustained local chains
often still fail to establish
([Senda et al., *Emerging Infectious Diseases* 2018](https://pmc.ncbi.nlm.nih.gov/articles/PMC6106439/)),
suggesting repeated introduction/viral load may matter beyond simple
suitability. The deliverable is a risk window and prioritization framework,
not a predicted outbreak.

## Two-year workplan (brief)

- **Year 1**: Aim 1 in full (nowcasting model, validation, subsampling
  experiment) — existing pipeline, no new data acquisition.
- **Year 2 Q1–Q2**: assemble travel-volume and vector-suitability layers for
  Japan; build the P(importation) and P(local transmission) components.
- **Year 2 Q3**: combine into the introduction-opportunity index; sanity-
  check against aggregate published import-case timing and the 2014/2019
  precedents.
- **Year 2 Q4**: finalize source-country prioritization, seasonal-window,
  and regional risk-map deliverables; draft the surveillance-trigger
  framework; write final report. Korea/China treated as comparative
  extensions only if travel/vector data assembly time allows.

## Priority reading list (15–20, all verified via live search/fetch)

**Direct methodological templates**
1. Sakamoto, Yamauchi & Kokaze, "Mathematical model estimation of dengue fever transmission risk from Southeast and South Asia into Japan between 2016 and 2018" — [PMC10495242](https://pmc.ncbi.nlm.nih.gov/articles/PMC10495242/) — read first; the direct precedent for Aim 2's model structure.
2. Manica et al., "Risk assessment and perspectives of local transmission of chikungunya and dengue in Italy, a European forerunner," *Nat. Commun.* 2025 — [PMC12234714](https://pmc.ncbi.nlm.nih.gov/articles/PMC12234714/)
3. Senda et al., "Estimating Frequency of Probable Autochthonous Cases of Dengue, Japan," *Emerging Infectious Diseases* 2018 — [PMC6106439](https://pmc.ncbi.nlm.nih.gov/articles/PMC6106439/) — surveillance-trigger/aberration-detection methodology, Tokyo/Osaka.

**Vector suitability in Japan (reuse, don't rebuild)**
4. Yang, Higa, Kasai et al., "Tiger prowling: Distribution modelling for northward-expanding *Aedes albopictus* in Japan," *PLOS ONE* 2024 — [article](https://journals.plos.org/plosone/article?id=10.1371%2Fjournal.pone.0303137) — the core Japan-specific vector layer for Aim 2, including future SSP-scenario projections to 2090.
5. Kraemer et al., "The global distribution of the arbovirus vectors *Aedes aegypti* and *Ae. albopictus*," *eLife* 2015 — [article](https://elifesciences.org/articles/08347) — for Korea/China if extended.
6. Kraemer et al., "Past and future spread of the arbovirus vectors *Aedes aegypti* and *Ae. albopictus*," *Nat. Microbiol.* 2019 — [article](https://www.nature.com/articles/s41564-019-0376-y)

**Proof-of-principle: Japan's actual local Aedes-borne transmission history**
7. "Autochthonous Dengue Fever, Tokyo, Japan, 2014" — [PMC4344289](https://pmc.ncbi.nlm.nih.gov/articles/PMC4344289/)
8. NIID IASR, "Dengue fever and dengue hemorrhagic fever, 2015–2019" (covers the 2019 recurrence) — [article](https://www.niid.go.jp/niid/en/basic-science/865-iasr/9916-484te.html)
9. "Ongoing local transmission of dengue in Japan, August to September 2014," *WPSAR* — [article](https://ojs.wpro.who.int/ojs/index.php/wpsar/article/download/285/419?inline=1)

**Nowcasting methodology (Aim 1 core)**
10. Höhle & an der Heiden, "Nowcasting the Number of New Symptomatic Cases During Infectious Disease Outbreaks Using Constrained P-spline Smoothing," *Epidemiology* 2019 — [PMC6684223](https://pmc.ncbi.nlm.nih.gov/articles/PMC6684223/)
11. Bastos et al., "Nowcasting by Bayesian Smoothing," *PLOS Comp Biol* — [article](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007735)

**Climate-driven transmission mechanics (temperature layer for Aim 2)**
12. Mordecai, Ryan et al., thermal biology of mosquito-borne disease, *Ecology Letters* 2019 — [PMC6744319](https://pmc.ncbi.nlm.nih.gov/articles/PMC6744319/)
13. "Reviewing estimates of the basic reproduction number for dengue, Zika and chikungunya across global climate zones" — [ScienceDirect](https://www.sciencedirect.com/science/article/pii/S0013935120300050)

**Importation-risk modelling frameworks (general precedent for Aim 2's structure)**
14. Wu et al., "Seasonal and interannual risks of dengue introduction from South-East Asia into China, 2005–2015" — [PMC6248995](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6248995/)
15. "Dengue virus importation risks in Africa: a modelling study," *Lancet Planetary Health* — [PMC11649930](https://pmc.ncbi.nlm.nih.gov/articles/PMC11649930/)

**Regional imported-case surveillance context**
16. Japan imported chikungunya cases 2006–2016 (source-country breakdown) — [PubMed 29394382](https://pubmed.ncbi.nlm.nih.gov/29394382/)
17. Retrospective analysis, 16 imported chikungunya cases, Japan — [PMC5827309](https://pmc.ncbi.nlm.nih.gov/articles/PMC5827309/)
18. South Korea imported dengue cases 2020–2024 — [PubMed 41276240](https://pubmed.ncbi.nlm.nih.gov/41276240/)
19. NIID/IASR imported chikungunya/dengue surveillance overview — [IASR](https://idsc.niid.go.jp/iasr/32/376/tpc376.html)
20. "Multiple early local transmissions of chikungunya virus, Mainland France, from May 2025" — [PMC12355907](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12355907/) — current comparative context (non-endemic Europe, repeated introductions).

## Verification notes

All URLs above were returned by live web search/fetch during this session
and cross-checked directly (title, authors, methods, findings confirmed via
WebFetch, not assumed from search snippets) — none fabricated. File paths
referenced (`build_reporting_triangle.R`, `chik_delay_dirmult.stan`,
`fit_chik_delay_nowcast.R`, `validate_nowcast.R`) match the pipeline already
scoped in this repo's `analysis_plan.md`.
