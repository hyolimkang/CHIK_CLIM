# Idea 2 — Detailed Analysis Plan
## Climate- and immunity-driven chikungunya epidemic dynamics, and importation risk to Japan, Korea, and China

This is a detailed, standalone plan for the idea selected in
[grant_topic_assessment.md](grant_topic_assessment.md): using short-term
climate variability and population immunity to explain chikungunya epidemic
timing/frequency/intensity, then using that to assess importation risk to
Japan, Korea, and China.

## The core design problem, and how it's resolved

Japan, Korea, and (outside Guangdong) China have no endemic chikungunya
transmission, so there is no local force-of-infection (FOI) or seroprevalence
estimate to plug into an "immunity" term for these countries — and there
shouldn't be one. The project is **not** one model fit twice (once on Brazil,
once on East Asia). It is two components with different jobs:

- **Part 1 (Brazil, or wherever data exists)** estimates how short-term
  climate variability and *time-varying population immunity* interact to
  drive epidemic timing/size — this needs real case + serology data, which
  exists for Brazil via the existing pipeline.
- **Part 2 (Japan/Korea/China)** does not re-estimate immunity locally.
  Because these countries have no meaningful chikungunya transmission
  history, the population there is treated as fully susceptible
  (S/N ≈ 1) — a standard, defensible simplifying assumption, not a data gap.
  Part 2 asks a *different* question: given a case is imported, could local
  climate + vector conditions sustain an outbreak? This only needs Part 1's
  **transferable climate-response transmission mechanism** (e.g. a
  temperature/rainfall → transmission-rate relationship), evaluated on East
  Asian climate data — not East Asian case data. This is the same design
  used in existing importation-risk literature (dengue → China,
  [Wu et al.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6248995/); Africa
  dengue importation, [Lancet Planetary Health](https://pmc.ncbi.nlm.nih.gov/articles/PMC11649930/)),
  neither of which requires destination-country immunity data.

**A second, separate correction**: Brazil is very unlikely to be a real
source of Japan/Korea/China's imported cases. Japan's actual imported
chikungunya cases are dominated by Indonesia, India, the Philippines, and
Thailand
([PubMed 29394382](https://pubmed.ncbi.nlm.nih.gov/29394382/);
[PMC5827309](https://pmc.ncbi.nlm.nih.gov/articles/PMC5827309/)). So Brazil's
role is as the **calibration setting** for the transferable climate-response
mechanism (because that's where rich FOI/case data and the existing pipeline
exist) — it is *not* treated as a real source-country for East Asian imports.
The "force of importation" component instead uses actual surveillance-reported
incidence from the real dominant source region (Southeast/South Asia), as a
simple incidence proxy, not a full mechanistic re-fit.

## Model structure

**Part 1 — Climate + immunity → epidemic timing (Brazil, calibration setting)**

- Extend the existing Stan/Bayesian SIR pipeline
  (`02_Script/fit_chik_ceara_stan_weekly.R`, v12f) with:
  - Short-term climate covariates on the transmission rate β(t): weekly
    precipitation (already available,
    `01_Data/ce_prcp_weekly.rds` / `02_Script/fetch_ce_prcp_weekly.R`, and
    extendable to other UFs the same way) and temperature (new — Open-Meteo's
    archive API used in `fetch_ce_prcp_weekly.R` also serves temperature;
    same fetch pattern, new variable).
  - A biologically motivated functional form for the climate → β
    relationship (unimodal thermal-response curve, following
    [Mordecai/Ryan et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC6744319/)),
    fit as a shared (partially pooled) curve across the 11-state UF panel
    (`02_Script/build_uf_weekly_panel.R`), so it is a general
    climate-transmission relationship, not a single-state idiosyncrasy.
  - Time-varying susceptibility already partially in place via the
    births-in demographic flow in v12f
    (`02_Script/stan/ceara_weekly_v12f.stan`) — extend with the
    cumulative-attack/immunity-depletion logic already used for the
    `s0_alpha/beta`, `cum_attack_ext_mean` priors.
- **Deliverable 1a**: posterior estimates of the climate-response
  transmission function (the transferable object Part 2 needs).
- **Deliverable 1b**: state-level fitted/hindcast epidemic timing, frequency,
  and intensity across Brazil 2015–2024, explained by climate variability and
  immunity — the standalone paper/chapter anchoring Aim 1.

**Part 2 — Importation-and-establishment risk index (Japan, Korea, China)**

Risk(country, month) = **Force of importation** × **Local establishment
potential given introduction**

1. **Force of importation** = travel volume (origin → destination) ×
   surveillance-reported case incidence at origin (proxy for infectious
   traveler prevalence), summed over the real dominant source
   countries/regions for each destination (empirically identified per
   destination from its own import-case line-list, not assumed):
   - Japan: NIID/JIHS import-case records (already show Indonesia > India >
     Philippines > Thailand dominance —
     [IASR](https://idsc.niid.go.jp/iasr/32/376/tpc376.html)).
   - Korea: KDCA import-case records (need to check equivalent breakdown;
     likely similar Southeast/South Asia dominance, to be confirmed early in
     Year 2 rather than assumed).
   - China: China CDC import-case records; note China also has domestic risk
     via the 2025 Foshan (Guangdong) outbreak, treated separately below.
   - Origin incidence: national dengue/chikungunya surveillance bulletins
     (WHO/PAHO-equivalent regional bulletins, ProMED, or published outbreak
     time series) for Indonesia/India/Philippines/Thailand — used as a
     simple incidence-proxy input, **not** a full mechanistic Stan re-fit;
     this keeps Part 2 tractable for PI+1RA.
   - Travel volume: destination immigration/visitor-arrival statistics by
     nationality and month (Japan Immigration Services Agency / JNTO publish
     this; Korea Immigration Service similarly; China's equivalent is less
     transparent — flag as a data-access risk to resolve early, with
     UNWTO aggregate statistics as a fallback).

2. **Local establishment potential given introduction** = Part 1's fitted
   climate-response transmission function, evaluated using **local**
   destination climate data (extend the existing weekly climate-fetch
   pattern to Tokyo/Osaka, Seoul/Busan, Guangzhou/Foshan/Shanghai — same
   Open-Meteo-based approach already used for Ceará), combined with local
   *Aedes albopictus* presence/density (published distribution maps —
   *Ae. albopictus* is established as far north as Tōhoku in Japan) and the
   S/N ≈ 1 assumption discussed above. This produces a local R0(t)-style
   quantity, converted to a probability of a sustained outbreak given one
   introduction via a standard branching-process calculation (offspring
   distribution from local R0 and an overdispersion parameter k — this is
   the same logic used for assessing outbreak potential from a single
   imported case in other importation-risk work).
   - **Future scenarios**: reuse the existing CMIP6 delta-downscaling
     pipeline (`02_Script/build_future_PRCP_pilot.R`,
     `02_Script/predict_future_PRCP_pilot.R`; MIROC6, SSP245/SSP585) to
     project how local establishment potential in Japan/Korea/China changes
     under future climate — a natural, low-marginal-cost extension since the
     downscaling machinery already exists.

3. **Validation anchor — the 2025 Foshan (Guangdong) outbreak.** This is a
   rare, directly usable out-of-sample test: a known introduction context in
   a specific climate/vector setting, with a known outcome (>10,000 cases,
   [China CDC Weekly](https://weekly.chinacdc.cn/en/article/doi/10.46234/ccdcw2025.172)).
   Early in Year 2, check whether Guangdong CDC or a related group published
   a post-outbreak serosurvey — if one exists, it becomes a genuine
   calibration/validation data point for the establishment-potential
   sub-model (predicted vs. observed outbreak size given the known
   introduction and local conditions). If no such serosurvey exists, Foshan
   is still usable as a qualitative validation check (did the model predict
   elevated local establishment potential for Guangdong in July 2025?).

## Two-year workplan

**Year 1 (Q1–Q4) — Part 1: Brazil climate-immunity transmission model**
- Q1: Extend climate-data fetch to temperature; assemble multi-UF weekly
  climate covariates alongside the existing 11-state case panel.
- Q2–Q3: Extend `ceara_weekly_v12f.stan` with the shared climate-response
  transmission function and immunity-depletion dynamics; fit across the
  11-UF panel; run diagnostics (extend
  `02_Script/check_v12e_mcmc_diagnostics.R` pattern).
- Q4: Finalize Deliverable 1a (climate-response function) and 1b
  (state-level epidemic timing/intensity results); draft paper/chapter 1.

**Year 2 (Q1–Q4) — Part 2: East Asia importation-and-establishment risk**
- Q1: Confirm real dominant source countries per destination from
  NIID/KDCA/China CDC import-case records (don't assume — verify, as done
  above for Japan); resolve travel-volume data access (Japan/Korea
  straightforward; China needs an early fallback decision).
- Q2: Assemble origin incidence proxy series (Indonesia/India/Philippines/
  Thailand) and compute the force-of-importation component.
- Q3: Extend climate-fetch to East Asian cities; evaluate Part 1's
  climate-response function locally; assemble vector-suitability layers;
  compute the local-establishment-potential component; run the Foshan
  validation check.
- Q4: Combine into the risk index (historical + CMIP6 future scenarios);
  finalize Deliverable 2 (risk index/dashboard for Japan/Korea/China);
  draft paper/chapter 2; write final report.

## Risk register

| Risk | Mitigation |
|---|---|
| China travel-volume data less transparent | UNWTO aggregate statistics as fallback; scope China's force-of-importation component as lower-confidence/sensitivity-tested if needed |
| Real dominant source countries for Korea not yet confirmed | Verify from KDCA import-case records in Year 2 Q1 before building the incidence-proxy pipeline, not assumed from Japan's pattern |
| Climate-response function may not transfer cleanly outside the tropics | Treat local establishment potential as one input to a probability (via branching process), not a deterministic outbreak prediction; validate against Foshan |
| Two-part coherence (reviewer concern) | Frame explicitly as calibration-setting → transferable-mechanism → risk-index, as in this document, not as "Brazil model + separate Asia model" |
| No post-Foshan serosurvey available | Fall back to qualitative validation (did the model flag elevated risk pre-outbreak) |

## Reference list (idea 2 only, all verified via live search)

- Mordecai, Ryan et al., thermal biology of mosquito-borne disease — [PMC6744319](https://pmc.ncbi.nlm.nih.gov/articles/PMC6744319/)
- Kikuti et al., localized CHIKV outbreak, Salvador — [PMC6396974](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6396974/)
- Spatiotemporal dynamics and recurrence of CHIKV in Brazil, *Lancet Microbe* — [full text](https://www.thelancet.com/journals/lanmic/article/PIIS2666-5247(23)00033-2/fulltext)
- Rising dengue risk with ENSO amplitude, *Nat. Commun.* 2025 — [article](https://www.nature.com/articles/s41467-025-63655-0)
- Seasonal/interannual dengue introduction risk, SE Asia → China — [PMC6248995](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6248995/)
- Dengue importation risk into Africa, *Lancet Planetary Health* — [PMC11649930](https://pmc.ncbi.nlm.nih.gov/articles/PMC11649930/)
- Japan imported chikungunya cases 2006–2016 (source-country breakdown) — [PubMed 29394382](https://pubmed.ncbi.nlm.nih.gov/29394382/)
- Retrospective analysis, 16 imported chikungunya cases, Japan — [PMC5827309](https://pmc.ncbi.nlm.nih.gov/articles/PMC5827309/)
- Japan 2014 Yoyogi Park autochthonous dengue outbreak — [PMC4318974](https://pmc.ncbi.nlm.nih.gov/articles/PMC4318974/)
- 2025 Foshan (Guangdong) chikungunya outbreak, China CDC Weekly — [article](https://weekly.chinacdc.cn/en/article/doi/10.46234/ccdcw2025.172)
- South Korea imported dengue cases 2020–2024 — [PubMed 41276240](https://pubmed.ncbi.nlm.nih.gov/41276240/)
- NIID/IASR imported chikungunya/dengue surveillance overview — [IASR](https://idsc.niid.go.jp/iasr/32/376/tpc376.html)

## Verification notes

All URLs above were returned by live web search during this session — none
fabricated. File paths referenced (`fetch_ce_prcp_weekly.R`,
`build_uf_weekly_panel.R`, `ceara_weekly_v12f.stan`,
`check_v12e_mcmc_diagnostics.R`, `build_future_PRCP_pilot.R`,
`predict_future_PRCP_pilot.R`) were confirmed to exist in this repo.
