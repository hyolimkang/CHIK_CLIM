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

- **Part 1 (Brazil + Indonesia)** estimates how short-term climate
  variability and *time-varying population immunity/transmission intensity*
  interact to drive epidemic timing/size — this needs real case, serology, or
  genomic data, which exists for Brazil (case/serology, via the existing
  pipeline) and, if you secure the collaboration described below, for
  Indonesia (genomic/phylogenetic).
- **Part 2 (Japan/Korea/China)** does not re-estimate immunity locally.
  Because these countries have no meaningful chikungunya transmission
  history, the population there is treated as fully susceptible (S/N ≈ 1) —
  a standard, defensible simplifying assumption, not a data gap. Part 2 asks
  a *different* question: given a case is imported, could local climate +
  vector conditions sustain an outbreak? This only needs Part 1's
  **transferable climate-response transmission mechanism**, evaluated on East
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
role is as one **calibration setting** for the transferable climate-response
mechanism (because that's where rich FOI/case data and the existing pipeline
exist) — it is *not* treated as a real source-country for East Asian imports.
The real force-of-importation component needs a real source-country signal —
which is exactly what an Indonesian collaborator provides.

## Value of the Indonesian phylogenetics/surveillance collaborator

**Yes — this is a substantial, arguably transformative asset**, and directly
resolves the weakest link in the original design (how to get real,
non-Brazil source-country data without building a full mechanistic model for
every possible source country). A collaborator who has already generated
chikungunya sequence data from Indonesia *and* phylogenetically confirmed at
least one Indonesia→Japan introduction event turns Part 2 from
"plausible risk index built on proxies" into an empirically grounded,
partially-validated framework. Concretely, their data can be used four ways:

1. **An independent, second test of the climate-transmission mechanism (Part
   1b).** Time-stamped viral sequences from Indonesia allow phylodynamic
   reconstruction of transmission intensity over time there (effective
   reproduction number or relative genetic diversity via skyline/skygrid or
   birth-death-skyline methods — see reading list) — a signal that doesn't
   depend on Indonesia's (likely incomplete) case-reporting system. Testing
   whether this phylodynamic signal correlates with Indonesian climate
   variability the same way Brazil's case-based signal does is a genuine,
   publishable cross-validation of the climate-response mechanism in a
   second, independent, and — critically — *actually relevant* setting
   (tropical Indonesia, the real source region, not just Brazil).
2. **A direct, empirical measurement of "force of importation" (Part
   2a).** If the collaborator has both Indonesian and Japan-imported-case
   sequences, a phylogeographic model (discrete-trait diffusion, à la
   [Lemey et al. 2009](https://journals.plos.org/ploscompbiol/article?id=10.1371%2Fjournal.pcbi.1000520))
   can estimate the Indonesia→Japan viral export/migration rate over time —
   this literally *is* the force-of-importation quantity Part 2 needs for
   this corridor, measured directly rather than inferred from
   incidence×travel proxies. It also lets you validate the cheaper
   proxy-based method (which you'll still need for Korea/China and for
   source countries without genomic collaborators) against genomic ground
   truth for at least one corridor.
3. **Real calibration events for the establishment-probability sub-model
   (Part 2b/validation).** Each phylogenetically confirmed introduction has
   an inferred date (via molecular-clock tip-dating) and a *known outcome* —
   in Japan's historical cases, no local establishment occurred. These become
   real negative-control data points: "given this date/season, the model
   should predict low local establishment probability" — directly
   testable, and complementary to the 2025 Foshan outbreak, which is a
   positive control (introduction that *did* establish and spread).
4. **A concrete, named collaboration for the grant's feasibility section.**
   Reviewers commonly ask "how will you actually get import-side data?" —
   having a named collaborator with real Indonesia–Japan transmission-linkage
   data is a direct, credible answer, not a promise to build a proxy method
   and hope it works.

**Practical implication for scope**: keep the proxy-based approach
(incidence × travel volume) as the general method that extends to Korea,
China, and other source countries (India, Philippines, Thailand) where you
won't have equivalent genomic linkage data — but use the Indonesia
collaboration to *validate* that method for the one corridor where you can
check it against ground truth. This is more defensible than either relying
on proxies alone (no ground truth check) or trying to get phylogenetic
linkage data for every corridor (unrealistic for PI+1RA in 2 years).

## Step-by-step modelling flow: how Part 1 feeds Part 2

Read this top to bottom — each step names its inputs, what you do with them,
what comes out, and exactly which later step consumes that output.

---

**Step 0 — Data inventory (Q1, Year 1)**
- *Do*: confirm what already exists vs. must be newly acquired.
- *Have already*: Brazil 11-UF weekly case panel + FOI/serology-informed
  priors (existing repo); weekly precipitation for Ceará + fetch pattern
  extendable to other UFs; CMIP6 downscaling pipeline.
- *To acquire*: Indonesian sequence data + metadata (dates, locations) from
  the collaborator; Japan-imported-case sequences/dates (from the same
  collaborator or NIID); temperature data (extend existing Open-Meteo fetch);
  East Asian city-level climate data; vector-suitability layers; travel/visa
  statistics; origin-country incidence bulletins.

**Step 1 — Part 1a: Brazil mechanistic climate-immunity transmission model**
- *Input*: 11-UF weekly case panel, weekly precipitation + temperature,
  existing FOI-informed immunity priors.
- *Method*: extend `02_Script/fit_chik_ceara_stan_weekly.R` (v12f) with a
  shared climate-response transmission function β(temperature, rainfall) —
  thermal-response-curve form ([Mordecai/Ryan et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC6744319/)) —
  fit across all 11 states, plus the existing immunity-depletion/cumulative-attack
  logic (`s0_alpha/beta`, `cum_attack_ext_mean`).
- *Output A*: the fitted **climate-response transmission function** — the
  key transferable object.
- *Output B*: Brazil epidemic timing/frequency/intensity results —
  standalone Deliverable 1a/1b, doesn't feed downstream, anchors the grant's
  Aim 1 paper.
- *Feeds into*: Step 2 (as a comparison target) and Step 5 (evaluated on East
  Asian climate).

**Step 2 — Part 1b: Indonesia phylodynamic cross-validation (NEW, if
collaboration secured)**
- *Input*: Indonesian CHIKV sequences with sample dates/locations (from
  collaborator), Indonesian climate data.
- *Method*: phylodynamic reconstruction of transmission intensity over time
  (birth-death skyline or skygrid — see reading list) from the sequence data;
  regress/correlate the resulting Reff(t) or diversity signal against local
  climate variability.
- *Output*: an independent test of whether Step 1's climate-response *shape*
  generalizes to a real, non-Brazil, actually-relevant tropical setting. If
  yes → strong justification for transferring Step 1's function to Part 2.
  If no → still valuable: derive an Indonesia-specific climate-response
  function instead, and use *that* in Step 5 for the Indonesia-linked risk
  pathway.
- *Feeds into*: Step 5 (which climate-response function to transfer/evaluate
  on East Asian climate) and Step 3 (transmission-intensity time series used
  there).

**Step 3 — Bridge: phylogeographic export-rate estimation
(Indonesia→Japan)**
- *Input*: combined Indonesian + Japan-imported-case sequences with dates.
- *Method*: discrete-trait phylogeography
  ([Lemey et al. 2009](https://journals.plos.org/ploscompbiol/article?id=10.1371%2Fjournal.pcbi.1000520))
  or structured-coalescent approach to estimate Indonesia→Japan viral
  migration/export events and rates over time; cross-reference inferred
  introduction dates against Step 2's Indonesian transmission-intensity peaks
  and against known travel seasonality.
- *Output*: (i) a genomically-grounded force-of-importation series for the
  Indonesia–Japan corridor specifically; (ii) a set of dated, confirmed
  introduction events with known outcome (no local establishment) — negative
  controls for Step 6.
- *Feeds into*: Step 4 (as the validation benchmark for the general proxy
  method) and Step 6 (as calibration data).

**Step 4 — Part 2a: Force of importation, generalized to all corridors**
- *Input*: origin-country incidence bulletins (Indonesia, India, Philippines,
  Thailand — surveillance reports/ProMED, not full mechanistic fits) ×
  destination travel/visitor-arrival statistics (Japan Immigration/JNTO,
  Korea Immigration Service, UNWTO as China fallback).
- *Method*: simple incidence × travel-volume proxy, summed across each
  destination's real dominant source countries/regions (confirm Korea's and
  China's breakdown from KDCA/China CDC import-case records — don't assume
  it mirrors Japan's).
- *Validation*: compare the proxy-implied Indonesia→Japan import pattern
  against Step 3's genomic ground truth; if they track reasonably well, that
  supports using the cheaper proxy for corridors without genomic data
  (Korea, China, and non-Indonesia sources).
- *Output*: **Force-of-importation(country, month)** series for Japan, Korea,
  China.
- *Feeds into*: Step 7 (combined risk index).

**Step 5 — Part 2b: Local establishment potential given introduction**
- *Input*: Step 1's (and, if it differs, Step 2's) climate-response
  transmission function; local East Asian city-level climate data (extend
  the Open-Meteo fetch pattern to Tokyo/Osaka, Seoul/Busan,
  Guangzhou/Foshan/Shanghai); published *Aedes albopictus* suitability/density
  layers; the S/N ≈ 1 assumption.
- *Method*: evaluate the transferred climate-response function on local
  climate to get a local R0(t)-style quantity; convert to a probability of a
  sustained outbreak given one introduction via a branching-process
  calculation (offspring distribution from local R0 + overdispersion
  parameter k).
- *Output*: **Local-establishment-probability(country, month)** series,
  plus future-scenario versions using the existing CMIP6 pipeline
  (`build_future_PRCP_pilot.R`, `predict_future_PRCP_pilot.R`).
- *Feeds into*: Step 6 (validation) and Step 7 (combined risk index).

**Step 6 — Validation**
- *Input*: Step 3's confirmed Japan introduction events (negative controls —
  known non-establishment); the 2025 Foshan, Guangdong outbreak (positive
  control — known establishment and large outbreak size, plus a
  post-outbreak serosurvey if one is published).
- *Method*: check whether Step 5's establishment-probability function
  correctly predicts low probability for the historical Japan
  introduction dates/seasons, and elevated probability for Guangdong,
  July 2025.
- *Output*: a calibration check on the whole Part 2 pipeline before treating
  its outputs as more than illustrative.

**Step 7 — Part 2c: Combined risk index**
- *Input*: Step 4's force-of-importation series × Step 5's
  local-establishment-probability series.
- *Method*: Risk(country, month) = Force of importation × Local
  establishment probability, historical and under CMIP6 future scenarios.
- *Output*: **Deliverable 2** — an importation-and-establishment risk
  index/dashboard for Japan, Korea, and China, historical + future, validated
  in Step 6.

---

## Two-year workplan

**Year 1 — Part 1 (Brazil + Indonesia)**
- Q1: Data inventory (Step 0); confirm Indonesia collaboration scope and
  data-sharing terms; extend climate-data fetch to temperature.
- Q2: Fit Part 1a (Step 1) on the Brazil 11-UF panel; run diagnostics.
- Q3: Phylodynamic analysis of Indonesian sequences (Step 2); compare against
  Part 1a's climate-response function.
- Q4: Phylogeographic Indonesia→Japan export-rate estimation (Step 3);
  finalize Deliverables 1a/1b; draft paper/chapter 1 (Brazil + Indonesia
  climate-transmission comparison).

**Year 2 — Part 2 (Japan/Korea/China)**
- Q1: Confirm real dominant source countries for Korea/China (don't assume);
  resolve travel-volume data access.
- Q2: Build the force-of-importation series (Step 4), validated against
  Step 3's genomic ground truth for Indonesia–Japan.
- Q3: Extend climate fetch to East Asian cities; build the
  local-establishment-probability series (Step 5); run the validation checks
  (Step 6).
- Q4: Combine into the risk index (Step 7), historical + CMIP6 scenarios;
  finalize Deliverable 2; draft paper/chapter 2; write final report.

## Risk register

| Risk | Mitigation |
|---|---|
| Indonesia collaboration data-sharing/authorship terms not finalized | Resolve in Year 1 Q1, before phylodynamic work is scheduled to start |
| China travel-volume data less transparent | UNWTO aggregate statistics as fallback; scope China's force-of-importation component as lower-confidence/sensitivity-tested if needed |
| Real dominant source countries for Korea not yet confirmed | Verify from KDCA import-case records in Year 2 Q1 before building the incidence-proxy pipeline, not assumed from Japan's pattern |
| Phylodynamic methods (BEAST/skyline) are new to PI/RA | Budget dedicated training time in Year 1 Q1–Q2 (see reading list); consider a short course or a methods co-supervisor from the Indonesia collaborator's team |
| Climate-response function may not transfer cleanly outside the tropics | Treat local establishment potential as one input to a probability (via branching process), not a deterministic outbreak prediction; validate against Foshan and against Step 3's Japan negative controls |
| Two-part coherence (reviewer concern) | Frame explicitly as calibration-setting → transferable-mechanism → risk-index, as in this document, not as "Brazil model + separate Asia model" |
| No post-Foshan serosurvey available | Fall back to qualitative validation (did the model flag elevated risk pre-outbreak) |

## Priority reading list

Organized by what each cluster is *for*, roughly in the order you'd want to
read them.

**1. Ground the source-country mechanism (Part 1a/1b core)**
- Mordecai, Ryan et al., thermal biology of mosquito-borne disease — [PMC6744319](https://pmc.ncbi.nlm.nih.gov/articles/PMC6744319/)
- Kikuti et al., localized CHIKV outbreak, Salvador — [PMC6396974](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6396974/)
- Spatiotemporal dynamics and recurrence of CHIKV in Brazil, *Lancet Microbe* — [full text](https://www.thelancet.com/journals/lanmic/article/PIIS2666-5247(23)00033-2/fulltext)
- **Increased interregional virus exchange and nucleotide diversity outline the expansion of chikungunya virus in Brazil**, *Nat. Commun.* — [article](https://www.nature.com/articles/s41467-023-40099-y) — directly relevant: a CHIKV phylogeography study *in Brazil*, read this first as the bridge between your existing case-based Part 1a and the phylodynamic methods used in Part 1b.
- Reviewing R0 estimates for dengue/Zika/chikungunya across global climate zones — [ScienceDirect](https://www.sciencedirect.com/science/article/pii/S0013935120300050)

**2. Learn the phylodynamics/phylogeography methods (new territory — budget real time here)**
- Lemey, Rambaut, Drummond, Suchard, "Bayesian Phylogeography Finds Its Roots," *PLOS Comp Biol* 2009 — [article](https://journals.plos.org/ploscompbiol/article?id=10.1371%2Fjournal.pcbi.1000520) — foundational discrete-trait phylogeography method (implemented in BEAST); read this before attempting Step 3.
- Faria et al., "Establishment and cryptic transmission of Zika virus in Brazil and the Americas," *Nature* 2017 — [PMC5722632](https://pmc.ncbi.nlm.nih.gov/articles/PMC5722632/) — the best applied template for exactly this kind of arbovirus genomics + epidemiology + spread analysis, from the same research tradition.
- Nguyen et al., "Multiple chikungunya virus introductions in Lao PDR from 2014 to 2020" — [PMC9286254](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9286254/) — closest direct methodological analog: using phylogenetics to identify and date CHIKV introduction events into a single country, exactly what Step 3 needs.

**3. Understand the real source country (Indonesia)**
- "Chikungunya virus infection in Indonesia: a systematic review and evolutionary analysis," *BMC Infect Dis* — [PMC6417237](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6417237/)
- "Circulation of two chikungunya variants during an outbreak in Bali, Indonesia, 2021–2022" — [ScienceDirect](https://www.sciencedirect.com/science/article/pii/S1876034126000924) — recent genomic data directly from the collaborator's likely area of work.
- An investigation into CHIKV epidemiology across neglected regions of Indonesia, *PLOS NTD* — [PMC7785224](https://pmc.ncbi.nlm.nih.gov/articles/PMC7785224/)

**4. Importation-risk modelling frameworks (Part 2 core)**
- Seasonal/interannual dengue introduction risk, SE Asia → China — [PMC6248995](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6248995/)
- Dengue importation risk into Africa, *Lancet Planetary Health* — [PMC11649930](https://pmc.ncbi.nlm.nih.gov/articles/PMC11649930/)
- "A combination of climatic conditions determines major within-season dengue outbreaks in Guangdong Province, China" — [PMC6341621](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6341621/) — directly relevant to the Foshan validation step (Step 6); same province.

**5. Japan/Korea/China surveillance context**
- Japan imported chikungunya cases 2006–2016 (source-country breakdown) — [PubMed 29394382](https://pubmed.ncbi.nlm.nih.gov/29394382/)
- Retrospective analysis, 16 imported chikungunya cases, Japan — [PMC5827309](https://pmc.ncbi.nlm.nih.gov/articles/PMC5827309/)
- Japan 2014 Yoyogi Park autochthonous dengue outbreak — [PMC4318974](https://pmc.ncbi.nlm.nih.gov/articles/PMC4318974/)
- 2025 Foshan (Guangdong) chikungunya outbreak, China CDC Weekly — [article](https://weekly.chinacdc.cn/en/article/doi/10.46234/ccdcw2025.172)
- South Korea imported dengue cases 2020–2024 — [PubMed 41276240](https://pubmed.ncbi.nlm.nih.gov/41276240/)
- NIID/IASR imported chikungunya/dengue surveillance overview — [IASR](https://idsc.niid.go.jp/iasr/32/376/tpc376.html)
- Rising dengue risk with ENSO amplitude, *Nat. Commun.* 2025 — [article](https://www.nature.com/articles/s41467-025-63655-0)

## Verification notes

All URLs above were returned by live web search during this session — none
fabricated. File paths referenced (`fetch_ce_prcp_weekly.R`,
`build_uf_weekly_panel.R`, `ceara_weekly_v12f.stan`,
`check_v12e_mcmc_diagnostics.R`, `build_future_PRCP_pilot.R`,
`predict_future_PRCP_pilot.R`) were confirmed to exist in this repo. The
specific paper describing the collaborator's own Indonesia→Japan
phylogenetic linkage finding was not independently located via public search
(it may be unpublished or in preparation) — send me the citation and I'll
fold it into this document and the reading list directly.
