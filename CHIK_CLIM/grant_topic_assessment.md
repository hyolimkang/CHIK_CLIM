# KAKENHI Wakate Grant Topic Assessment: Chikungunya Modelling Ideas

Comparative feasibility / policy-relevance / Japan-context assessment of three
candidate modelling topics for a 2-year KAKENHI Wakate (JSPS early-career
researcher grant — modest scope, expected to be doable by a PI + 1 RA).

## The three ideas

1. **Vaccine phase 4 / nowcasting**: chikungunya vaccine (live-attenuated)
   phase 4 vaccine-effectiveness (VE) trial design modelling using nowcasting
   simulation.
2. **Climate–immunity–importation**: how short-term climate variability and
   population immunity shape chikungunya epidemic timing/frequency/intensity,
   extended to importation risk modelling for Japan, Korea, China.
3. **Dengue vaccine phase 3 (Nagasaki)**: dengue vaccine (Nagasaki University
   candidate) phase 3 trial design modelling — site selection, sample size,
   feasibility.

## Comparison at a glance

| | Idea 1: Vaccine phase 4 nowcasting | Idea 2: Climate–immunity–importation | Idea 3: Dengue vaccine phase 3 (Nagasaki) |
|---|---|---|---|
| **Technical feasibility (PI+1RA, 2y)** | Feasible only if reframed (see below); self-contained simulation, no data partnership needed | Feasible; extends existing pipeline directly | Feasible only as a methods-adaptation exercise; needs external collaboration |
| **Existing repo infrastructure to reuse** | Age-structured vaccine-impact simulators (`postvacc_simulator_expanded.R`, `prevacc_simulator.R`, `sim_age_routine_v12f.R`) | Stan/Bayesian SIR pipeline, weekly precip auto-fetch, CMIP6 downscaling pipeline | None |
| **Policy relevance** | Real, but for outbreak-triggered emergency monitoring, not classic phase 4 | Real and current (NIID priority monitoring; sustained regional importation) | High, but candidate is still pre-Phase-II |
| **Japan-context relevance** | Weak/indirect (not endemic, no PMDA approval) | Strong and unusually timely (Foshan 2025, Yoyogi Park 2014) | Strongest institutional fit, but Nagasaki/NEKKEN-specific rather than national-policy-specific |
| **Novelty** | Moderate, if reframed | High — chikungunya-specific East Asia importation modelling is a real gap | Moderate — established methodology, new candidate |
| **Key risk** | "Phase 4 trial design" framing doesn't match how phase 4 studies are actually done | Two-part coherence: must be scoped as one project, not two | No dengue experience/collaborators; unsecured NEKKEN collaboration; candidate years from Phase 3 |

## Idea 1 — Chikungunya vaccine phase 4 / nowcasting VE trial design

**Technical feasibility.** Real phase 4/post-marketing VE studies for arboviral
vaccines use test-negative case-control or cohort designs — see the first
TAK-003/Qdenga post-marketing test-negative-design study (São Paulo, 2025)
([Lancet Infectious Diseases](https://www.thelancet.com/journals/laninf/article/PIIS1473-3099(25)00382-2/abstract)).
Nowcasting — correcting reporting-delay-truncated epidemic curves in real time
([Höhle & an der Heiden 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6684223/);
[Bastos et al. 2019](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007735))
— is an outbreak-situational-awareness tool. It is not currently used in the
VE/phase-4 literature; a search combining "nowcasting" with "vaccine
effectiveness" / "test-negative design" returned nothing. As stated, the idea
rests on a methodological mismatch between "phase 4 trial design" and what
nowcasting actually does.

A narrower, defensible reframe exists: *does real-time reporting-delay
correction improve rapid VE/safety monitoring during outbreak-triggered
chikungunya vaccination campaigns?* This is a self-contained simulation study
— synthetic epidemic + synthetic delayed reporting + synthetic vaccine
rollout, benchmarked against ground truth — that needs no manufacturer or
ministry data partnership, and can be built directly on this repo's existing
age-structured vaccine-impact simulators (`02_Script/postvacc_simulator_expanded.R`,
which already uses `ve_sus = 0.989`, the IXCHIQ trial efficacy estimate;
`02_Script/prevacc_simulator.R`; `02_Script/sim_age_routine_v12f.R`).

**Policy relevance.** Real and highly current, but for *emergency* monitoring
rather than classic phase 4. Both approved vaccines (Valneva's IXCHIQ,
Bavarian Nordic's VIMKUNYA) were deployed reactively during the 2024–25
La Réunion/Mayotte epidemic. IXCHIQ has had a turbulent regulatory history:
paused in adults ≥60y
([FDA, May 2025](https://www.fda.gov/vaccines-blood-biologics/safety-availability-biologics/fda-and-cdc-recommend-pause-use-ixchiq-chikungunya-vaccine-live-individuals-60-years-age-and-older)),
BLA suspended
([FDA, Aug 2025](https://www.fda.gov/safety/medical-product-safety-information/fda-update-safety-ixchiq-chikungunya-vaccine-live-fda-suspends-biologics-license-fda-safety)),
temporarily suspended in the UK for 65+
([MHRA](https://www.gov.uk/drug-safety-update/ixchiq-chikungunya-vaccine-temporary-suspension-in-people-aged-65-years-or-older)),
and voluntarily withdrawn (Jan 2026). VIMKUNYA was FDA-approved Feb 2025
([FDA approval letter](https://www.fda.gov/media/185478/download)). This is an
unsettled, high-attention regulatory space where rapid bias-corrected
monitoring is plausibly valuable.

**Japan-context relevance.** Weak/indirect. Chikungunya is not endemic in
Japan — cases are imported only
([JIHS fact page](https://id-info.jihs.go.jp/diseases/ta/chikungunya/index.html);
[NIID/IASR surveillance overview](https://idsc.niid.go.jp/iasr/32/376/tpc376.html)).
Neither IXCHIQ nor VIMKUNYA is PMDA-approved or confirmed under active review.
*Aedes albopictus* is established as far north as Tōhoku, so a vector-risk
link exists, but it is indirect.

## Idea 2 — Climate + immunity → epidemic timing, plus East Asia importation risk

**Technical feasibility.** Part 1 (climate/immunity → epidemic timing) is
well-precedented:

- Kikuti et al., localized CHIKV outbreak, Salvador — [PMC6396974](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6396974/)
- Spatiotemporal dynamics and recurrence of CHIKV in Brazil, *Lancet Microbe* — [full text](https://www.thelancet.com/journals/lanmic/article/PIIS2666-5247(23)00033-2/fulltext)
- Mordecai/Ryan et al., thermal biology of mosquito-borne disease, *Ecology Letters* 2019 — [PMC6744319](https://pmc.ncbi.nlm.nih.gov/articles/PMC6744319/)
- Rising dengue risk with ENSO amplitude, *Nat. Commun.* 2025 — [article](https://www.nature.com/articles/s41467-025-63655-0)

This extends the existing Stan/Bayesian SIR pipeline
(`02_Script/fit_chik_ceara_stan_weekly.R`) directly — incremental, not a new
methodological leap — especially given the climate-data infrastructure
already in this repo: a weekly precipitation auto-fetch from Open-Meteo
(`02_Script/fetch_ce_prcp_weekly.R`), a WorldClim/TerraClim spatial
FOI-covariate random-forest pipeline (`02_Script/current_covariates_baseline.R`),
and a CMIP6 future-climate delta-downscaling pipeline already built
(`02_Script/build_future_PRCP_pilot.R`, `02_Script/predict_future_PRCP_pilot.R`;
MIROC6, SSP245/SSP585). This substantially lowers the "steep learning curve in
climate" concern — the data engineering exists; what's missing is the
mechanistic climate+immunity link.

Part 2 (importation risk) is a distinct methodological domain — travel volume
× origin force-of-infection × destination vector suitability:

- Seasonal/interannual dengue introduction risk, South-East Asia → China — [PMC6248995](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6248995/)
- Dengue importation risk into Africa, *Lancet Planetary Health* — [PMC11649930](https://pmc.ncbi.nlm.nih.gov/articles/PMC11649930/)

A full mechanistic coupling of both parts is over-ambitious for 2 years with
PI+1RA, but is tractable if Year 2 uses a **lighter, largely independent**
importation-risk index (published vector-suitability layers + travel
statistics + Year-1-derived source-country FOI parameters) rather than a live
dynamical coupling. Framed this way it stays one coherent story: better
source-country epidemic models feeding better-parameterized importation risk.

**Policy relevance.** Genuine and current. NIID lists chikungunya among
priority imported-disease monitoring targets. Japan recorded 89 imported
chikungunya cases 2006–2016
([PubMed 29394382](https://pubmed.ncbi.nlm.nih.gov/29394382/)); South Korea
reported 551 imported dengue cases in 2020–2024 alone
([PubMed 41276240](https://pubmed.ncbi.nlm.nih.gov/41276240/)) — sustained
importation pressure region-wide with thin surveillance capacity for
translating source-country intelligence into local risk forecasts.

**Japan-context relevance.** Strong and unusually timely. Japan's 2014
autochthonous dengue outbreak (Yoyogi Park, 160 cases) proved local *Ae.
albopictus* transmission competence exists
([PMC4318974](https://pmc.ncbi.nlm.nih.gov/articles/PMC4318974/)). The
**July–August 2025 Foshan (Guangdong, China) chikungunya outbreak — over
10,000 cases**
([China CDC Weekly](https://weekly.chinacdc.cn/en/article/doi/10.46234/ccdcw2025.172))
is a directly relevant, very recent regional event demonstrating East Asia is
not immune to CHIKV establishment — a compelling, current hook.

## Idea 3 — Dengue vaccine phase 3 trial modelling (Nagasaki)

**Verified reality.** NEKKEN (Institute of Tropical Medicine, Nagasaki
University) does have a live-attenuated tetravalent dengue candidate,
**KD-382**, developed with KM Biologics under NEKKEN's Department of Tropical
Viral Vaccine Development (established 2023) with AMED funding
([NEKKEN TVVD](https://www.tm.nagasaki-u.ac.jp/nekken/en/departments/vaccine-development.html);
[AMED/DIDA announcement](https://dida.nagasaki-u.ac.jp/en/124/)). However
KD-382 has only completed **Phase I**
([npj Vaccines](https://www.nature.com/articles/s41541-025-01204-y);
[Phase I results, 2021](https://www.businesswire.com/news/home/20210324006033/en/KM-Biologics-Announces-Phase-I-Clinical-Study-Results-of-a-Live-Attenuated-Tetravalent-Dengue-Vaccine-KD-382)) —
no Phase II/III is scheduled or public. "Phase 3 trial modelling" is
therefore premature framing; a Phase 3 is plausibly still years away.

**Technical feasibility.** Site-selection/sample-size methodology for dengue
vaccine trials is a mature, well-published space:

- TAK-003 (Qdenga) 4.5-year efficacy — [Lancet Global Health](https://www.thelancet.com/journals/langlo/article/PIIS2214-109X(23)00522-3/fulltext)
- Butantan-DV, ~16,235 participants, NCT02406729 — [NEJM](https://www.nejm.org/doi/full/10.1056/NEJMoa2301790)
- Dengue force-of-infection mapping across countries — [PMC8730486](https://pmc.ncbi.nlm.nih.gov/articles/PMC8730486/)
- Serosurvey design for targeting vaccination — [PubMed 29944696](https://pubmed.ncbi.nlm.nih.gov/29944696/)
- WHO/TDR dengue vaccine trial guidelines — [PMC3056064](https://pmc.ncbi.nlm.nih.gov/articles/PMC3056064/)
- WHO position paper on dengue vaccines, May 2024 — [WHO](https://www.who.int/publications/i/item/who-wer-9918-203-224)

A Wakate project here would be *adapting* established methods to a new
candidate, not innovating methodologically — a defensible but modest scope.
The user has no dengue-specific modelling experience or collaborators, and no
dengue infrastructure exists in this repo, so credible site-selection
modelling would need either primary/secondary dengue seroprevalence data or an
actual working relationship with NEKKEN's TVVD/KD-382 team — an external
dependency not yet secured, and a real risk for a 2-year timeline.

**Policy relevance.** High once KD-382 nears Phase 2/3 — directly useful to
WHO/AMED/KM Biologics planning.

**Japan/Nagasaki-context relevance.** The strongest of the three —
direct institutional proximity to the actual vaccine developer, plus
NEKKEN's genuine overseas field stations (Vietnam, Kenya, Philippines) that
could plausibly host future trial-site work.

## Recommendation

**Idea 2** (climate + immunity → epidemic timing, scoped to feed a lighter
East Asia importation-risk index) is the strongest choice, on balance:

- **Best technical fit for PI+1RA** — more reusable existing infrastructure
  than either other idea: the Stan/Bayesian transmission pipeline, weekly
  precipitation auto-fetch, and an already-built CMIP6 delta-downscaling
  pipeline mean the "steep learning curve" concern is smaller than it first
  appears.
- **Strongest, most current Japan/policy hook** — the 2025 Foshan outbreak and
  Japan's own 2014 Yoyogi Park precedent give this proposal a compelling,
  timely narrative that neither other idea matches. Idea 1's Japan link is
  weak/indirect; Idea 3's is strong but hinges on an unsecured external
  collaboration and a candidate still pre-Phase-II.
- **Highest novelty** — most importation-risk literature is dengue-centric; a
  chikungunya-specific, climate-and-immunity-informed dynamic source model for
  East Asia importation risk is a genuine, defensible gap.
- **Main risk to manage: coherence.** This must be explicitly scoped as one
  project, not two (see workplan below).

Idea 1 remains viable as a smaller, self-contained methods paper (or a
future/parallel project) if honestly reframed away from "phase 4 trial
design" toward "outbreak-triggered rapid VE/safety monitoring" — it is the
most immediately buildable given existing simulator code, but its Japan
relevance is thin. Idea 3 has the best institutional/policy fit but is the
least feasible right now given no dengue expertise/infrastructure and an
unsecured collaboration dependency; worth revisiting once a TVVD collaboration
is secured and KD-382 nears Phase 2.

## Sketch 2-year workplan for Idea 2

**Year 1 — Climate + immunity → epidemic timing (source-country model)**
- Extend `02_Script/fit_chik_ceara_stan_weekly.R` (v12f) with short-term
  climate covariates (weekly precipitation already available via
  `01_Data/ce_prcp_weekly.rds` / `02_Script/fetch_ce_prcp_weekly.R`; extend to
  temperature) and a time-varying/immunity-dependent transmission term, fit
  across the existing 11-state UF panel (`02_Script/build_uf_weekly_panel.R`).
- Deliverable: posterior estimates of how short-term climate variability and
  population immunity jointly modulate epidemic timing, frequency, and
  intensity across Brazilian states — the paper/chapter that anchors the
  grant's Aim 1.

**Year 2 — East Asia importation-risk index**
- Build a literature/travel-data-based importation-risk index for Japan,
  Korea, China: origin-country force-of-infection/epidemic-timing parameters
  from Year 1 × travel volume × destination *Ae. albopictus* vector
  suitability (published layers) × destination climate suitability (reuse the
  existing CMIP6 pipeline in `02_Script/build_future_PRCP_pilot.R` /
  `02_Script/predict_future_PRCP_pilot.R` for future-scenario suitability, not
  a live dynamical coupling to Year 1's model).
- Deliverable: an importation-risk index/dashboard for the three countries,
  validated against known imported-case counts (Japan 2006–2016 series,
  Korea 2020–2024 dengue series as a validation analogue), explicitly framed
  as "informed by Year 1 source-country dynamics" rather than a fully coupled
  dynamical model.

## Verification notes

Every URL above was returned by live web research (WebSearch/WebFetch) during
this assessment — none were guessed or reconstructed from memory. File paths
referenced (e.g. `postvacc_simulator_expanded.R`, `fetch_ce_prcp_weekly.R`,
`build_future_PRCP_pilot.R`) were confirmed to exist in this repo at the time
of writing.
