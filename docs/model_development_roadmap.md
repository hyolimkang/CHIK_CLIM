# CHIK_CLIM: Model Development Roadmap and Decision Gates

**Updated:** 2026-09-07
**Scope:** Retrospective chikungunya susceptibility reconstruction in Ceara, Brazil, followed by conditional climate analysis and scenario projection.
**Immediate priority:** Correct population accounting and establish computationally credible full-period inference before adding climate.

> Status: this is a proposed development and evaluation plan, not a record that all listed models have been implemented or validated. The version history below is based on project updates and figures supplied in discussion. The current GitHub source could not be independently retrieved during this review. In particular, the exact v1 implementation must be checked against its archived code. Proposed version names are suggestions.

## 1. Current decision

The population-input problem has been identified: the old state series combined municipal estimates from different revisions with a 2022 census count and subsequent estimates. Treating those discontinuities as real demographic movements created artificial changes in susceptible and non-susceptible counts.

The immediate next comparison is:

```text
v2.1 full-period stress test, mixed population vintages
    -> v2.1 full-period population-data correction, otherwise unchanged
    -> explicit demographic accounting
    -> targeted improvements to identification and computation
    -> validated climate-free full-history model
```

Do not simultaneously change population inputs, births/deaths, reporting, temporal priors, and importation. A successful new fit should have an explainable source of improvement.

### Population correction reported by the project

These values were supplied by the project; they have not been independently re-extracted from the IBGE workbook in this review.

| Year | Previous population input | Replacement: 2024 revision |
|---|---:|---:|
| 2021 | 9,240,580 | 9,129,520 |
| 2022 | 8,794,957 | 9,162,955 |
| 2023 | 9,014,301 | 9,196,672 |
| 2024 | 9,233,656 | 9,233,656 |
| 2025 | 9,268,836 | 9,268,836 |

The revision, workbook, worksheet, reference date, geography and extraction procedure must be recorded. Use the IBGE 2024 revision archive as the source registry.[^ibge]

## 2. Version history: what has actually been learned

| Version | Purpose and reported changes | Interpretation and current status |
|---|---|---|
| **v1** | Initial renewal reconstruction prototype. | Archive and inspect the original code before describing its exact likelihood, depletion assumptions or priors. Do not infer its specification from its name. |
| **v2** | Short-period reconstruction using cases and a Juazeiro serology anchor; fixed population; fixed symptomatic probability; reporting component; flexible time-varying transmission structure. The later project description indicates a site-to-state offset was already present with a fixed geographic scale. | Established a useful reconstruction workflow. It should not be described as literally equating Juazeiro with Ceara if that offset was present. One reported run had max R-hat around 1.194 and minimum ESS around 16.5, so it was not a validated posterior benchmark. |
| **v2.1, approximately 2015-2019** | Two serosurveys with geographic offsets, estimated geographic scale, uncertain symptomatic probability, and time-varying population with a net-change adjustment. | Cases could still be fitted while latent infection burden and susceptibility changed. This demonstrates model sensitivity, not proof of a more accurate statewide immunity estimate. Reported diagnostics improved but some R-hat values remained above 1.01. |
| **v2.1 full-period stress test, 2015-2025** | Extended the fitting horizon without a planned redesign of the scientific model. | Reported max R-hat 1.663, min ESS about 9.5, median ESS about 61.3, and min E-BFMI 0.279. The combined-chain state summaries are not suitable for scientific inference. The population-vintage artifact was identified during this audit. |
| **v2.1 full-period population correction** | Proposed immediate checkpoint: replace mixed-vintage population input while preserving the same full-period likelihood, priors and demographic approximation. | Not yet an explicit birth/death model. Its purpose is to separate a data correction from a structural change. |
| **v3 demographic reconstruction** | Proposed next scientific development: explicit, reconciled demographic flows and unambiguous state definitions. | Proceed through the gates below. Do not label the model validated merely because it compiles or fits case peaks. |

A changed fitting horizon changes the conditioning data. If a time covariate is rescaled over the fitting window, it can also change the meaning of a prior. Record both effects rather than calling every date extension an otherwise identical model.

## 3. The central distinction: population size versus population composition

### 3.1 State definitions

Use one explicit timing convention:

- `N[t]`: total resident population at the start of week `t`.
- `S[t]`: susceptible residents at the start of week `t`.
- `U[t] = N[t] - S[t]`: residents no longer susceptible under the model's infection-derived protection assumption. This includes recently infected residents; it is not an explicit recovered compartment.
- `X[t]`: new infections occurring during week `t`.
- `B[t]`, `D[t]`: births and all-cause deaths during that week.
- `M_balance[t]`: signed net population reconciliation flow, defined below.

The model tracks infection incidence, not a separate infectious-person stock. Avoid calling `U` an infectious or recovered compartment unless those states are actually implemented.

**Infection moves people from S to U; it does not, by itself, remove them from N.** Disease deaths, if later represented explicitly, must not also be subtracted a second time through all-cause mortality.

### 3.2 Conservation identities

The demographic balance is

$$
N_{t+1}=N_t+B_t-D_t+M_t.
$$

A consistent state accounting is

$$
S_{t+1}=S_t-X_t+B_t-D_t^S+M_t^S,
$$

$$
U_{t+1}=U_t+X_t-D_t^U+M_t^U,
$$

with

$$
D_t^S+D_t^U=D_t,\qquad M_t^S+M_t^U=M_t.
$$

Summing the state equations recovers the population equation. The assumption that births enter S is an aggregate approximation; neonatal protection and age-dependent risks are not explicitly represented at this stage.

These equations are accounting requirements, not a demand to sample three independent population trajectories. `N` can remain an externally specified input and `U` can be derived.

### 3.3 Why net population change is insufficient

Consider a hypothetical population of 1,000,000, initially 40% susceptible, with no new infections. There are 10,000 births and 10,000 deaths. If deaths are allocated proportionally to the initial S/U composition, susceptible deaths are 4,000 and non-susceptible deaths are 6,000.

Then total population is unchanged, but:

```text
S: 400,000 - 4,000 + 10,000 = 406,000
U: 600,000 - 6,000          = 594,000
N:                           1,000,000
```

Thus a model driven only by `delta_N = 0` misses demographic turnover. Correcting the population vintage removes artificial jumps; it does not add this turnover mechanism.

## 4. Demographic input preparation

### 4.1 Population stocks

Use one revision for the complete retrospective series. Preserve the published population reference date rather than automatically assigning annual values to January 1.

For weekly interpolation:

1. Read the reference-date metadata; use July 1 only when that is the documented date for the selected table.
2. Interpolate between consistently defined annual stock estimates onto the model's weekly boundary dates.
3. Obtain boundary-year data needed to cover the first and final weekly states; do not silently extrapolate from incomplete annual coverage.
4. Keep the case pipeline's epidemiological-week convention unchanged.
5. Label the resulting values as interpolated population estimates, not weekly census observations.

Population interpolation must follow source harmonization. Smoothing a mixed-vintage discontinuity is not source harmonization.

### 4.2 Birth and death flows

Prefer a documented, internally coherent combination of population stocks and demographic components. Check the IBGE revision's component tables and their interval definitions first. An alternative is residence-based SINASC births and SIM deaths, reconciled with the population estimates.[^ibge][^sinasc][^sim]

Do not treat these alternatives as interchangeable:

- Projection-system demographic components are estimated quantities within that system.
- Registered births and deaths are surveillance/vital-registration quantities with their own coverage, revisions and reference periods.

As checked on 2026-09-07, the official SINASC 2025 CSV resource is labelled preliminary. Record the actual downloaded release and any later updates rather than treating a provisional year as finalized.[^sinasc2025]

If only annual or monthly totals are available, distribute them across covered days under an explicitly stated within-interval assumption, then aggregate to epidemiological weeks. Preserve totals, leap days and cross-year weeks. Do not divide every year by 52 or every week by an assumed constant month length.

If supplied data are rates, document their denominator and time units before converting to counts. Annual births and deaths are flows; annual population estimates are stocks at reference dates.

### 4.3 Reconcile population and flows without double counting

When the weekly population path is externally fixed, define

$$
M_t^{balance}=N_{t+1}-N_t-B_t+D_t.
$$

This ensures exact accounting. It is a **net migration/reconciliation residual**, not automatically a measurement of migration. It can also contain interpolation effects and discrepancies between data systems.

Do not implement

```text
S_next = S - X + delta_N + births - susceptible_deaths
```

because `delta_N` already includes births and deaths.

If directly estimated migration flows are available, check the population identity against them. Resolve discrepancies upstream or model the measurement uncertainty explicitly; do not silently count both migration and the complete residual.

### 4.4 Allocation across susceptible and non-susceptible residents

In the first age-aggregated implementation, proportional allocation of deaths across S/U can be a documented approximation. The age distributions of these groups may differ, so this is a structural assumption, not a biological fact.

Population totals alone do not identify the susceptibility of migrants. A signed proportional-composition adjustment can be used as a baseline for the reconciliation residual, but it is only an accounting closure. Net flows do not reveal gross inflows/outflows with different immune compositions.

Do not automatically classify all positive reconciliation as susceptible births. Births have already been accounted for separately.

Document the within-week operation order and use it consistently. Do not hide infeasible states with clipping; fail a test and investigate the cause.

### 4.5 Suggested weekly input contract

```text
week_start
week_end
N_start
N_end
births
all_cause_deaths
net_population_reconciliation
population_revision
population_reference_date_rule
birth_source_and_release
death_source_and_release
provisional_data_flag
```

The aggregate metadata can be stored once in a manifest rather than repeated in every row.

## 5. Separate the quantities currently conflated in plots

| Output | Definition | Can decrease over time? |
|---|---|---|
| `cumulative_infections_count[t]` | Sum of modelled infections from analysis start through week t | No |
| `current_susceptible_fraction[t]` | S[t] / N[t] | Yes |
| `current_non_susceptible_fraction[t]` | 1 - S[t] / N[t], under the current protection assumptions | Yes |
| `survey_seropositivity[j]` | Model probability of a positive test in survey j's sampled population/window | Yes across surveys/times |

A current resident immune/non-susceptible fraction is not necessarily the number of historical infections divided by current population. Births, deaths and movement break that equality.

Likewise, a sum of infections divided by a fixed baseline population is a cumulative burden measure, not automatically the probability that a current resident has ever been infected.

Keep the current assay assumptions unchanged in the next iteration, as requested, but describe the serology-to-immune-state link as an approximation. No new assay sensitivity/specificity parameters are required for the population correction.

## 6. Serology and reporting: retain information without overstating identification

### 6.1 Working serology registry

| Site | Collection window | Positives / tested | Intended role |
|---|---|---:|---|
| Juazeiro do Norte | June-December 2018 | 103 / 404 | Local community survey likelihood, not a direct statewide measurement.[^juazeiro] |
| Quixada | June 2018-December 2019 | 289 / 409 | Local survey likelihood; verify the exact assay endpoint, sampling design and collection dates in the final analysis registry.[^quixada] |

Do not infer state representativeness from the number of positive tests. Adding a second local survey provides information on local heterogeneity, but it does not produce a population-weighted statewide sample.

The existing Fortaleza subgroup evidence should not be presented as a clean validation sample for all residents without accounting for its different target population. Do not add it to the primary likelihood simply to narrow statewide uncertainty.

### 6.2 Geographic observation layer

The current site-to-state random-effects approach can be retained as a provisional transport model. It is not a full spatial transmission model.

Two sites provide limited information about a geographic variance. A zero-mean logit offset does not imply that its center equals the population-weighted state proportion. The latter is

$$
a_{state,t}=\frac{\sum_m N_{m,t}a_{m,t}}{\sum_m N_{m,t}}.
$$

Therefore statewide susceptibility can remain sensitive to geographic-offset priors and reporting assumptions even when both local observations are fitted closely.

For a time-constant site offset, a coherent survey-window probability is

$$
p_{j,t}=\operatorname{logit}^{-1}\{\operatorname{logit}(a_t)+u_j\},
\qquad
\bar p_j=\sum_{t\in W_j} w_{j,t}p_{j,t}.
$$

This generally differs from applying the offset after averaging `a[t]`. Audit the current convention; do not silently alter it during a population-only comparison. Use actual collection counts/dates when available. Uniform window weights are an assumption, not an observed sampling schedule. With date-stratified observations, use date-stratified likelihood contributions.

If site observations were used to fit the model, label matching plots **in-sample serology checks**, not external validation. For a genuinely held-out site, do not condition its predicted offset on its held-out serology.[^ppc]

### 6.3 Symptomatic probability and reporting

The observation scale includes

$$
q_t=p_{symp}\rho_t,\qquad E(Y_t)=q_t X_t.
$$

Report overall detection `q[t]` separately from symptomatic-case reporting `rho[t]`.

The previously suggested `Beta(30,28)` for `p_symp` is a working elicited prior, not a literature-validated estimate established by this document. Record its external justification before publication. Do not use the same evidence both to create an informative prior and again as an independent likelihood contribution without accounting for that reuse.

For each important scale parameter, compare prior and posterior and examine correlations with statewide immunity. A weak correlation is not proof of identification; strong priors can obscure a likelihood ridge.

## 7. Immediate execution plan

### Run A: population-data correction only

Suggested identifier: `v2_1_full_popfix`.

- Preserve the original mixed-vintage run and its input manifest.
- Replace the population series only. Keep the full 2015-2025 case window, serology, priors, initial-state convention and existing net-change demographic rule unchanged.
- Use the same interpolation convention in the first comparison where possible. If interpolation is also corrected, label that separately from the source-vintage replacement.
- Verify that no time rescaling or parameter-prior change has slipped into the comparison.
- Run the data tests, compile check and a modest multi-chain pilot before an expensive full fit.
- Inspect chain-specific summaries. If mixing remains severely poor, diagnose rather than indefinitely increasing iterations.

This is a controlled data-ablation experiment, not the final demographic model.

### Run B: explicit demographic accounting

Suggested identifier: `v3_0_demography`.

- Use the corrected population inputs and reconciled births/deaths.
- Replace the old net-change susceptible-entry rule; do not layer births on top of it.
- Keep the remaining scientific components unchanged for this comparison.
- Validate deterministic accounting on synthetic examples without fitting epidemic data.
- Refit the short calibration period first. Then evaluate successive longer prefixes, for example through 2021, 2023 and 2025.

A later-data fit can revise earlier latent states through Bayesian smoothing. To isolate a model change from a data-window change, also compare model versions fitted to the same horizon.

### Run C: targeted identification/computation work, only if needed

Before redesigning the temporal model, inspect where the prior is actually applied. A prior on a reconstructed transmission quantity is not necessarily a direct random walk on `log_hazard_relative`.

Distinguish:

1. **Equivalent reparameterization:** changes computational coordinates while preserving the probability model, including transformations and Jacobians where required.
2. **New scientific model:** changes priors, temporal covariance, reporting flexibility, migration assumptions or importation structure.

Centered and non-centered formulations are alternatives to test, not a universal ranking. An AR(1) process changes assumptions; it is not merely a sampler fix and does not necessarily reduce the number of weekly states.

Do not fix high reconstructed transmission values by imposing an arbitrary cap or by forcing susceptibility upward. Determine whether they are driven by the infection/reporting scale, low susceptible fractions, weak data support, coding errors or the process prior.

## 8. Decision gates for every version

A version can be archived even if it fails. Development can proceed from a documented failure; there is no need to force a known-bad data/model combination to converge first. Promotion to a scientific baseline requires evidence beyond attractive plots.

| Gate | Evaluation | Minimum decision rule |
|---|---|---|
| **G0: inputs and definitions** | Population vintage, epidemiological weeks, missing vs zero cases, provisional data, survey windows, variable definitions | No unexplained mismatch; all sources and assumptions recorded. |
| **G1: accounting and implementation** | State bounds, conservation, cumulative-count monotonicity, date/index mapping, deterministic tests | Tests pass without clipping or hidden rescaling. |
| **G2: computation** | Rank-normalized split R-hat; bulk/tail ESS; Monte Carlo SE; chain traces; divergences; E-BFMI; tree-depth limits | Stable sampling for the requested estimands and relevant latent blocks; no unaddressed warnings. |
| **G3: identification** | Prior/posterior comparison; reporting-infection-immunity trade-offs; geographic-prior dependence; recovery tests | Claims reflect actual information. Weakly identified quantities are labelled or the inferential target is narrowed. |
| **G4: data adequacy** | Replicated weekly observations, wave summaries, quiet periods, residual structure, serology checks | No major systematic discrepancy in the quantities central to the research question. |
| **G5: generalization** | Chronological held-out evaluation; geographically held-out evaluation when appropriate | Required before predictive or transfer claims, not replaced by in-sample fit. |

### Numerical sampling targets

Use these as working diagnostic targets, not proof of a correct scientific model:[^stan-diag][^cmdstan-diag]

- At least four independently initialized chains for a final fit.
- R-hat below 1.01 for relevant nonconstant parameters and derived quantities; audit nuisance blocks as well.
- Bulk ESS at least roughly 100 times the number of chains; with four chains this is 400. Check tail ESS separately, especially for interval endpoints. Apply a comparable tail-ESS working target when intervals are central to the analysis.
- Monte Carlo error small relative to the precision needed for the scientific statement; predefine tolerances for the main susceptibility and reporting summaries.
- Zero post-warmup divergences as the goal. Inspect E-BFMI by chain; values below about 0.3 warrant investigation.
- Report tree-depth hit counts and rates. A few hits are not equivalent to severe R-hat failure, but they must be investigated in context.

Remove warmup before posterior diagnostics. Do not discard inconvenient chains to improve metrics. Separate parameters, transformed parameters and generated quantities; thousands of bad derived elements are not thousands of independent failure mechanisms.

The supplied full-period stress test fails G2. Its 2025 susceptibility or 91% immune fraction should not be treated as established posterior estimates.

### Identification and recovery standards

Simulation-based calibration checks the inference implementation under its own generative model; it does not establish that the model is a correct description of Ceara.[^sbc]

Use repeated simulation experiments spanning plausible low/high detection, heterogeneous survey sites and different demographic turnover levels. A single attractive recovery plot is not sufficient. Record bias, interval coverage, computational failures and which quantities remain prior-dominated.

Do not require uncertainty to become narrower with every version. Better accounting or weaker unsupported assumptions can widen it. Do not require a preferred final immune fraction.

## 9. Standard output package

### Six-panel scientific summary

Maintain consistent labels and shared-axis comparisons where meaningful:

A. Observed cases and replicated-observation intervals. Distinguish these from uncertainty in the fitted expected case count.

B. Latent weekly infection incidence.

C. Current susceptible and non-susceptible fractions, with the state timing convention stated.

D. Effective and model-defined susceptibility-adjusted transmission quantities. Avoid presenting an aggregate inferred quantity as a directly measured biological R0.

E. Symptomatic probability, symptomatic-case reporting, and overall detection.

F. Survey-window statewide quantities, site-specific probabilities and observed local serology, explicitly labelled as calibration/in-sample checks.

### Required supplementary diagnostics

- Population stock and demographic-flow accounting, including the reconciliation residual.
- Separate cumulative infection counts and current immune fractions.
- Chain-specific state/scale summaries and rank/trace plots.
- Prior/posterior comparisons for symptom probability, reporting, geographic variation and transmission-process scales.
- Posterior dependence among infection scale, detection, state immunity and site offsets.
- Residual temporal checks and replicated summaries of peaks, annual totals and quiet periods.[^ppc]
- HMC diagnostics with parameter family, index, mapped date and chain information.
- A run manifest with code/data checksums, Git commit when available, software versions, seed, priors, sampler settings and observation window.

Summarize draws by `model x quantity x date` before drawing ribbons. Test for duplicated rows and accidental reuse of S summaries for reporting plots. Do not interpret plot artifacts as uncertainty.

## 10. Long-term scientific roadmap

### Phase I: trustworthy historical reconstruction

**Question:** What infection burden and susceptible trajectories are consistent with cases, local serology and coherent demography?

Complete the population correction, demographic accounting and inference gates. Assess reporting and serology transport uncertainty. Retain a simple climate-free reference model, but do not assume it uniquely identifies hidden states.

**Deliverable:** An uncertainty-qualified reconstruction with explicit limitations, not a single unquestioned susceptibility curve.

### Phase II: targeted full-history structure

**Question:** Which additional processes are required by the observed recurrence and by the model's failures?

Compare one change at a time when justified: long-term reporting structure, residual temporal covariance, geography or external introductions. Do not conclude that zero/low reported counts prove local extinction. Do not infer that a high model-derived transmission value alone proves missing importation.

Population movement in demographic accounting and pathogen introduction into transmission are different processes; do not reuse one as the other without an explicit model.

**Deliverable:** A defensible full-history baseline. Poorly identified quantities may remain conditional on priors and should be reported that way.

### Phase III: optional empirical climate analysis

**Question:** What descriptive nonlinear and delayed climate-case associations appear after appropriate temporal and spatial adjustment?

Use the municipality climate panel and, where helpful, DLNM/small-area time-series models. This is a complementary exploratory analysis, not a mandatory prerequisite for every renewal extension and not an independent validation when it uses the same cases.

Do not copy case-level DLNM relative risks or case-level lag coefficients directly into a transmission model. Do not turn an estimate from the same data into an independent informative prior without addressing data reuse and the mismatch of estimands.

### Phase IV: joint climate-renewal inference

**Question:** Does climate provide reproducible information about variation in transmission beyond susceptibility, observation processes and a residual temporal model?

Fit climate and latent epidemic states jointly where feasible. A posterior-mean susceptibility trajectory from Phase I is not error-free input. A modular alternative must propagate state uncertainty and clearly state what feedback is excluded.

Compare a climate-free model and a climate-informed model with comparable residual-process assumptions. Climate coefficients and residual temporal components can compete for the same seasonal variation.

Residual variance reduction is model-dependent evidence, not automatically causal attribution. Keep a residual component unless a restricted climate-only model is explicitly being tested as a comparator.

### Phase V: generalization

**Question:** Does the relationship improve calibrated prediction or transport to unobserved periods or regions?

Use rolling-origin/leave-future-out evaluation. A hindcast ending before an outbreak must not initialize from states smoothed using that outbreak's future cases or serology.[^lfo]

Also prevent leakage from future reporting estimates, fitted future random effects or climate transformations estimated with held-out observations. State whether realized future weather is supplied: this gives a weather-conditional hindcast, not an operational weather forecast.

Spatial validation needs a clear treatment of a new region's baseline and random effects. A two-site serology analysis alone cannot establish national transportability.

### Phase VI: conditional future scenarios

**Question:** How does future epidemic potential or burden change under specified climate, population, immunity, introduction and intervention assumptions?

Proceed only at the level supported by the validation evidence. If burden is not identifiable, report more limited quantities rather than presenting precise future case counts.

Propagate posterior states and recent infection history, demographic uncertainty, climate-model/scenario uncertainty and residual-process uncertainty. Clearly define future reporting when converting infections to reported cases.

Assess exposure-range support, bias adjustment and temporal/spatial scale compatibility before applying a historical climate function to future conditions. Separate climate-only contrasts from changing demography or intervention contrasts.

Call these **conditional scenario projections**, not unconditional predictions of a particular outbreak in a distant year.[^ipcc]

## 11. Suggested implementation milestones

```text
D0  Population source/vintage fix and audit
D1  Controlled v2.1 full-period population-only comparison
D2  Explicit demographic accounting with deterministic tests
D3  Matched short-period and progressively longer fits
D4  Targeted computation/identification improvements
D5  Full-history reconstruction accepted with stated limits
C0  Optional empirical climate analysis
C1  Joint climate-renewal comparison
C2  Chronological and spatial generalization tests
P0  Conditional future scenario analysis, where justified
```

The model name and data version should be separate in run metadata. A data correction is not automatically a new scientific model. Proposed labels such as `v3_0_demography` are organizational aids, not claims that the model has passed its gates.

## 12. Immediate agent task brief

```text
Preserve all existing model files and outputs.

1. Audit and freeze the single-vintage population input and its metadata.
2. Create a population-only full-period comparison with all other scientific
   assumptions held fixed. Report any unavoidable interpolation change.
3. Do not interpret pooled posterior trajectories unless sampling is adequate.
4. Prepare an explicit demographic input table: N_start, N_end, births,
   deaths, and net migration/reconciliation; do not double count delta_N.
5. Implement and test population/state accounting before epidemic fitting.
6. Separate cumulative infection counts from current immune fractions.
7. Refit matched short periods before expensive full-period fits.
8. Examine infection-detection-serology scale dependence by chain.
9. Propose one justified next model change at a time, with a testable purpose.
10. Do not automatically add climate, a new importation process, immunity
    waning, arbitrary R0 limits, or a more flexible reporting model.

Before any large run, report the exact run window, data revision, changed
components, priors, sampler configuration, expected outputs and resource use.
```

## References and provenance

Project-specific numerical results and version descriptions above come from user-supplied run summaries and figures in the project discussion, not from independent re-analysis of posterior draws. The corrected population table is likewise project-supplied. Link archived run manifests and exact commits when adding this document to the repository.

[^ibge]: IBGE. *Population Projections, 2024 Revision*, official release directory and methodology links. https://ftp.ibge.gov.br/Projecao_da_Populacao/Projecao_da_Populacao_2024/

[^sinasc]: Brazilian Ministry of Health. *Sistema de Informacao sobre Nascidos Vivos (SINASC)*. https://www.gov.br/saude/pt-br/composicao/svsa/sistemas-de-informacao/sinasc

[^sim]: Brazilian Ministry of Health. *Sistema de Informacao sobre Mortalidade (SIM)*, official data catalogue. https://dadosabertos.saude.gov.br/dataset/sim

[^sinasc2025]: Brazilian Ministry of Health. *Nascidos Vivos - 2025 - preliminar*, CSV resource; status checked 2026-09-07. https://dadosabertos.saude.gov.br/dataset/sistema-de-informacao-sobre-nascidos-vivos-sinasc/resource/7f63d984-97aa-44da-891d-72b1dd6c5122

[^juazeiro]: Barreto et al. (2020). *Seroprevalence, spatial dispersion and factors associated with flavivirus and chikungunya infection in a risk area: a population-based seroprevalence study in Brazil*. BMC Infectious Diseases. https://link.springer.com/article/10.1186/s12879-020-05611-5

[^quixada]: Braga. *Inquerito soroepidemiologico sobre Chikungunya, Dengue e Zika no municipio de Quixada, Ceara*. Universidade Federal do Ceara institutional repository. https://repositorio.ufc.br/handle/riufc/51582

[^stan-diag]: Stan Development Team. *How to Diagnose and Resolve Convergence Problems*. https://mc-stan.org/learn-stan/diagnostics-warnings.html

[^cmdstan-diag]: CmdStan User's Guide. *Diagnosing Biased Hamiltonian Monte Carlo Inferences*. https://mc-stan.org/docs/cmdstan-guide/diagnose_utility.html

[^ppc]: Stan User's Guide. *Posterior and Prior Predictive Checks*. https://mc-stan.org/docs/stan-users-guide/posterior-predictive-checks.html

[^sbc]: Stan User's Guide. *Simulation-Based Calibration Checking*. https://mc-stan.org/docs/stan-users-guide/simulation-based-calibration.html

[^lfo]: Stan/loo. *Approximate leave-future-out cross-validation for Bayesian time series models*. https://mc-stan.org/loo/articles/loo2-lfo.html

[^ipcc]: IPCC AR6 Working Group II. *Annex II: Glossary*, definitions of projections, predictions and scenarios. https://www.ipcc.ch/report/ar6/wg2/chapter/annex-ii/

### Recommended conceptual reading

Steyn, Parag, Thompson and Donnelly (2025). *A Primer on Inference and Prediction With Epidemic Renewal Models and Sequential Monte Carlo*. Statistics in Medicine. DOI: 10.1002/sim.70204. https://doi.org/10.1002/sim.70204
