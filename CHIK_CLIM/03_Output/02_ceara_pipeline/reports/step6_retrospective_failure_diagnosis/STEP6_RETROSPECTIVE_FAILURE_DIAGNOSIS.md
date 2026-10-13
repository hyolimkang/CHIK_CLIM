# STEP 6 -- Retrospective Failure Diagnosis: why does the Ceará chikungunya
# renewal model underpredict the observed 2022 resurgence?

Mechanistic diagnostic only -- no parameters were changed to improve the
2022 fit, no new Stan fitting occurred, no importation rate was selected
or tuned. All findings below are quantitative outputs of scripts in
`02_Script/02_ceara_pipeline/15_step6_retrospective_failure_diagnosis/`,
reusing the training-only posteriors and holdout simulations already
produced by `12_future_year_effect_validation` and
`13_retrospective_resurgence_diagnostic`. Nothing accepted was modified.
Not committed.

**PROVISIONAL**: the HMC reliability audit found low but nonzero
divergences and incomplete within-chain mixing in tau_year, most
pronounced in its upper tail (`HMC_AUDIT_model_level.csv`,
`HMC_RELIABILITY_AUDIT.md`). All quantitative findings below -- medians,
intervals, probabilities -- should be read as provisional until that
reliability is fully established, as explicitly flagged on every figure.

## Phase 0 -- Audit of the retrospective simulation design

`STEP6_design_and_implementation_audit.csv` (11-point check, extending
the audit already performed in `13_retrospective_resurgence_diagnostic/
PHASE0_design_audit.csv`):

| Item | Status | Finding |
|---|---|---|
| Training/serology | PASS | H1/H2 use only weeks through their own cutoffs; eligible 2018 serology included in both. |
| Observed infection seed | PASS | Training data have is_seed=X_seed=0; no observed holdout infection seed enters simulation. |
| Full-period posterior init. | PASS | Holdouts read only their own `ce_tau_year_training_{H1,H2}.rds` terminal states/history. |
| Terminal S/U/history | PASS | `advance_terminal_state()` supplies the one demographic step the Stan model's `t<N` guard skips; 8-week posterior X history carried forward. |
| tau_year draw | PASS-WITH-HMC-CAVEAT | G2 uses tau_year from the matching training-only posterior draw; quantitative results remain provisional. |
| **Future year-effect recentering** | **FAIL** | Future z_year values are re-centred (sum to zero) WITHIN each finite holdout horizon. This shrinks the marginal SD of each year's realised effect from tau_year to tau_year*sqrt((K-1)/K), where K=holdout years -- 0.913 for H1 (K=6), 0.866 for H2 (K=4). Quantified, disclosed, and small (median within-horizon SD 0.147/0.142 vs. median tau_year used 0.158/0.164 for H1/H2 respectively -- `STEP6_future_year_recentering_audit.csv`) relative to the ~170-fold entry-infection gap identified below. |
| R0/Reff/seasonal phase | PASS | Reconstructed log-R0 components (baseline+year+seasonal+climate) reproduce the saved simulated R0 to within 4.4e-16 (floating-point only) for both H1 (n=31,200) and H2 (n=20,800) -- `STEP6_R0_component_reconstruction_check.csv`. |
| Demography/population | PASS | S+U=N exactly (max abs error = 0); all states nonnegative. |
| Importation accounting | PASS | M_t enters once as M_t/N_t in the hazard -- **introduction PRESSURE, not realised imported infections or added cases**, per the task's own clarification. |
| Observation process | PASS | Reported cases ~ NegBinomial(mean=q*X, phi_obs), q=0.05 fixed; aggregate simulated reported-to-latent ratio = 0.05004 (matches q to 4 significant figures) -- `STEP6_numerical_consistency.csv`. |
| Climate information | PASS-WITH-CAVEAT | Hindcasts use observed holdout climate standardised with full-2015-2025 scaling constants -- observed-climate-conditional with disclosed preprocessing leakage (same finding, independently reproduced, as `13_retrospective_resurgence_diagnostic/PHASE0_climatology_leakage_magnitude.csv`). |

No serious leakage requiring a fix before scientific interpretation was
found. The one FAIL (recentering) is real but small in magnitude and is
carried forward as a disclosed limitation, not corrected here (correcting
it would change the generator's own definition and requires a
design decision out of Step 6's scope).

## Phase 1-2 -- Susceptibility entering 2022

`STEP6_training_cutoff_susceptibility.csv`, `STEP6_Jan2022_susceptibility_G2.csv`:

| Quantity | H1 | H2 |
|---|---|---|
| S/N at own training cutoff | 0.739 [0.714, 0.763] | 0.744 [0.722, 0.770] |
| S/N, first week of 2022 (G2) | **0.741 [0.597, 0.767]** | **0.744 [0.722, 0.770]** |

**Finding: susceptibility entering 2022 is essentially indistinguishable
between H1 and H2** -- medians differ by 0.003 (0.4%), well within each
other's 95% ranges. H1's interval is wider (because it is a genuinely
simulated 2-year-ahead quantity, not a direct posterior read-out like
H2's), but **not shifted low or high** -- H1 generates neither excessive
nor insufficient pre-2022 infections relative to H2 on average. This
directly answers Phase 2's question: **outcome (C) -- susceptibility
looks reasonable; the main problem lies elsewhere.**

## Phase 3 -- Transmission potential (R0, Reff) in 2022

`STEP6_2022_headline_mechanistic_metrics.csv` (G2, estimated tau_year):

| Quantity | H1 | H2 |
|---|---|---|
| Mean R0 | 1.26 [0.82, 1.72] | 1.38 [1.04, 1.91] |
| Peak R0 | 1.90 [1.19, 3.47] | 1.91 [1.23, 3.20] |
| Mean Reff | 0.93 [0.57, 1.30] | 1.01 [0.78, 1.24] |
| Peak Reff | 1.41 [0.79, 2.58] | 1.40 [0.91, 2.38] |
| Weeks with Reff>1 | 17 [0, 33] | 23 [0, 29] |
| Longest consecutive Reff>1 run | 17 [0, 24] | 19 [0, 29] |

**Finding: peak transmission potential is essentially identical between
H1 and H2** (peak R0 1.90 vs. 1.91; peak Reff 1.41 vs. 1.40). H1 shows
modestly *fewer* weeks with Reff>1 (17 vs. 23) and a shorter typical
favourable run, a real but much smaller difference than the entry-infection
gap found below. FIG21 shows the full weekly R0(t)/Reff(t) trajectories
are visually near-overlapping between holdouts. Reff>1 is confirmed,
consistent with the task's framing, as evidence of POTENTIAL growth, not
a guarantee of an outbreak -- which is exactly what Phase 4 investigates.

## Phase 4 -- Introduction, standing infections, and establishment

This is where H1 and H2 diverge sharply. `STEP6_entry2022_infection_history.csv`
(reusing the generation-weighted infectiousness definition `sum(w*X)`
over the 8 weeks preceding 2022-01-02):

| Source | Median recent X (8wk) | Median weighted infectiousness | P(=0) | Median X, first week of 2022 |
|---|---|---|---|---|
| H1: simulated final 8 weeks of 2021 | **37** | **4.77** [0.43, 306.5] | 0.2% | **6** [0, 277] |
| H2: training-posterior final 8 weeks of 2021 | **6,385** | **886.4** [671, 1,119] | 0% | **974** [643, 1,531] |

**H2 enters 2022 with ~170-190x more standing infectiousness than H1**,
because H1's 2020-2021 infection history is the product of its OWN
2-year stochastic holdout simulation (no real surveillance information
at all, generated purely from the training-only posterior's parameters
plus `lambda_import=1`), whereas H2's is read directly from the
training-only posterior's OWN fitted state, which was estimated against
the REAL observed 2015-2021 case series via the likelihood. This is not
a susceptibility difference (Phase 2) and not a peak-transmission-potential
difference (Phase 3) -- it is a difference in how much infection is
already circulating when the 2022 season's favourable conditions arrive.

`STEP6_2022_establishment_summary.csv` and `STEP6_2022_mechanism_classification.csv`
(G2, estimated tau_year; established := >=1 week with Reff>1 AND M_t>0 in
2022, followed by >=4 consecutive weeks where local transmission hazard
exceeds the introduction hazard -- a mechanistic, not calibrated, rule
fixed before any outcome was inspected):

| | H1/G2 | H2/G2 |
|---|---|---|
| P(any latent infection in 2022) | 100% | 100% |
| P(favourable pressure: Reff>1 & M_t>0) | 90.9% | 94.0% |
| P(established, per the rule above) | 90.9% | 94.0% |
| **P(major annual epidemic)** | **10.55%** | **99.55%** |
| Mechanism class "established; no major epidemic" | **80.35%** | 0% |

**Finding: establishment itself is common and similar in both holdouts
(90.9% vs. 94.0%) -- introductions are neither absent nor badly timed
relative to favourable Reff in either holdout.** The gap opens entirely
AFTER establishment: 80.35% of H1/G2 trajectories establish local
transmission but still fail to cross the major-epidemic threshold within
the calendar year, versus essentially 0% of H2/G2 trajectories (FIG22).
This is exactly what the ~170-fold entry-infectiousness gap predicts: a
chain that establishes from ~6 standing infections needs many more weeks
of sustained growth to reach epidemic scale than one that establishes
from ~974, and the 2022 transmission season is finite.

Within H1/G2 specifically, trajectories that DO become major epidemics
differ descriptively from those that don't (`STEP6_G2_2022_outcome_contrasts.csv`,
FIG23A; conditioning on simulated outcome, NOT a causal claim): higher
cumulative latent infections (median 20,617 vs. 797), higher peak Reff
(1.93 vs. 1.36), higher realised year_effect[2022] (+0.068 vs. -0.025)
and higher A_year[2022] (1.23 vs. 0.98). This is consistent with -- not
independent proof of -- the entry-infection mechanism: a major epidemic
in H1 requires BOTH a favourable stochastic year draw AND enough standing
infections to exploit it within one season.

## Phase 5 -- Dominant mechanism (descriptive, not causal)

**Of the four candidate explanations (A: susceptibility, B: transmission
potential, C: introduction/establishment, D: implementation), the evidence
most strongly and specifically implicates a facet of (C) that is better
described as "entry infection-history" than "introduction timing" per
se**: introductions are not absent and establishment is not rare in
either holdout; what differs by ~2 orders of magnitude is how much
infection has already accumulated by the time 2022 begins, which is a
direct consequence of H1 having to simulate TWO FULL YEARS (2020-2021)
from the generative model alone (no real data) under a constant
`lambda_import=1`, during which Reff hovers at or below 1 for much of
the period (FIG21) -- a regime where stochastic fadeout is a real and
frequent outcome. (A) susceptibility and (B) peak transmission potential
are both ruled out as primary drivers by direct, matched comparison
(Phase 2-3 above). (D) implementation problems are largely ruled out by
the numerical consistency checks, with the one disclosed exception
(recentering) being too small in magnitude to explain a ~170-fold gap.

**Important scope limitation, stated explicitly per the task's own
caution**: H1 and H2 were fitted on different training windows with
different posteriors and are NOT a paired causal experiment -- the
comparison above is a reconstruction of two different retrospective
designs, not a controlled contrast. Where this report DOES use a valid
common-random-number pairing is the G1-vs-G2 generator substitution
WITHIN each holdout (same training posterior, same z_year/z_season/
climate/import draws per matched posterior_draw_id and
stochastic_replicate_id -- `STEP6_controlled_G1_G2_contrasts.csv`); this
pairing is NOT valid across H1 and H2, and is not used that way anywhere
in this report.

## Phase 6 -- Predictive performance

`STEP6_annual_predictive_performance.csv`, `STEP6_multiyear_cumulative_performance.csv`,
`STEP6_2022_exceedance_probability_mc.csv`, `STEP6_G1_vs_G2_interval_widths.csv`:

- **2022 exceedance probability** (simulated annual reported cases >=
  observed 50,986), with draw-cluster bootstrap 95% uncertainty: H1/G1
  2.9% [0.1%, 6.6%]; **H1/G2 0.3% [0%, 0.9%]**; H2/G1 30.4% [21.8%, 39.4%];
  **H2/G2 11.9% [6.0%, 18.5%]**.
- **Multi-year cumulative holdout performance**: H1/G1 predictive median
  53,629 (close to observed 56,903; WIS=19,452); **H1/G2 median 1,600 (34x
  low; WIS=34,199 -- WORSE than G1 despite G2's narrower, more "precise"
  interval)**; H2/G1 median 69,508 (WIS=21,214); **H2/G2 median 16,103
  (3.4x low; WIS=16,303 -- slightly BETTER than G1 for this holdout)**.
- **G1's intervals are 2-66x wider than G2's** across all holdout-years
  (`STEP6_G1_vs_G2_interval_widths.csv`) -- G1's apparent "success" at
  covering 2022 is substantially a product of its much larger, less
  informative predictive spread, not of a better-specified mechanism;
  this is exactly the caution the task asked not to overlook.
- **Containing the observed value in a 95% interval is not, by itself,
  evidence of good prediction here**: both G1 and G2 technically "cover"
  2022 for H1, but WIS -- a proper score that penalises both bias and
  sharpness -- ranks G1 better for H1 and G2 (modestly) better for H2,
  a holdout-dependent result that would be invisible from coverage alone.

## Answers to the seven questions

**1. Why does H1 underpredict the observed 2022 resurgence?**
Because H1 must simulate two full years (2020-2021) of epidemic dynamics
purely from the generative model, with no real case information, under
a constant importation-pressure assumption, and in the median trajectory
this produces a near-fadeout state (median ~6 standing infections,
~170x below H2) by the start of 2022 -- not because susceptibility or
peak transmission potential (R0/Reff) are wrong; both are statistically
indistinguishable from H2's.

**2. Why does H2 generate a different predictive distribution?**
H2's holdout begins exactly at 2022-01-02, so its initial state is read
directly from the training-only posterior fitted to the REAL observed
2015-2021 surveillance (via the Bayesian likelihood), carrying forward
the actual, apparently substantial, 2020-2021 epidemic signal rather
than reconstructing it from scratch stochastically.

**3. Is the difference primarily susceptibility, transmission, importation,
or their interaction?**
Primarily a consequence of entry infection-history, a facet closest to
introduction/establishment (C) but distinct from the task's default
notion of "introduction timing" -- it is about standing infection SCALE
at a fixed, favourable calendar point, not about missed introduction
windows (introductions and establishment are common and similar between
H1 and H2). Susceptibility and peak R0/Reff are not implicated.

**4. Does estimated tau_year appear too restrictive, or is this conclusion
unsupported?**
**Unsupported as the primary explanation.** Realised year_effect[2022]
values ARE somewhat lower, on average, among H1/G2's non-major
trajectories than among its major ones (-0.025 vs. +0.068, FIG23A) --
but this is a descriptive, outcome-conditioned association, not a
demonstration that tau_year itself is too small: susceptibility and
transmission potential match H2 closely, and the entry-infection gap
(170x) is far larger in both magnitude and mechanistic specificity than
any plausible effect of tau_year's scale. The recentering FAIL (Phase 0)
is a real, quantified, but small (<15%) contributor, not a dominant one.

**5. Is lambda_import=1 a plausible working assumption, or is its
adequacy unresolved?**
**Unresolved by Step 6's design, by instruction** -- no importation
sensitivity grid was run here. What Step 6 DOES show is that lambda_import=1,
combined with a Reff(t) trajectory that spends much of 2020-2021 at or
below 1 (FIG21), produces a MEDIAN outcome close to fadeout by the start
of 2022 in H1's unconditioned 2-year simulation. Whether a different
(not yet selected) importation rate or a structurally different
introduction process would change this is explicitly a STEP 7 question.

**6. Can the current G2 generator be retained for further investigation?**
Yes, with the caveats disclosed throughout: it is mechanistically
well-behaved (numerically consistent R0/Reff/demography/observation
layers; PASS on 9 of 11 audit items), its tau_year posterior location is
robust (per the separate HMC audit), and its qualitative failure mode
(entry-infection-history-driven underprediction in H1, much closer
agreement in H2) is now mechanistically understood rather than a black
box. It should not yet be used for unconditional point forecasts given
the WIS/exceedance results above and the still-provisional HMC status.

**7. What is the smallest scientifically justified next step?**
A **STEP 7 importation-pressure sensitivity analysis** (a prespecified
lambda_import grid, NOT tuned to 2022 outcomes, exactly as this task's
scope already anticipates) is the most direct way to test whether the
entry-infection-history gap identified here is sensitive to the
constant-rate introduction assumption -- paired with an explicit test of
whether a constant-rate Poisson process is the right STRUCTURE (vs. a
discrete reintroduction-event process), since Step 6's finding describes
a magnitude gap consistent with either a rate problem or a structural
one. The task's own permitted "optional zero-introduction diagnostic"
was NOT run in this Step 6 pass (it would have required a new simulation
beyond what was already available) and is a natural, minimal first
probe before a full STEP 7 grid.

## Final decision

**B -- Introduction pressure is an important unresolved mechanism; proceed
to STEP 7 importation sensitivity and validation.**

This is an evidence-driven conclusion, not a default: (A) susceptibility
and the transmission-potential component of (C)/(B) -- peak R0/Reff --
were directly, quantitatively ruled out as primary drivers by matched
comparison (Phase 2-3); (D) implementation problems were ruled out by
numerical consistency checks, with one small (<15%), disclosed exception
too minor to explain a ~170-fold entry-infection gap. The evidence
converges specifically and repeatedly (entry-infection-history table,
establishment-vs-major gap, mechanism classification, outcome contrasts)
on introduction-pressure-adjacent dynamics during the two years H1 must
simulate without any real data. (C) "transmission or susceptibility
structure needs reconsideration" is explicitly NOT supported, since
neither was found to differ materially between H1 and H2. (D) "multiple
explanations remain indistinguishable" does not apply either: the
evidence discriminates cleanly.

## Files, scripts, and reused fits

```
02_Script/02_ceara_pipeline/15_step6_retrospective_failure_diagnosis/
  00_config_step6.R, 00_run_step6.R, 01_build_diagnostics.R, 02_make_figures_and_report.R

Reused (NOT refit, NOT modified):
  03_Output/.../model_fits/future_year_effect_validation/ce_tau_year_training_{H1,H2}.rds
  03_Output/.../results/future_year_effect_validation/PHASE3_4_holdout_simulation_results.rds, PHASE3_4_observed_holdout_summaries.rds
  03_Output/.../results/retrospective_resurgence_diagnostic/PHASE1_diagnostic_simulation_results.rds (bit-for-bit-verified extended trajectories: S,U,N,X_latent,Y_reported,R0,R_eff,M_import,lambda_t)
  03_Output/.../tables/retrospective_resurgence_diagnostic/PHASE0_design_audit.csv, PHASE0_climatology_leakage_magnitude.csv, PHASE1_integrity_check_vs_original.csv
  03_Output/.../tables/hmc_reliability_audit/HMC_AUDIT_model_level.csv

New outputs (this report's own, isolated):
  03_Output/02_ceara_pipeline/{tables,results}/step6_retrospective_failure_diagnosis/STEP6_*.csv (23 tables/results)
  03_Output/02_ceara_pipeline/figures/step6_retrospective_failure_diagnosis/FIG21-23.{pdf,png@600dpi}
  03_Output/02_ceara_pipeline/reports/step6_retrospective_failure_diagnosis/STEP6_RETROSPECTIVE_FAILURE_DIAGNOSIS.md (this file)
```

**No new Stan fitting was performed.** If STEP 7 proceeds, it will need
no additional Stan fitting either (it reuses the same training-only
posteriors); it requires only new forward stochastic simulations across
a prespecified lambda_import grid, at a compute cost comparable to the
diagnostic re-simulations already run in `13_retrospective_resurgence_diagnostic`
and this module (minutes, not a new HMC fit).

## Uncertainty and limitations

- All quantitative findings are **PROVISIONAL** pending full resolution
  of the tau_year HMC reliability caveats (mild divergences, upper-tail
  mixing; `HMC_RELIABILITY_AUDIT.md`).
- H1 and H2 are not a paired experiment (different training posteriors);
  all H1-vs-H2 comparisons here are reconstructions of two designs, not
  a controlled contrast. Only the within-holdout G1-vs-G2 substitution is
  a valid common-random-number pairing.
- The climate-anomaly standardisation leakage (PASS-WITH-CAVEAT) and the
  future year-effect recentering (FAIL) are both real, disclosed, and
  judged too small to be the dominant mechanism, but neither has been
  corrected, and their combined contribution has not been formally
  bounded.
- The entry-infection-history finding is strong and consistent across
  multiple independent tables (entry history, establishment summary,
  mechanism classification, outcome contrasts) but remains a descriptive,
  not causal, account of why H1 and H2 differ -- it explains the
  structural asymmetry between the two retrospective designs, not a
  manipulable causal pathway within a single model.
- Major-epidemic classification reuses the project's already-accepted
  absolute/incidence thresholds (100 cases, 5/100k), not re-derived or
  tuned here.

## Recommended next action

Run a **STEP 7 prespecified lambda_import sensitivity grid** (values
fixed before inspecting 2022 outcomes, as in `13_retrospective_resurgence_diagnostic`'s
own Experiment A) on H1 specifically, to test whether raising introduction
pressure materially closes the entry-infection-history gap identified
here -- and, in parallel, resolve the HMC reliability caveats on tau_year
(adapt_delta/parameterisation work already scoped in `HMC_RELIABILITY_AUDIT.md`)
before treating any STEP 7 result as more than provisional.
