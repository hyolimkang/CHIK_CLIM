# Pernambuco Spatial-0 -- U14 Recife Serology Likelihood Audit

**Spatial-0 is NOT yet frozen because Recife U14 validation failed.** No
production run, no vaccine projection, and no model freeze has occurred.
This audit uses ONLY the 3 chains that already reached post-warmup
sampling (chains 2, 4, 5 -- LOW/INTERMEDIATE/HIGH branches, confirmed to
agree with each other in script 38); no new HMC run was performed.

## 1. Is U14 actually in the likelihood? YES -- implementation is correct.

Exact values passed to Stan (`PE_5strata_stan_data.rds`):

| Field | Value |
|---|---|
| `RECIFE_INDEX` | 1 (stratum `"1_Recife"` -- confirmed correct, not Metropolitana_remainder or any other stratum) |
| `J_sero` | 1 |
| `sero_n_positive` | 770 |
| `sero_n_tested` | 2070 |
| survey prevalence | 0.3720 |
| model weeks used | 2018-08-05 to 2019-02-24 (30 weeks; reported survey window 2018-08-01 to 2019-02-28 -- matches to within weekly-binning resolution) |
| `kappa_sero` | 50 (fixed) |

Exact Stan statements (`renewal_pernambuco_v4_9_spatial0.stan`):

```stan
// transformed parameters -- links p_site_window to the RECIFE stratum only
for (j in 1:J_sero) {
  int start = sero_window_start_idx[j];
  int nwk = sero_window_n_weeks[j];
  real acc = 0;
  for (k in 1:nwk) {
    int t = start + k - 1;
    acc += fmin(1 - 1e-9, fmax(1e-9, immune_prop[t, RECIFE_INDEX]));
  }
  p_site_window[j] = acc / nwk;
  alpha_sero_site[j] = p_site_window[j] * kappa_sero;
  beta_sero_site[j] = (1 - p_site_window[j]) * kappa_sero;
}

// model block -- the likelihood contribution itself
if (J_sero > 0) {
  for (j in 1:J_sero) {
    sero_n_positive[j] ~ beta_binomial(sero_n_tested[j], alpha_sero_site[j], beta_sero_site[j]);
  }
}
```

This is a genuine `~` sampling statement inside the `model` block (not
merely computed in `generated quantities`), it indexes `immune_prop[t,
RECIFE_INDEX]` with `RECIFE_INDEX=1` (the Recife stratum, not statewide,
not Metropolitana_remainder), and the window indices land on the correct
calendar dates. **U14 is genuinely and correctly wired into the
likelihood.**

## 2. Does the likelihood use the same latent prevalence shown in panel F? YES.

Independently recomputing `p_Recife_U14` outside Stan, directly from the
stored `immune_prop[iter, t, RECIFE_INDEX]` draws averaged over the same
30-week window:

| Quantity | Median |
|---|---|
| Stan's `p_site_window` (used in the likelihood) | 0.025731 |
| Independently recomputed `p_Recife_U14` | 0.025731 |
| Max absolute difference | 6.7e-8 (floating-point noise only) |

**No indexing or computational error.** Panel F's plotted model value IS
the value the likelihood actually sees.

## 3. Beta-binomial weight audit

`kappa_sero = 50` (fixed, not estimated) is a substantially DOWN-WEIGHTED
observation model relative to a simple binomial survey of n=2070:

| Quantity | Value |
|---|---|
| Implied intraclass correlation, rho = 1/(kappa+1) | 0.0196 |
| Implied EFFECTIVE sample size at n=2070 | 49.8 (2.4% of the raw 2070) |
| log-lik of U14 at model median p=0.0257 | -24.71 |
| log-lik of U14 at observed p=0.372 | -5.89 |
| Log-lik ratio (observed favoured by) | **18.82 nats** |
| Posterior predictive y_rep/n: median [95% PI] | 0.0198 [0.0010, 0.0877] |
| Max of 1500 posterior predictive draws | 0.198 |
| Is observed 0.372 within the 95% PI? | **NO** |

Even though the likelihood is already heavily down-weighted (effective
n approx 50, not 2070), it STILL prefers the observed value by 18.8 nats
-- and the posterior predictive check fails outright: **not even the
single most extreme of 1500 posterior predictive draws reaches half of
the observed prevalence.** A weak likelihood should still visibly pull
the posterior toward the data if the conflict were mild; here it visibly
cannot, because it is comprehensively outvoted by the 4 other regions'
case likelihoods sharing the same `q_PE`. This point counts against
"likelihood merely too weak" (B) as a *sufficient* explanation and toward
a *structural* conflict (C) -- see Section 6.

## 4. Original survey uncertainty

In-repo documentation (`02_build_pe_serology_window.R`) records only: household-based
stratified random survey by socioeconomic level, IgM or IgG, n=2070,
positive=770 (raw prevalence 0.372), reported 95% CI [0.340, 0.404]. No
design weights, no age-stratified breakdown, and no age range are recorded
in-repo beyond this. This is a genuine documentation gap -- flagged, not
filled in with assumptions.

What CAN be derived from the reported CI itself:

| Quantity | Value |
|---|---|
| Naive binomial 95% half-width (n=2070, p=0.372) | 0.0208 |
| Reported 95% half-width | 0.0320 |
| Implied design effect | **2.36x** |

A design effect of approximately 2.4 (consistent with household
clustering in a stratified survey) is unremarkable and does NOT
materially change the conclusion: even inflating the survey's
uncertainty by 2.36x, the reported point estimate (37.2%) remains
enormously far from the model's median (2.6%) and from the maximum of
1500 posterior predictive draws (19.8%). No plausible design-based
correction closes this gap.

## 5. Age / geographic compatibility

Age range and precise geographic sampling frame are NOT recorded in-repo
beyond "Recife household residents." This cannot be fully verified from
project data alone and is flagged as an open provenance question for the
original P10 source. However, quantitatively: to reconcile an all-age
model prevalence of 2.6% with an observed 37.2% via an age-restriction
argument alone would require the surveyed age band to have a prevalence
roughly 14x the model's all-age Recife average -- i.e., a level of
age-concentration for which there is no supporting evidence and which
would be biologically extreme for an arbovirus with broadly age-uniform
exposure risk. Geographic frame mismatch (survey drawn from a
sub-area of Recife with atypically high transmission) is possible but
cannot explain the gap either without implying the REST of Recife had
substantially negative-tending correction, which is not plausible. No
eta_geo or other ad hoc absorption term was used to paper over this (per
instruction).

## 6. Case-serology conflict diagnostic

Recife's OWN cumulative reported case series (not the shared model
q_PE) directly implies what q would be needed to reconcile it with the
observed 37.2% prevalence, via the same crude case-implied-attack
convention used throughout this project's diagnostics:

| Quantity | Value |
|---|---|
| Recife cumulative reported cases through the U14 window (2015-01-04 to 2019-02-24) | 6,293 |
| Recife population at that time | 1,645,727 |
| Crude cumulative reported incidence | 0.382% |
| **q_Recife needed to reproduce 37.2% prevalence** (simple division) | **0.0103** |
| Spatial-0's fitted SHARED q_PE (median) | 0.0623 |
| Ratio | shared q_PE is **6.0x higher** than the q that would reconcile Recife alone |

Fixed-q profile (crude, descriptive -- what Recife's prevalence-by-U14-window
would be under a range of q values, holding Recife's own case series fixed):

| q | Implied Recife prevalence |
|---|---|
| 0.010 | 38.2% (matches observed) |
| 0.020 | 19.1% |
| 0.050 | 7.6% |
| **0.062 (fitted)** | **6.2%** |
| 0.100 | 3.8% |
| 0.200 | 1.9% |

This crude linear profile (6.2% at q=0.062) is in the same order of
magnitude as, though somewhat higher than, the full nonlinear
renewal-model estimate (2.6%) -- the difference is expected, since the
full model's susceptible-depletion feedback and the beta-binomial
likelihood's (weak) pull toward the data both act on the exact figure,
while the crude profile is a simple linear cross-check. Both agree
Recife's own case trajectory, under ANY q remotely compatible with the
other 4 regions' case series (approx 0.04-0.10), cannot reach anywhere
near 37% -- **Recife would need a REGION-SPECIFIC q approximately 6x
lower than the shared estimate**, which the current one-common-q
architecture (Section G, by design) cannot supply.

This is diagnostic only. No region-specific q was introduced into the
inferential model.

## 7. Decision

**Primary diagnosis: C. GENUINE CASE-SEROLOGY STRUCTURAL CONFLICT.**

- (A) Implementation error: **RULED OUT.** Data indexing, window dates,
  and the latent-prevalence computation used by the likelihood are all
  verified correct (Sections 1-2).
- (B) Likelihood too weak: **contributing factor, but not sufficient by
  itself.** kappa_sero=50 does heavily down-weight U14 (effective n approx
  50), but even this weak likelihood assigns the observed value 18.8 more
  nats of log-likelihood than the model's fitted value, and the posterior
  predictive check fails completely (observed value outside the full
  range of 1500 draws). A merely-weak-but-compatible likelihood would show
  some visible pull toward the data and a border-line-passing posterior
  predictive check; this shows neither.
- (C) **Genuine structural conflict: PRIMARY DIAGNOSIS.** Recife's own
  case series, under any q compatible with the other 4 regions
  (approx 0.04-0.10), implies a cumulative attack fraction an order of
  magnitude below the observed serosurvey. The ONE-shared-q_PE
  architecture (deliberately, per the pre-registered design) cannot let
  Recife's ascertainment differ from the state, so the shared parameter
  is pulled by the other 4 regions' case likelihoods (which vastly
  outweigh one serology point in total likelihood mass) toward a value
  that is far too high to be consistent with Recife's own reported
  cases and observed serology simultaneously.
- (D) Mixed/unclear: not selected -- the evidence in Sections 1-3 and 6 is
  internally consistent and points to one primary mechanism.

## 8. Recommended minimum next model change

Do not implement yet (per instruction) -- reporting only. The minimum
change suggested by this diagnostic is to allow ascertainment (or an
equivalent local-intensity adjustment) to differ for the Recife stratum
specifically, since the conflict is structurally located there (i.e. some
form of the geographic-offset mechanism the homogeneous model already
used, or a Recife-specific q, would be the natural fix) -- but this is a
recommendation for review, not an action taken here.

## 9. Outputs

- `PE_U14_LIKELIHOOD_AUDIT.md` (this file)
- `03_Output/tables/pernambuco_v4_9_replication/spatial_model/PE_U14_likelihood_diagnostics.csv`
- `03_Output/figures/pernambuco_v4_9_replication/spatial0/PE_U14_posterior_predictive.png`
- Supporting: `PE_U14_case_serology_conflict.rds`, `PE_U14_fixed_q_profile.csv`, `PE_U14_audit_draws.rds`

## 10. Status

**Spatial-0 remains UNFROZEN.** No 4-chain production run, no vaccine
projection, and no model change has been made. Awaiting review of this
diagnosis before any further model development.
