# Pernambuco Spatial-1 (5-strata, sigma_region hierarchical) -- Geometry Diagnostic

**Status: Spatial-1 STOPPED as a computational diagnostic, not a candidate
final posterior. No q, susceptibility, or epidemiological point estimate
from this run is reported as an inferential result.**

## 1. What happened

Spatial-1 (5-strata model with a non-centred, sum-to-zero hierarchical
regional intercept `delta_region[r] = sigma_region * (region_contrast_basis
%*% z_region_contrast)`) was launched as a 4-chain, 1000-warmup +
500-sampling pilot. Chain 1 exhibited a severe, persistent warmup pathology
with no recovery trend:

| Metric | Chain 1 (stopped) | Chains 2-4 |
|---|---|---|
| Iterations completed | 78 / 1500 | 1500 / 1500 (all complete) |
| Step size (median) | 0.000132 | 0.033-0.078 |
| Treedepth-14 hits | 26 (all of chain 1's td14 hits) | 0 |
| n_leapfrog (typical) | up to 16,383 | 63-511 |
| Divergences so far | 5 | 23, 24, 27 |

Projected time to complete Chain 1's warmup at its observed rate: >10 hours.
Per the pre-registered live-warmup stop rule, the run was terminated. Chains
2-4 had already completed all 1500 iterations normally by the time Chain 1
was killed, so their post-warmup draws were preserved (raw `sample_file`
CSVs) and used ONLY for the geometry diagnosis below.

## 2. Is sigma_region the new pathology? NO.

Comparing sigma_region itself across chains:

| | Chain 1 | Chains 2-4 pooled |
|---|---|---|
| sigma_region median | 0.0499 | 0.0285 |
| sigma_region range | [0.015, 0.103] | [0.010, 0.142] |

Chain 1's sigma_region is NOT pinned near zero and its range fully overlaps
the "normal" chains' range. **Sigma_region shows no separation whatsoever
between the pathological and healthy chains** (see
`PE_spatial1_sigma_region_pairs.png`, all four panels: the x-axis position
of Chain 1's black points is indistinguishable from the grey "normal"
cloud). This rules out classification A (hierarchical-scale/funnel
pathology dominated by sigma_region).

## 3. What IS different about Chain 1: the SAME q/depletion/R0 mode seen throughout this project

Plotting logit_q, alpha_global, min(S_prop), and max(R0) against
sigma_region (`PE_spatial1_sigma_region_pairs.png`) shows Chain 1 (black)
occupying an **entirely separate cluster** from the healthy chains'
posterior cloud (grey), on every axis EXCEPT sigma_region:

| Quantity | Chain 1 (stuck mode) | Chains 2-4 (healthy mode) |
|---|---|---|
| logit_q | -3.9 to -2.9 (q approx 0.02-0.05) | -0.6 to 0.5 (q approx 0.35-0.62) |
| alpha_global | 0.05-0.46 | approx -0.01 |
| min S_prop (any stratum/week) | 0.28-0.90 | approx 0.94-0.95 |
| max R0 (any stratum/week) | 2.8-5.3 | approx 1.8-2.1 |

This is, quantitatively and visually, the **exact same low-q /
high-depletion / high-R0 mode** that has produced HMC pathology throughout
this entire Pernambuco modelling effort (homogeneous Models A-D). Chain 1
was initialised at the most extreme of the four dispersed q-starting values
(`q_starts = c(0.02, 0.05, 0.10, 0.20)`, Chain 1 -> 0.02) and appears to
have been captured by this mode before warmup adaptation could pull it
toward the other, better-behaved mode that Chains 2-4 (started at
q=0.05/0.10/0.20) found instead.

Divergent draws from the HEALTHY chains (red points) sit WITHIN the normal
posterior cloud on every axis -- they are scattered, mild warmup noise, not
instances of the low-q/high-depletion mode. This further isolates Chain 1's
problem as mode-specific, not a general property of the model.

## 4. Classification

**FAIL -- unresolved cross-chain branch separation / posterior geometry
(dominant mechanism: B, continuing q/depletion/high-R0 pathology, NOT A,
sigma_region funnel).**

Adding regional susceptible pools and a hierarchical regional intercept did
NOT introduce a new source of geometric difficulty (sigma_region is
innocent -- it overlaps fully between Chain 1 and Chains 2-4, and
treedepth-14 pathology is not concentrated at sigma_region -> 0). Instead,
Chain 1 and Chains 2-4 occupy two DISTINCT, well-separated epidemiological
branches of the SAME posterior:

- **Chain 1 (stopped)**: lower-q (approx 0.02-0.05), lower-susceptibility
  (min S/N down to approx 0.6-0.9 in some strata/weeks), higher-alpha /
  higher-R0 (up to approx 5.3) -- the same low-q/high-depletion/high-R0
  mode that has produced HMC pathology throughout every homogeneous PE
  model attempt (A-D).
- **Chains 2-4 (completed)**: higher-q (approx 0.35-0.62), near-fully
  susceptible (S/N approx 0.94-0.95 throughout), lower-alpha / lower-R0
  (approx 1.8-2.1).

Because Spatial-1 is a single joint model, these are not "4 independent
opinions" -- they are 2 chains apparently stuck in different modes of one
multimodal (or near-non-identifiable) posterior surface that Chain
initialisation (dispersed q-starts 0.02/0.05/0.10/0.20) was able to expose.
**Spatial-1 is therefore recorded as a FAILED candidate model.** Its draws
are diagnostic warmup/posterior STATES only, not usable posterior samples:
per instruction, Spatial-1's epidemiological summaries (q, susceptibility,
etc.) are NOT to be reported as inferential results, Chain 1 is not to be
resumed, and the posterior is NOT to be salvaged by retaining only Chains
2-4 (that would silently discard the low-q branch, which is exactly the
kind of unacknowledged mode-selection this diagnostic exists to prevent).

This directly motivates Spatial-0's branch-connectivity design: rather than
assume either branch is correct, or that avoiding extreme initial values
alone would resolve it, Spatial-0 deliberately starts chains in BOTH
branches (plus an intermediate/central class) to test whether removing
sigma_region/delta_region entirely lets all chains converge to one
posterior region.

## 5. Outputs

Table: `03_Output/tables/pernambuco_v4_9_replication/spatial_model/PE_spatial1_geometry_summary.csv`
(per-iteration diagnostics, all 4 chains), `PE_spatial1_chain_summary.csv`
(by-chain summary).

Figures: `03_Output/figures/pernambuco_v4_9_replication/spatial5/PE_spatial1_sigma_region_geometry.png`
(sigma_region vs iteration/treedepth/leapfrog/stepsize, by chain),
`PE_spatial1_sigma_region_pairs.png` (sigma_region vs logit_q/alpha_global/min-S/max-R0,
divergent and Chain-1 points flagged).

Preserved raw chain CSVs (warmup diagnostics + all draws, including
generated quantities): `outputs/spatial5/chains/chain_{1,2,3,4}.csv`
(Chain 1: 78 rows, warmup only; Chains 2-4: 1500 rows each, complete).

## 6. Next step

Per instruction, Spatial-1 is not modified further. Proceeding to
Spatial-0: the same 5 strata and separate S/U/X pools, but with NO
sigma_region/delta_region (all strata share one common R0(t) trajectory) --
testing whether spatially asynchronous depletion ALONE, without any
regional transmission heterogeneity, is sufficient to resolve the
homogeneous-model pathology. Given the branch separation documented above,
Spatial-0 uses 6 chains started in 3 explicit classes -- LOW (near Chain
1's mode: logit_q=-3.199, alpha=0.185), INTERMEDIATE (central prior
region), HIGH (near Chains 2-4's mode: logit_q=-0.462, alpha=-0.009) -- so
that branch-connectivity (do all chains converge to the same posterior
region?) is tested directly and explicitly, rather than assumed away by a
single central initialisation.
