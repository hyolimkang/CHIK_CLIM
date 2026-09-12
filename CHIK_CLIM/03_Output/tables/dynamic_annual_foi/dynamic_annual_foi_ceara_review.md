# Ceará annual dynamic FOI pilot: scientific review

## Computational gate

**HMC FAIL — DO NOT INTERPRET SUSCEPTIBILITY OR ANNUAL FOI TRAJECTORIES.**

For the primary M1-Q1 run (4 chains, 2,000 warmup and 2,000 post-warmup
iterations per chain, `adapt_delta = 0.999`), the diagnostic gate was not met:

- 1 divergent transition;
- no maximum-treedepth hits;
- maximum R-hat 1.0353;
- minimum bulk ESS 179.7;
- minimum BFMI 0.123 (target at least 0.3).

The M0 static comparator had no divergences and BFMI 0.917, whereas M1 had
BFMI 0.123. This localises the computational difficulty to the dynamic annual
FOI process, rather than the population accounting or constant-q observation
model alone.

## Inputs retained for the pilot

- Years: 2015–2025; no 2014 surveillance year was used.
- Outcome: SINAN `CLASSI_FIN == 13`, confirmed chikungunya.
- Long-term anchor: Ceará `foi_equiv` ensemble median 0.01230; valid-population
  coverage 95.4%.
- Initial state: S(2015)=N(2015), U(2015)=0.
- Detection: a single time-constant overall probability q, with Q1–Q3
  pre-specified sensitivities based on the prior product of symptomatic
  probability and Brazil-wide symptomatic reporting.

## Non-interpretable posterior quantities

The M1-Q1 file contains posterior summaries and figures for debugging only.
For example, its median q was 0.109 and median cumulative infections were
1.88 million, but these are not scientific estimates because the HMC gate
failed. Posterior correlations showed q with cumulative infections = -0.40 and
q with end-2025 S/N = 0.40; this is material scale dependence, although it is
not by itself evidence that the external q anchor is absent.

## Recommendation

**STOP for scientific review.** Do not derive weekly S trajectories or
S-at-outbreak-onset quantities. Before any scientific interpretation, review
whether the annual AR(1) deviation process can be retained with a more stable
identified parameterisation, or whether the pre-specified independent,
strongly regularised M2 comparison should be evaluated. Do not add additional
latent reporting, climate, renewal, R0, or weekly components to address this
failure.
