# Convergence assessment

## What passed

- Divergences: 0.
- Maximum R-hat: 1.0099.
- Minimum bulk ESS: 424.964.
- BFMI was above the prespecified 0.30 threshold in every chain.

## Why this run is not accepted

- Maximum-treedepth hits: 64 (1.60% of retained transitions).
- Hits by chain: 51, 0, 5, 8.
- Tree-depth quantiles: median 10, 95th percentile 11, 99th percentile 12.

## Interpretation

This is not a catastrophic failure: there are no divergences, R-hat is below 1.01, and effective sample sizes are adequate. However, the sampler repeatedly reached its numerical trajectory cap, concentrated in one chain. Posterior intervals and the six-panel figure may be used for exploratory model development only; they are not accepted as the baseline for vaccine-impact inference.

## Prespecified next steps

1. Run the identical model with a larger maximum treedepth in a separate output folder. If all other diagnostics remain satisfactory and no trajectories hit the new cap, retain that numerical run.
2. If cap hits remain, retain the epidemiological assumptions but reparameterize the annual effects as an explicit sum-to-zero basis, then reassess HMC. This removes a redundant annual-effect direction without adding biological complexity.
3. If reparameterization still fails, reduce annual transmission flexibility before adding climate, reporting, serology, age, or vaccination components. Do not solve this by loosening convergence thresholds.
