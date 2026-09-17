# FAILED_STAGE_B2_ADAPT099_TREEDEPTH -- diagnostic summary

Preserved copy of the Stage B2 canary (adapt_delta=0.99, diag_e metric,
warmup=400, sampling=200, 4 chains, q inits at prior 20th/40th/60th/80th
percentiles, transmission-parameter inits drawn from 4 distinct v4.0/td14
posterior draws). Terminated by user instruction once live warmup
diagnostics already showed the strict gate could not be satisfied --
**no chain had completed warmup** when stopped (all still within the first
17% of the 400-iteration warmup phase), so this is a warmup-phase-only
diagnostic; no post-warmup (sampling-phase) rows exist.

## Warmup-so-far summary (all 4 chains, before Stage-A-terminated point)

| chain | rows so far | treedepth median | treedepth max | eq14 hits | n_leapfrog median | n_leapfrog max | stepsize (last) |
|---|---|---|---|---|---|---|---|
| 1 | 63 | 12.0 | 14 | 8 (12.7%) | 5477 | 16383 | 0.000047 |
| 2 | 47 | 13.0 | 14 | 18 (38.3%) | 8191 | 16383 | 0.000067 |
| 3 | 68 | 12.0 | 14 | 10 (14.7%) | 4095 | 16383 | 0.000039 |
| 4 | 66 | 12.0 | 14 | 11 (16.7%) | 4095 | 16383 | 0.000052 |

## Interpretation

Raising adapt_delta from 0.95 to 0.99 made treedepth **worse**, not better
-- median treedepth roughly doubled (7-8.5 in Stage B's completed chains vs
12-13 here), and treedepth==14 hits are already occurring in 13-38% of
warmup iterations (vs 0 in Stage B's post-warmup chains). Stepsize
collapsed to ~0.00004-0.00007 (roughly 100-250x smaller than Stage B's
~0.01), consistent with the classic adapt_delta/treedepth trade-off: a
smaller step size (fewer/no divergences, in principle) requires far more
leapfrog steps to traverse the same trajectory length, which is expensive
in exactly the region where Stage B was already struggling. This confirms
Stage B2 was correctly judged un-viable without waiting for it to reach
sampling.

## Next step

Per instruction: do NOT increase adapt_delta further, do NOT increase
max_treedepth. Instead, isolate the mass-matrix metric as a single factor:
return to adapt_delta=0.95 (Stage B's original value) and use `dense_e`
instead of the default `diag_e` metric, to test whether the difficulty is
driven by strong posterior correlations a diagonal metric cannot represent.
