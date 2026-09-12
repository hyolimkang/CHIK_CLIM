# FAILED_STAGE_B_ADAPT095_DIVERGENCES -- diagnostic summary

Preserved copy of the Stage B canary (adapt_delta=0.95, warmup=300,
sampling=150, 4 chains), terminated early once post-warmup divergences were
confirmed (strict gate requires divergences=0, already unrecoverable).
Chains 1 and 4 were still in warmup at termination; chains 2 and 3 had
completed all 450 iterations (300 warmup + 150 sampling).

## Post-warmup summary (chains 2 and 3 only -- 1 and 4 did not reach sampling)

| chain | post-warmup divergent | % | treedepth median/p90/p95/p99/max | eq14 hits | n_leapfrog median/p95/max | stepsize (post-adapt) | accept_stat median |
|---|---|---|---|---|---|---|---|
| 2 | 8/150 | 5.3% | 8.5 / 11.1 / 12.0 / 13.0 / 13 | 0 | 495 / 7179 / 14719 | 0.01162 | 0.961 |
| 3 | 25/150 | 16.7% | 7.0 / 9.0 / 10.0 / 11.0 / 11 | 0 | 181 / 1561 / 3583 | 0.01055 | 0.408 |

Chains 1 and 4 (warmup only, not yet at sampling when stopped): treedepth
staying low throughout warmup as well (median ~4, no eq14 hits observed).

## Interpretation

- **Treedepth is dramatically better than v4.1 q-only** (median 7-8.5 here
  vs 10-13 there; max treedepth 11-13 here vs the cap of 14 hit repeatedly
  there). No treedepth==14 hits at all in either completed chain. This
  supports the working hypothesis that the Juazeiro serology anchor
  substantially shortened the q<->alpha_R ridge.
- **But real, recurring post-warmup divergences appeared** -- not a single
  isolated event. Chain 3's divergent iterations are spread throughout the
  post-warmup range (2, 22, 24, 50, 54, 57, 61, 62, 66, 68, 70, 83, 85, 93,
  95, 98, 101, 109, 111, 114, 115, 128, 138, 144, 150 -- of 150), and its
  accept_stat median (0.408) is far below the adapt_delta target (0.95),
  indicating the sampler is persistently struggling in some region, not
  just occasionally clipping a single hard corner.
- Chain 2 is comparatively much healthier (8/150 divergent, accept_stat
  median 0.961 -- close to target) but still fails the divergences=0
  criterion.
- This pattern (clean treedepth + genuine divergences, especially
  concentrated/recurring in one chain) is classically associated with a
  *local* high-curvature region (a "funnel"-like feature) rather than the
  long, flat, hard-to-navigate ridge seen in the pure q-only model. This is
  consistent with the task's own working interpretation: the serology
  anchor shortened the ridge but left a sharp local feature, plausibly
  where the case-data-preferred region and the serology-preferred region
  are reconciled.

## Preserved for the record

- `chains/chain_{1,2,3,4}.csv` -- raw CmdStan-style output, unmodified.
- This run used: v4.0/td14 backbone, q ~ Beta(16.2,108.7), Juazeiro do Norte
  2018 serology (103/404) via beta-binomial(kappa_sero=50), adapt_delta=0.95,
  max_treedepth=14, seed=20260911, q inits at prior 10th/35th/65th/90th
  percentiles by chain.

## Next step

Stage B2: identical scientific model, adapt_delta raised to 0.99 only (per
Section 3 of the task spec), to test whether the divergences are resolvable
by a smaller step size before concluding anything about the model/data
geometry itself.
