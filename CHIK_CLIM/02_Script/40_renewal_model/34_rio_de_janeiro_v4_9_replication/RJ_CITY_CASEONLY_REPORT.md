# RJ CITY — Model A (case-only, global-q) report

## Purpose

Auxiliary geography-matched analysis fitting the frozen v4.9 renewal
structure to **Rio de Janeiro CITY** (muni6 = 330455) weekly case data alone,
with the same broad `logit_q ~ Normal(logit(0.10), 1)` prior family used at
state level, and dispersed inits (q0 = 0.05/0.10/0.15/0.20, NOT centred at
0.010). This model exists to (a) provide a city-scale case-only baseline
that Model B (city + U10 serology) can be compared against, and (b) test
whether Rio CITY's case dynamics alone imply a different q than the RJ
STATE case-only fit.

## Data

- 573 weeks (2015–2025), city population 6.2M–6.8M, muni-level births and
  population genuine; deaths are a disclosed population-share proxy of
  state deaths (see `RJ_CITY_DATA_AUDIT.md`).
- Total city cases = 53,957 (42.7% of the state total of 126,428).

## HMC gate: PASS

| chain | divergent | max_treedepth_hits | mean_stepsize | mean_accept | bfmi |
|---|---|---|---|---|---|
| 1 | 0 | 0 | 0.01575 | 0.9805 | 0.9066 |
| 2 | 0 | 0 | 0.01339 | 0.9851 | 0.9259 |
| 3 | 0 | 0 | 0.01549 | 0.9786 | 0.8558 |
| 4 | 0 | 0 | 0.01793 | 0.9762 | 0.8008 |

Overall: divergences = 0, max_treedepth_hits = 0, max Rhat = 1.0033,
min bulk ESS = 1051.4. **hmc_pass = TRUE**, clean on the first run (adapt_delta
= 0.98, 1200 warmup / 750 sampling, dense_e metric), no rescue needed.

## Identification / scientific results

- **q_city (case-only): median = 0.0150, 95% CrI = [0.0129, 0.0180]**
- cor(logit_q, alpha_R) = **-0.702** — still a q–transmission-level ridge, of
  similar character to (slightly less severe than) the state fit's -0.839,
  consistent with a shorter, more geographically concentrated series
  carrying somewhat more identifying information.
- Model-implied immune fraction over the U10 survey window (2018-07-01 to
  2018-10-31, case-only, i.e. *not* informed by serology): median = 9.57%,
  95% CrI = [7.80%, 11.70%].
- immune_2025 (end of series): median = 48.81%, 95% CrI = [42.65%, 53.47%].
- min S/N (population susceptible, minimum over the series): median = 50.94%.
- max R0(t): median = 2.60.

## City vs. state comparison (case-only)

| quantity | RJ STATE (case-only) | RJ CITY (case-only) |
|---|---|---|
| q median [95% CrI] | 0.0150 [0.0130, 0.0179] | 0.0150 [0.0129, 0.0180] |
| cor(logit_q, alpha_R) | -0.839 | -0.702 |

The city and state case-only q estimates are **essentially identical**
(medians match to 3 decimal places, CrIs almost fully overlapping). This is
the first indication that Rio CITY's case-reporting/ascertainment behaves
similarly to the state as a whole, and that the state-level q is not being
distorted by a city/non-city ascertainment split.

See `RJ_CITY_caseonly_trajectories.png` for the 3-panel case PPC / S-U
immunity / q posterior figure.
