# Dynamic annual FOI priors

The Ceará FOI anchor is based on the 100 `foi_equiv` random-forest ensemble
members. It is described as **ensemble uncertainty**, not Bayesian posterior
uncertainty. Its log-scale mean is -4.38704 and SD is 0.16264 (the preparation
script recomputes these directly from the CSV each run).

Overall detection is parameterised directly as `logit(q)`. Q1 is a fixed-seed
moment-matched logit-normal approximation to the product of the existing
renewal priors `p_symp ~ Beta(30,28)` and `rho_sym_Brazil ~ Beta(20,60)`:

| Scenario | logit(q) mean | logit(q) SD |
|---|---:|---:|
| Q1 primary | -1.93304 | 0.26919 |
| Q2 broader | -1.93304 | 0.53838 |
| Q3 lower detection | -2.68969 | 0.26919 |

Q3 centres q at one half of the simulated Q1 median. These priors are used
only as a transparent external anchor for the confirmed-case observation
probability; q is estimated in Stan and remains constant over 2015–2025.

For M1: `rho_delta ~ Normal(0,0.45)` restricted to [-0.95,0.95],
`sigma_delta ~ Normal(0,0.75)` restricted positive, and standard-normal
non-centred innovations. `phi_obs ~ Gamma(2,0.1)`, retaining the regularising
negative-binomial dispersion form used in the previous surveillance models.
Q1's FOI prior SD is also multiplied by 1.5 and 2.0 in planned anchor-scale
sensitivity reruns; these do not replace the required Q1–Q3 detection check.
