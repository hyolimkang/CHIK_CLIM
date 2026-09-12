# Ceará retrospective annual FOI / susceptibility reconstruction

This is a Ceará-only pilot for calendar years 2015–2025. It is not a climate,
renewal, reproduction-number, weekly-transmission, recurrence, or forecasting
model. Its latent outputs are annual FOI `lambda[y]`, annual infections, and
the deterministic annual susceptible trajectory.

The externally aggregated Ceará `foi_equiv` ensemble anchors the long-run
geometric/central FOI scale:

\[
\log(\mu_\lambda) \sim N(\overline{\log(foi\_equiv)}, SD[\log(foi\_equiv)]).
\]

It is a prior, not a constraint on the realised 2015–2025 mean. For M1,
`log(lambda[y]) = log(mu_lambda) + delta[y]`, with a stationary AR(1) prior
for `delta`; no observed-period centring or multiplier normalisation is done.
`mu_lambda` is therefore a geometric/central long-run scale. The generated
quantity `implied_long_run_arithmetic_mean_lambda` reports the lognormal AR(1)
implication separately.

At the start of 2015, `S=N` and `U=0`, reflecting assumed negligible prior
CHIKV immunity before establishment. For each annual period,
`attack_prob = 1-exp(-lambda)` and `infections = S_start * attack_prob`.
Births enter S. Deaths and the exact population reconciliation residual are
allocated between S and U using the post-infection composition, preserving
`S+U=N_end` within numerical tolerance.

The outcome is SINAN confirmed chikungunya (`CLASSI_FIN == 13`), aggregated by
the cleaned event date (onset preferred; notification fallback). The observation
model is `cases[y] ~ NegBinomial2(q * infections[y], phi)`, with one constant
overall confirmed-case detection probability `q` for all years.

M0 fixes `delta=0`. M1 is the primary AR(1) model. Q1–Q3 are pre-specified
detection-prior sensitivities. Results must not be interpreted if the HMC gate
or the q–infection-scale identifiability checks fail.
