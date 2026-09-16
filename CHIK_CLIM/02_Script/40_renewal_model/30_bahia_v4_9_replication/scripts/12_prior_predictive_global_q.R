# Prior-predictive check for logit_q ~ normal(logit(0.10), 1.0), BEFORE
# fitting the global-q model, per the project's standing requirement that
# any new prior is checked prior-predictively before use.

q_prior_draws <- plogis(rnorm(20000, qlogis(0.10), 1.0))
message(sprintf("q prior: median=%.4f | 50%% CrI=[%.4f, %.4f] | 95%% CrI=[%.4f, %.4f]",
                 median(q_prior_draws), quantile(q_prior_draws, .25), quantile(q_prior_draws, .75),
                 quantile(q_prior_draws, .025), quantile(q_prior_draws, .975)))
message(sprintf("Coverage check vs the fixed-q sensitivity grid (0.05-0.30): P(q in [0.05,0.30]) = %.3f",
                 mean(q_prior_draws >= 0.05 & q_prior_draws <= 0.30)))
message(sprintf("P(q < 0.02) = %.3f | P(q > 0.50) = %.3f (checking the prior does not collapse onto an implausible corner)",
                 mean(q_prior_draws < 0.02), mean(q_prior_draws > 0.50)))
