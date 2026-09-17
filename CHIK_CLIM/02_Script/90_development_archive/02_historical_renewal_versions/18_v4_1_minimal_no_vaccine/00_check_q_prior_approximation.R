# Monte Carlo check of the informative prior for the single constant
# reporting fraction q in v4.1.
#
# External epidemiological information: p_symptomatic ~ Beta(30,28),
# p_detect_given_symptomatic ~ Beta(20,60). The case likelihood identifies
# mainly their product, so v4.1 estimates ONE constant q with a direct
# informative prior approximating the induced distribution of
# Beta(30,28) * Beta(20,60), rather than the two factors separately.
# This script draws >=1e6 Monte Carlo samples of that product, reports its
# moments/quantiles, compares them against the previously suggested
# Beta(16.2, 108.7) approximation, and -- if that approximation is not
# adequate -- moment-matches a better Beta prior and documents it.

set.seed(20260911L)
n_mc <- 2e6

p_symptomatic <- rbeta(n_mc, 30, 28)
p_detect <- rbeta(n_mc, 20, 60)
q_induced <- p_symptomatic * p_detect

summarise_dist <- function(x) {
  c(mean = mean(x), sd = sd(x), median = median(x),
    q025 = quantile(x, .025, names = FALSE), q975 = quantile(x, .975, names = FALSE))
}

induced_summary <- summarise_dist(q_induced)

# Moment-match a Beta(a,b) to the INDUCED distribution's own mean/variance
# (method of moments), for comparison against the previously suggested
# Beta(16.2, 108.7).
mean_induced <- induced_summary["mean"]
var_induced <- induced_summary["sd"]^2
common <- mean_induced * (1 - mean_induced) / var_induced - 1
a_matched <- mean_induced * common
b_matched <- (1 - mean_induced) * common

candidates <- list(
  previously_suggested = c(a = 16.2, b = 108.7),
  moment_matched = c(a = unname(a_matched), b = unname(b_matched))
)

set.seed(20260911L)
candidate_draws <- lapply(candidates, function(ab) rbeta(n_mc, ab["a"], ab["b"]))
candidate_summaries <- lapply(candidate_draws, summarise_dist)

report <- rbind(
  induced_product = induced_summary,
  previously_suggested_Beta_16.2_108.7 = candidate_summaries$previously_suggested,
  moment_matched_Beta = candidate_summaries$moment_matched
)
report <- as.data.frame(report)
report$parameters <- c(
  "p_symptomatic ~ Beta(30,28) * p_detect ~ Beta(20,60) [Monte Carlo, n=2e6]",
  sprintf("Beta(%.4g, %.4g)", candidates$previously_suggested["a"], candidates$previously_suggested["b"]),
  sprintf("Beta(%.4g, %.4g)", candidates$moment_matched["a"], candidates$moment_matched["b"])
)

# KS distance as a distributional-adequacy check (not just moment match).
ks_previously_suggested <- suppressWarnings(ks.test(q_induced, candidate_draws$previously_suggested)$statistic)
ks_moment_matched <- suppressWarnings(ks.test(q_induced, candidate_draws$moment_matched)$statistic)

cat("=== q prior Monte Carlo check (n =", n_mc, ") ===\n\n")
print(report, digits = 4)
cat("\nKS distance (induced product vs previously-suggested Beta(16.2,108.7)):", signif(ks_previously_suggested, 4), "\n")
cat("KS distance (induced product vs moment-matched Beta):", signif(ks_moment_matched, 4), "\n")

adequate_threshold <- 0.02
decision <- if (ks_previously_suggested <= adequate_threshold) {
  "ADEQUATE: use the previously suggested Beta(16.2, 108.7) as q's prior in v4.1."
} else {
  sprintf("NOT adequate (KS=%.4f > %.2f threshold): use the moment-matched Beta(%.4g, %.4g) as q's prior in v4.1 instead.",
          ks_previously_suggested, adequate_threshold, a_matched, b_matched)
}
cat("\nDECISION:", decision, "\n")

root <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/Documents/GitHub/CHIK_CLIM/CHIK_CLIM"
out_dir <- file.path(root, "03_Output", "tables", "renewal_v4_1_minimal_no_vaccine")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
report_out <- report
report_out$ks_vs_induced <- c(NA, ks_previously_suggested, ks_moment_matched)
write.csv(report_out, file.path(out_dir, "q_prior_monte_carlo_check.csv"), row.names = TRUE)
writeLines(decision, file.path(out_dir, "q_prior_decision.txt"))
cat("\nSaved:", file.path(out_dir, "q_prior_monte_carlo_check.csv"), "\n")
