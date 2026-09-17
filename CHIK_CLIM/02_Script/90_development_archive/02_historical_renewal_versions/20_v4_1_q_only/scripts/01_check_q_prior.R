# Reproduces the q prior adequacy check (>=1e6 Monte Carlo draws) for this
# experiment's own record. q = p_symptomatic * p_detect_given_symptomatic,
# p_symptomatic ~ Beta(30,28), p_detect_given_symptomatic ~ Beta(20,60).
# Compares against Beta(16.2, 108.7) (the prior actually used in
# renewal_ceara_v4_1_q_only.stan). Do NOT weaken/narrow/truncate the prior
# based on this check or to help HMC converge -- report only.

required_packages <- c("here", "tibble", "readr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(tibble); library(readr) })

q_prior_root <- function() {
  root <- here::here()
  if (file.exists(file.path(root, "CHIK_CLIM.Rproj"))) return(root)
  candidate <- file.path(root, "CHIK_CLIM")
  if (file.exists(file.path(candidate, "CHIK_CLIM.Rproj"))) return(candidate)
  stop("Could not identify the inner CHIK_CLIM project root.")
}

run_check_q_prior <- function() {
  root <- q_prior_root()
  out_dir <- file.path(root, "02_Script", "90_development_archive", "02_historical_renewal_versions", "20_v4_1_q_only")
  n_mc <- 2e6
  set.seed(20260911L)
  q_induced <- rbeta(n_mc, 30, 28) * rbeta(n_mc, 20, 60)
  set.seed(20260911L)
  q_beta_approx <- rbeta(n_mc, 16.2, 108.7)

  probs <- c(.025, .05, .25, .5, .75, .95, .975)
  summarise_dist <- function(x) tibble(
    statistic = c("mean", "sd", paste0("p", probs * 100)),
    value = c(mean(x), sd(x), quantile(x, probs, names = FALSE))
  )
  induced_summary <- summarise_dist(q_induced)
  beta_summary <- summarise_dist(q_beta_approx)
  comparison <- merge(induced_summary, beta_summary, by = "statistic", suffixes = c("_induced_product", "_Beta_16.2_108.7"))
  ks <- suppressWarnings(ks.test(q_induced, q_beta_approx)$statistic)

  write_csv(comparison, file.path(out_dir, "q_prior_check.csv"))
  message("[q-prior] induced product (Beta(30,28) x Beta(20,60), n=", n_mc, ") vs Beta(16.2, 108.7):")
  print(comparison)
  message(sprintf("[q-prior] KS distance = %.4f (adequate if <= 0.02) -> %s", ks, if (ks <= 0.02) "ADEQUATE, using q ~ Beta(16.2, 108.7) as-is" else "NOT adequate -- see note below"))
  if (ks > 0.02) message("[q-prior] NOTE: per instructions this prior is used regardless (do not weaken/narrow/truncate for computational convenience); report the discrepancy rather than silently substituting.")
  invisible(comparison)
}

if (sys.nframe() == 0L) run_check_q_prior()
