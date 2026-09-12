# =============================================================================
# Diagnose the primary Ceará annual dynamic FOI fit (M1_Q1).
# HMC gates are reported before any scientific interpretation.
# =============================================================================

required_packages <- c("here", "rstan", "posterior", "dplyr", "tidyr", "readr", "ggplot2", "patchwork", "scales")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(rstan); library(posterior); library(dplyr); library(tidyr); library(readr); library(ggplot2); library(patchwork) })

project_root <- function() { root <- here::here(); if (file.exists(file.path(root, "CHIK_CLIM.Rproj"))) return(root); file.path(root, "CHIK_CLIM") }
project_path <- function(...) file.path(project_root(), ...)
qfun <- function(x, p) unname(stats::quantile(x, p, na.rm = TRUE))
interval <- function(x) c(q025 = qfun(x, .025), median = qfun(x, .5), q975 = qfun(x, .975))

hmc_diagnostics <- function(fit, elapsed_seconds, max_treedepth = 12L) {
  draws_array <- posterior::as_draws_array(rstan::extract(fit, permuted = FALSE))
  draw_summary <- posterior::summarise_draws(draws_array, posterior::rhat, posterior::ess_bulk, posterior::ess_tail) |>
    rename(
      rhat = `posterior::rhat`,
      ess_bulk = `posterior::ess_bulk`,
      ess_tail = `posterior::ess_tail`
    ) |>
    filter(!grepl("^(C_pred|log_lik_cases)", variable))
  sampler <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
  divergences <- sum(vapply(sampler, function(x) sum(x[, "divergent__"]), numeric(1)))
  treedepth_hits <- sum(vapply(sampler, function(x) sum(x[, "treedepth__"] >= max_treedepth), numeric(1)))
  bfmi <- vapply(sampler, function(x) mean(diff(x[, "energy__"])^2) / stats::var(x[, "energy__"]), numeric(1))
  gate <- divergences == 0 && treedepth_hits == 0 && max(draw_summary$rhat, na.rm = TRUE) <= 1.01 && min(draw_summary$ess_bulk, na.rm = TRUE) >= 400 && min(draw_summary$ess_tail, na.rm = TRUE) >= 400 && all(bfmi >= .3)
  list(
    summary = draw_summary,
    headline = tibble(
      hmc_gate = ifelse(gate, "PASS", "FAIL"), divergences = divergences,
      max_treedepth_hits = treedepth_hits, max_rhat = max(draw_summary$rhat, na.rm = TRUE),
      n_rhat_gt_1_01 = sum(draw_summary$rhat > 1.01, na.rm = TRUE),
      min_bulk_ess = min(draw_summary$ess_bulk, na.rm = TRUE), min_tail_ess = min(draw_summary$ess_tail, na.rm = TRUE),
      min_bfmi = min(bfmi), elapsed_seconds = elapsed_seconds
    ),
    bfmi = tibble(chain = seq_along(bfmi), bfmi = bfmi)
  )
}

make_trace_data <- function(fit, parameters) {
  array <- rstan::extract(fit, pars = parameters, permuted = FALSE, inc_warmup = FALSE)
  out <- lapply(seq_along(parameters), function(i) {
    expand_grid(iteration = seq_len(dim(array)[1]), chain = seq_len(dim(array)[2])) |>
      mutate(parameter = parameters[i], value = as.vector(array[, , i]))
  })
  bind_rows(out)
}

diagnose_dynamic_annual_foi <- function(fit_id = "M1_Q1") {
  fit_path <- project_path("02_Script", "stan", paste0("dynamic_annual_foi_ceara_", fit_id, ".rds"))
  if (!file.exists(fit_path)) stop("Fit not found: ", fit_path, ". Run fit_dynamic_annual_foi.R first.")
  bundle <- readRDS(fit_path); fit <- bundle$fit; prepared <- bundle$prepared; annual <- prepared$annual
  table_dir <- project_path("03_Output", "tables", "dynamic_annual_foi")
  figure_dir <- project_path("03_Output", "figures")
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE); dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  dr <- rstan::extract(fit, permuted = TRUE)
  years <- annual$year; index <- seq_along(years)

  annual_summary <- tibble(
    year = years, observed_cases = annual$observed_cases,
    baseline_foi_equiv_mid = prepared$foi_anchor$median,
    lambda_median = apply(dr$lambda, 2, qfun, .5), lambda_q025 = apply(dr$lambda, 2, qfun, .025), lambda_q25 = apply(dr$lambda, 2, qfun, .25), lambda_q75 = apply(dr$lambda, 2, qfun, .75), lambda_q975 = apply(dr$lambda, 2, qfun, .975),
    attack_prob_median = apply(dr$attack_prob, 2, qfun, .5), attack_prob_q025 = apply(dr$attack_prob, 2, qfun, .025), attack_prob_q975 = apply(dr$attack_prob, 2, qfun, .975),
    latent_infections_median = apply(dr$latent_infections, 2, qfun, .5), latent_infections_q025 = apply(dr$latent_infections, 2, qfun, .025), latent_infections_q975 = apply(dr$latent_infections, 2, qfun, .975),
    S_start_prop_median = apply(dr$S_start_prop, 2, qfun, .5), S_start_prop_q025 = apply(dr$S_start_prop, 2, qfun, .025), S_start_prop_q975 = apply(dr$S_start_prop, 2, qfun, .975),
    S_end_prop_median = apply(dr$S_end_prop, 2, qfun, .5), S_end_prop_q025 = apply(dr$S_end_prop, 2, qfun, .025), S_end_prop_q975 = apply(dr$S_end_prop, 2, qfun, .975),
    expected_cases_median = apply(dr$expected_cases, 2, qfun, .5), predicted_cases_median = apply(dr$C_pred, 2, qfun, .5), predicted_cases_q025 = apply(dr$C_pred, 2, qfun, .025), predicted_cases_q975 = apply(dr$C_pred, 2, qfun, .975)
  )
  write_csv(annual_summary, project_path("03_Output", "tables", "dynamic_annual_foi_ceara_summary.csv"))

  hmc <- hmc_diagnostics(fit, bundle$elapsed_seconds)
  write_csv(hmc$headline, file.path(table_dir, "dynamic_annual_foi_hmc_headline.csv"))
  write_csv(hmc$bfmi, file.path(table_dir, "dynamic_annual_foi_bfmi.csv"))
  write_csv(arrange(hmc$summary, desc(rhat)) |> slice_head(n = 30), file.path(table_dir, "dynamic_annual_foi_worst_rhat.csv"))
  write_csv(arrange(hmc$summary, ess_bulk) |> slice_head(n = 30), file.path(table_dir, "dynamic_annual_foi_worst_ess.csv"))

  correlation_data <- cbind(
    q = dr$q, mu_lambda = dr$mu_lambda, cumulative_latent_infections = dr$cumulative_latent_infections,
    S_end_2025 = dr$S_end_prop[, ncol(dr$S_end_prop)], dr$lambda
  )
  colnames(correlation_data)[5:ncol(correlation_data)] <- paste0("lambda_", years)
  correlation_matrix <- stats::cor(correlation_data)
  write_csv(as.data.frame(correlation_matrix) |> mutate(parameter = rownames(correlation_matrix), .before = 1), file.path(table_dir, "dynamic_annual_foi_correlation_matrix.csv"))

  q_prior <- filter(prepared$q_scenarios, scenario == bundle$config$scenario)
  posterior_scalars <- tibble(
    parameter = c("mu_lambda", "q", "phi_obs", "rho_delta", "sigma_delta", "realized_mean_lambda", "realized_mean_to_mu_ratio", "cumulative_latent_infections", "S_end_2025"),
    median = c(qfun(dr$mu_lambda,.5), qfun(dr$q,.5), qfun(dr$phi_obs,.5), qfun(dr$rho_delta,.5), qfun(dr$sigma_delta,.5), qfun(dr$realized_mean_lambda_2015_2025,.5), qfun(dr$realized_mean_to_mu_ratio,.5), qfun(dr$cumulative_latent_infections,.5), qfun(dr$S_end_prop[,ncol(dr$S_end_prop)],.5)),
    q025 = c(qfun(dr$mu_lambda,.025), qfun(dr$q,.025), qfun(dr$phi_obs,.025), qfun(dr$rho_delta,.025), qfun(dr$sigma_delta,.025), qfun(dr$realized_mean_lambda_2015_2025,.025), qfun(dr$realized_mean_to_mu_ratio,.025), qfun(dr$cumulative_latent_infections,.025), qfun(dr$S_end_prop[,ncol(dr$S_end_prop)],.025)),
    q975 = c(qfun(dr$mu_lambda,.975), qfun(dr$q,.975), qfun(dr$phi_obs,.975), qfun(dr$rho_delta,.975), qfun(dr$sigma_delta,.975), qfun(dr$realized_mean_lambda_2015_2025,.975), qfun(dr$realized_mean_to_mu_ratio,.975), qfun(dr$cumulative_latent_infections,.975), qfun(dr$S_end_prop[,ncol(dr$S_end_prop)],.975))
  )
  write_csv(posterior_scalars, file.path(table_dir, "dynamic_annual_foi_primary_posterior_scalars.csv"))

  p_week <- ggplot(prepared$weekly_cases, aes(week_start, cases_confirmed)) + geom_line(linewidth=.25) + labs(title="A. Weekly confirmed cases (context)", x=NULL, y="Cases") + theme_classic()
  p_cases <- ggplot(annual_summary, aes(year, observed_cases)) + geom_col(fill="#0072B2") + labs(title="B. Annual observed confirmed cases", y="Cases") + theme_classic()
  p_lambda <- ggplot(annual_summary, aes(year, lambda_median)) + geom_ribbon(aes(ymin=lambda_q025,ymax=lambda_q975),fill="#56B4E9",alpha=.35) + geom_line() + geom_hline(yintercept=prepared$foi_anchor$median,linetype=2,colour="#D55E00") + scale_y_continuous(trans="log10") + labs(title="C. Annual FOI: posterior and long-term anchor",y="FOI (log scale)") + theme_classic()
  p_attack <- ggplot(annual_summary,aes(year,attack_prob_median))+geom_ribbon(aes(ymin=attack_prob_q025,ymax=attack_prob_q975),fill="#CC79A7",alpha=.35)+geom_line()+labs(title="D. Annual attack probability",y="Probability")+theme_classic()
  p_s <- ggplot(annual_summary,aes(year,S_end_prop_median))+geom_ribbon(aes(ymin=S_end_prop_q025,ymax=S_end_prop_q975),fill="#009E73",alpha=.35)+geom_line()+labs(title="E. Susceptible proportion at annual end",y="S/N")+theme_classic()
  p_ppc <- ggplot(annual_summary,aes(year,observed_cases))+geom_point(size=2)+geom_ribbon(aes(ymin=predicted_cases_q025,ymax=predicted_cases_q975),fill="grey70",alpha=.5)+geom_line(aes(y=predicted_cases_median),colour="#D55E00")+labs(title="F. Annual posterior predictive check",y="Cases")+theme_classic()
  p_q <- ggplot(tibble(q=dr$q),aes(q))+geom_histogram(bins=35,fill="#E69F00")+labs(title="G. Posterior overall detection q",x="q",y="Draws")+theme_classic()
  prior_draws <- tibble(parameter=rep(c("mu_lambda","q","rho_delta","sigma_delta"),each=3000), value=c(exp(rnorm(3000,prepared$foi_anchor$log_mean,prepared$foi_anchor$log_sd)),plogis(rnorm(3000,q_prior$q_logit_prior_mean,q_prior$q_logit_prior_sd)),pmax(-.95,pmin(.95,rnorm(3000,0,.45))),abs(rnorm(3000,0,.75))))
  post_draws <- tibble(parameter=rep(c("mu_lambda","q","rho_delta","sigma_delta"),each=length(dr$q)),value=c(dr$mu_lambda,dr$q,dr$rho_delta,dr$sigma_delta),source="Posterior") |> bind_rows(mutate(prior_draws,source="Prior"))
  p_prior <- ggplot(post_draws,aes(value,fill=source))+geom_density(alpha=.35)+facet_wrap(~parameter,scales="free",ncol=2)+labs(title="H. Prior versus posterior")+theme_classic()+theme(legend.position="bottom")
  ggsave(file.path(figure_dir,"dynamic_annual_foi_ceara_diagnostics.pdf"),(p_week|p_cases)/(p_lambda|p_attack)/(p_s|p_ppc)/(p_q|p_prior),width=250,height=300,units="mm",device=cairo_pdf)

  trace_parameters <- c("mu_lambda","q","phi_obs","rho_delta","sigma_delta",paste0("delta[",c(2,3,8),"]"),paste0("lambda[",c(2,3,8),"]"))
  trace_data <- make_trace_data(fit, trace_parameters)
  trace_plot <- ggplot(trace_data,aes(iteration,value,colour=factor(chain)))+geom_line(linewidth=.18)+facet_wrap(~parameter,scales="free_y",ncol=3)+labs(colour="Chain",title="Ceará annual FOI traceplots")+theme_classic(base_size=8)
  ggsave(file.path(figure_dir,"dynamic_annual_foi_ceara_traceplots.pdf"),trace_plot,width=250,height=220,units="mm",device=cairo_pdf)
  invisible(list(hmc=hmc, annual_summary=annual_summary, correlation_matrix=correlation_matrix))
}

if (sys.nframe() == 0L) diagnose_dynamic_annual_foi()
