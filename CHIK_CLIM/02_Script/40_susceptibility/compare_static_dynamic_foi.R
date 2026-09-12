# Compare pre-specified M0 (static) and M1 (dynamic AR(1)) Ceará fits.

required_packages <- c("here", "rstan", "dplyr", "readr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(rstan); library(dplyr); library(readr) })
project_root <- function() { root <- here::here(); if (file.exists(file.path(root,"CHIK_CLIM.Rproj"))) return(root); file.path(root,"CHIK_CLIM") }
project_path <- function(...) file.path(project_root(), ...)
bfmi <- function(fit) vapply(rstan::get_sampler_params(fit, inc_warmup=FALSE), function(x) mean(diff(x[,"energy__"])^2)/var(x[,"energy__"]), numeric(1))

compare_static_dynamic_foi <- function() {
  ids <- c("M0_Q1", "M1_Q1")
  paths <- project_path("02_Script","stan",paste0("dynamic_annual_foi_ceara_",ids,".rds"))
  if (any(!file.exists(paths))) stop("M0_Q1 and M1_Q1 fits are required before comparison.")
  rows <- lapply(seq_along(paths), function(i) {
    b <- readRDS(paths[i]); fit <- b$fit; d <- rstan::extract(fit)
    sampler <- rstan::get_sampler_params(fit, inc_warmup=FALSE)
    summary <- summary(fit)$summary
    summary <- summary[!grepl("^(C_pred|log_lik_cases)", rownames(summary)), , drop = FALSE]
    tibble(
      fit_id=ids[i], elapsed_seconds=b$elapsed_seconds,
      divergences=sum(vapply(sampler,function(x)sum(x[,"divergent__"]),numeric(1))),
      treedepth_hits=sum(vapply(sampler,function(x)sum(x[,"treedepth__"]>=12),numeric(1))),
      max_rhat=max(summary[,"Rhat"],na.rm=TRUE), min_n_eff=min(summary[,"n_eff"],na.rm=TRUE), min_bfmi=min(bfmi(fit)),
      q_median=median(d$q), mu_lambda_median=median(d$mu_lambda),
      cumulative_infections_median=median(d$cumulative_latent_infections),
      S_end_2025_median=median(d$S_end_prop[,ncol(d$S_end_prop)]),
      ppc_rmse=sqrt(mean((apply(d$C_pred, 2, median)-b$prepared$annual$observed_cases)^2))
    )
  }) |> bind_rows()
  out <- project_path("03_Output","tables","dynamic_annual_foi","dynamic_annual_foi_M0_M1_comparison.csv")
  write_csv(rows,out)
  rows
}
if (sys.nframe() == 0L) print(compare_static_dynamic_foi())
