# Pernambuco Spatial-1 (5-strata, sigma_region hierarchical) -- Part B.
#
# Spatial-1 was STOPPED as a computational diagnostic (Chain 1: 78/1500
# iterations, stepsize~1e-4, treedepth pegged at 14, ~16000 leapfrogs/iter,
# no recovery trend -- projected >10h to finish warmup). Chains 2-4
# completed all 1500 iterations normally and are used HERE ONLY to
# characterise sampler geometry -- NOT as a final posterior. No q,
# susceptibility, or epidemiological point estimate from this script is an
# inferential result.

required_packages <- c("data.table", "dplyr", "ggplot2", "tidyr", "patchwork")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(data.table); library(dplyr); library(ggplot2); library(tidyr); library(patchwork) })

root <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/Documents/GitHub/CHIK_CLIM/CHIK_CLIM"
pe_root <- file.path(root, "02_Script/40_renewal_model/33_pernambuco_v4_9_replication")
chains_dir <- file.path(pe_root, "outputs/spatial5/chains")
table_dir <- file.path(root, "03_Output/tables/pernambuco_v4_9_replication/spatial_model")
figure_dir <- file.path(root, "03_Output/figures/pernambuco_v4_9_replication/spatial5")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

R_STRATA <- 5

read_chain_diagnostics <- function(chain_id) {
  path <- file.path(chains_dir, sprintf("chain_%d.csv", chain_id))
  raw_head <- readLines(path, n = 200)
  skip_n <- which(!startsWith(raw_head, "#"))[1] - 1L
  header_line <- raw_head[!startsWith(raw_head, "#")][1]
  all_cols <- strsplit(header_line, ",")[[1]]

  # Sampler diagnostics + scalar parameters (cheap).
  scalar_pattern <- "^(lp__|accept_stat__|stepsize__|treedepth__|n_leapfrog__|divergent__|energy__|alpha_global|sigma_region|logit_q|z_region_contrast\\.[0-9]+)$"
  scalar_cols <- all_cols[grepl(scalar_pattern, all_cols)]

  # delta_region (transformed parameter, vector[R]).
  delta_cols <- all_cols[grepl("^delta_region\\.[0-9]+$", all_cols)]

  # fill=TRUE guards against Chain 1's final row being cut off mid-write
  # when the run was stopped; any resulting incomplete row is dropped below.
  read_cols <- unique(c(scalar_cols, delta_cols))
  dt <- fread(path, skip = skip_n, select = read_cols, showProgress = FALSE, fill = TRUE)
  dt <- dt[stats::complete.cases(dt[, ..scalar_cols])]
  dt[, chain := chain_id]
  dt[, iter := .I]
  n_good <- nrow(dt)

  # Min S_prop / max R0 / max R_eff across all (t, r): read these wide
  # matrices separately (still large, but only 3 of them, one at a time).
  for (mat_name in c("S_prop", "R0_t", "R_eff_t")) {
    mat_cols <- all_cols[grepl(paste0("^", mat_name, "\\."), all_cols)]
    if (length(mat_cols) == 0L) next
    mat_dt <- fread(path, skip = skip_n, select = mat_cols, showProgress = FALSE, fill = TRUE)
    mat_dt <- mat_dt[seq_len(n_good)]
    dt[[paste0(mat_name, "_extreme")]] <- if (mat_name == "S_prop") {
      apply(mat_dt, 1, min)
    } else {
      apply(mat_dt, 1, max)
    }
  }

  # Statewide S_PE_prop / U_PE_prop at the FINAL week only (column N).
  n_weeks <- max(as.integer(sub("^S_PE_prop\\.", "", all_cols[grepl("^S_PE_prop\\.", all_cols)])))
  for (statewide_name in c("S_PE_prop", "U_PE_prop")) {
    col <- paste0(statewide_name, ".", n_weeks)
    if (col %in% all_cols) {
      v <- fread(path, skip = skip_n, select = col, showProgress = FALSE, fill = TRUE)[[1]]
      dt[[paste0(statewide_name, "_final")]] <- v[seq_len(n_good)]
    }
  }
  dt
}

message("Reading chain diagnostics (this scans several large matrix blocks per chain)...")
all_chains <- rbindlist(lapply(1:4, read_chain_diagnostics), fill = TRUE)
message("Loaded ", nrow(all_chains), " total rows across 4 chains (chain 1 warmup-only, chains 2-4 complete).")

all_chains[, max_abs_delta_region := apply(.SD, 1, function(x) max(abs(x), na.rm = TRUE)),
           .SDcols = patterns("^delta_region\\.")]
all_chains[, is_divergent := divergent__ == 1]
all_chains[, is_td14 := treedepth__ >= 14]
all_chains[, chain_label := factor(chain, levels = 1:4, labels = c("Chain 1 (STOPPED, warmup-only)", "Chain 2", "Chain 3", "Chain 4"))]

fwrite(all_chains[, .(chain, iter, lp__, accept_stat__, stepsize__, treedepth__, n_leapfrog__, divergent__,
                       alpha_global, sigma_region, logit_q, max_abs_delta_region,
                       S_prop_extreme, R0_t_extreme, R_eff_t_extreme, S_PE_prop_final, U_PE_prop_final)],
       file.path(table_dir, "PE_spatial1_geometry_summary.csv"))
message("[saved] ", file.path(table_dir, "PE_spatial1_geometry_summary.csv"))

# ---- Summary: pathological (divergent / treedepth==14) vs healthy draws ---
summarise_by_flag <- function(df, flag_col, label) {
  df |> group_by(across(all_of(flag_col))) |>
    summarise(n = n(), median_sigma_region = median(sigma_region), median_logit_q = median(logit_q),
              median_alpha_global = median(alpha_global), median_max_abs_delta_region = median(max_abs_delta_region),
              median_min_S = median(S_prop_extreme, na.rm = TRUE), median_max_R0 = median(R0_t_extreme, na.rm = TRUE),
              median_max_Reff = median(R_eff_t_extreme, na.rm = TRUE), .groups = "drop") |>
    mutate(comparison = label)
}
df <- as.data.frame(all_chains)
message("\n=== Divergent vs non-divergent (chains 2-4 pooled; chain 1 has 0 recorded divergences so far) ===")
div_summary <- summarise_by_flag(df, "is_divergent", "divergent_vs_not")
print(as.data.frame(div_summary), digits = 3)

message("\n=== treedepth==14 vs <14 ===")
td_summary <- summarise_by_flag(df, "is_td14", "td14_vs_not")
print(as.data.frame(td_summary), digits = 3)

message("\n=== By chain ===")
chain_summary <- df |> group_by(chain_label) |>
  summarise(n = n(), n_divergent = sum(is_divergent), n_td14 = sum(is_td14),
            median_sigma_region = median(sigma_region), sigma_region_range = paste(round(range(sigma_region), 4), collapse = "-"),
            median_stepsize = median(stepsize__), median_logit_q = median(logit_q), .groups = "drop")
print(as.data.frame(chain_summary), digits = 3)
write.csv(chain_summary, file.path(table_dir, "PE_spatial1_chain_summary.csv"), row.names = FALSE)

# ---- Is Chain 1 concentrated near sigma_region -> 0? ----------------------
message(sprintf("\nChain 1 sigma_region: median=%.5f, range=[%.5f, %.5f]",
                 median(df$sigma_region[df$chain == 1]), min(df$sigma_region[df$chain == 1]), max(df$sigma_region[df$chain == 1])))
message(sprintf("Chains 2-4 pooled sigma_region: median=%.5f, range=[%.5f, %.5f]",
                 median(df$sigma_region[df$chain != 1]), min(df$sigma_region[df$chain != 1]), max(df$sigma_region[df$chain != 1])))

# ---- Figures ----------------------------------------------------------------
theme_v4 <- theme_classic(base_size = 9) + theme(panel.grid.major.y = element_line(colour = "grey90"))
chain_colours <- c("Chain 1 (STOPPED, warmup-only)" = "black", "Chain 2" = "#0072B2", "Chain 3" = "#009E73", "Chain 4" = "#D55E00")

p1 <- ggplot(df, aes(iter, sigma_region, colour = chain_label)) + geom_line(linewidth = 0.3) +
  scale_colour_manual(values = chain_colours, name = NULL) +
  labs(title = "sigma_region vs iteration, by chain", x = "Iteration (incl. warmup)", y = "sigma_region") +
  theme_v4 + theme(legend.position = "top")

p2 <- ggplot(df, aes(treedepth__, sigma_region, colour = chain_label)) +
  geom_jitter(size = 0.5, alpha = 0.3, width = 0.15) +
  scale_colour_manual(values = chain_colours, guide = "none") +
  labs(title = "sigma_region vs treedepth", x = "treedepth", y = "sigma_region") + theme_v4

p3 <- ggplot(df, aes(n_leapfrog__, sigma_region, colour = chain_label)) + geom_point(size = 0.5, alpha = 0.3) +
  scale_x_log10() + scale_colour_manual(values = chain_colours, guide = "none") +
  labs(title = "sigma_region vs leapfrog count", x = "n_leapfrog (log scale)", y = "sigma_region") + theme_v4

p4 <- ggplot(df, aes(stepsize__, sigma_region, colour = chain_label)) + geom_point(size = 0.5, alpha = 0.3) +
  scale_x_log10() + scale_colour_manual(values = chain_colours, guide = "none") +
  labs(title = "sigma_region vs step size", x = "stepsize (log scale)", y = "sigma_region") + theme_v4

combined1 <- p1 / (p2 | p3 | p4)
ggsave(file.path(figure_dir, "PE_spatial1_sigma_region_geometry.png"), combined1, width = 260, height = 220, units = "mm", dpi = 300, bg = "white")
message("\n[saved] ", file.path(figure_dir, "PE_spatial1_sigma_region_geometry.png"))

# ---- Pair plots with divergent / treedepth-14 / chain-1 flagged -----------
df$point_flag <- with(df, ifelse(chain == 1, "Chain 1 (stopped)",
                            ifelse(is_divergent, "Divergent", ifelse(is_td14, "Treedepth=14", "Normal (chains 2-4)"))))
flag_colours <- c("Normal (chains 2-4)" = "grey70", "Treedepth=14" = "#E69F00", "Divergent" = "red", "Chain 1 (stopped)" = "black")
flag_sizes <- c("Normal (chains 2-4)" = 0.3, "Treedepth=14" = 0.8, "Divergent" = 1.2, "Chain 1 (stopped)" = 1.0)
flag_alpha <- c("Normal (chains 2-4)" = 0.15, "Treedepth=14" = 0.6, "Divergent" = 0.9, "Chain 1 (stopped)" = 0.7)

df$point_flag <- factor(df$point_flag, levels = c("Normal (chains 2-4)", "Treedepth=14", "Divergent", "Chain 1 (stopped)"))
df_ordered <- df[order(df$point_flag), ] # plot normal first, flagged points on top

pair_plot <- function(xvar, yvar, xlab, ylab, logx = FALSE) {
  p <- ggplot(df_ordered, aes(.data[[xvar]], .data[[yvar]], colour = point_flag, size = point_flag, alpha = point_flag)) +
    geom_point() +
    scale_colour_manual(values = flag_colours, name = NULL) + scale_size_manual(values = flag_sizes, guide = "none") +
    scale_alpha_manual(values = flag_alpha, guide = "none") +
    labs(x = xlab, y = ylab) + theme_v4
  if (logx) p <- p + scale_x_log10()
  p
}

pp1 <- pair_plot("sigma_region", "logit_q", "sigma_region", "logit_q")
pp2 <- pair_plot("sigma_region", "alpha_global", "sigma_region", "alpha_global")
pp3 <- pair_plot("sigma_region", "S_prop_extreme", "sigma_region", "min S_prop (any t,r)")
pp4 <- pair_plot("sigma_region", "R0_t_extreme", "sigma_region", "max R0 (any t,r)")

combined2 <- (pp1 | pp2) / (pp3 | pp4) + plot_layout(guides = "collect") & theme(legend.position = "top")
ggsave(file.path(figure_dir, "PE_spatial1_sigma_region_pairs.png"), combined2, width = 240, height = 220, units = "mm", dpi = 300, bg = "white")
message("[saved] ", file.path(figure_dir, "PE_spatial1_sigma_region_pairs.png"))
