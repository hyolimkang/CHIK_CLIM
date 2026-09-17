# Pernambuco Spatial-0 -- Sections 9-10: PRIMARY gate is cross-chain branch
# agreement (do LOW/INTERMEDIATE/HIGH-initialised chains converge to the
# SAME posterior region?), evaluated BEFORE the standard HMC gate. Reads
# directly from the raw per-chain sample_file CSVs so it works whether the
# run completed normally or was stopped early (same robust approach as the
# Spatial-1 geometry diagnostic, script 35).

required_packages <- c("data.table", "dplyr", "ggplot2", "patchwork")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(data.table); library(dplyr); library(ggplot2); library(patchwork) })

root <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/Documents/GitHub/CHIK_CLIM/CHIK_CLIM"
pe_root <- file.path(root, "03_Output/model_fits/pernambuco/v4_9_replication")
chains_dir <- file.path(pe_root, "outputs/spatial0/chains")
table_dir <- file.path(root, "03_Output/tables/pernambuco_v4_9_replication/spatial_model")
figure_dir <- file.path(root, "03_Output/figures/pernambuco_v4_9_replication/spatial0")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

CHAIN_BRANCH <- c("LOW", "LOW", "INTERMEDIATE", "INTERMEDIATE", "HIGH", "HIGH")
WARMUP_N <- 1000L # rows 1..WARMUP_N in each chain's CSV are warmup

read_chain <- function(chain_id) {
  path <- file.path(chains_dir, sprintf("chain_%d.csv", chain_id))
  if (!file.exists(path)) return(NULL)
  raw_head <- readLines(path, n = 200)
  skip_n <- which(!startsWith(raw_head, "#"))[1] - 1L
  header_line <- raw_head[!startsWith(raw_head, "#")][1]
  all_cols <- strsplit(header_line, ",")[[1]]

  scalar_cols <- all_cols[grepl("^(lp__|accept_stat__|stepsize__|treedepth__|n_leapfrog__|divergent__|energy__|alpha_R|logit_q)$", all_cols)]
  dt <- fread(path, skip = skip_n, select = scalar_cols, showProgress = FALSE, fill = TRUE)
  dt <- dt[stats::complete.cases(dt)]
  n_good <- nrow(dt)
  dt[, chain := chain_id]
  dt[, iter := .I]
  dt[, branch := CHAIN_BRANCH[chain_id]]
  dt[, phase := ifelse(iter <= WARMUP_N, "warmup", "sampling")]

  for (mat_name in c("S_prop", "R0_t")) {
    mat_cols <- all_cols[grepl(paste0("^", mat_name, "\\."), all_cols)]
    if (length(mat_cols) == 0L) next
    mat_dt <- fread(path, skip = skip_n, select = mat_cols, showProgress = FALSE, fill = TRUE)[seq_len(n_good)]
    dt[[paste0(mat_name, "_extreme")]] <- if (mat_name == "S_prop") apply(mat_dt, 1, min) else apply(mat_dt, 1, max)
  }
  n_weeks <- max(as.integer(sub("^U_PE_prop\\.", "", all_cols[grepl("^U_PE_prop\\.", all_cols)])))
  col <- paste0("U_PE_prop.", n_weeks)
  if (col %in% all_cols) {
    dt$immune_2025 <- fread(path, skip = skip_n, select = col, showProgress = FALSE, fill = TRUE)[[1]][seq_len(n_good)]
  }
  dt
}

message("Reading Spatial-0 chain diagnostics...")
chain_files <- list.files(chains_dir, pattern = "^chain_[0-9]+\\.csv$")
chain_ids <- as.integer(gsub("[^0-9]", "", chain_files))
all_chains <- rbindlist(lapply(sort(chain_ids), read_chain), fill = TRUE)
all_chains[, q := plogis(logit_q)]
message("Loaded ", nrow(all_chains), " rows across chains: ", paste(sort(chain_ids), collapse = ", "))

by_chain_status <- all_chains[, .(n_rows = .N, n_warmup = sum(phase == "warmup"), n_sampling = sum(phase == "sampling"),
                                   branch = branch[1]), by = chain]
message("\n=== Per-chain row counts (branch) ===")
print(as.data.frame(by_chain_status))

# ---- Section 9: PRIMARY gate -- cross-chain branch agreement (post-warmup only, if any) ----
post_warmup <- all_chains[phase == "sampling"]
if (nrow(post_warmup) == 0L) {
  message("\n*** No chain has reached post-warmup sampling yet. Branch-agreement check uses LATE-WARMUP draws (last 200 rows per chain) as a provisional signal only. ***")
  post_warmup <- all_chains[, .SD[max(1, .N - 199):.N], by = chain]
}

branch_summary <- post_warmup[, .(
  n = .N, median_q = median(q), lo95_q = quantile(q, .025), hi95_q = quantile(q, .975),
  median_alpha_R = median(alpha_R), median_min_S = median(S_prop_extreme, na.rm = TRUE),
  median_max_R0 = median(R0_t_extreme, na.rm = TRUE), median_immune_2025 = median(immune_2025, na.rm = TRUE)
), by = .(chain, branch)]
message("\n=== Section 9: per-chain posterior summary (post-warmup, or late-warmup if unavailable) ===")
print(as.data.frame(branch_summary), digits = 3)
write.csv(branch_summary, file.path(table_dir, "PE_spatial0_branch_agreement.csv"), row.names = FALSE)
message("[saved] ", file.path(table_dir, "PE_spatial0_branch_agreement.csv"))

branch_ranges <- post_warmup[, .(q_lo = quantile(q, .025), q_hi = quantile(q, .975)), by = branch]
message("\n=== Branch q ranges (95% of post-warmup/late-warmup draws) ===")
print(as.data.frame(branch_ranges))

# Overlap check: do the three branches' 95% q intervals overlap?
ranges_list <- split(branch_ranges, branch_ranges$branch)
overlap <- function(a, b) a$q_lo <= b$q_hi && b$q_lo <= a$q_hi
if (all(c("LOW", "INTERMEDIATE", "HIGH") %in% names(ranges_list))) {
  ov_low_mid <- overlap(ranges_list$LOW, ranges_list$INTERMEDIATE)
  ov_mid_high <- overlap(ranges_list$INTERMEDIATE, ranges_list$HIGH)
  ov_low_high <- overlap(ranges_list$LOW, ranges_list$HIGH)
  message(sprintf("\nq-interval overlap: LOW-INTERMEDIATE=%s, INTERMEDIATE-HIGH=%s, LOW-HIGH=%s",
                   ov_low_mid, ov_mid_high, ov_low_high))
  branch_agreement_pass <- ov_low_mid && ov_mid_high
} else {
  message("\nNot all three branches have post-warmup/late-warmup draws yet -- branch agreement not yet assessable.")
  branch_agreement_pass <- NA
}

# ---- Figure: q trajectory by chain/branch over the whole run ----
theme_v4 <- theme_classic(base_size = 9) + theme(panel.grid.major.y = element_line(colour = "grey90"))
branch_colours <- c("LOW" = "#D55E00", "INTERMEDIATE" = "#009E73", "HIGH" = "#0072B2")
p_q_traj <- ggplot(all_chains, aes(iter, q, colour = branch, group = chain)) +
  geom_line(linewidth = 0.3, alpha = 0.8) +
  geom_vline(xintercept = WARMUP_N, linetype = "dashed", colour = "grey50") +
  scale_colour_manual(values = branch_colours, name = "Init branch") +
  labs(title = "Pernambuco Spatial-0: q trajectory by chain (dashed line = end of warmup)",
       subtitle = "Do LOW/INTERMEDIATE/HIGH-initialised chains converge to the same q?", x = "Iteration", y = "q") +
  theme_v4

p_q_density <- ggplot(post_warmup, aes(q, fill = branch)) + geom_density(alpha = 0.5) +
  scale_fill_manual(values = branch_colours, name = "Init branch") +
  labs(title = "q posterior density by branch (post-warmup / late-warmup)", x = "q", y = "Density") + theme_v4

ggsave(file.path(figure_dir, "PE_spatial0_branch_agreement.png"), p_q_traj / p_q_density,
       width = 220, height = 200, units = "mm", dpi = 300, bg = "white")
message("[saved] ", file.path(figure_dir, "PE_spatial0_branch_agreement.png"))

# ---- Section 10: standard HMC gate (only meaningful if branch agreement passes) ----
hmc_by_chain <- all_chains[phase == "sampling", .(
  n_sampling = .N, divergent = sum(divergent__), max_treedepth_hits = sum(treedepth__ >= 14),
  mean_stepsize = mean(stepsize__), mean_leapfrog = mean(n_leapfrog__), mean_accept = mean(accept_stat__)
), by = .(chain, branch)]
message("\n=== Section 10: HMC diagnostics by chain (post-warmup only) ===")
print(as.data.frame(hmc_by_chain), digits = 4)
write.csv(hmc_by_chain, file.path(table_dir, "PE_spatial0_HMC.csv"), row.names = FALSE)
message("[saved] ", file.path(table_dir, "PE_spatial0_HMC.csv"))

message("\n=== VERDICT ===")
if (is.na(branch_agreement_pass)) {
  message("Branch agreement: NOT YET ASSESSABLE (run still in progress or incomplete).")
} else if (!branch_agreement_pass) {
  message("Branch agreement: FAIL -- LOW/INTERMEDIATE/HIGH branches remain separated in q.")
  message("SPATIAL-0 CLASSIFICATION: FAIL -- multimodality / unresolved cross-chain branch separation.")
} else {
  total_divergent <- sum(hmc_by_chain$divergent)
  total_td14 <- sum(hmc_by_chain$max_treedepth_hits)
  message(sprintf("Branch agreement: PASS. Total post-warmup divergences=%d, treedepth-14 hits=%d.", total_divergent, total_td14))
  if (total_divergent == 0 && total_td14 == 0) {
    message("SPATIAL-0 CLASSIFICATION (provisional): PASS on branch agreement + standard HMC gate -- proceed to identification checks (Section 11).")
  } else {
    message("SPATIAL-0 CLASSIFICATION: FAIL -- HMC geometry (branches agree, but standard HMC gate not clean).")
  }
}
