# Bahia global-q local-consistency audit -- remaining Section 15 figures.

required_packages <- c("dplyr", "readr", "ggplot2", "tidyr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(dplyr); library(readr); library(ggplot2); library(tidyr) })

root <- ROOT # from 00_project_setup.R (sourced by the stage runner / .Rprofile) -- no scientific change
table_dir <- file.path(root, "03_Output/03_bahia_pipeline/tables/bahia_global_q_local_consistency_audit")
figure_dir <- file.path(root, "03_Output/03_bahia_pipeline/figures/bahia_global_q_local_consistency_audit")
theme_v4 <- theme_classic(base_size = 9) + theme(panel.grid.major.y = element_line(colour = "grey90"))

case_implied_wide <- read_csv(file.path(table_dir, "bahia_serology_case_implied_attack_by_q.csv"), show_col_types = FALSE)
implied_q_tbl <- read_csv(file.path(table_dir, "bahia_serology_implied_q_audit.csv"), show_col_types = FALSE)
posterior_local <- read_csv(file.path(table_dir, "bahia_global_q_posterior_local_attack_summary.csv"), show_col_types = FALSE)

# ---- Figure: observed seroprevalence vs case-implied attack fraction by q --
long_df <- case_implied_wide |>
  pivot_longer(starts_with("p_case_implied_q"), names_to = "q", values_to = "p_case_implied") |>
  mutate(q = as.numeric(gsub("p_case_implied_q", "", q)))

p1 <- ggplot(long_df, aes(q, p_case_implied)) +
  geom_line(colour = "#315A7D") + geom_point(colour = "#315A7D", size = 1.3) +
  geom_hline(aes(yintercept = observed_prevalence), colour = "#D55E00", linetype = 2) +
  facet_wrap(~sero_id, scales = "free_y") +
  labs(title = "Bahia: observed seroprevalence (dashed) vs crude case-implied attack fraction by q",
       subtitle = "Case-implied values use municipality-level cumulative reported incidence aligned to each survey window -- NOT a fitted prevalence",
       x = "q", y = "fraction") +
  theme_v4 + theme(strip.background = element_blank())
ggsave(file.path(figure_dir, "bahia_serology_observed_vs_case_implied_by_q.png"), p1, width = 220, height = 160, units = "mm", dpi = 300)
message("[saved] ", file.path(figure_dir, "bahia_serology_observed_vs_case_implied_by_q.png"))

# ---- Figure: crude serology-implied q with untruncated CI ------------------
p2 <- ggplot(implied_q_tbl, aes(sero_id, q_implied)) +
  geom_pointrange(aes(ymin = pmax(q_implied_low_untruncated, 0), ymax = q_implied_high_untruncated), colour = "#76558F", fatten = 2) +
  geom_hline(yintercept = c(0.05, 0.136, 0.359), linetype = 3, colour = "grey50") +
  annotate("text", x = 0.6, y = c(0.05, 0.136, 0.359), label = c("q=0.05", "q=0.136 (current)", "q=0.359 (sero-off)"),
           hjust = 0, size = 2.5, colour = "grey40") +
  scale_y_continuous(labels = scales::label_number(accuracy = 0.001)) +
  coord_cartesian(clip = "off") +
  labs(title = "Bahia: crude serology-implied q per survey (municipality reported incidence / observed seroprevalence)",
       subtitle = "Error bars = untruncated Wilson-CI-propagated range. All six surveys imply q orders of magnitude BELOW any value in the fixed-q sensitivity grid.",
       x = "survey", y = "q_implied") +
  theme_v4 + theme(plot.margin = margin(5, 60, 5, 5))
ggsave(file.path(figure_dir, "bahia_serology_crude_implied_q.png"), p2, width = 220, height = 130, units = "mm", dpi = 300)
message("[saved] ", file.path(figure_dir, "bahia_serology_crude_implied_q.png"))

# ---- Figure: full q-posterior-propagated case-implied vs observed ----------
p3 <- ggplot(posterior_local, aes(sero_id)) +
  geom_pointrange(aes(y = median, ymin = lo95, ymax = hi95), colour = "#76558F", fatten = 2) +
  geom_point(aes(y = observed_prevalence), colour = "#D55E00", size = 2.5, shape = 18) +
  scale_y_log10(labels = scales::label_percent(accuracy = 0.01)) +
  labs(title = "Bahia: q-posterior-propagated case-implied attack fraction (purple) vs observed seroprevalence (orange diamond)",
       subtitle = "Log scale. Full global-q posterior (median 0.136, 95% CrI 0.070-0.261) propagated through each survey's municipality-level R_survey.",
       x = "survey", y = "fraction (log scale)") +
  theme_v4
ggsave(file.path(figure_dir, "bahia_global_q_posterior_vs_local_serology.png"), p3, width = 200, height = 130, units = "mm", dpi = 300)
message("[saved] ", file.path(figure_dir, "bahia_global_q_posterior_vs_local_serology.png"))
