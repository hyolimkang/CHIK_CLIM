# Orchestration: source all modular scripts and run the full pipeline
# through the three-arm execution step. Figures/tables/report are run as
# separate later steps once results/three_arm_results.rds exists.

root <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/Documents/GitHub/CHIK_CLIM/CHIK_CLIM"
script_dir <- file.path(root, "02_Script/40_renewal_model/37_ceara_historical_age12_vaccine_counterfactual/scripts")
setwd(root)

source(file.path(script_dir, "00_config.R"))
source(file.path(script_dir, "02_generic_age_vaccine_simulator.R"))
source(file.path(script_dir, "01_load_accepted_inputs.R"))
source(file.path(script_dir, "03_run_no_vaccine_arm.R"))
source(file.path(script_dir, "04_run_full_feedback_vaccine_arm.R"))
source(file.path(script_dir, "05_run_direct_only_diagnostic_arm.R"))
source(file.path(script_dir, "06_execute_three_arms_all_draws.R"))

inputs <- load_accepted_inputs()
execute_three_arms_all_draws(inputs)
message("\n[99_run_all] Three-arm execution complete.")
