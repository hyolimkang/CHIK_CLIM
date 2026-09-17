# Pernambuco Spatial-0 -- data builder. Reuses the exact same 5-strata data
# construction as Spatial-1 (script 30: cases, demography, population-scaled
# imports, Recife-direct serology) since NONE of that changed -- only the
# transmission block (no sigma_region/delta_region/region_contrast_basis)
# differs, which is a Stan-side change only. Strips the region_contrast_basis
# field since Spatial-0's Stan file does not declare it.

root <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/Documents/GitHub/CHIK_CLIM/CHIK_CLIM"
source(file.path(root, "02_Script", "04_pernambuco_pipeline", "05_spatial_validation", "11_build_pe_5strata_stan_data.R"))

build_pe_spatial0_stan_data <- function() {
  built <- build_pe_5strata_stan_data()
  built$stan_data$region_contrast_basis <- NULL
  built
}

if (sys.nframe() == 0L) {
  built <- build_pe_spatial0_stan_data()
  message("Spatial-0 data ready: N=", built$stan_data$N, " R=", built$stan_data$R, " Y=", built$stan_data$Y)
}
