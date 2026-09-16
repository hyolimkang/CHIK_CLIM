# Rio de Janeiro CITY Auxiliary Analysis -- Data Audit

Auxiliary validation analysis only. Does NOT modify or overwrite the RJ
STATE model (`RJ_STATE_GLOBAL_Q_CASE_ONLY`, HMC PASS, q_state ~ 0.015),
which is preserved unchanged.

## Municipality identity

Rio de Janeiro **CITY** (municipality) = `muni6 = "330455"`, `muni7 =
3304557`, confirmed against `01_Data/ibge_muni_name_lookup.rds`
(`name_muni == "Rio de Janeiro"`, `uf == "RJ"`).

## Weekly case series

Source: `01_Data/chik_dlnm_panel_muni_week_2015_2025.rds` (the SAME
validated national muni-week panel used for the state model and every
other state in this project), filtered to `muni6 == "330455"`, same
2015-01-04 to 2025-12-21 window as the state model.

| Check | Result |
|---|---|
| Weeks | 573 (matches state window exactly) |
| Missing weeks | 0 |
| Duplicate weeks | 0 |
| NA / negative cases | 0 |
| Total cases (2015-2025) | 53,957 (42.7% of the state total, 126,428) |
| Peak week | 2019-05-05 (3,058 cases) -- IDENTICAL peak week to the state series |

The city accounts for the large majority of the state's single dominant
2018-2020 epidemic in both timing and share, consistent with Rio de
Janeiro city being the historical point of introduction and largest
population centre in the state.

## Demographic inputs

| Component | Source | Status |
|---|---|---|
| Population | `01_Data/ibge_pop_muni_year_2015_2025.rds` (municipality-level IBGE estimates), linearly interpolated to weekly (1 July reference), 2014 back-filled with the 2015 value (population series starts 2015; births start late Dec 2014) | **Genuine municipality data** |
| Births | SINASC individual records filtered to `CODMUNRES == "330455"` (script `09e_fetch_riodejaneiro_city_sinasc_weekly_births.R`), same national SINASC cache as every other state/city build in this project | **Genuine municipality data** |
| Deaths | **DISCLOSED PROXY**: no municipality-level annual all-cause-death count is cached anywhere in this project (only UF-level, via `ibge_population_projection_uf_2024revision.rds`). City deaths are obtained by allocating the RJ STATE's annual all-cause deaths proportional to the city's share of state population each year: `deaths_city[y] = deaths_state[y] * pop_city[y] / pop_state[y]`. This is the same population-scaled-allocation convention already used and disclosed elsewhere in this project (e.g. the Pernambuco spatial-model importation allocation). | **Proxy, not genuine city-level data** |

Demographic accounting closure (`N_end - N_start - births + deaths -
reconciliation`): **max abs error = 0** (exact, by construction --
`net_population_reconciliation` absorbs any residual, exactly as in every
other state-level build in this project).

Population range: 6,211,223 to 6,775,561 (2015-2025), consistent with Rio
de Janeiro city's known population (~6.7M).

## Assessment of the death-allocation proxy

Because deaths are allocated by population SHARE (not independently
measured), the city's death rate is implicitly assumed identical to the
state average death rate. This could bias the city's net demographic
reconciliation term if Rio city's true mortality rate differs materially
from the state average (plausible, given urban/rural mortality
differentials), but any such bias is absorbed into
`net_population_reconciliation`, not into births or population (both
genuine), and does not affect case counts. Given deaths are a small
correction relative to the population scale (~1,100/week against a ~6.7M
population, i.e. an annual rate of ~0.85%, in line with Brazil's national
crude death rate), this proxy is judged acceptable for an AUXILIARY
validation analysis. It would need to be revisited before any city-level
result is used for a primary (non-auxiliary) claim.

## Outputs

- `RJ_CITY_DATA_AUDIT.md` (this file)
- `03_Output/tables/rio_de_janeiro_v4_9_replication/city/RJ_CITY_weekly_cases_clean.csv`
- `03_Output/tables/rio_de_janeiro_v4_9_replication/city/RJ_CITY_case_audit.csv`
- `03_Output/tables/rio_de_janeiro_v4_9_replication/city/RJ_CITY_demographic_inputs.csv`
- `03_Output/tables/rio_de_janeiro_v4_9_replication/city/riodejaneiro_city_weekly_input.rds` (audited Stan input)
- `03_Output/figures/rio_de_janeiro_v4_9_replication/city/RJ_CITY_cases_2015_2025.png`
