# v3.2 Coarse-Regional Mapping — Validation Report

**Status: COMPLETE for the mapping/final_model_unit stage.** The 184-row
municipality → health-macroregion crosswalk and the 8-unit `final_model_unit`
derivation have been built and validated below. Per instructions, the Stan
model has not been touched and no HMC has been run. Regional-panel
aggregation, reconciliation, conservation checks, and serology-mapping
(Sections 5–10 of the original task spec) have **not** been attempted yet —
see "Next step / open blocker" at the end of this document.

Outputs produced:

- `01_Data/reference/ceara_municipality_health_macroregion_2018.csv` (184 rows)
- `01_Data/reference/ceara_municipality_final_model_unit_2018.csv` (184 rows, adds `final_model_unit`)

---

## 1. Source used for the official macroregion classification

The repository contains no authoritative Ceará health-macroregion mapping
(confirmed by exhaustive repo search — see prior version of this report /
conversation history). Per instruction, an external authoritative source was
required. This is what was found and used, and what could **not** be
directly retrieved:

### What could not be retrieved directly

`saude.ce.gov.br` (the SESA-CE domain hosting the PDR document and the
resolution archive, including the specific **Resolução CIB/CE nº 03/2018**
"Alteração do PDR — com atualização populacional 2017") returned either an
expired-TLS-certificate error or an F5 WAF 403 ("The requested URL was
rejected") for every access path attempted, including via a text-extraction
proxy. `cosemsce.org.br` (which mirrors CIB resolutions) was blocked by
Cloudflare. **The primary legal resolution text itself was not read.**

### What was used instead (retrieved and cross-validated)

1. **Municipality → 22 health-region assignment**: the DATASUS/IBGE
   municipality-to-health-region territorial table (2019 vintage;
   `regiao_id, estado_abrev, regional_id, regional_nome, municipio_id,
   municipio_id_sdv, municipio` columns), obtained via the public GitHub
   mirror `github.com/lansaviniec/shapefile_das_regionais_de_saude_sus`
   (`DTB.csv`), which states it integrates official IBGE + DATASUS territorial
   sources. Filtered to Ceará (`estado_abrev == "CE"`): exactly **184 rows,
   184 unique municipalities**, grouped into the 22 named health regions
   ("1ª Região Fortaleza" … "22ª Região Cascavel").

2. **22-region → 5-macroregion grouping**: Ceará's own official
   **Plano Estadual de Saúde (PES) 2020–2023**, retrieved via its CONASS
   mirror (`conass.org.br/wp-content/uploads/2022/02/PES-20-23-Atualizado.pdf`,
   279 pages, text-extracted with `pdftools`). Page 159–160 states verbatim:

   > "O Estado está dividido em 05 (cinco) Regiões de Saúde... 1ª Região de
   > Saúde de Fortaleza com 05 áreas descentralizadas Caucaia, Maracanaú,
   > Baturité, Itapipoca e Cascavel; 2ª Região de Saúde de Sobral com 04
   > áreas descentralizadas Acaraú, Tianguá, Crateús e Camocim; 3ª Região de
   > Saúde de Cariri com 04 áreas descentralizadas Icó, Iguatu, Brejo Santo e
   > Crato; 4ª Região de Saúde de Sertão Central com 02 áreas descentralizadas
   > Canindé e Tauá; e a 5ª Região de Saúde de Litoral Leste Jaguaribe com 02
   > áreas descentralizadas Aracati e Russas."

   This is the same 5-macroregion / 22-microregion administrative structure
   named in the task (Fortaleza, Sobral, Cariri, Sertão Central, Litoral
   Leste/Jaguaribe), from Ceará's own official state health-planning
   document, not a third-party summary.

3. **Independent corroboration** (used to cross-check, not as primary
   input): a Ministry of Health / Hospital Alemão Oswaldo Cruz technical note
   on the Cariri macroregion (states 45 municipalities in 5 microregions:
   Icó, Iguatu, Brejo Santo, Crato, Juazeiro do Norte with Barbalha); a
   SESA-CE/SEADE Sertão Central regional profile giving the exact 20
   municipalities under ADS Canindé/Quixadá/Tauá; and a news report on the
   Litoral Leste/Jaguaribe macroregion's creation (2014) giving the exact 20
   constituent municipalities.

### Cross-validation performed

- Every one of the 3 independently-sourced municipality lists (Fortaleza
  health region: 4 munis; Caucaia health region: 10 munis; Acaraú health
  region: 7 munis; Tianguá health region: 8 munis; Sertão Central's 20
  munis; Litoral Leste/Jaguaribe's 20 munis) was checked **name-for-name**
  against the DTB.csv assignment for the corresponding `regional_id` — **100%
  match, zero discrepancies**.
- The macroregion municipality totals computed from DTB.csv using the PES
  page-160 grouping (Fortaleza=44, Sobral=55, Cariri=45, Sertão
  Central=20, Litoral Leste/Jaguaribe=20; sum=184) match the PES text's own
  stated total for the Fortaleza macroregion ("44 municípios... 53% do
  Estado") and the LAPUR/UFC academic source's stated total for Cariri ("45
  municípios... 5 microrregiões") exactly.
- `DTB.csv`'s 184 Ceará `muni6` values and municipality names were checked
  against the repository's existing `01_Data/ibge_muni_name_lookup.rds`
  (filtered to `uf == "CE"`) and against the distinct `muni6` values in
  `01_Data/chik_dlnm_panel_muni_week_2015_2025.rds`: **zero missing on
  either side, zero name mismatches**.

**Caveat, stated plainly:** this crosswalk was not built by reading
Resolução CIB/CE nº 03/2018 itself — that document could not be retrieved
because the hosting site blocked all access attempts. It was built from (a)
an independently-maintained but source-attributed replication of the
national DATASUS/IBGE health-region territorial table, which is the same
data source SINAN's own `ID_REGIONA`/`ID_RG_RESI` fields draw on, and (b)
Ceará's own current official State Health Plan describing the identical
5-macroregion/22-region structure by name, with matching totals. The
resulting structure is internally consistent, matches every independent
academic/administrative source found, and reproduces the exact macroregion
names and pole-city composition given in the task. If bit-exact provenance
to the specific 2018 resolution number is required (e.g. for a publication
methods section), the underlying PDF should be obtained directly from SESA-CE
through a channel that isn't subject to the WAF/TLS blocks encountered here.

---

## 2. Exact original macroregion/region labels (verbatim from source)

`health_region_22` (verbatim `regional_nome` from the DATASUS/IBGE table):

| id | health_region_22 | n_municipalities | official_macroregion_5 |
|----|-------------------|------------------:|--------------------------|
| 23001 | 1ª Região Fortaleza | 4 | Fortaleza |
| 23002 | 2ª Região Caucaia | 10 | Fortaleza |
| 23003 | 3ª Região Maracanaú | 8 | Fortaleza |
| 23004 | 4ª Região Baturité | 8 | Fortaleza |
| 23006 | 6ª Região Itapipoca | 7 | Fortaleza |
| 23022 | 22ª Região Cascavel | 7 | Fortaleza |
| 23011 | 11ª Região Sobral | 24 | Norte_Sobral |
| 23012 | 12ª Região Acaraú | 7 | Norte_Sobral |
| 23013 | 13ª Região Tianguá | 8 | Norte_Sobral |
| 23015 | 15ª Região Crateús | 11 | Norte_Sobral |
| 23016 | 16ª Região Camocim | 5 | Norte_Sobral |
| 23017 | 17ª Região Icó | 7 | Cariri |
| 23018 | 18ª Região Iguatú | 10 | Cariri |
| 23019 | 19ª Região Brejo Santo | 9 | Cariri |
| 23020 | 20ª Região Crato | 13 | Cariri |
| 23021 | 21ª Região Juazeiro do Norte | 6 | Cariri |
| 23005 | 5ª Região Canindé | 6 | Sertao_Central |
| 23008 | 8ª Região Quixadá | 10 | Sertao_Central |
| 23014 | 14ª Região Tauá | 4 | Sertao_Central |
| 23007 | 7ª Região Aracati | 4 | Litoral_Leste_Jaguaribe |
| 23009 | 9ª Região Russas | 5 | Litoral_Leste_Jaguaribe |
| 23010 | 10ª Região Limoeiro do Norte | 11 | Litoral_Leste_Jaguaribe |

## 3. Normalization rules applied

`official_macroregion_5` values are the ASCII-normalized macroregion names,
chosen to match the vocabulary already used by the `final_model_unit`
derivation rules given in the task:

| Official macroregion (PES 2020-2023 wording) | official_macroregion_5 |
|---|---|
| Região de Saúde de Fortaleza | `Fortaleza` |
| Região de Saúde de Sobral | `Norte_Sobral` |
| Região de Saúde de Cariri | `Cariri` |
| Região de Saúde de Sertão Central | `Sertao_Central` |
| Região de Saúde de Litoral Leste Jaguaribe | `Litoral_Leste_Jaguaribe` |

("Sobral" is normalized to `Norte_Sobral` — the SESA structure and the PES
text call it both "Sobral" and, in the Superintendências listing, "Norte";
`Norte_Sobral` is used so the value matches the `final_model_unit` rule set,
which checks for that exact string.)

`health_region_22` is kept **verbatim** as `regional_nome` from the source
table (no normalization) per the instruction to preserve the original label.

## 4. Municipality counts per official macroregion

| official_macroregion_5 | n_municipalities |
|---|---:|
| Fortaleza | 44 |
| Norte_Sobral | 55 |
| Cariri | 45 |
| Sertao_Central | 20 |
| Litoral_Leste_Jaguaribe | 20 |
| **Total** | **184** |

## 5. Special-municipality membership validation (performed before deriving final_model_unit)

| Municipality | muni6 | health_region_22 | official_macroregion_5 | Required | Result |
|---|---|---|---|---|---|
| Fortaleza | 230440 | 1ª Região Fortaleza | Fortaleza | Fortaleza | **PASS** |
| Juazeiro do Norte | 230730 | 21ª Região Juazeiro do Norte | Cariri | Cariri | **PASS** |
| Quixadá | 231130 | 8ª Região Quixadá | Sertao_Central | Sertao_Central | **PASS** |

All three checks passed; no discrepancy required stopping before
`final_model_unit` construction.

## 6. final_model_unit counts (8 modelling units)

| final_model_unit | n_municipalities |
|---|---:|
| Fortaleza | 1 |
| Fortaleza_remainder | 43 |
| Juazeiro_do_Norte | 1 |
| Cariri_remainder | 44 |
| Quixada | 1 |
| Sertao_Central_remainder | 19 |
| Norte_Sobral | 55 |
| Litoral_Leste_Jaguaribe | 20 |
| **Total** | **184** |

(44 − 1 = 43 Fortaleza_remainder; 45 − 1 = 44 Cariri_remainder; 20 − 1 = 19
Sertao_Central_remainder — arithmetic checks out.)

## 7. Completeness / duplicate checks

- 184/184 unique `muni6` values; **0 duplicates**.
- **0** missing `municipality_name`.
- **0** missing `official_macroregion_5`.
- **0** unmatched municipalities against the modelling panel
  (`01_Data/chik_dlnm_panel_muni_week_2015_2025.rds`) — every `muni6` in the
  crosswalk is present in the panel and vice versa.
- **0** name mismatches between the DATASUS/IBGE source names and the
  repository's existing `ibge_muni_name_lookup.rds` names.
- Every municipality belongs to exactly one of the 5 official macroregions
  and exactly one of the 8 `final_model_unit` values (verified by
  `stopifnot` checks in the build script; no municipality appears twice).

---

## Next step / open blocker for regional aggregation (Sections 5–10, not yet attempted)

Building the actual regional weekly panel (`final_model_unit × week_start`,
184-municipality → 8-unit aggregation of cases/`N_start`/`N_end`/births/
deaths, reconciliation, and conservation checks) was **not attempted in this
pass** — the task that produced this report was scoped to the crosswalk and
`final_model_unit` derivation only. Note for when that step is picked up: a
prior repository search found that municipality-level `N_start`/`N_end`/
`births`/`deaths` fields **do not currently exist** — those fields exist only
at **state level** in `01_Data/ceara_weekly_demography.csv`. Producing the
8-unit regional demographic panel will require either (a) a
municipality-level demographic build-out first, or (b) an explicit,
documented method for allocating state-level demographic totals down to the
8 `final_model_unit` regions (e.g. by population share) — this should be
confirmed with the user before proceeding, since it is a modelling choice,
not a data-lookup.
