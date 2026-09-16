# Rio de Janeiro Serology Audit

**Status: audit only. NO serology observation is included in the Phase-1
global-q case-only model.**

Source: `Ref/Brazil_CHIKV_serology_inventory_2026-09-10.xlsx` (the same
master inventory used to source Bahia's 6 multisite surveys and
Pernambuco's U14/P10 survey), sheets `UF_coverage`, `Survey_units`,
`Raw_observations`, `Sources`.

## UF_coverage summary

State-level representativeness for RJ: **`NOT ESTABLISHED`**. Note:
*"Contains local / selected samples only; no whole-state estimate
established."* Two survey units on file (`U10`/`P07`, `U35`/`P21`).

## Survey inventory

### U10 / P07 — Perisse et al. 2020 ("Zika, dengue and chikungunya
population prevalence in Rio de Janeiro city, Brazil...", PLOS ONE, doi
10.1371/journal.pone.0243239)

| Field | Value |
|---|---|
| Location | **Rio de Janeiro city** (municipality), NOT Rio de Janeiro state |
| Survey dates | July-October 2018 |
| Sampling design | Administrative-region stratified household survey, all ages, design/nonresponse weights |
| n tested | 2,120 |
| n positive | 371 |
| Raw prevalence | 0.175 |
| Reported survey-weighted prevalence | 0.180 (95% CI 0.148-0.212) |
| Assay | Rapid IgM/IgG test (BahiaFarma) |
| Extraction tier | `PRIMARY_FULL_TEXT` |
| Flag | `LOCAL_POPULATION_PRIORITY`; explicit note: "City population-based; NOT Rio de Janeiro state" |
| Note on raw vs. weighted | "Raw 371/2120 differs from weighted 18.0%; retain both." |

**Classification: B — municipality/local-population survey.** This is
the direct RJ analogue of Pernambuco's U14 (Recife). If used at all, it
would require the SAME local-to-state transport framework (geographic
offset) already used for Bahia/PE, or a state-specific-stratum treatment
if RJ is later spatially disaggregated -- explicitly NOT to be treated as
a state-representative anchor. Not used in Phase 1.

### U35 / P21 — Sant'Anna et al. (HEMORIO blood-centre donors, Rio de
Janeiro, 2019-2022 pooled)

| Field | Value |
|---|---|
| Location | HEMORIO blood centre, Rio de Janeiro |
| Survey dates | 2019-2022, pooled (no within-period breakdown) |
| Sampling design | Blood donors (donation-eligible population only, NOT general population) |
| n tested | 778 |
| n positive | **NOT retrieved** -- only rounded percentages available (21.3% unweighted "CHIKV prevalence, exact marker not recovered"; 4.7% unweighted IgM-only) |
| Assay | Serological/molecular screening; full definition not retrieved for the primary 21.3% endpoint; IgM-only definition given for the 4.7% figure |
| Extraction tier | `PRIMARY_ABSTRACT` only (weaker than U10's full-text extraction) |
| Flag | `HOLD_ENDPOINT_DEFINITION`; explicit instruction on file: "Numerator NOT retrieved; do not reconstruct from rounded %" |

**Classification: D — duplicate/unusable (for quantitative fitting).**
Blood donors are a doubly-selected population (healthy, donation-eligible,
self-selected) and, independent of that, the source data explicitly
prohibits reconstructing a raw numerator from the rounded percentage. Two
candidate endpoints (21.3% all-CHIKV, 4.7% IgM-only) are inconsistent with
each other and neither has a usable `n_positive`. Not usable as a
quantitative likelihood observation under any framework without further
primary-source retrieval.

## Conclusion

No RJ survey is state-representative. One survey (U10/Perisse 2020,
Rio de Janeiro city) is a defensible LOCAL anchor with usable raw counts
(371/2120) and a full-text extraction, directly analogous to Pernambuco's
U14 -- it is a serious candidate for a LATER serology-augmented model, but
is explicitly excluded from Phase 1 per instruction. The second survey
(HEMORIO donors) is not usable quantitatively regardless of phase, due to
missing numerator data and a doubly-selected sampling frame.

Phase 1 (this report) proceeds as **case-only**, consistent with: *"No
independent serological anchor is available for RJ [as a
state-representative source]"* -- continuing with case-only
identifiability analysis, per instruction.

## Outputs

- `RJ_SEROLOGY_AUDIT.md` (this file)
- `03_Output/tables/rio_de_janeiro_v4_9_replication/RJ_serology_inventory.csv`
