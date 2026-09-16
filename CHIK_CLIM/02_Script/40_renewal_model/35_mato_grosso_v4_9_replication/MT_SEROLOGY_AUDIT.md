# Mato Grosso (MT) Serology Audit — Information Only

Per Phase-1 Section 6: this audit is for the record only. **No serology is
used in the Phase-1 global-q case-only fit.**

## Source

`Ref/Brazil_CHIKV_serology_inventory_2026-09-10.xlsx` (search cutoff
2026-09-10, v1.0) — the same master inventory used to source the RJ (U10),
PE (U14), and Bahia serology anchors elsewhere in this project.

## UF_coverage sheet, MT row

| field | value |
|---|---|
| uf | MT |
| state_name | Mato Grosso |
| macroregion | Central-West |
| source_papers | 0 |
| sample_site_time_records | 0 |
| local_population_priority_records | 0 |
| has_state_representative_survey | NOT ESTABLISHED |
| source_ids | (none) |
| coverage_note | "No extractable observation identified in this search; NOT evidence that no survey exists." |

## Raw_observations sheet

Cross-checked directly: the `uf` column of `Raw_observations` contains
observations only for **AC, AP, BA, CE, DF, MG, PE, RJ, SP** (9 UFs total,
matching the "UFs with at least one extracted observation" count on the
inventory's summary page). **Mato Grosso does not appear.**

## Classification

**No independent MT serological anchor identified.**

This is an absence-of-evidence finding from a single systematic search
pass, not a claim that no MT chikungunya serosurvey has ever been
conducted — the inventory's own coverage note is explicit on this point.
If a defensible MT-specific serosurvey is identified later (state- or
municipality-representative, with reported numerator/denominator or a
reported prevalence + CI), it can be considered for a serology-augmented
follow-up model, analogous to the RJ CITY U10 analysis. Until then, MT
absolute-scale inference can only be assessed via case-only
identification diagnostics (q-alpha_R ridge, prior sensitivity, targeted
fixed-q profiling if warranted) — never via HEMORIO/U14/U10-style
external validation, since none of those apply to MT geographically.
