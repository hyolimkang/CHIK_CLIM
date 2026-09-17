# Pernambuco 5-Strata Input Audit

Descriptive pre-model audit only. No Stan fitting.

## Stratum population (final observed week, 2025-12-21)

| Stratum | Population |
|---|---|
| 1_Recife | 1,588,376 |
| 2_Metropolitana_remainder | 4,038,134 |
| 3_Agreste | 1,978,503 |
| 4_Sertao |   881,487 |
| 5_Vale_Sao_Francisco_Araripe | 1,072,166 |

Confirmed: summing the 5 strata's weekly cases reproduces the exact PE state-level series (max abs diff = 0 across 573 weeks).

## Dominant stratum per wave

| Wave | Dominant stratum | Cases in wave | Share of wave |
|---|---|---|---|
| PE_wave_03 | 1_Recife | 5862 | 37.0% |
| PE_wave_04 | 1_Recife | 312 | 47.6% |
| PE_wave_05 | 1_Recife | 772 | 56.8% |
| PE_wave_06 | 1_Recife | 2066 | 68.8% |
| PE_wave_07 | 1_Recife | 12448 | 62.7% |
| PE_wave_08 | 5_Vale_Sao_Francisco_Araripe | 8189 | 63.0% |
| PE_wave_09 | 1_Recife | 361 | 46.5% |
| PE_wave_10 | 1_Recife | 996 | 45.1% |
| PE_wave_11 | 1_Recife | 570 | 47.3% |

## Key-wave stratum breakdown

### PE_wave_03 (2016 index epidemic)

| Stratum | Cases | Share |
|---|---|---|
| 1_Recife | 5862 | 37.0% |
| 2_Metropolitana_remainder | 5225 | 33.0% |
| 3_Agreste | 3652 | 23.0% |
| 5_Vale_Sao_Francisco_Araripe | 831 | 5.2% |
| 4_Sertao | 279 | 1.8% |

### PE_wave_07 (2021 epidemic)

| Stratum | Cases | Share |
|---|---|---|
| 1_Recife | 12448 | 62.7% |
| 2_Metropolitana_remainder | 6765 | 34.1% |
| 3_Agreste | 472 | 2.4% |
| 5_Vale_Sao_Francisco_Araripe | 92 | 0.5% |
| 4_Sertao | 69 | 0.3% |

### PE_wave_08 (2022 epidemic)

| Stratum | Cases | Share |
|---|---|---|
| 5_Vale_Sao_Francisco_Araripe | 8189 | 63.0% |
| 2_Metropolitana_remainder | 1506 | 11.6% |
| 3_Agreste | 1451 | 11.2% |
| 1_Recife | 1037 | 8.0% |
| 4_Sertao | 811 | 6.2% |

### PE_wave_09 (2022/23 recurrence)

| Stratum | Cases | Share |
|---|---|---|
| 1_Recife | 361 | 46.5% |
| 2_Metropolitana_remainder | 249 | 32.0% |
| 3_Agreste | 96 | 12.4% |
| 5_Vale_Sao_Francisco_Araripe | 39 | 5.0% |
| 4_Sertao | 32 | 4.1% |

### PE_wave_10 (2023/24 recurrence)

| Stratum | Cases | Share |
|---|---|---|
| 1_Recife | 996 | 45.1% |
| 2_Metropolitana_remainder | 898 | 40.7% |
| 5_Vale_Sao_Francisco_Araripe | 156 | 7.1% |
| 3_Agreste | 90 | 4.1% |
| 4_Sertao | 69 | 3.1% |

### PE_wave_11 (2024/25 recurrence)

| Stratum | Cases | Share |
|---|---|---|
| 1_Recife | 570 | 47.3% |
| 2_Metropolitana_remainder | 513 | 42.6% |
| 5_Vale_Sao_Francisco_Araripe | 55 | 4.6% |
| 3_Agreste | 45 | 3.7% |
| 4_Sertao | 21 | 1.7% |

## Interpretation

See PE_5strata_weekly_cases.png: the two-panel figure (common absolute-case
scale on top, common incidence scale below) shows whether PE's recurrent
state-level trajectory is visibly composed of asynchronous regional
epidemics, or whether all 5 strata rise and fall together. See the
dominant-stratum and key-wave tables above for the quantitative answer.
