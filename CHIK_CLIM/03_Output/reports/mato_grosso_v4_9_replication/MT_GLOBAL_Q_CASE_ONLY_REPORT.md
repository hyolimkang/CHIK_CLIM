# Mato Grosso (MT) — Phase 1: Frozen v4.9 Global-q Case-Only Report

## Model

Frozen v4.9 state-level homogeneous renewal model, reused byte-identical
from CE/BA/PE/RJ (`renewal_bahia_v4_9_global_q_multisite_serology.stan`,
state-agnostic, `J_sero=0`). One estimated global `q_MT`. No serology
(none identified for MT anyway — see `MT_SEROLOGY_AUDIT.md`). No seed (no
MT inter-wave gap qualifies under the 399-week Ceara-2022 standard; max
gap = 184 weeks). No spatial structure, regional q/R0, mobility, or
Pernambuco-style reparameterisation.

## HMC run history

| run | adapt_delta | divergences | max treedepth hits | max Rhat | min bulk ESS | verdict |
|---|---|---|---|---|---|---|
| 1st pass | 0.95 | 2 (both chain 2) | 0 | 1.0089 | 932 | NEAR-PASS |
| rescue (Section 10 rule) | 0.98 | 0 | 0 | 1.0056 (all params) / 1.0023 (gate pars) | 1035 | **PASS** |

The 1st-pass run showed exactly the qualifying "isolated numerical noise"
signature (Section 10): 1-2 divergences, zero treedepth-14 hits, chains
clearly overlapping, healthy BFMI (0.84-0.93) and step sizes. Per the
predefined rule this authorized exactly ONE rescue at adapt_delta=0.98,
same warmup/sampling (1200/750), same model/priors/inits. The rescue
passed cleanly and is retained as the MT computational fit (canonical
path: `outputs/caseonly/mt_global_q_case_only.rds`; the ad095 near-pass
run is archived as `..._ad095_nearpass_ARCHIVED.rds`). No further
adapt_delta escalation or reparameterisation was needed or attempted.

## q / S identification diagnostics (final PASS fit)

- q prior: broad, `logit_q ~ Normal(logit(0.10), 1)`, 95% width = 0.424
- **q_MT posterior: median = 0.0407, 95% CrI = [0.0336, 0.0524]** (95%
  width = 0.019, ratio to prior = 0.044 — posterior meaningfully narrower
  than the prior)
- **cor(logit_q, alpha_R) = -0.793** — a STRONG ridge, of the same
  character and magnitude as the RJ q-alpha_R ridge (|cor| 0.70-0.84
  across RJ's case-only/city models). A narrow marginal q posterior
  alongside this |cor| does **not** constitute clean identification, per
  the standing convention established for RJ.
- cor(q, immune_2025) = -0.873, cor(q, min S/N) = 0.873, cor(q, max R0) = -0.097

## Checkpoints (S/N, cumulative immune fraction, R0)

| checkpoint | S/N median [95% CrI] | immune fraction median [95% CrI] | R0 median |
|---|---|---|---|
| 2016-12-31 | 0.999 [0.999,1.000] | 0.0006 [0.0004,0.0008] | 2.62 |
| 2018-12-31 | 0.904 [0.878,0.927] | 0.096 [0.073,0.122] | 1.48 |
| 2020-12-31 | 0.903 [0.877,0.927] | 0.097 [0.073,0.123] | 1.54 |
| 2022-12-31 | 0.904 [0.879,0.928] | 0.096 [0.072,0.121] | 1.78 |
| 2025-12-21 | 0.472 [0.394,0.581] | 0.528 [0.419,0.606] | 1.97 |

The 2018-2022 period is essentially flat (S/N pinned at ~90%, consistent
with the genuinely quiet intermediate epidemic period). Almost all
model-implied cumulative infection (an S/N drop from ~90% to ~47%, i.e.
~43 percentage points) is concentrated in the 2024-2025 recurrence alone.

## Case posterior predictive check

- Annual and weekly PPC figures: `MT_global_q_case_PPC.png`,
  `MT_global_q_annual_PPC.png`.
- Per-wave check (early 2017-2019 epidemic vs. late 2024-2025 recurrence):

| wave | era | observed total | predicted median [95% CrI] | within 95% PPC? |
|---|---|---|---|---|
| MT_wave_05 | early | 2,929 | 3,648 [2,756, 4,911] | yes |
| MT_wave_06 | early | 13,419 | 9,955 [7,781, **12,708**] | **no — observed exceeds upper bound** |
| MT_wave_07 | early | 247 | 231 [173, 300] | yes |
| MT_wave_13 | late | 19,449 | 20,849 [15,710, 28,569] | yes |
| MT_wave_14 | late | 47,315 | 46,249 [36,093, 58,394] | yes |

**One wave (MT_wave_06, the largest early-period wave) is mildly
underpredicted** — observed 13,419 vs. a 95% predictive upper bound of
12,708 (~5.6% overshoot beyond the 97.5th percentile). This is a modest,
single-wave miss, not a systematic failure of one era: the late
recurrence (both 2024 and 2025 waves) is well reproduced, and two of the
three early-period waves are also well reproduced. **One shared
transmission/susceptibility history can reproduce both epidemic eras
essentially, with one modest early-wave underprediction** — the model
does not achieve its late-recurrence fit by sacrificing the early period,
or vice versa.

## Biological plausibility

- Immune fraction immediately after the early epidemic (end MT_wave_06,
  2019-08-18): median = 9.75% [7.34%, 12.34%]
- S/N immediately before the late recurrence (onset MT_wave_13,
  2024-01-07): median = 90.30% [87.81%, 92.67%] — still overwhelmingly
  susceptible going into 2024, consistent with the ~4x larger 2024-2025
  recurrence.
- Max R0(t) overall: median = 3.13 [2.84, 3.46] — occurs during the
  2024/2025 peaks, i.e. exactly the period of steepest depletion. This is
  the highest seasonal R0 peak in the whole 2015-2025 series but is not
  biologically implausible for chikungunya (R0 estimates elsewhere in
  this project's other states reach comparable ranges).
- Min S/N overall: 47.18% [39.43%, 58.13%]; cumulative infection
  (immune_2025): 52.82% [41.87%, 60.57%].
- No evidence of "extreme depletion + implausibly high R0" pathology: the
  recurrence is large, but the model does not need to invoke an
  extreme/unphysical R0 to reproduce it.

## Phase-1 Final Classification

**B. MT COMPUTATIONALLY CLEAN BUT q/S IDENTIFICATION UNCERTAIN**

Rationale: the HMC gate is cleanly PASSED (0 divergences, max Rhat=1.006,
healthy ESS/BFMI) and the case PPC is essentially good (4/5 waves fully
within the 95% predictive interval, one wave modestly underpredicted).
However, the q-alpha_R ridge (|cor|=0.79) is as strong as the one that
has driven weak absolute-scale identification in every other state
modeled this way (RJ in particular). A narrow q posterior sitting on top
of a ridge this strong is not, by the project's own standing convention,
evidence of genuine identification — it reflects how tightly the prior
and case-curve shape jointly constrain the transmission-level/observation
scale, not an independently well-identified ascertainment fraction. MT is
NOT yet an accepted absolute-susceptibility reconstruction on the basis
of this Phase-1 run alone.

## MT GLOBAL-q CASE-ONLY PHASE COMPLETE — STOP FOR REVIEW
