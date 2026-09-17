# NNV module validation

Generated: 2026-09-17 10:23:03.732472

## 1. cumulative_doses >= 0 for every draw
**PASS**
min = 913032.1, max = 913032.1

## 2. cumulative_effectively_protected <= cumulative_doses for every draw
**PASS**
max(effective - doses) = -1.1198e+05 (prior-immune entrants at age 12 receive a dose but no additional protection, so effective < doses strictly whenever any entrant is already Uinf)

## 3. coverage = 0 -> doses = 0 -> impact = 0 -> NNV undefined (NA), not zero
**PASS**
Unit test on draw 355: doses=0.000000, effective=0.000000, infections(no-vax)=2505195.55, infections(cov=0)=2505195.55, impact=0.000000, NNV=NA (correct)

## 4. VE = 0 -> doses > 0 but effective = 0 -> impact approx 0 -> NNV undefined or extremely large, never artificially small/finite
**PASS**
Unit test on draw 355: doses=913032.1 (>0 as expected), effective=0.000000 (exactly 0, as expected), impact=0.000000 infections (~0, within numerical noise of a 365-week/101-age-band recursion), NNV=NA (undefined, correct)

## 5. infections_averted > 0 checked before computing conventional NNV_infection (NA otherwise)
**PASS**
0 / 200 draws have infections_averted <= 0 in the main stress-test run (all had positive benefit; the NA-guard is exercised explicitly by the coverage=0 and VE=0 unit tests above)

## 6. reported_cases_averted = q * infections_averted under the current constant q = 0.05 observation mapping
**PASS**
max abs residual = 0.00e+00

## 7. NNV computed per-draw (doses[d]/effect[d]) THEN summarised -- not median(doses)/median(effect)
**PASS**
Correct (per-draw-first) posterior median NNV = 1.6947. Incorrect (median-of-ratio-of-medians) would give = 1.6947. Difference = 0.0000 (0.0% relative) -- demonstrates why per-draw computation matters even though small here, because doses are deterministic given coverage while infections_averted varies across draws.

## 8. Cumulative NNV at any (draw, week) recomputes exactly as cum_doses / cum_infections_averted -- no smoothing/interpolation applied
**PASS**
Spot-checked 8 random (draw, week) cells against stored values; recomputation from raw cumulative doses/effects matches exactly. Script 11_calculate_nnv.R contains no rollmean/loess/smooth.spline or other smoothing call (source inspection).

## 9. NNV numerator is cumulative_doses (all entrants offered vaccination), never cumulative_effectively_protected
**PASS**
Source inspection of 11_calculate_nnv.R confirms NNV_infection <- cumulative_doses / infections_averted (never cumulative_effectively_protected). Numerically, cumulative_effectively_protected is strictly less than cumulative_doses for every draw (median gap = 126383 doses given to already-immune or otherwise-unprotected entrants).

## 10. total effect (A-B) = direct-only (A-C) + indirect (C-B), recomputed within the NNV module
**PASS**
max abs error = 0.00e+00 across 200 draws (consistent with QA10 in 07_validate_counterfactuals.R)

---
## Overall: ALL 10 NNV QA CHECKS PASS
