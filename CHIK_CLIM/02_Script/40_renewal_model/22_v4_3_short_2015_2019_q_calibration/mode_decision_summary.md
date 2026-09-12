# Mode Decision Summary — v4.3/v4.4 q-identifiability audit

Full evidence base: `mode_validity_audit.md`. This file applies the
pre-registered four-outcome decision rule.

## Applying the decision rule

- **Outcome C** (high-q remains genuinely plausible even under 4-week
  aggregation → stop the continuous renewal branch): **ruled out.** The
  chains initialized in the high-q region (q≈0.70–0.90) did not remain
  there in either Fit A4 or Fit B4 — they moved down and merged with the
  low-q-initialized chains in both fits. High-q is not robust to reducing
  the weekly likelihood's pseudo-replication.

- **Outcome A** (4-week model recommended as next production candidate,
  q kept fully uncertain, not fixed at any single value): **partially
  supported but not cleanly met.** Fit B4 (+ serology) does converge to a
  single, epidemiologically coherent, serology-consistent region (q
  median 0.113, 95% CrI [0.043, 0.226]; predicted immune fraction and
  serology posterior-predictive both bracket the observed Juazeiro
  result). But its HMC gate does not pass (max Rhat 1.034, just above the
  1.01 threshold), and Fit A4 (no serology) has **real divergences** (19),
  not just a slow-mixing but healthy geometry.

- **Outcome D** (4-week residual autocorrelation remains severe → do not
  add an AR observation process yet; prepare, do not implement, a
  two-part epidemic-burden + within-epidemic-shape observation model):
  **directly supported.** Block-level lag-1 residual ACF is 0.79–0.82 in
  A4/B4 — essentially unchanged from the weekly analysis (0.81–0.86).
  4-week aggregation reduced the *number* of pseudo-independent
  observation units (which is why the extreme modes merged) but did
  **not** fix the underlying temporal dependence in the observation
  error.

- **Outcome B** (pseudo-replication contributed but did not fully
  resolve identification; do not choose a mode by initialization): **best
  overall fit to the evidence.** The extreme bimodality substantially
  dissolved once pseudo-replication was reduced, which confirms
  pseudo-replication was a major driver of it — but neither fit achieves
  a clean HMC pass, and the observation model's core mis-specification
  (temporal residual dependence) persists at the block level.

**Combined verdict: Outcome B + Outcome D.** Reducing weekly
pseudo-replication (4-week blocks) collapsed the pathological high-q/low-q
separation and produced a single, more defensible region when serology is
included — but this is not yet a validated production result, because (i)
neither fit passes the strict HMC gate, and (ii) the residual
autocorrelation that motivated this diagnostic in the first place is still
present at the block level, just as strongly as at the weekly level.

## Recommended next action (exactly one)

**Do not adopt a weekly or 4-week block NB observation likelihood as the
production model, and do not select or fix a single q value from either.**
Instead, design (but do not yet implement or fit) a two-part
epidemic-burden + within-epidemic-shape observation model, in the spirit
of Outcome D: one component modeling total case burden per (longer) period
at a resolution where the conditional-independence assumption is
defensible (informed by the observed effective-sample-size reduction of
roughly 10×, based on the AR(1)-implied effective N of ~25–30 out of 313
weekly points), and a separate Dirichlet-multinomial-style component
modeling the within-period temporal *shape* of cases, which does not
require treating each week as an independent NB draw. This design should
be reviewed before any fitting begins, and should preserve — not
pre-empt — the full q uncertainty seen in Fit B4 (q ≈ 0.04–0.23), rather
than centering a new prior on 0.05 or on 0.11 merely because those regions
appeared in this audit's diagnostic fits.
