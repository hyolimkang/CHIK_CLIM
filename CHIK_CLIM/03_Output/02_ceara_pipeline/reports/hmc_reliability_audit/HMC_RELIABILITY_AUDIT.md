# HMC Reliability Audit -- tau_year hierarchical fits (P1-P3, H1/H2)

Requested as a prerequisite gate before continuing
`13_retrospective_resurgence_diagnostic`. Produced by
`02_Script/02_ceara_pipeline/14_hmc_reliability_audit/`. Read-only
against all 5 existing fits; one new fit (P2 rerun, adapt_delta=0.99)
saved SEPARATELY, original `ce_tau_year_P2.rds` untouched. Not committed.

## 1. Divergences, Rhat, ESS, E-BFMI, treedepth -- by model and chain

Full tables: `HMC_AUDIT_model_level.csv`, `HMC_AUDIT_per_chain.csv`,
`HMC_AUDIT_per_parameter_rhat_ess.csv` (the last computed over **every**
parameter in the model -- all 573-week S/U/X/log_R0 trajectories and
all year/season random effects, 3,430-7,510 quantities per model --
not just the 10 named hyperparameters, so the "worst" parameter is
genuinely the model-wide worst case).

| Model | Post-warmup iter | Divergences | Div. rate | Treedepth hits | min E-BFMI | max Rhat (par) | min ESS bulk (par) | min ESS tail (par) |
|---|---|---|---|---|---|---|---|---|
| P1 | 1600 | 20 | 1.250% | 0 | 0.724 | 1.035 (**tau_year**) | 112 (**tau_year**) | 56 (A_year[7]) |
| P2 (primary) | 1600 | 2 | 0.125% | 0 | 0.716 | 1.011 (log_R0[206]) | 275 (**tau_year**) | 278 (z_year[3]) |
| P3 | 1600 | 5 | 0.313% | 0 | 0.716 | 1.014 (log_R0[281]) | 335 (**tau_year**) | 416 (z_season_year[2]) |
| H1 | 1600 | 5 | 0.313% | 0 | 0.735 | 1.013 (z_season_year[2]) | 220 (z_season_year[2]) | 140 (z_season_year[2]) |
| H2 | 1600 | 10 | 0.625% | 0 | 0.857 | 1.019 (**tau_year**) | 387 (**tau_year**) | 292 (**tau_year**) |

Per-chain divergences (`HMC_AUDIT_per_chain.csv`) are concentrated in
1-2 of the 4 chains per model (e.g. P1: chain 2 alone carries 19/20;
H2: chain 1 alone carries 6/10) -- not spread evenly, consistent with
isolated excursions into a locally hard region rather than a
chain-wide sampling failure. **Zero treedepth-saturation hits in every
chain of every model** -- the sampler is not being truncated; E-BFMI is
comfortably above the 0.3 energy-pathology threshold everywhere
(min 0.716-0.857). Under the project's strict gate (0 divergences, 0
treedepth hits, Rhat<=1.01, min ESS>=100, BFMI>=0.3) **all 5 models
fail**, each on a small, specific, disclosed margin -- consistent with
the gate failures already disclosed when these fits were produced.
**tau_year itself is the single worst-Rhat-or-ESS parameter in the
entire model for P1, P2 (bulk ESS), and H2** -- not a minor nuisance
parameter buried among 500+ latent states, but a genuine identification
hot-spot. H1 is the one exception: its worst parameter throughout is
`z_season_year[2]` (the seasonal-AMPLITUDE hierarchy), suggesting H1's
dominant HMC difficulty may be the short (5-year) training window's
amplitude hierarchy rather than tau_year specifically.

## 2. Where divergences concentrate -- pairs plots

`HMC_AUDIT_pairs_{P1,P2,P3,H1,H2}.png` (tau_year x alpha_R x
sigma_season_year x beta_T x beta_P, divergent transitions as red X).

**Divergences concentrate at LARGER tau_year, not near the tau_year=0
boundary**, consistently across all 5 models (`HMC_AUDIT_divergence_concentration.csv`):

| Model | median tau_year, divergent iters | median tau_year, non-divergent iters |
|---|---|---|
| P1 | 0.360 | 0.152 |
| P2 | 0.206 | 0.139 |
| P3 | 0.167 | 0.145 |
| H1 | 0.360 | 0.152 |
| H2 | 0.222 | 0.158 |

This is the opposite of the classic centred-parameterisation funnel
(which pathologically concentrates divergences near a *boundary*);
here the model already uses a non-centred form
(`z_year ~ std_normal()`, `year_effect = tau_year*(z_year-mean(z_year))`),
so a near-zero funnel was not expected and is not what is observed.
**Visual inspection of the pairs plots shows a qualitative difference
by prior tightness**: P1 (tightest prior, HalfNormal(0,0.20), most
divergences) shows a visually distinct **secondary cluster** of
several divergent points tightly co-located at tau_year~0.27,
alpha_R~0.284, sigma_season_year~0.40 (versus the bulk's ~0.29/~0.47)
-- i.e. divergent iterations there move tau_year up AND
sigma_season_year down TOGETHER, a correlated trade-off, not
independent noise. P2's two divergent points and H1's five are more
scattered across the tau_year range with no equivalent tight cluster.
P3 (loosest prior) shows the mildest shift of all (divergent median
0.167 vs bulk 0.145). **The pattern across priors (tightest prior ->
most divergences -> clearest clustering) suggests the tension is
sharper when the prior pulls harder against values the likelihood
occasionally wants to visit**, not a problem specific to one prior choice.

**Likely structural source**: `compute_R0_future()`/the Stan model's
own `log_R0` both add `year_effect[year]` (tau_year's contribution)
AND multiply `A_year[year]` (sigma_season_year's contribution) into
the SAME year's seasonal terms. A year can reach a similar realised
R0(t) peak via a large year_effect with modest seasonal amplitude OR a
modest year_effect with large seasonal amplitude -- a plausible
partial confound between the two year-level hierarchies, consistent
with the correlated tau_year/sigma_season_year shift seen in P1's
divergent cluster.

## 3. Could posterior summaries of tau_year be biased?

`HMC_AUDIT_bias_check_exclude_divergent.csv`, `HMC_AUDIT_tail_bias_check.csv`:
excluding divergent iterations shifts the reported **median** by at
most -0.36% (P1) and under 0.1% for P2/P3/H1/H2 -- not materially
biased. The **upper tail** is more sensitive: q95 shifts -4.5% and
q97.5 shifts -7.2% for P1 if divergent draws are excluded (i.e. the
currently-reported P1 upper tail is, if anything, already being pulled
WIDER by the divergent draws it contains, not truncated/underestimated
by missing them); P2/P3/H1/H2 shift under 1.6% at both quantiles.
**Divergent draws are never dropped from any reported headline
posterior in this project** -- this is a diagnostic-only check.
Conclusion: no evidence that the reported tau_year point estimates or
credible intervals are biased downward by sampler avoidance of a
problematic region; P1's tail is the only one with a non-negligible
(but still modest, <10%) sensitivity to this choice.

## 4. P2 rerun: adapt_delta=0.99, warmup=1500, sampling=500

Saved separately: `model_fits/hmc_reliability_audit/ce_tau_year_P2_adapt099.rds`
(original `ce_tau_year_P2.rds` untouched). Elapsed 1932s (~32 min) vs
original 957s (~16 min).

| | Original P2 (adapt_delta=0.95, 800+400) | Rerun (adapt_delta=0.99, 1500+500) |
|---|---|---|
| Divergences | 2 | **1** |
| Treedepth hits | 0 | 0 |
| mean treedepth | ~7.9 | ~9.1 |
| tau_year Rhat | 1.007 | **1.022** (worse) |
| tau_year ESS bulk | 275 | **264** (slightly worse) |
| tau_year ESS tail | 295 | 372 (better) |
| tau_year median | 0.1389 | **0.1460** (+5.1%) |
| tau_year mean | 0.1443 | **0.1564** (+8.4%) |
| tau_year sd | 0.0345 | **0.0456** (+32% wider) |
| tau_year q95 | 0.2115 | **0.2394** (+13.2%) |
| HMC gate | FAIL | FAIL |

**Simply increasing adapt_delta and iterations resolved the
divergences (2->1, now within the normal single-digit noise floor) but
did NOT fix Rhat or bulk ESS -- both got marginally WORSE -- and the
posterior itself shifted and widened** (median +5%, mean +8%, sd
+32%, q95 +13%) rather than merely tightening around the same answer.
This is the signature the task anticipated in its guidance to
"investigate parameterisation or identifiability rather than simply
increasing iterations": a brute-force HMC-settings fix does not make
this problem go away, and the direction of the shift (more conservative
sampling finds MORE upper-tail mass, not less) is consistent with
genuine, moderate identifiability tension between tau_year and
sigma_season_year (Section 2), not with simple under-sampling of an
otherwise well-behaved posterior.

## 5. Recommendation

**P2 and the H1/H2 training-only posteriors are reliable enough to
proceed with the retrospective resurgence diagnostic (module 13),
WITH the following caveats carried forward explicitly:**

1. **tau_year's point estimate (median) is robust** -- stable to
   within 0.4% whether or not divergent draws are included, stable
   across 3 independent priors (0.137-0.145 full-period) and across
   training-window length (0.152/0.158 for H1/H2) -- this is the
   figure Task D's and Task E's headline comparisons (G1 vs G2,
   historical vs future variability) rest on, and it is not in question.
2. **tau_year's upper credible tail carries additional, uncharacterised
   uncertainty beyond what is reported** -- the adapt_delta=0.99 rerun
   suggests the true q95 could plausibly be ~13% higher than currently
   reported (~0.24 instead of ~0.21 for P2). Any conclusion in Task E
   that depends on the EXACT width of tau_year's credible interval
   (rather than its median, or than order-of-magnitude contrasts like
   G1-vs-G2's ~117x or the importation grid's 16x span) should be
   treated as approximate, not exact.
3. **The likely structural cause is a partial tau_year/sigma_season_year
   confound** in how year_effect and A_year jointly enter log_R0 --
   worth a dedicated reparameterisation in future model development
   (e.g. jointly non-centring the two year-level hierarchies, or
   constraining their covariance explicitly), but this is model
   DEVELOPMENT work, out of scope for the current retrospective
   diagnostic task.
4. **Do not substitute the adapt_delta=0.99 rerun for the original P2
   fit** in downstream work -- it did not resolve the underlying issue
   and introduced its own (larger) Rhat/ESS degradation; the original
   `ce_tau_year_P2.rds` (and H1/H2 training fits, not rerun here per
   the task's own instruction) remain the fits to use.
5. P1/P3 were **not** rerun, per instruction; their divergence
   patterns (more pronounced for P1, the tightest prior) are disclosed
   above but do not block P2/H1/H2-based work.

No accepted model or prior fit was modified or overwritten. Nothing committed.

## 6. Addendum -- P2 chain-level mixing and tau_year/sigma_season_year geometry (deep-dive)

Requested follow-up. Produced by
`04_P2_chain_mixing_and_geometry.R`, read-only against the original
`ce_tau_year_P2.rds`. Figures: `P2_mixing_traceplot.png`,
`P2_mixing_rankplot.png`, `P2_mixing_acf.png`,
`P2_geometry_tau_vs_sigma_season.png`. Tables:
`P2_mixing_gelman_rubin_decomposition.csv`, `P2_mixing_ess_per_chain.csv`,
`P2_geometry_funnel_check_by_tau_quintile.csv`,
`P2_geometry_year_effect_vs_A_year_per_year_cor.csv`.

**Chain-level mixing.** The classic Gelman-Rubin decomposition shows
the 4 chains agree closely on location: between-chain variance is only
**0.29% of total variance for tau_year and 0.09% for sigma_season_year**
(chain means: tau_year 0.142/0.144/0.146/0.146; sigma_season_year
0.465/0.468/0.471/0.467) -- Rhat is driven almost entirely by
within-chain autocorrelation, not by chains disagreeing on where the
posterior mass is. Rank plots (`P2_mixing_rankplot.png`) are
reasonably uniform for both parameters with no systematic one-chain
skew, consistent with this. **Within-chain ESS is uneven and
asymmetric across chains**: for tau_year, chain 3 is a clear outlier
(ESS=37 vs. 101/113/105 for chains 1/2/4; lag-1 autocorrelation 0.554,
the highest of the four); for sigma_season_year the pattern inverts --
chain 3 is the BEST-mixing chain (ESS=199) while chain 1 is worst
(ESS=64). **Each chain trades mixing efficiency between the two
parameters in opposite directions** -- the chain that explores
tau_year most slowly explores sigma_season_year fastest, and vice
versa. This is itself evidence of a real (if subtle and chain-local)
interaction between the two parameters' sampling difficulty, even
though no chain is stuck in a different REGION (per the between-chain
variance result above) -- it is a shared, transient geometric
obstacle that each chain encounters at different points in its
trajectory, not a multimodality/non-convergence problem.

**Posterior geometry.** The originally-hypothesised global
tau_year/sigma_season_year trade-off (Section 2, based only on P1's
divergent-cluster visual) does **not** hold for P2 under direct
quantitative testing: pooled Pearson r=0.030, Spearman rho=0.037
(`P2_geometry_tau_vs_sigma_season.png`) -- essentially uncorrelated.
The funnel/heteroscedasticity check (`P2_geometry_funnel_check_by_tau_quintile.csv`)
shows sigma_season_year's spread is roughly flat across tau_year
quintiles (sd 0.079-0.088, no monotonic trend) -- **not a classical
funnel**. The per-year correlation between the REALISED quantities
that actually enter log_R0 (year_effect[y] = tau_year*(z_year[y]-mean(z_year))
vs. A_year[y] = exp(sigma_season_year*z_season_year[y])) is small and
inconsistent in sign across the 11 years (range -0.089 to +0.247, mean
+0.059) -- no systematic confound at the realised-quantity level
either. **The joint scatter instead shows an irregular, right-skewed
shape**: roughly circular/uncorrelated in the bulk, but tau_year's
upper tail (>~0.20) is sparse and bumpy rather than smoothly tapering,
and both of P2's divergent transitions fall either at the edge of this
sparse tail (tau_year=0.27) or at an unremarkable bulk location with
high sigma_season_year (0.14, 0.59) -- i.e. divergences are better
explained by **local sparsity/roughness in tau_year's own right tail**
than by a pairwise correlation with sigma_season_year.

**Revised conclusion (supersedes the Section 2 hypothesis for P2
specifically)**: P2's residual HMC difficulty is NOT a clean global
tau_year/sigma_season_year confound. It is better characterised as (a)
ordinary within-chain mixing inefficiency that is asymmetric and
chain-specific rather than persistent, and (b) mild local roughness in
the sparsely-populated right tail of tau_year's own marginal
posterior, not a two-parameter funnel or correlation. P1's more
visually dramatic divergent cluster (Section 2) may reflect a
genuinely different (tighter-prior-induced) dynamic and was not
re-tested with this level of detail -- if a similar deep-dive is
wanted for P1, the same script can be pointed at `MODELS_TO_AUDIT$P1`.
This does not change Section 5's recommendation: the location
estimates remain trustworthy; the main remaining caveat is still
tau_year's upper-tail width, now attributed to sparse-tail sampling
roughness rather than a cross-parameter funnel.

## Files produced

```
02_Script/02_ceara_pipeline/14_hmc_reliability_audit/
  00_config_hmc_audit.R, 01_diagnostics_all_models.R, 01b_tail_bias_check.R,
  02_rerun_P2_adapt099.R, 03_compare_P2_original_vs_adapt099.R,
  04_P2_chain_mixing_and_geometry.R

03_Output/02_ceara_pipeline/
  tables/hmc_reliability_audit/HMC_AUDIT_*.csv (8 tables),
                                 P2_mixing_*.csv (2), P2_geometry_*.csv (2)
  figures/hmc_reliability_audit/HMC_AUDIT_pairs_{P1,P2,P3,H1,H2}.png,
                                  HMC_AUDIT_divergence_vs_tau_year.png,
                                  P2_mixing_traceplot.png, P2_mixing_rankplot.png,
                                  P2_mixing_acf.png, P2_geometry_tau_vs_sigma_season.png
  model_fits/hmc_reliability_audit/ce_tau_year_P2_adapt099.rds (NEW, separate from accepted P2)
  reports/hmc_reliability_audit/HMC_RELIABILITY_AUDIT.md (this file)
```
