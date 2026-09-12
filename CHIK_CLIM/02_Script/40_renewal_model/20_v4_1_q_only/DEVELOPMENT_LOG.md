# Ceará renewal-model development log -- v4.1 q-only ablation

## Purpose

A clean, single-factor ablation of the HMC-passing v4.0/td14 backbone
(archived + checksum-verified in `archive/`, matched byte-for-byte by
`scripts/00_verify_ablation_identity.R`). The ONLY scientific change: q_fixed
(data, =0.10) becomes q (parameter, `q ~ beta(16.2, 108.7)`). Annual-effect
parameterisation, year_effect_prior_sd (=0.40, fixed), all other priors,
generation interval, imports_per_week, S/U bookkeeping, negative-binomial
likelihood, and chain initialisation strategy are UNCHANGED from v4.0/td14.

This isolates the earlier confounded experiments: the aborted joint-q
attempt (`../18_v4_1_minimal_no_vaccine/`) changed q AND sigma_year AND the
annual-effect basis simultaneously; the modular-q attempt
(`../19_v4_1_modular_q/`) kept q fixed (not jointly sampled) but varied it
across separate fits. Neither isolates q's own effect on HMC geometry in a
single joint fit. This experiment does.

## q prior check (reproduced for this experiment's own record)

Monte Carlo (n=2e6) of the induced Beta(30,28) x Beta(20,60) product vs
Beta(16.2, 108.7): KS distance 0.0096 (adequate, <=0.02 threshold). Full
table: `q_prior_check.csv`. Prior used as-is: `q ~ beta(16.2, 108.7)` (mean
~0.129, SD ~0.030, 95% interval ~[0.077, 0.194]).

## Stage A -- implementation validation (2026-09-11)

5 warmup + 5 sampling, all 4 chains. Purpose: confirm compile, data
dimensions (N=573 weeks confirmed), q initialised at dispersed prior
quantiles (0.0927/0.1165/0.1395/0.1692), X/S/U recursion valid (no
`reject()` triggered), sample_file CSV writing confirmed growing correctly
during the run. NOT assessed for convergence (divergences/Rhat from this
run are meaningless and were not used for any decision). PASS as an
implementation check.

## Stage B -- geometry canary (2026-09-11): FAIL, unambiguous and severe

4 chains, 300 warmup + 150 sampling, adapt_delta=0.95, max_treedepth=14,
dispersed q inits, otherwise identical settings to v4.0/td14.

**STRICT HMC GATE: FAIL.**

| metric | v4.0/td14 (1000+1000) | v4.1 q-only Stage B (300+150) | ratio |
|---|---|---|---|
| runtime (s) | 1478.0 | 1774.7 | 1.20x (despite 4.4x FEWER iterations -- true per-iteration cost ~5.3x higher) |
| median treedepth | 10 | 13 | 1.30x |
| p95 treedepth | 11 | 14 (= cap) | 1.27x |
| **treedepth==14 hits (post-warmup)** | **0** | **106** (of 600 post-warmup draws, 17.7%) | infinite (qualitative failure that did not exist at all in td14) |
| median n_leapfrog | 1023 | 8191 | 8.01x |
| p95 n_leapfrog | 3071 | 16383 (= max possible) | 5.33x |
| divergences | 0 | 0 | -- (this criterion alone still passes) |
| BFMI range | [0.926, 1.083] | [0.700, 1.118] | degraded but still >=0.30 |
| max Rhat | 1.0090 | 1.0481 | degraded (short run inflates this somewhat) |
| min bulk ESS | 480.8 | 77.6 | degraded (short run reduces this somewhat) |

Per-chain treedepth==14 hits (out of 150 post-warmup draws each): chain
1=20, chain 2=8, chain 3=47, chain 4=31 -- **every chain** hits the cap
repeatedly, unlike a single-chain outlier pattern.

### Mechanism identified: this is NOT a chain-split (unlike the earlier modular-q q50 failure) -- it is a uniform prior-likelihood tension

Per-chain posterior means for alpha_R, beta_cos, phi_obs, and q are tightly
consistent ACROSS ALL 4 CHAINS (e.g. alpha_R: 0.687-0.690 in every chain; q:
0.0276-0.0281 in every chain) -- so, unlike q50's dramatic 3-vs-1 split into
two different basins, all four chains here agree with each other. But they
agree on a value that is in severe tension with the informative prior: **q
posterior mean ~0.028 is roughly 4.6x smaller than the prior mean (0.129),
and well below even the prior's 2.5th percentile (0.077)**. alpha_R has
correspondingly risen far above its td14 value (~0.69 vs ~0.146) --
consistent with the anticipated q <-> alpha_R compensation (a smaller q
requires proportionally more latent infections X, hence higher R0/alpha_R,
to explain the same observed case counts). This creates a persistent,
difficult ridge between q's informative prior and the case likelihood's
preferred region, which NUTS must repeatedly fight through at high
treedepth on nearly every iteration -- explaining both the elevated
treedepth/leapfrog cost and its uniform presence across all four
independently-initialised chains (they all get pulled to the same
difficult corner, rather than getting stuck in different local optima).

**Important caveat**: Stage B is short (150 post-warmup draws); this is not
yet a fully-equilibrated posterior. But the consistency across all 4
independently-initialised chains, converging to the same far-from-prior
region and hitting the treedepth cap there repeatedly, is a strong signal
that this is a real, stable feature of the joint (q, alpha_R, ...)
posterior surface -- not transient warmup noise.

### Decision, per Sections 7/8/16 of the task spec

Per Section 7 ("do not proceed automatically" if Stage B is unhealthy) and
Section 16 ("if the q-only model fails despite preserving the exact
successful v4.0 structure, do NOT immediately try another sampler tuning
trick... conclude there is strong evidence for a structural q <-> latent
infections <-> susceptibility <-> R0 identification problem. STOP and
report."): **stopping here.** Not proceeding to Stage C or D. Not
increasing adapt_delta or max_treedepth. Not weakening/narrowing the q
prior. Output preserved at `outputs/stageB/`.

This result, combined with the earlier modular-q q50 finding (median q node
also failed, via a different chain-split mechanism), now gives two
independent pieces of evidence -- from two structurally different
experiments -- that estimating q (whether jointly in one fit, as here, or
even just fixing it at a value 27% away from 0.10 across separate fits, as
in modular-q) disrupts the HMC geometry that the fixed q_fixed=0.10 backbone
enjoys. Per Section 16, the next scientific option under consideration is
introducing genuinely independent cumulative-infection information (most
likely a carefully selected serology likelihood) -- NOT implemented here;
see `serology_audit_candidates.csv` for the preparation-only audit (Section
15), which additionally found that the "Quixada" (409 tested/289 positive)
value used throughout every existing renewal-model diagnose script since
v3.0 has **no traceable source anywhere in this repository** (flagged, not
resolved) -- a separate, pre-existing data-provenance gap surfaced during
this audit, unrelated to the q-only ablation result itself but relevant to
any future serology-informed model.

## Generation-interval audit

Unchanged from v4.0 (confirmed byte-identical). See
`../19_v4_1_modular_q/generation_interval_audit_summary.csv` /
`generation_interval_audit_weights.csv` for the full audit (Gamma(shape=4,
rate=2) in weeks, mean 2 weeks/14 days, SD 1 week/7 days, truncated at G=8,
99.991% probability retained before renormalisation) -- reused here
unchanged since the generation interval was not touched in this experiment
either.
