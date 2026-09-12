# Dynamic annual FOI diagnostic rules

Primary HMC acceptance requires zero divergences, zero maximum-treedepth hits,
maximum R-hat <= 1.01, bulk and tail ESS >= 400, and BFMI >= 0.3 for each
chain. Otherwise outputs are labelled **UNCONVERGED — DO NOT INTERPRET
SUSCEPTIBILITY**.

The diagnosis writes the HMC headline, BFMI, worst R-hat/ESS parameters,
annual posterior summary, posterior correlation matrix, scalar summary,
diagnostic PDF, traceplot PDF, and prior-predictive PDF. The key
identifiability checks are correlations of q with mu_lambda, cumulative
infections, S at end-2025, and each annual lambda. Strong q/infection-scale
dependence is a reason to stop rather than add more latent structure.

M0 and M1 are compared on HMC, posterior predictive RMSE, q, infection burden,
and end susceptibility. In-sample fit alone is not a model-selection rule.
