// M1v2-NB: same normalized independent annual FOI shape, absolute NB observation model.
data {
  int<lower=2> Y; array[Y] int<lower=0> C;
  vector<lower=0>[Y] N_start; vector<lower=0>[Y] N_end; vector<lower=0>[Y] births; vector<lower=0>[Y] deaths; vector[Y] reconciliation;
  real log_foi_prior_mean; real<lower=0> log_foi_prior_sd; real<lower=0> sigma_year_prior_scale;
  real q_logit_prior_mean; real<lower=0> q_logit_prior_sd;
}
parameters { real log_lambda_bar; real log_sigma_year; vector[Y] z_year; real logit_q; real log_phi; }
transformed parameters {
  real<lower=0> lambda_bar = exp(log_lambda_bar); real<lower=0> sigma_year = exp(log_sigma_year); real<lower=0,upper=1> q = inv_logit(logit_q); real<lower=0> phi = exp(log_phi);
  vector[Y] eta = sigma_year * z_year; real log_mean_weight = log_sum_exp(eta) - log(Y);
  vector<lower=0>[Y] weight = exp(eta-log_mean_weight); vector<lower=0>[Y] lambda=lambda_bar*weight; vector<lower=0,upper=1>[Y] attack_prob=-expm1(-lambda);
  vector<lower=0>[Y] X; vector<lower=0>[Y] S_start; vector<lower=0>[Y] S_end; vector<lower=0>[Y] U_start; vector<lower=0>[Y] U_end; vector<lower=0,upper=1>[Y] S_end_prop; vector<lower=0>[Y] expected_cases;
  S_start[1]=N_start[1]; U_start[1]=0;
  for(y in 1:Y) { real s_after; real u_after; real frac_s; X[y]=S_start[y]*attack_prob[y]; s_after=S_start[y]-X[y]; u_after=U_start[y]+X[y]; frac_s=s_after/N_start[y]; S_end[y]=s_after+births[y]-deaths[y]*frac_s+reconciliation[y]*frac_s; U_end[y]=u_after-deaths[y]*(1-frac_s)+reconciliation[y]*(1-frac_s); S_end_prop[y]=S_end[y]/N_end[y]; expected_cases[y]=q*X[y]+1e-9; if(fabs(S_end[y]+U_end[y]-N_end[y])>1e-5*N_end[y]) reject("Population accounting failure."); if(y<Y){S_start[y+1]=S_end[y];U_start[y+1]=U_end[y];} }
}
model { log_lambda_bar~normal(log_foi_prior_mean,log_foi_prior_sd); target+=normal_lpdf(sigma_year|0,sigma_year_prior_scale)+log_sigma_year; z_year~normal(0,1); logit_q~normal(q_logit_prior_mean,q_logit_prior_sd); target+=gamma_lpdf(phi|2,.1)+log_phi; C~neg_binomial_2(expected_cases,phi); }
generated quantities { real q_implied=sum(C)/sum(X); int q_implied_gt_one=q_implied>1; real max_annual_multiplier=max(weight); real cumulative_infections=sum(X); array[Y] int C_pred; for(y in 1:Y) C_pred[y]=neg_binomial_2_rng(expected_cases[y],phi); }
