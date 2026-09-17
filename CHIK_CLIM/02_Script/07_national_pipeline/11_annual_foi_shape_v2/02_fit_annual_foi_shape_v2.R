# Fit AR(1)-free normalized annual-FOI revision in its own output namespace.
required <- c("here", "rstan", "dplyr", "readr", "ggplot2", "patchwork")
if (length(m <- required[!vapply(required, requireNamespace, logical(1), quietly=TRUE)])) stop("Missing: ", paste(m, collapse=", "))
suppressPackageStartupMessages({library(rstan);library(dplyr);library(readr);library(ggplot2);library(patchwork)})
rstan_options(auto_write=TRUE); options(mc.cores=min(4L,parallel::detectCores()))
project_root_v2 <- function(){r<-here::here();if(file.exists(file.path(r,"CHIK_CLIM.Rproj")))r else file.path(r,"CHIK_CLIM")}
root <- project_root_v2(); source(file.path(root,"02_Script","90_development_archive","06_legacy_susceptibility","prepare_dynamic_annual_foi_data.R"))
env_int <- function(x,d){z<-Sys.getenv(x,"");if(!nzchar(z))d else as.integer(z)}

prior_shape_simulation <- function(prepared, sigma_scale, n=1500L) {
  set.seed(20260913 + as.integer(100*sigma_scale)); y<-nrow(prepared$annual)
  sims<-lapply(seq_len(n),function(i){eta<-abs(rnorm(1,0,sigma_scale))*rnorm(y);w<-exp(eta-(log(sum(exp(eta)))-log(y)));lambda<-exp(rnorm(1,prepared$foi_anchor$log_mean,prepared$foi_anchor$log_sd))*w;tibble::tibble(draw=i,year=prepared$annual$year,lambda=lambda,weight=w,attack_prob=-expm1(-lambda))})|>bind_rows()
  sims |> group_by(year) |> summarise(lambda_mid=median(lambda),lambda_lo=quantile(lambda,.025),lambda_hi=quantile(lambda,.975),attack_mid=median(attack_prob),max_multiplier=median(max(weight)),.groups="drop")
}

fit_annual_foi_shape_v2 <- function(){
  prepared<-build_dynamic_annual_foi_data(write_outputs=FALSE); out_t<-file.path(root,"03_Output", "07_national_pipeline", "tables", "annual_foi_shape_v2");out_f<-file.path(root,"03_Output", "07_national_pipeline", "figures", "annual_foi_shape_v2");dir.create(out_t,recursive=TRUE,showWarnings=FALSE);dir.create(out_f,recursive=TRUE,showWarnings=FALSE)
  prior<-bind_rows(lapply(c(.5,.75,1),function(s)prior_shape_simulation(prepared,s)|>mutate(sigma_prior_scale=s)))
  write_csv(prior,file.path(out_t,"annual_foi_shape_v2_prior_predictive_summary.csv"))
  pprior<-ggplot(prior,aes(year,lambda_mid,colour=factor(sigma_prior_scale),fill=factor(sigma_prior_scale)))+geom_ribbon(aes(ymin=lambda_lo,ymax=lambda_hi),alpha=.16,colour=NA)+geom_line()+scale_y_continuous(trans="log10")+labs(title="Prior predictive normalized annual FOI",colour="sigma_year prior SD",fill="sigma_year prior SD",y="Annual FOI (log scale)")+theme_classic()
  ggsave(file.path(out_f,"annual_foi_shape_v2_prior_predictive.pdf"),pprior,width=190,height=120,units="mm",device=cairo_pdf)
  q1<-dplyr::filter(prepared$q_scenarios,scenario=="Q1_primary_product_anchor")
  configs<-tibble::tribble(~fit_id,~model,~sigma_scale,~foi_sd_multiplier,
    "M1v2_SHAPE_s0_5","SHAPE",.5,1,"M1v2_SHAPE_s0_75","SHAPE",.75,1,"M1v2_SHAPE_s1_0","SHAPE",1,1,
    "M1v2_SHAPE_foisd1_5","SHAPE",.75,1.5,"M1v2_SHAPE_foisd2","SHAPE",.75,2,"M1v2_NB_s0_75","NB",.75,1)
  selected<-strsplit(Sys.getenv("ANNUAL_FOI_SHAPE_V2_MODELS",paste(configs$fit_id,collapse=",")),",",fixed=TRUE)[[1]];configs<-dplyr::filter(configs,fit_id%in%selected)
  shape_model<-stan_model(file.path(root,"02_Script", "00_shared", "stan", "current", "national", "chik_dynamic_annual_foi_shape_v2.stan"));nb_model<-stan_model(file.path(root,"02_Script", "00_shared", "stan", "current", "national", "chik_dynamic_annual_foi_nb_v2.stan"))
  iter<-env_int("ANNUAL_FOI_SHAPE_V2_ITER",2000L);warmup<-env_int("ANNUAL_FOI_SHAPE_V2_WARMUP",1000L);chains<-env_int("ANNUAL_FOI_SHAPE_V2_CHAINS",4L); manifest<-vector("list",nrow(configs))
  for(i in seq_len(nrow(configs))){d<-prepared$stan_base_data;d$log_foi_prior_sd<-d$log_foi_prior_sd*configs$foi_sd_multiplier[i];d$sigma_year_prior_scale<-configs$sigma_scale[i]; if(configs$model[i]=="SHAPE"){d$C_total<-sum(d$C);model<-shape_model}else{d$q_logit_prior_mean<-q1$q_logit_prior_mean;d$q_logit_prior_sd<-q1$q_logit_prior_sd;model<-nb_model};message("[fit-v2] ",configs$fit_id[i]);start<-Sys.time();fit<-sampling(model,data=d,chains=chains,iter=iter,warmup=warmup,seed=20260920+i,control=list(adapt_delta=.99,max_treedepth=12),refresh=max(1,floor(iter/10)));elapsed<-as.numeric(difftime(Sys.time(),start,units="secs"));path<-file.path(root,"03_Output","07_national_pipeline","model_fits","annual_foi_shape_v2",paste0("annual_foi_shape_v2_",configs$fit_id[i],".rds"));saveRDS(list(fit=fit,prepared=prepared,config=as.list(configs[i,]),elapsed_seconds=elapsed),path);manifest[[i]]<-mutate(configs[i,],elapsed_seconds=elapsed,path=path)}
  write_csv(bind_rows(manifest),file.path(out_t,"annual_foi_shape_v2_fit_manifest.csv"))
}
if(sys.nframe()==0L)fit_annual_foi_shape_v2()
