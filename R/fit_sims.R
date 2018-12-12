fit_sims <- function(n_sites, n_visits, n_years, annual_r, survey_interval, spp, 
                     base_C, negbin, theta, sigma_site, sigma_r_site, sigma_site_cor,
                     sigma_yr, site_cov_b, obs_cov_b, omit_site_cov, omit_obs_cov) {
  
  stopifnot(require(dplyr))
  stopifnot(require(glmmTMB))
  # Need full path here since we'll be using multiple cores
  # You may need to modify for your system...
  source(normalizePath("~/FWS_Projects/MABM/Power_analysis/R/sim_MABM_C.R"))
  
  # Get species level settings
  bC <- base_C[[spp]]
  nb <- negbin[[spp]]
  th <- theta[[spp]]
  sigsite <- sigma_site[[spp]]
  sigsitecor <- sigma_site_cor[[spp]]
  sigyr <- sigma_yr[[spp]]
  
  dat <- sim_MABM_C(n_sites, n_visits, n_years, delta_r = annual_r, survey_interval, 
                    base_C = bC, negbin = nb, theta = th, site_cov_b, obs_cov_b, 
                    sigma_site = sigsite, sigma_r_site, sigma_site_cor = sigsitecor,
                    sigma_yr = sigyr)
  if (nb) fm <- nbinom2 else fm <- poisson
  
  os <- omit_site_cov
  oo <- omit_obs_cov
  
  # Parameter estimates of interest...
  all_parms <- c("base_C_est", "ann_r_est", "site_cov_est", "obs_cov_est", 
                 "theta_est", "site_sd_est", "site_r_sd_est", 
                 "site_re_cor", "yr_sd_est", "ann_r_p")
  
  # Modify base formula, if necessary
  form <- C ~ year + site_cov + obs_cov + (1|year) + (year|site)
  if (os) {
    if (oo) form <- update.formula(form, . ~ . - site_cov - obs_cov)
    else form <- update.formula(form, . ~ . - site_cov)
  } else {
    if (oo) form <- update.formula(form, . ~ . - obs_cov)
  }
  
  # Catch warnings to check for convergence/Hessian matrix issues (overspecification) 
  m <- tryCatch(glmmTMB(form, data = dat, family = nbinom2),
                warning = function(w) w)
  # Check for over-specified model
  if (is(m, "warning")) {
    simp_re <- grepl("Model convergence problem", m$message)
    if (!simp_re) warning("Non-convergence/Hessian-related warning")
  } else simp_re <- FALSE
  if (simp_re) {
    form <- update(form, . ~ . - (year|site) + (year - 1|site) + (1|site))
    m <- glmmTMB(form, data = dat, family = nbinom2)
  }
  
  parm_nms <- c("base_C_est", "ann_r_est",
                if (os) NULL else "site_cov_est",
                if (oo) NULL else "obs_cov_est",
                if (nb) "theta_est" else NULL,
                "site_sd_est", "site_r_sd_est", 
                if (!simp_re) "DISCARD",
                "yr_sd_est")
  parms <- m$sdr$par.fixed
  names(parms) <- parm_nms 
  sds <- sqrt(diag(m$sdr$cov.fixed)); names(sds) <- parm_nms
  parms_sd <- sqrt(diag(m$sdr$cov.fixed))
  ann_r_se <- sds["ann_r_est"]
  ann_r_p <- 2 * pnorm(abs(parms["ann_r_est"]/ann_r_se), lower.tail = FALSE)
  to_exp <- which(!grepl("site_cov|obs_cov", parm_nms))
  parms[to_exp] <- exp(parms[to_exp])
  parms["ann_r_est"] <- parms["ann_r_est"] - 1
  parms["ann_r_p"] <- ann_r_p
  # Hoops to get correlation among REs
  if (simp_re)
    parms["site_re_cor"] <- NA
  else
    parms["site_re_cor"] <- m$obj$env$report(m$fit$parfull)[["corr"]][[1]][2]
  if (os) parms["site_cov_est"] <- NA
  if (oo) parms["obs_cov_est"] <- NA
  # Reorder to standardize
  parms[all_parms]
}
