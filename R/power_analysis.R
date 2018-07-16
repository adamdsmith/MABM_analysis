if (!requireNamespace("pacman", quietly = TRUE))
  install.packages("pacman")
pacman::p_load(dplyr, ggplot2, purrr, tidyr, simsalapar)
source("./R/powerplot.R")

var_list <- varlist(
  n.sim = list(type = "N", expr = quote(N[sim]), value = 250),
  n_years = list(type = "grid", expr = quote(n[yrs]), value = c(10, 20)),
  survey_interval = list(type = "grid", expr = quote(n[yrs]), value = c(1, 2)), 
  n_sites = list(type = "grid", expr = quote(n[sites]),value = c(50, 100, 200)),
  n_visits = list(type = "frozen", expr = quote(n[visits]),value = 2),
  # 25% and 50% declines in 25 years, resp, plus 'catastrophic' 5% annual decline
  annual_r = list(type = "grid", expr = quote(r[ann]), value = c(-0.011441, -0.027345, -0.05)),  
  spp = list(type = "grid", value = c("LABO", "PESU/NYHU", "EPFU", "MYLU")),
  base_C = list(type = "frozen", expr = quote(C[0]),
                value = list(LABO = 12, `PESU/NYHU` = 12, EPFU = 1, MYLU = 2)),
  negbin = list(type = "frozen", expr = quote(nb),
                value = list(LABO = TRUE, `PESU/NYHU` = TRUE, EPFU = TRUE, MYLU = TRUE)),
  theta = list(type = "frozen", 
               value = list(LABO = 2.8, `PESU/NYHU` = 3.6, EPFU = 2.2, MYLU = 2.2)),
  sigma_site = list(type = "frozen", expr = quote(sigma[site]),
                    value = list(LABO = 0.45, `PESU/NYHU` = 0.65, EPFU = 1.1, MYLU = 0.2)),
  sigma_r_site = list(type = "frozen", expr = quote(sigma[r]), value = 0.15),
  sigma_site_cor = list(type = "frozen", expr = quote(r[site[re]]),
                        value = list(LABO = -0.3, `PESU/NYHU` = -0.3, EPFU = 0.25, MYLU = -0.3)),
  sigma_yr = list(type = "frozen", expr = quote(sigma[yr]), 
                  value = list(LABO = 0.2, `PESU/NYHU` = 0.1, EPFU = 0.05, MYLU = 0.1)),
  # Could modify the following if there was interest in exploring failure to capture
  # important site/observation covariates explaining relative abundance
  site_cov_b = list(type = "frozen", expr = quote(b[site]), value = 0),
  obs_cov_b = list(type = "frozen", expr = quote(b[obs]), value = 0),
  omit_site_cov = list(type = "frozen", expr = quote(omit[b[site]]), value = FALSE),
  omit_obs_cov = list(type = "frozen", expr = quote(omit[b[obs]]), value = FALSE))

# Create PDF document of power analysis set-up
tmp <- tempfile(fileext = ".Rmd")
writeLines(c(c("---", "title: 'MABM power analysis parameters'", 
               "output: pdf_document", "tables: true", "---"),
             "$$\\dimendef\\prevdepth=0",
             paste0("$$", toLatex(var_list)[1]),
             toLatex(var_list)[-1]), tmp)
rmarkdown::render(tmp, 
                  output_file = "C:/Users/adsmith/Documents/FWS_Projects/MABM/Analysis/power_analysis_setup.pdf")

fit_sims <- function(n_sites, n_visits, n_years, annual_r, survey_interval, spp, 
                     base_C, negbin, theta, sigma_site, sigma_r_site, sigma_site_cor,
                     sigma_yr, site_cov_b, obs_cov_b, omit_site_cov, omit_obs_cov) {

  stopifnot(require(dplyr))
  stopifnot(require(glmmTMB))
  source("C:/Users/adsmith/Documents/FWS_Projects/MABM/Analysis/R/sim_MABM_C.R")

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

# Set simulation seeds for reproducible results
# seeds <- .Random.seed[seq(var_list$n.sim$value)]
# saveRDS(seeds, file = "./Output/power_seeds.rds")
seeds <- readRDS("./Output/power_seeds.rds")
res <- doClusterApply(var_list, cluster = parallel::makeCluster(4L), 
                      sfile = "./Output/power_results.rds", 
                      doOne = fit_sims, seed = seeds,
                      monitor = interactive())

val <- getArray(res)
# Check dimensions and names
dim(val)
# Give returned vector a name
names(dimnames(val))[1] <- "parm" 

# Retrieve error, warning, and timing information
err <- getArray(res, "error")
warn <- getArray(res, "warning")
time <- getArray(res, "time") # in milliseconds

# Check for warnings/errors
ftable(err, col.vars = c("spp", "n_years", "n_sites"), row.vars = c("survey_interval"))
ftable(warn, col.vars = c("spp", "n_years", "n_sites"), row.vars = c("survey_interval"))
# Convert timing to minutes
ftable(round(time/1000/60, 1), col.vars = c("n_years", "n_sites"), row.vars = c("spp", "survey_interval"))

