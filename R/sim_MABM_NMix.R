# This function (and the power analysis that uses it) is making the following assumptions:
#   1. "true" (simulated) population change is linear
#   2. Detection probability is constant or modelled with relevant covariates
#   3. We include the most important covariates that influence site-level relative abundance (e.g., habitat)
#   4. We have reasonable estimates for:
#        - baseline detections per average length route
#        - detection probability (maybe not a strict requirement)
#        - site-level variation in baseline detections
#        - site-level variation in annual change estimates
#        - year-to-year variation in baseline detections
#        - the negative binomial dispersion parameter
#   5. Similar length MABM transects (~ 20 miles on average)

sim_MABM_NMix <- function(n_sites = 50,          # number of routes, assumed independent and with non-zero abundance in principle
                          n_visits = 2,          # number of replicates of each route per year
                          n_years = 10,          # survey period in years
                          delta_r = -0.01144,    # population change (proportion, from -0.99 to 0.99) to detect over n_years
                          survey_interval = 1,   # years between route surveys (e.g., yearly = 1, every third year = 3, etc.)
                          base_C = 2,            # average number of bats DETECTED (counted) per survey in year 0
                          base_p = 0.10,         # assumed (or estimated) average detection probability
                          negbin = TRUE,         # Negative binomial model (TRUE) or Poisson model (FALSE)?
                          theta = 2,             # overdispersion (NB size parameter); ignored if negbin = FALSE
                          site_cov_b = 0,        # coefficient (link scale) of scaled site-level covariate
                          obs_cov_b = 0,         # coefficient (link scale) of scaled observation-level covariate
                          sigma_site = 0.5,      # random intercept SD across sites in log(baseline bats DETECTED)
                          sigma_r_site = 0.15,   # random slope SD across sites in log(annual rate of change)  
                          sigma_yr = 0.15,       # random intercept SD across years in log(baseline bats DETECTED)
                          store_meta = FALSE) {  # Keep simulation parameters as metadata?
  
  # Calculate baseline N (actual abund) from C (detected) and detection probability
  base_N <- base_C / base_p

  # Variable indicating years from initial year
  year_cov <- rep(seq_len(n_years) - 1, each = n_visits)
  
  # Filter data to years surveyed
  if (survey_interval > 1) {
    keep_yrs <- seq(1, n_years, by = survey_interval)
    year_cov <- year_cov[year_cov %in% (keep_yrs - 1)]
  }
  n_yrs_surv <- n_distinct(year_cov)

  # Site-level relative abundance covariate...
  site_cov <- as.numeric(scale(rep(rnorm(n_sites), each = n_yrs_surv * n_visits)))

  ### THIS WILL BE INCORPORATED INTO EFFECT ON P BELOW
  # Observation-level (e.g., detection) covariate...
  # We try to make it look like a sequential variable (e.g., day of year)
  obs_cov <- as.vector(replicate(n_sites * n_yrs_surv, sort(runif(n_visits, -30, 30))))
  obs_cov <- obs_cov_v <- as.numeric(scale(obs_cov))
  obs_cov <- matrix(obs_cov, ncol = 2, byrow = TRUE)
  obs_cov <- aperm(`dim<-`(t(obs_cov), list(n_visits, n_sites, n_yrs_surv)), c(2,1,3))

  # Random site intercept and annual change slope
  # Assumed independent
  site_eta <- rnorm(n_sites, 0, sigma_site)
  site_eta <- rep(as.numeric(scale(site_eta, scale = FALSE)), 
                  each = n_yrs_surv * n_visits)
  site_r_eta <- rnorm(n_sites, 0, sigma_r_site)
  site_r_eta <- rep(as.numeric(scale(site_r_eta, scale = FALSE)), 
                    each = n_yrs_surv * n_visits)
  yr_eta <- rnorm(n_yrs_surv, 0, sigma_yr)
  yr_eta <- as.numeric(scale(rep(rep(yr_eta, each = n_visits), n_sites), scale = FALSE))
  
  # Calculate site x year "actual" abundances
  log_lam <- log(base_N) + site_eta + yr_eta +   # Intercepts
    site_cov_b * site_cov +                      # Site-level covariate
    log(1 + delta_r) * year_cov +                # Annual rate of change           
    site_r_eta * year_cov                        # Random slope in annual rate of change
  h <- rgamma(n = n_sites * n_yrs_surv, shape = theta, rate = theta)
  N_P <- rpois(n = n_sites * n_yrs_surv, lambda = exp(log_lam))
  N_NB <- rpois(n = n_sites * n_yrs_surv, lambda = h * exp(log_lam))
  if (!negbin) N_vec <- N_P else N_vec <- N_NB

  # Populate actual abundance array
  N <- array(N_vec, dim = c(n_sites, n_yrs_surv))

  # Populate the count array and the logit detection probability array
  C <- array(NA_integer_, dim = c(n_sites, n_visits, n_yrs_surv))
  lp <- qlogis(base_p) + obs_cov_b * obs_cov 

  for (i in 1:n_sites){
    for (k in 1:n_yrs_surv){
      # If assuming constant p, could have estimated all j without loop
      # C[i, ,k] <- rbinom(n_visits, size = N[i,k], prob = base_p)
      for (j in 1:n_visits){ 
        C[i,j,k] <- rbinom(1, size = N[i,k], prob = plogis(lp[i,j,k]))
  }}}

  if (store_meta)
    meta <- data.frame(variable = c("n_sites", "n_visits", "n_years",
                                    "survey_interval", "n_years_w_surveys",
                                    "delta_r", "base_C", "base_p", "negbin", "theta",
                                    "sigma_site", "sigma_yr", "sigma_r_site"),
                       value = c(n_sites, n_visits, n_years,
                                 survey_interval = paste("every", survey_interval, "years"),
                                 n_yrs_surv,
                                 paste0(round(delta_r * 100, 3), "%"),
                                 base_C, base_p, negbin, 
                                 ifelse(negbin, theta, NA_character_),
                                 sigma_site, sigma_yr, sigma_r_site),
                       stringsAsFactors = FALSE)

  dat <- data.frame(C = as.vector(C),
                    visit = 1:n_visits,
                    year = year_cov,
                    yf = as.factor(year_cov),
                    site = rep(1:n_sites, each = n_yrs_surv * n_visits),
                    site_cov = site_cov,
                    obs_cov = obs_cov_v) %>%
    arrange(site, year, visit)
  
  if (!store_meta)
    return(dat)
  else
    return(list(meta = meta, data = dat))
}