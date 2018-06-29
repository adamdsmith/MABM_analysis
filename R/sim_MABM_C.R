# This function (and the power analysis that uses it) is making the following assumptions:
#   1. "true" (simulated) population change is linear
#   2. Detection probability is constant or modelled with relevant covariates
#   3. We include the most important covariates that influence site-level relative abundance (e.g., habitat).
#      This is relaxed a bit because we'll actually explore how missing important covariates influences power
#   4. We have reasonable estimates for:
#        - baseline detections per average length route
#        - detection probability (maybe not a strict requirement)
#        - site-level variation in baseline detections
#        - site-level variation in annual change estimates
#        - year-to-year variation in baseline detections
#        - the negative binomial dispersion parameter
#   5. Similar length MABM transects (~ 20 miles on average)

sim_MABM_C <- function(n_sites = 50,          # number of routes, assumed independent and with non-zero abundance in principle
                       n_visits = 2,          # number of replicates of each route per year
                       n_years = 10,          # survey period in years
                       delta_r = -0.01144,    # population change (proportion, from -0.99 to 0.99) to detect over n_years
                       survey_interval = 1,   # years between route surveys (e.g., yearly = 1, every third year = 3, etc.)
                       base_C = 2,            # average number of bats DETECTED (counted) per survey in year 0
                       negbin = TRUE,         # Negative binomial model (TRUE) or Poisson model (FALSE)?
                       theta = 2,             # overdispersion (NB size parameter); ignored if negbin = FALSE
                       site_cov_b = 0,        # coefficient (link scale) of scaled site-level covariate
                       obs_cov_b = 0,         # coefficient (link scale) of scaled observation-level covariate
                       sigma_site = 0.5,      # random intercept SD across sites in log(baseline bats DETECTED)
                       sigma_r_site = 0.15,   # random slope SD across sites in log(annual rate of change)  
                       sigma_yr = 0.1,        # random intercept SD across years in log(baseline bats DETECTED)
                       store_meta = FALSE) {  # Keep simulation parameters as metadata?
  
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
  
  # Observation-level (e.g., detection) relative abundance covariate...
  # We try to make it look like a sequential variable (e.g., day of year)
  obs_cov <- as.vector(replicate(n_sites * n_yrs_surv, sort(runif(n_visits, -30, 30))))
  obs_cov <- as.numeric(scale(obs_cov))
  
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

  # Calculate relative abundances
  log_lam <- log(base_C) + site_eta + yr_eta +   # Intercepts
    site_cov_b * site_cov +             # Site-level covariate
    obs_cov_b * obs_cov +               # Observation-level covariate
    log(1 + delta_r) * year_cov +       # Annual rate of change           
    site_r_eta * year_cov               # Random slope in annual rate of chanegmap
  h <- rgamma(n = n_sites * n_yrs_surv, shape = theta, rate = theta)
  C_P <- rpois(n = n_sites * n_yrs_surv * n_visits, lambda = exp(log_lam))
  C_NB <- rpois(n = n_sites * n_yrs_surv * n_visits, lambda = h * exp(log_lam))
  if (!negbin) C_vec <- C_P else C_vec <- C_NB
  
  if (store_meta)
    meta <- data.frame(variable = c("n_sites", "n_visits", "n_years",
                                    "survey_interval", "n_years_w_surveys",
                                    "delta_r", "base_C", "negbin", "theta",
                                    "sigma_site", "sigma_yr", "sigma_r_site"),
                       value = c(n_sites, n_visits, n_years,
                                 survey_interval = paste("every", survey_interval, "years"),
                                 n_yrs_surv,
                                 paste0(round(delta_r * 100, 3), "%"),
                                 base_C, negbin, 
                                 ifelse(negbin, theta, NA_character_),
                                 sigma_site, sigma_yr, sigma_r_site),
                       stringsAsFactors = FALSE)
  
  dat <- data.frame(C = C_vec,
                    visit = 1:n_visits,
                    year = year_cov,
                    yf = as.factor(year_cov),
                    site = rep(1:n_sites, each = n_yrs_surv * n_visits),
                    site_cov = site_cov,
                    obs_cov = obs_cov) %>%
    arrange(site, year, visit)

  if (!store_meta)
    return(dat)
  else
    return(list(meta = meta, data = dat))
}
