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
                          sigma_site_cor = -0.2, # correlation between random intercept and slope
                          sigma_yr = 0.15) {     # random intercept SD across years in log(baseline bats DETECTED)

  # Calculate baseline N (actual abund) from C (detected) and detection probability
  base_N <- base_C / base_p

  ### Year-related effects
  year_cov <- seq(n_years) - 1
  
  # Filter data to years surveyed
  if (survey_interval > 1) {
    keep_yrs <- seq(0, n_years - 1, by = survey_interval)
    year_cov <- year_cov[year_cov %in% keep_yrs]
  }
  n_yrs_surv <- n_distinct(year_cov)
  
  # Random intercept for each year
  yr_eta <- rnorm(n_yrs_surv, 0, sigma_yr)
  year_df <- data.frame(year = year_cov, yr_eta)
  ###

  ### With years to be surveyed decided, create base data.frame
  dat <- expand.grid(visit = seq(n_visits),
                     site = seq(n_sites),
                     year = year_cov)
  
  ### Site-level effects
  # Site-level relative abundance covariate...
  site_cov <- runif(n_sites)
  
  # Potentially correlated random site intercept and annual change slope
  re_covar <- sigma_site_cor * sigma_site * sigma_r_site
  site_vcov <- matrix(c(sigma_site^2, rep(re_covar, 2), sigma_r_site^2), 2, 2)
  site_re <- MASS::mvrnorm(n_sites, mu = c(0,0), Sigma = site_vcov)
  
  site_df <- data.frame(site = seq(n_sites), site_cov,
                        site_eta = site_re[, 1], site_r_eta = site_re[, 2])
  ###

  ### Observation-level (e.g., detection) effects
  # Make it appear as a sequential variable (e.g., day of year)
  obs_cov <- as.vector(replicate(n_sites * n_yrs_surv, sort(runif(n_visits))))
  ###

  # Paste it all together, scale relevant variables, and calculate linear predictor
  # of for (latent) actual abundance...
  scale <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  dat <- left_join(dat, site_df, by = "site") %>%
    left_join(year_df, by = "year") %>%
    cbind(obs_cov) %>%
    mutate_at(vars(c("site_cov", "obs_cov")), scale) %>%
    mutate(lp_N = log(base_N) + site_eta + yr_eta + # Intercepts
             site_cov_b * site_cov +              # Site-level covariate
             log(1 + delta_r) * year +           # Annual rate of change
             site_r_eta * year)                  # Random slope in annual rate of change
  
  # Add relative abundances
  h <- rgamma(n = n_sites * n_yrs_surv, shape = theta, rate = theta)
  N_P <- rpois(n = n_sites * n_yrs_surv * n_visits, lambda = exp(dat$lp_N))
  N_NB <- rpois(n = n_sites * n_yrs_surv * n_visits, lambda = h * exp(dat$lp_N))
  if (!negbin) N <- N_P else N <- N_NB
  dat$N <- N

  # Calculate detection probability from the observation model,
  # and simulate observed count of individuals
  dat <- mutate(dat,
                p = plogis(qlogis(base_p) + obs_cov_b * obs_cov)) %>%
    rowwise() %>%
    mutate(C = rbinom(1, N, p)) %>%
    ungroup()

  dat <- dat %>%
    select(site, year, visit, N, p, C, site_cov, obs_cov) %>%
    arrange(site, year, visit)
  dat
}