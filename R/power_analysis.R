pacman::p_load(dplyr, ggplot2, purrr, tidyr, simsalapar)

# Define the power analysis structure and output PDF of this structure
var_list <- varlist(
  n.sim = list(type = "N", expr = quote(N[sim]), value = 1000),
  n_years = list(type = "grid", expr = quote(n[yrs]), value = c(10, 20)),
  survey_interval = list(type = "grid", expr = quote(si[yrs]), value = c(1, 2)), 
  n_sites = list(type = "grid", expr = quote(n[sites]),value = c(50, 100, 200)),
  n_visits = list(type = "grid", expr = quote(n[visits]),value = c(2, 3)),
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
# (nrow(mkGrid(var_list)))

# Create PDF document of power analysis set-up
tmp <- tempfile(fileext = ".Rmd")
writeLines(c(c("---", "title: 'MABM power analysis parameters'", 
               "output: pdf_document", "tables: true", "---"),
             "$$\\dimendef\\prevdepth=0",
             paste0("$$", toLatex(var_list)[1]),
             toLatex(var_list)[-1]), tmp)
rmarkdown::render(tmp, 
                  output_file = "Output/power_analysis_setup.pdf")
message("Power analysis setup summarized in Output/power_analysis_setup.pdf")

# Load function to simulate the data and perform the power analysis
source("./R/fit_sims.R")

# Set simulation seeds for reproducible results
# seeds <- as.integer(runif(var_list$n.sim$value, -2e9, 2e9))
# all.equal(var_list$n.sim$value, length(unique(seeds)))
# saveRDS(seeds, file = "./Output/power_seeds.rds")
seeds <- readRDS("./Output/power_seeds.rds")
n_cores <- min(parallel::detectCores(), 24)
if (n_cores < 24) stop("This power analysis took nearly 5 days to complete ",
                       "on a machine using 24 cores. If you're OK with continuing, ",
                       "comment out lines 54-56 of R/power_analysis.R and carry on.")
res <- doClusterApply(var_list, cluster = parallel::makeCluster(n_cores), 
                      sfile = "./Output/MABM_power_sim.rds", 
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

# Create figures of power analysis and open it...
source("./R/powerplot.R")
powerplot(val, output = "pdf")

# Create figure of estimated population trend under various scenarios
# fix example to 100 sites, 10 years, annual surveys x2
trends <- simsalapar::array2df(sim_vals) %>%
  filter(parm %in% c("ann_r_p", "ann_r_est")) %>%
  spread(parm, value) %>%
  filter(n_sites == "50",
         n_years == "10",
         survey_interval == "1",
         n_visits == "2") %>%
  mutate(detected = ann_r_p <= alpha,
         wrong_sign = ifelse(detected, ann_r_est > 0, FALSE)) %>%
  group_by(spp, annual_r, wrong_sign) %>%
  mutate(prop_detected = sum(detected)/1000) %>%
            # min_dec = min(ann_r_est_decline, na.rm = TRUE),
            # max_dec = max(ann_r_est_decline, na.rm = TRUE),
            # min_inc = min(ann_r_est_increase, na.rm = TRUE),
            # max_inc = max(ann_r_est_increase, na.rm = TRUE)) %>%
  ungroup() %>%
  # Prettier organization for plots
  mutate(spp_label = factor(spp, levels = c("EPFU", "LABO", "MYLU", "PESU/NYHU"),
                            labels = c("EPFU\n", "LABO\n", "MYLU\n", "NYHU/\nPESU")),
         act_decline = round(as.numeric(as.character(annual_r)), 3),
         annual_r = factor(annual_r, labels = c("1.14% annual decline (25% over 25 years)",
                                                "2.73% annual decline (50% over 25 years)",
                                                "5% annual decline (~ 72% over 25 years)"))) %>%
  filter(detected)
    
pd <- position_dodge(width = 0.45, preserve = "single")
ggplot(trends, aes(spp_label, ann_r_est, color = wrong_sign, fill = prop_detected)) +
  geom_hline(aes(yintercept = act_decline), lty = "dashed", color = "gray50", lwd = 1) +
  geom_boxplot(position = pd) +
  facet_wrap(~ annual_r, ncol = 1) +
  scale_fill_viridis_c("Proportion of trends detected", end = 0.8) +
  scale_color_manual(values = c("black", "black")) + 
  xlab("Species") +
  scale_y_continuous("Estimated annual population change") +
  # ggtitle(paste0("MABM Power Analysis: ", pretty_r, "% annual change")) +
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5, barwidth = unit(2, "inches")),
         color = "none") +
  theme_bw() +
  theme(legend.position = "top")
ggsave("Output/MABM_trend_detection.png", width = 6.5, height = 9)
