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

# Create graphs of power analysis and open it...
source("./R/powerplot.R")
powerplot(val, output = "pdf")
