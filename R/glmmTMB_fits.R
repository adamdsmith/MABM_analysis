# Package requirements
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes", quiet = TRUE)
if (!requireNamespace("nrsmisc", quietly = TRUE)) remotes::install_github("adamdsmith/nrsmisc")
if (!requireNamespace("broom.mixed", quietly = TRUE)) remotes::install_github("bbolker/broom.mixed")
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman", quiet = TRUE)
pacman::p_load(dplyr, glmmTMB, bbmle, broom, broom.mixed, ggplot2, DHARMa)
options(scipen = 10, mc.cores = 4)
source("./R/utils.R")

# -------------------------- Prepare data for analysis ----------------------------
# These input files were created by the code in MABM_analysis_overview.Rmd
bats <- readRDS("./Output/MABM_spp_counts.rds")
site_meta <- readRDS("./Output/MABM_routes.rds") %>%
  dplyr::select(site, len_mi:lon)
site_habitat <- readRDS("./Output/MABM_site_env_data.rds") %>%
  dplyr::select(-contains("_1500"), -contains("_100"))
site_data <- left_join(site_meta, site_habitat, by = "site") %>%
  dplyr::select(site, len_mi, wood_250, urban_250)
survey_meta <- readRDS("./Output/MABM_survey_info.rds") %>%
  dplyr::select(site:surv_date, doy)
# Some MABM averages 
avg_doy <- mean(survey_meta$doy)
avg_wood <- mean(site_data$wood_250)
avg_urban <- mean(site_data$urban_250)

# ---------------- Join it with count and survey level data -----------------------
bats <- lapply(bats, function(i) {
  left_join(i, survey_meta, by = c("site", "surv_date")) %>%
    left_join(site_data, by = "site")
})

# --------------------------------- Pick a species --------------------------------
# Bat species of interest
boi <- c("LABO", "EPFU", "NYHU", "PESU", "MYLU")

# Set up data.frame to catch coefficients
coefs <- data.frame()

for (sp in boi) {
message(sp, " at the plate...")
# ------------------------ Some exploratory data analysis -------------------------
spp_dat <- bats[[sp]]
nrsmisc::pairs_plus(spp_dat[, 4:ncol(spp_dat)])

# ---------------------------- Some tweaks for analysis ---------------------------
spp_dat <- spp_dat %>%
  mutate(site = as.factor(site),
         # Created scale versions for covariate effect magnitude evaluation
         wood_std = as.numeric(scale(wood_250)),
         urban_std = as.numeric(scale(urban_250)),
         doy_std = as.numeric(scale(doy)),
         wood_250 = wood_250 / 0.1, # interpret as change in count per 10% increase
         urban_250 = urban_250 / 0.1, # ditto
         wk_jun1 = (doy - 152) / 7, # interpret as change in count per week since 1 June
         # offset; set so intercept interpretation is for mean route length
         l_len = log(len_mi/mean(site_data$len_mi)), 
         year = year - 2012)

# ------------------------ Set up some candidate models ---------------------------
# Need different RE structure for MYLU due to data sparsity
form <- count ~ year + wk_jun1 + wood_250 + urban_250 + offset(l_len) + (1|year) + (year|site)
if (sp == "MYLU")
  form <- update(form, . ~ . - (year|site) + (year - 1|site) + (1|site))

# Using interpretable versions of variables
m_pois <- glmmTMB(form, data = spp_dat, family = poisson)
m_nb1 <- glmmTMB(form, data = spp_dat, family = nbinom1)
m_nb2 <- glmmTMB(form, data = spp_dat, family = nbinom2)
mods <- list(m_pois, m_nb1, m_nb2); names(mods) <- c("Poisson", "NB1", "NB2")
print(aictab <- AICtab(mods, base = TRUE, weights = TRUE))
saveRDS(mods, file = file.path("./Output/Models", paste0(sp, "_glmmTMB.rds")))
best_mod <- names(mods)[which.min(sapply(mods, AIC))]
final_mod <- mods[[best_mod]]

# Time for some diagnostics
sim_out <- simulateResiduals(final_mod)
plot(sim_out)
Sys.sleep(3)
plotResiduals(spp_dat$year, sim_out$scaledResiduals)
Sys.sleep(3)
plotResiduals(spp_dat$wk_jun1, sim_out$scaledResiduals)
Sys.sleep(3)
plotResiduals(spp_dat$wood_250, sim_out$scaledResiduals)
Sys.sleep(3)
plotResiduals(spp_dat$urban_250, sim_out$scaledResiduals)
Sys.sleep(3)

# Get estimated # detections for average length route at MABM-wide averages in start year (2012)
nd <- expand.grid(year = 0, wk_jun1 = (avg_doy - 152) / 7, 
                  wood_250 = avg_wood / 0.1, urban_250 = avg_urban / 0.1, 
                  l_len = 0, site = "new")
fit <- predict(final_mod, newdata = nd, allow.new.levels = TRUE, type = "link", se = TRUE)

# Using scaled versions of variables
final_mod_sc <- update(final_mod, . ~ . - wood_250 - urban_250 - wk_jun1 + doy_std + wood_std + urban_std)
fixed_coef <- tidy(final_mod) %>%
  filter(effect == "fixed") %>%
  bind_rows(filter(tidy(final_mod_sc), grepl("doy_std|wood_std|urban_std", term))) %>%
  add_row(group = "NB", term = "theta", estimate = final_mod$sdr$par.fixed["betad"],
          std.error = sqrt(diag(final_mod$sdr$cov.fixed)["betad"])) %>%
  add_row(group = "fixed", term = "baseC", estimate = fit$fit,
          std.error = fit$se.fit) 
random_coef <- tidy(final_mod) %>%
  filter(effect == "ran_pars") %>%
  mutate(term = 
           ifelse(grepl("\\(Intercept\\)$", term), paste("sd", group, "int", sep = "_"),
                  ifelse(grepl("^cor_", term), paste("cor", group, "int", "slope", sep = "_"),
                         paste("sd", sub("^.*_", "", term), "slope", sep = "_"))))
all_coef <- bind_rows(fixed_coef, random_coef) %>%
  mutate(spp = sp, model = best_mod,
         lcl = estimate + qnorm(0.025) * std.error,
         hcl = estimate + qnorm(0.975) * std.error) %>%
  dplyr::select(spp, model, term, estimate, lcl, hcl) %>%
  filter(!grepl("Intercept", term))

coefs <- bind_rows(coefs, all_coef)

}

p <- ggplot(coefs, aes(x = term, y = estimate, color = spp)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(position = position_dodge(width = 0.6), size = 2) +
  geom_errorbar(aes(ymin = lcl, ymax = hcl), lwd = 1.25,
                 position = position_dodge(0.6), width = 0) + 
  xlab("Model parameter") + ylab("Coefficient (95% CI, if appropriate), on link scale") + 
  coord_flip() + 
  guides(colour = guide_legend("Species", reverse=TRUE)) +
  viridis::scale_color_viridis(discrete = TRUE) + 
  theme_black() +
  theme(axis.text.y = element_text(hjust = 1))
p
message("Coefficient plot saved at ./Output/MABM_analysis.pdf")
ggsave("./Output/MABM_analysis.pdf", width = 6.5, height = 9)
