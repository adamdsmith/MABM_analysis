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

# Set up list to catch AIC tables
aictabs <- vector(mode = "list", length = length(boi))
names(aictabs) <- boi

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
  print(aictabs[[sp]] <- AICtab(mods, base = TRUE, weights = TRUE, logLik = TRUE))
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
  saveRDS(final_mod, file = file.path("./Output/Models", paste0(sp, "_final_glmmTMB.rds")))
  saveRDS(final_mod_sc, file = file.path("./Output/Models", paste0(sp, "_final_glmmTMB_scaled.rds")))
  
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
saveRDS(coefs, file = "Output/MABM_GLMM_coefficients.rds")

# AIC table
aictabs <- lapply(aictabs, function(i) {
  class(i) <- "data.frame"
  i$model <- rownames(i)
  i
})
aictabs <- bind_rows(aictabs, .id = "Species") %>%
  mutate_if(is.numeric, round, digits = 3)
tmp <- tempfile(fileext = ".Rmd")
writeLines(c(c("---", "title: 'MABM analyis AIC model comparison'", 
               "output: pdf_document", "tables: true", "---"),
             kableExtra::kable(aictabs, "latex", booktabs = TRUE, digits = 3,
                               linesep = c('', '', '\\addlinespace'))),
           tmp)
rmarkdown::render(tmp, 
                  output_file = "~/FWS_Projects/MABM/Analysis/Output/MABM_analysis_AIC_tables.pdf")

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

# Create figure of standardized variables for manuscript
make_par_fig <- FALSE
if (make_par_fig) {
  # Rotate the error bar legend key to match figure
  GeomErrorbar$draw_key <-  function (data, params, size)     {
    draw_key_vpath <- function (data, params, size) {
      grid::segmentsGrob(0.5, 0.1, 0.5, 0.9, 
                   gp = grid::gpar(col = alpha(data$colour, data$alpha), 
                             lwd = data$size * .pt, lty = data$linetype, 
                             lineend = "butt"), arrow = params$arrow)
    }
    grid::grobTree(draw_key_vpath(data, params, size), 
             draw_key_point(transform(data, size = data$size), params))
  }
  coefs <- coefs %>%
    mutate(term_label = 
             factor(term,
                    levels = c("year", "wk_jun1", "doy_std", "wood_250", "wood_std", 
                               "urban_250", "urban_std"),
                    labels = c("Annual trend\n(linear)", 
                               "Survey date\n(weeks since 1 June)",
                               "Survey day of year\n(scaled)",
                               "Weighted forest cover\n(per 10% cover increment)",
                               "Weighted forest cover\n(sigma = 250 m)\n(scaled)",
                               "Weighted urban cover\n(per 10% cover increment)",
                               "Weighted urban cover\n(sigma = 250 m)\n(scaled)")))

  p <- ggplot(filter(coefs, term %in% c("year", "wk_jun1", "wood_250", "urban_250")),
              aes(x = spp, y = 100 * (exp(estimate) - 1))) +
    geom_hline(yintercept = 0, linetype = "dashed", lwd = 1, color = "gray50") +
    geom_point(position = position_dodge(width = 0.6), size = 3) +
    geom_errorbar(aes(ymin = 100 * (exp(lcl) - 1), 
                      ymax = 100 * (exp(hcl) - 1)), lwd = 1.5, width = 0) + 
    facet_wrap(~ term_label, nrow = 2, ncol = 2, scales = "free_y") +
    scale_y_continuous("Percent change (and 95% CI) in bat relative abundance") +
                       # limits = c(0, NA), minor_breaks = seq(-0.4, 0.9, by = 0.1)) +
    labs(x = "Species", color = NULL) + 
    
    # scale_color_viridis_d(end = 0.8) + 
    theme_bw() + theme(legend.position = "top")
  p
  ggsave("Output/MABM_scaled_parameter_estimates.png", width = 6.5, height = 4.5)
}

# Parameter estimate table
resp_scale <- function(x) exp(x) - 1
parms <- filter(coefs, term %in% c("year", "wk_jun1", "wood_250", "urban_250")) %>%
  select(-term_label) %>%
  mutate_if(is.numeric, resp_scale) %>%
  left_join(tibble(spp = c("EPFU", "LABO", "MYLU", "NYHU", "PESU")), .)
names(parms) <- c("Species", "Count model", "Parameter", "Estimate", "LCL", "UCL")
tmp <- tempfile(fileext = ".Rmd")
writeLines(c(c("---", "title: 'MABM parameter estimates (response scale)'", 
               "output: pdf_document", "tables: true", "---"),
             kableExtra::kable(parms, "latex", booktabs = TRUE, digits = 3,
                               linesep = c('', '', '', '\\addlinespace'))),
             tmp)
rmarkdown::render(tmp, 
                  output_file = "~/FWS_Projects/MABM/Analysis/Output/MABM_parameter_estimates.pdf")
message("Power analysis setup summarized in Output/MABM_parameter_estimates.pdf")

make_hab_fig <- FALSE
if (make_hab_fig) {
  pacman::p_load(cowplot)
  get_fits <- pbapply::pblapply(boi, function(sp) {
    mods <- readRDS(file = file.path("./Output/Models", paste0(sp, "_glmmTMB.rds")))
    best_mod <- names(mods)[which.min(sapply(mods, AIC))]
    final_mod <- mods[[best_mod]]

    # Get estimated # detections for varying main variables of interest across their range
    # while holding others at their MABM-wide averages in start year (2012)
    nd <- bind_rows(
      # Vary day of year first
      expand.grid(year = 0, wk_jun1 = (150:210 - 152) / 7, 
                  wood_250 = avg_wood / 0.1, urban_250 = avg_urban / 0.1,
                  l_len = 0, site = "new", var = "wk_jun1", stringsAsFactors = FALSE),
      # Then wooded cover
      expand.grid(year = 0, wk_jun1 = (avg_doy - 152) / 7, 
                  wood_250 = seq(0.1, 0.9, by = 0.01) / 0.1, urban_250 = avg_urban / 0.1,
                  l_len = 0, site = "new", var = "wood_250", stringsAsFactors = FALSE),
      # Then urban cover
      expand.grid(year = 0, wk_jun1 = (avg_doy - 152) / 7, 
                  wood_250 = avg_wood / 0.1, urban_250 = seq(0, 0.15, by = 0.005) / 0.1,
                  l_len = 0, site = "new", var = "urban_250", stringsAsFactors = FALSE)) %>%
      mutate(doy = wk_jun1 * 7 + 152,
             wood_250_orig = wood_250 / 10,
             urban_250_orig = urban_250 / 10)
    fit <- predict(final_mod, newdata = nd, allow.new.levels = TRUE, type = "link", se = TRUE)
    nd <- bind_cols(nd, fit) %>%
      mutate(spp = sp,
             lcl = exp(fit + qnorm(0.025) * se.fit),
             hcl = exp(fit + qnorm(0.975) * se.fit),
             fit = exp(fit))
  })
  get_fits <- bind_rows(get_fits)
  # df to fix scales of plots
  dummy <- get_fits %>% 
    group_by(spp, var) %>% 
    summarise(min = min(lcl), max = max(hcl)) %>% 
    tidyr::gather(metric, value, -spp, -var) %>%
    mutate(wk_jun1 = avg_doy - 1 + as.Date("2012-01-01"), 
           wood_250 = avg_wood, 
           urban_250 = avg_urban)
  p_doy <- ggplot(filter(get_fits, var == "wk_jun1"),
                  aes(x = doy - 1 + as.Date("2012-01-01"), y = fit)) + 
    geom_ribbon(aes(ymin = lcl, ymax = hcl), fill = "grey70") +
    geom_line() +
    geom_blank(data = dummy, aes(wk_jun1, value)) + 
    facet_grid(spp ~ ., scales = "free") + 
    theme_bw() + 
    theme(strip.background = element_blank(),
          strip.text = element_blank()) +
    labs(x = "Survey date", y = "Expected number of detections")
  p_wood <- ggplot(filter(get_fits, var == "wood_250"),
                  aes(x = wood_250_orig * 100, y = fit)) + 
    geom_ribbon(aes(ymin = lcl, ymax = hcl), fill = "grey70") +
    geom_line() +
    geom_blank(data = dummy, aes(wood_250 * 100, value)) + 
    facet_grid(spp ~ ., scales = "free") + 
    theme_bw() + 
    theme(strip.background = element_blank(),
          strip.text = element_blank()) +
    labs(x = "Weighted % forest cover", y = NULL)
  p_urban <- ggplot(filter(get_fits, var == "urban_250"),
                   aes(x = urban_250_orig * 100, y = fit)) + 
    geom_ribbon(aes(ymin = lcl, ymax = hcl), fill = "grey70") +
    geom_line() +
    geom_blank(data = dummy, aes(urban_250 * 100, value)) + 
    facet_grid(spp ~ ., scales = "free") + 
    theme_bw() + 
    labs(x = "Weighted % urban cover", y = NULL)
  plot_grid(p_doy, p_wood, p_urban, nrow = 1, rel_widths = c(1.07, 1.00, 1.09),
            labels = "AUTO", label_x = c(0.24, 0.14, 0.13), label_y = 0.99)
  ggsave("Output/MABM_habitat_associations.png", width = 8.5, height = 5.5, dpi = 300)
  
}
