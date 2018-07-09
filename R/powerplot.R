powerplot <- function(sim_vals,
                      bat = c("LABO", "EPFU", "NYHU", "PESU", "MYLU"),
                      alpha = 0.10) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    install.packages("ggplot2")
  binCI <- function(x, n) {
    bt <- binom.test(x, n)
    out <- as.list(c(bt$estimate, bt$conf.int))
    names(out) <- c("power", "lci", "uci")
    out
  }
  simdf <- simsalapar::array2df(sim_vals) %>%
    filter(spp == bat,
                  parm %in% c("ann_r_p", "ann_r_est")) %>%
    spread(parm, value) %>%
    mutate(detected = ann_r_p <= alpha,
                  ann_r_est = ifelse(detected, ann_r_est, NA)) %>%
    group_by(n_years, survey_interval, n_sites, annual_r) %>%
    summarize(detections = sum(detected),
                     n = n(),
                     tmp = list(binCI(detections, n)),
                     ann_r_est = mean(ann_r_est, na.rm = TRUE)) %>%
    mutate(power = map_dbl(tmp, 1),
           lci = map_dbl(tmp, 2),
           uci = map_dbl(tmp, 3)) %>%
    dplyr::select(-detections, -n, -tmp) %>%
    ungroup()
                     

  p <- ggplot(simdf) +
    geom_hline(yintercept = 0.8, lty = "dashed", color = "gray50", lwd = 1.5) +
    geom_pointrange(aes(x = annual_r, y = power, ymin = lci, ymax = uci, 
                        fill = survey_interval),
                    pch = 21,
                    position = position_dodge(width = 0.35)) +
    viridis::scale_fill_viridis("Survey interval (yrs)", 
                                begin = 0.2, end = 0.8, discrete = TRUE) +
    facet_grid(n_years ~ n_sites) +
    xlab("Annual rate of population change") +
    ylab(bquote(Power~"("*alpha~"="~.(alpha)*")")) +
    ggtitle(paste("MABM Power Analysis:", bat)) +
    theme_bw() +
    theme(legend.position = "top")
  print(p)
  arrange(simdf, n_sites, n_years, survey_interval)
}
