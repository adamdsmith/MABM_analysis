powerplot <- function(sim_vals,
                      alpha = 0.10,
                      ref_line = 0.80,
                      pdf = FALSE) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    install.packages("ggplot2")
  binCI <- function(x, n) {
    bt <- binom.test(x, n)
    out <- as.list(c(bt$estimate, bt$conf.int))
    names(out) <- c("power", "lci", "uci")
    out
  }
  simdf <- simsalapar::array2df(sim_vals) %>%
    filter(parm %in% c("ann_r_p", "ann_r_est")) %>%
    spread(parm, value) %>%
    mutate(detected = ann_r_p <= alpha,
           ann_r_est = ifelse(detected, ann_r_est, NA)) %>%
    group_by(n_years, survey_interval, n_sites, annual_r, spp) %>%
    summarize(detections = sum(detected),
                     n = n(),
                     tmp = list(binCI(detections, n)),
                     ann_r_est = mean(ann_r_est, na.rm = TRUE)) %>%
    mutate(power = map_dbl(tmp, 1),
           lci = map_dbl(tmp, 2),
           uci = map_dbl(tmp, 3)) %>%
    dplyr::select(-detections, -n, -tmp) %>%
    ungroup()
                 
  rs <- levels(simdf$annual_r)
  
  plots <- lapply(rs, function(r) {
    pretty_r <- round(as.numeric(r) * 100, 3)
    ggplot(filter(simdf, annual_r == r)) +
      geom_hline(yintercept = ref_line, lty = "dashed", color = "gray50", lwd = 1.5) +
      geom_pointrange(aes(x = spp, y = power, ymin = lci, ymax = uci, 
                          fill = survey_interval),
                      pch = 21,
                      position = position_dodge(width = 0.35)) +
      scale_fill_manual("Survey interval (yrs)",
                        values = c("#414487FF", "#7AD151FF")) +
      facet_grid(n_years ~ n_sites) +
      xlab("Species") +
      scale_y_continuous(bquote(Power~"("*alpha~"="~.(alpha)*")"),
                         limits = c(0, 1)) +
      ggtitle(paste0("MABM Power Analysis: ", pretty_r, "% annual change")) +
      theme_bw() +
      theme(legend.position = "top")
  })

  if (pdf) {
    tmp_fig <- tempfile(fileext = ".pdf")
    pdf(tmp_fig, width = 10, height = 7.5)
  }
  for (i in seq_along(plots)) print(plots[[i]])
  if (pdf) {
    invisible(dev.off())
    system(paste0('open "', normalizePath(tmp_fig), '"'))
  }
  arrange(simdf, n_sites, n_years, survey_interval)
}
