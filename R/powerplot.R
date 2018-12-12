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
    filter(!is.na(ann_r_p)) %>%
    mutate(detected = ann_r_p <= alpha,
           ann_r_est = ifelse(detected, ann_r_est, NA)) %>%
    group_by(spp, n_sites, n_years, survey_interval, n_visits, annual_r) %>%
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
    rdat <- filter(simdf, annual_r == r) %>%
      mutate(si_nv = factor(paste(survey_interval, n_visits, sep = ":"),
                            labels = c("Annual (2 visits)", "Annual (3 visits)",
                                       "Biennial (2 visits)", "Biennial (3 visits)")))
    pd <- position_dodge(width = 0.45)
    ggplot(rdat) +
      geom_hline(yintercept = ref_line, lty = "dashed", color = "gray50", lwd = 1.5) +
      geom_linerange(aes(x = spp, ymin = lci, ymax = uci, group = si_nv), position = pd) +
      geom_point(aes(x = spp, y = power, fill = si_nv, shape = si_nv), position = pd) +
      # geom_pointrange(aes(x = spp, y = power, ymin = lci, ymax = uci, 
      #                     fill = interaction(n_visits, survey_interval)),
      #                 pch = 21,
      #                 position = pd)) +
      scale_fill_manual("Survey interval (# visits/yr)",
                        values = c("#d8b365", "#8c510a", "#5ab4ac", "#01665e")) +
      scale_shape_manual("Survey interval (# visits/yr)", values = c(22, 23, 22, 23)) +
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
