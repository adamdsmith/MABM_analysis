powerplot <- function(sim_vals,
                      alpha = 0.10,
                      ref_line = 0.80,
                      output = c("pdf", "png", "none")) {
  output <- match.arg(output)
  pdf <- identical(output, "pdf")
  png <- identical(output, "png")
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
           wrong_sign = ifelse(detected, ann_r_est > 0, FALSE),
           ann_r_est_decline = ifelse(detected & !wrong_sign, ann_r_est, NA),
           ann_r_est_increase = ifelse(detected & wrong_sign, ann_r_est, NA)) %>%
    group_by(spp, n_sites, n_years, survey_interval, n_visits, annual_r) %>%
    summarize(detections = sum(detected),
              n = n(),
              tmp = list(binCI(detections, n)),
              ann_r_est_decline = mean(ann_r_est_decline, na.rm = TRUE),
              ann_r_est_increase = mean(ann_r_est_increase, na.rm = TRUE),
              wrong_sign = sum(wrong_sign) / detections) %>%
    mutate(power = map_dbl(tmp, 1),
           lci = map_dbl(tmp, 2),
           uci = map_dbl(tmp, 3)) %>%
    dplyr::select(-detections, -n, -tmp) %>%
    ungroup() %>%
    # Prettier organization for plots
    mutate(spp_label = factor(spp, levels = c("EPFU", "LABO", "MYLU", "PESU/NYHU"),
                              labels = c("EPFU", "LABO/LASE", "MYLU", "NYHU/PESU")),
           n_years = factor(n_years, labels = c("10 years\nof surveys", "20 years\nof surveys")))
                 
  rs <- levels(simdf$annual_r)
  
  plots <- lapply(rs, function(r) {
    pretty_r <- round(as.numeric(r) * 100, 2)
    rdat <- filter(simdf, annual_r == r) %>%
      mutate(si_nv = factor(paste(survey_interval, n_visits, sep = ":"),
                            labels = c("Annual (2 visits)", "Annual (3 visits)",
                                       "Biennial (2 visits)", "Biennial (3 visits)")))
    pd <- position_dodge(width = 0.45)
    ggplot(rdat) +
      geom_hline(yintercept = ref_line, lty = "dashed", color = "gray50", lwd = ifelse(pdf, 1.5, 1)) +
      geom_linerange(aes(x = n_sites, ymin = lci, ymax = uci, group = si_nv), position = pd) +
      geom_point(aes(x = n_sites, y = power, fill = si_nv, shape = si_nv), position = pd,
                 size = ifelse(pdf, 2, 1.5)) +
      scale_fill_manual("Survey interval (# visits/yr)",
                        values = viridis::viridis(4, end = 0.8)) +
      scale_shape_manual("Survey interval (# visits/yr)", values = c(22, 23, 22, 23)) +
      facet_grid(n_years ~ spp_label) +
      xlab("# sites surveyed") +
      scale_y_continuous(bquote(Power~"("*alpha~"="~.(alpha)*")"),
                         limits = c(0, 1)) +
      ggtitle(paste0(pretty_r, "% annual change")) +
      guides(fill = guide_legend(nrow = 1, title.position = "top", title.hjust = 0.5), 
             shape = guide_legend(nrow = 1), title.position = "top", title.hjust = 0.5) +
      theme_bw() +
      theme(legend.position = "top",
            legend.key.width = unit(0.65, "lines"))
  })

  if (!requireNamespace("cowplot", quietly = TRUE))
    pacman::p_load("cowplot")
  legend <- cowplot::get_legend(plots[[1]])
  plots_png <- lapply(seq_along(plots), function(i) {
    p <- plots[[i]] + theme(legend.position = "none")
    if (i < 3) p <- p + theme(axis.text.x = element_blank(), axis.title.x=element_blank())
    p
  })
  plots_png <- c(list(legend), plots_png)
  cowplot::plot_grid(plotlist = plots_png, ncol = 1, rel_heights = c(0.2, 0.84, 0.84, 1),
                     labels = c("", LETTERS[1:3])) #label_y = c(1, 115, 1, 1))
  
  if (pdf) {
    ggsave("Output/FIG3.pdf", width = 190, height = 240, units = "mm")
    system(paste0('open "', normalizePath("Output/FIG3.pdf"), '"'))
  } else if (png) {
    ggsave("Output/FIG3.png", width = 6.5, height = 9)
  } 
  arrange(simdf, n_sites, n_years, survey_interval)
}
