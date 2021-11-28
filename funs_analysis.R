### ------------------------------------------------------------------------ ###
### wormplot of projection ####
### ------------------------------------------------------------------------ ###


plot_worm <- function(stk, stk_hist, refpts,
                      input, res,
                      history = TRUE, its = 1:5,
                      n_years = 20, yr_end = 2040, yr_start = 2000,
                      xintercept = 2020, ymax_catch = NA, ymax_rec = NA, 
                      ymax_ssb = NA, ymax_fbar = NA, 
                      title_rec = "Recruitment [1000s]",
                      title_catch = "Catch [1000t]",
                      title_ssb = "SSB [1000t]",
                      title_fbar = paste0("Mean F (ages ", 
                                          range(stk)[["minfbar"]], "-", 
                                          range(stk)[["maxfbar"]], ")")
) {
  #browser()
  ### load projection
  stk_plot <- stk_hist
  yrs_res <- dimnames(stk)$year
  stk_plot[, ac(yrs_res)] <- stk
  stk <- stk_plot
  stk <- window(stk, max(yr_start - 10, dims(stk)$minyear), end = yr_end)
  ### load reference points
  refpts <- iterMedians(refpts)
  ### get metrics
  qnts <- FLQuants(catch = catch(stk)/1000, rec = rec(stk)/1000,
                   ssb = ssb(stk)/1000, fbar = fbar(stk))
  ### percentiles
  qnts_perc <- lapply(qnts, quantile, probs = c(0.05, 0.25, 0.5, 0.75, 0.95),
                      na.rm = TRUE)
  qnts_perc <- FLQuants(qnts_perc)
  qnts_perc <- as.data.frame(qnts_perc)
  qnts_perc <- qnts_perc %>% select(year, iter, data, qname) %>%
    pivot_wider(names_from = iter, values_from = data)
  ### individual iterations
  qnts_iter <- as.data.frame(iter(qnts, its))
  ### plot
  p_catch <- ggplot() +
    geom_vline(xintercept = xintercept, colour = "grey", size = 0.5) +
    geom_ribbon(data = qnts_perc %>% filter(qname == "catch"),
                aes(x = year, ymin = `5%`, ymax = `95%`), alpha = 0.1,
                show.legend = FALSE) +
    geom_ribbon(data = qnts_perc %>% filter(qname == "catch"),
                aes(x = year, ymin = `25%`, ymax = `75%`), alpha = 0.1,
                show.legend = FALSE) +
    geom_line(data = qnts_iter %>% filter(qname == "catch"),
              aes(x = year, y = data, colour = iter), 
              size = 0.1, show.legend = FALSE) + 
    geom_line(data = qnts_perc %>% filter(qname == "catch"),
              aes(x = year, y = `50%`), size = 0.4) +
    geom_hline(yintercept = c(median(refpts["Cmsy"], na.rm = TRUE))/1000,
               colour = "black", size = 0.5, linetype = "dashed") +
    coord_cartesian(xlim = c(yr_start, yr_end), ylim = c(0, ymax_catch)) + 
    labs(y = title_catch) +
    theme_bw(base_size = 8) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank())
  p_rec <- ggplot() +
    geom_vline(xintercept = xintercept, colour = "grey", size = 0.5) +
    geom_ribbon(data = qnts_perc %>% filter(qname == "rec"),
                aes(x = year, ymin = `5%`, ymax = `95%`), alpha = 0.1,
                show.legend = FALSE) +
    geom_ribbon(data = qnts_perc %>% filter(qname == "rec"),
                aes(x = year, ymin = `25%`, ymax = `75%`), alpha = 0.1,
                show.legend = FALSE) +
    geom_line(data = qnts_iter %>% filter(qname == "rec"),
              aes(x = year, y = data, colour = iter), 
              size = 0.1, show.legend = FALSE) + 
    geom_line(data = qnts_perc %>% filter(qname == "rec"),
              aes(x = year, y = `50%`), size = 0.4) +
    coord_cartesian(xlim = c(yr_start, yr_end), ylim = c(0, ymax_rec)) + 
    labs(y = title_rec) +
    theme_bw(base_size = 8) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank())
  p_ssb <- ggplot() +
    geom_vline(xintercept = xintercept, colour = "grey", size = 0.5) +
    geom_ribbon(data = qnts_perc %>% filter(qname == "ssb"),
                aes(x = year, ymin = `5%`, ymax = `95%`), alpha = 0.1,
                show.legend = FALSE) +
    geom_ribbon(data = qnts_perc %>% filter(qname == "ssb"),
                aes(x = year, ymin = `25%`, ymax = `75%`), alpha = 0.1,
                show.legend = FALSE) +
    geom_line(data = qnts_iter %>% filter(qname == "ssb"),
              aes(x = year, y = data, colour = iter), 
              size = 0.1, show.legend = FALSE) + 
    geom_line(data = qnts_perc %>% filter(qname == "ssb"),
              aes(x = year, y = `50%`), size = 0.4) +
    geom_hline(yintercept = c(median(refpts["Bmsy"], na.rm = TRUE))/1000,
               colour = "black", size = 0.5, linetype = "dashed") +
    geom_hline(yintercept = c(median(refpts["Blim"], na.rm = TRUE))/1000,
               colour = "black", size = 0.5, linetype = "dotted") +
    coord_cartesian(xlim = c(yr_start, yr_end), ylim = c(0, ymax_ssb)) + 
    labs(y = title_ssb) +
    theme_bw(base_size = 8)
  p_fbar <- ggplot() +
    geom_vline(xintercept = xintercept, colour = "grey", size = 0.5) +
    geom_ribbon(data = qnts_perc %>% filter(qname == "fbar"),
                aes(x = year, ymin = `5%`, ymax = `95%`), alpha = 0.1,
                show.legend = FALSE) +
    geom_ribbon(data = qnts_perc %>% filter(qname == "fbar"),
                aes(x = year, ymin = `25%`, ymax = `75%`), alpha = 0.1,
                show.legend = FALSE) +
    geom_line(data = qnts_iter %>% filter(qname == "fbar"),
              aes(x = year, y = data, colour = iter), 
              size = 0.1, show.legend = FALSE) + 
    geom_line(data = qnts_perc %>% filter(qname == "fbar"),
              aes(x = year, y = `50%`), size = 0.4) +
    geom_hline(yintercept = c(median(refpts["Fmsy"], na.rm = TRUE)),
               colour = "black", size = 0.5, linetype = "dashed") +
    coord_cartesian(xlim = c(yr_start, yr_end), ylim = c(0, ymax_fbar)) + 
    labs(y = title_fbar) +
    theme_bw(base_size = 8)
  p <-  plot_grid(p_catch, p_rec, p_fbar, p_ssb, align = "v", 
                  rel_heights = c(1, 1.1))
  return(p)
}
