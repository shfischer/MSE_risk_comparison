### ------------------------------------------------------------------------ ###
### ple.27.7e - analyse MSE results ####
### ------------------------------------------------------------------------ ###
library(mse)
library(GA)
library(tidyr)
library(dplyr)
library(cowplot)
library(ggplot2)
library(foreach)
source("funs.R")
source("funs_GA.R")

### ------------------------------------------------------------------------ ###
### collate results ####
### ------------------------------------------------------------------------ ###
OM = c("baseline", "M_low", "M_high", "M_Gislason", "M_dd",
       "M_no_migration", "no_discards", "rec_no_AC", 
       "rec_failure")
res <-   foreach(MP = c("rfb", "2over3", "2over3_XSA", "hr", "ICES_SAM")) %:%
  foreach(stock = c("ple.27.7e", "cod.27.47d20"), .combine = bind_rows) %:%
  foreach(OM = c("baseline"), .combine = bind_rows) %:%
  foreach(optimised = c("default", "multiplier", "all"), 
          .combine = bind_rows) %:%
  foreach(period = c("1-20", "11-20"), .combine = bind_rows) %do% {
    #browser()
    res_i <- data.frame(stock = stock, OM = OM, MP = MP, optimised = optimised, 
                        period = period)
    if (MP %in% c("rfb", "hr")) {
      ga_path <- paste0("output/", stock, "/", OM, "/1000_20/multiplier/", MP, 
                        "/")
      if (isTRUE(MP %in% c("rfb", "hr")) & 
          optimised %in% c("default", "multiplier")) {
        ga_prefix <- "multiplier-upper_constraint1.2-lower_constraint0.7"
      } else if (identical(MP, "rfb") & identical(optimised, "all")) {
        ga_prefix <- paste0("lag_idx-range_idx_1-range_idx_2-exp_r-exp_f-exp_b-",
                            "interval-multiplier-upper_constraint1.2-lower_",
                            "constraint0.7")
      } else if (identical(MP, "hr") & identical(optimised, "all")) {
        ga_prefix <- paste0("idxB_lag-idxB_range_3-comp_b_multiplier-interval-",
                            "multiplier-upper_constraint1.2-lower_constraint0.7")
      }
      ga_res_file <- paste0(ga_prefix, "--obj_ICES_res_", period, ".rds")
      ga_runs_file <- paste0(ga_prefix, "--obj_ICES_runs.rds")
      
      ### load file
      if (!file.exists(paste0(ga_path, ga_res_file))) return(NULL)
      ga_res <- readRDS(paste0(ga_path, ga_res_file))
      ga_solution <- as.data.frame(as.list(ga_res@solution[1, ]))
      ### round parameters to significance used in optimisation
      if (identical(MP, "rfb")) {
        ga_solution[c(1:4, 8)] <- round(ga_solution[c(1:4, 8)], 0)
        ga_solution[c(5:7)] <- round(ga_solution[c(5:7)], 1)
        ga_solution[c(9:11)] <- round(ga_solution[c(9:11)], 2)
      } else if (identical(MP, "hr")) {
        ga_solution[c(1:2, 5)] <- round(ga_solution[c(1:2, 5)], 0)
        ga_solution[c(3:4)] <- round(ga_solution[c(3:4)], 1)
        ga_solution[c(6:8)] <- round(ga_solution[c(6:8)], 2)
      }
      if (identical(optimised, "default") & identical(MP, "rfb")) {
        ga_solution$multiplier <- 0.95
      } else if (identical(optimised, "default") & identical(MP, "hr")) {
        ga_solution$multiplier <- 1
      }
      
      ### load stats
      ga_runs <- readRDS(paste0(ga_path, ga_runs_file))
      ga_run_i <- ga_runs[[paste0(ga_solution, collapse = "_")]]
      ga_stats_i <- bind_rows(ga_run_i$stats)
      
      res_i$generations <- ga_res@iter
      res_i$fitness <- ga_res@fitnessValue
      if (identical(optimised, "default")) {
        res_i$generations <- NA
        res_i$fitness <- NA
      }
      ### combine definition and stats
      res_i <- bind_cols(res_i, ga_solution, ga_stats_i)
      ### check if MP results exist
      path_i <- paste0("output/", stock, "/", OM, "/1000_20/", MP, "/")
      file_i <- paste0("mp_", paste0(ga_solution, collapse = "_"), ".rds")
      file_i <- paste0(path_i, file_i)
      res_i$file <- ifelse(file.exists(file_i), file_i, NA)
    ### other MPs - only default exists
    } else if (identical(optimised, "default")) {
      ### load stats
      path_i <- paste0("output/", stock, "/", OM, "/1000_20/", MP, "/")
      file_i <- paste0(path_i, "stats.rds")
      if (!file.exists(file_i)) return(NULL)
      stats_i <- readRDS(file_i)
      stats_i <- bind_rows(stats_i)
      ### check if MP results exist
      path_mp_i <- paste0(path_i, "mp.rds")
      stats_i$file <- ifelse(file.exists(path_mp_i), path_mp_i, NA)
      res_i <- bind_cols(res_i, stats_i)
    } else {
      return(NULL)
    }
    ### calculate fitness
    if (isTRUE(is.na(res_i$fitness)) | is.null(res_i$fitness)) {
      tmp_catch <- unlist(res_i[, paste0(gsub(period, pattern = "-",
                                              replacement = ":"), 
                                         "_Catch_rel") ])
      tmp_risk <- unlist(res_i[, paste0(gsub(period, pattern = "-",
                                             replacement = ":"),
                                        "_risk_Blim_max") ])
      res_i$fitness <- tmp_catch - penalty(x = tmp_risk, 
                                           negative = FALSE, max = 1, 
                                           inflection = 0.06, 
                                           steepness = 1000)
    }
    return(res_i)
}
res <- do.call(bind_rows, res)
res <- res[, c(1:14, 92:94, 15:91)]

View(res)
res$file

write.csv(res, file = "output/MPs_baseline.csv", row.names = FALSE)
saveRDS(res, file = "output/MPs_baseline.rds")

### ------------------------------------------------------------------------ ###
### rfb-rule - multiplier ####
### ------------------------------------------------------------------------ ###

### GA summary
ga_mult <- readRDS(paste0("output/ple.27.7e/baseline/1000_20/multiplier/rfb/",
                          "multiplier-upper_constraint1.2-lower_constraint0.7",
                          "--obj_ICES_res_1-20.rds"))
ga_mult@solution
ga_mult@fitnessValue

### all runs
runs_mult <- readRDS(paste0("output/ple.27.7e/baseline/1000_20/multiplier/rfb/",
                            "multiplier-upper_constraint1.2-lower_constraint0.7",
                            "--obj_ICES_runs.rds"))
runs_mult <- lapply(runs_mult, function(x) {
  #browser()
  c(x$pars, x$stats)
})
runs_mult <- do.call(rbind, runs_mult)
row.names(runs_mult) <- NULL
runs_mult <- as.data.frame(runs_mult)

names(runs_mult) <- gsub(x = names(runs_mult), pattern = "1:5_", "1:5.")
names(runs_mult) <- gsub(x = names(runs_mult), pattern = "1:10_", "1:10.")
names(runs_mult) <- gsub(x = names(runs_mult), pattern = "6:10_", "6:10.")
names(runs_mult) <- gsub(x = names(runs_mult), pattern = "11:20_", "11:20.")
names(runs_mult) <- gsub(x = names(runs_mult), pattern = "1:20_", "1:20.")
runs_mult2 <- runs_mult %>%
  select(lag_idx:lower_constraint, "1:5.risk_Blim":"1:20.ICV") %>%
  pivot_longer(cols = "1:5.risk_Blim":"1:20.ICV",
               names_sep = "\\.", names_to = c("years", "stat")) %>%
  filter(stat %in% c("risk_Blim", "risk_Blim_max", "SSB_rel", #"Fbar_rel",
                     "Catch_rel", "ICV")) %>%
  filter(years %in% c("1:5", "6:10", "11:20", "1:20")) %>%
  mutate(type = "average") %>%
  mutate(type = ifelse(stat == "risk_Blim_max", "max", type)) %>%
  mutate(stat = ifelse(stat == "risk_Blim_max", "risk_Blim", stat)) %>%
  mutate(years = factor(years, levels = c("1:5", "6:10", "11:20", "1:20"),
                        labels = c("'short-term\n(years 1-5)'",
                                   "'medium-term\n(years 6-10)'",
                                   "'long-term\n(years 11-20)'",
                                   "'full period\n(years 1-20)'")),
         stat = factor(stat, levels = c("SSB_rel", "Catch_rel",
                                        "risk_Blim", "risk_Blim_max", "ICV"),
                       labels = c("SSB/italic(B)[MSY]", "Catch/MSY",
                                  "italic(B)[lim]*' risk'", "'max '*B[lim]*' risk'",
                                  "'ICV'"))) %>%
  mutate(value = as.numeric(value),
         multiplier = as.numeric(multiplier))

runs_mult2 %>%
  ggplot(aes(x = multiplier, y = value, linetype = type)) +
  geom_hline(data = data.frame(stat = c("italic(B)[lim]*' risk'"),
                               value = 0.05,
                               type = NA),
             aes(yintercept = value),
             colour = "red", linetype = "dashed", size = 0.3) +
  geom_vline(xintercept = 1.14, color = "blue", linetype = "solid", size = 0.3) +
  geom_line(size = 0.4) +
  facet_grid(stat ~ years, scales = "free", labeller = label_parsed, 
             switch = "y") +
  geom_blank(data = data.frame(multiplier = 1,
                               value = c(1, 1),
                               stat = c("italic(B)[lim]*' risk'", "'ICV'"),
                               type = NA)) +
  scale_linetype("") + 
  ylim(c(0, NA)) +
  labs(x = "multiplier") +
  theme_bw() +
  theme(axis.title.y = element_blank(), 
        strip.placement.y = "outside",
        strip.background.y = element_blank(),
        strip.text.x = element_text(vjust = -1, 
                                    margin = margin(t = 8, b = 4)),
        legend.position = c(0.07, 0.45),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.key.height = unit(0.6, "lines"))
ggsave(filename = "output/ple.27.7e/plots/rfb_mult_stats.png", 
       width = 17, height = 12, units = "cm", dpi = 600, type = "cairo")
# ggsave(filename = "output/ple.27.7e/plots/rfb_mult_stats.pdf", 
#        width = 17, height = 12, units = "cm", dpi = 600)

### plot penalty function
pen <- data.frame(risk = seq(0, 0.1, 0.0001))
pen$penalty <- penalty(x = pen$risk, 
                       negative = FALSE, max = 1, 
                       inflection = 0.05 + 0.01, 
                       steepness = 1000)
pen %>% 
  ggplot(aes(x = risk, y = penalty)) +
  geom_line(size = 0.4) +
  geom_vline(xintercept = 0.05, colour = "red", size = 0.3, 
             linetype = "dashed") +
  labs(x = expression(italic(B)[lim]*' risk')) +
  scale_x_continuous(breaks = seq(0, 0.1, 0.02)) +
  theme_bw()
ggsave(filename = "output/ple.27.7e/plots/penalty_function.png", 
       width = 8.5, height = 5, units = "cm", dpi = 600, type = "cairo")


### explore fitness function formulation
res_fit <- runs_mult %>%
  select(multiplier, Catch_rel, risk_Blim_max) %>%
  mutate(multiplier = unlist(multiplier),
         Catch_rel = unlist(Catch_rel),
         risk_Blim_max = unlist(risk_Blim_max))
head(res_fit)
res_fit$fit1 <- unlist(res_fit$Catch_rel) -
  penalty(x = unlist(res_fit$risk_Blim_max), 
              negative = FALSE, max = 1, 
              inflection = 0.05 + 0.01, 
              steepness = 1000)
#View(res_fit)
res_fit[114:120, ]
res_fit %>%
  pivot_longer(c("Catch_rel", "risk_Blim_max", "fit1")) %>%
  ggplot(aes(x = multiplier, y = value)) +
  geom_line() +
  facet_grid(name ~ 1, scales = "free") +
  theme_bw() +
  geom_vline(xintercept = 1.165)

ggplot() +
  geom_ribbon(data = res_fit,
              aes(ymin = -10, ymax = fit1, x = multiplier),
              fill = "grey30", alpha = 0.5) +
  geom_line(data = res_fit %>%
              pivot_longer(c("Catch_rel", "risk_Blim_max", "fit1")) %>%
              mutate(name = factor(name, 
                                   levels = c("fit1", "Catch_rel", 
                                              "risk_Blim_max"),
                                   labels = c("fitness", "Catch/MSY", 
                                              "'max '*italic(B)[lim]*' risk'"))),
            aes(x = multiplier, y = value, colour = name),
            size = 0.4) +
  scale_colour_manual("", values = c("#377EB8", "#4DAF4A", "#E41A1C"),
                        labels = expression(fitness, Catch/MSY,
                                            'max '*italic(B)[lim]*' risk')) +
  labs(y = "") +
  coord_cartesian(ylim = c(-1.005, 1.005), xlim = c(-0.05, 2.05),
                  expand = FALSE) +
  theme_bw() +
  theme(legend.background = element_blank(),
        legend.key = element_blank(),
        legend.key.height = unit(0.6, "lines"),
        legend.position = c(0.22, 0.3),
        axis.title.y = element_blank(), 
        legend.text.align = 0)
ggsave(filename = "output/ple.27.7e/plots/rfb_mult_stats_fitness.png", 
       width = 8.5, height = 5, units = "cm", dpi = 600, type = "cairo")

### ------------------------------------------------------------------------ ###
### load projections: rfb-rule (PA, mult, all), 2 over 3, SAM ####
### ------------------------------------------------------------------------ ###

### rfb default PA parameterisation: multiplier = 0.95
mp_PA <- readRDS(paste0("output/ple.27.7e/baseline/1000_20/rfb/",
                          "mp_1_2_3_1_1_1_1_2_0.95_1.2_0.7.rds"))
# plot(mp_PA)

### rfb multiplier
ga_mult <- readRDS(paste0("output/ple.27.7e/baseline/1000_20/multiplier/rfb/",
                          "multiplier-upper_constraint1.2-lower_constraint0.7",
                          "--obj_ICES_res_1-20.rds"))
# ga_mult@solution
mp_mult <- readRDS(paste0("output/ple.27.7e/baseline/1000_20/rfb/",
                          "mp_", paste0(ga_mult@solution, collapse = "_"),
                          ".rds"))
# plot(mp_mult)

### rfb all parameters (excluding caps)
ga_all <- readRDS(paste0("output/ple.27.7e/baseline/1000_20/multiplier/rfb/",
                          "range_idx_1-range_idx_2-exp_r-exp_f-exp_b-interval",
                         "-multiplier-upper_constraint1.2-lower_constraint0.7",
                         "--obj_ICES_res_1-20.rds"))
ga_all@solution
par_all <- ga_all@solution
par_all[c(1, 2, 3, 4, 8)] <- round(par_all[c(1, 2, 3, 4, 8)])
par_all[c(5, 6, 7)] <- round(par_all[c(5, 6, 7)], 1)
par_all[c(9)] <- round(par_all[c(9)], 2)
mp_all <- readRDS(paste0("output/ple.27.7e/baseline/1000_20/rfb/",
                          "mp_", paste0(par_all, collapse = "_"), ".rds"))
# plot(mp_mult)

### Category 1 with SAM
mp_SAM <- readRDS("output/ple.27.7e/baseline/1000_20/ICES_SAM/mp.rds")

### 2 over 3 (with XSA)
mp_2over3_XSA <- readRDS("output/ple.27.7e/baseline/1000_20/2over3_XSA/mp.rds")
### 2 over 3 (index only)
mp_2over3 <- readRDS("output/ple.27.7e/baseline/1000_20/2over3/mp.rds")

### load historical data
stk_hist <- readRDS("input/ple.27.7e/baseline/1000_100/stk.rds")
stk_hist <- window(stk_hist, end = 2020)
# plot(stk_hist)

### ------------------------------------------------------------------------ ###
### plot rfb-rule (PA, mult, all), 2 over 3, SAM ####
### ------------------------------------------------------------------------ ###

### extract quantiles for plotting
mp_list <- list(hist = stk_hist, rfb_PA = mp_PA@stock, rfb_mult = mp_mult@stock,
                rfb_all = mp_all@stock, `2over3_XSA` = mp_2over3_XSA@stock,
                `2over3` = mp_2over3@stock, SAM = mp_SAM@stock)
### get quantiles
mp_df <- lapply(seq_along(mp_list), function(x) {#browser()
  qnts <- collapse_correction(mp_list[[x]])
  qnts <- lapply(qnts, quantile, probs = c(0.05, 0.25, 0.5, 0.75, 0.95),
                 na.rm = TRUE)
  qnts <- FLQuants(qnts)
  qnts <- as.data.frame(qnts)
  qnts <- qnts %>%
    mutate(data = ifelse(qname == "catch", data/1000, data),
           data = ifelse(qname == "ssb", data/1000, data),
           data = ifelse(qname == "rec", data/1000, data)) %>%
    #mutate(data = ifelse(qname == "fbar", pmin(data, 1), data)) %>%
    select(year, iter, data, qname) %>%
    pivot_wider(names_from = iter, values_from = data) %>%
    mutate(source = names(mp_list)[[x]])
  return(qnts)
})
mp_df <- do.call(rbind, mp_df)
mp_df <- mp_df %>%
  mutate(source = factor(source, levels = c("hist", "rfb_PA", "rfb_mult", 
                                            "rfb_all", "2over3_XSA", "2over3",
                                            "SAM"),
                         labels = c("hist", "rfb: PA", "rfb: multiplier", 
                                    "rfb: all parameters", "2 over 3 (XSA)",
                                    "2 over 3", "SAM")),
         qname = factor(qname, levels = c("catch", "rec", "fbar", "ssb"),
                        labels = c("Catch [1000t]", "Recruitment [millions]",
                                   "F (ages 3-6)", "SSB [1000 t]")))
### distribution in last simulation year
mp_distr <- lapply(seq_along(mp_list)[-1], function(x) {#browser()
  qnts <- collapse_correction(mp_list[[x]])
  qnts$catch <- qnts$catch/1000
  qnts$ssb <- qnts$ssb/1000
  qnts$rec <- qnts$rec/1000
  qnts <- lapply(seq_along(qnts), function(y) { 
    data.frame(val = c(qnts[[y]][, ac(2040)]), 
               qname = names(qnts)[y],
               source = names(mp_list)[x])
  })
  do.call(rbind, qnts)
})
mp_distr <- do.call(rbind, mp_distr)
mp_distr <- mp_distr %>%
  mutate(source = factor(source, levels = c("rfb_PA", "rfb_mult", 
                                            "rfb_all", "2over3_XSA", "2over3",
                                            "SAM"),
                         labels = c("rfb: PA", "rfb: multiplier", 
                                    "rfb: all parameters", "2 over 3 (XSA)",
                                    "2 over 3", "SAM")),
         qname = factor(qname, levels = c("catch", "rec", "fbar", "ssb"),
                        labels = c("Catch [1000t]", "Recruitment [millions]",
                                   "F (ages 3-6)", "SSB [1000 t]")))

### MSY reference levels
refpts <- readRDS("input/ple.27.7e/baseline/1000_100/refpts_mse.rds")
df_refs <- data.frame(qname = c("Catch [1000t]", "F (ages 3-6)", 
                                "SSB [1000 t]", "Recruitment [millions]"),
                      value = c(c(iterMedians(refpts["Cmsy"]))/1000, 
                                c(iterMedians(refpts["Fmsy"])), 
                                c(iterMedians(refpts["Bmsy"]))/1000, 
                                NA),
                      source = NA)


# debugonce(plot_projection)
plot_projection <- function(mp_df, mp_distr, source, df_refs,
                            ylim = c(0, 4.5), xlim = c(2008, 2040),
                            x_max = 2040, legend.position,
                            plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "pt"),
                            plot_distr = FALSE, distr_limit = 1.05, 
                            distr_yr_limit = 2045) {
  p <- ggplot() +
  geom_ribbon(data = mp_df %>% filter(source == "hist" & 
                                        qname == !!source),
              aes(x = year, ymin = `5%`, ymax = `95%`), alpha = 0.1,
              show.legend = FALSE) +
  geom_ribbon(data = mp_df %>% filter(source == "hist" & 
                                        qname == !!source),
              aes(x = year, ymin = `25%`, ymax = `75%`), alpha = 0.1,
              show.legend = FALSE) +
  geom_ribbon(data = mp_df %>% filter(source != "hist" & 
                                        qname == !!source),
              aes(x = year, ymin = `5%`, ymax = `95%`, fill = source),
              alpha = 0.1, show.legend = FALSE) +
  geom_ribbon(data = mp_df %>% filter(source != "hist" & 
                                        qname == !!source),
              aes(x = year, ymin = `25%`, ymax = `75%`, fill = source),
              alpha = 0.1, show.legend = FALSE) +
  geom_vline(xintercept = 2020, colour = "grey") +
  scale_alpha(guide = guide_legend(label = FALSE)) +
  geom_hline(data = df_refs %>% filter(qname == !!source), 
             aes(yintercept = value, linetype = "MSY"), 
             size = 0.4) +
  geom_line(data = mp_df %>% filter(source == "hist" & 
                                      qname == !!source),
            aes(x = year, y = `50%`), show.legend = FALSE, size = 0.4) +
  geom_line(data = mp_df %>% filter(source != "hist" & 
                                      qname == !!source),
            aes(x = year, y = `50%`, colour = source), show.legend = TRUE,
            size = 0.4) +
  scale_linetype_manual("", values = c("dashed"), 
    guide = guide_legend(override.aes = list(linetype = c("dashed"),
                                             colour = c("black")),
                                             order = 2)) +
  scale_colour_brewer("", palette = "Set1", guide = guide_legend(order = 1)) +
  scale_fill_brewer("", palette = "Set1", 
                    guide = guide_legend(order = 1)) +
  labs(y = source) +
  xlim(c(NA, x_max)) + 
  coord_cartesian(ylim = ylim, xlim = xlim, expand = FALSE) +
  theme_bw(base_size = 8) +
  theme(legend.position = legend.position,
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.title = element_blank(), 
        legend.spacing.y = unit(0, "lines"),
        legend.key = element_blank(),
        legend.key.height = unit(0.6, "lines"),
        legend.key.width = unit(1, "lines"),
        plot.margin = plot.margin)
  if (isTRUE(plot_distr)) {
    p_distr <- mp_distr %>% filter(qname == !!source) %>%
      ggplot(aes(x = val, colour = source)) +
      geom_density(aes(y = ..scaled..), show.legend = FALSE, size = 0.3) +
      scale_colour_brewer("", palette = "Set1") +
      theme_bw(base_size = 8) +
      theme_void() +
      coord_flip(xlim = ylim, ylim = c(0, distr_limit), expand = FALSE) +
      theme(panel.background = element_rect(fill = "white", color = "white"))
    p <- p +
      annotation_custom(grob = ggplotGrob(p_distr),
                        xmin = x_max, xmax = distr_yr_limit,
                        ymin = ylim[1], ymax = ylim[2])
  }
  return(p)
}

p_catch <- plot_projection(mp_df = mp_df, mp_distr = mp_distr,
                df_refs = df_refs, 
                source = "Catch [1000t]", ylim = c(0, 4.5), 
                xlim = c(2008, 2040), x_max = 2040,
                legend.position = c(0.6, 0.7),
                plot.margin = unit(c(5.5, 40, 5.5, 5.5), "pt"),
                plot_distr = TRUE, distr_limit = 1.05, 
                distr_yr_limit = 2045)
p_ssb <- plot_projection(mp_df = mp_df, mp_distr = mp_distr,
                         df_refs = df_refs, 
                         source = "SSB [1000 t]", ylim = c(0, 35), 
                         xlim = c(2008, 2040), x_max = 2040,
                         legend.position = "none",
                         plot.margin = unit(c(5.5, 40, 5.5, 5.5), "pt"),
                         plot_distr = TRUE, distr_limit = 1.05, 
                         distr_yr_limit = 2045)
p_fbar <- plot_projection(mp_df = mp_df, mp_distr = mp_distr,
                          df_refs = df_refs, 
                          source = "F (ages 3-6)", ylim = c(0, 0.8), 
                          xlim = c(2008, 2040), x_max = 2040,
                          legend.position = "none",
                          plot.margin = unit(c(5.5, 40, 5.5, 5.5), "pt"),
                          plot_distr = TRUE, distr_limit = 1.05, 
                          distr_yr_limit = 2045)
p_rec <- plot_projection(mp_df = mp_df, mp_distr = mp_distr,
                         df_refs = df_refs, 
                         source = "Recruitment [millions]", ylim = c(0, 22), 
                         xlim = c(2008, 2040), x_max = 2040,
                         legend.position = "none",
                         plot.margin = unit(c(5.5, 40, 5.5, 5.5), "pt"),
                         plot_distr = TRUE, distr_limit = 1.05, 
                         distr_yr_limit = 2045)

plot_grid(p_catch, p_rec, p_fbar, p_ssb,
          ncol = 2, nrow = 2,
          align = "vh")
ggsave(filename = "output/ple.27.7e/plots/all_MPs_trajectories_combined.png", 
       width = 17, height = 12, units = "cm", dpi = 600, type = "cairo")



### ------------------------------------------------------------------------ ###
### wormplots ####
### ------------------------------------------------------------------------ ###

plot_worm <- function(stock_id = "ple.27.7e", OM, MP, scenario = "", 
                      history = TRUE, its = 1:5,
                      refpts = TRUE, file_name = NULL, 
                      n_years = 20, yr_end = 2040, yr_start = 2000,
                      n_iter = 1000, 
                      xintercept = 2020, ymax_catch = NA, ymax_rec = NA, 
                      ymax_ssb = NA, ymax_fbar = NA, 
                      title_rec = "Recruitment [1000t]",
                      title_catch = "Catch [1000t]",
                      title_ssb = "SSB [1000t]",
                      title_fbar = "Mean F (ages 3-6)"
                      ) {
  #browser()
  ### load projection
  if (is.null(file_name)) {
    file_name <- switch(MP,
                        "rfb: PA" = "mp_1_2_3_1_1_1_1_2_0.95_1.2_0.7.rds",
                        "rfb: mult" = "mp_1_2_3_1_1_1_1_2_1.16_1.2_0.7.rds",
                        "rfb: all" = "mp_1_3_2_1_1_0.3_1.1_2_1.1_1.2_0.7.rds",
                        "2 over 3 (XSA)" = "mp.rds",
                        "2 over 3" = "mp.rds",
                        "SAM" = "mp.rds")
  }
  dir_MP <- switch(MP,
                   "rfb: PA" = "rfb",
                   "rfb: mult" = "rfb",
                   "rfb: all" = "rfb",
                   "2 over 3 (XSA)" = "2over3_XSA",
                   "2 over 3" = "2over3",
                   "SAM" = "ICES_SAM")
  dir_OM <- OM
  if (identical(OM, "rec_failure")) {
    scenario <- OM
    dir_OM <- "baseline"
  }
  file_path <- paste0("output/", stock_id, "/", dir_OM, "/1000_", n_years, 
                      "/", scenario, "/", dir_MP, "/", file_name)
  if (!file.exists(file_path)) return(FALSE)
  res <- readRDS(file_path)
  if (isTRUE(history)) {
    stk <- readRDS(paste0("input/", stock_id, "/", dir_OM, "/1000_100/stk.rds"))
    yrs_res <- dimnames(stock(res))$year
    stk[, ac(yrs_res)] <- stock(res)
  } else {
    stk <- stock(res)
  }
  stk <- window(stk, end = yr_end)
  ### load reference points
  refpts <- readRDS(paste0("input/", stock_id, "/", dir_OM, 
                           "/1000_100/refpts_mse.rds"))
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
  ### iterations
  qnts_iter <- as.data.frame(iter(qnts, its))
  ### plot
  p_catch <- ggplot() +
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
    geom_vline(xintercept = xintercept, colour = "grey") +
    geom_hline(yintercept = c(median(refpts["Cmsy"], na.rm = TRUE))/1000,
               colour = "black", size = 0.5, linetype = "dashed") +
    coord_cartesian(xlim = c(yr_start, yr_end), ylim = c(0, ymax_catch)) + 
    labs(y = title_catch) +
    theme_bw(base_size = 8) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank())
  p_rec <- ggplot() +
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
    geom_vline(xintercept = xintercept, colour = "grey") +
    coord_cartesian(xlim = c(yr_start, yr_end), ylim = c(0, ymax_rec)) + 
    labs(y = title_rec) +
    theme_bw(base_size = 8) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank())
  p_ssb <- ggplot() +
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
    geom_vline(xintercept = xintercept, colour = "grey") +
    geom_hline(yintercept = c(median(refpts["Bmsy"], na.rm = TRUE))/1000,
               colour = "black", size = 0.5, linetype = "dashed") +
    geom_hline(yintercept = c(median(refpts["Blim"], na.rm = TRUE))/1000,
               colour = "black", size = 0.5, linetype = "dotted") +
    coord_cartesian(xlim = c(yr_start, yr_end), ylim = c(0, ymax_ssb)) + 
    labs(y = title_ssb) +
    theme_bw(base_size = 8)
  p_fbar <- ggplot() +
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
    geom_vline(xintercept = xintercept, colour = "grey") +
    geom_hline(yintercept = c(median(refpts["Fmsy"], na.rm = TRUE)),
               colour = "black", size = 0.5, linetype = "dashed") +
    coord_cartesian(xlim = c(yr_start, yr_end), ylim = c(0, ymax_fbar)) + 
    labs(y = title_fbar) +
    theme_bw(base_size = 8)
  p <-  plot_grid(p_catch, p_rec, p_fbar, p_ssb, align = "v", 
                  rel_heights = c(1, 1.1))
  return(p)
}

### baseline OM - rfb: PA
plot_worm(OM = "baseline", MP = "rfb: PA")
ggsave(filename = "output/ple.27.7e/plots/wormplots/baseline_rfb_PA.png", 
       width = 17, height = 12, units = "cm", dpi = 600, type = "cairo")
plot_worm(OM = "baseline", MP = "rfb: mult")

### plot all MPs and all OMs
MPs <- c("rfb: PA", "rfb: mult", "rfb: all", "2 over 3", "2 over 3 (XSA)",
         "SAM")
MPs_dir <- c("rfb_PA", "rfb_mult", "rfb_all", "2over3", "2over3_XSA", "SAM")
OMs <- c("baseline", "M_low", "M_high", "M_Gislason", "rec_no_AC", 
         "no_discards", "rec_failure")
. <- foreach(MP = MPs, MP_dir = MPs_dir) %:% foreach(OM = OMs) %do% {
  p <- plot_worm(OM = OM, MP = MP, ymax_fbar = 0.8)
  if (isFALSE(p)) return(NULL)
  ggsave(filename = paste0("output/ple.27.7e/plots/wormplots/",
                           OM, "_", MP_dir, ".png"), 
         plot = p,
         width = 17, height = 12, units = "cm", dpi = 600, type = "cairo")
}



### ------------------------------------------------------------------------ ###
### alternative OMs - stats ####
### ------------------------------------------------------------------------ ###

OM_stats <- foreach(OM = c("baseline", "M_low", "M_high", "M_Gislason",
                           "rec_no_AC", "no_discards", "rec_failure"),
                    stock_id = rep("ple.27.7e", 7),
                    .combine = bind_rows) %:%
  foreach(MP = c("rfb: PA", "rfb: mult", "rfb: all", "2 over 3 (XSA)", 
                 "2 over 3", "SAM"), 
          interval = c(2, 2, 2, 1, 2, 1),
          .combine = bind_rows) %do% {
    #browser()
    
    ### load projection
    file_name <- switch(MP,
                        "rfb: PA" = "mp_1_2_3_1_1_1_1_2_0.95_1.2_0.7.rds",
                        "rfb: mult" = "mp_1_2_3_1_1_1_1_2_1.16_1.2_0.7.rds",
                        "rfb: all" = "mp_1_3_2_1_1_0.3_1.1_2_1.1_1.2_0.7.rds",
                        "2 over 3 (XSA)" = "mp.rds",
                        "2 over 3" = "mp.rds",
                        "SAM" = "mp.rds")
    dir_MP <- switch(MP,
                     "rfb: PA" = "rfb",
                     "rfb: mult" = "rfb",
                     "rfb: all" = "rfb",
                     "2 over 3 (XSA)" = "2over3_XSA",
                     "2 over 3" = "2over3",
                     "SAM" = "ICES_SAM")
    dir_OM <- OM
    # n_years <- ifelse(MP == "SAM" & OM == "baseline", 100, 20)
    scenario <- ""
    if (identical(OM, "rec_failure")) {
      scenario <- OM
      dir_OM <- "baseline"
    }
    
    file_path <- paste0("output/", stock_id, "/", dir_OM, "/1000_", 20, 
                        "/", scenario, "/", dir_MP, "/", file_name)
    if (isFALSE(file.exists(file_path))) return(NULL)
    res <- readRDS(file_path)
    ### load reference points
    refpts <- readRDS(paste0("input/", stock_id, "/", dir_OM, 
                             "/1000_", 100, "/refpts_mse.rds"))
    
    ### extract metrics
    stk <- window(res@stock, start = 2031, end = 2040)
    stk_icv <- window(res@stock, start = 2030, end = 2040)
    ssb_i <- c(ssb(stk)/refpts["Bmsy"])
    ssb20_i <- c(ssb(stk)[, ac(2040)]/refpts["Bmsy"])
    catch_i <- c(catch(stk)/refpts["Cmsy"])
    fbar_i <- c(fbar(stk)/refpts["Fmsy"])
    risk_i <- c(apply(ssb(stk) < c(refpts["Blim"]), 2, mean))
    icv_i <- c(iav(catch(stk_icv), period = interval))
    if (identical(OM, "no_discards")) catch_i <- c(landings(stk)/refpts["Cmsy"])
    ### combine
    df <- do.call(rbind, list(data.frame(val = ssb_i, metric = "SSB"),
                              data.frame(val = ssb20_i, metric = "SSB20"),
                              data.frame(val = catch_i, metric = "catch"),
                              data.frame(val = fbar_i, metric = "Fbar"),
                              data.frame(val = icv_i, metric = "ICV"),
                              data.frame(val = risk_i, metric = "risk")
    ))
    df$OM <- OM
    df$MP <- MP
    df$stock_id <- stock_id
    return(df)
    
}
### combine all
stats_combined <- OM_stats %>%
  mutate(OM_group = 
           case_when(OM %in% c("baseline") ~ "baseline",
                     OM %in% c("M_low", "M_high", "M_Gislason") ~ "M",
                     OM %in% c("rec_no_AC", "rec_failure") ~ "Rec",
                     OM %in% c("no_discards") ~ "Catch"),
         OM_group = factor(OM_group, levels = c("baseline", "M", "Rec", 
                                                "Catch")),
         OM_label = factor(OM, levels = c("baseline", "M_low", "M_high", 
                                          "M_Gislason", "rec_no_AC", 
                                          "rec_failure", "no_discards"),
                           labels = c("baseline", "low", "high", "Gislason",
                                      "no AC", "failure", "no discards")),
         MP = factor(MP, levels = c("rfb: PA", "rfb: mult", "rfb: all", 
                                    "2 over 3 (XSA)", "2 over 3", "SAM")))

### baseline OM - all MPs
p_catch <- stats_combined %>%
  filter(OM == "baseline" & metric == "catch") %>%
  ggplot(aes(x = MP, y = val)) +
  geom_hline(yintercept = 1, colour = "grey") +
  geom_violin(aes(fill = MP), size = 0.2, show.legend = FALSE,
              position = position_dodge(width = 0.5)) +
  scale_fill_brewer(palette = "Set1") +
  geom_boxplot(aes(group = interaction(OM, MP)), 
               position = position_dodge(width = 0.5),
               fill = "white", width = 0.1, size = 0.2,
               outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
               outlier.fill = "transparent") +
  labs(y = expression(catch/MSY)) +
  coord_cartesian(ylim = c(0, 2.5)) +
  theme_bw(base_size = 8) +
  theme(panel.spacing.x = unit(0, "lines"),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
p_fbar <- stats_combined %>%
  filter(OM == "baseline" &  metric == "Fbar") %>%
  ggplot(aes(x = MP, y = val)) +
  geom_hline(yintercept = 1, colour = "grey") +
  geom_violin(aes(fill = MP), size = 0.2, show.legend = FALSE,
              position = position_dodge(width = 0.5)) +
  scale_fill_brewer(palette = "Set1") +
  geom_boxplot(aes(group = interaction(OM, MP)), 
               position = position_dodge(width = 0.5),
               fill = "white", width = 0.1, size = 0.2,
               outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
               outlier.fill = "transparent") +
  labs(y = expression(F/F[MSY])) +
  coord_cartesian(ylim = c(0, 3.5)) +
  theme_bw(base_size = 8) +
  theme(panel.spacing.x = unit(0, "lines"),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
p_SSB <- stats_combined %>%
  filter(OM == "baseline" &  metric == "SSB") %>%
  ggplot(aes(x = MP, y = val)) +
  geom_hline(yintercept = 1, colour = "grey") +
  geom_violin(aes(fill = MP), size = 0.2, show.legend = FALSE,
              position = position_dodge(width = 0.5)) +
  scale_fill_brewer(palette = "Set1") +
  geom_boxplot(aes(group = interaction(OM, MP)), 
               position = position_dodge(width = 0.5),
               fill = "white", width = 0.1, size = 0.2,
               outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
               outlier.fill = "transparent") +
  labs(y = expression(SSB/B[MSY])) +
  coord_cartesian(ylim = c(0, 3.5)) +
  theme_bw(base_size = 8) +
  theme(panel.spacing.x = unit(0, "lines"),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
p_SSB20 <- stats_combined %>%
  filter(OM == "baseline" & metric == "SSB20") %>%
  ggplot(aes(x = MP, y = val)) +
  geom_hline(yintercept = 1, colour = "grey") +
  geom_violin(aes(fill = MP), size = 0.2, show.legend = FALSE,
              position = position_dodge(width = 0.5)) +
  scale_fill_brewer(palette = "Set1") +
  geom_boxplot(aes(group = interaction(OM, MP)), 
               position = position_dodge(width = 0.5),
               fill = "white", width = 0.1, size = 0.2,
               outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
               outlier.fill = "transparent") +
  labs(y = expression(SSB20/B[MSY])) +
  theme_bw(base_size = 8) +
  theme(panel.spacing.x = unit(0, "lines"),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
p_ICV <- stats_combined %>%
  filter(OM == "baseline" & metric == "ICV") %>%
  ggplot(aes(x = MP, y = val)) +
  geom_violin(aes(fill = MP), size = 0.2, show.legend = FALSE,
              position = position_dodge(width = 0.5)) +
  scale_fill_brewer(palette = "Set1") +
  geom_boxplot(aes(group = interaction(OM, MP)), 
               position = position_dodge(width = 0.5),
               fill = "white", width = 0.1, size = 0.2,
               outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
               outlier.fill = "transparent") +
  labs(y = expression(ICV)) +
  coord_cartesian(ylim = c(0, 0.5)) +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.spacing.x = unit(0, "lines"),
        axis.title.x = element_blank())
p_risk <- stats_combined %>%
  filter(OM == "baseline" & metric == "risk") %>%
  ggplot() +
  geom_hline(yintercept = 0.05, colour = "red") +
  geom_col(data = stats_combined %>%
             filter(OM == "baseline" & metric == "risk") %>%
             group_by(MP) %>%
             summarise(val = max(val)),
           aes(x = MP, y = val, fill = MP), 
           show.legend = FALSE, width = 0.5, colour = "black", size = 0.2,
           position = position_dodge(width = 0.5)) +
  scale_fill_brewer("", palette = "Set1") +
  geom_boxplot(aes(x = MP, y = val, group = interaction(OM, MP)),
               position = position_dodge(width = 0.5),
               fill = "white", width = 0.1, size = 0.2,
               outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
               outlier.fill = "transparent") +
  labs(y = expression(risk)) +
  ylim(c(0, NA)) +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.spacing.x = unit(0, "lines"),
        axis.title.x = element_blank(),
        legend.position = c(0.1, 0.7),
        legend.key = element_blank(),
        legend.key.height = unit(0.4, "lines"),
        legend.key.width = unit(0.3, "lines"),
        legend.background = element_blank())

plot_grid(p_catch, p_fbar, p_SSB, p_SSB20,
          p_ICV, p_risk, ncol = 2, align = "v", 
          rel_heights = c(1, 1, 1.3))
ggsave(filename = "output/ple.27.7e/plots/baseline_MPs_stats_11-20_log.png", 
       width = 17, height = 9, units = "cm", dpi = 600, type = "cairo")


### alternative OMs - rfb-rule (PA & mult) & SAM
p_catch <- stats_combined %>%
  filter(MP %in% c("rfb: PA", "rfb: mult", "SAM") & 
           metric == "catch") %>%
  ggplot(aes(x = OM_label, y = val)) +
  geom_hline(yintercept = 1, colour = "grey") +
  geom_violin(aes(fill = MP), size = 0.2, show.legend = FALSE,
              position = position_dodge(width = 0.5)) +
  geom_boxplot(aes(group = interaction(OM, MP)), 
               position = position_dodge(width = 0.5),
               fill = "white", width = 0.1, size = 0.2,
               outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
               outlier.fill = "transparent") +
  facet_grid(~ OM_group, scales = "free_x", space = "free_x") +
  labs(y = expression(catch/MSY)) +
  coord_cartesian(ylim = c(0, 2.5)) +
  theme_bw(base_size = 8) +
  theme(panel.spacing.x = unit(0, "lines"),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
p_fbar <- stats_combined %>%
  filter(MP %in% c("rfb: PA", "rfb: mult", "SAM") & metric == "Fbar") %>%
  ggplot(aes(x = OM_label, y = val)) +
  geom_hline(yintercept = 1, colour = "grey") +
  geom_violin(aes(fill = MP), size = 0.2, show.legend = FALSE,
              position = position_dodge(width = 0.5)) +
  geom_boxplot(aes(group = interaction(OM, MP)), 
               position = position_dodge(width = 0.5),
               fill = "white", width = 0.1, size = 0.2,
               outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
               outlier.fill = "transparent") +
  facet_grid(~ OM_group, scales = "free_x", space = "free_x") +
  labs(y = expression(F/F[MSY])) +
  coord_cartesian(ylim = c(0, 3.5)) +
  theme_bw(base_size = 8) +
  theme(panel.spacing.x = unit(0, "lines"),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
p_SSB <- stats_combined %>%
  filter(MP %in% c("rfb: PA", "rfb: mult", "SAM") & metric == "SSB") %>%
  ggplot(aes(x = OM_label, y = val)) +
  geom_hline(yintercept = 1, colour = "grey") +
  geom_violin(aes(fill = MP), size = 0.2, show.legend = FALSE,
              position = position_dodge(width = 0.5)) +
  geom_boxplot(aes(group = interaction(OM, MP)), 
               position = position_dodge(width = 0.5),
               fill = "white", width = 0.1, size = 0.2,
               outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
               outlier.fill = "transparent") +
  facet_grid(~ OM_group, scales = "free_x", space = "free_x") +
  labs(y = expression(SSB/B[MSY])) +
  coord_cartesian(ylim = c(0, 3.5)) +
  theme_bw(base_size = 8) +
  theme(panel.spacing.x = unit(0, "lines"),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
p_SSB20 <- stats_combined %>%
  filter(MP %in% c("rfb: PA", "rfb: mult", "SAM") & metric == "SSB20") %>%
  ggplot(aes(x = OM_label, y = val)) +
  geom_hline(yintercept = 1, colour = "grey") +
  geom_violin(aes(fill = MP), size = 0.2, show.legend = FALSE,
              position = position_dodge(width = 0.5)) +
  geom_boxplot(aes(group = interaction(OM, MP)), 
               position = position_dodge(width = 0.5),
               fill = "white", width = 0.1, size = 0.2,
               outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
               outlier.fill = "transparent") +
  facet_grid(~ OM_group, scales = "free_x", space = "free_x") +
  labs(y = expression(SSB20/B[MSY])) +
  coord_cartesian(ylim = c(0, 3.5)) +
  theme_bw(base_size = 8) +
  theme(panel.spacing.x = unit(0, "lines"),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
p_ICV <- stats_combined %>%
  filter(MP %in% c("rfb: PA", "rfb: mult", "SAM") & metric == "ICV") %>%
  ggplot(aes(x = OM_label, y = val)) +
  geom_violin(aes(fill = MP), size = 0.2, show.legend = FALSE,
              position = position_dodge(width = 0.5)) +
  geom_boxplot(aes(group = interaction(OM, MP)), 
               position = position_dodge(width = 0.5),
               fill = "white", width = 0.1, size = 0.2,
               outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
               outlier.fill = "transparent") +
  facet_grid(~ OM_group, scales = "free_x", space = "free_x") +
  labs(y = expression(ICV)) +
  coord_cartesian(ylim = c(0, 0.5)) +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.spacing.x = unit(0, "lines"),
        axis.title.x = element_blank())
p_risk <- stats_combined %>%
  filter(MP %in% c("rfb: PA", "rfb: mult", "SAM") & metric == "risk") %>%
  ggplot() +
  geom_hline(yintercept = 0.05, colour = "red") +
  geom_col(data = stats_combined %>%
             filter(MP %in% c("rfb: PA", "rfb: mult", "SAM") & 
                      metric == "risk") %>%
             group_by(MP, OM, OM_label, OM_group) %>%
             summarise(val = max(val)),
           aes(x = OM_label, y = val, fill = MP), 
           show.legend = TRUE, width = 0.5, colour = "black", size = 0.2,
           position = position_dodge(width = 0.5)) +
  geom_boxplot(aes(x = OM_label, y = val, group = interaction(OM, MP)),
               position = position_dodge(width = 0.5),
               fill = "white", width = 0.1, size = 0.2,
               outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
               outlier.fill = "transparent") +
  scale_fill_discrete("") +
  facet_grid(~ OM_group, scales = "free_x", space = "free_x") +
  labs(y = expression(risk)) +
  ylim(c(0, NA)) +
  # scale_y_continuous(trans = "sqrt", 
  #                    breaks = c(0, 0.001, 0.01, 0.025, 0.05, 0.1, 0.15, 0.25), 
  #                    limits = c(0, NA)) +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.spacing.x = unit(0, "lines"),
        axis.title.x = element_blank(),
        legend.position = c(0.25, 0.85),
        legend.key = element_blank(),
        legend.key.height = unit(0.4, "lines"),
        legend.key.width = unit(0.3, "lines"),
        legend.background = element_blank())

plot_grid(p_catch, p_fbar, p_SSB, p_SSB20,
          p_ICV, p_risk, ncol = 2, align = "v", 
          rel_heights = c(1, 1, 1.3))
ggsave(filename = "output/ple.27.7e/plots/robustness_altOMs_stats_11-20_log.png", 
       width = 17, height = 9, units = "cm", dpi = 600, type = "cairo")





