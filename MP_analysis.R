### ------------------------------------------------------------------------ ###
### ple.27.7e - analyse MSE results ####
### ------------------------------------------------------------------------ ###
library(mse)
library(GA)
library(tidyr)
library(dplyr)
library(cowplot)
library(patchwork)
library(ggplot2)
library(foreach)
source("funs.R")
source("funs_GA.R")
source("funs_analysis.R")

### ------------------------------------------------------------------------ ###
### collate results - baseline MPs ####
### ------------------------------------------------------------------------ ###

### baseline OM runs for all MPs (including optimisations)
res <- foreach(MP = c("rfb", "2over3", "2over3_XSA", "hr", "ICES_SAM"), 
               .combine = bind_rows) %:%
  foreach(stock = c("ple.27.7e", "cod.27.47d20", "her.27.3a47d"), 
          .combine = bind_rows) %:%
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
        if (stock %in% c("ple.27.7e", "cod.27.47d20")) {
          ga_solution$multiplier <- 0.95
        } else if (stock %in% c("her.27.3a47d")) {
          ga_solution$multiplier <- 0.9
        }
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
      ### catch interval
      res_i$interval <- case_when(
        MP %in% c("2over3") ~ 2,
        MP %in% c("2over3_XSA", "ICES_SAM") ~ 1
      )
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

write.csv(res, file = "output/MPs_baseline.csv", row.names = FALSE)
saveRDS(res, file = "output/MPs_baseline.rds")
res <- readRDS("output/MPs_baseline.rds")

### ------------------------------------------------------------------------ ###
### summary table for baseline OM: optimisation results ####
### ------------------------------------------------------------------------ ###
res <- readRDS("output/MPs_baseline.rds")
res %>% 
  filter(period == "11-20" & MP %in% c("rfb", "hr")) %>%
  select(`stock`:lower_constraint, idxB_lag, idxB_range_3, comp_b_multiplier, 
         -OM, -period) %>%
  group_by(stock, MP) %>%
  ### relative growth -> used because some fitness values might be negative
  mutate(fitness_improvement = ((fitness - min(fitness))/abs(min(fitness)))*100,
         .after = fitness) %>%
  mutate(fitness_improvement = round(fitness_improvement)) %>%
  write.csv("tmp.csv", row.names = FALSE)


### ------------------------------------------------------------------------ ###
### collate results - alternative OMs ####
### ------------------------------------------------------------------------ ###

### include alternative OMs
OMs_ple <- c("baseline", "M_low", "M_high", "M_Gislason", 
             "no_discards", "rec_no_AC", "rec_failure")
OMs_cod <- c("baseline", "rec_higher", "M_dd", "M_no_migration", "rec_failure")
OMs_her <- c("baseline", "rec_higher", "rec_failure")
res_alt <- res %>%
  filter(period == "11-20") %>%
  select(stock:period, fitness, lag_idx:lower_constraint, file) %>%
  mutate(id = seq(n()))
res_alt <- foreach(OM = unique(c(OMs_ple, OMs_cod, OMs_her)),
                   .combine = bind_rows) %:% 
  foreach(i = split(res_alt, f = seq(nrow(res_alt))), 
          .combine = bind_rows) %do% {
    #browser()
    i$OM <- OM
    if (identical(OM, "baseline")) return(i)
    if ((identical(i$stock, "ple.27.7e") & !OM %in% OMs_ple) |
        (identical(i$stock, "cod.27.47d20") & !OM %in% OMs_cod) |
        (identical(i$stock, "her.27.3a47d") & !OM %in% OMs_her)) 
      return(NULL)
    file_i <- i$file
    if (identical(OM, "rec_failure")) {
      file_i <- gsub(x = file_i, pattern = "baseline/1000_20", 
                     replacement = "baseline/1000_20/rec_failure")
    } else {
      file_i <- gsub(x = file_i, pattern = "baseline", replacement = OM)
    }
    i$file <- ifelse(file.exists(file_i), file_i, NA)
    i$fitness <- NA
    i$OM_group = case_when(
      OM == "baseline" ~ "baseline",
      OM %in% c("M_low", "M_high", "M_Gislason", "M_dd", 
                        "M_no_migration") ~ "M",
      OM %in% c("rec_higher", "rec_no_AC", "rec_failure") ~ "Rec",
      OM %in% c("no_discards") ~ "Catch")
    return(i)
}
res_alt$OM_group[res_alt$OM == "baseline"] <- "baseline"
#sum(table(res_alt$OM_group))
#nrow(res_alt)
#table(res_alt$id)
#any(is.na(res_alt$file))
#View(res_alt)

saveRDS(res_alt, file = "output/MPs_alternative_OMs.rds")


### ------------------------------------------------------------------------ ###
### collate results - alternative OMs - stats ####
### ------------------------------------------------------------------------ ###
res_alt <- readRDS("output/MPs_alternative_OMs.rds")
stats_alt <- foreach(i = split(res_alt, f = seq(nrow(res_alt))), 
        .combine = bind_rows) %do% {#browser()
  ### get projections and reference points
  mp_i <- readRDS(i$file)
  path_OM <- paste0("input/", i$stock, "/", 
                    ifelse(identical(i$OM, "rec_failure"), "baseline", i$OM), 
                    "/1000_100/")
  refpts_i <- iterMedians(readRDS(paste0(path_OM, "refpts_mse.rds")))
  ### find OM stock - slot depends on version of mse package
  . <- try(stk <- mp_i@om@stock, silent = TRUE)
  if (is(., "try-error")) stk <- mp_i@stock
  ### extract metrics
  stk_icv <- window(stk, start = 2030, end = 2040)
  stk <- window(stk, start = 2031, end = 2040)
  ssb_i <- c(ssb(stk)/refpts_i["Bmsy"])
  ssb20_i <- c(ssb(stk)[, ac(2040)]/refpts_i["Bmsy"])
  catch_i <- c(catch(stk)/refpts_i["Cmsy"])
  fbar_i <- c(fbar(stk)/refpts_i["Fmsy"])
  risk_i <- c(apply(ssb(stk) < c(refpts_i["Blim"]), 2, mean))
  icv_i <- c(iav(catch(stk_icv), period = i$interval))
  if (identical(i$OM, "no_discards")) 
    catch_i <- c(landings(stk)/refpts_i["Cmsy"])
  ### combine
  df <- do.call(rbind, list(data.frame(val = ssb_i, metric = "SSB"),
                            data.frame(val = ssb20_i, metric = "SSB20"),
                            data.frame(val = catch_i, metric = "catch"),
                            data.frame(val = fbar_i, metric = "Fbar"),
                            data.frame(val = icv_i, metric = "ICV"),
                            data.frame(val = risk_i, metric = "risk")
  ))
  return(bind_cols(i, df))
}
stats_alt <- stats_alt %>%
  mutate(
    OM_label = factor(OM, 
      levels = c("baseline", "M_low", "M_high", "M_Gislason",
                 "M_dd", "M_no_migration", "no_discards",
                 "rec_higher", "rec_no_AC", "rec_failure"),
      labels = c("baseline", "low", "high", "Gislason",
                 "dens. dep.", "no migration", "no discards",
                 "higher", "no AC", "failure")), .after = "OM") %>%
  mutate(OM_group = factor(OM_group, c("baseline", "M", "Catch", "Rec"))) %>%
  mutate(
    MP_label = factor(paste0(MP, "_", optimised), 
      levels = c("2over3_default", "2over3_XSA_default", 
                 "rfb_default", "rfb_multiplier", "rfb_all", 
                 "hr_default", "hr_multiplier", "hr_all", 
                 "ICES_SAM_default"),
      labels = c("2 over 3", "2 over 3 (XSA)", 
                 "rfb (generic)", "rfb (multiplier)", "rfb (all)",
                 "hr (generic)", "hr (multiplier)", "hr (all)",
                 "ICES MSY (SAM)")), .after = "MP") %>%
  mutate(stock_label = factor(stock, 
                              levels = c("ple.27.7e", "cod.27.47d20",
                                         "her.27.3a47d"),
                              labels = c("Plaice", "Cod", "Herring")))

saveRDS(stats_alt, file = "output/MPs_alternative_OMs_stats.rds")
# stats_alt <- readRDS("output/MPs_alternative_OMs_stats.rds")

### ------------------------------------------------------------------------ ###
### wormplots for all MPs/OMs ####
### ------------------------------------------------------------------------ ###
res_alt <- readRDS("output/MPs_alternative_OMs.rds")

### go through all runs
for (i in split(res_alt, f = seq(nrow(res_alt)))) {#browser()
  ### find OM stock - slot depends on version of mse package
  res_mp <- readRDS(i$file)
  . <- try(stk_res <- res_mp@om@stock, silent = TRUE)
  if (is(., "try-error")) stk_res <- res_mp@stock
  path_OM <- paste0("input/", i$stock, "/", 
                    ifelse(identical(i$OM, "rec_failure"), "baseline", i$OM), 
                    "/1000_100/")
  stk_hist <- readRDS(paste0(path_OM, "stk.rds"))
  refpts <- readRDS(paste0(path_OM, "refpts_mse.rds"))
  p <- plot_worm(stk = stk_res, stk_hist = stk_hist, refpts = refpts)
  ggsave(filename = paste0("output/plots/wormplots/", i$stock, "_", i$OM, "_",
                           i$MP, "_", i$optimised, ".png"), plot = p, 
         width = 17, height = 8, units = "cm", dpi = 600, type = "cairo")
  ggsave(filename = paste0("output/plots/wormplots/", i$stock, "_", i$OM, "_",
                           i$MP, "_", i$optimised, ".pdf"), plot = p, 
         width = 17, height = 8, units = "cm")
}

### ------------------------------------------------------------------------ ###
### violin plots - all stocks/OMs/MPs  ####
### ------------------------------------------------------------------------ ###
stats_alt <- readRDS("output/MPs_alternative_OMs_stats.rds")

col_vals <- c("2 over 3" = "#9e9ac8", 
              "2 over 3 (XSA)" = "#6a51a3", 
              "rfb (generic)" = "#bdd7e7", 
              "rfb (multiplier)" = "#6baed6", 
              "rfb (all)" = "#2171b5", 
              "hr (generic)" = "#fcae91", 
              "hr (multiplier)" = "#fb6a4a", 
              "hr (all)" = "#cb181d", 
              "ICES MSY" = "#ffff00")
p_ple_catch <- stats_alt %>%
  filter(stock == "ple.27.7e" &
           metric == "catch") %>%
  ggplot(aes(x = OM_label, y = val)) +
  geom_hline(yintercept = 1, colour = "grey") +
  geom_violin(aes(fill = MP_label), size = 0.2, show.legend = FALSE,
              position = position_dodge(width = 0.8), scale = "width") +
  geom_boxplot(aes(group = interaction(OM, MP_label)), 
               position = position_dodge(width = 0.8),
               fill = "white", width = 0.1, size = 0.2,
               outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
               outlier.fill = "transparent") +
  facet_grid(~ OM_group, scales = "free_x", space = "free_x") +
  scale_fill_manual("", values = col_vals) +
  labs(y = expression(Catch/MSY)) +
  coord_cartesian(ylim = c(0, 2.5)) +
  theme_bw(base_size = 8) +
  theme(panel.spacing.x = unit(0, "lines"),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
p_ple_SSB <- stats_alt %>%
  filter(stock == "ple.27.7e" &
           metric == "SSB") %>%
  ggplot(aes(x = OM_label, y = val)) +
  geom_hline(yintercept = 1, colour = "grey") +
  geom_violin(aes(fill = MP_label), size = 0.2, show.legend = FALSE,
              position = position_dodge(width = 0.8), scale = "width") +
  geom_boxplot(aes(group = interaction(OM, MP_label)), 
               position = position_dodge(width = 0.8),
               fill = "white", width = 0.1, size = 0.2,
               outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
               outlier.fill = "transparent") +
  facet_grid(~ OM_group, scales = "free_x", space = "free_x") +
  scale_fill_manual("", values = col_vals) +
  labs(y = expression(SSB/B[MSY])) +
  coord_cartesian(ylim = c(0, 3.5)) +
  theme_bw(base_size = 8) +
  theme(panel.spacing.x = unit(0, "lines"),
        strip.text.x = element_blank(),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
p_ple_risk <- stats_alt %>%
  filter(stock == "ple.27.7e" & metric == "risk") %>%
  ggplot() +
  geom_hline(yintercept = 0.055, colour = "red") +
  geom_col(data = stats_alt %>%
             filter(stock == "ple.27.7e" & metric == "risk") %>%
             group_by(MP_label, OM, OM_label, OM_group) %>%
             summarise(val = max(val)),
           aes(x = OM_label, y = val, fill = MP_label), 
           show.legend = TRUE, width = 0.8, colour = "black", size = 0.2,
           position = position_dodge(width = 0.8)) +
  geom_boxplot(aes(x = OM_label, y = val, group = interaction(OM, MP_label)),
               position = position_dodge(width = 0.8),
               fill = "white", width = 0.1, size = 0.2,
               outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
               outlier.fill = "transparent") +
  scale_fill_manual("", values = col_vals) +
  facet_grid(~ OM_group, scales = "free_x", space = "free_x") +
  labs(y = expression(max.~B[lim]~risk)) +
  ylim(c(0, NA)) +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.spacing.x = unit(0, "lines"),
        strip.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = c(0.25, 0.75),
        legend.key = element_blank(),
        legend.key.height = unit(0.4, "lines"),
        legend.key.width = unit(0.3, "lines"),
        legend.background = element_blank())
p_cod_catch <- stats_alt %>%
  filter(stock == "cod.27.47d20" &
           metric == "catch") %>%
  ggplot(aes(x = OM_label, y = val)) +
  geom_hline(yintercept = 1, colour = "grey") +
  geom_violin(aes(fill = MP_label), size = 0.2, show.legend = FALSE,
              position = position_dodge(width = 0.8), scale = "width") +
  geom_boxplot(aes(group = interaction(OM, MP_label)), 
               position = position_dodge(width = 0.8),
               fill = "white", width = 0.1, size = 0.2,
               outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
               outlier.fill = "transparent") +
  facet_grid(~ OM_group, scales = "free_x", space = "free_x") +
  scale_fill_manual("", values = col_vals) +
  labs(y = expression(Catch/MSY)) +
  coord_cartesian(ylim = c(0, 2.5)) +
  theme_bw(base_size = 8) +
  theme(panel.spacing.x = unit(0, "lines"),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
p_cod_SSB <- stats_alt %>%
  filter(stock == "cod.27.47d20" &
           metric == "SSB") %>%
  ggplot(aes(x = OM_label, y = val)) +
  geom_hline(yintercept = 1, colour = "grey") +
  geom_violin(aes(fill = MP_label), size = 0.2, show.legend = FALSE,
              position = position_dodge(width = 0.8), scale = "width") +
  geom_boxplot(aes(group = interaction(OM, MP_label)), 
               position = position_dodge(width = 0.8),
               fill = "white", width = 0.1, size = 0.2,
               outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
               outlier.fill = "transparent") +
  facet_grid(~ OM_group, scales = "free_x", space = "free_x") +
  scale_fill_manual("", values = col_vals) +
  labs(y = expression(SSB/B[MSY])) +
  coord_cartesian(ylim = c(0, 6)) +
  theme_bw(base_size = 8) +
  theme(panel.spacing.x = unit(0, "lines"),
        strip.text.x = element_blank(),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
p_cod_risk <- stats_alt %>%
  filter(stock == "cod.27.47d20" & metric == "risk") %>%
  ggplot() +
  geom_hline(yintercept = 0.055, colour = "red") +
  geom_col(data = stats_alt %>%
             filter(stock == "cod.27.47d20" & metric == "risk") %>%
             group_by(MP_label, OM, OM_label, OM_group) %>%
             summarise(val = max(val)),
           aes(x = OM_label, y = val, fill = MP_label), 
           show.legend = TRUE, width = 0.8, colour = "black", size = 0.2,
           position = position_dodge(width = 0.8)) +
  geom_boxplot(aes(x = OM_label, y = val, group = interaction(OM, MP_label)),
               position = position_dodge(width = 0.8),
               fill = "white", width = 0.1, size = 0.2,
               outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
               outlier.fill = "transparent") +
  facet_grid(~ OM_group, scales = "free_x", space = "free_x") +
  scale_fill_manual("", values = col_vals) +
  labs(y = expression(max.~B[lim]~risk)) +
  ylim(c(0, NA)) +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.spacing.x = unit(0, "lines"),
        strip.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.key = element_blank(),
        legend.key.height = unit(0.4, "lines"),
        legend.key.width = unit(0.3, "lines"),
        legend.background = element_blank())
p_her_catch <- stats_alt %>%
  filter(stock == "her.27.3a47d" &
           metric == "catch") %>%
  ggplot(aes(x = OM_label, y = val)) +
  geom_hline(yintercept = 1, colour = "grey") +
  geom_violin(aes(fill = MP_label), size = 0.2, show.legend = FALSE,
              position = position_dodge(width = 0.8), scale = "width") +
  geom_boxplot(aes(group = interaction(OM, MP_label)), 
               position = position_dodge(width = 0.8),
               fill = "white", width = 0.1, size = 0.2,
               outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
               outlier.fill = "transparent") +
  facet_grid(~ OM_group, scales = "free_x", space = "free_x") +
  scale_fill_manual("", values = col_vals) +
  labs(y = expression(Catch/MSY)) +
  coord_cartesian(ylim = c(0, 2.5)) +
  theme_bw(base_size = 8) +
  theme(panel.spacing.x = unit(0, "lines"),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
p_her_SSB <- stats_alt %>%
  filter(stock == "her.27.3a47d" &
           metric == "SSB") %>%
  ggplot(aes(x = OM_label, y = val)) +
  geom_hline(yintercept = 1, colour = "grey") +
  geom_violin(aes(fill = MP_label), size = 0.2, show.legend = FALSE,
              position = position_dodge(width = 0.8), scale = "width") +
  geom_boxplot(aes(group = interaction(OM, MP_label)), 
               position = position_dodge(width = 0.8),
               fill = "white", width = 0.1, size = 0.2,
               outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
               outlier.fill = "transparent") +
  facet_grid(~ OM_group, scales = "free_x", space = "free_x") +
  scale_fill_manual("", values = col_vals) +
  labs(y = expression(SSB/B[MSY])) +
  coord_cartesian(ylim = c(0, 6)) +
  theme_bw(base_size = 8) +
  theme(panel.spacing.x = unit(0, "lines"),
        strip.text.x = element_blank(),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
p_her_risk <- stats_alt %>%
  filter(stock == "her.27.3a47d" & metric == "risk") %>%
  ggplot() +
  geom_hline(yintercept = 0.055, colour = "red") +
  geom_col(data = stats_alt %>%
             filter(stock == "her.27.3a47d" & metric == "risk") %>%
             group_by(MP_label, OM, OM_label, OM_group) %>%
             summarise(val = max(val)),
           aes(x = OM_label, y = val, fill = MP_label), 
           show.legend = TRUE, width = 0.8, colour = "black", size = 0.2,
           position = position_dodge(width = 0.8)) +
  geom_boxplot(aes(x = OM_label, y = val, group = interaction(OM, MP_label)),
               position = position_dodge(width = 0.8),
               fill = "white", width = 0.1, size = 0.2,
               outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
               outlier.fill = "transparent") +
  facet_grid(~ OM_group, scales = "free_x", space = "free_x") +
  scale_fill_manual("", values = col_vals) +
  labs(y = expression(max.~B[lim]~risk)) +
  ylim(c(0, 1)) +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.spacing.x = unit(0, "lines"),
        strip.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none")

### combine catch, SSB & risk plots for each stock
p_ple <- plot_grid(p_ple_catch, p_ple_SSB, 
                   p_ple_risk + theme(legend.position = "none"),
                   ncol = 1, align = "v", 
                   rel_heights = c(1.2, 1, 1.6))
p_cod <- plot_grid(p_cod_catch, p_cod_SSB, 
                   p_cod_risk + theme(legend.position = "none"),
                   ncol = 1, align = "v", 
                   rel_heights = c(1.2, 1, 1.6))
p_her <- plot_grid(p_her_catch, p_her_SSB, 
                   p_her_risk + theme(legend.position = "none"),
                   NULL,
                   ncol = 1, align = "v", 
                   rel_heights = c(1.2, 1, 1.4, 0.2))

p <- plot_grid(NULL,
               plot_grid(p_ple, get_legend(p_cod_risk), 
                         nrow = 1, rel_widths = c(1, 0.22)),
               plot_grid(NULL, NULL, nrow = 1, rel_widths = c(1, 0.6),
                         labels = c("(b) Cod", "(c) Herring"),
                         label_size = 9, label_x = c(-0.025)),
               plot_grid(p_cod, p_her, align = "vh", axis = "tb",
                         nrow = 1, rel_widths = c(1, 0.6)),
               labels = c("(a) Plaice", "", "", ""),
               label_size = 9, label_x = c(-0.025),
               ncol = 1, rel_heights = c(0.07, 1, 0.07, 1), align = "v")
p
ggsave(filename = "output/plots/risk_MPs_all_OMs_all_violin.png", plot = p,
       width = 17, height = 15, units = "cm", dpi = 600, type = "cairo",
       bg = "white")
ggsave(filename = "output/plots/risk_MPs_all_OMs_all_violin.pdf", plot = p,
       width = 17, height = 15, units = "cm", 
       bg = "white")

### ------------------------------------------------------------------------ ###
### violin plots - baseline OM & all stocksMPs  ####
### ------------------------------------------------------------------------ ###
stats_alt <- readRDS("output/MPs_alternative_OMs_stats.rds")
stats_baseline <- readRDS("output/MPs_baseline.rds")
stats_baseline <- stats_baseline %>%
  filter(period == "11-20")
col_vals <- c("2 over 3" = "#9e9ac8", 
              "2 over 3 (XSA)" = "#6a51a3", 
              "rfb (generic)" = "#bdd7e7", 
              "rfb (multiplier)" = "#6baed6", 
              "rfb (all)" = "#2171b5", 
              "hr (generic)" = "#fcae91", 
              "hr (multiplier)" = "#fb6a4a", 
              "hr (all)" = "#cb181d", 
              "ICES MSY" = "#ffff00")
p_catch <- stats_alt %>%
  filter(stock %in% c("ple.27.7e", "cod.27.47d20", "her.27.3a47d") & 
           metric == "catch" & OM == "baseline") %>%
  ggplot(aes(x = MP_label, y = val)) +
  geom_hline(yintercept = 1, colour = "grey") +
  geom_violin(aes(fill = MP_label), size = 0.2, show.legend = FALSE,
              scale = "width") +
  geom_boxplot(fill = "white", width = 0.1, size = 0.2,
               outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
               outlier.fill = "transparent") +
  facet_grid(~ stock_label, scales = "free_x", space = "free_x") +
  scale_fill_manual("", values = col_vals) +
  labs(y = expression(Catch/MSY)) +
  coord_cartesian(ylim = c(0, 2.5)) +
  theme_bw(base_size = 8) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
p_ssb <- stats_alt %>%
  filter(stock %in% c("ple.27.7e", "cod.27.47d20", "her.27.3a47d") & 
           metric == "SSB" & OM == "baseline") %>%
  ggplot(aes(x = MP_label, y = val)) +
  geom_hline(yintercept = 1, colour = "grey") +
  geom_violin(aes(fill = MP_label), size = 0.2, show.legend = FALSE,
              scale = "width") +
  geom_boxplot(fill = "white", width = 0.1, size = 0.2,
               outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
               outlier.fill = "transparent") +
  facet_grid(~ stock_label, scales = "free_x", space = "free_x") +
  scale_fill_manual("", values = col_vals) +
  labs(y = expression(SSB/B[MSY])) +
  coord_cartesian(ylim = c(0, 7)) +
  theme_bw(base_size = 8) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text.x = element_blank())
p_risk <- stats_alt %>%
  filter(stock %in% c("ple.27.7e", "cod.27.47d20", "her.27.3a47d") & 
           metric == "risk" & OM == "baseline") %>%
  ggplot() +
  geom_hline(yintercept = 0.0525, colour = "red", size = 0.75) +
  geom_col(data = stats_alt %>%
             filter(stock %in% c("ple.27.7e", "cod.27.47d20", "her.27.3a47d") &
                      metric == "risk" & OM == "baseline") %>%
             group_by(stock_label, MP_label) %>%
             summarise(val = max(val)),
           aes(x = MP_label, y = val, fill = MP_label),
           show.legend = TRUE, width = 0.8, colour = "black", size = 0.2) +
  geom_boxplot(aes(x = MP_label, y = val),
               position = position_dodge(width = 0.8),
               fill = "white", width = 0.1, size = 0.2,
               outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
               outlier.fill = "transparent") +
  scale_fill_manual("", values = col_vals) +
  facet_grid(~ stock_label, scales = "free_x", space = "free_x") +
  labs(y = expression(max.~B[lim]~risk)) +
  ylim(c(0, NA)) +
  theme_bw(base_size = 8) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        strip.text.x = element_blank())
p_fitness <- stats_alt %>%
  filter(stock %in% c("ple.27.7e", "cod.27.47d20", "her.27.3a47d") &
           OM == "baseline") %>%
  select(stock_label, MP_label, fitness) %>%
  unique() %>%
  ggplot() +
  geom_hline(yintercept = 1, colour = "grey") +
  geom_hline(yintercept = 0, colour = "grey") +
  geom_col(aes(x = MP_label, y = fitness, fill = MP_label), 
           show.legend = TRUE, width = 0.8, colour = "black", size = 0.2) +
  scale_fill_manual("", values = col_vals) +
  facet_grid(~ stock_label, scales = "free_x", space = "free_x") +
  labs(y = expression("fitness "*italic(phi)["MSY-PA"])) +
  scale_y_continuous(#limits = c(-0.27, 1.05), 
                     breaks = c(-0.25, 0, 0.25, 0.5, 0.75, 1)) +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        strip.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none")

p <- plot_grid(p_catch, p_ssb, p_risk, p_fitness,
               ncol = 1, align = "v", 
               rel_heights = c(1.15, 1, 1, 1.45))
p
ggsave(filename = "output/plots/risk_MPs_all_baseline_violin.png", plot = p,
       width = 17, height = 15, units = "cm", dpi = 600, type = "cairo",
       bg = "white")
ggsave(filename = "output/plots/risk_MPs_all_baseline_violin.pdf", plot = p,
       width = 17, height = 15, units = "cm",
       bg = "white")

### ------------------------------------------------------------------------ ###
### plot OM trajectories vs. ICES assessment ####
### ------------------------------------------------------------------------ ###
refpts_cod <- readRDS("input/cod.27.47d20/baseline/1000_100/refpts_mse.rds")
refpts_cod <- iterMedians(refpts_cod)
refpts_ple <- readRDS("input/ple.27.7e/baseline/1000_100/refpts_mse.rds")
refpts_ple <- iterMedians(refpts_ple)
refpts_her <- readRDS("input/her.27.3a47d/baseline/1000_100/refpts_mse.rds")
refpts_her <- iterMedians(refpts_her)

### values from operating model
df_OM <- foreach(stock = c("Plaice", "Cod", "Herring"),
                stock_id = c("ple.27.7e", "cod.27.47d20", "her.27.3a47d"),
                .combine = bind_rows) %do% {
  #browser()
  stk_i <- readRDS(paste0("input/", stock_id, "/baseline/1000_100/stk.rds"))
  ### get metrics
  qnts <- FLQuants(ssb = ssb(stk_i)/1000, fbar = fbar(stk_i))
  ### percentiles
  qnts_perc <- lapply(qnts, quantile, probs = c(0.05, 0.25, 0.5, 0.75, 0.95),
                      na.rm = TRUE)
  qnts_perc <- FLQuants(qnts_perc)
  qnts_perc <- as.data.frame(qnts_perc)
  qnts_perc <- qnts_perc %>% select(year, iter, data, qname) %>%
    pivot_wider(names_from = iter, values_from = data) %>%
    mutate(stock = stock, source = "OM")
  return(qnts_perc)
}
### get ICES assessment summary
df_ICES <- foreach(stock = c("Plaice", "Cod", "Herring"),
                 stock_id = c("ple.27.7e", "cod.27.47d20", "her.27.3a47d"),
                 .combine = bind_rows) %do% {
  #browser()
  smry <- read.csv(paste0("input/", stock_id, "/preparation/",
                          "ices_assessment_summary.csv"))
  names(smry)[1] <- "year"
  smry <- smry %>% 
    select(year, SSB, F) %>%
    mutate(SSB = SSB/1000) %>%
    rename(fbar = F, ssb = SSB) %>%
    pivot_longer(c(ssb, fbar), names_to = "qname") %>%
    mutate(stock = stock, source = "ICES")
  return(smry)
}
df_ICES <- df_ICES %>% filter(year <= 2020)
df_combined <- 
  bind_rows(df_ICES,
            df_OM %>% 
              select(year, qname, `50%`, stock, source) %>%
              rename(value = `50%`)) %>%
  mutate(source = factor(source, levels = c("OM", "ICES"),
                         labels = c("Operating model",
                                    "ICES assessment")))


p_ple_ssb <- ggplot() +
  geom_ribbon(data = df_OM %>%
                filter(stock == "Plaice" & qname == "ssb"),
              aes(x = year, ymin = `5%`, ymax = `95%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_ribbon(data = df_OM %>%
                filter(stock == "Plaice" & qname == "ssb"),
              aes(x = year, ymin = `25%`, ymax = `75%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_line(data = df_combined %>% filter(stock == "Plaice" & qname == "ssb"),
            aes(x = year, y = value, linetype = source, colour = source),
            size = 0.3) +
  scale_color_manual("", values = c("Operating model" = "black",
                                    "ICES assessment" = "red")) +
  scale_linetype_manual("", values = c("Operating model" = "solid",
                                       "ICES assessment" = "2121")) +
  geom_hline(yintercept = c(refpts_ple["Bmsy"])/1000, 
             size = 0.3, linetype = "dashed") + 
  geom_hline(yintercept = c(refpts_ple["Blim"])/1000, 
             size = 0.3, linetype = "dotted") + 
  coord_cartesian(xlim = c(1979, 2021), ylim = c(0, 10.1), expand = FALSE) +
  facet_wrap(~ "Plaice") + 
  labs(y = "SSB [1000t]", x = "Year") +
  theme_bw(base_size = 8) +
  theme(legend.position = c(0.5, 0.75),
        legend.key = element_blank(),
        legend.key.height = unit(0.5, "lines"),
        legend.key.width = unit(0.8, "lines"),
        legend.background = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
p_cod_ssb <- ggplot() +
  geom_ribbon(data = df_OM %>%
                filter(stock == "Cod" & qname == "ssb" & year <= 2020),
              aes(x = year, ymin = `5%`, ymax = `95%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_ribbon(data = df_OM %>%
                filter(stock == "Cod" & qname == "ssb" & year <= 2020),
              aes(x = year, ymin = `25%`, ymax = `75%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_line(data = df_combined %>% 
              filter(stock == "Cod" & qname == "ssb" & year <= 2020),
            aes(x = year, y = value, linetype = source, colour = source),
            size = 0.3) +
  scale_color_manual("", values = c("Operating model" = "black",
                                    "ICES assessment" = "red")) +
  scale_linetype_manual("", values = c("Operating model" = "solid",
                                       "ICES assessment" = "2121")) +
  geom_hline(yintercept = c(refpts_cod["Bmsy"])/1000, 
             size = 0.3, linetype = "dashed") + 
  geom_hline(yintercept = c(refpts_cod["Blim"])/1000, 
             size = 0.3, linetype = "dotted") + 
  coord_cartesian(xlim = c(1962, 2022), ylim = c(0, 270), expand = FALSE) +
  facet_wrap(~ "Cod") + 
  labs(y = "SSB [1000t]", x = "Year") +
  theme_bw(base_size = 8) +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
p_her_ssb <- ggplot() +
  geom_ribbon(data = df_OM %>%
                filter(stock == "Herring" & qname == "ssb" & year <= 2020),
              aes(x = year, ymin = `5%`, ymax = `95%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_ribbon(data = df_OM %>%
                filter(stock == "Herring" & qname == "ssb" & year <= 2020),
              aes(x = year, ymin = `25%`, ymax = `75%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_line(data = df_combined %>% 
              filter(stock == "Herring" & qname == "ssb" & year <= 2020),
            aes(x = year, y = value, linetype = source, colour = source),
            size = 0.3) +
  scale_color_manual("", values = c("Operating model" = "black",
                                    "ICES assessment" = "red")) +
  scale_linetype_manual("", values = c("Operating model" = "solid",
                                       "ICES assessment" = "2121")) +
  geom_hline(yintercept = c(refpts_her["Bmsy"])/1000, 
             size = 0.3, linetype = "dashed") + 
  geom_hline(yintercept = c(refpts_her["Blim"])/1000, 
             size = 0.3, linetype = "dotted") + 
  coord_cartesian(xlim = c(1946, 2022), ylim = c(0, 7000), expand = FALSE) +
  facet_wrap(~ "Herring") + 
  labs(y = "SSB [1000t]", x = "Year") +
  theme_bw(base_size = 8) +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

p_ple_fbar <- ggplot() +
  geom_ribbon(data = df_OM %>%
                filter(stock == "Plaice" & qname == "fbar" & year <= 2020),
              aes(x = year, ymin = `5%`, ymax = `95%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_ribbon(data = df_OM %>%
                filter(stock == "Plaice" & qname == "fbar" & year <= 2020),
              aes(x = year, ymin = `25%`, ymax = `75%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_line(data = df_combined %>% 
              filter(stock == "Plaice" & qname == "fbar" & year <= 2020),
            aes(x = year, y = value, linetype = source, colour = source),
            size = 0.3) +
  scale_color_manual("", values = c("Operating model" = "black",
                                    "ICES assessment" = "red")) +
  scale_linetype_manual("", values = c("Operating model" = "solid",
                                       "ICES assessment" = "2121")) +
  geom_hline(yintercept = c(refpts_ple["Fmsy"]), 
             size = 0.3, linetype = "dashed") + 
  coord_cartesian(xlim = c(1979, 2021), ylim = c(0, 1), expand = FALSE) +
  labs(y = "mean F (ages 3-6)", x = "Year") +
  theme_bw(base_size = 8) +
  theme(legend.position = "none")
p_cod_fbar <- ggplot() +
  geom_ribbon(data = df_OM %>%
                filter(stock == "Cod" & qname == "fbar" & year <= 2020),
              aes(x = year, ymin = `5%`, ymax = `95%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_ribbon(data = df_OM %>%
                filter(stock == "Cod" & qname == "fbar" & year <= 2020),
              aes(x = year, ymin = `25%`, ymax = `75%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_line(data = df_combined %>% 
              filter(stock == "Cod" & qname == "fbar" & year <= 2020),
            aes(x = year, y = value, linetype = source, colour = source),
            size = 0.3) +
  scale_color_manual("", values = c("Operating model" = "black",
                                    "ICES assessment" = "red")) +
  scale_linetype_manual("", values = c("Operating model" = "solid",
                                       "ICES assessment" = "2121")) +
  geom_hline(yintercept = c(refpts_cod["Fmsy"]), 
             size = 0.3, linetype = "dashed") + 
  coord_cartesian(xlim = c(1962, 2021), ylim = c(0, 1.3), expand = FALSE) +
  labs(y = "mean F (ages 2-4)", x = "Year") +
  theme_bw(base_size = 8) +
  theme(legend.position = "none")
p_her_fbar <- ggplot() +
  geom_ribbon(data = df_OM %>%
                filter(stock == "Herring" & qname == "fbar" & year <= 2020),
              aes(x = year, ymin = `5%`, ymax = `95%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_ribbon(data = df_OM %>%
                filter(stock == "Herring" & qname == "fbar" & year <= 2020),
              aes(x = year, ymin = `25%`, ymax = `75%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_line(data = df_combined %>% 
              filter(stock == "Herring" & qname == "fbar" & year <= 2020),
            aes(x = year, y = value, linetype = source, colour = source),
            size = 0.3) +
  scale_color_manual("", values = c("Operating model" = "black",
                                    "ICES assessment" = "red")) +
  scale_linetype_manual("", values = c("Operating model" = "solid",
                                       "ICES assessment" = "2121")) +
  geom_hline(yintercept = c(refpts_her["Fmsy"]), 
             size = 0.3, linetype = "dashed") + 
  coord_cartesian(xlim = c(1946, 2022), ylim = c(0, 1.6), expand = FALSE) +
  labs(y = "mean F (ages 2-6)", x = "Year") +
  theme_bw(base_size = 8) +
  theme(legend.position = "none")

p <- p_ple_ssb + p_cod_ssb + p_her_ssb + 
  p_ple_fbar + p_cod_fbar + p_her_fbar + plot_layout(ncol = 3, widths = 1)
p
ggsave(filename = "output/plots/risk_OM_vs_ICES.png", plot = p, 
       width = 17, height = 8, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/risk_OM_vs_ICES.pdf", plot = p, 
       width = 17, height = 8, units = "cm")


### ------------------------------------------------------------------------ ###
### projections - baseline OM with all MPs and stocks ####
### ------------------------------------------------------------------------ ###
res_alt <- readRDS("output/MPs_alternative_OMs.rds")
res_lookup <- res_alt %>%
  filter(OM == "baseline") %>%
  bind_rows(data.frame(stock = c("ple.27.7e", "cod.27.47d20", "her.27.3a47d"),
                       MP = "history"))
refpts_cod <- readRDS("input/cod.27.47d20/baseline/1000_100/refpts_mse.rds")
refpts_cod <- iterMedians(refpts_cod)
refpts_ple <- readRDS("input/ple.27.7e/baseline/1000_100/refpts_mse.rds")
refpts_ple <- iterMedians(refpts_ple)
refpts_her <- readRDS("input/her.27.3a47d/baseline/1000_100/refpts_mse.rds")
refpts_her <- iterMedians(refpts_her)


proj <- foreach(i = split(res_lookup, f = seq(nrow(res_lookup))),
                .combine = bind_rows) %do% {
  # browser()
  if (identical(i$MP, "history")) {
    stk_i <- readRDS(paste0("input/", i$stock, "/baseline/1000_100/stk.rds"))
    stk_i <- window(stk_i, end = 2020)
  } else {
    ### find OM stock - slot depends on version of mse package
    res_mp <- readRDS(i$file)
    . <- try(stk_i <- res_mp@om@stock, silent = TRUE)
    if (is(., "try-error")) stk_i <- res_mp@stock
  }
  ### get metrics
  qnts <- FLQuants(catch = catch(stk_i)/1000, rec = rec(stk_i)/1000,
                   ssb = ssb(stk_i)/1000, fbar = fbar(stk_i))
  ### percentiles
  qnts_perc <- lapply(qnts, quantile, probs = c(0.05, 0.25, 0.5, 0.75, 0.95),
                      na.rm = TRUE)
  qnts_perc <- FLQuants(qnts_perc)
  qnts_perc <- as.data.frame(qnts_perc)
  qnts_perc <- qnts_perc %>% select(year, iter, data, qname) %>%
    pivot_wider(names_from = iter, values_from = data)
  return(bind_cols(i, qnts_perc))
}
proj <- proj %>%
  mutate(
    MP_label = factor(paste0(MP, "_", optimised), 
                      levels = c("2over3_default", "2over3_XSA_default", 
                                 "rfb_default", "rfb_multiplier", "rfb_all", 
                                 "hr_default", "hr_multiplier", "hr_all", 
                                 "ICES_SAM_default"),
                      labels = c("2 over 3", "2 over 3 (XSA)", 
                                 "rfb (generic)", "rfb (multiplier)", "rfb (all)",
                                 "hr (generic)", "hr (multiplier)", "hr (all)",
                                 "ICES MSY")), .after = "MP") %>%
  mutate(stock_label = factor(stock, 
                              levels = c("ple.27.7e", "cod.27.47d20",
                                         "her.27.3a47d"),
                              labels = c("Plaice", "Cod", "Herring")))
proj_distr <- foreach(i = split(res_lookup, f = seq(nrow(res_lookup))),
                .combine = bind_rows) %do% {
  # browser()
  if (identical(i$MP, "history")) {
    return(NULL)
  } else {
    ### find OM stock - slot depends on version of mse package
    res_mp <- readRDS(i$file)
    . <- try(stk_i <- res_mp@om@stock, silent = TRUE)
    if (is(., "try-error")) stk_i <- res_mp@stock
  }
  ### get metrics & distribution in last year
  qnts <- FLQuants(catch = catch(stk_i)[, ac(2040)]/1000, 
                   rec = rec(stk_i)[, ac(2040)]/1000,
                   ssb = ssb(stk_i)[, ac(2040)]/1000, 
                   fbar = fbar(stk_i)[, ac(2040)])
  qnts <- as.data.frame(qnts)
  qnts <- qnts %>% select(year, iter, data, qname)
  return(bind_cols(i, qnts))
}
proj_distr <- proj_distr %>%
  mutate(
    MP_label = factor(paste0(MP, "_", optimised), 
                      levels = c("2over3_default", "2over3_XSA_default", 
                                 "rfb_default", "rfb_multiplier", "rfb_all", 
                                 "hr_default", "hr_multiplier", "hr_all", 
                                 "ICES_SAM_default"),
                      labels = c("2 over 3", "2 over 3 (XSA)", 
                                 "rfb (generic)", "rfb (multiplier)", "rfb (all)",
                                 "hr (generic)", "hr (multiplier)", "hr (all)",
                                 "ICES MSY")), .after = "MP")


col_vals <- c("2 over 3" = "#8c6bb1", 
              "2 over 3 (XSA)" = "#811b7c", 
              "rfb (generic)" = "#6baed6", 
              "rfb (multiplier)" = "#2271b5", 
              "rfb (all)" = "#08306b", 
              "hr (generic)" = "#66c2a4", 
              "hr (multiplier)" = "#248b45", 
              "hr (all)" = "#00441b", 
              "ICES MSY" = "#fc8d59")
col_vals <- c("2 over 3" = "#9e9ac8", 
              "2 over 3 (XSA)" = "#6a51a3", 
              "rfb (generic)" = "#bdd7e7", 
              "rfb (multiplier)" = "#6baed6", 
              "rfb (all)" = "#2171b5", 
              "hr (generic)" = "#fcae91", 
              "hr (multiplier)" = "#fb6a4a", 
              "hr (all)" = "#cb181d", 
              "ICES MSY" = "#ffff00")
lty_vals <- c("2 over 3" = "solid", 
              "2 over 3 (XSA)" = "1111", 
              "rfb (generic)" = "solid", 
              "rfb (multiplier)" = "3131", 
              "rfb (all)" = "1111", 
              "hr (generic)" = "solid", 
              "hr (multiplier)" = "3131", 
              "hr (all)" = "1111", 
              "ICES MSY" = "solid")

p_ple_catch <- ggplot() +
  geom_ribbon(data = proj %>%
                filter(stock_label == "Plaice" & MP == "history" & 
                         qname == "catch"),
              aes(x = year, ymin = `5%`, ymax = `95%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_ribbon(data = proj %>%
                filter(stock_label == "Plaice" & MP == "history" & 
                         qname == "catch"),
              aes(x = year, ymin = `25%`, ymax = `75%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_ribbon(data = proj %>%
                filter(stock_label == "Plaice" & MP != "history" & 
                         qname == "catch"),
              aes(x = year, ymin = `5%`, ymax = `95%`, fill = MP_label), 
              alpha = 0.05,
              show.legend = FALSE) +
  geom_line(data = proj %>%
              filter(stock_label == "Plaice" & MP == "history" & 
                       qname == "catch"),
            aes(x = year, y = `5%`), 
            alpha = 0.5, size = 0.04,
            show.legend = FALSE) +
  geom_line(data = proj %>%
              filter(stock_label == "Plaice" & MP == "history" & 
                       qname == "catch"),
            aes(x = year, y = `95%`), 
            alpha = 0.5, size = 0.04,
            show.legend = FALSE) +
  geom_line(data = proj %>%
              filter(stock_label == "Plaice" & MP != "history" & 
                       qname == "catch"),
            aes(x = year, y = `5%`, colour = MP_label), 
            alpha = 0.5, size = 0.04,
            show.legend = FALSE) +
  geom_line(data = proj %>%
              filter(stock_label == "Plaice" & MP != "history" & 
                       qname == "catch"),
            aes(x = year, y = `95%`, colour = MP_label), 
            alpha = 0.5, size = 0.04,
            show.legend = FALSE) +
  geom_ribbon(data = proj %>%
                filter(stock_label == "Plaice" & MP != "history" & 
                         qname == "catch"),
              aes(x = year, ymin = `25%`, ymax = `75%`, fill = MP_label), 
              alpha = 0.05,
              show.legend = FALSE) +
  geom_line(data = proj %>%
              filter(stock_label == "Plaice" & MP == "history" & 
                       qname == "catch"),
            aes(x = year, y = `50%`),
            colour = "black", size = 0.4) + 
  geom_line(data = proj %>%
              filter(stock_label == "Plaice" & MP != "history" & 
                       qname == "catch"),
            aes(x = year, y = `50%`, colour = MP_label, linetype = MP_label),
            size = 0.4) +
  geom_hline(yintercept = c(refpts_ple["Cmsy"])/1000, 
             size = 0.4, linetype = "dashed") + 
  scale_alpha(guide = guide_legend(label = FALSE)) +
  scale_colour_manual("", values = col_vals) + 
  scale_fill_manual("", values = col_vals) + 
  scale_linetype_manual("", values = lty_vals) +
  coord_cartesian(xlim = c(2009, 2040), ylim = c(0, 4), expand = FALSE) +
  facet_wrap(~ "Plaice") + 
  labs(y = "Catch [1000t]", x = "Year") +
  theme_bw(base_size = 8) +
  theme(legend.position = "none",
        plot.margin = unit(c(5, 25, 2, 5), "pt"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
p_ple_catch_distr <- proj_distr %>% 
  filter(stock == "ple.27.7e" & MP != "history" & qname == "catch") %>%
  ggplot(aes(x = data, linetype = MP_label, colour = MP_label)) +
  geom_density(aes(y = ..scaled..), show.legend = FALSE, size = 0.2) +
  scale_colour_manual("", values = col_vals) + 
  scale_linetype_manual("", values = lty_vals) +
  theme_bw(base_size = 8) +
  theme_void() +
  coord_flip(xlim = c(0, 4), ylim = c(0, 1.02), expand = FALSE) +
  theme(panel.background = element_rect(fill = "white", color = "white"))
p_ple_catch <- p_ple_catch +
  annotation_custom(grob = ggplotGrob(p_ple_catch_distr),
                    xmin = 2040, xmax = 2043,
                    ymin = 0, ymax = 4)
p_cod_catch <- ggplot() +
  geom_ribbon(data = proj %>%
                filter(stock_label == "Cod" & MP == "history" & 
                         qname == "catch"),
              aes(x = year, ymin = `5%`, ymax = `95%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_ribbon(data = proj %>%
                filter(stock_label == "Cod" & MP == "history" & 
                         qname == "catch"),
              aes(x = year, ymin = `25%`, ymax = `75%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_ribbon(data = proj %>%
                filter(stock_label == "Cod" & MP != "history" & 
                         qname == "catch"),
              aes(x = year, ymin = `5%`, ymax = `95%`, fill = MP_label), 
              alpha = 0.05,
              show.legend = FALSE) +
  geom_line(data = proj %>%
              filter(stock_label == "Cod" & MP == "history" & 
                       qname == "catch"),
            aes(x = year, y = `5%`), 
            alpha = 0.5, size = 0.04,
            show.legend = FALSE) +
  geom_line(data = proj %>%
              filter(stock_label == "Cod" & MP == "history" & 
                       qname == "catch"),
            aes(x = year, y = `95%`), 
            alpha = 0.5, size = 0.04,
            show.legend = FALSE) +
  geom_line(data = proj %>%
              filter(stock_label == "Cod" & MP != "history" & 
                       qname == "catch"),
            aes(x = year, y = `5%`, colour = MP_label), 
            alpha = 0.5, size = 0.04,
            show.legend = FALSE) +
  geom_line(data = proj %>%
              filter(stock_label == "Cod" & MP != "history" & 
                       qname == "catch"),
            aes(x = year, y = `95%`, colour = MP_label), 
            alpha = 0.5, size = 0.04,
            show.legend = FALSE) +
  geom_ribbon(data = proj %>%
                filter(stock_label == "Cod" & MP != "history" & 
                         qname == "catch"),
              aes(x = year, ymin = `25%`, ymax = `75%`, fill = MP_label), 
              alpha = 0.05,
              show.legend = FALSE) +
  geom_line(data = proj %>%
              filter(stock_label == "Cod" & MP == "history" & 
                       qname == "catch"),
            aes(x = year, y = `50%`),
            colour = "black", size = 0.4) + 
  geom_line(data = proj %>%
              filter(stock_label == "Cod" & MP != "history" & 
                       qname == "catch"),
            aes(x = year, y = `50%`, colour = MP_label, linetype = MP_label),
            size = 0.4) +
  geom_hline(yintercept = c(refpts_cod["Cmsy"])/1000, 
             size = 0.4, linetype = "dashed") +
  scale_alpha(guide = guide_legend(label = FALSE)) +
  scale_colour_manual("", values = col_vals) + 
  scale_fill_manual("", values = col_vals) + 
  scale_linetype_manual("", values = lty_vals) +
  xlim(c(2008, NA)) +
  coord_cartesian(xlim = c(2009, 2040), ylim = c(0, 150), expand = FALSE) +
  facet_wrap(~ "Cod") + 
  labs(y = "Catch [1000t]", x = "Year") +
  theme_bw(base_size = 8) +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        plot.margin = unit(c(5, 25, 2, 5), "pt"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
p_cod_catch_distr <- proj_distr %>% 
  filter(stock == "cod.27.47d20" & MP != "history" & qname == "catch") %>%
  ggplot(aes(x = data, colour = MP_label, linetype = MP_label)) +
  geom_density(aes(y = ..scaled..), show.legend = FALSE, size = 0.2) +
  scale_colour_manual("", values = col_vals) + 
  scale_linetype_manual("", values = lty_vals) +
  theme_bw(base_size = 8) +
  theme_void() +
  coord_flip(xlim = c(0, 150), ylim = c(0, 1.02), expand = FALSE) +
  theme(panel.background = element_rect(fill = "white", color = "white"))
p_cod_catch <- p_cod_catch +
  annotation_custom(grob = ggplotGrob(p_cod_catch_distr),
                    xmin = 2040, xmax = 2043,
                    ymin = 0, ymax = 150)
p_her_catch <- ggplot() +
  geom_ribbon(data = proj %>%
                filter(stock_label == "Herring" & MP == "history" & 
                         qname == "catch"),
              aes(x = year, ymin = `5%`, ymax = `95%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_ribbon(data = proj %>%
                filter(stock_label == "Herring" & MP == "history" & 
                         qname == "catch"),
              aes(x = year, ymin = `25%`, ymax = `75%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_ribbon(data = proj %>%
                filter(stock_label == "Herring" & MP != "history" & 
                         qname == "catch"),
              aes(x = year, ymin = `5%`, ymax = `95%`, fill = MP_label), 
              alpha = 0.05,
              show.legend = FALSE) +
  geom_line(data = proj %>%
              filter(stock_label == "Herring" & MP == "history" & 
                       qname == "catch"),
            aes(x = year, y = `5%`), 
            alpha = 0.5, size = 0.04,
            show.legend = FALSE) +
  geom_line(data = proj %>%
              filter(stock_label == "Herring" & MP == "history" & 
                       qname == "catch"),
            aes(x = year, y = `95%`), 
            alpha = 0.5, size = 0.04,
            show.legend = FALSE) +
  geom_line(data = proj %>%
              filter(stock_label == "Herring" & MP != "history" & 
                       qname == "catch"),
            aes(x = year, y = `5%`, colour = MP_label), 
            alpha = 0.5, size = 0.04,
            show.legend = FALSE) +
  geom_line(data = proj %>%
              filter(stock_label == "Herring" & MP != "history" & 
                       qname == "catch"),
            aes(x = year, y = `95%`, colour = MP_label), 
            alpha = 0.5, size = 0.04,
            show.legend = FALSE) +
  geom_ribbon(data = proj %>%
                filter(stock_label == "Herring" & MP != "history" & 
                         qname == "catch"),
              aes(x = year, ymin = `25%`, ymax = `75%`, fill = MP_label), 
              alpha = 0.05,
              show.legend = FALSE) +
  geom_line(data = proj %>%
              filter(stock_label == "Herring" & MP == "history" & 
                       qname == "catch"),
            aes(x = year, y = `50%`),
            colour = "black", size = 0.4) + 
  geom_line(data = proj %>%
              filter(stock_label == "Herring" & MP != "history" & 
                       qname == "catch"),
            aes(x = year, y = `50%`, colour = MP_label, linetype = MP_label),
            size = 0.4) +
  geom_hline(yintercept = c(refpts_her["Cmsy"])/1000, 
             size = 0.4, linetype = "dashed") + 
  scale_alpha(guide = guide_legend(label = FALSE)) +
  scale_colour_manual("", values = col_vals) + 
  scale_fill_manual("", values = col_vals) + 
  scale_linetype_manual("", values = lty_vals) +
  coord_cartesian(xlim = c(2009, 2040), ylim = c(0, 700), expand = FALSE) +
  facet_wrap(~ "Herring") + 
  labs(y = "Catch [1000t]", x = "Year") +
  theme_bw(base_size = 8) +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        plot.margin = unit(c(5, 25, 2, 5), "pt"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
p_her_catch_distr <- proj_distr %>% 
  filter(stock == "her.27.3a47d" & MP != "history" & qname == "catch") %>%
  ggplot(aes(x = data, linetype = MP_label, colour = MP_label)) +
  geom_density(aes(y = ..scaled..), show.legend = FALSE, size = 0.2) +
  scale_colour_manual("", values = col_vals) + 
  scale_linetype_manual("", values = lty_vals) +
  theme_bw(base_size = 8) +
  theme_void() +
  coord_flip(xlim = c(0, 700), ylim = c(0, 1.02), expand = FALSE) +
  theme(panel.background = element_rect(fill = "white", color = "white"))
p_her_catch <- p_her_catch +
  annotation_custom(grob = ggplotGrob(p_her_catch_distr),
                    xmin = 2040, xmax = 2043,
                    ymin = 0, ymax = 700)

p_ple_ssb <- ggplot() +
  geom_ribbon(data = proj %>%
                filter(stock_label == "Plaice" & MP == "history" & 
                         qname == "ssb"),
              aes(x = year, ymin = `5%`, ymax = `95%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_ribbon(data = proj %>%
                filter(stock_label == "Plaice" & MP == "history" & 
                         qname == "ssb"),
              aes(x = year, ymin = `25%`, ymax = `75%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_ribbon(data = proj %>%
                filter(stock_label == "Plaice" & MP != "history" & 
                         qname == "ssb"),
              aes(x = year, ymin = `5%`, ymax = `95%`, fill = MP_label), 
              alpha = 0.05,
              show.legend = FALSE) +
  geom_line(data = proj %>%
              filter(stock_label == "Plaice" & MP == "history" & 
                       qname == "ssb"),
            aes(x = year, y = `5%`), 
            alpha = 0.5, size = 0.04,
            show.legend = FALSE) +
  geom_line(data = proj %>%
              filter(stock_label == "Plaice" & MP == "history" & 
                       qname == "ssb"),
            aes(x = year, y = `95%`), 
            alpha = 0.5, size = 0.04,
            show.legend = FALSE) +
  geom_line(data = proj %>%
              filter(stock_label == "Plaice" & MP != "history" & 
                       qname == "ssb"),
            aes(x = year, y = `5%`, colour = MP_label), 
            alpha = 0.5, size = 0.04,
            show.legend = FALSE) +
  geom_line(data = proj %>%
              filter(stock_label == "Plaice" & MP != "history" & 
                       qname == "ssb"),
            aes(x = year, y = `95%`, colour = MP_label), 
            alpha = 0.5, size = 0.04,
            show.legend = FALSE) +
  geom_ribbon(data = proj %>%
                filter(stock_label == "Plaice" & MP != "history" & 
                         qname == "ssb"),
              aes(x = year, ymin = `25%`, ymax = `75%`, fill = MP_label), 
              alpha = 0.05,
              show.legend = FALSE) +
  geom_line(data = proj %>%
              filter(stock_label == "Plaice" & MP == "history" & 
                       qname == "ssb"),
            aes(x = year, y = `50%`),
            colour = "black", size = 0.4) + 
  geom_line(data = proj %>%
              filter(stock_label == "Plaice" & MP != "history" & 
                       qname == "ssb"),
            aes(x = year, y = `50%`, colour = MP_label, linetype = MP_label),
            size = 0.4) +
  geom_hline(yintercept = c(refpts_ple["Bmsy"])/1000, 
             size = 0.4, linetype = "dashed") +
  geom_hline(yintercept = c(refpts_ple["Blim"])/1000, 
             size = 0.4, linetype = "dotted") +
  scale_alpha(guide = guide_legend(label = FALSE)) +
  scale_colour_manual("", values = col_vals) + 
  scale_fill_manual("", values = col_vals) + 
  scale_linetype_manual("", values = lty_vals) +
  coord_cartesian(xlim = c(2009, 2040), ylim = c(0, 33), expand = FALSE) +
  #facet_wrap(~ "Plaice") + 
  labs(y = "SSB [1000t]", x = "Year") +
  theme_bw(base_size = 8) +
  theme(legend.position = "none",
        plot.margin = unit(c(2, 25, 5, 5), "pt"))
p_ple_ssb_distr <- proj_distr %>% 
  filter(stock == "ple.27.7e" & MP != "history" & qname == "ssb") %>%
  ggplot(aes(x = data, colour = MP_label, linetype = MP_label)) +
  geom_density(aes(y = ..scaled..), show.legend = FALSE, size = 0.2) +
  scale_colour_manual("", values = col_vals) + 
  scale_linetype_manual("", values = lty_vals) +
  theme_bw(base_size = 8) +
  theme_void() +
  coord_flip(xlim = c(0, 33), ylim = c(0, 1.02), expand = FALSE) +
  theme(panel.background = element_rect(fill = "white", color = "white"))
p_ple_ssb <- p_ple_ssb +
  annotation_custom(grob = ggplotGrob(p_ple_ssb_distr),
                    xmin = 2040, xmax = 2043,
                    ymin = 0, ymax = 33)
p_cod_ssb <- ggplot() +
  geom_ribbon(data = proj %>%
                filter(stock_label == "Cod" & MP == "history" & 
                         qname == "ssb"),
              aes(x = year, ymin = `5%`, ymax = `95%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_ribbon(data = proj %>%
                filter(stock_label == "Cod" & MP == "history" & 
                         qname == "ssb"),
              aes(x = year, ymin = `25%`, ymax = `75%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_ribbon(data = proj %>%
                filter(stock_label == "Cod" & MP != "history" & 
                         qname == "ssb"),
              aes(x = year, ymin = `5%`, ymax = `95%`, fill = MP_label), 
              alpha = 0.05,
              show.legend = FALSE) +
  geom_line(data = proj %>%
              filter(stock_label == "Cod" & MP == "history" & 
                       qname == "ssb"),
            aes(x = year, y = `5%`), 
            alpha = 0.5, size = 0.04,
            show.legend = FALSE) +
  geom_line(data = proj %>%
              filter(stock_label == "Cod" & MP == "history" & 
                       qname == "ssb"),
            aes(x = year, y = `95%`), 
            alpha = 0.5, size = 0.04,
            show.legend = FALSE) +
  geom_line(data = proj %>%
              filter(stock_label == "Cod" & MP != "history" & 
                       qname == "ssb"),
            aes(x = year, y = `5%`, colour = MP_label), 
            alpha = 0.5, size = 0.04,
            show.legend = FALSE) +
  geom_line(data = proj %>%
              filter(stock_label == "Cod" & MP != "history" & 
                       qname == "ssb"),
            aes(x = year, y = `95%`, colour = MP_label), 
            alpha = 0.5, size = 0.04,
            show.legend = FALSE) +
  geom_ribbon(data = proj %>%
                filter(stock_label == "Cod" & MP != "history" & 
                         qname == "ssb"),
              aes(x = year, ymin = `25%`, ymax = `75%`, fill = MP_label), 
              alpha = 0.05,
              show.legend = FALSE) +
  geom_line(data = proj %>%
              filter(stock_label == "Cod" & MP == "history" & 
                       qname == "ssb"),
            aes(x = year, y = `50%`),
            colour = "black", size = 0.4) + 
  geom_line(data = proj %>%
              filter(stock_label == "Cod" & MP != "history" & 
                       qname == "ssb"),
            aes(x = year, y = `50%`, colour = MP_label, linetype = MP_label),
            size = 0.4) +
  geom_hline(yintercept = c(refpts_cod["Bmsy"])/1000, 
             size = 0.4, linetype = "dashed") +
  geom_hline(yintercept = c(refpts_cod["Blim"])/1000, 
             size = 0.4, linetype = "dotted") +
  scale_alpha(guide = guide_legend(label = FALSE)) +
  scale_colour_manual("", values = col_vals) + 
  scale_fill_manual("", values = col_vals) + 
  scale_linetype_manual("", values = lty_vals) +
  xlim(c(2008, NA)) +
  coord_cartesian(xlim = c(2009, 2040), ylim = c(0, 620), expand = FALSE) +
  #facet_wrap(~ "Cod") + 
  labs(y = "SSB [1000t]", x = "Year") +
  theme_bw(base_size = 8) +
  theme(legend.position = c(0.3, 0.65),
        legend.key.height = unit(0.5, "lines"),
        legend.background = element_blank(),
        legend.key.width = unit(0.7, "lines"),
        axis.title.y = element_blank(),
        plot.margin = unit(c(2, 25, 5, 5), "pt"))
p_cod_ssb_distr <- proj_distr %>% 
  filter(stock == "cod.27.47d20" & MP != "history" & qname == "ssb") %>%
  ggplot(aes(x = data, colour = MP_label, linetype = MP_label)) +
  geom_density(aes(y = ..scaled..), show.legend = FALSE, size = 0.2) +
  scale_colour_manual("", values = col_vals) + 
  scale_linetype_manual("", values = lty_vals) +
  theme_bw(base_size = 8) +
  theme_void() +
  coord_flip(xlim = c(0, 620), ylim = c(0, 1.02), expand = FALSE) +
  theme(panel.background = element_rect(fill = "white", color = "white"))
p_cod_ssb <- p_cod_ssb +
  annotation_custom(grob = ggplotGrob(p_cod_ssb_distr),
                    xmin = 2040, xmax = 2043,
                    ymin = 0, ymax = 620)
p_her_ssb <- ggplot() +
  geom_ribbon(data = proj %>%
                filter(stock_label == "Herring" & MP == "history" & 
                         qname == "ssb"),
              aes(x = year, ymin = `5%`, ymax = `95%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_ribbon(data = proj %>%
                filter(stock_label == "Herring" & MP == "history" & 
                         qname == "ssb"),
              aes(x = year, ymin = `25%`, ymax = `75%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_ribbon(data = proj %>%
                filter(stock_label == "Herring" & MP != "history" & 
                         qname == "ssb"),
              aes(x = year, ymin = `5%`, ymax = `95%`, fill = MP_label), 
              alpha = 0.05,
              show.legend = FALSE) +
  geom_line(data = proj %>%
              filter(stock_label == "Herring" & MP == "history" & 
                       qname == "ssb"),
            aes(x = year, y = `5%`), 
            alpha = 0.5, size = 0.04,
            show.legend = FALSE) +
  geom_line(data = proj %>%
              filter(stock_label == "Herring" & MP == "history" & 
                       qname == "ssb"),
            aes(x = year, y = `95%`), 
            alpha = 0.5, size = 0.04,
            show.legend = FALSE) +
  geom_line(data = proj %>%
              filter(stock_label == "Herring" & MP != "history" & 
                       qname == "ssb"),
            aes(x = year, y = `5%`, colour = MP_label), 
            alpha = 0.5, size = 0.04,
            show.legend = FALSE) +
  geom_line(data = proj %>%
              filter(stock_label == "Herring" & MP != "history" & 
                       qname == "ssb"),
            aes(x = year, y = `95%`, colour = MP_label), 
            alpha = 0.5, size = 0.04,
            show.legend = FALSE) +
  geom_ribbon(data = proj %>%
                filter(stock_label == "Herring" & MP != "history" & 
                         qname == "ssb"),
              aes(x = year, ymin = `25%`, ymax = `75%`, fill = MP_label), 
              alpha = 0.05,
              show.legend = FALSE) +
  geom_line(data = proj %>%
              filter(stock_label == "Herring" & MP == "history" & 
                       qname == "ssb"),
            aes(x = year, y = `50%`),
            colour = "black", size = 0.4) + 
  geom_line(data = proj %>%
              filter(stock_label == "Herring" & MP != "history" & 
                       qname == "ssb"),
            aes(x = year, y = `50%`, colour = MP_label, linetype = MP_label),
            size = 0.4) +
  geom_hline(yintercept = c(refpts_her["Bmsy"])/1000, 
             size = 0.4, linetype = "dashed") +
  geom_hline(yintercept = c(refpts_her["Blim"])/1000, 
             size = 0.4, linetype = "dotted") +
  scale_alpha(guide = guide_legend(label = FALSE)) +
  scale_colour_manual("", values = col_vals) + 
  scale_fill_manual("", values = col_vals) + 
  scale_linetype_manual("", values = lty_vals) +
  xlim(c(2008, NA)) +
  coord_cartesian(xlim = c(2009, 2040), ylim = c(0, 3500), expand = FALSE) +
  #facet_wrap(~ "Cod") + 
  labs(y = "SSB [1000t]", x = "Year") +
  theme_bw(base_size = 8) +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        plot.margin = unit(c(2, 25, 5, 5), "pt"))
p_her_ssb_distr <- proj_distr %>% 
  filter(stock == "her.27.3a47d" & MP != "history" & qname == "ssb") %>%
  ggplot(aes(x = data, colour = MP_label, linetype = MP_label)) +
  geom_density(aes(y = ..scaled..), show.legend = FALSE, size = 0.2) +
  scale_colour_manual("", values = col_vals) + 
  scale_linetype_manual("", values = lty_vals) +
  theme_bw(base_size = 8) +
  theme_void() +
  coord_flip(xlim = c(0, 3500), ylim = c(0, 1.02), expand = FALSE) +
  theme(panel.background = element_rect(fill = "white", color = "white"))
p_her_ssb <- p_her_ssb +
  annotation_custom(grob = ggplotGrob(p_her_ssb_distr),
                    xmin = 2040, xmax = 2043,
                    ymin = 0, ymax = 3500)


p <- plot_grid(p_ple_catch, p_cod_catch, p_her_catch,
               p_ple_ssb, p_cod_ssb, p_her_ssb,
          ncol = 3, align = "v", axis = "l")
library(patchwork)
p <- p_ple_catch + p_cod_catch + p_her_catch + 
  p_ple_ssb + p_cod_ssb + p_her_ssb + plot_layout(ncol = 3, widths = 1)
p
ggsave(filename = "output/plots/risk_baseline_MPs_projection.png", plot = p, 
       width = 17, height = 8, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/risk_baseline_MPs_projection.pdf", plot = p, 
       width = 17, height = 8, units = "cm")

### ------------------------------------------------------------------------ ###
### projections - rec_failure OM for cod ####
### ------------------------------------------------------------------------ ###
res_alt <- readRDS("output/MPs_alternative_OMs.rds")
res_lookup <- res_alt %>%
  filter(OM == "rec_failure" & stock == "cod.27.47d20" &
         ((MP == "rfb" & optimised == "multiplier") | 
            (MP == "hr" & optimised == "multiplier") |
           MP == "ICES_SAM")) %>%
  bind_rows(data.frame(stock = "cod.27.47d20", MP = "history",
                       optimised = "default"))
refpts_cod <- readRDS("input/cod.27.47d20/baseline/1000_100/refpts_mse.rds")
refpts_cod <- iterMedians(refpts_cod)

### load projections
proj <- foreach(i = split(res_lookup, f = seq(nrow(res_lookup))),
                .combine = bind_rows) %do% {
  # browser()
  if (identical(i$MP, "history")) {
    stk_i <- readRDS(paste0("input/", i$stock, "/baseline/1000_100/stk.rds"))
    stk_i <- window(stk_i, end = 2020)
  } else {
    ### find OM stock - slot depends on version of mse package
    res_mp <- readRDS(i$file)
    . <- try(stk_i <- res_mp@om@stock, silent = TRUE)
    if (is(., "try-error")) stk_i <- res_mp@stock
  }
  ### get metrics
  qnts <- FLQuants(catch = catch(stk_i)/1000, rec = rec(stk_i)/1000,
                   ssb = ssb(stk_i)/1000, fbar = fbar(stk_i))
  ### percentiles
  qnts_perc <- lapply(qnts, quantile, probs = c(0.05, 0.25, 0.5, 0.75, 0.95),
                      na.rm = TRUE)
  qnts_perc <- FLQuants(qnts_perc)
  qnts_perc <- as.data.frame(qnts_perc)
  qnts_perc <- qnts_perc %>% select(year, iter, data, qname) %>%
    pivot_wider(names_from = iter, values_from = data)
  return(bind_cols(i, qnts_perc))
}
proj <- proj %>%
  mutate(
    MP_label = factor(paste0(MP, "_", optimised), 
                      levels = c("rfb_multiplier", "hr_multiplier",
                                 "ICES_SAM_default"),
                      labels = c("rfb (multiplier)", "hr (multiplier)", 
                                 "ICES MSY")), 
    .after = "MP") %>%
  mutate(stock_label = "cod.27.47d20") %>%
  mutate(qname_label = factor(qname,
                              levels = c("rec", "ssb", "catch", "fbar"),
                              labels = c("Recruitment [millions]", 
                                         "SSB [1000 t]", "Catch [1000 t]",
                                         "F (ages 2-4)")))
proj_plot <- proj %>% 
  filter(MP_label %in% c("rfb (multiplier)", "hr (multiplier)", 
                         "ICES MSY") & 
           qname %in% c("catch", "ssb", "rec"))
col_vals <- c("rfb (multiplier)" = "#6baed6", 
              "hr (multiplier)" = "#fb6a4a", 
              "ICES MSY" = "#ffff00")
lty_vals <- c("rfb (multiplier)" = "3131", 
              "hr (multiplier)" = "3131", 
              "ICES MSY" = "solid")

df_area <- data.frame(
  xmin = rep(2021, 3), 
  xmax = rep(2025, 3), 
  ymin = rep(0, 3), 
  ymax = c(300, 390, 47),
  qname_label = factor(c("Recruitment [millions]", 
                         "SSB [1000 t]", "Catch [1000 t]"),
                       levels = c("Recruitment [millions]",
                                  "SSB [1000 t]", "Catch [1000 t]")))
df_text <- data.frame(
  x = rep(2023, 3), y = c(200, 200, 20),
  text = c("recruitment\nfailure:\n2021-2025", "", ""),
  qname_label = factor(c("Recruitment [millions]", 
                         "SSB [1000 t]", "Catch [1000 t]"),
                       levels = c("Recruitment [millions]",
                                  "SSB [1000 t]", "Catch [1000 t]")))
p <- ggplot() +
  geom_rect(data = df_area,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            alpha = 0.2, linetype = 0) +
  geom_text(data = df_text, aes(x = x, y = y, label = text),
            size = 7 * 0.35) +
  geom_line(data = proj_plot,
            aes(x = year, y = `50%`, colour = MP_label, linetype = MP_label)) +
  # facet_wrap(~ qname_label, scales = "free_y", strip.position = "left",
  #            ncol = 1) +
  facet_grid(qname_label ~ "Cod", scales = "free_y", switch = "y") +
  scale_colour_manual("", values = col_vals) +
  scale_linetype_manual("", values = lty_vals) +
  coord_cartesian(xlim = c(2020, 2040.5), ylim = c(0, NA), expand = FALSE) +
  #ylim(c(0, NA)) +
  labs(x = "Year") +
  theme_bw(base_size = 8) +
  theme(strip.placement = "outside",
        strip.background.y = element_blank(),
        strip.text = element_text(size = 8),
        strip.switch.pad.grid = unit(0, "pt"),
        axis.title.y = element_blank(),
        legend.position = c(0.8, 0.8),
        legend.key.height = unit(0.5, "lines"),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.key.width = unit(0.7, "lines"))
p
ggsave(filename = "output/plots/risk_rec_failure_cod_projection.png", plot = p, 
       width = 8, height = 10, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/risk_rec_failure_cod_projection.pdf", plot = p, 
       width = 8, height = 10, units = "cm")


### ------------------------------------------------------------------------ ###
### visualisation of OM MSY values ####
### ------------------------------------------------------------------------ ###

MSY_runs <- foreach(stock = c("ple.27.7e", "cod.27.47d20", "her.27.3a47d"),
                    .combine = bind_rows) %:% 
  foreach(OM = c("baseline", "M_low", "M_high", "M_Gislason",
                 "no_discards_not_hidden",
                 "rec_no_AC", "rec_higher", "M_dd", "M_no_migration"),
          .combine = bind_rows) %do% {#browser()
    file_i <- paste0("input/", stock, "/", OM, "/1000_100/MSY_trace.rds")
    if (!file.exists(file_i)) return(NULL)
    runs_i <- readRDS(file_i)
    runs_i <- bind_rows(runs_i)
    runs_i$stock <- stock
    runs_i$OM <- OM
    return(runs_i)
}
MSY_runs <- MSY_runs %>%
  group_by(stock, OM) %>%
  mutate(MSY = ifelse(catch == max(catch), TRUE, FALSE))


MSY_runs_plot <- MSY_runs %>%
  pivot_longer(c(catch, ssb)) %>%
  mutate(value = value/1000) %>%
  mutate(label = factor(name, levels = c("catch", "ssb"),
                        labels = c("Catch [1000t]", "SSB [1000t]")),
         OM_label = factor(OM, 
                           levels = c("baseline", 
                                      "rec_no_AC", "rec_higher",
                                      "M_low", "M_high", "M_Gislason",
                                      "M_dd", "M_no_migration",
                                      "no_discards_not_hidden"),
                           labels = c("baseline",
                                      "R: no AC", "R: higher",
                                      "M: low", "M: high", "M: Gislason",
                                      "M: dens. dep.", "M: no migration",
                                      "Catch: no discards")))

p_MSY_ple <- MSY_runs_plot %>%
  filter(stock == "ple.27.7e") %>%
  ggplot(aes(x = Ftrgt, y = value)) +
  geom_point(size = 0.4) +
  geom_smooth(span = 0.4,
              se = FALSE, size = 0.2, colour = "blue") +
  geom_vline(data = MSY_runs_plot %>%
               filter(stock == "ple.27.7e" & MSY == TRUE),
             aes(xintercept = Ftrgt), size = 0.4, colour = "red") +
  facet_grid(label ~ OM_label, scales = "free_y", switch = "y") +
  labs(x = "target F", title = "(a) Plaice") +
  scale_x_continuous(breaks = seq(0, 1, 0.2)) +
  theme_bw(base_size = 8) +
  theme(strip.background.y = element_blank(),
        strip.placement = "outside",
        strip.switch.pad.grid = unit(0, "pt"),
        axis.title.y = element_blank(), 
        plot.title = element_text(size = 9, face = "bold"))
p_MSY_cod <- MSY_runs_plot %>%
  filter(stock == "cod.27.47d20") %>%
  ggplot(aes(x = Ftrgt, y = value)) +
  geom_point(size = 0.4) +
  stat_smooth(span = 0.2, 
              se = FALSE, size = 0.2, colour = "blue") +
  geom_vline(data = MSY_runs_plot %>%
               filter(stock == "cod.27.47d20" & MSY == TRUE),
             aes(xintercept = Ftrgt), size = 0.4, colour = "red") +
  facet_grid(label ~ OM_label, scales = "free_y", switch = "y") +
  labs(x = "target F", title = "(b) Cod") +
  scale_x_continuous(breaks = seq(0, 1, 0.2)) +
  theme_bw(base_size = 8) +
  theme(strip.background.y = element_blank(),
        strip.placement = "outside",
        strip.switch.pad.grid = unit(0, "pt"),
        axis.title.y = element_blank(), 
        plot.title = element_text(size = 9, face = "bold"))
p_MSY_her <- MSY_runs_plot %>%
  filter(stock == "her.27.3a47d") %>%
  ggplot(aes(x = Ftrgt, y = value)) +
  geom_point(size = 0.4) +
  stat_smooth(span = 0.25, 
              se = FALSE, size = 0.2, colour = "blue") +
  geom_vline(data = MSY_runs_plot %>%
               filter(stock == "her.27.3a47d" & MSY == TRUE),
             aes(xintercept = Ftrgt), size = 0.4, colour = "red") +
  facet_grid(label ~ OM_label, scales = "free_y", switch = "y") +
  labs(x = "target F", title = "(c) Herring") +
  scale_x_continuous(breaks = seq(0, 1, 0.2)) +
  theme_bw(base_size = 8) +
  theme(strip.background.y = element_blank(),
        strip.placement = "outside",
        strip.switch.pad.grid = unit(0, "pt"),
        axis.title.y = element_blank(), 
        plot.title = element_text(size = 9, face = "bold"))

p_MSY <- plot_grid(p_MSY_ple, 
          plot_grid(p_MSY_cod, p_MSY_her, nrow = 1, rel_widths = c(1, 0.5)),
          ncol = 1)

ggsave(filename = "output/plots/risk_MSY_all.png", plot = p_MSY, 
       width = 17, height = 10, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/risk_MSY_all.pdf", plot = p_MSY, 
       width = 17, height = 10, units = "cm")

### get Blim
Blims <- foreach(stock = c("ple.27.7e", "cod.27.47d20", "her.27.3a47d"),
                    .combine = bind_rows) %:% 
  foreach(OM = c("baseline", "M_low", "M_high", "M_Gislason",
                 "no_discards_not_hidden",
                 "rec_no_AC", "rec_higher", "M_dd", "M_no_migration"),
          .combine = bind_rows) %do% {#browser()
  file_i <- paste0("input/", stock, "/", OM, "/1000_100/refpts_mse.rds")
  if (!file.exists(file_i)) return(NULL)
  i <- readRDS(file_i)
  return(data.frame(stock = stock, OM = OM ,Blim = median(c(i["Blim"]))))
}

### table
### MSY values
MSY_runs %>%
  filter(MSY == TRUE) %>%
  select(-tsb) %>% unique() %>%
  ### unfished values
  full_join(MSY_runs %>%
              filter(Ftrgt == 0) %>%
              select(ssb, rec, stock, OM) %>%
              rename(B0 = ssb, R0 = rec)) %>%
  ### Blim
  full_join(Blims) %>%
  select(stock, OM, B0, R0, Ftrgt, catch, ssb, rec, Blim) %>%
  rename(FMSY = Ftrgt, MSY = catch, BMSY = ssb, RMSY = rec) %>%
  write.csv(file = "input/OM_refpts.csv", row.names = FALSE)

























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





