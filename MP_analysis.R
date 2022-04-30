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
# res_alt <- readRDS("output/MPs_alternative_OMs.rds")

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
                 "ICES MSY")), .after = "MP") %>%
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
### violin plots - all stocks/OMs/MPs - grouped by OM ####
### ------------------------------------------------------------------------ ###
stats_alt <- readRDS("output/MPs_alternative_OMs_stats.rds")

col_vals <- c("2 over 3" = "#A0A0A0",#"#9e9ac8", 
              "2 over 3 (XSA)" = "#606060",#"#6a51a3", 
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
### violin plots - all stocks/OMs/MPs - grouped by MP ####
### ------------------------------------------------------------------------ ###
stats_alt <- readRDS("output/MPs_alternative_OMs_stats.rds")
### create OM labels
stats_alt_MP <- stats_alt %>%
  mutate(OM_label_axis = paste0(OM_group, ": ", OM_label)) %>%
  mutate(OM_label_axis = factor(OM_label_axis, 
    levels = c("baseline: baseline", "M: low", "M: high", "M: Gislason",
               "M: dens. dep.", "M: no migration",
               "Catch: no discards",
               "Rec: no AC", "Rec: higher", "Rec: failure"),
    labels = c("baseline", "M: low", "M:high", "M: Gislason",
               "M: dens. dep.", "M: no migr.",
               "Catch: no disc.",
               "R: no AC", "R: higher", "R: failure")))

col_vals <- c("2 over 3" = "#A0A0A0",#"#9e9ac8", 
              "2 over 3 (XSA)" = "#606060",#"#6a51a3", 
              "rfb (generic)" = "#bdd7e7", 
              "rfb (multiplier)" = "#6baed6", 
              "rfb (all)" = "#2171b5", 
              "hr (generic)" = "#fcae91", 
              "hr (multiplier)" = "#fb6a4a", 
              "hr (all)" = "#cb181d", 
              "ICES MSY" = "#ffff00")


p_ple_risk <- stats_alt_MP %>%
  filter(stock == "ple.27.7e" & metric == "risk") %>%
  ggplot() +
  geom_hline(yintercept = 0.055, colour = "red") +
  geom_col(data = stats_alt_MP %>%
             filter(stock == "ple.27.7e" & metric == "risk") %>%
             group_by(MP_label, OM, OM_label_axis, OM_group) %>%
             summarise(val = max(val)),
           aes(x = OM_label_axis, y = val, fill = MP_label), 
           show.legend = TRUE, width = 0.8, colour = "black", size = 0.2,
           position = position_dodge(width = 0.8)) +
  geom_boxplot(aes(x = OM_label_axis, y = val, group = interaction(OM, MP_label)),
               position = position_dodge(width = 0.8),
               fill = "white", width = 0.1, size = 0.2,
               outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
               outlier.fill = "transparent") +
  scale_fill_manual("", values = col_vals) +
  facet_grid(~ MP_label, scales = "free_x", space = "free_x") +
  labs(y = expression(max.~B[lim]~risk), title = "(a) Plaice") +
  ylim(c(0, NA)) +
  theme_bw(base_size = 8) +
  theme(panel.spacing.x = unit(0, "lines"),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(face = "bold"),
        legend.position = "none")
p_ple_catch <- stats_alt_MP %>%
  filter(stock == "ple.27.7e" &
           metric == "catch") %>%
  ggplot(aes(x = OM_label_axis, y = val)) +
  geom_hline(yintercept = 1, colour = "grey") +
  geom_violin(aes(fill = MP_label), size = 0.2, show.legend = FALSE,
              position = position_dodge(width = 0.8), scale = "width") +
  geom_boxplot(aes(group = interaction(OM, MP_label)), 
               position = position_dodge(width = 0.8),
               fill = "white", width = 0.1, size = 0.2,
               outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
               outlier.fill = "transparent") +
  facet_grid(~ MP_label, scales = "free_x", space = "free_x") +
  scale_fill_manual("", values = col_vals) +
  labs(y = expression(Catch/MSY)) +
  coord_cartesian(ylim = c(0, 2.5)) +
  theme_bw(base_size = 8) +
  theme(panel.spacing.x = unit(0, "lines"),
        strip.text.x = element_blank(),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
p_ple_SSB <- stats_alt_MP %>%
  filter(stock == "ple.27.7e" &
           metric == "SSB") %>%
  ggplot(aes(x = OM_label_axis, y = val)) +
  geom_hline(yintercept = 1, colour = "grey") +
  geom_violin(aes(fill = MP_label), size = 0.2, show.legend = FALSE,
              position = position_dodge(width = 0.8), scale = "width") +
  geom_boxplot(aes(group = interaction(OM, MP_label)), 
               position = position_dodge(width = 0.8),
               fill = "white", width = 0.1, size = 0.2,
               outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
               outlier.fill = "transparent") +
  facet_grid(~ MP_label, scales = "free_x", space = "free_x") +
  scale_fill_manual("", values = col_vals) +
  labs(y = expression(SSB/B[MSY])) +
  coord_cartesian(ylim = c(0, 3.5)) +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5),
        panel.spacing.x = unit(0, "lines"),
        strip.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none",#c(0.25, 0.75),
        legend.key = element_blank(),
        legend.key.height = unit(0.4, "lines"),
        legend.key.width = unit(0.3, "lines"),
        legend.background = element_blank())

p_cod_risk <- stats_alt_MP %>%
  filter(stock == "cod.27.47d" & metric == "risk") %>%
  ggplot() +
  geom_hline(yintercept = 0.055, colour = "red") +
  geom_col(data = stats_alt_MP %>%
             filter(stock == "cod.27.47d20" & metric == "risk") %>%
             group_by(MP_label, OM, OM_label_axis, OM_group) %>%
             summarise(val = max(val)),
           aes(x = OM_label_axis, y = val, fill = MP_label), 
           show.legend = TRUE, width = 0.8, colour = "black", size = 0.2,
           position = position_dodge(width = 0.8)) +
  geom_boxplot(aes(x = OM_label_axis, y = val, group = interaction(OM, MP_label)),
               position = position_dodge(width = 0.8),
               fill = "white", width = 0.1, size = 0.2,
               outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
               outlier.fill = "transparent") +
  scale_fill_manual("", values = col_vals) +
  facet_grid(~ MP_label, scales = "free_x", space = "free_x") +
  labs(y = expression(max.~B[lim]~risk), title = "(b) Cod") +
  ylim(c(0, NA)) +
  theme_bw(base_size = 8) +
  theme(panel.spacing.x = unit(0, "lines"),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(face = "bold"),
        legend.position = "none")
p_cod_catch <- stats_alt_MP %>%
  filter(stock == "cod.27.47d20" &
           metric == "catch") %>%
  ggplot(aes(x = OM_label_axis, y = val)) +
  geom_hline(yintercept = 1, colour = "grey") +
  geom_violin(aes(fill = MP_label), size = 0.2, show.legend = FALSE,
              position = position_dodge(width = 0.8), scale = "width") +
  geom_boxplot(aes(group = interaction(OM, MP_label)), 
               position = position_dodge(width = 0.8),
               fill = "white", width = 0.1, size = 0.2,
               outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
               outlier.fill = "transparent") +
  facet_grid(~ MP_label, scales = "free_x", space = "free_x") +
  scale_fill_manual("", values = col_vals) +
  labs(y = expression(Catch/MSY)) +
  coord_cartesian(ylim = c(0, 2.5)) +
  theme_bw(base_size = 8) +
  theme(panel.spacing.x = unit(0, "lines"),
        strip.text.x = element_blank(),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
p_cod_SSB <- stats_alt_MP %>%
  filter(stock == "cod.27.47d20" &
           metric == "SSB") %>%
  ggplot(aes(x = OM_label_axis, y = val)) +
  geom_hline(yintercept = 1, colour = "grey") +
  geom_violin(aes(fill = MP_label), size = 0.2, show.legend = FALSE,
              position = position_dodge(width = 0.8), scale = "width") +
  geom_boxplot(aes(group = interaction(OM, MP_label)), 
               position = position_dodge(width = 0.8),
               fill = "white", width = 0.1, size = 0.2,
               outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
               outlier.fill = "transparent") +
  facet_grid(~ MP_label, scales = "free_x", space = "free_x") +
  scale_fill_manual("", values = col_vals) +
  labs(y = expression(SSB/B[MSY])) +
  coord_cartesian(ylim = c(0, 6)) +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5),
        panel.spacing.x = unit(0, "lines"),
        strip.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none",#c(0.25, 0.75),
        legend.key = element_blank(),
        legend.key.height = unit(0.4, "lines"),
        legend.key.width = unit(0.3, "lines"),
        legend.background = element_blank())

p_her_risk <- stats_alt_MP %>%
  filter(stock == "her.27.3a47d" & metric == "risk") %>%
  ggplot() +
  geom_hline(yintercept = 0.055, colour = "red") +
  geom_col(data = stats_alt_MP %>%
             filter(stock == "her.27.3a47d" & metric == "risk") %>%
             group_by(MP_label, OM, OM_label_axis, OM_group) %>%
             summarise(val = max(val)),
           aes(x = OM_label_axis, y = val, fill = MP_label), 
           show.legend = TRUE, width = 0.8, colour = "black", size = 0.2,
           position = position_dodge(width = 0.8)) +
  geom_boxplot(aes(x = OM_label_axis, y = val, group = interaction(OM, MP_label)),
               position = position_dodge(width = 0.8),
               fill = "white", width = 0.1, size = 0.2,
               outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
               outlier.fill = "transparent") +
  scale_fill_manual("", values = col_vals) +
  facet_grid(~ MP_label, scales = "free_x", space = "free_x") +
  labs(y = expression(max.~B[lim]~risk), title = "(c) Herring") +
  ylim(c(0, NA)) +
  theme_bw(base_size = 8) +
  theme(panel.spacing.x = unit(0, "lines"),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(face = "bold"),
        legend.position = "none")
p_her_catch <- stats_alt_MP %>%
  filter(stock == "her.27.3a47d" &
           metric == "catch") %>%
  ggplot(aes(x = OM_label_axis, y = val)) +
  geom_hline(yintercept = 1, colour = "grey") +
  geom_violin(aes(fill = MP_label), size = 0.2, show.legend = FALSE,
              position = position_dodge(width = 0.8), scale = "width") +
  geom_boxplot(aes(group = interaction(OM, MP_label)), 
               position = position_dodge(width = 0.8),
               fill = "white", width = 0.1, size = 0.2,
               outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
               outlier.fill = "transparent") +
  facet_grid(~ MP_label, scales = "free_x", space = "free_x") +
  scale_fill_manual("", values = col_vals) +
  labs(y = expression(Catch/MSY)) +
  coord_cartesian(ylim = c(0, 2.5)) +
  theme_bw(base_size = 8) +
  theme(panel.spacing.x = unit(0, "lines"),
        strip.text.x = element_blank(),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
p_her_SSB <- stats_alt_MP %>%
  filter(stock == "her.27.3a47d" &
           metric == "SSB") %>%
  ggplot(aes(x = OM_label_axis, y = val)) +
  geom_hline(yintercept = 1, colour = "grey") +
  geom_violin(aes(fill = MP_label), size = 0.2, show.legend = FALSE,
              position = position_dodge(width = 0.8), scale = "width") +
  geom_boxplot(aes(group = interaction(OM, MP_label)), 
               position = position_dodge(width = 0.8),
               fill = "white", width = 0.1, size = 0.2,
               outlier.size = 0.35, outlier.shape = 21, outlier.stroke = 0.2,
               outlier.fill = "transparent") +
  facet_grid(~ MP_label, scales = "free_x", space = "free_x") +
  scale_fill_manual("", values = col_vals) +
  labs(y = expression(SSB/B[MSY])) +
  coord_cartesian(ylim = c(0, 6)) +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5),
        panel.spacing.x = unit(0, "lines"),
        strip.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none",#c(0.25, 0.75),
        legend.key = element_blank(),
        legend.key.height = unit(0.4, "lines"),
        legend.key.width = unit(0.3, "lines"),
        legend.background = element_blank())

layout <- "AA\nBB\nCC\n#D\n#E\n#F\n#G\n#H\n#I\n"
p <- p_ple_risk + p_ple_catch + p_ple_SSB +
  p_cod_risk + p_cod_catch + p_cod_SSB +
  p_her_risk + p_her_catch + p_her_SSB +
  plot_layout(design = layout, widths = c(0.195, 1))
p
ggsave(filename = "output/plots/risk_MPs_all_OMs_all_violin_MP.png", plot = p,
       width = 18, height = 23, units = "cm", dpi = 600, type = "cairo",
       bg = "white")
ggsave(filename = "output/plots/risk_MPs_all_OMs_all_violin_MP.pdf", plot = p,
       width = 18, height = 23, units = "cm", 
       bg = "white")


### ------------------------------------------------------------------------ ###
### violin plots - baseline OM & all stocksMPs  ####
### ------------------------------------------------------------------------ ###
stats_alt <- readRDS("output/MPs_alternative_OMs_stats.rds")
stats_baseline <- readRDS("output/MPs_baseline.rds")
stats_baseline <- stats_baseline %>%
  filter(period == "11-20")
col_vals <- c("2 over 3" = "#A0A0A0",#"#9e9ac8", 
              "2 over 3 (XSA)" = "#606060",#"#6a51a3", 
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
        axis.ticks.x = element_blank(),
        legend.position = "none",
        strip.text.x = element_blank())
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
        strip.text.x = element_blank(),
        legend.position = "none")
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
        legend.position = "none")
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

p <- p_risk + p_catch + p_ssb + p_fitness + plot_layout(ncol = 1)

p
ggsave(filename = "output/plots/risk_MPs_all_baseline_violin.png", plot = p,
       width = 18, height = 16, units = "cm", dpi = 600, type = "cairo",
       bg = "white")
ggsave(filename = "output/plots/risk_MPs_all_baseline_violin.pdf", plot = p,
       width = 18, height = 16, units = "cm",
       bg = "white")

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


col_vals <- c("2 over 3" = "#A0A0A0",#"#9e9ac8", 
              "2 over 3 (XSA)" = "#606060",#"#6a51a3", 
              "rfb (generic)" = "#bdd7e7", 
              "rfb (multiplier)" = "#6baed6", 
              "rfb (all)" = "#2171b5", 
              "hr (generic)" = "#fcae91", 
              "hr (multiplier)" = "#fb6a4a", 
              "hr (all)" = "#cb181d", 
              "ICES MSY" = "#ffff00")
lty_vals <- c("2 over 3" = "solid", 
              "2 over 3 (XSA)" = "2121", 
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


# p <- plot_grid(p_ple_catch, p_cod_catch, p_her_catch,
#                p_ple_ssb, p_cod_ssb, p_her_ssb,
#           ncol = 3, align = "v", axis = "l")
p <- p_ple_catch + p_cod_catch + p_her_catch + 
  p_ple_ssb + p_cod_ssb + p_her_ssb + plot_layout(ncol = 3, widths = 1)
p
ggsave(filename = "output/plots/risk_baseline_MPs_projection.png", plot = p, 
       width = 18, height = 8, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/risk_baseline_MPs_projection.pdf", plot = p, 
       width = 18, height = 8, units = "cm")

### ------------------------------------------------------------------------ ###
### projections - rec_failure OM for cod ####
### ------------------------------------------------------------------------ ###
res_alt <- readRDS("output/MPs_alternative_OMs.rds")
res_lookup <- res_alt %>%
  filter(OM == "rec_failure" & stock == "cod.27.47d20" &
         ((MP == "rfb" & optimised == "multiplier") | 
            (MP == "hr" & optimised == "multiplier") |
           MP == "ICES_SAM"))
refpts_cod <- readRDS("input/cod.27.47d20/baseline/1000_100/refpts_mse.rds")
refpts_cod <- iterMedians(refpts_cod)

### load projections
proj <- foreach(i = split(res_lookup, f = seq(nrow(res_lookup))),
                .combine = bind_rows) %do% {
  # browser()
  ### find OM stock - slot depends on version of mse package
  res_mp <- readRDS(i$file)
  . <- try(stk_i <- res_mp@om@stock, silent = TRUE)
  if (is(., "try-error")) stk_i <- res_mp@stock
  ### get metrics
  qnts <- FLQuants(catch = catch(stk_i)/1000, rec = rec(stk_i)/1000,
                   ssb = ssb(stk_i)/1000, fbar = fbar(stk_i),
                   risk = ssb(stk_i) %=% NA_real_)
  ### percentiles
  qnts_perc <- lapply(qnts, quantile, probs = c(0.05, 0.25, 0.5, 0.75, 0.95),
                      na.rm = TRUE)
  qnts_perc <- FLQuants(qnts_perc)
  ### annual risks
  qnts_perc$risk[,,,,, "50%"] <- apply(qnts$ssb < c(refpts_cod["Blim"]/1000), 2,
                                       mean, na.rm = TRUE)
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
                              levels = c("rec", "ssb", "risk", "catch", "fbar"),
                              labels = c("'Recruitment\n  [millions]'", 
                                         "'SSB [1000 t]'", "B[lim]~'risk'",
                                         "'Catch [1000 t]'",
                                         "'F (ages 2-4)'")))
proj_plot <- proj %>% 
  filter(MP_label %in% c("rfb (multiplier)", "hr (multiplier)", 
                         "ICES MSY") & 
           qname %in% c("catch", "ssb", "rec", "risk"))
col_vals <- c("rfb (multiplier)" = "#6baed6", 
              "hr (multiplier)" = "#fb6a4a", 
              "ICES MSY" = "#ffff00")
lty_vals <- c("rfb (multiplier)" = "3131", 
              "hr (multiplier)" = "3131", 
              "ICES MSY" = "solid")

df_area <- data.frame(
  xmin = rep(2021, 4), 
  xmax = rep(2025, 4), 
  ymin = rep(0, 4), 
  ymax = c(300, 390, 1, 47),
  qname_label = factor(c("'Recruitment\n  [millions]'", 
                         "'SSB [1000 t]'", "B[lim]~'risk'",
                         "'Catch [1000 t]'"),
                       levels = c("'Recruitment\n  [millions]'",
                                  "'SSB [1000 t]'", "B[lim]~'risk'",
                                  "'Catch [1000 t]'")))
df_text <- data.frame(
  x = rep(2023, 4), y = c(200, 200, 1, 20),
  text = c("recruitment\nfailure:\n2021-2025", "", "", ""),
  qname_label = factor(c("'Recruitment\n  [millions]'", 
                         "'SSB [1000 t]'", "B[lim]~'risk'",
                         "'Catch [1000 t]'"),
                       levels = c("'Recruitment\n  [millions]'",
                                  "'SSB [1000 t]'", "B[lim]~'risk'",
                                  "'Catch [1000 t]'")))
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
  facet_grid(qname_label ~ "Cod", scales = "free_y", switch = "y", 
             labeller = label_parsed) +
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
        legend.position = c(0.8, 0.85),
        legend.key.height = unit(0.5, "lines"),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.key.width = unit(0.7, "lines"),
        strip.text.y = element_text(margin = unit(c(0, 0, 0, 8), "pt"))
        )
        #strip.text.y = unit(c(4, 4, 4, 20), "pt"))
p
ggsave(filename = "output/plots/risk_rec_failure_cod_projection.png", plot = p, 
       width = 8, height = 10, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/risk_rec_failure_cod_projection.pdf", plot = p, 
       width = 8, height = 10, units = "cm")


