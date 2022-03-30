### ------------------------------------------------------------------------ ###
### Operating model figures ####
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
### plot OM trajectories vs. ICES assessment - baseline OMs ####
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
  coord_cartesian(xlim = c(1979, 2021), ylim = c(0, 10.5), expand = FALSE) +
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
  coord_cartesian(xlim = c(1962, 2021.5), ylim = c(0, 270), expand = FALSE) +
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
  coord_cartesian(xlim = c(1962, 2021.5), ylim = c(0, 1.3), expand = FALSE) +
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
### plot OM trajectories vs. ICES assessment - alternative OMs ####
### ------------------------------------------------------------------------ ###
### find alternative OMs
res_alt <- readRDS("output/MPs_alternative_OMs.rds")
res_alt_OMs <- res_alt %>%
  select(stock, OM_group, OM) %>%
  unique() %>%
  mutate(
    OM_label = factor(OM, 
      levels = c("baseline", "M_low", "M_high", "M_Gislason",
                 "M_dd", "M_no_migration", "no_discards",
                 "rec_higher", "rec_no_AC", "rec_failure"),
      labels = c("baseline", "low", "high", "Gislason",
                 "dens. dep.", "no migration", "no discards",
                 "higher", "no AC", "failure")), .after = "OM") %>%
  mutate(
    OM_label2 = factor(OM_label, 
       levels = c("baseline", "low", "high", "Gislason",
                  "dens. dep.", "no migration", "no discards",
                  "higher", "no AC", "failure"),
       labels = c("baseline", 
                  "M: low", "M: high", "M: Gislason",
                  "M: dens. dep.", 
                  "M: no migration", "Catch: no discards",
                  "Rec: higher", "Rec: no AC", "Rec: failure")),
    .after = "OM_label") %>%
  mutate(OM_group = factor(OM_group, c("baseline", "M", "Catch", "Rec"))) %>%
  mutate(stock_label = factor(stock, 
                              levels = c("ple.27.7e", "cod.27.47d20",
                                         "her.27.3a47d"),
                              labels = c("Plaice", "Cod", "Herring")))

### load reference points for all stocks and OMs
df_refpts <- foreach(i = seq(nrow(res_alt_OMs)), .combine = bind_rows) %do% {
  OM_i <- res_alt_OMs$OM[i]
  if (identical(OM_i, "rec_failure")) OM_i <- "baseline"
  stock_i <- res_alt_OMs$stock[i]
  refpts_i <- readRDS(paste0("input/", stock_i, "/",
                             OM_i, "/1000_100/refpts_mse.rds"))
  refpts_i <- iterMedians(refpts_i)
  res_alt_OMs[i, ] %>%
    mutate(Fmsy = c(refpts_i["Fmsy"]),
           Bmsy = c(refpts_i["Bmsy"]),
           Cmsy = c(refpts_i["Cmsy"]),
           Blim = c(refpts_i["Blim"])) %>%
    rename(stock_id = stock, stock = stock_label)
}

### get values from operating models
df_OM <- foreach(stock_id = res_alt_OMs$stock,
                 OM = res_alt_OMs$OM,
                 OM_label2 = res_alt_OMs$OM_label2,
                 stock = res_alt_OMs$stock_label,
                 .combine = bind_rows) %do% {
  #browser()
  stk_i <- readRDS(paste0("input/", stock_id, "/", 
                          ifelse(identical(OM, "rec_failure"),
                                 "baseline", OM), 
                          "/1000_100/stk.rds"))
  ### get metrics
  qnts <- FLQuants(ssb = ssb(stk_i)/1000, fbar = fbar(stk_i),
                   catch = catch(stk_i)/1000, rec = rec(stk_i)/1000)
  if (isTRUE(OM_label2 == "Catch: no discards"))
    qnts$catch <- landings(stk_i)/1000
  ### percentiles
  qnts_perc <- lapply(qnts, quantile, probs = c(0.05, 0.25, 0.5, 0.75, 0.95),
                       na.rm = TRUE)
  qnts_perc <- FLQuants(qnts_perc)
  qnts_perc <- as.data.frame(qnts_perc)
  qnts_perc <- qnts_perc %>% select(year, iter, data, qname) %>%
    pivot_wider(names_from = iter, values_from = data) %>%
    mutate(stock = stock, stock_id = stock_id,
           OM = OM, OM_label2 = OM_label2,
           source = "OM")
  return(qnts_perc)
}
df_OM <- df_OM %>% filter(year <= 2020)
### get ICES assessment summary
df_ICES <- foreach(stock = c("Plaice", "Cod", "Herring"),
                   stock_id = c("ple.27.7e", "cod.27.47d20", "her.27.3a47d"),
                   .combine = bind_rows) %do% {
  #browser()
  smry <- read.csv(paste0("input/", stock_id, "/preparation/",
                          "ices_assessment_summary.csv"))
  names(smry)[1] <- "year"
  smry <- smry %>% 
    select(year, SSB, F, Catches, Landings, Recruitment) %>%
    mutate(SSB = SSB/1000, Catches = Catches/1000, Landings = Landings/1000,
           Recruitment = Recruitment/1000) %>%
    rename(fbar = F, ssb = SSB, catch = Catches, landings = Landings,
           rec = Recruitment) %>%
    pivot_longer(c(ssb, fbar, catch, landings, rec), names_to = "qname") %>%
    mutate(stock = stock, stock_id = stock_id, source = "ICES")
  return(smry)
}
df_ICES <- df_ICES %>% filter(year <= 2020)
### add OM levels
df_ICES <- df_ICES %>%
  full_join(df_OM %>% 
              select(stock, stock_id, OM, OM_label2) %>% 
              unique())
# ### use landings for plaice scenario
# df_ICES$value[df_ICES$stock == "Plaice" & 
#                 df_ICES$OM_label2 == "Catch: no discards" &
#                 df_ICES$qname == "catch"] <- 
#   df_ICES$value[df_ICES$stock == "Plaice" & 
#                   df_ICES$OM_label2 == "Catch: no discards" &
#                   df_ICES$qname == "landings"]

df_combined <- 
  bind_rows(df_ICES %>%
              select(year, qname, value, source, stock, stock_id,
                     OM, OM_label2),
            df_OM %>% 
              select(year, qname, value = `50%`, stock, stock_id, source, 
                     OM, OM_label2)) %>%
  mutate(source = factor(source, levels = c("OM", "ICES"),
                         labels = c("Operating model",
                                    "ICES assessment")))


p_ple_catch <- ggplot() +
  geom_ribbon(data = df_OM %>%
                filter(stock == "Plaice" & qname == "catch"),
              aes(x = year, ymin = `5%`, ymax = `95%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_ribbon(data = df_OM %>%
                filter(stock == "Plaice" & qname == "catch"),
              aes(x = year, ymin = `25%`, ymax = `75%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_line(data = df_combined %>% 
              filter(stock == "Plaice" & qname == "catch"),
            aes(x = year, y = value, linetype = source, colour = source),
            size = 0.3) +
  geom_hline(data = df_refpts %>%
               filter(stock == "Plaice"),
             aes(yintercept = Cmsy/1000),
             size = 0.3, linetype = "dashed") +
  scale_color_manual("", values = c("Operating model" = "black",
                                    "ICES assessment" = "red")) +
  scale_linetype_manual("", values = c("Operating model" = "solid",
                                       "ICES assessment" = "2121")) +
  coord_cartesian(xlim = c(1979, 2021), ylim = c(0, 3.2), expand = FALSE) +
  facet_wrap(~ OM_label2, ncol = 1, strip.position = "right") + 
  labs(y = "Catch [1000t]", x = "Year") +
  theme_bw(base_size = 8) +
  theme(legend.position = "none",
        legend.key = element_blank(),
        legend.key.height = unit(0.5, "lines"),
        legend.key.width = unit(0.8, "lines"),
        legend.background = element_blank(),
        strip.text.y = element_blank())
p_ple_rec <- ggplot() +
  geom_ribbon(data = df_OM %>%
                filter(stock == "Plaice" & qname == "rec"),
              aes(x = year, ymin = `5%`, ymax = `95%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_ribbon(data = df_OM %>%
                filter(stock == "Plaice" & qname == "rec"),
              aes(x = year, ymin = `25%`, ymax = `75%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_line(data = df_combined %>% 
              filter(stock == "Plaice" & qname == "rec"),
            aes(x = year, y = value, linetype = source, colour = source),
            size = 0.3) +
  scale_color_manual("", values = c("Operating model" = "black",
                                    "ICES assessment" = "red")) +
  scale_linetype_manual("", values = c("Operating model" = "solid",
                                       "ICES assessment" = "2121")) +
  coord_cartesian(xlim = c(1979, 2021), ylim = c(0, 39), expand = FALSE) +
  facet_wrap(~ OM_label2, ncol = 1, strip.position = "right") + 
  labs(y = "Recruitment [1000s]", x = "Year") +
  theme_bw(base_size = 8) +
  theme(legend.position = c(0.5, 0.98),
        legend.key = element_blank(),
        legend.key.height = unit(0.5, "lines"),
        legend.key.width = unit(0.8, "lines"),
        legend.background = element_blank(),
        strip.text.y = element_blank())
p_ple_ssb <- ggplot() +
  geom_ribbon(data = df_OM %>%
                filter(stock == "Plaice" & qname == "ssb"),
              aes(x = year, ymin = `5%`, ymax = `95%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_ribbon(data = df_OM %>%
                filter(stock == "Plaice" & qname == "ssb"),
              aes(x = year, ymin = `25%`, ymax = `75%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_line(data = df_combined %>% 
              filter(stock == "Plaice" & qname == "ssb"),
            aes(x = year, y = value, linetype = source, colour = source),
            size = 0.3) +
  scale_color_manual("", values = c("Operating model" = "black",
                                    "ICES assessment" = "red")) +
  scale_linetype_manual("", values = c("Operating model" = "solid",
                                       "ICES assessment" = "2121")) +
  geom_hline(data = df_refpts %>%
               filter(stock == "Plaice"),
             aes(yintercept = Bmsy/1000),
             size = 0.3, linetype = "dashed") +
  geom_hline(data = df_refpts %>%
               filter(stock == "Plaice"),
             aes(yintercept = Blim/1000),
             size = 0.3, linetype = "dotted") +
  coord_cartesian(xlim = c(1979, 2021), ylim = c(0, 23), expand = FALSE) +
  facet_wrap(~ OM_label2, ncol = 1, strip.position = "right") + 
  labs(y = "SSB [1000t]", x = "Year") +
  theme_bw(base_size = 8) +
  theme(legend.position = "none",
        legend.key = element_blank(),
        legend.key.height = unit(0.5, "lines"),
        legend.key.width = unit(0.8, "lines"),
        legend.background = element_blank(),
        strip.text.y = element_blank())
p_ple_fbar <- ggplot() +
  geom_ribbon(data = df_OM %>%
                filter(stock == "Plaice" & qname == "fbar"),
              aes(x = year, ymin = `5%`, ymax = `95%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_ribbon(data = df_OM %>%
                filter(stock == "Plaice" & qname == "fbar"),
              aes(x = year, ymin = `25%`, ymax = `75%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_line(data = df_combined %>% 
              filter(stock == "Plaice" & qname == "fbar"),
            aes(x = year, y = value, linetype = source, colour = source),
            size = 0.3) +
  scale_color_manual("", values = c("Operating model" = "black",
                                    "ICES assessment" = "red")) +
  scale_linetype_manual("", values = c("Operating model" = "solid",
                                       "ICES assessment" = "2121")) +
  geom_hline(data = df_refpts %>%
               filter(stock == "Plaice"),
             aes(yintercept = Fmsy),
             size = 0.3, linetype = "dashed") +
  coord_cartesian(xlim = c(1979, 2021), ylim = c(0, 0.9), expand = FALSE) +
  facet_wrap(~ OM_label2, ncol = 1, strip.position = "right") + 
  labs(y = "mean F (ages 3-6)", x = "Year") +
  theme_bw(base_size = 8) +
  theme(legend.position = "none",
        legend.key = element_blank(),
        legend.key.height = unit(0.5, "lines"),
        legend.key.width = unit(0.8, "lines"),
        legend.background = element_blank())

p_cod_catch <- ggplot() +
  geom_ribbon(data = df_OM %>%
                filter(stock == "Cod" & qname == "catch"),
              aes(x = year, ymin = `5%`, ymax = `95%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_ribbon(data = df_OM %>%
                filter(stock == "Cod" & qname == "catch"),
              aes(x = year, ymin = `25%`, ymax = `75%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_line(data = df_combined %>% 
              filter(stock == "Cod" & qname == "catch"),
            aes(x = year, y = value, linetype = source, colour = source),
            size = 0.3) +
  geom_hline(data = df_refpts %>%
               filter(stock == "Cod"),
             aes(yintercept = Cmsy/1000),
             size = 0.3, linetype = "dashed") +
  scale_color_manual("", values = c("Operating model" = "black",
                                    "ICES assessment" = "red")) +
  scale_linetype_manual("", values = c("Operating model" = "solid",
                                       "ICES assessment" = "2121")) +
  coord_cartesian(xlim = c(1962, 2021.5), ylim = c(0, 610), expand = FALSE) +
  facet_wrap(~ OM_label2, ncol = 1, strip.position = "right") + 
  labs(y = "Catch [1000t]", x = "Year") +
  theme_bw(base_size = 8) +
  theme(legend.position = "none",
        legend.key = element_blank(),
        legend.key.height = unit(0.5, "lines"),
        legend.key.width = unit(0.8, "lines"),
        legend.background = element_blank(),
        strip.text.y = element_blank())
p_cod_rec <- ggplot() +
  geom_ribbon(data = df_OM %>%
                filter(stock == "Cod" & qname == "rec"),
              aes(x = year, ymin = `5%`/1000, ymax = `95%`/1000), alpha = 0.15,
              show.legend = FALSE) +
  geom_ribbon(data = df_OM %>%
                filter(stock == "Cod" & qname == "rec"),
              aes(x = year, ymin = `25%`/1000, ymax = `75%`/1000), alpha = 0.15,
              show.legend = FALSE) +
  geom_line(data = df_combined %>% 
              filter(stock == "Cod" & qname == "rec"),
            aes(x = year, y = value/1000, linetype = source, colour = source),
            size = 0.3) +
  scale_color_manual("", values = c("Operating model" = "black",
                                    "ICES assessment" = "red")) +
  scale_linetype_manual("", values = c("Operating model" = "solid",
                                       "ICES assessment" = "2121")) +
  coord_cartesian(xlim = c(1962, 2021.5), ylim = c(0, 3.4), expand = FALSE) +
  facet_wrap(~ OM_label2, ncol = 1, strip.position = "right") + 
  labs(y = "Recruitment [millions]", x = "Year") +
  theme_bw(base_size = 8) +
  theme(legend.position = "none",
        legend.key = element_blank(),
        legend.key.height = unit(0.5, "lines"),
        legend.key.width = unit(0.8, "lines"),
        legend.background = element_blank(),
        strip.text.y = element_blank())
p_cod_ssb <- ggplot() +
  geom_ribbon(data = df_OM %>%
                filter(stock == "Cod" & qname == "ssb"),
              aes(x = year, ymin = `5%`, ymax = `95%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_ribbon(data = df_OM %>%
                filter(stock == "Cod" & qname == "ssb"),
              aes(x = year, ymin = `25%`, ymax = `75%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_line(data = df_combined %>% 
              filter(stock == "Cod" & qname == "ssb"),
            aes(x = year, y = value, linetype = source, colour = source),
            size = 0.3) +
  scale_color_manual("", values = c("Operating model" = "black",
                                    "ICES assessment" = "red")) +
  scale_linetype_manual("", values = c("Operating model" = "solid",
                                       "ICES assessment" = "2121")) +
  geom_hline(data = df_refpts %>%
               filter(stock == "Cod"),
             aes(yintercept = Bmsy/1000),
             size = 0.3, linetype = "dashed") +
  geom_hline(data = df_refpts %>%
               filter(stock == "Cod"),
             aes(yintercept = Blim/1000),
             size = 0.3, linetype = "dotted") +
  coord_cartesian(xlim = c(1962, 2021.5), ylim = c(0, 290), expand = FALSE) +
  facet_wrap(~ OM_label2, ncol = 1, strip.position = "right") + 
  labs(y = "SSB [1000t]", x = "Year") +
  theme_bw(base_size = 8) +
  theme(legend.position = "none",
        legend.key = element_blank(),
        legend.key.height = unit(0.5, "lines"),
        legend.key.width = unit(0.8, "lines"),
        legend.background = element_blank(),
        strip.text.y = element_blank())
p_cod_fbar <- ggplot() +
  geom_ribbon(data = df_OM %>%
                filter(stock == "Cod" & qname == "fbar"),
              aes(x = year, ymin = `5%`, ymax = `95%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_ribbon(data = df_OM %>%
                filter(stock == "Cod" & qname == "fbar"),
              aes(x = year, ymin = `25%`, ymax = `75%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_line(data = df_combined %>% 
              filter(stock == "Cod" & qname == "fbar"),
            aes(x = year, y = value, linetype = source, colour = source),
            size = 0.3) +
  scale_color_manual("", values = c("Operating model" = "black",
                                    "ICES assessment" = "red")) +
  scale_linetype_manual("", values = c("Operating model" = "solid",
                                       "ICES assessment" = "2121")) +
  geom_hline(data = df_refpts %>%
               filter(stock == "Cod"),
             aes(yintercept = Fmsy),
             size = 0.3, linetype = "dashed") +
  coord_cartesian(xlim = c(1962, 2021.5), ylim = c(0, 1.3), expand = FALSE) +
  facet_wrap(~ OM_label2, ncol = 1, strip.position = "right") + 
  labs(y = "mean F (ages 2-4)", x = "Year") +
  theme_bw(base_size = 8) +
  theme(legend.position = "none",
        legend.key = element_blank(),
        legend.key.height = unit(0.5, "lines"),
        legend.key.width = unit(0.8, "lines"),
        legend.background = element_blank())

p_her_catch <- ggplot() +
  geom_ribbon(data = df_OM %>%
                filter(stock == "Herring" & qname == "catch"),
              aes(x = year, ymin = `5%`, ymax = `95%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_ribbon(data = df_OM %>%
                filter(stock == "Herring" & qname == "catch"),
              aes(x = year, ymin = `25%`, ymax = `75%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_line(data = df_combined %>% 
              filter(stock == "Herring" & qname == "catch"),
            aes(x = year, y = value, linetype = source, colour = source),
            size = 0.3) +
  geom_hline(data = df_refpts %>%
               filter(stock == "Herring"),
             aes(yintercept = Cmsy/1000),
             size = 0.3, linetype = "dashed") +
  scale_color_manual("", values = c("Operating model" = "black",
                                    "ICES assessment" = "red")) +
  scale_linetype_manual("", values = c("Operating model" = "solid",
                                       "ICES assessment" = "2121")) +
  coord_cartesian(xlim = c(1946, 2022), ylim = c(0, 1450), expand = FALSE) +
  facet_wrap(~ OM_label2, ncol = 1, strip.position = "right") + 
  labs(y = "Catch [1000t]", x = "Year") +
  theme_bw(base_size = 8) +
  theme(legend.position = "none",
        legend.key = element_blank(),
        legend.key.height = unit(0.5, "lines"),
        legend.key.width = unit(0.8, "lines"),
        legend.background = element_blank(),
        strip.text.y = element_blank())
p_her_rec <- ggplot() +
  geom_ribbon(data = df_OM %>%
                filter(stock == "Herring" & qname == "rec"),
              aes(x = year, ymin = `5%`/1000, ymax = `95%`/1000), alpha = 0.15,
              show.legend = FALSE) +
  geom_ribbon(data = df_OM %>%
                filter(stock == "Herring" & qname == "rec"),
              aes(x = year, ymin = `25%`/1000, ymax = `75%`/1000), alpha = 0.15,
              show.legend = FALSE) +
  geom_line(data = df_combined %>% 
              filter(stock == "Herring" & qname == "rec"),
            aes(x = year, y = value/1000, linetype = source, colour = source),
            size = 0.3) +
  scale_color_manual("", values = c("Operating model" = "black",
                                    "ICES assessment" = "red")) +
  scale_linetype_manual("", values = c("Operating model" = "solid",
                                       "ICES assessment" = "2121")) +
  coord_cartesian(xlim = c(1946, 2022), ylim = c(0, 95), expand = FALSE) +
  facet_wrap(~ OM_label2, ncol = 1, strip.position = "right") + 
  labs(y = "Recruitment [1000s]", x = "Year") +
  theme_bw(base_size = 8) +
  theme(legend.position = "none",
        legend.key = element_blank(),
        legend.key.height = unit(0.5, "lines"),
        legend.key.width = unit(0.8, "lines"),
        legend.background = element_blank(),
        strip.text.y = element_blank())
p_her_ssb <- ggplot() +
  geom_ribbon(data = df_OM %>%
                filter(stock == "Herring" & qname == "ssb"),
              aes(x = year, ymin = `5%`, ymax = `95%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_ribbon(data = df_OM %>%
                filter(stock == "Herring" & qname == "ssb"),
              aes(x = year, ymin = `25%`, ymax = `75%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_line(data = df_combined %>% 
              filter(stock == "Herring" & qname == "ssb"),
            aes(x = year, y = value, linetype = source, colour = source),
            size = 0.3) +
  scale_color_manual("", values = c("Operating model" = "black",
                                    "ICES assessment" = "red")) +
  scale_linetype_manual("", values = c("Operating model" = "solid",
                                       "ICES assessment" = "2121")) +
  geom_hline(data = df_refpts %>%
               filter(stock == "Herring"),
             aes(yintercept = Bmsy/1000),
             size = 0.3, linetype = "dashed") +
  geom_hline(data = df_refpts %>%
               filter(stock == "Herring"),
             aes(yintercept = Blim/1000),
             size = 0.3, linetype = "dotted") +
  coord_cartesian(xlim = c(1946, 2022), ylim = c(0, 7200), expand = FALSE) +
  facet_wrap(~ OM_label2, ncol = 1, strip.position = "right") + 
  labs(y = "SSB [1000t]", x = "Year") +
  theme_bw(base_size = 8) +
  theme(legend.position = "none",
        legend.key = element_blank(),
        legend.key.height = unit(0.5, "lines"),
        legend.key.width = unit(0.8, "lines"),
        legend.background = element_blank(),
        strip.text.y = element_blank())
p_her_fbar <- ggplot() +
  geom_ribbon(data = df_OM %>%
                filter(stock == "Herring" & qname == "fbar"),
              aes(x = year, ymin = `5%`, ymax = `95%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_ribbon(data = df_OM %>%
                filter(stock == "Herring" & qname == "fbar"),
              aes(x = year, ymin = `25%`, ymax = `75%`), alpha = 0.15,
              show.legend = FALSE) +
  geom_line(data = df_combined %>% 
              filter(stock == "Herring" & qname == "fbar"),
            aes(x = year, y = value, linetype = source, colour = source),
            size = 0.3) +
  scale_color_manual("", values = c("Operating model" = "black",
                                    "ICES assessment" = "red")) +
  scale_linetype_manual("", values = c("Operating model" = "solid",
                                       "ICES assessment" = "2121")) +
  geom_hline(data = df_refpts %>%
               filter(stock == "Herring"),
             aes(yintercept = Fmsy),
             size = 0.3, linetype = "dashed") +
  coord_cartesian(xlim = c(1946, 2022), ylim = c(0, 1.59), expand = FALSE) +
  facet_wrap(~ OM_label2, ncol = 1, strip.position = "right") + 
  labs(y = "mean F (ages 2-6)", x = "Year") +
  theme_bw(base_size = 8) +
  theme(legend.position = "none",
        legend.key = element_blank(),
        legend.key.height = unit(0.5, "lines"),
        legend.key.width = unit(0.8, "lines"),
        legend.background = element_blank())


p_ple <- p_ple_catch + 
  ggtitle(label = "(a) Plaice") + 
  theme(plot.title = element_text(face = "bold")) +
  p_ple_rec + p_ple_ssb + p_ple_fbar + 
  plot_layout(nrow = 1)
ggsave(filename = "output/plots/risk_OM_vs_ICES_ple.png", plot = p_ple, 
       width = 16, height = 16, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/risk_OM_vs_ICES_ple.pdf", plot = p_ple, 
       width = 16, height = 16, units = "cm")

p_cod_her <- p_cod_catch + 
  ggtitle(label = "(b) Cod") + 
  theme(plot.title = element_text(face = "bold")) +
  p_cod_rec + p_cod_ssb + p_cod_fbar + 
  p_her_catch + 
  ggtitle(label = "(c) Herring") + 
  theme(plot.title = element_text(face = "bold")) +
  p_her_rec + p_her_ssb + p_her_fbar + 
  plot_layout(nrow = 2, heights = c(5, 3))
ggsave(filename = "output/plots/risk_OM_vs_ICES_cod_her.png", plot = p_cod_her, 
       width = 16, height = 18, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/risk_OM_vs_ICES_cod_her.pdf", plot = p_cod_her, 
       width = 16, height = 18, units = "cm")

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
                   plot_grid(p_MSY_cod, p_MSY_her, nrow = 1, 
                             rel_widths = c(1, 0.5)),
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







