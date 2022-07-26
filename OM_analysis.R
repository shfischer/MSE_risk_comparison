### ------------------------------------------------------------------------ ###
### Operating model figures ####
### ------------------------------------------------------------------------ ###

library(mse)
library(GA)
library(tidyr)
library(dplyr)
library(tibble)
library(cowplot)
library(patchwork)
library(ggplot2)
library(foreach)
library(stockassessment)
library(FLfse)
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
       width = 18, height = 8.5, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/risk_OM_vs_ICES.pdf", plot = p, 
       width = 18, height = 8.5, units = "cm")

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
  coord_cartesian(xlim = c(1946, 2022), ylim = c(0, 245), expand = FALSE) +
  facet_wrap(~ OM_label2, ncol = 1, strip.position = "right") + 
  labs(y = "Recruitment [millions]", x = "Year") +
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
  coord_cartesian(xlim = c(1946, 2022), ylim = c(0, 11000), expand = FALSE) +
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
  plot_layout(nrow = 2, heights = c(5, 5))
ggsave(filename = "output/plots/risk_OM_vs_ICES_cod_her.png", plot = p_cod_her, 
       width = 16, height = 22, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/risk_OM_vs_ICES_cod_her.pdf", plot = p_cod_her, 
       width = 16, height = 22, units = "cm")

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

p_MSY <- p_MSY_ple + p_MSY_cod + p_MSY_her +
  plot_layout(design = "AA\nB#\nC#\n", widths = c(4, 2))
ggsave(filename = "output/plots/risk_MSY_all.png", plot = p_MSY, 
       width = 17, height = 15, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/risk_MSY_all.pdf", plot = p_MSY, 
       width = 17, height = 15, units = "cm")

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
### Plaice: recruitment model and residual visualisation ####
### ------------------------------------------------------------------------ ###

### input data, including discard estimates
stk_data <- readRDS("input/ple.27.7e/preparation/model_input_stk_d.RDS")
idx_data <- readRDS("input/ple.27.7e/preparation/model_input_idx.RDS")
### use configuration similar to accepted XSA assessment
SAM_conf <- list(keyLogFpar = 
                   matrix(data = c(rep(-1, 9),
                                   0:5, 5, -1, -1,
                                   6:11, 11, 11, -1),
                          ncol = 9, nrow = 3, byrow = TRUE))

fit <- FLR_SAM(stk_data, idx_data, conf = SAM_conf)
fit_stk <- SAM2FLStock(object = fit, stk = stk_data)
sr <- as.FLSR(fit_stk, model = "bevholtSV")
sr <- fmle(sr, method = 'L-BFGS-B', fixed = list(), 
           control = list(trace = 0))
sr_params <- abPars("bevholt", s = params(sr)["s"], v = params(sr)["v"], 
                    spr0 = params(sr)["spr0"])
sr_params <- list(a = c(sr_params$a), b = c(sr_params$b))
sr_model <- function(ssb, a, b) {(a*ssb)/(b + ssb)}

### data.frame with data
df_sr <- data.frame(year = as.numeric(dimnames(sr)$year),
                    ssb = c(ssb(sr)), 
                    rec = c(rec(sr)), 
                    fitted = c(fitted(sr)),
                    residuals = c(residuals(sr)))
### get density
dens <- density(x = df_sr$residuals)
df_dens <- data.frame(x = dens$x, y = dens$y)

p_res <- df_sr %>%
  ggplot(aes(x = ssb, y = rec)) +
  geom_linerange(aes(x = ssb, ymin = rec, ymax = fitted), 
                 colour = "blue", size = 0.2, 
                 linetype = "2121") +
  geom_point(size = 0.5) +
  geom_function(fun = sr_model, args = sr_params, size = 0.4) +
  scale_colour_manual(values = c(observation = "black", model = "black")) +
  scale_y_continuous(breaks = c(0, 5000, 10000, 15000),
                     labels = c(0, 5, 10, 15),
                     limits = c(0, 18000), expand = c(0, 0)) +
  scale_x_continuous(breaks = c(0, 2000, 4000, 6000), 
                     labels = c(0, 2, 4, 6), 
                     limits = c(0, 7500), expand = c(0, 0)) +
  labs(x = "SSB [1000t]", y = "Recruitment [1000s]") +
  theme_bw(base_size = 8)
p_hist <- df_sr %>%
  ggplot(aes(residuals)) +
  geom_histogram(bins = 15, colour = "black", fill = "white", size = 0.4) +
  geom_line(aes(x = x, y = y * 10, colour = "kernel"), 
            data = df_dens, show.legend = TRUE, size = 0.4) + 
  scale_y_continuous(sec.axis = sec_axis(~./10, name = "Density"),
                     limits = c(-0.05, 10), expand = c(0, 0)) +
  scale_colour_manual("", values = c(kernel = "red"), 
                      labels = c("kernel density")) +
  labs(y = "Residual count", x = "Log residuals") +
  theme_bw(base_size = 8) +
  theme(axis.title.y.right = element_text(angle = 90),
        legend.position = c(0.8, 0.8),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.key.width = unit(0.8, "lines"))

p_res / p_hist + plot_annotation(tag_levels = "a", 
                                 tag_prefix = "(", tag_suffix = ")")  &
  theme(plot.tag = element_text(face = "bold"))
ggsave(filename = "output/plots/OM/OM_ple_rec.png", 
       width = 9, height = 9, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/OM/OM_ple_rec.pdf", 
       width = 9, height = 9, units = "cm", dpi = 600)


### ------------------------------------------------------------------------ ###
### Recruitment models of baseline OMs (all stocks) ####
### ------------------------------------------------------------------------ ###

### plaice
sr_ple <- readRDS("input/ple.27.7e/baseline/1000_100/sr.rds")
df_sr_ple <- data.frame(year = as.numeric(dimnames(ssb(sr_ple))$year),
                        ssb = c(iterMedians(ssb(sr_ple))), 
                        rec = c(iterMedians(rec(sr_ple))))
sr_pars_ple <- abPars("bevholt", s = params(sr_ple)["s"], 
                      v = params(sr_ple)["v"], 
                      spr0 = params(sr_ple)["spr0"])
sr_pars_ple_med <- list(a = c(iterMedians(sr_pars_ple$a)), 
                    b = c(iterMedians(sr_pars_ple$b)))
sr_model_bevholt <- function(ssb, a, b) {(a*ssb)/(b + ssb)}
ssbs_ple <- seq(from = 0, to = 7500, by = 10)
sr_ple_df <- lapply(1:1000, function(x) {
  data.frame(ssb = ssbs_ple,
             rec = sr_model_bevholt(ssb = ssbs_ple, 
                                    a = c(sr_pars_ple$a[, x]),
                                    b = c(sr_pars_ple$b[, x])),
             iter = x)
})
sr_ple_df <- bind_rows(sr_ple_df)

df_sr_ple %>%
  ggplot(aes(x = ssb, y = rec)) +
  geom_point() +
  xlim(c(0, NA)) + ylim(c(0, NA)) +
  geom_function(fun = sr_model_bevholt, args = sr_pars_ple_med, size = 0.4)

p_ple <- df_sr_ple %>%
  ggplot(aes(x = ssb, y = rec)) +
  geom_line(data = sr_ple_df %>% filter(iter %in% 1:1000),
            aes(x = ssb, y = rec, group = iter), size = 0.1, alpha = 0.05) +
  geom_function(fun = sr_model_bevholt, args = sr_pars_ple_med, size = 0.4,
                colour = "red") +
  geom_point(size = 0.3, colour = "red") +
  scale_x_continuous("SSB [1000t]", breaks = c(0, 2000, 4000, 6000),
                     labels = c(0, 2, 4, 6)) +
  scale_y_continuous("Recruitment [1000s]", 
                     breaks = c(0, 5000, 10000, 15000),
                     labels = c(0, 5, 10, 15)) +
  coord_cartesian(xlim = c(0, 7200), ylim = c(0, 17000), expand = FALSE) +
  facet_wrap(~ "Plaice") +
  theme_bw(base_size = 8)

### cod
sr_cod <- readRDS("input/cod.27.47d20/baseline/1000_100/sr.rds")
df_sr_cod <- data.frame(year = as.numeric(dimnames(ssb(sr_cod))$year),
                        ssb = c(iterMedians(ssb(sr_cod))), 
                        rec = c(iterMedians(rec(sr_cod))))
sr_pars_cod <- params(sr_cod)
sr_pars_cod_med <- list(a = c(iterMedians(sr_pars_cod$a)), 
                        b = c(iterMedians(sr_pars_cod$b)))
sr_model_segreg <- function(ssb, a, b) {ifelse(ssb <= b, a * ssb, a * b)}
#ssbs_cod <- seq(from = 0, to = 110000, by = 100)
ssbs_cod <- c(0, c(sr_pars_cod$b), 110000) ## only breakpoint and min/max
sr_cod_df <- lapply(1:1000, function(x) {
  data.frame(ssb = ssbs_cod,
             rec = sr_model_segreg(ssb = ssbs_cod, 
                                   a = c(sr_pars_cod$a[, x]),
                                   b = c(sr_pars_cod$b[, x])),
             iter = x)
})
sr_cod_df <- bind_rows(sr_cod_df)

df_sr_cod %>%
  ggplot(aes(x = ssb, y = rec)) +
  geom_point() +
  xlim(c(0, NA)) + ylim(c(0, NA)) +
  geom_function(fun = sr_model_segreg, args = sr_pars_cod_med, size = 0.4)

p_cod <- df_sr_cod %>%
  ggplot(aes(x = ssb, y = rec)) +
  geom_line(data = sr_cod_df %>% filter(iter %in% 1:1000),
            aes(x = ssb, y = rec, group = iter), size = 0.1, alpha = 0.05) +
  geom_function(fun = sr_model_segreg, args = sr_pars_cod_med, size = 0.4,
                colour = "red") +
  geom_point(size = 0.3, colour = "red") +
  scale_x_continuous("SSB [1000t]", breaks = c(0, 25000, 50000, 75000, 100000),
                     labels = c(0, 25, 50, 75, 100)) +
  scale_y_continuous("Recruitment [1000s]",
                     breaks = c(0, 100000, 200000, 300000, 400000),
                     labels = c(0, 100, 200, 300, 400)) +
  coord_cartesian(xlim = c(0, 105000), ylim = c(0, 500000), expand = FALSE) +
  facet_wrap(~ "Cod") +
  theme_bw(base_size = 8)


### herring
sr_her <- readRDS("input/her.27.3a47d/baseline/1000_100/sr.rds")
df_sr_her <- data.frame(year = as.numeric(dimnames(ssb(sr_her))$year),
                        ssb = c(iterMedians(ssb(sr_her))), 
                        rec = c(iterMedians(rec(sr_her))))
sr_pars_her <- params(sr_her)
sr_pars_her_med <- list(a = c(iterMedians(sr_pars_her$a)), 
                        b = c(iterMedians(sr_pars_her$b)))
sr_model_segreg <- function(ssb, a, b) {ifelse(ssb <= b, a * ssb, a * b)}
#ssbs_her <- seq(from = 0, to = 2600000, by = 1000)
ssbs_her <- c(0, unique(c(sr_pars_her$b)), 26000000)
sr_her_df <- lapply(1:1000, function(x) {
  data.frame(ssb = ssbs_her,
             rec = sr_model_segreg(ssb = ssbs_her, 
                                   a = c(sr_pars_her$a[, x]),
                                   b = c(sr_pars_her$b[, x])),
             iter = x)
})
sr_her_df <- bind_rows(sr_her_df)

df_sr_her %>%
  ggplot(aes(x = ssb, y = rec)) +
  geom_point() +
  xlim(c(0, NA)) + ylim(c(0, NA)) +
  geom_function(fun = sr_model_segreg, args = sr_pars_her_med, size = 0.4,
                n = 100)

p_her <- df_sr_her %>%
  ggplot(aes(x = ssb, y = rec)) +
  geom_line(data = sr_her_df %>% filter(iter %in% 1:1000),
            aes(x = ssb, y = rec, group = iter), size = 0.1, alpha = 0.05) +
  geom_function(fun = sr_model_segreg, args = sr_pars_her_med, size = 0.4,
                colour = "red", n = 10000) +
  geom_point(size = 0.3, colour = "red") +
  scale_x_continuous("SSB [million t]", 
                     breaks = c(0, 500000, 1000000, 1500000, 2000000),
                     labels = c("0", "0.5", "1.0", "1.5", "2.0")) +
  scale_y_continuous("Recruitment [millions]",
                     breaks = c(0, 10000000, 20000000, 30000000, 40000000),
                     labels = c(0, 10, 20, 30, 40)) +
  coord_cartesian(xlim = c(0, 2400000), ylim = c(0, 50000000), expand = FALSE) +
  facet_wrap(~ "Herring") +
  theme_bw(base_size = 8)


p <- p_ple + p_cod + p_her + plot_layout(nrow = 1)
ggsave(filename = "output/plots/OM/OM_rec_baseline.png", plot = p,
       width = 16, height = 6, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/OM/OM_rec_baseline.pdf", plot = p,
       width = 16, height = 6, units = "cm", dpi = 600)

### ------------------------------------------------------------------------ ###
### Recruitment models of alternative OMs (all stocks) ####
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


### load recruitment values and stock-recruit pairs - use medians
df_OM <- foreach(stock_id = res_alt_OMs$stock,
                 OM = res_alt_OMs$OM,
                 OM_label2 = res_alt_OMs$OM_label2,
                 stock = res_alt_OMs$stock_label) %do% {
  #browser()
  sr_i <- readRDS(paste0("input/", stock_id, "/", 
                         ifelse(identical(OM, "rec_failure"),
                                "baseline", OM), 
                         "/1000_100/sr.rds"))
  ### SSB-recruitment pairs used for model fitting
  df_pairs <- data.frame(year = as.numeric(dimnames(ssb(sr_i))$year),
                         ssb = c(iterMedians(ssb(sr_i))), 
                         rec = c(iterMedians(rec(sr_i)))) %>%
    filter(!is.na(ssb)) %>%
    mutate(stock = stock, stock_id = stock_id,
           OM = OM, OM_label2 = OM_label2)
  ### table with modelled recruitment
  sr_i_med <- iter(sr_i, 1)
  params(sr_i_med) <- iterMedians(params(sr_i))
  ssb_max <- max(c(iterMedians(ssb(sr_i))), na.rm = TRUE)
  
  if (isTRUE(model(FLSR(model = bevholt)) == model(sr_i))) {
    ssb_vals <- seq(from = 0, to = ssb_max * 3, length.out = 1500)
    rec_vals <- c(predict(sr_i_med, ssb = FLQuant(ssb_vals)))
  } else if (isTRUE(model(FLSR(model = segreg)) == model(sr_i))) {
    ssb_vals <- c(0, c(params(sr_i_med)$b), ssb_max * 3)
    rec_vals <- c(0, c(params(sr_i_med)$a * params(sr_i_med)$b),
                  c(params(sr_i_med)$a * params(sr_i_med)$b))
  } else {
    stop("unknown model")
  }
  df_pairs_model <- data.frame(ssb = ssb_vals,
                               rec = rec_vals) %>%
    mutate(stock = stock, stock_id = stock_id,
           OM = OM, OM_label2 = OM_label2, rec_failure = FALSE)
  if (isTRUE(OM == "rec_failure")) {
    df_pairs_model <- bind_rows(
      df_pairs_model %>%
        mutate(rec_failure = FALSE),
      df_pairs_model %>%
        mutate(rec_failure = TRUE)) %>%
      mutate(rec = ifelse(rec_failure, rec/10, rec))
  }
  ### recruitment model parameters
  pars_i <- data.frame(a = c(iterMedians(params(sr_i)$a)),
                       b = c(iterMedians(params(sr_i)$b))) %>%
    mutate(stock = stock, stock_id = stock_id,
           OM = OM, OM_label2 = OM_label2)
  return(list(pairs = df_pairs, pairs_add = df_pairs_model, pars = pars_i))
}
df_pairs <- bind_rows(lapply(df_OM, "[[", 1))
df_pairs_add <- bind_rows(lapply(df_OM, "[[", 2))
df_pars <- bind_rows(lapply(df_OM, "[[", 3))


p_ple <- ggplot() +
  geom_point(data = df_pairs %>%
               filter(stock == "Plaice"),
             aes(x = ssb, y = rec),
             size = 0.3) +
  geom_line(data = df_pairs_add %>% 
              filter(stock == "Plaice"),
            aes(x = ssb, y = rec, linetype = rec_failure),
            size = 0.4, show.legend = FALSE) + 
  facet_wrap(~ OM_label2, ncol = 5) +
  scale_x_continuous("SSB [1000t]",
                     breaks = c(0, 2000, 4000, 6000, 8000),
                     labels = c(0, 2, 4, 6, 8)) +
  scale_y_continuous("Recruitment [1000s]",
                     breaks = c(0, 10000, 20000, 30000),
                     labels = c(0, 10, 20, 30)) +
  coord_cartesian(xlim = c(0, 9500), ylim = c(0, 35000), expand = FALSE) +
  labs(title = "(a) Plaice") +
  theme_bw(base_size = 8) +
  theme(plot.title = element_text(face = "bold"))

p_cod <- ggplot() +
  geom_point(data = df_pairs %>%
               filter(stock == "Cod"),
             aes(x = ssb, y = rec),
             size = 0.3) +
  geom_line(data = df_pairs_add %>% 
              filter(stock == "Cod"),
            aes(x = ssb, y = rec, linetype = rec_failure),
            size = 0.4, show.legend = FALSE) + 
  facet_wrap(~ OM_label2, ncol = 5) +
  scale_x_continuous("SSB [1000t]",
                     breaks = c(0, 25000, 50000, 75000, 100000),
                     labels = c(0, 25, 50, 75, 100)) +
  scale_y_continuous("Recruitment [millions]",
                     breaks = c(0, 500000, 1000000, 1500000),
                     labels = c(0, "0.5", "1.0", "1.5")) +
  coord_cartesian(xlim = c(0, 125000), ylim = c(0, 1500000), expand = FALSE) +
  labs(title = "(b) Cod") +
  theme_bw(base_size = 8) +
  theme(plot.title = element_text(face = "bold"))

p_her <- ggplot() +
  geom_point(data = df_pairs %>%
               filter(stock == "Herring"),
             aes(x = ssb, y = rec),
             size = 0.3) +
  geom_line(data = df_pairs_add %>% 
              filter(stock == "Herring"),
            aes(x = ssb, y = rec, linetype = rec_failure),
            size = 0.4, show.legend = FALSE) + 
  facet_wrap(~ OM_label2, ncol = 5) +
  scale_x_continuous("SSB [million t]",
                     breaks = c(0, 2000000, 4000000, 6000000),
                     labels = c(0, 2, 4, 6)) +
  scale_y_continuous("Recruitment [millions]",
                     breaks = c(0, 50000000, 100000000, 150000000),
                     labels = c(0, 50, 100, 150)) +
  coord_cartesian(xlim = c(0, 6200000), ylim = c(0, 180000000), expand = FALSE) +
  labs(title = "(c) Herring") +
  theme_bw(base_size = 8) +
  theme(plot.title = element_text(face = "bold"))

p <- (p_ple + p_cod + p_her + plot_layout(heights = c(2.4, 1, 1)))
p
ggsave(filename = "output/plots/OM/OM_rec_OMs.png", plot = p,
       width = 16, height = 16, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/OM/OM_rec_OMs.pdf", plot = p,
       width = 16, height = 16, units = "cm", dpi = 600)

### ------------------------------------------------------------------------ ###
### fishery selectivity ####
### ------------------------------------------------------------------------ ###

fit_ple <- readRDS("input/ple.27.7e/baseline/1000_100/SAM_fit.rds")
fit_her <- readRDS("input/her.27.3a47d/baseline/1000_100/SAM_fit.rds")
fit_cod <- readRDS("input/cod.27.47d20/baseline/1000_100/SAM_fit.rds")

bind_rows(as.data.frame(apply(faytable(fit_ple), 1, function(x) x/max(x))) %>%
            rownames_to_column("age") %>%
            mutate(age = as.numeric(age)) %>%
            pivot_longer(-age, names_to = c("year"), 
                         names_transform = c(year = "as.numeric")) %>%
            mutate(stock = "Plaice"),
          as.data.frame(apply(faytable(fit_cod), 1, function(x) x/max(x))) %>%
            rownames_to_column("age") %>%
            mutate(age = as.numeric(age)) %>%
            pivot_longer(-age, names_to = c("year"), 
                         names_transform = c(year = "as.numeric")) %>%
            mutate(stock = "Cod"),
          as.data.frame(apply(faytable(fit_her), 1, function(x) x/max(x))) %>%
            rownames_to_column("age") %>%
            mutate(age = as.numeric(age)) %>%
            pivot_longer(-age, names_to = c("year"), 
                         names_transform = c(year = "as.numeric")) %>%
            mutate(stock = "Herring")
          ) %>%
  mutate(stock = factor(stock, levels = c("Plaice", "Cod", "Herring"))) %>%
  filter(year >= 2010) %>%
  ggplot(aes(x = age, y = value, colour = as.factor(year))) +
  geom_line() +
  scale_color_discrete("year") +
  scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10)) +
  facet_wrap(~ stock) + 
  labs(y = "Fishery selectivity") +
  theme_bw(base_size = 8)
  
  
  
  
as.data.frame(apply(faytable(fit_ple), 1, function(x) x/max(x))) %>%
  rownames_to_column("age") %>%
  mutate(age = as.numeric(age)) %>%
  pivot_longer(-age, names_to = c("year"), 
               names_transform = c(year = "as.numeric")) %>%
  filter(year >= 2010) %>%
  ggplot(aes(x = age, y = value, colour = as.factor(year))) +
  geom_line() +
  scale_color_discrete("year") +
  labs(y = "Fishery selectivity") +
  theme_bw(base_size = 8)
as.data.frame(apply(faytable(fit_her), 1, function(x) x/max(x))) %>%
  rownames_to_column("age") %>%
  mutate(age = as.numeric(age)) %>%
  pivot_longer(-age, names_to = c("year"), 
               names_transform = c(year = "as.numeric")) %>%
  filter(year >= 2010) %>%
  ggplot(aes(x = age, y = value, colour = as.factor(year))) +
  geom_line() +
  scale_color_discrete("year") +
  theme_bw(base_size = 8)
as.data.frame(apply(faytable(fit_cod), 1, function(x) x/max(x))) %>%
  rownames_to_column("age") %>%
  mutate(age = as.numeric(age)) %>%
  pivot_longer(-age, names_to = c("year"), 
               names_transform = c(year = "as.numeric")) %>%
  filter(year >= 2010) %>%
  ggplot(aes(x = age, y = value, colour = as.factor(year))) +
  geom_line() +
  scale_color_discrete("year") +
  theme_bw(base_size = 8)
