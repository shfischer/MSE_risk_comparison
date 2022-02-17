### ------------------------------------------------------------------------ ###
### figures for visualising empirical MP parameters ####
### ------------------------------------------------------------------------ ###

library(ggplot2)
library(FLCore)
library(ggplotFL)
library(foreach)
library(dplyr)
library(tidyr)
library(patchwork)

source("funs.R")
source("funs_WKNSMSE.R")
source("funs_OM.R")

### ------------------------------------------------------------------------ ###
### hr rule: target harvest rate ####
### ------------------------------------------------------------------------ ###
### first: check mean catch length ####

### plaice
ple_LFeM <- 36
ple_lmean <- read.csv("input/ple.27.7e/preparation/lmean.csv")
### always below LFeM
ple_lmean <- ple_lmean %>%
  rename(year = Year) %>%
  mutate(stock = "Plaice",
         LFeM = ple_LFeM)
p_ple_lmean <- ggplot() +
  geom_line(data = ple_lmean, 
            aes(x = year, y = Lmean), size = 0.4) +
  geom_point(data = ple_lmean, 
            aes(x = year, y = Lmean), size = 0.4) +
  geom_hline(yintercept = ple_LFeM, colour = "black", size = 0.4) +
  ylim(c(0, NA)) + xlim(c(2003, 2020)) +
  facet_wrap(~ "Plaice") + 
  labs(y = "mean catch length [cm]") +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())

### cod
cod_Lc <- 20
cod_LFeM <- 0.75*20 + 0.25*117
cod_stk <- readRDS("input/cod.27.47d20/preparation/stk.rds")
cod_ALKs <- readRDS("input/cod.27.47d20/preparation/ALK_MSE.rds")
cod_lmean <- left_join(
  ### observed catch numbers at age
  x = as.data.frame(catch.n(cod_stk)[, ac(1992:2020)]) %>%
    select(year, age, data) %>%
    rename("caa" = "data") %>%
    mutate(caa = caa + .Machine$double.eps),
  ### merge with ALKs
  y = cod_ALKs,
  by = c("year", "age")) %>%
  ### calculate numbers at length
  mutate(cal = caa * freq) %>%
  ### keep only numbers where length >= Lc
  filter(length >= cod_Lc) %>% 
  ### mean catch length above Lc
  group_by(year) %>%
  summarise(Lmean = weighted.mean(x = length, w = cal)) %>%
  mutate(stock = "Cod",
         LFeM = cod_LFeM,
         use = ifelse(Lmean >= LFeM, TRUE, FALSE))
p_cod_lmean <- ggplot() +
  geom_line(data = cod_lmean, 
            aes(x = year, y = Lmean), size = 0.4) +
  geom_point(data = cod_lmean, 
             aes(x = year, y = Lmean), size = 0.4) +
  geom_hline(yintercept = cod_LFeM, colour = "black", size = 0.4) +
  geom_point(data = cod_lmean %>% filter(use == TRUE), 
             aes(x = year, y = Lmean), size = 0.4, colour = "red") +
  ylim(c(0, NA)) + xlim(c(1992, 2020)) +
  facet_wrap(~ "Cod") + 
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

### herring
her_Lc <- 25
her_LFeM <- 0.75*25 + 0.25*31
her_stk <- readRDS("input/her.27.3a47d/preparation/stk.rds")
her_ALKs <- readRDS("input/her.27.3a47d/preparation/ALK_MSE.rds")
### combine all ALKs for years <2010 without ALK data
### equal weighting of years
her_ALKs_hist <- her_ALKs %>%
  select(-year) %>%
  group_by(age) %>%
  mutate(freq = freq/sum(freq)) %>%
  group_by(age, length) %>%
  summarise(freq = sum(freq))
### estimate mean length in catch
her_lmean <- left_join(
  ### observed catch numbers at age
  x = as.data.frame(catch.n(her_stk)[, ac(2010:2020)]) %>%
    select(year, age, data) %>%
    rename("caa" = "data") %>%
    mutate(caa = caa + .Machine$double.eps),
  ### merge with ALKs
  y = her_ALKs,
  by = c("year", "age")) %>%
  ### calculate numbers at length
  mutate(cal = caa * freq) %>%
  ### keep only numbers where length >= Lc
  filter(length >= her_Lc) %>% 
  ### mean catch length above Lc
  group_by(year) %>%
  summarise(Lmean = weighted.mean(x = length, w = cal)) %>%
  mutate(stock = "Herring",
         LFeM = her_LFeM)
### mean length in catch for years < 2010 (with pooled ALK)
her_lmean_hist <- left_join(
  ### observed catch numbers at age
  x = as.data.frame(catch.n(her_stk)[, ac(1989:2009)]) %>%
    select(year, age, data) %>%
    rename("caa" = "data") %>%
    mutate(caa = caa + .Machine$double.eps),
  ### merge with ALKs
  y = her_ALKs_hist,
  by = c("age")) %>%
  ### calculate numbers at length
  mutate(cal = caa * freq) %>%
  ### keep only numbers where length >= Lc
  filter(length >= her_Lc) %>% 
  ### mean catch length above Lc
  group_by(year) %>%
  summarise(Lmean = weighted.mean(x = length, w = cal)) %>%
  mutate(stock = "Herring",
         LFeM = her_LFeM)
### combine both
her_lmean <- bind_rows(her_lmean_hist, her_lmean) %>%
  mutate(use = ifelse(Lmean >= LFeM, TRUE, FALSE))
### plot
p_her_lmean <- ggplot() +
  geom_line(data = her_lmean, 
            aes(x = year, y = Lmean), size = 0.4) +
  geom_point(data = her_lmean, 
             aes(x = year, y = Lmean), size = 0.4) +
  geom_point(data = her_lmean %>% filter(use == TRUE), 
             aes(x = year, y = Lmean), size = 0.4, colour = "red") +
  geom_hline(yintercept = her_LFeM, colour = "black", size = 0.4) +
  ylim(c(0, NA)) + xlim(c(1989, 2020)) +
  facet_wrap(~ "Herring") + 
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

### second: harvest rate target 

### plaice
ple_stk <- readRDS("input/ple.27.7e/preparation/model_input_stk_d.RDS")
### indices 
### use observed values - equivalent to simulated plus added uncertainty
ple_idx <- readRDS("input/ple.27.7e/preparation/model_input_idx.RDS")
ple_idx$`FSP-7e`@index ### 2003-2020
### aggregated biomass index
ple_idxB <- quantSums(ple_idx$`FSP-7e`@index * 
                        catch.wt(ple_stk)[ac(2:8), ac(2003:2020)])
### corresponding catch
ple_idxC <- catch(ple_stk)[, ac(2003:2020)]
### harvest rate
ple_hr <- as.data.frame(ple_idxC/ple_idxB) %>%
  select(year, hr = data) %>%
  mutate(stock = "Plaice") %>%
  mutate(use = ifelse(year == 2014, TRUE, FALSE))

p_ple_hr <- ggplot() +
  geom_line(data = ple_hr, 
            aes(x = year, y = hr), size = 0.4) +
  geom_point(data = ple_hr, 
             aes(x = year, y = hr), size = 0.4) +
  geom_point(data = ple_hr %>% filter(use), 
             aes(x = year, y = hr), size = 0.4, colour = "red") +
  geom_hline(yintercept = mean(ple_hr$hr[ple_hr$use]),
             linetype = "dashed", colour = "red", size = 0.4) +
  ylim(c(0, NA)) + xlim(c(2003, 2020)) +
  labs(y = "Harvest rate (catch/index)", x = "") +
  theme_bw(base_size = 8) +
  theme()



### cod
cod_stk <- readRDS("input/cod.27.47d20/preparation/stk.rds")
cod_idx <- readRDS("input/cod.27.47d20/preparation/idx.rds")
cod_idx$IBTS_Q3_gam@index ### 1992-2020
cod_idx_weights_Q3 <- readRDS("input/cod.27.47d20/preparation/idx_weights_Q3.rds")
cod_idxB <- quantSums(cod_idx$IBTS_Q3_gam@index * cod_idx_weights_Q3)
plot(cod_idxB) + ylim(c(0, NA))
cod_idxC <- catch(cod_stk)[, ac(1992:2020)]
### harvest rate
cod_hr <- as.data.frame(cod_idxC/cod_idxB) %>%
  select(year, hr = data) %>%
  mutate(stock = "Cod") %>%
  mutate(use = ifelse(year %in% cod_lmean$year[cod_lmean$use], TRUE, FALSE))
p_cod_hr <- ggplot() +
  geom_line(data = cod_hr, 
            aes(x = year, y = hr), size = 0.4) +
  geom_point(data = cod_hr, 
             aes(x = year, y = hr), size = 0.4) +
  geom_point(data = cod_hr %>% filter(use), 
             aes(x = year, y = hr), size = 0.4, colour = "red") +
  geom_hline(yintercept = mean(cod_hr$hr[cod_hr$use]),
             linetype = "dashed", colour = "red", size = 0.4) +
  ylim(c(0, NA)) + xlim(c(1992, 2020)) +
  theme_bw(base_size = 8) +
  theme(axis.title.y = element_blank())



### herring
her_stk <- readRDS("input/her.27.3a47d/preparation/stk.rds")
her_idx <- readRDS("input/her.27.3a47d/preparation/idx_wts.rds")
her_idx$HERAS@index ### 1989-2020
her_idxB <- quantSums(her_idx$HERAS@index * her_idx$HERAS@catch.wt)
plot(her_idxB) + ylim(c(0, NA))
her_idxC <- catch(her_stk)[, ac(1989:2020)]
### harvest rate
her_hr <- as.data.frame(her_idxC/her_idxB) %>%
  select(year, hr = data) %>%
  mutate(stock = "Herring") %>%
  mutate(use = ifelse(year %in% her_lmean$year[her_lmean$use], TRUE, FALSE))

p_her_hr <- ggplot() +
  geom_line(data = her_hr, 
            aes(x = year, y = hr), size = 0.4) +
  geom_point(data = her_hr, 
             aes(x = year, y = hr), size = 0.4) +
  geom_point(data = her_hr %>% filter(use), 
             aes(x = year, y = hr), size = 0.4, colour = "red") +
  geom_hline(yintercept = mean(her_hr$hr[her_hr$use]),
             linetype = "dashed", colour = "red", size = 0.4) +
  ylim(c(0, NA)) + xlim(c(1989, 2020)) +
  labs(x = "") +
  theme_bw(base_size = 8) +
  theme(axis.title.y = element_blank())

### plot harvest rate target ####
p <- p_ple_lmean + p_cod_lmean + p_her_lmean +
  p_ple_hr + p_cod_hr + p_her_hr + plot_layout(ncol = 3, widths = 1)
p
ggsave(filename = "output/plots/risk_MP_hr_target.png", plot = p, 
       width = 16, height = 10, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/risk_MP_hr_target.pdf", plot = p, 
       width = 16, height = 10, units = "cm")

### ------------------------------------------------------------------------ ###
### biomass indices & safeguard ####
### ------------------------------------------------------------------------ ###

### plaice
### plaice
ple_stk <- readRDS("input/ple.27.7e/preparation/model_input_stk_d.RDS")
ple_idx <- readRDS("input/ple.27.7e/preparation/model_input_idx.RDS")
ple_idxB <- quantSums(ple_idx$`FSP-7e`@index * 
                        catch.wt(ple_stk)[ac(2:8), ac(2003:2020)])

### cod
cod_stk <- readRDS("input/cod.27.47d20/preparation/stk.rds")
cod_idx <- readRDS("input/cod.27.47d20/preparation/idx.rds")
cod_idx_weights_Q3 <- readRDS("input/cod.27.47d20/preparation/idx_weights_Q3.rds")
cod_idxB <- quantSums(cod_idx$IBTS_Q3_gam@index * cod_idx_weights_Q3)

### herring
her_stk <- readRDS("input/her.27.3a47d/preparation/stk.rds")
her_idx <- readRDS("input/her.27.3a47d/preparation/idx_wts.rds")
her_idxB <- quantSums(her_idx$HERAS@index * her_idx$HERAS@catch.wt)

### combine all stocks
idxB_all <- FLQuants(Plaice = ple_idxB, Cod = cod_idxB, Herring = her_idxB) %>%
  as.data.frame() %>%
  select(year, index = data, stock = qname) %>%
  mutate(stock = factor(stock, 
                        levels = c("Plaice", "Cod", "Herring"),
                        labels = c("Plaice - UK-FSP Q3",
                                   "Cod - IBTS Q3",
                                   "Herring - HERAS"))) %>%
  group_by(stock) %>%
  mutate(lowest = ifelse(index == min(index, na.rm = TRUE), TRUE, FALSE))
### reference levels
idxB_all_refs <- idxB_all %>%
  filter(lowest) %>%
  select(index, stock) %>%
  mutate(Iloss = index, Itrigger = index * 1.4) %>%
  pivot_longer(c(Iloss, Itrigger), names_to = "level") %>%
  mutate(level = factor(level, levels = c("Itrigger", "Iloss"),
                        labels = c("I[trigger]", "I[loss]")))
  

p <- idxB_all %>%
  ggplot(aes(x = year, y = index)) +
  geom_line(size = 0.4) +
  geom_point(data = idxB_all %>% filter(lowest), 
             aes(x = year, y = index), size = 0.4, colour = "red") +
  geom_hline(data = idxB_all_refs,
             aes(yintercept = value, linetype = level), 
             colour = "red", size = 0.4) +
  facet_wrap(~ stock, scales = "free") +
  scale_linetype_manual("", values = c("I[trigger]" = "dashed", 
                                       "I[loss]" = "solid"),
                        labels = scales::parse_format()) +
  ylim(c(0, NA)) +
  labs(x = "year", y = "biomass index") +
  theme_bw(base_size = 8) +
  theme(legend.key.height = unit(0.6, "lines"),
        legend.key.width = unit(1, "lines"))
p
ggsave(filename = "output/plots/risk_MP_idxB.png", plot = p, 
       width = 16, height = 6, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/risk_MP_idxB.pdf", plot = p, 
       width = 16, height = 6, units = "cm")

