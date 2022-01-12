### ------------------------------------------------------------------------ ###
### create OM for North Sea herring her.27.3a47d ####
### ------------------------------------------------------------------------ ###
### base OM on SAM model fit
### follow OM routines developed for ICES WKNSMSE 2018
### based on 2021 ICES assessment
### simplifications:
### - use SAM R package instead of FLSAM
### - exclude the four partial LAI indices
###   (same as during WKNSMSE, surveys have negligible impact,
###    but require different SAM version and substantially slow down SAM)


library(ggplot2)
library(FLCore)
library(FLAssess)
library(FLXSA)
library(FLash)
library(FLfse)
library(ggplotFL)
library(stockassessment)
library(foreach)
library(dplyr)
library(tidyr)
library(doParallel)
library(mse)

source("funs.R")
source("funs_WKNSMSE.R")
source("funs_OM.R")

### load input data
stk <- readRDS("input/her.27.3a47d/preparation/stk.rds")
idx <- readRDS("input/her.27.3a47d/preparation/idx.rds")
ctrl <- readRDS("input/her.27.3a47d/preparation/conf.rds")

### index
### use HERAS for biomass index
### stock.wt is taken from the HERAS survey
### but uses 3-year average 
### -> use raw weights from index for HERAS biomass index
idx.wt <- readVPAFile("input/her.27.3a47d/preparation/west_raw.txt")
idx.wt <- idx.wt[, ac(1989:2020)]
### includes age 9, need to create plusgroup 8 (last age in stock assessment)
idx.n <- FLQuant(NA, dimnames = list(age = 1:9, year = 1989:2020))
idx.n[ac(1:8)] <- catch.n(idx$HERAS)
idx.n[ac(9)] <- index(idx$HERAS)[ac(8)] - catch.n(idx$HERAS)[ac(8)]
### weighting of index number in age 8 and 9+
idx.n_wt <- idx.n[ac(8:9)]/rep(c(apply(idx.n[ac(8:9)], 2, sum)), each = 2)
### calculate mean weight, weighted by numbers
idx.wt[ac(8)] <- apply(idx.wt[ac(8:9)] * idx.n_wt, 2, sum)
idx.wt <- idx.wt[ac(1:8)]
### insert into index
idx$HERAS@catch.wt <- idx.wt

### age length key
ALK_MSE <- readRDS("input/her.27.3a47d/preparation/ALK_MSE.rds")

refpts <- list(
  ### ICES style EqSim reference points (run with SAM fit)
  EqSim_Btrigger = 1232828, EqSim_Fmsy = 0.31, EqSim_Fpa = 0.31, 
  EqSim_Bpa = 956483, EqSim_Blim = 874198,
  ### real OM MSY values
  Fmsy = NA, Bmsy = NA, Cmsy = NA, Blim = NA,
  ### length reference points
  Lc = 25, Lref = 0.75*25 + 0.25*31
)


### ------------------------------------------------------------------------ ###
### create OMs ####
### ------------------------------------------------------------------------ ###

### default OM
#debugonce(create_OM)
create_OM(stk_data = stk, idx_data = idx, n = 1000, n_years = 100,
          yr_data = 2020, int_yr = TRUE,
          SAM_conf = ctrl, SAM_newtonsteps = 0, SAM_rel.tol = 0.001,
          SAM_NA_rm = FALSE,
          n_sample_yrs = 10, sr_model = "segreg", sr_start = 2002,
          sr_parallel = 10, sr_ar_check = TRUE, 
          process_error = TRUE, catch_oem_error = TRUE,
          idx_weights = c("index.wt", "none", "none", "none"), idxB = "HERAS", 
          idxL = TRUE, ALKs = ALK_MSE,
          ALK_yrs = 2016:2020, length_samples = 2000, PA_status = TRUE,
          refpts = refpts, stock_id = "her.27.3a47d", OM = "baseline", 
          save = TRUE, return = FALSE, M_alternative = NULL)


### higher recruitment - from 1981 instead of 2002
create_OM(stk_data = stk, idx_data = idx, n = 1000, n_years = 100,
          yr_data = 2020, int_yr = TRUE,
          SAM_conf = ctrl, SAM_newtonsteps = 0, SAM_rel.tol = 0.001,
          SAM_NA_rm = FALSE,
          n_sample_yrs = 10, sr_model = "segreg", sr_start = 1981,
          sr_parallel = 10, sr_ar_check = TRUE, 
          process_error = TRUE, catch_oem_error = TRUE,
          idx_weights = c("index.wt", "none", "none", "none"), idxB = "HERAS", 
          idxL = TRUE, ALKs = ALK_MSE,
          ALK_yrs = 2016:2020, length_samples = 2000, PA_status = TRUE,
          refpts = refpts, stock_id = "her.27.3a47d", OM = "rec_higher", 
          save = TRUE, return = FALSE, M_alternative = NULL)



### ------------------------------------------------------------------------ ###
### MSY reference points ####
### ------------------------------------------------------------------------ ###
### called from OM_MSY.pbs -> OM_MSY.R

### ------------------------------------------------------------------------ ###
### update MSY reference points for alternative OMs ####
### ------------------------------------------------------------------------ ###

stk_baseline <- readRDS("input/her.27.3a47d/baseline/1000_100/stk.rds")
Blim <- 874198 # from ICES Advice Sheet 2021
sr_baseline <- readRDS("input/her.27.3a47d/baseline/1000_100/sr.rds")
### find ratio of R(SSB=Blim)/R0 -> definition of Blim
Blim/iterMedians(params(sr_baseline))["b"]
### alternative: use position of Blim relative to hockey-stick breakpoint
Blim_ratio <- Blim / c(iterMedians(params(sr_baseline))["b"])
### 0.4701833
### Blim is breakpoint * 0.4701833


refpts <- FLPar(refpts, iter = 1000, unit = "")
update_refpts <- function(stock_id = "her.27.3a47d", OM, refpts, 
                          Blim_ratio = FALSE) {
  ### get MSY levels 
  refpts_MSY <- readRDS(paste0("input/", stock_id, "/", OM,
                               "/1000_100/MSY_trace.rds"))
  refpts_MSY <- refpts_MSY[[which.max(sapply(refpts_MSY, function(x) x$catch))]]
  ### update
  refpts["Fmsy"] <- refpts_MSY$Ftrgt
  refpts["Bmsy"] <- refpts_MSY$ssb
  refpts["Cmsy"] <- refpts_MSY$catch
  ### load recruitment model and estimate Blim
  sr_mse <- readRDS(paste0("input/", stock_id, "/", OM, "/1000_100/sr.rds"))
  pars <- iterMedians(params(sr_mse))
  if (!isFALSE(Blim_ratio))
    refpts["Blim"] <- c(pars["b"])*Blim_ratio
  print(refpts)
  ### save updated values
  saveRDS(refpts, file = paste0("input/", stock_id, "/", OM,
                                "/1000_100/refpts_mse.rds"))
}

### baseline
update_refpts(OM = "baseline", refpts = refpts, Blim_ratio = Blim_ratio)
### higher recruitment
update_refpts(OM = "rec_higher", refpts = refpts, Blim_ratio = Blim_ratio)



### ------------------------------------------------------------------------ ###
### for harvest rate: check mean catch length history ####
### ------------------------------------------------------------------------ ###

### load stock
stk <- readRDS("input/cod.27.47d20/preparation/stk.rds")

### load full ALK history from IBTS
ALKs <- readRDS("input/cod.27.47d20/preparation/ALK_MSE.rds")

### indices 
### use observed values - equivalent to simulated plus added uncertainty
idx <- readRDS("input/cod.27.47d20/preparation/idx.rds")
idx$IBTS_Q3_gam@index ### 1992-2020
idx_weights_Q3 <- readRDS("input/cod.27.47d20/preparation/idx_weights_Q3.rds")

### aggregated biomass index
idxB <- quantSums(idx$IBTS_Q3_gam@index * idx_weights_Q3)
plot(idxB) + ylim(c(0, NA))
### corresponding catch
idxC <- catch(stk)[, ac(1992:2020)]

### harvest rate
plot(idxC/idxB) + ylim(c(0, NA))

### calculate mean catch length
Lc <- 20
LFeM <- 0.75*20 + 0.25*113
lmean <- left_join(
  ### observed catch numbers at age
  x = as.data.frame(catch.n(stk)[, ac(1992:2020)]) %>%
    select(year, age, data) %>%
    rename("caa" = "data") %>%
    mutate(caa = caa + .Machine$double.eps),
  ### merge with ALKs
  y = ALKs,
  by = c("year", "age")) %>%
  ### calculate numbers at length
  mutate(cal = caa * freq) %>%
  ### keep only numbers where length >= Lc
  filter(length >= Lc) %>% 
  ### mean catch length above Lc
  group_by(year) %>%
  summarise(mean = weighted.mean(x = length, w = cal))
lmean_above <- lmean %>% filter(mean >= LFeM)
lmean_above$year
# [1] 2008 2009 2010 2011 2012 2013 2014 2015 2016 2017 2018 2019
### 2008-2019

### plot mean length
ggplot() +
  geom_hline(yintercept = LFeM, size = 0.4, colour = "red") +
  geom_line(data = lmean, aes(x = year, y = mean),
            size = 0.3) +
  geom_point(data = lmean_above, 
            aes(x = year, y = mean),
            size = 0.5) +
  ylim(c(0, NA)) + xlim(c(1990, 2020)) +
  labs(y = "mean catch length [cm]") +
  theme_bw(base_size = 8)
ggsave(filename = "output/plots/OM/OM_cod_mean_length.png", 
       width = 8.5, height = 5, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/OM/OM_cod_mean_length.pdf", 
       width = 8.5, height = 5, units = "cm", dpi = 600)

### plot harvest rate
df_hr <- as.data.frame(idxC/idxB)
ggplot() +
  geom_line(data = df_hr, aes(x = year, y = data),
            size = 0.3) +
  geom_point(data = df_hr %>% filter(year %in% lmean_above$year), 
             aes(x = year, y = data),
             size = 0.5, colour = "red") +
  geom_line(data = df_hr %>% 
              filter(year %in% lmean_above$year) %>%
              mutate(data = mean(data)), 
            aes(x = year, y = data),
            size = 0.5, colour = "red") +
  ylim(c(0, NA)) + xlim(c(1990, 2020)) +
  labs(y = "harvest rate (catch/index)") +
  theme_bw(base_size = 8)
ggsave(filename = "output/plots/OM/OM_cod_mean_length_hr_target.png", 
       width = 8.5, height = 5, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/OM/OM_cod_mean_length_hr_target.pdf", 
       width = 8.5, height = 5, units = "cm", dpi = 600)
