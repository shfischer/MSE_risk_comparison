### ------------------------------------------------------------------------ ###
### create OM for western English Channel plaice ple.27.7e ####
### ------------------------------------------------------------------------ ###
### base OM on SAM model fit
### follow OM routines developed for ICES WKNSMSE 2018

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
library(patchwork)

source("funs.R")
source("funs_WKNSMSE.R")
source("funs_OM.R")

### input data, including discard estimates
stk_data <- readRDS("input/ple.27.7e/preparation/model_input_stk_d.RDS")
idx_data <- readRDS("input/ple.27.7e/preparation/model_input_idx.RDS")
### use configuration similar to accepted XSA assessment
SAM_conf <- list(keyLogFpar = 
                   matrix(data = c(rep(-1, 9),
                                   0:5, 5, -1, -1,
                                   6:11, 11, 11, -1),
                          ncol = 9, nrow = 3, byrow = TRUE))
ALKs <- readRDS("input/ple.27.7e/preparation/ALK_MSE.rds")
refpts <- list(
  ### ICES style EqSim reference points (run with SAM fit)
  EqSim_Btrigger = 2954, EqSim_Fmsy = 0.241, EqSim_Fpa = 0.392, 
  EqSim_Bpa = 2954, EqSim_Blim = 2110,
  ### ICES reference points from WGCSE (run with XSA)
  ICES_Btrigger = 2443, ICES_Fmsy = 0.238,
  ### real OM MSY values
  Fmsy = 0.164, Bmsy = 9536, Cmsy = 1703, Blim = 2110,
  ### length reference points
  Lc = 26, Lref = 0.75*26 + 0.25*66
)
# round(min(ssbtable(fit)[, "Estimate"]))

### default OM
create_OM(stk_data = stk_data, idx_data = idx_data, n = 1000, n_years = 100,
          yr_data = 2020, 
          SAM_conf = SAM_conf, SAM_newtonsteps = 0, SAM_rel.tol = 0.001,
          n_sample_yrs = 5, sr_model = "bevholtSV", sr_parallel = 10,
          sr_ar_check = TRUE, process_error = TRUE, catch_oem_error = TRUE,
          idx_weights = c("catch.wt", "stock.wt"), idxB = "FSP-7e", 
          idxL = TRUE, ALKs = ALKs, ALK_yrs = 2016:2020, length_samples = 2000,
          PA_status = TRUE,
          refpts = refpts, stock_id = "ple.27.7e", OM = "baseline", save = TRUE,
          return = FALSE, M_alternative = NULL)
### high M: +50%
create_OM(stk_data = stk_data, idx_data = idx_data, n = 1000, n_years = 100,
          yr_data = 2020, 
          SAM_conf = SAM_conf, SAM_newtonsteps = 0, SAM_rel.tol = 0.001,
          n_sample_yrs = 5, sr_model = "bevholt", sr_parallel = 10,
          sr_ar_check = TRUE, process_error = TRUE, catch_oem_error = TRUE,
          idx_weights = c("catch.wt", "stock.wt"), idxB = "FSP-7e", 
          idxL = TRUE, ALKs = ALKs, ALK_yrs = 2016:2020, length_samples = 2000,
          PA_status = TRUE,
          refpts = refpts, stock_id = "ple.27.7e", OM = "M_high", save = TRUE,
          return = FALSE, M_alternative = 0.18)
### low M: -50%
create_OM(stk_data = stk_data, idx_data = idx_data, n = 1000, n_years = 100,
          yr_data = 2020, 
          SAM_conf = SAM_conf, SAM_newtonsteps = 0, SAM_rel.tol = 0.001,
          n_sample_yrs = 5, sr_model = "bevholt", sr_parallel = 10,
          sr_ar_check = TRUE, process_error = TRUE, catch_oem_error = TRUE,
          idx_weights = c("catch.wt", "stock.wt"), idxB = "FSP-7e", 
          idxL = TRUE, ALKs = ALKs, ALK_yrs = 2016:2020, length_samples = 2000,
          PA_status = TRUE,
          refpts = refpts, stock_id = "ple.27.7e", OM = "M_low", save = TRUE,
          return = FALSE, M_alternative = 0.06)
### Gislason age-dependent M 
Linf = 66 ### from WGCSE 2021
k = 0.1   ### from WGCSE 2021
t0 = -2   ### from WGCSE 2021
ages <- (2:10) + 0.5
lengths <- Linf * (1 - exp(-k * (ages - t0)))
M_age <- exp(0.55 - 1.61*log(lengths) + 1.44*log(Linf) + log(k))
create_OM(stk_data = stk_data, idx_data = idx_data, n = 1000, n_years = 100,
          yr_data = 2020, 
          SAM_conf = SAM_conf, SAM_newtonsteps = 0, SAM_rel.tol = 0.001,
          n_sample_yrs = 5, sr_model = "bevholt", sr_parallel = 10,
          sr_ar_check = TRUE, process_error = TRUE, catch_oem_error = TRUE,
          idx_weights = c("catch.wt", "stock.wt"), idxB = "FSP-7e", 
          idxL = TRUE, ALKs = ALKs, ALK_yrs = 2016:2020, length_samples = 2000,
          PA_status = TRUE,
          refpts = refpts, stock_id = "ple.27.7e", OM = "M_Gislason", 
          save = TRUE, return = FALSE, M_alternative = M_age)
### default OM but no auto-correlated recruitment
create_OM(stk_data = stk_data, idx_data = idx_data, n = 1000, n_years = 100,
          yr_data = 2020, 
          SAM_conf = SAM_conf, SAM_newtonsteps = 0, SAM_rel.tol = 0.001,
          n_sample_yrs = 5, sr_model = "bevholt", sr_parallel = 10,
          sr_ar_check = FALSE, process_error = TRUE, catch_oem_error = TRUE,
          idx_weights = c("catch.wt", "stock.wt"), idxB = "FSP-7e", 
          idxL = TRUE, ALKs = ALKs, ALK_yrs = 2016:2020, length_samples = 2000,
          PA_status = TRUE,
          refpts = refpts, stock_id = "ple.27.7e", OM = "rec_no_AC", save = TRUE,
          return = FALSE, M_alternative = NULL)

# plot(FLStocks(default = stk, M_low = stk_M_low, M_high = stk_M_high, M_Gislason = stk_M_Gislason), probs = c(0.05, 0.25, 0.5, 0.75, 0.95)) + xlim(c(NA, 2020))

### no discards - 100% discard survival
### - assume 0 discards in OM, but MP assumes discards do happen
### (discards appear in OM stk for easier conditioning of OM, but do not 
###  affect stock numbers or F)
# debugonce(create_OM)
create_OM(stk_data = stk_data, idx_data = idx_data, n = 1000, n_years = 100,
          yr_data = 2020,
          SAM_conf = SAM_conf, SAM_newtonsteps = 0, SAM_rel.tol = 0.001,
          n_sample_yrs = 5, sr_model = "bevholt", sr_parallel = 10, 
          sr_start = 1990,
          sr_ar_check = TRUE, process_error = TRUE, catch_oem_error = TRUE,
          idx_weights = c("catch.wt", "stock.wt"), idxB = "FSP-7e", 
          idxL = TRUE, ALKs = ALKs, ALK_yrs = 2016:2020, length_samples = 2000,
          PA_status = TRUE,
          refpts = refpts, stock_id = "ple.27.7e", OM = "no_discards",
          save = TRUE, return = FALSE, M_alternative = NULL,
          disc_survival = 1, disc_survival_hidden = TRUE)
### repeat but exclude dead discards in OM
### only used for MSY calculation,
create_OM(stk_data = stk_data, idx_data = idx_data, n = 1000, n_years = 100,
          yr_data = 2020,
          SAM_conf = SAM_conf, SAM_newtonsteps = 0, SAM_rel.tol = 0.001,
          n_sample_yrs = 5, sr_model = "bevholt", sr_parallel = 10,
          sr_start = 1990,
          sr_ar_check = TRUE, process_error = TRUE, catch_oem_error = TRUE,
          idx_weights = c("catch.wt", "stock.wt"), idxB = "FSP-7e", 
          idxL = TRUE, ALKs = ALKs, ALK_yrs = 2016:2020, length_samples = 2000,
          PA_status = TRUE,
          refpts = refpts, stock_id = "ple.27.7e", 
          OM = "no_discards_not_hidden",
          save = TRUE, return = FALSE, M_alternative = NULL,
          disc_survival = 1, disc_survival_hidden = FALSE)


### ------------------------------------------------------------------------ ###
### MSY reference points ####
### ------------------------------------------------------------------------ ###
### usually called from OM_MSY.pbs -> OM_MSY.R

if (FALSE) {
  ### set up parallel processing
  req_pckgs <- c("FLCore", "FLash", "FLBRP", "mse", "FLfse", "FLXSA",
                 "GA", "doParallel", "doRNG",
                 "tidyr", "dplyr", "stockassessment")
  for (i in req_pckgs) library(package = i, character.only = TRUE)
  cl <- makeCluster(10)
  registerDoParallel(cl)
  cl_length <- length(cl)
  ### load packages and functions into parallel workers
  . <- foreach(i = seq(cl_length)) %dopar% {
    for (i in req_pckgs) library(package = i, character.only = TRUE,
                                 warn.conflicts = FALSE, verbose = FALSE,
                                 quietly = TRUE)
    source("funs.R", echo = FALSE)
    source("funs_GA.R", echo = FALSE)
    source("funs_WKNSMSE.R", echo = FALSE)
  }
  
  
  ### baseline OM
  res <- est_MSY(OM = "baseline")
  res$result[which.max(res$result$catch), ]
  #        Ftrgt   catch      ssb      tsb      rec
  # 15 0.1638845 1702.92 9536.053 11136.77 6542.729
  
}

### ------------------------------------------------------------------------ ###
### update MSY reference points for alternative OMs ####
### ------------------------------------------------------------------------ ###

### find ratio of R(SSB=Blim)/R0 -> definition of Blim
stk_baseline <- readRDS("input/ple.27.7e/baseline/1000_100/stk.rds")
Blim <- min(iterMedians(ssb(stk_baseline)), na.rm = TRUE)
dimnames(stk_baseline)$year[which.min(iterMedians(ssb(stk_baseline)))]
### Blim is SSB in 2007
sr_baseline <- readRDS("input/ple.27.7e/baseline/1000_100/sr.rds")
RR0 <- c(((iterMedians(params(sr_baseline)["a"])*Blim) /
            (iterMedians(params(sr_baseline))["b"] + Blim)) /
  iterMedians(params(sr_baseline))["a"])
RR0
### Blim corresponds to SSB at ~ 77% of R0

refpts <- FLPar(refpts, iter = 1000, unit = "")
update_refpts <- function(stock_id = "ple.27.7e", OM, refpts, RR0) {
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
  refpts["Blim"] <- c(pars["b"])*(RR0/(1 - RR0))
  print(refpts)
  ### save updated values
  saveRDS(refpts, file = paste0("input/", stock_id, "/", OM,
                                "/1000_100/refpts_mse.rds"))
}

### baseline
update_refpts(OM = "baseline", refpts = refpts, RR0 = RR0)
### M low
update_refpts(OM = "M_low", refpts = refpts, RR0 = RR0)
### M high
update_refpts(OM = "M_high", refpts = refpts, RR0 = RR0)
### M Gislason
update_refpts(OM = "M_Gislason", refpts = refpts, RR0 = RR0)
### no recruitment AC
update_refpts(OM = "rec_no_ac", refpts = refpts, RR0 = RR0)
### 100% discards survival
### - use no_discards_not_hidden OM where discards are not included in OM stk
update_refpts(OM = "no_discards_not_hidden", refpts = refpts, RR0 = RR0)
file.copy(from = "input/ple.27.7e/no_discards_not_hidden/1000_100/refpts_mse.rds",
          to = "input/ple.27.7e/no_discards/1000_100/refpts_mse.rds", 
          overwrite = TRUE)


### ------------------------------------------------------------------------ ###
### for harvest rate: check mean catch length history ####
### ------------------------------------------------------------------------ ###

### load stock
stk <- readRDS("input/ple.27.7e/preparation/model_input_stk_d.RDS")

### indices 
### use observed values - equivalent to simulated plus added uncertainty
idx <- readRDS("input/ple.27.7e/preparation/model_input_idx.RDS")
idx$`FSP-7e`@index ### 2003-2020

### aggregated biomass index
idxB <- quantSums(idx$`FSP-7e`@index * catch.wt(stk)[ac(2:8), ac(2003:2020)])
plot(idxB) + ylim(c(0, NA))
### corresponding catch
idxC <- catch(stk)[, ac(2003:2020)]

### harvest rate
plot(idxC/idxB) + ylim(c(0, NA))

### get mean catch length from WGCSE 2021
Lc <- 26
LFeM <- 36
lmean <- read.csv("input/ple.27.7e/preparation/lmean.csv")
### always below LFeM

### plot mean length
ggplot() +
  geom_hline(yintercept = LFeM, size = 0.4, colour = "red") +
  geom_line(data = lmean, aes(x = Year, y = Lmean),
            size = 0.3) +
  ylim(c(0, NA)) + xlim(c(2010, 2020)) +
  labs(y = "mean catch length [cm]") +
  theme_bw(base_size = 8)
ggsave(filename = "output/plots/OM/OM_ple_mean_length.png", 
       width = 8.5, height = 5, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/OM/OM_ple_mean_length.pdf", 
       width = 8.5, height = 5, units = "cm", dpi = 600)

### plot harvest rate
df_hr <- as.data.frame(idxC/idxB)
ggplot() +
  geom_line(data = df_hr, aes(x = year, y = data),
            size = 0.3) +
  geom_point(data = df_hr %>% filter(year %in% 2014), 
             aes(x = year, y = data),
             size = 0.5, colour = "red") +
  ylim(c(0, NA)) + xlim(c(2010, 2020)) +
  labs(y = "harvest rate (catch/index)") +
  theme_bw(base_size = 8)
ggsave(filename = "output/plots/OM/OM_ple_mean_length_hr_target.png", 
       width = 8.5, height = 5, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/OM/OM_ple_mean_length_hr_target.pdf", 
       width = 8.5, height = 5, units = "cm", dpi = 600)


