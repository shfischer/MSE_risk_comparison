### ------------------------------------------------------------------------ ###
### create OM for North Sea cod cod.27.47d20 ####
### ------------------------------------------------------------------------ ###
### base OM on SAM model fit
### follow OM routines developed for ICES WKNSMSE 2018

# remotes::install_github("fishfollower/SAM/stockassessment@bioparprocess", 
#                         INSTALL_opts = "--no-multiarch")
# remotes::install_github("fishfollower/SAM/stockassessment", 
#                         INSTALL_opts = "--no-multiarch")


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
stk <- readRDS("input/cod.27.47d20/preparation/stk.rds")
idx <- readRDS("input/cod.27.47d20/preparation/idx.rds")
ctrl <- readRDS("input/cod.27.47d20/preparation/ctrl.rds")
fit <- readRDS("input/cod.27.47d20/preparation/fit.rds")

### age length key
ALK_MSE <- readRDS("input/cod.27.47d20/preparation/ALK_MSE.rds")

### survey weights
idx_Q1_wts <- readRDS("input/cod.27.47d20/preparation/idx_weights_Q1.rds")
idx_Q3_wts <- readRDS("input/cod.27.47d20/preparation/idx_weights_Q3.rds")
### add to survey
catch.wt(idx$IBTS_Q1_gam) <- idx_Q1_wts
catch.wt(idx$IBTS_Q3_gam) <- idx_Q3_wts

### density dependent M
dd_M_relation <- read.csv("input/cod.27.47d20/preparation/predation.csv")

refpts <- list(
  ### ICES style EqSim reference points (run with SAM fit)
  EqSim_Btrigger = 97777, EqSim_Fmsy = 0.28, EqSim_Fpa = 0.49, 
  EqSim_Bpa = 97777, EqSim_Blim = 69841,
  ### real OM MSY values
  Fmsy = NA, Bmsy = NA, Cmsy = NA, Blim = NA,
  ### length reference points
  Lc = 20, Lref = 0.75*20 + 0.25*113
)


### include maturity estimates from model 
stk_input <- stk
mat(stk_input) <- SAM2FLStock(fit, mat_est = TRUE)@mat

### simplify model configuration
ctrl_input <- ctrl
ctrl_input$matureModel <- 0
ctrl_input$keyMatureMean[] <- rep(NA, 6)
ctrl_input$keyCorObs[] <- NA

### check configuration
if (FALSE) {
  ctrl_tmp <- ctrl
  ctrl_tmp$matureModel <- 0
  ctrl_tmp$keyMatureMean[] <- rep(NA, 6)
  ctrl_tmp$keyCorObs[c(2, 4), ] <- NA
  ctrl_tmp$keyCorObs[c(3), ] <- c(0, 1, 1, -1, -1)
  fit_tmp <- FLR_SAM(stk = stk_input, idx = idx, conf = ctrl_tmp, 
                     conf_full = TRUE,
                     idx_weight = "index.var")
  set.seed(1)
  unc_tmp <- SAM_uncertainty(fit = fit_tmp, n = 10)
}
# fit_tmp <- FLR_SAM(stk = stk_input, idx = idx, conf = ctrl_input, conf_full = TRUE,
#                    idx_weight = "index.var")
# 
# set.seed(1)
# unc_light <- SAM_uncertainty2(fit = fit_tmp, n = 1000)
# ### all on
# set.seed(1)
# unc_orig <- SAM_uncertainty2(fit = fit, n = 1000)
# set.seed(1)
# unc_orig_1 <- SAM_uncertainty2(fit = fit, n = 1)
# 
# plot(c(original = fit, simplified = fit_tmp))
# fbarplot(fit)
# fbarplot(fit_tmp)
# fbarplot(c(original = fit, simplified = fit_tmp))

### changes since WKNSMSE
### survey ages for Q1 & Q3 correlated (excluding Q3 recruitment index)
### -> year effects included by default

### errors from survey covariance
### errors from catch numbers

### ------------------------------------------------------------------------ ###
### create OMs ####
### ------------------------------------------------------------------------ ###

### default OM
#debugonce(create_OM)
create_OM(stk_data = stk_input, idx_data = idx, n = 1000, n_years = 100,
          yr_data = 2020, int_yr = TRUE,
          SAM_conf = ctrl_input, SAM_conf_full = TRUE, 
          SAM_idx_weight = "index.var",
          SAM_newtonsteps = 0, SAM_rel.tol = 0.001,
          n_sample_yrs = 5, 
          sr_model = "segreg", sr_start = 1998, sr_parallel = 10,
          sr_ar_check = TRUE, process_error = TRUE, catch_oem_error = TRUE,
          idx_weights = c("index.wt", "index.wt", "none"), idxB = "IBTS_Q1_gam", 
          idxL = TRUE, ALKs = ALK_MSE,
          ALK_yrs = 2016:2020, length_samples = 2000, PA_status = TRUE,
          refpts = refpts, stock_id = "cod.27.47d20", OM = "baseline", 
          save = TRUE, return = FALSE, M_alternative = NULL)


### higher recruitment - from 1988 instead of 1998
create_OM(stk_data = stk_input, idx_data = idx, n = 1000, n_years = 100,
          yr_data = 2020, int_yr = TRUE,
          SAM_conf = ctrl_input, SAM_conf_full = TRUE, 
          SAM_idx_weight = "index.var",
          SAM_newtonsteps = 0, SAM_rel.tol = 0.001,
          n_sample_yrs = 5, 
          sr_model = "segreg", sr_start = 1988, sr_parallel = 10,
          sr_ar_check = TRUE, process_error = TRUE, catch_oem_error = TRUE,
          idx_weights = c("index.wt", "index.wt", "none"), idxB = "IBTS_Q1_gam", 
          idxL = TRUE, ALKs = ALK_MSE,
          ALK_yrs = 2016:2020, length_samples = 2000, PA_status = TRUE,
          refpts = refpts, stock_id = "cod.27.47d20", OM = "rec_higher", 
          save = TRUE,
          return = FALSE, M_alternative = NULL)

### density dependent M
create_OM(stk_data = stk_input, idx_data = idx, n = 1000, n_years = 100,
          yr_data = 2020, int_yr = TRUE,
          SAM_conf = ctrl_input, SAM_conf_full = TRUE, 
          SAM_idx_weight = "index.var",
          SAM_newtonsteps = 0, SAM_rel.tol = 0.001,
          n_sample_yrs = 5, 
          sr_model = "segreg", sr_start = 1998, sr_parallel = 10,
          sr_ar_check = TRUE, process_error = TRUE, catch_oem_error = TRUE,
          idx_weights = c("index.wt", "index.wt", "none"), idxB = "IBTS_Q1_gam", 
          idxL = TRUE, ALKs = ALK_MSE,
          ALK_yrs = 2016:2020, length_samples = 2000, PA_status = TRUE,
          refpts = refpts, stock_id = "cod.27.47d20", OM = "M_dd", 
          save = TRUE, return = FALSE, M_alternative = NULL, 
          M_dd = TRUE, M_dd_relation = dd_M_relation, M_dd_yr = 2020,
          M_dd_migration = TRUE)

### M alternative: no migration correction for age 3+
m_alt <- m(stk_input)
m_alt[ac(3:6), ac(2011:2021)] <- m_alt[ac(3:6), ac(2011:2021)] + log(1 - 0.15)
#debugonce(create_OM)
create_OM(stk_data = stk_input, idx_data = idx, n = 1000, n_years = 100,
          yr_data = 2020, int_yr = TRUE,
          SAM_conf = ctrl_input, SAM_conf_full = TRUE, 
          SAM_idx_weight = "index.var",
          SAM_newtonsteps = 0, SAM_rel.tol = 0.001,
          n_sample_yrs = 5, 
          sr_model = "segreg", sr_start = 1998, sr_parallel = 10,
          sr_ar_check = TRUE, process_error = TRUE, catch_oem_error = TRUE,
          idx_weights = c("index.wt", "index.wt", "none"), idxB = "IBTS_Q1_gam", 
          idxL = TRUE, ALKs = ALK_MSE,
          ALK_yrs = 2016:2020, length_samples = 2000, PA_status = TRUE,
          refpts = refpts, stock_id = "cod.27.47d20", OM = "M_no_migration", 
          save = TRUE, return = FALSE, M_alternative = m_alt)


# debugonce(input_mp)
# tmp <- input_mp(stock_id = "cod.27.47d20", OM = "M_dd", n_iter = 10, 
#                 MP = "ICES_SAM")
# debugonce(tmp$oem@method)
# tmp$args$fy <- 2025
# res <- do.call(mp, tmp)

### ------------------------------------------------------------------------ ###
### MSY reference points ####
### ------------------------------------------------------------------------ ###
### called from OM_MSY.pbs -> OM_MSY.R

### ------------------------------------------------------------------------ ###
### update MSY reference points for alternative OMs ####
### ------------------------------------------------------------------------ ###

stk_baseline <- readRDS("input/cod.27.47d20/baseline/1000_100/stk.rds")
Blim <- 69841 # from ICES Advice Sheet 2021
sr_baseline <- readRDS("input/cod.27.47d20/baseline/1000_100/sr.rds")
### find ratio of R(SSB=Blim)/R0 -> definition of Blim
iterMedians(params(sr_baseline))["b"]
### R/R0 approach does not work because Blim is on plateau of hockey-stick model
### alternative: use position of Blim relative to hockey-stick breakpoint
Blim_ratio <- Blim / c(iterMedians(params(sr_baseline))["b"])
### 1.11329
### Blim is breakpoint * 1.11329


refpts <- FLPar(refpts, iter = 1000, unit = "")
update_refpts <- function(stock_id = "cod.27.47d20", OM, refpts, 
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
### density dependent M
update_refpts(OM = "M_dd", refpts = refpts, Blim_ratio = Blim_ratio)
### M alternative: no migration correction for age 3+
update_refpts(OM = "M_no_migration", refpts = refpts, Blim_ratio = Blim_ratio)











