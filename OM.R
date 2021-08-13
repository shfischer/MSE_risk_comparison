verbose <- TRUE
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

source("funs.R")
source("funs_WKNSMSE.R")

### ------------------------------------------------------------------------ ###
### load stock ####
### ------------------------------------------------------------------------ ###

### input data, including discard estimates
stk_data <- readRDS("input/model_input_stk_d.RDS")
idx_data <- readRDS("input/model_input_idx.RDS")
### use configuration similar to accepted XSA assessment
conf <- list(keyLogFpar = 
               matrix(data = c(rep(-1, 9),
                               0:5, 5, -1, -1,
                               6:11, 11, 11, -1),
                      ncol = 9, nrow = 3, byrow = TRUE))
fit <- FLR_SAM(stk_data, idx_data, conf = conf)

### check fitting time
system.time({
  fit_i <- FLR_SAM(stk_data, idx_data, conf = conf)
})
### ~6s
### relax model convergence
newtonsteps <- 0
rel.tol <- 0.001
system.time({
  fit_i <- FLR_SAM(stk_data, idx_data, conf = conf, 
                   newtonsteps = newtonsteps, rel.tol = rel.tol)
})
### ~2.7s
pars_ini <- getpars(fit)
system.time({
  fit_i <- FLR_SAM(stk_data, idx_data, conf = conf, 
                   newtonsteps = newtonsteps, rel.tol = rel.tol, 
                   par_ini = pars_ini)
})
### ~ 1.7s
### update initial parameters with relaxed model convergence
pars_ini <- getpars(fit_i)

### ------------------------------------------------------------------------ ###
### simulation specifications ####
### ------------------------------------------------------------------------ ###

### number of iterations/replicates
n <- 10
### number of projection years
n_years <- 20
### last data year
yr_data <- 2020

### ------------------------------------------------------------------------ ###
### create FLStock ####
### ------------------------------------------------------------------------ ###
### create template with 1 iteration

stk <- SAM2FLStock(object = fit, stk = stk_data)

if (isTRUE(verbose)) summary(stk)

if (isTRUE(verbose)) plot(stk)

### save for later comparison
stk_orig <- stk

### ------------------------------------------------------------------------ ###
### projection years ####
### ------------------------------------------------------------------------ ###
yrs_hist <- as.numeric(dimnames(stk)$year)
yrs_proj <- seq(from = dims(stk)$maxyear + 1, length.out = n_years)
yrs_mse <- sort(unique(c(yrs_hist, yrs_proj)))

### ------------------------------------------------------------------------ ###
### add uncertainty ####
### ------------------------------------------------------------------------ ###
### first approach: use variance-covariance


### add iteration dimension
stk <- FLCore::propagate(stk, n)
if (isTRUE(verbose)) dim(stk)

### add uncertainty estimated by SAM as iterations
set.seed(1)
uncertainty <- SAM_uncertainty(fit = fit, n = n)
### add noise to stock
stock.n(stk)[] <- uncertainty$stock.n
stock(stk)[] <- computeStock(stk)
### add noise to F
harvest(stk)[] <- uncertainty$harvest

### add noise to catch numbers
catch.n(stk) <- uncertainty$catch.n
catch(stk) <- computeCatch(stk)

if (isTRUE(verbose)) plot(stk, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))

### maximum observed F
if (isTRUE(verbose)) max(fbar(stk))
if (isTRUE(verbose)) max(harvest(stk))

### ------------------------------------------------------------------------ ###
### extend stock for MSE simulation ####
### ------------------------------------------------------------------------ ###

stk_stf <- stf(stk, n_years)

### ------------------------------------------------------------------------ ###
### biological data for OM ####
### ------------------------------------------------------------------------ ###
### Resample weights, maturity and natural mortality from the last 5 years 
### (2016-2020)
### set up an array with one resampled year for each projection year 
### (including intermediate year) and replicate
### use the same resampled year for all biological parameters
set.seed(2)

### use last five data years to sample biological parameters
sample_yrs <- 2016:2020
### get year position of sample years
sample_yrs_pos <- which(dimnames(stk_stf)$year %in% sample_yrs)

### create samples for biological data (weights, etc.)
### the historical biological parameters are identical for all iterations
### and consequently do not need to be treated individually
### (but keep age structure)
### create vector with resampled years
bio_samples <- sample(x = sample_yrs_pos, 
                      size = (n_years) * n, replace = TRUE)
### do the same for selectivity
sel_samples <- sample(x = sample_yrs_pos, 
                      size = (n_years) * n, replace = TRUE)
### years to be populated
bio_yrs <- which(dimnames(stk_stf)$year %in% 2021:dims(stk_stf)$maxyear)


### insert values
catch.wt(stk_stf)[, bio_yrs] <- c(catch.wt(stk)[, bio_samples,,,, 1])
stock.wt(stk_stf)[, bio_yrs] <- c(stock.wt(stk)[, bio_samples,,,, 1])
landings.wt(stk_stf)[, bio_yrs] <- c(landings.wt(stk)[, bio_samples,,,, 1])
discards.wt(stk_stf)[, bio_yrs] <- c(discards.wt(stk)[, bio_samples,,,, 1])
# m(stk_stf)[, bio_yrs] <- c(m(stk)[, bio_samples,,,, 1]) ### constant
# mat(stk_stf)[, bio_yrs] <- c(mat(stk)[, bio_samples,,,, 1]) ### constant
### use different samples for selectivity
harvest(stk_stf)[, bio_yrs] <- c(harvest(stk)[, sel_samples,,,, 1])

if (isTRUE(verbose)) plot(stk_stf)

### ------------------------------------------------------------------------ ###
### stock recruitment ####
### ------------------------------------------------------------------------ ###
### fit hockey-stick model
### get residuals from smoothed residuals

### create FLSR object
sr <- as.FLSR(stk_stf, model = "segreg")
### fit model individually to each iteration and suppress output to screen
if (n <= 10) {
  suppressWarnings(. <- capture.output(sr <- fmle(sr)))
} else {
  ### run in parallel
  cl <- makeCluster(10)
  registerDoParallel(cl)
  sr <- fmle_parallel(sr, cl)
  stopCluster(cl)
  registerDoSEQ()
  ### run again for failed iterations - not needed
  pos_error <- which(is.na(params(sr)["a"]))
  # sr_corrected <- fmle(FLCore::iter(sr, pos_error))
  # sr[,,,,, pos_error] <- sr_corrected[]
  # params(sr)[, pos_error] <- params(sr_corrected)
}

### check autocorrelation of residuals for SAM median perception
sr_med <- as.FLSR(stk_orig, model = "segreg")
suppressWarnings(. <- capture.output(sr_med <- fmle(sr_med)))
sr_acf <- acf(residuals(sr_med))
sr_rho <- sr_acf$acf[2]
if (isTRUE(verbose)) sr_rho
### lag-1 auto-correlation of recruitment residuals -> adapt residuals in MSE

### generate residuals for MSE
### years with missing residuals
yrs_res <- colnames(rec(sr))[which(is.na(iterMeans(rec(sr))))]

### try creating residuals for median
if (isTRUE(verbose)) {
  set.seed(1)
  ### get residuals
  res <- c(residuals(sr_med))
  res <- res[!is.na(res)]
  ### calculate kernel density of residuals
  density <- density(x = res)
  plot(density)
  ### sample residuals
  mu <- sample(x = res, size = length(yrs_res), replace = TRUE)
  plot(mu)
  hist(mu)
  ### "smooth", i.e. sample from density distribution
  res_new <- rnorm(n = length(yrs_res), mean = mu, sd = density$bw)
  plot(mu)
  plot(res_new)
  ### "add" autocorrelation
  res_ac <- rep(0, length(yrs_res))
  res_ac[1] <- sr_rho * tail(res, 1) + sqrt(1 - sr_rho^2) * res_new[1]
  for (r in 2:length(res_ac)) {
    res_ac[r] <- sr_rho * res_ac[r - 1] + sqrt(1 - sr_rho^2) * res_new[r]
  }
  plot(res_ac, type = "l")
  lines(res_new, col = "red")
  acf(res_ac, plot = FALSE, lag.max = 1)$acf[2]
  acf(res, plot = FALSE, lag.max = 1)$acf[2]
  acf(res_new, plot = FALSE, lag.max = 1)$acf[2]
  hist(res_ac)
  hist(res_new)
}

### go through iterations and create residuals
### use kernel density to create smooth distribution of residuals
### and sample from this distribution
res_new <- foreach(iter_i = seq(dim(sr)[6])) %do% {
  
  set.seed(iter_i)
  
  ### get residuals for current iteration
  res_i <- c(FLCore::iter(residuals(sr), iter_i))
  res_i <- res_i[!is.na(res_i)]
  
  ### calculate kernel density of residuals
  density <- density(x = res_i)
  ### sample residuals
  mu <- sample(x = res_i, size = length(yrs_res), replace = TRUE)
  ### "smooth", i.e. sample from density distribution
  res_new <- rnorm(n = length(yrs_res), mean = mu, sd = density$bw)
  
  # return(res_new)
  
  ### "add" autocorrelation
  sr_acf_i <- acf(res_i, lag.max = 1, plot = FALSE)
  sr_rho_i <- sr_acf_i$acf[2]
  res_ac <- rep(0, length(yrs_res))
  res_ac[1] <- sr_rho * tail(res_i, 1) + sqrt(1 - sr_rho^2) * res_new[1]
  for (r in 2:length(res_ac)) {
    res_ac[r] <- sr_rho * res_ac[r - 1] + sqrt(1 - sr_rho^2) * res_new[r]
  }
  return(res_ac)
  
}
if (isTRUE(verbose)) summary(exp(unlist(res_new)))
### insert into model
residuals(sr)[, yrs_res] <- unlist(res_new)
### exponentiate residuals to get factor
residuals(sr) <- exp(residuals(sr))
sr_res <- residuals(sr)

if (isTRUE(verbose)) plot(sr_res)

### ------------------------------------------------------------------------ ###
### process noise ####
### ------------------------------------------------------------------------ ###
### create FLQuant with process noise
### this will be added to the values obtained from fwd() in the MSE

### create noise for process error
set.seed(3)
proc_res <- stock.n(stk_stf) %=% 0 ### template FLQuant
proc_res[] <- stats::rnorm(n = length(proc_res), mean = 0, 
                           sd = uncertainty$proc_error)
### the proc_res values follow a normal distribution,
### exponentiate to get log-normal residuals
proc_res <- exp(proc_res)
### proc_res is a factor by which the numbers at age are multiplied

### for historical period, numbers already include process error from SAM
### -> remove deviation
proc_res[, dimnames(proc_res)$year <= 2020] <- 1

### remove deviation for first age class (recruits)
proc_res[1, ] <- 1

### try saving in stock recruitment model ... 
### this gets passed on to the projection module
fitted(sr) <- proc_res

if (isTRUE(verbose)) plot(proc_res)

### ------------------------------------------------------------------------ ###
### stf ####
### ------------------------------------------------------------------------ ###

### not done for plaice
stk_fwd <- stk_stf


### ------------------------------------------------------------------------ ###
### biological data for OEM ####
### ------------------------------------------------------------------------ ###

### base on OM stock
stk_oem <- stk_fwd

### projection years
proj_yrs <- 2021:range(stk_oem)[["maxyear"]]

### use means of sampled values for projection period
catch.wt(stk_oem)[, ac(proj_yrs)] <- 
  yearMeans(catch.wt(stk_oem)[, ac(sample_yrs)])
landings.wt(stk_oem)[, ac(proj_yrs)] <- 
  yearMeans(landings.wt(stk_oem)[, ac(sample_yrs)])
discards.wt(stk_oem)[, ac(proj_yrs)] <- 
  yearMeans(discards.wt(stk_oem)[, ac(sample_yrs)])
stock.wt(stk_oem)[, ac(proj_yrs)] <- 
  yearMeans(stock.wt(stk_oem)[, ac(sample_yrs)])
m(stk_oem)[, ac(proj_yrs)] <- yearMeans(m(stk_oem)[, ac(sample_yrs)])
mat(stk_oem)[, ac(proj_yrs)] <- yearMeans(mat(stk_oem)[, ac(sample_yrs)])

### remove stock assessment results
stock.n(stk_oem)[] <- stock(stk_oem)[] <- harvest(stk_oem)[] <- NA

### ------------------------------------------------------------------------ ###
### catch noise ####
### ------------------------------------------------------------------------ ###
### take estimates from sam: uncertainty$catch_sd is "logSdLogObs"
### assume catch observed by SAM in projection is log-normally distributed
### around operating model catch

### create noise for catch
set.seed(5)
catch_res <- catch.n(stk_fwd) %=% 0 ### template FLQuant
catch_res[] <- stats::rnorm(n = length(catch_res), mean = 0, 
                            sd = uncertainty$catch_sd)
### the catch_res values are on a normale scale,
### exponentiate to get log-normal 
catch_res <- exp(catch_res)
### catch_res is a factor by which the numbers at age are multiplied

### for historical period, pass on real observed catch
### -> remove deviation
catch_res[, dimnames(catch_res)$year <= 2020] <- 
  window(catch.n(stk_orig), end = 2020) / window(catch.n(stk_fwd), end = 2020)

if (isTRUE(verbose)) plot(catch_res, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))

### ------------------------------------------------------------------------ ###
### indices ####
### ------------------------------------------------------------------------ ###
### use real FLIndices object as template
idx <- idx_data
### use both surveys
names(idx)
### extend for simulation period
idx <- window(idx, end = yr_data + n_years)
### add iterations
idx <- lapply(idx, propagate, n)


idx_Q1_yrs <- as.numeric(dimnames(idx$Q1SWBeam)$year)
idx_FSP_yrs <- as.numeric(dimnames(idx$`FSP-7e`)$year)
idx_Q1_yrs_hist <- setdiff(idx_Q1_yrs, yrs_proj)
idx_FSP_yrs_hist <- setdiff(idx_FSP_yrs, yrs_proj)
idx_Q1_ages <- dimnames(idx$Q1SWBeam)$age
idx_FSP_ages <- dimnames(idx$`FSP-7e`)$age

### insert catchability
### set catchability for projection
for (idx_i in seq_along(idx)) {
  
  ### set catchability for projection
  index.q(idx[[idx_i]])[] <- uncertainty$survey_catchability[[idx_i]]
  
}

### index weights
### use observed weights
### Q1SWBeam - use stock weights (beginning of year)
idx$Q1SWBeam@catch.wt <- stock.wt(stk_oem)[ac(idx_Q1_ages), ac(idx_Q1_yrs)]
### FSP - use catch weights (mid-year)
idx$`FSP-7e`@catch.wt <- catch.wt(stk_oem)[ac(idx_FSP_ages), ac(idx_FSP_yrs)]

### create copy of index with original values
idx_raw <-  lapply(idx ,index)
### calculate index values
idx <- calc_survey(stk = stk_fwd, idx = idx, use_q = TRUE, use_time = TRUE, 
                   use_wt = FALSE)

### create deviances for indices
### first, get template
idx_dev <- lapply(idx, index)
### create random noise based on sd
set.seed(4)
for (idx_i in seq_along(idx_dev)) {
  ### insert sd
  idx_dev[[idx_i]][] <- uncertainty$survey_sd[[idx_i]]
  ### noise
  idx_dev[[idx_i]][] <- stats::rnorm(n = length(idx_dev[[idx_i]]),
                                     mean = 0, sd = idx_dev[[idx_i]])
  ### exponentiate to get from normal to log-normal scale
  idx_dev[[idx_i]] <- exp(idx_dev[[idx_i]])
}



### modify residuals for historical period so that index values passed to 
### stock assessment are the ones observed in reality
idx_dev$Q1SWBeam[, ac(idx_Q1_yrs_hist)] <-
  idx_raw$Q1SWBeam[, ac(idx_Q1_yrs_hist)] /
  index(idx$Q1SWBeam)[, ac(idx_Q1_yrs_hist)]
idx_dev$`FSP-7e`[, ac(idx_FSP_yrs_hist)] <-
  idx_raw$`FSP-7e`[, ac(idx_FSP_yrs_hist)] /
  index(idx$`FSP-7e`)[, ac(idx_FSP_yrs_hist)]

### check biomass index
if (isTRUE(verbose)) plot(quantSums(idx$Q1SWBeam@catch.wt * idx$Q1SWBeam@index))
### including uncertainty
if (isTRUE(verbose)) 
  plot(quantSums(idx$Q1SWBeam@catch.wt * idx$Q1SWBeam@index * idx_dev$Q1SWBeam))
if (isTRUE(verbose)) 
  plot(quantSums(idx$`FSP-7e`@catch.wt * idx$`FSP-7e`@index * idx_dev$`FSP-7e`))

### add template for biomass index
idxB <- quantSums(idx$`FSP-7e`@catch.wt * idx$`FSP-7e`@index) %=% NA_real_
idx <- FLIndices("Q1SWBeam" = idx$Q1SWBeam,
                 "FSP-7e" = idx$`FSP-7e`,
                 idxB = FLIndex(idxB))
idx_dev[length(idx_dev) + 1] <- list(idxB)
names(idx_dev)[length(idx_dev)] <- "idxB"

### ------------------------------------------------------------------------ ###
### length index ####
### ------------------------------------------------------------------------ ###

### load age-length keys
ALKs <- readRDS("input/ALK_MSE.rds")

### keep only last 5 years
ALKs <- ALKs %>%
  filter(year %in% 2016:2020)
### sort
ALKs <- ALKs %>%
  arrange(year, age, length)

### scale up
### define which ALKs are used
alk_yrs <- 2016:2020
### random samples
set.seed(89)
alk_samples <- catch(stk_fwd) %=% NA_real_
alk_samples[] <- sample(x = alk_yrs, size = length(yrs_mse) * n, replace = TRUE)
### use existing ALKs for historical years
alk_samples[, ac(alk_yrs)] <- alk_yrs

### length at first capture: from ICES WGCSE 2021
Lc <- 26

### pre-populate index
set.seed(91)
data_yr <- 1980:2020
### catch numbers at age
cn <- as.data.frame(catch.n(stk_fwd)[, ac(data_yr)]) %>%
  select(year, age, iter, data) %>%
  rename("caa" = "data") %>%
  mutate(iter = as.numeric(as.character(iter)))
### select ALK years
alks <- as.data.frame(alk_samples[, ac(data_yr)]) %>%
  select(year, iter, data) %>%
  rename("alk_year" = "data") %>%
  mutate(iter = as.numeric(as.character(iter)))
### match ALK year with data year
cn <- full_join(cn, alks)
### merge ALKs
cn <- left_join(cn, 
                 ALKs %>% rename("alk_year" = "year"))
### calculate numbers at length
cn <- cn %>%
  mutate(cal = caa * freq)
### mean catch length above Lc with sampling
length_samples <- 2000
means <- cn %>%
  filter(length >= Lc) %>%
  group_by(year, iter) %>%
  summarise(data = mean(sample(x = length, prob = cal, 
                                 size = length_samples, replace = TRUE))) %>%
  arrange(as.numeric(as.character(iter)))
### convert into index
idxL <- FLIndex(index = as.FLQuant(means))
### extend for projection
idxL <- window(idxL, end = dims(stk_fwd)$maxyear)

### add length index to index object
idx <- FLIndices(c(idx, idxL = idxL))
idx_dev[length(idx_dev) + 1] <- list(index(idxL) %=% 1)
names(idx_dev)[length(idx_dev)] <- "idxL"
idx_dev[length(idx_dev) + 1] <- list(alk_samples)
names(idx_dev)[length(idx_dev)] <- "alk_yrs"

### ------------------------------------------------------------------------ ###
### PA buffer for 2 over 3 rule ####
### ------------------------------------------------------------------------ ###
### SPiCT performance based on 
### Fischer et al. 2021 https://doi.org/10.1093/icesjms/fsab018

### index deviation
PA_status_dev <- FLQuant(NA, dimnames = list(age = c("positive", "negative"), 
                                             year = dimnames(stk_fwd)$year, 
                                             iter = dimnames(stk_fwd)$iter))
set.seed(1)
PA_status_dev["positive"] <- rbinom(n = PA_status_dev["positive"], 
                                    size = 1, prob = 0.9886215)
set.seed(2)
PA_status_dev["negative"] <- rbinom(n = PA_status_dev["negative"], 
                                    size = 1, prob = 1 - 0.4216946)

### PA status index template
PA_status_template <- FLIndex(index = ssb(stk_fwd) %=% NA_integer_)

### add to index object
idx <- FLIndices(c(idx, PA_status = PA_status_template))
idx_dev[length(idx_dev) + 1] <- list(PA_status_dev)
names(idx_dev)[length(idx_dev)] <- "PA_status"



### mimic reference points
### MSYBtrigger ~ Bmsy/2
### Basis: Bloss * 1.4
#PA_Bmsy <- median(ssb(stk_fwd)[, ac(1980)]) * 2
PA_Bmsy <- min(ssb(stk_orig)) * 1.4 * 2
### Fmsy ~ minimum Fbar (2013)
PA_Fmsy <- min(fbar(stk_orig))

### ------------------------------------------------------------------------ ###
### save ####
### ------------------------------------------------------------------------ ###

# save.image(paste0("input/image_", n, ".RData"))
# load(paste0("input/image_", n, ".RData"))

### path
input_path <- paste0("input/ple.27.7e/baseline/", n, "_", n_years, "/")
dir.create(input_path, recursive = TRUE)
### stock
saveRDS(stk_fwd, file = paste0(input_path, "stk.rds"))
# stk_fwd <- readRDS(paste0(input_path, "stk.rds"))
### stock recruitment
saveRDS(sr, file = paste0(input_path, "sr.rds"))
# sr <- readRDS(paste0(input_path, "sr.rds"))
### surveys
saveRDS(idx, file = paste0(input_path, "idx.rds"))
# idx <- readRDS(paste0(input_path, "idx.rds"))
saveRDS(idx_dev, file = paste0(input_path, "idx_dev.rds"))
# idx_dev <- readRDS(paste0(input_path, "idx_dev.rds"))
### catch noise
saveRDS(catch_res, file = paste0(input_path, "catch_res.rds"))
# catch_res <- readRDS(paste0(input_path, "catch_res.rds"))
### process error
saveRDS(proc_res, file = paste0(input_path, "proc_res.rds"))
# proc_res <- readRDS(paste0(input_path, "proc_res.rds"))
### observed stock
saveRDS(stk_oem, file = paste0(input_path, "stk_oem.rds"))
# stk_oem <- readRDS(paste0(input_path, "stk_oem.rds"))
### sam initial parameters
saveRDS(pars_ini, file = paste0(input_path, "sam_initial.rds"))
# pars_ini <- readRDS(paste0(input_path, "sam_initial.rds"))
### sam configuration
saveRDS(conf, file = paste0(input_path, "conf.rds"))
# conf <- readRDS(paste0(input_path, "conf.rds"))

#save.image(file = paste0(input_path, "image.RData"))



### ------------------------------------------------------------------------ ###
### prepare OM for MSE ####
### ------------------------------------------------------------------------ ###

### length target
Lref <- rep(0.75*26 + 0.25*66, n)
### I_trigger = 1.4 * I_loss
I_trigger = apply(quantSums(index(idx$`FSP-7e`) * idx_dev$`FSP-7e` * 
                              catch.wt(idx$`FSP-7e`)),
                  6, min, na.rm = TRUE) * 1.4

### some arguments (passed to mp())
args <- list(fy = dims(stk_fwd)$maxyear, ### final simulation year
             y0 = dims(stk_fwd)$minyear, ### first data year
             iy = yr_data, ### first simulation (intermediate) year
             nsqy = 3, ### not used, but has to provided
             nblocks = 1, ### block for parallel processing
             seed = 1 ### random number seed before starting MSE
)

### operating model
om <- FLom(stock = stk_fwd, ### stock 
           sr = sr, ### stock recruitment and precompiled residuals
           projection = mseCtrl(method = fwd_attr, 
                                args = list(maxF = 5,
                                            ### process noise on stock.n
                                            proc_res = "fitted",
                                            dupl_trgt = FALSE
                                ))
)

### observation (error) model
oem <- FLoem(method = obs_generic,
  observations = list(
    stk = stk_oem, 
    idx = idx), 
  deviances = list(
    stk = FLQuants(catch.dev = catch_res), 
    idx = idx_dev),
  args = list(use_catch_residuals = TRUE, 
              use_idx_residuals = TRUE,
              use_stk_oem = TRUE,
              use_wt = TRUE,
              PA_status = FALSE, PA_status_dev = FALSE,
              PA_Bmsy = PA_Bmsy,
              PA_Fmsy = PA_Fmsy,
              alks = as.data.frame(ALKs),
              use_age_idcs = c("Q1SWBeam", "FSP-7e"),
              biomass_index = "FSP-7e",
              length_idx = TRUE,
              Lc = 26,
              lngth_samples = 2000))
### implementation error model (banking and borrowing)
# iem <- FLiem(method = iem_WKNSMSE, 
#              args = list(BB = TRUE))
ctrl <- mpCtrl(list(
  est = mseCtrl(method = est_comps,
                args = list(comp_r = TRUE, comp_f = TRUE, comp_b = TRUE,
                            comp_c = TRUE, comp_m = 1,
                            idxB_lag = 1, idxB_range_1 = 2, idxB_range_2 = 3,
                            idxB_range_3 = 1,
                            catch_lag = 0, ### 0 to mimic advice
                            catch_range = 1,
                            Lref = Lref, 
                            I_trigger = c(I_trigger),
                            idxL_lag = 1, idxL_range = 1,
                            pa_buffer = FALSE, pa_size = 0.8, pa_duration = 3,
                            FLXSA = FALSE,
                            FLXSA.control = NULL)),
  phcr = mseCtrl(method = phcr_comps,
                 args = list(exp_r = 1, exp_f = 1, exp_b = 1)),
  hcr = mseCtrl(method = hcr_comps,
                args = list(interval = 2)),
  isys = mseCtrl(method = is_comps,
                 args = list(interval = 2, 
                             upper_constraint = 1.2, lower_constraint = 0.7, 
                             cap_below_b = FALSE))
))
### additional tracking metrics
tracking <- c("comp_c", "comp_i", "comp_r", "comp_f", "comp_b",
              "multiplier", "exp_r", "exp_f", "exp_b")

### reference points
refpts_mse <- FLPar(ICES_Btrigger = 2954, ICES_Ftrgt = 0.241, 
                    ICES_Fpa = 0.392, ICES_Bpa = 2954,
                    Blim = 2110, 
                    Fmsy = 0.18, Bmsy = 8543, Cmsy = 1666,
                    iter = seq(n))


### save mse objects
input <- list(om = om, oem = oem, ctrl = ctrl,
              args = args, tracking = tracking, refpts = refpts_mse,
              cut_hist = TRUE)

saveRDS(input, file = paste0(input_path, "input_rfb.rds"))

if (FALSE) {
  #input$args$nblocks <- 200
  #debugonce(input$oem@method)
  #debugonce(goFish)
  #debugonce(input$ctrl$isys@method)
  set.seed(1)
  res <- do.call(mp, input)
  saveRDS(res, file = paste0("output/default/res_", n, "_rfb.rds"))
  plot(res, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
  plot(res, probs = c(0.05, 0.25, 0.5, 0.75, 0.95), iter = 1:5)
  plot(res@oem@observations$idx$idxB@index, 
       probs = c(0.05, 0.25, 0.5, 0.75, 0.95)) +
    ylim(c(0, NA))
  plot(res@oem@observations$idx$idxL@index, 
       probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
}
### ------------------------------------------------------------------------ ###
### 2 over 3 rule with PA buffer ####
### ------------------------------------------------------------------------ ###

input_2over3 <- input
input_2over3$oem@args$length_idx <- FALSE
input_2over3$oem@args$PA_status <- TRUE
input_2over3$oem@args$PA_status_dev <- TRUE
input_2over3$ctrl$est@args$pa_buffer <- TRUE
input_2over3$ctrl$est@args$comp_f <- FALSE
input_2over3$ctrl$est@args$comp_b <- FALSE
input_2over3$ctrl$isys@args$upper_constraint <- 1.2
input_2over3$ctrl$isys@args$lower_constraint <- 0.8
input_2over3$ctrl$isys@args$cap_below_b <- TRUE

saveRDS(input_2over3, file = paste0("input/input_", n, "_2over3.rds"))

if (FALSE) {
  #debugonce(goFish)
  set.seed(1)
  res_2over3 <- do.call(mp, input_2over3)
  saveRDS(res_2over3, file = paste0("output/res_", n, "_2over3.rds"))
  plot(res_2over3, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
}

### ------------------------------------------------------------------------ ###
### FLXSA ####
### ------------------------------------------------------------------------ ###

FLXSA_control <- FLXSA.control(fse = 1.0, rage = 0, qage = 6, 
                               shk.n = FALSE, 
                               shk.ages = 3, shk.yrs = 3, 
                               min.nse = 0.3, 
                               tspower = 0, tsrange = 100,
                               maxit = 100)


input_FLXSA <- input
input_FLXSA$ctrl$est@args$FLXSA <- TRUE
input_FLXSA$ctrl$est@args$FLXSA_control <- FLXSA_control
input_FLXSA$ctrl$est@args$FLXSA_landings <- TRUE ### landings only?
input_FLXSA$ctrl$est@args$FLXSA_idcs <- c("Q1SWBeam", "FSP-7e")
input_FLXSA$ctrl$est@args$FLXSA_stf <- TRUE
input_FLXSA$ctrl$est@args$FLXSA_Btrigger <- 2443
input_FLXSA$ctrl$est@args$FLXSA_Ftrigger <- 0.238
input_FLXSA$oem@args$length_idx <- FALSE
input_FLXSA$oem@args$PA_status <- TRUE
input_FLXSA$oem@args$PA_status_dev <- TRUE
input_FLXSA$ctrl$est@args$pa_buffer <- TRUE
input_FLXSA$ctrl$est@args$comp_f <- FALSE
input_FLXSA$ctrl$est@args$comp_b <- FALSE
input_FLXSA$ctrl$isys@args$upper_constraint <- 1.2
input_FLXSA$ctrl$isys@args$lower_constraint <- 0.8
input_FLXSA$ctrl$isys@args$cap_below_b <- TRUE
input_FLXSA$ctrl$hcr@args$interval <- 1
input_FLXSA$ctrl$isys@args$interval <- 1

saveRDS(input_FLXSA, file = paste0("input/input_", n, "_FLXSA.rds"))

if (FALSE) {
  debugonce(input_FLXSA$ctrl$est@method)
  set.seed(1)
  res_FLXSA <- do.call(mp, input_FLXSA)
  plot(res_FLXSA, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
}


### ------------------------------------------------------------------------ ###
### prepare OM for category 1 assessment ####
### ------------------------------------------------------------------------ ###
newtonsteps <- 0
rel.tol = 0.001

### reference points
refpts_mse <- list(Btrigger = 2954,
                   Ftrgt = 0.241,
                   Fpa = 0.392,
                   Bpa = 2954,
                   Blim = 2110)
### some specifications for short term forecast with SAM
stf_def <- list(fwd_yrs_rec_start = 1980,
                fwd_splitLD = TRUE,
                fwd_yrs_average = -4:0,
                fwd_yrs_sel = -4:0)

### some arguments (passed to mp())
args <- list(fy = dims(stk_fwd)$maxyear, ### final simulation year
             y0 = dims(stk_fwd)$minyear, ### first data year
             iy = 2020, ### first simulation (intermediate) year
             nsqy = 3, ### not used, but has to provided
             nblocks = 1, ### block for parallel processing
             seed = 1 ### random number seed before starting MSE
)

### operating model
om <- FLom(stock = stk_fwd, ### stock 
           sr = sr, ### stock recruitment and pre-compiled residuals
           projection = mseCtrl(method = fwd_attr, 
                                args = list(maxF = 5,
                                          proc_res = "fitted" ### process noise
                                ))
)

### observation (error) model
oem <- FLoem(method = obs_generic,
             observations = list(stk = stk_oem, 
                                 idx = idx[c("Q1SWBeam", "FSP-7e")]), 
             deviances = list(stk = FLQuants(catch.dev = catch_res), 
                              idx = idx_dev[c("Q1SWBeam", "FSP-7e")]),
             args = list(cut_idx = TRUE,
                         idx_timing = c(-1, -1),
                         catch_timing = -1,
                         use_catch_residuals = TRUE, 
                         use_idx_residuals = TRUE,
                         use_stk_oem = TRUE))
### default management
ctrl_obj <- mpCtrl(list(
  est = mseCtrl(method = SAM_wrapper,
                args = c(### short term forecast specifications
                  forecast = TRUE, 
                  fwd_trgt = list(c("fsq", "fsq", "fsq")), fwd_yrs = 2, 
                  stf_def,
                  newtonsteps = newtonsteps, rel.tol = rel.tol,
                  par_ini = list(pars_ini),
                  track_ini = TRUE, 
                  conf = list(conf)
                )),
  phcr = mseCtrl(method = phcr_WKNSMSE,
                 args = refpts_mse),
  hcr = mseCtrl(method = hcr_WKNSME, args = list(option = "A")),
  isys = mseCtrl(method = is_WKNSMSE, 
                 args = c(hcrpars = list(refpts_mse),
                          fwd_trgt = list(c("fsq", "fsq", "hcr")), fwd_yrs = 3,
                          stf_def
                 ))
))
### additional tracking metrics - not used but required
tracking_add <- c("BB_return", "BB_bank_use", "BB_bank", "BB_borrow")

### save mse objects
input <- list(om = om, oem = oem, ctrl = ctrl_obj, args = args, 
              tracking = tracking_add)
saveRDS(object = input, 
        file = paste0(input_path, "input_SAM.rds"))

if (FALSE) {
  input$args$nblocks <- 1
  # debugonce(input$oem@method)
  # debugonce(goFish)
  #debugonce(input$ctrl$isys@method)
  set.seed(1)
  res_SAM <- do.call(mp, input)
}
