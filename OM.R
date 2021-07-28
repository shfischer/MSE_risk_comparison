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

source("funs.R")

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
n <- 1000
### number of projection years
n_years <- 100
### last data year
yr_data <- 2020

### ------------------------------------------------------------------------ ###
### create FLStock ####
### ------------------------------------------------------------------------ ###
### create template with 1 iteration

stk <- SAM2FLStock(object = fit, stk = stk_data)
summary(stk)

plot(stk)

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
dim(stk)

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

plot(stk, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))

### maximum observed F
max(fbar(stk))
max(harvest(stk))

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

plot(stk_stf)

### ------------------------------------------------------------------------ ###
### stock recruitment ####
### ------------------------------------------------------------------------ ###
### fit hockey-stick model
### get residuals from smoothed residuals

### create FLSR object
sr <- as.FLSR(stk_stf, model = "segreg")
### fit model individually to each iteration and suppress output to screen
# suppressWarnings(. <- capture.output(sr <- fmle(sr)))

### run in parallel
library(doParallel)
cl <- makeCluster(10)
registerDoParallel(cl)
sr <- fmle_parallel(sr, cl)
stopCluster(cl)
registerDoSEQ()
### run again for failed iterations
pos_error <- which(is.na(params(sr)["a"]))
# sr_corrected <- fmle(FLCore::iter(sr, pos_error))
# sr[,,,,, pos_error] <- sr_corrected[]
# params(sr)[, pos_error] <- params(sr_corrected)

### generate residuals for MSE
### years with missing residuals
yrs_res <- colnames(rec(sr))[which(is.na(iterMeans(rec(sr))))]

### go through iterations and create residuals
### use kernel density to create smooth distribution of residuals
### and sample from this distribution
res_new <- foreach(iter_i = seq(dim(sr)[6]), .packages = "FLCore", 
                   .errorhandling = "pass") %do% {
  
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
  
  return(res_new)
  
}
summary(exp(unlist(res_new)))
### insert into model
residuals(sr)[, yrs_res] <- unlist(res_new)
### exponeniate residuals to get factor
residuals(sr) <- exp(residuals(sr))
sr_res <- residuals(sr)

plot(sr_res)

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
### the proc_res values are on a normale scale,
### exponentiate to get log-normal 
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

plot(proc_res)

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

plot(catch_res, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))

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

### index weights - use stock weigths
idx$Q1SWBeam@catch.wt <- stock.wt(stk_fwd)[ac(idx_Q1_ages), ac(idx_Q1_yrs)]
idx$`FSP-7e`@catch.wt <- stock.wt(stk_fwd)[ac(idx_FSP_ages), ac(idx_FSP_yrs)]

### create copy of index with original values
idx_raw <-  lapply(idx ,index)
### calculate index values
idx <- calc_survey(stk = stk_fwd, idx = idx, use_q = TRUE, use_time = TRUE, 
                   use_biomass = TRUE)

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
# idx_dev$Q1SWBeam[, ac(idx_Q1_yrs_hist)] <- 
#   idx_raw$Q1SWBeam[, ac(idx_Q1_yrs_hist)]/
#   index(idx$Q1SWBeam)[, ac(idx_Q1_yrs_hist)]
# idx_dev$`FSP-7e`[, ac(idx_FSP_yrs_hist)] <- 
#   idx_raw$`FSP-7e`[, ac(idx_FSP_yrs_hist)]/
#   index(idx$`FSP-7e`)[, ac(idx_FSP_yrs_hist)]

### check biomass index
plot(quantSums(idx$Q1SWBeam@catch.wt * idx$Q1SWBeam@index))
### including uncertainty
plot(quantSums(idx$Q1SWBeam@catch.wt * idx$Q1SWBeam@index * idx_dev$Q1SWBeam))

### add template for biomass index
idx <- FLIndices("Q1SWBeam" = idx$Q1SWBeam,
                 "FSP-7e" = idx$`FSP-7e`,
                 "idxB" = FLIndex(index = quantSums(idx$Q1SWBeam@catch.wt * 
                                                      idx$Q1SWBeam@index)))


### ------------------------------------------------------------------------ ###
### length index ####
### ------------------------------------------------------------------------ ###

### load age-length keys
ALKs <- read.csv("input/ALKs.csv", as.is = TRUE)
names(ALKs)[1] <- "year"
ALKs <- ALKs %>% 
  pivot_longer(c(X1:X26), names_to = "age", 
               values_to = "value", names_prefix = "X", 
               values_ptypes = list(length = as.numeric()),
               values_drop_na = TRUE) %>%
  mutate(age = as.integer(age)) %>%
  ### remove data, quarter, sex
  group_by(year, length, age) %>%
  summarise(value = sum(value)) %>%
  ungroup() %>%
  ### keep only lengths with >= 1 fish
  filter(value > 0)
#str(ALKs)

dims(stk)$min
dims(stk)$max

### use age 10 as plusgroup
### and combine ages 1 & 2
# ALKs <- ALKs %>%
#   filter(age <= 10)
ALKs <- ALKs %>% 
  mutate(age = ifelse(age <= 10, age, 10)) %>%
  mutate(age = ifelse(age <= 2, 2, age)) %>%
  group_by(year, length, age) %>%
  summarise(value = sum(value))

### standardise frequencies per age
ALKs <- ALKs %>%
  group_by(year, age) %>%
  ### total numbers per age
  mutate(total = sum(value)) %>%
  ### get frequency
  ungroup() %>%
  mutate(freq = value/total, value = NULL, total = NULL)

# ALKs %>%
#   group_by(year, age) %>%
#   summarise(total = sum(freq)) %>%
#   summary()

### OM starts at age 2
ALKs <- ALKs %>%
  filter(age >= 2)

### keep only last 5 years
ALKs <- ALKs %>%
  filter(year %in% 2016:2020)

### sort
ALKs <- ALKs %>%
  arrange(year, age, length)

# ### check number of length samples from InterCatch
# table1 <- read.csv("input/table1_hist.txt", as.is = TRUE)
# 
# ### extract necessary data
# samples <- table1 %>%
#   filter(CATONRaisedOrImported == "Imported_Data" &
#            SampledOrEstimated == "Sampled_Distribution") %>%
#   select(Year, Country, CatchCategory, Fleet, Season, 
#          No..of.Length.Samples, No..of.Length.Measured) %>% 
#   filter(No..of.Length.Measured > 0)
# ### length readings per year
# samples %>%
#   group_by(Year) %>%
#   summarise(readings = sum(No..of.Length.Measured),
#             samples = sum(No..of.Length.Samples))


### scale up
### define which ALKs are used
alk_yrs <- 2016:2020
### random samples
set.seed(89)
alk_samples <- catch(stk_fwd)
alk_samples[] <- sample(x = alk_yrs, size = length(yrs_mse) * n, replace = TRUE)
### use existing ALKs for historical years
alk_samples[, ac(alk_yrs)] <- alk_yrs

Lc <- 26

### pre-populate index
set.seed(91)
data_yr <- 1980:2020
stk <- stk_fwd
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
### keep only numbers where length >= Lc
cn <- cn %>%
  filter(length >= Lc)
### mean catch length above Lc
means <- cn %>%
  group_by(year, iter) %>%
  summarise(data = mean(sample(x = length, prob = cal, 
                                 size = 100, replace = TRUE))) %>%
  arrange(as.numeric(as.character(iter)))
### convert into index
idxL <- FLIndex(index = as.FLQuant(means))
### extend for projection
idxL <- window(idxL, end = dims(stk_fwd)$maxyear)

### ------------------------------------------------------------------------ ###
### PA buffer for 2 over 3 rule ####
### ------------------------------------------------------------------------ ###

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

save.image(paste0("input/image_", n, ".RData"))
load(paste0("input/image_", n, ".RData"))

### ------------------------------------------------------------------------ ###
### prepare OM for MSE ####
### ------------------------------------------------------------------------ ###
library(mse)
source("funs.R")

### length target
Lref <- rep(0.75*26 + 0.25*66, n)
### I_trigger = 1.4 * I_loss
I_trigger = apply(quantSums(index(idx$Q1SWBeam) * idx_dev$Q1SWBeam), 
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
    idx = FLIndices(Q1SWBeam = idx$Q1SWBeam, 
                    "FSP-7e" = idx$`FSP-7e`,
                    idxB = idx$idxB,
                    idxL = idxL,
                    PA_status = PA_status_template)), 
  deviances = list(
    stk = FLQuants(catch.dev = catch_res), 
    idx = FLQuants(Q1SWBeam = idx_dev$Q1SWBeam,
                   "FSP-7e" = idx_dev$`FSP-7e`,
                   alk_yrs = alk_samples,
                   PA_status = PA_status_dev)),
  args = list(use_catch_residuals = TRUE, 
              use_idx_residuals = TRUE,
              use_stk_oem = TRUE,
              use_biomass = TRUE,
              PA_status = FALSE, PA_status_dev = FALSE,
              PA_Bmsy = PA_Bmsy,
              PA_Fmsy = PA_Fmsy,
              alks = as.data.frame(ALKs),
              use_age_idcs = c("Q1SWBeam", "FSP-7e"),
              biomass_index = "Q1SWBeam",
              length_idx = TRUE,
              Lc = 26,
              lngth_samples = 100))
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

### save mse objects
input <- list(om = om, oem = oem, ctrl = ctrl,
              args = args, tracking = tracking, cut_hist = FALSE)

saveRDS(input, file = paste0("input/input_", n, ".rds"))


#input$args$nblocks <- 100
#debugonce(input$oem@method)
#debugonce(goFish)
#debugonce(input$ctrl$isys@method)
set.seed(1)
res <- do.call(mp, input)
saveRDS(res, file = paste0("output/res_", n, "_rfb.rds"))
plot(res, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))


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

#debugonce(goFish)
set.seed(1)
res_2over3 <- do.call(mp, input_2over3)
saveRDS(res_2over3, file = paste0("output/res_", n, "_2over3.rds"))
plot(res_2over3, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))

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


debugonce(input_FLXSA$ctrl$est@method)
set.seed(1)
res_FLXSA <- do.call(mp, input_FLXSA)
plot(res_FLXSA, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
