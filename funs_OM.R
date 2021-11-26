### ------------------------------------------------------------------------ ###
### functions for creating operating models (OMs) ####
### ------------------------------------------------------------------------ ###

### ------------------------------------------------------------------------ ###
### convert Beverton-Holt stock-recruitment model ####
### ------------------------------------------------------------------------ ###
### formulation with steepness, virgin biomass -> a & b

bevholtSV_to_bevholt <- function(sr) {
  sr_new <- sr ### duplicate model
  model(sr_new) <- "bevholt" ### change model type to bevholt
  ### convert parameters into normal bevholt parameters
  sr_pars <- abPars(model = "bevholt", spr0 = (params(sr)["spr0"]),
                    s = (params(sr)["s"]), v = (params(sr)["v"]))
  dimnames(sr_pars$a)$params <- "a"
  dimnames(sr_pars$b)$params <- "b"
  ### combine all parameters and insert them
  sr_pars <- rbind(rbind(sr_pars$a, sr_pars$b), params(sr))
  params(sr_new) <- sr_pars
  ### also insert some more slots
  ssb(sr_new) <- ssb(sr)
  rec(sr_new) <- rec(sr)
  logLik(sr_new) <- logLik(sr)
  details(sr_new) <- details(sr)
  residuals(sr_new) <- residuals(sr)
  fitted(sr_new) <- fitted(sr)
  
  return(sr_new)
}

### ------------------------------------------------------------------------ ###
### create operating model ####
### ------------------------------------------------------------------------ ###
### this function creates the elements required for an OM to run an MSE,
### e.g. OM stock, stock-recruitment model, survey indices, etc.
create_OM <- function(stk_data, idx_data, 
                      n = 1000, n_years = 100, yr_data = 2020, int_yr = FALSE,
                      SAM_conf, SAM_conf_full = FALSE, 
                      SAM_idx_weight = FALSE,
                      SAM_newtonsteps = 0, SAM_rel.tol = 0.001,
                      n_sample_yrs = 5, 
                      sr_model = "bevholtSV", sr_start = NULL,
                      sr_parallel = 10, sr_ar_check = TRUE,
                      process_error = TRUE, catch_oem_error = TRUE,
                      idx_weights = c("none"), ### "none"/"catch.wt"/"stock.wt"
                      idxB = 1, ### FALSE, index name or numeric index
                      idxL = TRUE, ALKs, ALK_yrs = NULL,
                      length_samples = 2000,
                      PA_status = TRUE, ### status evaluation success rate
                      refpts = list(),
                      stock_id = "ple.27.7e",
                      OM = "baseline",
                      save = TRUE,
                      return = FALSE,
                      M_alternative = NULL, ### 1 value (const.) or M@age
                      M_dd = FALSE, ### density dependent M
                      M_dd_relation = NULL,
                      M_dd_yr = NULL, ### last key run
                      M_dd_migration = NULL, ### "correct" ages 3+ for migration
                      disc_survival = 0, ### discard survival proportion
                      disc_survival_hidden = TRUE
                      ) {
  
  ### ---------------------------------------------------------------------- ###
  ### preparation for alternative OMs ####
  stk_data_input <- stk_data
  
  ### ---------------------------------------------------------------------- ###
  ### alternative M scenario? ####
  if (!is.null(M_alternative)) {
    message("using alternative M scenario")
    m(stk_data)[] <- M_alternative
  }
  
  ### ---------------------------------------------------------------------- ###
  ### alternative discard survival scenario? ####
  if (isTRUE(disc_survival > 0)) {
    message("using alternative discard survival scenario")
    discards.n(stk_data)[is.na(discards.n(stk_data))] <- 0
    discards.wt(stk_data)[is.na(discards.wt(stk_data))] <- 0
    discards.n(stk_data)[] <- discards.n(stk_data) * (1 - disc_survival)
    discards(stk_data) <- computeDiscards(stk_data)
    catch(stk_data) <- computeCatch(stk_data, slot = "all")
  }
  
  ### ---------------------------------------------------------------------- ###
  ### fit SAM ####
  message("fitting SAM")
  fit <- FLR_SAM(stk_data, idx_data, conf = SAM_conf, conf_full = SAM_conf_full,
                 idx_weight = SAM_idx_weight)
  ### fit SAM with relaxed convergence
  fit_mse <- FLR_SAM(stk_data, idx_data, conf = SAM_conf,
                     conf_full = SAM_conf_full, idx_weight = SAM_idx_weight,
                     newtonsteps = SAM_newtonsteps, rel.tol = SAM_rel.tol)
  ### get initial parameters
  pars_ini <- getpars(fit_mse)
  
  ### ---------------------------------------------------------------------- ###
  ### create FLStock ####
  message("create OM FLStock")
  stk <- SAM2FLStock(object = fit, stk = stk_data)
  stk_orig <- stk
  
  ### ---------------------------------------------------------------------- ###
  ### projection years ####
  yrs_hist <- as.numeric(dimnames(stk)$year)
  yrs_proj <- seq(from = dims(stk)$maxyear + 1, length.out = n_years)
  if (isTRUE(int_yr)) {
    yrs_hist <- yrs_hist[-length(yrs_hist)]
    yrs_proj <- seq(from = dims(stk)$maxyear + 0, length.out = n_years)
    n_years_project <- n_years - 1
  }
  yrs_mse <- sort(unique(c(yrs_hist, yrs_proj)))
  
  ### ---------------------------------------------------------------------- ###
  ### add uncertainty with variance-covariance matrix ####
  message("add uncertainty to OM with variance-covariance matrix")
  
  ### add iteration dimension
  stk <- FLCore::propagate(stk, n)
  ### add uncertainty estimated by SAM as iterations
  set.seed(1)
  uncertainty <- SAM_uncertainty(fit = fit, n = n)
  ### add noise to stock
  stock.n(stk)[] <- uncertainty$stock.n
  stock(stk)[] <- computeStock(stk)
  ### add noise to F
  harvest(stk)[] <- uncertainty$harvest
  ### add noise to catch numbers
  catch.n(stk)[, ac(yrs_hist)] <- uncertainty$catch.n
  catch(stk) <- computeCatch(stk)
  
  ### ---------------------------------------------------------------------- ###
  ### alternative discard survival scenario ####
  ### by default, the surviving discards are included in discards.n slot
  if (isTRUE(disc_survival > 0) & isTRUE(disc_survival_hidden)) {
    ### use catch.n uncertainty from SAM for landings and discards
    ### get catch.n uncertainty (SAM estimate/input data)
    c_noise <- catch.n(stk)/catch.n(stk_data)
    ### discards ratio in data
    dratio <- discards.n(stk_data_input)/catch.n(stk_data_input)
    dratio[is.na(dratio)] <- 0
    ### landings ratio in data
    lratio <- landings.n(stk_data_input)/catch.n(stk_data_input)
    ### total (after accounting for discard survival)
    ctotal <- dratio*(1 - disc_survival) + lratio
    ### update landings
    landings.n(stk)[] <- catch.n(stk) * lratio/ctotal ### catch.n includes noise
    ### discards relative to landings -> used to reproduce original discards
    dlratio <- discards.n(stk_data_input)/landings.n(stk_data_input)
    ### update discards (including surviving discards)
    discards.n(stk)[] <- landings.n(stk) * dlratio
    ### update total catch (including surviving discards)
    catch.n(stk)[] <- catch.n(stk_data_input) * c_noise
    ### get original catch weights (these were changed for SAM)
    catch.wt(stk)[] <- catch.wt(stk_data_input)
    ### update totals
    catch(stk) <- computeCatch(stk)
    landings(stk) <- computeLandings(stk)
    discards(stk) <- computeDiscards(stk)
  }
  
  ### ---------------------------------------------------------------------- ###
  ### extend stock for MSE simulation ####
  message("extend OM stock for projection")
  stk_stf <- stf(stk, n_years_project)
  
  ### ---------------------------------------------------------------------- ###
  ### biological data for OM ####
  message("create biological and fishery data for projection")
  ### Resample weights, maturity and natural mortality from the last 5 years 
  ### set up an array with one resampled year for each projection year 
  ### (including intermediate year) and replicate
  ### use the same resampled year for all biological parameters
  set.seed(2)
  ### use last five data years to sample biological parameters
  sample_yrs <- seq(to = yr_data, length.out = n_sample_yrs)
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
  bio_yrs <- which(dimnames(stk_stf)$year %in% 
                     (yr_data + 1):dims(stk_stf)$maxyear)
  ### insert values
  catch.wt(stk_stf)[, bio_yrs] <- c(catch.wt(stk)[, bio_samples,,,, 1])
  stock.wt(stk_stf)[, bio_yrs] <- c(stock.wt(stk)[, bio_samples,,,, 1])
  landings.wt(stk_stf)[, bio_yrs] <- c(landings.wt(stk)[, bio_samples,,,, 1])
  discards.wt(stk_stf)[, bio_yrs] <- c(discards.wt(stk)[, bio_samples,,,, 1])
  m(stk_stf)[, bio_yrs] <- c(m(stk)[, bio_samples,,,, 1])
  mat(stk_stf)[, bio_yrs] <- c(mat(stk)[, bio_samples,,,, 1])
  ### use different samples for selectivity
  harvest(stk_stf)[, bio_yrs] <- c(harvest(stk)[, sel_samples,,,, 1])
  
  ### density dependent M
  if (isTRUE(M_dd) & isTRUE(int_yr)) {
    ### update M in intermediate year 
    m(stk_stf)[, ac(yr_data + 1)] <- 
      calculate_ddM(stk_stf, yr_data + 1, relation = M_dd_relation,
                    migration = M_dd_migration)
  }
  
  ### ---------------------------------------------------------------------- ###
  ### stock recruitment ####
  message("creating stock-recruitment model")
  ### fit stock-recruitment model and get residuals from smoothed residuals
  
  ### create FLSR object
  sr <- as.FLSR(stk_stf, model = sr_model)
  if (!is.null(sr_start)) sr <- window(sr, start = sr_start)
  ### fit model individually to each iteration and suppress output to screen
  if (isFALSE(sr_parallel) | isTRUE(sr_parallel == 0)) {
    sr <- fmle(sr, method = 'L-BFGS-B', control = list(trace = 0))
  } else {
    ### run in parallel
    message("fitting stock-recruitment model in parallel")
    cl_tmp <- makeCluster(as.numeric(sr_parallel))
    registerDoParallel(cl_tmp)
    sr <- fmle_parallel(sr, cl_tmp, method = 'L-BFGS-B')
    stopCluster(cl_tmp)
  }
  ### run again for failed iterations - if needed
  pos_error <- which(is.na(params(sr)[1]))
  if (isTRUE(length(pos_error) > 0)) {
    message("repeating failed iterations")
    sr_corrected <- fmle(FLCore::iter(sr, pos_error), method = 'L-BFGS-B', 
                         control = list(trace = 0))
    sr[,,,,, pos_error] <- sr_corrected[]
    params(sr)[, pos_error] <- params(sr_corrected)
  }
  if (identical(sr_model, "bevholtSV")) {
    sr <- bevholtSV_to_bevholt(sr)
  }
  
  ### check autocorrelation of residuals for SAM median perception
  sr_med <- as.FLSR(stk_orig, model = sr_model)
  if (!is.null(sr_start)) sr_med <- window(sr_med, start = sr_start)
  sr_med <- fmle(sr_med, method = 'L-BFGS-B', control = list(trace = 0))
  sr_acf <- acf(residuals(sr_med), plot = FALSE, na.action = na.exclude)
  sr_rho <- sr_acf$acf[2]
  ### only include if lag-1 auto-correlation is above threshold
  ci <- qnorm((1 + 0.95)/2)/sqrt(sr_acf$n.used)
  if (isTRUE(sr_rho >= ci)) {
    message(paste0("residual lag-1 autocorrelation ",
                   round(sr_rho, 2), " above threshold of ", round(ci, 2)))
  } else {
    message(paste0("residual lag-1 autocorrelation ",
                   round(sr_rho, 2), " below threshold of ", round(ci, 2)))
  }
  if (isFALSE(sr_ar_check) | sr_rho < ci) {
    message("- NOT including auto-correlation")
  } else {
    message("- including auto-correlation")
  }

  ### generate residuals for MSE
  ### years with missing residuals
  yrs_res <- colnames(rec(sr))[which(is.na(iterMeans(rec(sr))))]
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
    ### "add" autocorrelation
    if (isTRUE(sr_ar_check) & isTRUE(sr_rho >= ci)) {
      sr_acf_i <- acf(res_i, lag.max = 1, plot = FALSE, na.action = na.exclude)
      sr_rho_i <- sr_acf_i$acf[2]
      res_ac <- rep(0, length(yrs_res))
      res_ac[1] <- sr_rho * tail(res_i, 1) + sqrt(1 - sr_rho^2) * res_new[1]
      for (r in 2:length(res_ac)) {
        res_ac[r] <- sr_rho * res_ac[r - 1] + sqrt(1 - sr_rho^2) * res_new[r]
      }
      res_new <- res_ac
    }
    return(res_new)
  }
  ### insert into model
  residuals(sr)[, yrs_res] <- unlist(res_new)
  ### exponentiate residuals to get factor
  residuals(sr) <- exp(residuals(sr))
  sr_res <- residuals(sr)
  
  ### ---------------------------------------------------------------------- ###
  ### process noise ####
  ### create FLQuant with process noise
  ### this will be added to the values obtained from fwd() in the MSE
  if (isTRUE(process_error)) {
    message("including survival process error")
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
    proc_res[, dimnames(proc_res)$year <= yr_data] <- 1
    ### remove deviation for first age class (recruits)
    proc_res[1, ] <- 1
  } else {
    message("NOT including survival process error")
    proc_res <- 1
  }
  ### try saving in stock recruitment model ... 
  ### this gets passed on to the projection module
  fitted(sr) <- proc_res
  
  ### ---------------------------------------------------------------------- ###
  ### stf ####
  ### for intermediate year
  ### not done for plaice
  stk_fwd <- stk_stf
  
  ### ---------------------------------------------------------------------- ###
  ### biological data for OEM ####
  message("OEM biological data")
  ### base on OM stock
  stk_oem <- stk_fwd
  
  ### if alternative M scenario, MP does not know this
  if (!is.null(M_alternative)) {
    M_yrs <- dimnames(stk_data_input)$year
    m(stk_oem)[, M_yrs] <- m(stk_data_input)
  }
  
  ### projection years
  proj_yrs <- (yr_data + 1):range(stk_oem)[["maxyear"]]
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
  ### density dependent M
  if (isTRUE(M_dd)) {
    ### last key run was in M_dd_yr, with the last data year M_dd_yr - 1
    ### so M M_dd_yr:M_dd_yr+2 will be the mean of M_dd_yr-1:M_dd_yr-3 of OM
    ### and M M_dd_yr+3:M_dd_yr+5 mean of M_dd_yr:M_dd_yr+2
    m(stk_oem)[, ac(M_dd_yr:(M_dd_yr + 2))] <- 
      yearMeans(m(stk_oem)[, ac((M_dd_yr - 1):(M_dd_yr - 3))])
  }
  
  ### ---------------------------------------------------------------------- ###
  ### catch noise ####
  ### take estimates from sam: uncertainty$catch_sd is "logSdLogObs"
  ### assume catch observed by SAM in projection is log-normally distributed
  ### around operating model catch
  if (isTRUE(catch_oem_error)) {
    message("including catch observation error")
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
    if (isFALSE(disc_survival > 0)) {
      catch_res[, dimnames(catch_res)$year <= yr_data] <- 
        window(catch.n(stk_orig), end = yr_data) / 
        window(catch.n(stk_fwd), end = yr_data)
    } else {
      catch_res[, dimnames(catch_res)$year <= yr_data] <- 
        window(catch.n(stk_data_input), end = yr_data) / 
        window(catch.n(stk_fwd), end = yr_data)
    }
  } else {
    message("NOT including catch observation error")
    catch_res <- catch.n(stk_fwd) %=% 1
  }
  
  ### ---------------------------------------------------------------------- ###
  ### indices ####
  message("creating OM indices")
  ### use real FLIndices object as template (use all)
  idx <- idx_data
  ### extend for simulation period
  idx <- window(idx, end = yr_data + n_years)
  ### add iterations
  idx <- lapply(idx, propagate, n)
  ### extract some dimension names
  idx_yrs <- lapply(idx, function(x) as.numeric(dimnames(x)$year))
  idx_yrs_hist <- lapply(idx_yrs, function(x) setdiff(x, yrs_proj))
  idx_ages <- lapply(idx, function(x) dimnames(x)$age)
  ### set catchability for projection
  for (idx_i in seq_along(idx)) {
    index.q(idx[[idx_i]])[] <- uncertainty$survey_catchability[[idx_i]]
  }
  ### index weights
  if (isTRUE(length(idx_weights) < length(idx))) {
    idx_weights <- rep(idx_weights, length(idx))[seq(length(idx))]
  }
  for (idx_i in seq_along(idx)) {
    if (isTRUE(idx_weights[idx_i] == "none")) {
      next
    ### take weights from OM stock
    } else if (isTRUE(idx_weights[idx_i] %in% c("stock.wt", "catch.wt"))) {
      get.weights_i <- get(x = idx_weights[idx_i])
      catch.wt(idx[[idx_i]]) <- 
        get.weights_i(stk_oem)[ac(idx_ages[[idx_i]]), ac(idx_yrs[[idx_i]])]
    ### use weights from index - resample
    } else if (isTRUE(idx_weights[idx_i] %in% c("index.wt"))) {
      catch.wt(idx[[idx_i]])[, ac(proj_yrs)] <- 
        yearMeans(catch.wt(idx[[idx_i]])[, ac(sample_yrs)])
    } 
  }
  ### create copy of index with original values
  idx_raw <-  lapply(idx ,index)
  ### calculate index numbers
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
  for (idx_i in seq_along(idx_dev)) {
    idx_dev[[idx_i]][, ac(idx_yrs_hist[[idx_i]])] <-
      idx_raw[[idx_i]][, ac(idx_yrs_hist[[idx_i]])] /
      index(idx[[idx_i]])[, ac(idx_yrs_hist[[idx_i]])]
  }

  ### add template for biomass index
  if (!isFALSE(idxB)) {
    message("adding biomass index template")
    if (is.numeric(idxB)) idxB <- names(idx)[idxB]
    idxB_template <- quantSums(idx[[idxB]]@catch.wt * 
                                 idx[[idxB]]@index) %=% NA_real_
    idx <- FLIndices(c(idx, idxB = FLIndex(index = idxB_template)))
    idx_dev[["idxB"]] <- idxB_template
  }
  
  ### ---------------------------------------------------------------------- ###
  ### length index ####
  if (isTRUE(idxL)) {
    message("adding length index")
    if (!is.null(ALK_yrs))
      ALKs <- ALKs %>%
        filter(year %in% ALK_yrs)
    ALKs <- ALKs %>%
      arrange(year, age, length)
    ### randomly select ALK years
    set.seed(89)
    alk_samples <- catch(stk_fwd) %=% NA_real_
    alk_samples[] <- sample(x = ALK_yrs, size = length(yrs_mse) * n, 
                            replace = TRUE)
    ### use existing ALKs for historical years
    alk_samples[, ac(ALK_yrs)] <- ALK_yrs
    ### pre-populate index
    set.seed(91)
    data_yrs <- range(stk_fwd)[["minyear"]]:yr_data
    ### match ALK year with data year
    lmean <- full_join(
      ### use observed catch
      x = as.data.frame(catch.n(stk_fwd)[, ac(data_yrs)]) %>%
        select(year, age, iter, data) %>%
        mutate(data = ifelse(data < 0, 0, data)) %>% ### shouldn't happen...
        rename("caa" = "data") %>%
        mutate(iter = as.numeric(as.character(iter))) %>%
        mutate(caa = caa + .Machine$double.eps), ### avoid 0s
      y = as.data.frame(alk_samples[, ac(data_yrs)]) %>%
        select(year, iter, data) %>%
        rename("alk_year" = "data") %>%
        mutate(iter = as.numeric(as.character(iter))), 
      by = c("year", "iter")) %>%
      ### merge with ALKs
      left_join(ALKs %>% rename("alk_year" = "year"),
                by = c("age", "alk_year")) %>%
      ### calculate numbers at length
      mutate(cal = caa * freq) %>%
      ### keep only numbers where length >= Lc
      filter(length >= refpts$Lc) %>% 
      ### mean catch length above Lc
      group_by(year, iter) %>%
      summarise(data = mean(sample(x = length, prob = cal, 
                                   size = length_samples, replace = TRUE)),
                .groups = "keep") %>%
      arrange(as.numeric(as.character(iter)))
    ### convert into FLQuant
    idxL <-  window(as.FLQuant(lmean), end = dims(stk_fwd)$maxyear)
    ### add to index
    idx <- FLIndices(c(idx, idxL = FLIndex(index = idxL)))
    idx_dev[["idxL"]] <- idxL %=% 1
    idx_dev[["alk_yrs"]] <- alk_samples
  }
  
  ### ------------------------------------------------------------------------ ###
  ### PA buffer for 2 over 3 rule ####
  ### SPiCT performance based on 
  ### Fischer et al. 2021 https://doi.org/10.1093/icesjms/fsab018
  ### index deviation
  if (isTRUE(PA_status)) {
    message("adding template for PA status")
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
    idx_dev[["PA_status"]] <- PA_status_dev
  }
  
  ### ---------------------------------------------------------------------- ###
  ### format reference points ####
  refpts <- FLPar(refpts, unit = "")
  
  ### ---------------------------------------------------------------------- ###
  ### save ####
  if (isTRUE(save)) {
    
    ### path
    input_path <- paste0("input/", stock_id, "/", OM, "/", n, "_", n_years, "/")
    dir.create(input_path, recursive = TRUE)
    
    ### SAM model fit
    saveRDS(fit, file = paste0(input_path, "SAM_fit.rds"))
    ### stock
    saveRDS(stk_fwd, file = paste0(input_path, "stk.rds"))
    ### stock recruitment
    saveRDS(sr, file = paste0(input_path, "sr.rds"))
    ### surveys
    saveRDS(idx, file = paste0(input_path, "idx.rds"))
    saveRDS(idx_dev, file = paste0(input_path, "idx_dev.rds"))
    ### catch noise
    saveRDS(catch_res, file = paste0(input_path, "catch_res.rds"))
    ### process error
    saveRDS(proc_res, file = paste0(input_path, "proc_res.rds"))
    ### observed stock
    saveRDS(stk_oem, file = paste0(input_path, "stk_oem.rds"))
    ### sam initial parameters
    saveRDS(pars_ini, file = paste0(input_path, "SAM_initial.rds"))
    ### sam configuration
    saveRDS(SAM_conf, file = paste0(input_path, "SAM_conf.rds"))
    ### reference values
    saveRDS(refpts, file = paste0(input_path, "refpts_mse.rds"))
    ### age-length keys
    saveRDS(ALKs, file = paste0(input_path, "ALKs.rds"))
    ### density dependent M
    if (isTRUE(M_dd)) {
      saveRDS(list(dd_M_relation = dd_M_relation, M_dd_yr = M_dd_yr,
                   M_dd_migration = M_dd_migration), 
              file = paste0(input_path, "dd_M.rds"))
    }
  }
  if (isTRUE(return)) {
    return(list(stk_fwd = stk_fwd, sr = sr, idx = idx, idx_dev = idx_dev,
                catch_res = catch_res, proc_res = proc_res, stk_oem = stk_oem,
                pars_ini = pars_ini, SAM_conf = SAM_conf, refpts = refpts,
                ALKs = ALKs))
  }
  
}

### ------------------------------------------------------------------------ ###
### create input for mp() ####
### ------------------------------------------------------------------------ ###
### this function loads the elements required for running mse::mp(),
### adapts them if necessary (dimensions), sets the MP and 
### creates the input object for mp()

input_mp <- function(stock_id = "ple.27.7e", OM = "baseline", n_iter = 1000,
                     n_yrs = 100, yr_start = 2021, iy = yr_start - 1,
                     n_blocks = 1, seed = 1, cut_hist = TRUE, MP = "rfb",
                     migration = NULL,
                     disc_survival = 0, rec_failure = FALSE,
                     use_age_idcs = NULL, biomass_index = NULL,
                     idx_timing = NULL, catch_timing = NULL,
                     fwd_yrs_rec_start = NULL,
                     fwd_splitLD = NULL,
                     fwd_yrs_average = NULL,
                     fwd_yrs_sel = NULL,
                     fwd_trgt = NULL, fwd_yrs = NULL
                     ) {
  
  ### path to input objects
  path_input <- paste0("input/", stock_id, "/", OM, "/1000_100/")
  
  ### load objects
  ### use full dimensions (100 years, 1000 iterations) - reduced later
  stk_fwd <- readRDS(paste0(path_input, "stk.rds"))
  sr <- readRDS(paste0(path_input, "sr.rds"))
  idx <- readRDS(paste0(path_input, "idx.rds"))
  idx_dev <- readRDS(paste0(path_input, "idx_dev.rds"))
  catch_res <- readRDS(paste0(path_input, "catch_res.rds"))
  proc_res <- readRDS(paste0(path_input, "proc_res.rds"))
  stk_oem <- readRDS(paste0(path_input, "stk_oem.rds"))
  SAM_pars_ini <- readRDS(paste0(path_input, "SAM_initial.rds"))
  SAM_conf <- readRDS(paste0(path_input, "SAM_conf.rds"))
  ALKs <- readRDS(paste0(path_input, "ALKs.rds"))
  refpts_mse <- readRDS(paste0(path_input, "refpts_mse.rds"))
  if (identical(OM, "M_dd")) 
    dd_M <- readRDS(paste0(path_input, "dd_M.rds"))
  
  ### find last year
  yr_end <- yr_start + n_yrs - 1
  ### reduce dimensions, if requested
  if (isTRUE(n_yrs < 100)) {
    ### OM stock
    stk_fwd <- window(stk_fwd, end = yr_end)
    ### OEM stock
    stk_oem <- window(stk_oem, end = yr_end)
    ### stock-recruitment model - need to go through elements individually
    rec(sr) <- window(rec(sr), end = yr_end)
    ssb(sr) <- window(ssb(sr), end = yr_end)
    residuals(sr) <- window(residuals(sr), end = yr_end)
    fitted(sr) <- window(fitted(sr), end = yr_end)
    range(sr)[["maxyear"]] <- yr_end
    ### index
    idx <- window(idx, end = yr_end)
    ### residuals
    idx_dev <- window(idx_dev, end = yr_end)
    catch_res <- window(catch_res, end = yr_end)
    proc_res <- window(proc_res, end = yr_end)
  }
  if (isTRUE(n_iter < 1000)) {
    ### OM stock
    stk_fwd <- FLCore::iter(stk_fwd, seq(n_iter))
    ### OEM stock
    stk_oem <- FLCore::iter(stk_oem, seq(n_iter))
    ### stock-recruitment model 
    sr <- FLCore::iter(sr, seq(n_iter))
    ### index
    idx <- FLCore::iter(idx, seq(n_iter))
    ### residuals
    idx_dev <- FLCore::iter(idx_dev, seq(n_iter))
    catch_res <- FLCore::iter(catch_res, seq(n_iter))
    proc_res <- FLCore::iter(proc_res, seq(n_iter))
    ### reference points
    if (isTRUE(n_iter < dims(refpts_mse)$iter))
      refpts_mse <- iter(refpts_mse, seq(n_iter))
  }
  
  ### ---------------------------------------------------------------------- ###
  ### recruitment failure? ####
  ### ---------------------------------------------------------------------- ###
  if (!isFALSE(rec_failure)) {
    
    residuals(sr)[, ac(rec_failure)] <- residuals(sr)[, ac(rec_failure)] * 0.1
    
  }
  
  ### ---------------------------------------------------------------------- ###
  ### generic arguments ####
  args <- list(fy = yr_end, ### final simulation year
               y0 = dims(stk_fwd)$minyear, ### first data year
               iy = iy, ### first simulation (intermediate) year
               nsqy = 3, ### not used, but has to provided
               nblocks = n_blocks, ### block for parallel processing
               seed = seed ### random number seed before starting MSE
  )
  
  ### ---------------------------------------------------------------------- ###
  ### reference values ####
  # refpts_mse <- FLPar(refpts_mse, iter = n_iter)
  
  ### ---------------------------------------------------------------------- ###
  ### Operating model (OM) ####
  om <- FLom(stock = stk_fwd, ### stock 
             sr = sr, ### stock recruitment and precompiled residuals
             projection = mseCtrl(method = fwd_attr, 
                                  args = list(maxF = 5,
                                              ### process noise on stock.n
                                              proc_res = "fitted",
                                              dupl_trgt = FALSE,
                                              disc_survival = disc_survival
                                  ))
  )
  
  ### ---------------------------------------------------------------------- ###
  ### migration ####
  if (!is.null(migration)) {
    if (isTRUE(length(migration) == 1)) {
      migr_fct <- seq(from = 1, to = migration, length.out = n_yrs + 1)
      migr_fct <- migr_fct[-1]/migr_fct[-length(migr_fct)]
      #prod(migr_fct)
    } else {
      migr_fct <- migration
    }
    migr_qnt <- ssb(stk_fwd) %=% 1
    migr_qnt[, ac(yr_start:dims(migr_qnt)$maxyear)] <- migr_fct
    ### store in sr ssb slot
    om@sr@ssb <- migr_qnt
    om@projection@args$migration = "ssb"
  }
  
  ### ---------------------------------------------------------------------- ###
  ### Observation (error) model OEM ####
  if (identical(stock_id, "ple.27.7e")) {
    if (is.null(use_age_idcs)) use_age_idcs <- c("Q1SWBeam", "FSP-7e")
    if (is.null(biomass_index)) biomass_index <- "FSP-7e"
    if (is.null(idx_timing)) idx_timing <- c(-1, -1)
    if (is.null(catch_timing)) catch_timing <- -1
  } else if (identical(stock_id, "cod.27.47d20")) {
    if (is.null(use_age_idcs)) 
      use_age_idcs <- c("IBTS_Q1_gam", "IBTS_Q3_gam", "IBTS_Q3_gam_age0")
    if (is.null(biomass_index)) biomass_index <- "IBTS_Q3_gam"
    if (is.null(idx_timing)) idx_timing <- c(0, -1, 0)
    if (is.null(catch_timing)) catch_timing <- -1
  }
  
  ### default oem for rfb rule
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
                           alks = as.data.frame(ALKs),
                           use_age_idcs = use_age_idcs,
                           biomass_index = biomass_index,
                           length_idx = TRUE,
                           Lc = unique(c(refpts_mse["Lc"])),
                           lngth_samples = 2000))
  
  ### 2 over 3 rule: biomass index and PA status relative to OM + error
  if (isTRUE(MP == "2over3")) {
    oem@args$length_idx <- FALSE ### no length index
    oem@args$PA_status <- TRUE 
    oem@args$PA_status_dev <- TRUE
    oem@args$PA_Bmsy <- unique(c(refpts_mse["Bmsy"])) ### real MSY from OM
    oem@args$PA_Fmsy <- unique(c(refpts_mse["Fmsy"]))
    
    ### 2 over 3 based on XSA 
  } else if (isTRUE(MP == "2over3_XSA")) {
    oem@args$length_idx <- FALSE ### no length index
    oem@args$PA_status <- TRUE 
    oem@args$PA_status_dev <- TRUE
  } else if (isTRUE(MP == "ICES_SAM")) {
    oem@observations$idx <- oem@observations$idx[use_age_idcs] ### age indices only
    oem@deviances$idx <- oem@deviances$idx[use_age_idcs]
    oem@args <- list(cut_idx = TRUE, idx_timing = idx_timing,
                     catch_timing = catch_timing, use_catch_residuals = TRUE,
                     use_idx_residuals = TRUE, use_stk_oem = TRUE)
  } else if (isTRUE(MP == "constF")) {
    oem <- FLoem(observations = list(stk = FLQuant(0),
                                     idx = FLQuant()),
                 deviances = list(stk = FLQuant(0),
                                  idx = FLQuant()))
  }
  
  ### ---------------------------------------------------------------------- ###
  ### Management procedure MP ####
  
  if (isTRUE(MP == "rfb")) {
    ### I_trigger = 1.4 * I_loss
    I_trigger = apply(quantSums(index(idx[[biomass_index]]) *
                                  idx_dev[[biomass_index]] * 
                                  catch.wt(idx[[biomass_index]])),
                      6, min, na.rm = TRUE) * 1.4
    ctrl <- mpCtrl(list(
      est = mseCtrl(method = est_comps,
                    args = list(comp_r = TRUE, comp_f = TRUE, comp_b = TRUE,
                                comp_c = TRUE, comp_m = 1,
                                idxB_lag = 1, idxB_range_1 = 2, idxB_range_2 = 3,
                                idxB_range_3 = 1,
                                catch_lag = 0, ### 0 to mimic advice
                                catch_range = 1,
                                Lref = unique(c(refpts_mse["Lref"])), 
                                I_trigger = c(I_trigger),
                                idxL_lag = 1, idxL_range = 1)),
      phcr = mseCtrl(method = phcr_comps,
                     args = list(exp_r = 1, exp_f = 1, exp_b = 1)),
      hcr = mseCtrl(method = hcr_comps,
                    args = list(interval = 2)),
      isys = mseCtrl(method = is_comps,
                     args = list(interval = 2, 
                                 upper_constraint = 1.2, lower_constraint = 0.7, 
                                 cap_below_b = FALSE))
    ))
  } else if (isTRUE(MP == "2over3")) {
    ctrl <- mpCtrl(list(
      est = mseCtrl(method = est_comps,
                    args = list(comp_r = TRUE, comp_f = FALSE, comp_b = FALSE,
                                comp_c = TRUE, comp_m = 1, 
                                idxB_lag = 1, idxB_range_1 = 2, idxB_range_2 = 3,
                                catch_lag = 0, ### 0 to mimic advice
                                catch_range = 1, 
                                pa_buffer = TRUE, 
                                pa_size = 0.8, pa_duration = 3)),
      phcr = mseCtrl(method = phcr_comps),
      hcr = mseCtrl(method = hcr_comps,
                    args = list(interval = 2)),
      isys = mseCtrl(method = is_comps,
                     args = list(interval = 2, 
                                 upper_constraint = 1.2, lower_constraint = 0.8, 
                                 cap_below_b = TRUE))))
  } else if (isTRUE(MP == "2over3_XSA")) {
    FLXSA_control <- FLXSA.control(fse = 1.0, rage = 0, qage = 6, 
                                   shk.n = FALSE, 
                                   shk.ages = 3, shk.yrs = 3, 
                                   min.nse = 0.3, 
                                   tspower = 0, tsrange = 100,
                                   maxit = 100)
    ctrl <- mpCtrl(list(
      est = mseCtrl(method = est_comps,
                    args = list(comp_r = TRUE, comp_f = FALSE, comp_b = FALSE,
                                comp_c = TRUE, comp_m = 1, pa_buffer = TRUE,
                                idxB_lag = 1, 
                                idxB_range_1 = 2, idxB_range_2 = 3,
                                catch_lag = 0, ### 0 to mimic advice
                                catch_range = 1,
                                FLXSA = TRUE,
                                FLXSA_control = FLXSA_control,
                                FLXSA.control = NULL,
                                FLXSA_landings = TRUE, ### landings only?
                                FLXSA_idcs = c("Q1SWBeam", "FSP-7e"),
                                FLXSA_stf = TRUE,
                                FLXSA_Btrigger = unique(c(refpts_mse["ICES_Btrigger"])), ### from WGCSE
                                FLXSA_Ftrigger = unique(c(refpts_mse["ICES_Fmsy"]))
                    )),
      phcr = mseCtrl(method = phcr_comps),
      hcr = mseCtrl(method = hcr_comps,
                    args = list(interval = 1)), ### annual
      isys = mseCtrl(method = is_comps,
                     args = list(interval = 2, ### annual
                                 upper_constraint = 1.2, lower_constraint = 0.8, 
                                 cap_below_b = TRUE))))
  } else if (isTRUE(MP == "ICES_SAM")) {
    ### some specifications for short term forecast with SAM
    if (identical(stock_id, "ple.27.7e")) {
      if (is.null(fwd_yrs_rec_start)) fwd_yrs_rec_start <- 1980
      if (is.null(fwd_splitLD)) fwd_splitLD <- TRUE
      if (is.null(fwd_yrs_average)) fwd_yrs_average <- -4:0
      if (is.null(fwd_yrs_sel)) fwd_yrs_sel <- -4:0
      if (is.null(fwd_trgt)) fwd_trgt <- c("fsq", "fsq", "fsq")
      if (is.null(fwd_yrs)) fwd_yrs <- 2
    } else if (identical(stock_id, "cod.27.47d20")) {
      if (is.null(fwd_yrs_rec_start)) fwd_yrs_rec_start <- 1998
      if (is.null(fwd_splitLD)) fwd_splitLD <- TRUE
      if (is.null(fwd_yrs_average)) fwd_yrs_average <- -2:0
      if (is.null(fwd_yrs_sel)) fwd_yrs_sel <- NULL
      if (is.null(fwd_trgt)) fwd_trgt <- c("fsq")
      if (is.null(fwd_yrs)) fwd_yrs <- 1
    }
    SAM_stf_def <- list(fwd_yrs_rec_start = fwd_yrs_rec_start,
                        fwd_splitLD = fwd_splitLD,
                        fwd_yrs_average = fwd_yrs_average,
                        fwd_yrs_sel = fwd_yrs_sel)
    ctrl <- mpCtrl(list(
      est = mseCtrl(method = SAM_wrapper,
                    args = c(### short term forecast specifications
                      forecast = TRUE, 
                      fwd_trgt = list(fwd_trgt), fwd_yrs = fwd_yrs, 
                      SAM_stf_def, ### without list structure
                      newtonsteps = 0, rel.tol = 0.001,
                      par_ini = list(SAM_pars_ini),
                      track_ini = TRUE, 
                      conf = list(SAM_conf)
                    )),
      phcr = mseCtrl(method = phcr_WKNSMSE,
                     args = list(
                       Btrigger = unique(c(refpts_mse["EqSim_Btrigger"])), 
                       Ftrgt = unique(c(refpts_mse["EqSim_Fmsy"])), 
                       Blim = unique(c(refpts_mse["EqSim_Blim"])))),
      hcr = mseCtrl(method = hcr_WKNSME, args = list(option = "A")),
      isys = mseCtrl(method = is_WKNSMSE, 
                     args = c(hcrpars = list(
                       Btrigger = unique(c(refpts_mse["EqSim_Btrigger"])), 
                       Ftrgt = unique(c(refpts_mse["EqSim_Fmsy"])), 
                       Blim = unique(c(refpts_mse["EqSim_Blim"]))),
                       fwd_trgt = list(c(fwd_trgt, "hcr")), 
                       fwd_yrs = fwd_yrs + 1, SAM_stf_def
                     ))))
  } else if (isTRUE(MP == "constF")) {
    ctrl <- mpCtrl(list(hcr = mseCtrl(method = fixedF.hcr,
                                      args = list(ftrg = 0))))
  }
  
  ### ---------------------------------------------------------------------- ###
  ### density dependent M ####
  if (identical(OM, "M_dd")) {
    ### adapt projection argument of OM
    om@projection@args$dd_M <- TRUE
    om@projection@args$dd_M_relation <- dd_M$dd_M_relation
    om@projection@args$dd_M_fun <- calculate_ddM
    ### adapt observations - key runs
    ### only applicable when using SAM
    if (isTRUE(MP == "ICES_SAM")) {
      oem@args$dd_M <- TRUE
      oem@args$dd_M_relation <- dd_M$dd_M_relation
      oem@args$dd_M_yr <- dd_M$M_dd_yr
    }
  }
  
  ### ---------------------------------------------------------------------- ###
  ### additional tracking metrics ####
  tracking <- c("comp_c", "comp_i", "comp_r", "comp_f", "comp_b",
                "multiplier", "exp_r", "exp_f", "exp_b")
  if (isTRUE(MP == "ICES_SAM")) {
    tracking <- c("BB_return", "BB_bank_use", "BB_bank", "BB_borrow")
  }
  
  ### ---------------------------------------------------------------------- ###
  ### create input object ####
  
  input <- list(om = om, oem = oem, ctrl = ctrl,
                args = args, tracking = tracking, refpts = refpts_mse,
                cut_hist = cut_hist)
  
  ### within scenario parallelisation?
  if (isTRUE(n_blocks > 1)) {
    input$args$nblocks <- n_blocks
  }
  
  return(input)
  
}

### ------------------------------------------------------------------------ ###
### function for estimating MSY reference values ####
### ------------------------------------------------------------------------ ###
# res <- est_MSY(n_iter = 10, n_blocks = 1, vals_ini = c(0, 0.5, 1), tol = 0.1)
est_MSY <- function(stock_id = "ple.27.7e", OM = "baseline",
                    yr_start = 2021, n_blocks = 10, n_iter = 1000,
                    vals_ini = seq(0, 1, 0.1),
                    lower = 0, upper = 0.3, tol = 0.001,
                    plot = TRUE, x_label = "F (ages 3-6)",
                    save = TRUE) {
  #browser()
  path <- path_input <- paste0("input/", stock_id, "/", OM, "/", n_iter, 
                               "_100/")
  
  ### load mp input
  input <- input_mp(n_iter = n_iter, stock_id = stock_id, OM = OM,
                    yr_start = yr_start, MP = "constF", n_blocks = n_blocks)
  
  ### check if some results already exist
  if (isTRUE(save)) {
    if (isTRUE(file.exists(paste0(path, "MSY_trace.rds")))) {
      res_trace_ini <- readRDS(paste0(path, "MSY_trace.rds"))
    } else {
      res_trace_ini <- list()
    }
  }
  ### create object in new environment for storing results
  trace_env <- new.env()
  assign(x = "res_trace", value = res_trace_ini, envir = trace_env)
  
  ### define function for running projection and returning stats
  mp_catch <- function(input, Ftrgt, minimise = FALSE) {#browser()
    res_trace_i <- get("res_trace", envir = trace_env)
    ### if run already exists, do not run again
    if (isTRUE(Ftrgt %in% sapply(res_trace_i, function(x) x$Ftrgt))) {
      res_i <- res_trace_i[[which(Ftrgt == sapply(res_trace_i, 
                                                  function(x) x$Ftrgt))[[1]]]]
      catch_i <- res_i$catch
      ssb_i <- res_i$ssb
      tsb_i <- res_i$ssb
      rec_i <- res_i$rec
    } else {
      ### run projection
      input$ctrl$hcr@args$ftrg <- Ftrgt
      res_i <- do.call(mp, input)
      catch_i <- median(tail(catch(res_i@stock), 10), na.rm = TRUE)
      ssb_i <- median(tail(ssb(res_i@stock), 10), na.rm = TRUE)
      tsb_i <- median(tail(tsb(res_i@stock), 10), na.rm = TRUE)
      rec_i <- median(tail(rec(res_i@stock), 10), na.rm = TRUE)
    }
    ### print results of current run
    cat(paste0("Ftrgt=", Ftrgt, "; C=", catch_i, "; SSB=", ssb_i, "; R=", rec_i,
               "\n"))
    ### save results in res_trace
    res_add <- list(list(Ftrgt = Ftrgt, catch = catch_i,
                         ssb = ssb_i, tsb = tsb_i, rec = rec_i))
    assign(value = append(get("res_trace", envir = trace_env), res_add), 
           x = "res_trace", envir = trace_env)
    if (isTRUE(minimise)) catch_i <- -catch_i
    return(catch_i)
  }
  
  ### first, check some values
  names(vals_ini) <- vals_ini
  res_ini <- lapply(vals_ini, mp_catch, input = input)
  
  ### use optimise - 1D golden-section search
  res_optimise_MSY <- optimise(f = mp_catch, input = input,
                               interval = c(lower, upper), 
                               lower = lower, upper = upper,
                               maximum = TRUE,
                               tol = tol)
  
  ### get results and format
  res_trace_list <- unique(get("res_trace", envir = trace_env))
  res_trace <- as.data.frame(do.call(rbind, res_trace_list))
  res_trace <- as.data.frame(apply(res_trace, 2, unlist))
  res_trace <- unique(res_trace)
  
  
  if (isTRUE(plot)) {
    ### plot 
    p <- res_trace %>%
      select(Ftrgt, catch, ssb, rec) %>%
      pivot_longer(cols = c("catch", "ssb", "rec")) %>%
      mutate(value = value/1000,
             name = factor(name, levels = c("catch", "ssb", "rec"), 
                           labels = c("Catch [1000t]", "SSB [1000t]",
                                      "Recruitment [millions]"))) %>%
      ggplot(aes(x = Ftrgt, y = value)) +
      geom_point(size = 0.8) +
      stat_smooth(aes(alpha = "loess smoother"), size = 0.5,
                  se = FALSE, span = 0.3, n = 100, show.legend = TRUE) + 
      scale_alpha_manual("", values = 1) +
      facet_wrap(~ name, scales = "free_y", strip.position = "left") +
      labs(x = x_label, y = "") +
      ylim(c(0, NA)) +
      scale_x_continuous(breaks = seq(0, 1, 0.2)) +
      theme_bw() +
      theme(legend.position = c(0.85, 0.2),
            legend.background = element_blank(),
            legend.key = element_blank(),
            strip.placement = "outside",
            strip.background = element_blank(),
            axis.title.y = element_blank())
    ggsave(paste0(path, "MSY_search.png"), plot = p,
           width = 17, height = 6, units = "cm", dpi = 300, type = "cairo")
    ggsave(paste0(path, "MSY_search.pdf"), plot = p,
           width = 17, height = 6, units = "cm")
  } else {
    p <- NULL
  }
  if (isTRUE(save)) {
    saveRDS(res_trace_list, file = paste0(path, "MSY_trace.rds"))
    write.csv(res_trace, file = paste0(path, "MSY_trace.csv"), 
              row.names = FALSE)
  }
  
  return(list(plot = p, result = res_trace))
  
}




