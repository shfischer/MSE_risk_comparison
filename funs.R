### ------------------------------------------------------------------------ ###
### run fmle for FLSR in parallel ####
### ------------------------------------------------------------------------ ###

fmle_parallel <- function(sr, cl, seed = NULL, ...) {
  
  ### split into parts
  it_parts <- split(1:dim(sr)[6], cut(1:dim(sr)[6], length(cl)))
  names(it_parts) <- NULL
  
  ### split sr into junks
  sr_chunks <- lapply(it_parts, function(x) {
    FLCore::iter(sr, x)
  })
  
  ### set RNG seed
  if (!is.null(seed)) {
    clusterEvalQ(cl = cl, expr = {set.seed(seed)})
  }
  
  ### fit model to each junk
  fits <- foreach(sr_chunk = sr_chunks, .packages = "FLCore") %dopar% {
    
    fmle(sr_chunk, ...)
    
  }
  
  ### combine into single object
  ### slots rec, ssb, covar, logerrror, model, logl, gr, distribution, inital,
  ### name, desc, range
  ### do not need to be changed
  logLik(sr) <- logLik(fits[[1]])
  sr@logLik[seq(dims(sr)$iter)] <- c(sapply(fits, logLik))
  
  fitted(sr) <- do.call(FLCore::combine, lapply(fits, "fitted"))
  
  residuals(sr) <- do.call(FLCore::combine, lapply(fits, "residuals"))
  
  vcov(sr) <- array(data = unlist(lapply(fits, vcov)), 
                    dim = c(dim(fits[[1]]@vcov)[1], 
                            dim(fits[[1]]@vcov)[2], 
                            dims(sr)$iter), 
                    dimnames = list(dimnames(fits[[1]]@vcov)[[1]], 
                                    dimnames(fits[[1]]@vcov)[[2]], 
                                    iter = ac(dimnames(sr)$iter)))
  
  params(sr) <- propagate(params(sr), dims(sr)$iter)
  params(sr)[] <- unlist(lapply(fits, function(x) {
    split(x@params@.Data, slice.index(x@params, MARGIN = 2))
  }))
  
  hessian(sr) <- array(data = unlist(lapply(fits, hessian)), 
                       dim = c(dim(fits[[1]]@hessian)[1], 
                               dim(fits[[1]]@hessian)[2], 
                               dims(sr)$iter), 
                       dimnames = list(dimnames(fits[[1]]@hessian)[[1]], 
                                       dimnames(fits[[1]]@hessian)[[2]], 
                                       iter = ac(dimnames(sr)$iter)))
  
  details(sr) <- details(fits[[1]])
  
  return(sr)
  
}

### ------------------------------------------------------------------------ ###
### calculate survey index/indices from stock ####
### ------------------------------------------------------------------------ ###
#' calculate survey index/indices from FLStock
#'
#' This function calculates survey indices from the numbers at age of an 
#' FLStock object
#'
#' @param stk Object of class \linkS4class{FLStock} with stock and fishery data.
#' @param idx Object of class \linkS4class{FLIndices} or \linkS4class{FLIndex}.
#' @param use_q If \code{TRUE} index numbers are calculated by multiplying
#'   numbers at age in the stock with the catchability.
#' @param use_time If \code{TRUE} observed numbers in the survey are corrected
#'   for fishing and natural mortality.
#'
#' @return An object of class \code{FLIndex} or \code{FLIndex} with the survey
#'   index stored in the \code{index} slot.
#'
#' @export

setGeneric("calc_survey", function(stk, idx, use_q = TRUE, use_time = TRUE,
                                   use_wt = FALSE) {
  standardGeneric("calc_survey")
})

### stk = FLStock, idx = FLIndices
#' @rdname calc_survey
setMethod(f = "calc_survey",
          signature = signature(stk = "FLStock", idx = "FLIndices"),
          definition = function(stk, idx, use_q = TRUE, use_time = TRUE,
                                use_wt = FALSE) {
            
            ### apply function to every element of FLIndices
            lapply(X = idx, FUN = calc_survey_ind, stk = stk, use_q = use_q, 
                   use_time = use_time, use_wt = use_wt)
            
          })
### stk = FLStock, idx = FLIndex
#' @rdname calc_survey
setMethod(f = "calc_survey",
          signature = signature(stk = "FLStock", idx = "FLIndex"),
          definition = function(stk, idx, use_q = TRUE, use_time = TRUE,
                                use_wt = FALSE) {
            
            calc_survey_ind(stk = stk, idx = idx, use_q = use_q, 
                            use_time = use_time, use_wt = use_wt)
            
          })



calc_survey_ind <- function(stk, idx, 
                            use_q = TRUE, ### catchability
                            use_time = TRUE, ### timing of survey
                            use_wt = FALSE ### use weights?
) {
  
  ### find ranges for years, ages & iters
  ages <- intersect(dimnames(index(idx))$age, dimnames(stock.n(stk))$age)
  years <- intersect(dimnames(index(idx))$year, dimnames(stock.n(stk))$year)
  iter <- intersect(dimnames(index(idx))$iter, dimnames(stock.n(stk))$iter)
  
  ### timing of survey
  if (isTRUE(use_time)) {
    ### use mean of fishing period
    time <- mean(range(idx)[c("startf", "endf")])
  } else {
    ### otherwise assume beginning of year
    time <- 0
  }
  
  ### extract stock numbers for requested/available dimensions
  index.n <- stock.n(stk)[ac(ages), ac(years),,,, ac(iter)]
  ### get Z = M & F
  Z <- m(stk)[ac(ages), ac(years),,,, ac(iter)] +
    harvest(stk)[ac(ages), ac(years),,,, ac(iter)]
  
  ### estimate stock numbers at time of survey
  index.n <- index.n * exp(-time * Z)
  catch.n(idx)[ac(ages), ac(years),,,, ac(iter)] <- index.n
  
  ### add catchability, if requested
  if (isTRUE(use_q)) {
    index.n <- index.n * index.q(idx)[ac(ages), ac(years),,,, ac(iter)]
  }
  
  ### use biomass (i.e. include weights)?
  if (isTRUE(use_wt)) {
    index.n <- index.n * catch.wt(idx)[ac(ages), ac(years),,,, ac(iter)]
  }
  
  ### insert values into index
  index(idx)[ac(ages), ac(years),,,, ac(iter)] <- index.n
  
  return(idx)
  
}

### ------------------------------------------------------------------------ ###
### observations ####
### ------------------------------------------------------------------------ ###
obs_generic <- function(stk, observations, deviances, args, tracking,
                        biomass = FALSE,
                        simple_idx = FALSE, ### simple index calculation
                        ssb_idx = FALSE, tsb_idx = FALSE, ### use SSB idx
                        idx_dev = FALSE, ### add index deviation?
                        PA_status = FALSE, ### for PA buffer
                        PA_status_dev = FALSE,
                        PA_Bmsy = FALSE, PA_Fmsy = FALSE,
                        use_stk_oem = FALSE, ### biological parameters, wts etc
                        use_catch_residuals = FALSE,
                        use_idx_residuals = FALSE,
                        cut_idx = FALSE, ### cut off indices years after ay
                        catch_timing = -1,
                        idx_timing = -1,
                        use_wt = FALSE,
                        use_age_idcs = NULL, ### survey to use calc_survey()
                        alks = NULL, ### age-length keys
                        biomass_index = NULL, ### which survey used for biomass?
                        length_idx = FALSE, ### estimate length index?
                        lngth = FALSE, ### simpler length index
                        lngth_dev = FALSE, ### deviation for lngth
                        Lc = 0, ### length at first capture
                        lngth_samples = 100, ### number of length samples
                        dd_M = FALSE, ### density dependent M -> keyruns
                        dd_M_relation = NULL,
                        dd_M_yr = NULL,
                        dd_M_period = 3,
                        ...) {
  
  ay <- args$ay
  
  ### Density-dependent M
  ### Calculate 3-year means of M from the OM on key-run years
  ### to simulate the SAM process
  if (isTRUE(dd_M)) {
    if (isTRUE(((ay - dd_M_yr) %% dd_M_period) == 0)) {
      m(observations$stk)[, ac(ay:(ay + 2))] <- 
        yearMeans(m(stk)[, ac((ay - 3):(ay - 1))])
    }
  }
  
  ### create object for observed stock
  if (!isTRUE(use_stk_oem)) {
    
    ### use OM as template
    stk0 <- stk
    
  } else {
    ### otherwise use observed stock provided in observations object
    ### this can include different biological data, e.g. weights at age
    stk0 <- observations$stk
    
    ### update fishery data
    catch.n(stk0) <- catch.n(stk)
    catch(stk0) <- catch(stk)
    discards.n(stk0) <- discards.n(stk)
    discards(stk0) <- discards(stk)
    landings.n(stk0) <- landings.n(stk)
    landings(stk0) <- landings(stk)
    
  }
  
  ### add uncertainty to catch
  if (isTRUE(use_catch_residuals)) {
    
    ### implement for catch at age
    catch.n(stk0) <- catch.n(stk) * deviances$stk$catch.dev
    
    ### split catch into discards and landings, based on landing fraction
    landings.n(stk0) <- catch.n(stk0) * (landings.n(stk) / catch.n(stk))
    discards.n(stk0) <- catch.n(stk0) * (1 - landings.n(stk) / catch.n(stk))
    
    ### update total catch/discards/landings
    catch(stk0) <- computeCatch(stk0)
    landings(stk0) <- computeLandings(stk0)
    discards(stk0) <- computeDiscards(stk0)
    
  }
  
  ### calculate age indices
  if (is.null(use_age_idcs)) use_age_idcs <- names(observations$idx)
  idcs_new <- calc_survey(stk = stk, 
                          idx = observations$idx[use_age_idcs],
                          use_wt = use_wt)
  ### bug in FLIndices, need to iterate through them when replacing...
  for (i in use_age_idcs) {
    observations$idx[[i]] <- idcs_new[[i]]
  }
  
  ### observed survey
  idx0 <- observations$idx 
  
  ### add uncertainty to observed age indices
  if (isTRUE(use_idx_residuals)) {
    
    for (idx_i in use_age_idcs) {
      index(idx0[[idx_i]]) <- index(idx0[[idx_i]]) * deviances$idx[[idx_i]]
    }
    
  }
  
  ### cut off trailing years? used when running stock assessment
  if (isTRUE(cut_idx)) {
  
    ### cut off years
    ### workaround for NS cod: survey until intermediate year, but catch stops
    ### 1 year earlier
    ### slots such as natural mortality, maturity need to be kept, otherwise
    ### SAM will fall over
    if (any(idx_timing > catch_timing)) {
    
      ### keep stock until last survey data year
      stk0 <- window(stk0, end = ay + max(idx_timing))
      ### find years to remove
      yrs_remove <- (ay + catch_timing + 1):ay
      ### remove catch data
      catch(stk0)[, ac(yrs_remove)] <- NA
      catch.n(stk0)[, ac(yrs_remove)] <- NA
      catch.wt(stk0)[, ac(yrs_remove)] <- NA
      landings(stk0)[, ac(yrs_remove)] <- NA
      landings.n(stk0)[, ac(yrs_remove)] <- NA
      landings.wt(stk0)[, ac(yrs_remove)] <- NA
      discards(stk0)[, ac(yrs_remove)] <- NA
      discards.n(stk0)[, ac(yrs_remove)] <- NA
      discards.wt(stk0)[, ac(yrs_remove)] <- NA
    
    } else {
    
      stk0 <- window(stk0, end = ay + catch_timing)
    
    }
    
    ### timing of survey
    ### 0: intermediate/assessment year
    ### <0: fewer years available & vice versa
    if (length(idx_timing) < length(idx0)) { 
      idx_timing <- rep(idx_timing, length(idx0))
    }
    ### restrict years for indices based on timing
    idx0 <- lapply(seq_along(idx0), function(x) {
      window(idx0[[x]], end = ay + idx_timing[x])
    })
    idx0 <- FLIndices(idx0) ### restore class
    names(idx0) <- names(observations$idx)

  }
  
  ### create biomass index
  if (!is.null(biomass_index)) {
    idxB_yrs <- intersect(dimnames(index(observations$idx$idxB))$year,
                          dimnames(index(observations$idx[[biomass_index]]))$year)
    index(observations$idx$idxB)[, idxB_yrs] <- quantSums(index(observations$idx[[biomass_index]]))[, idxB_yrs]
    index(idx0$idxB)[, idxB_yrs] <- quantSums(index(idx0[[biomass_index]]))[, idxB_yrs]
  }
  
  ### catch length index
  if (isTRUE(length_idx)) {
    
    ### find years which need updating
    ### do not update beyond year ay
    ### always include/overwrite ay (in case MSE start earlier)
    yrs_update <- which(is.na(iterMeans(window(index(idx0$idxL), end = ay))))
    yrs_update <- unique(c(dimnames(idx0$idxL)$year[yrs_update], ay))
    ### match ALK year with data year
    lmean <- full_join(
      ### use observed catch
      x = as.data.frame(catch.n(stk0)[, ac(yrs_update)]) %>%
        select(year, age, iter, data) %>%
        mutate(data = ifelse(data < 0, 0, data)) %>% ### shouldn't happen...
        rename("caa" = "data") %>%
        mutate(iter = as.numeric(as.character(iter))) %>%
        mutate(caa = caa + .Machine$double.eps), ### avoid 0s
      y = as.data.frame(deviances$idx$alk_yrs[, ac(yrs_update)]) %>%
        select(year, iter, data) %>%
        rename("alk_year" = "data") %>%
        mutate(iter = as.numeric(as.character(iter))), 
      by = c("year", "iter")) %>%
      ### merge with ALKs
      left_join(alks %>% rename("alk_year" = "year"),
                by = c("age", "alk_year")) %>%
      ### calculate numbers at length
      mutate(cal = caa * freq) %>%
      ### keep only numbers where length >= Lc
      filter(length >= Lc) %>% 
      ### mean catch length above Lc
      group_by(year, iter) %>%
      summarise(data = mean(sample(x = length, prob = cal, 
                                   size = lngth_samples, replace = TRUE)),
                .groups = "keep") %>%
      arrange(as.numeric(as.character(iter)))
    
    ### save
    index(observations$idx$idxL[, ac(yrs_update)]) <- as.FLQuant(lmean)
    index(idx0$idxL) <- index(observations$idx$idxL)
    
  }
  
  ### use SSB as index?
  if (isTRUE(ssb_idx)) {
    index(observations$idx$idxB) <- ssb(observations$stk)
    ### TSB?
  } else  if (isTRUE(tsb_idx)) {
    index(observations$idx$idxB) <- tsb(observations$stk)
    ### otherwise calculate biomass index
  } else if (isTRUE(simple_idx)) {
    index(observations$idx$idxB) <- quantSums(stk@stock.n * stk@stock.wt * 
                                              index(observations$idx$sel))
  }
  
  ### simpler length index (for compatibility)
  if (isTRUE(lngth)) {
    index(idx0$idxL) <- index(observations$idx$idxL) <- 
      lmean(stk = stk, params = lngth_par)
  }
  
  ### stock status for PA buffer?
  if (isTRUE(PA_status)) {
    index(observations$idx$PA_status)[] <- ssb(stk) > 0.5*PA_Bmsy & 
      fbar(stk) < PA_Fmsy
  }

  ### add more deviances to index?
  if (isTRUE(idx_dev)) {
    if (isTRUE(ssb_idx) | isTRUE(tsb_idx)) {
      index(idx0$idxB) <- index(observations$idx$idxB) * deviances$idx$idxB
    } else {
      index(idx0$idxB) <- quantSums(stk@stock.n * stk@stock.wt * 
                               observations$idx$sel * deviances$idx$sel)
      if (isTRUE("idxB" %in% names(deviances$idx)) & 
          all.equal(dim(deviances$idx$idxB), dim(index(idx0$idxB))))
        index(idx0$idxB) <- index(idx0$idxB) * deviances$idx$idxB
    }
  }
  ### uncertainty for catch length
  if (isTRUE(lngth) & isTRUE(lngth_dev)) {
    index(idx0$idxL) <- index(observations$idx$idxL) * deviances$idx$idxL
  }
  ### uncertainty for stock status for PA buffer
  if (isTRUE(PA_status) & isTRUE(PA_status_dev)) {
    index(idx0$PA_status) <- ifelse(index(observations$idx$PA_status) == TRUE, 
                                    deviances$idx$PA_status["positive", ],
                                    deviances$idx$PA_status["negative", ])
  }
  
  return(list(stk = stk0, idx = idx0, observations = observations,
              tracking = tracking))
  
}

### ------------------------------------------------------------------------ ###
### estimator ####
### ------------------------------------------------------------------------ ###

est_comps <- function(stk, idx, tracking, args,
                      comp_r = FALSE, comp_f = FALSE, comp_b = FALSE,
                      comp_i = FALSE, comp_c = TRUE, comp_m = FALSE,
                      comp_hr = FALSE,
                      idxB_lag = 1, idxB_range_1 = 2, idxB_range_2 = 3,
                      idxB_range_3 = 1,
                      catch_lag = 1, catch_range = 1,
                      Lref, I_trigger,
                      idxL_lag = 1, idxL_range = 1,
                      pa_buffer = FALSE, pa_size = 0.8, pa_duration = 3,
                      Bmsy = NA,
                      FLXSA = FALSE, ### run FLXSA?
                      FLXSA_control = NULL, 
                      FLXSA_landings = FALSE, ### landings only?
                      FLXSA_idcs = NULL, ### indices to use for XSA
                      FLXSA_stf = FALSE, ### add short-term forecast?
                      FLXSA_Btrigger = 0, ### SSB trigger for PA buffer
                      FLXSA_Ftrigger = Inf, ### F trigger for PA buffer
                      ...) {
  
  ay <- args$ay
  
  ### prepare biomass index
  if (isTRUE(is(idx$idxB, "FLIndex"))) {
    idxB <- index(idx$idxB)
  } else {
    idxB <- idx$idxB
  }
  ### prepare length index
  if (isTRUE(is(idx$idxL, "FLIndex"))) {
    idxL <- index(idx$idxL)
  } else {
    idxL <- idx$idxL
  }
  
  if (isTRUE(FLXSA)) {
    ### prepare stock
    stk_FLXSA <- stk
    ### use only landings in assessment?
    if (isTRUE(FLXSA_landings)) {
      discards.n(stk_FLXSA) <- 0
      discards(stk_FLXSA) <- computeDiscards(stk_FLXSA)
      catch.n(stk_FLXSA) <- landings.n(stk_FLXSA)
      catch(stk_FLXSA) <- landings(stk_FLXSA)
    }
    ### prepare indices
    if (is.null(FLXSA_idcs)) FLXSA_idcs <- names(idx)
    idx_FLXSA <- idx[FLXSA_idcs]
    ### crop years
    stk_FLXSA <- window(stk_FLXSA, end = ay - 1)
    idx_FLXSA <- window(idx_FLXSA, end = ay - 1)
    ### run FLXSA
    xsa_res <- FLXSA(stock = stk_FLXSA, indices = idx_FLXSA, 
                     control = FLXSA_control)
    stk_FLXSA@stock.n[] <- xsa_res@stock.n
    stock(stk_FLXSA) <- computeStock(stk_FLXSA)
    stk_FLXSA@harvest[] <- xsa_res@harvest
    ### add short-term forecast to get SSB in following year?
    if (isTRUE(FLXSA_stf)) {
      ### set up control object
      stf_ctrl <- fwdControl(data.frame(year = ay,
                                        val = 0, ### dummy value
                                        quantity = "f"))
      ### define stf recruitment
      stf_rec <- apply(rec(stk_FLXSA), 6, function(x) exp(mean(log(x))))
      stf_rec <- FLPar(a = c(stf_rec), iter = dims(stk_FLXSA)$iter)
      ### extend stock by 1 year
      stk_FLXSA <- stf(stk_FLXSA, 1)
      ### forecast
      stk_FLXSA <- fwd(object = stk_FLXSA, control = stf_ctrl,
                     sr = list(model = "mean", 
                               params = stf_rec))
      dimnames(stk_FLXSA) <- list(iter = dimnames(stk)$iter) ### fix iter names
    }
    ### use SSB as stock index
    idxB <- ssb(stk_FLXSA)
    
    ### overwrite PA buffer indicator, base on FLXSA results
    ### compare SSB to MSYBtrigger and Fbar to Fmsy
    index(idx$PA_status)[] <- 0
    ### SSB status
    PA_status_SSB <- ssb(stk_FLXSA) >= FLXSA_Btrigger
    ### Fbar status
    PA_status_F <- fbar(stk_FLXSA) <= FLXSA_Ftrigger
    ### shift by 1 year because of stf
    if (isTRUE(FLXSA_stf)) {
      ### remove first year because of shift
      PA_status_SSB <- PA_status_SSB[, -1] 
      dimnames(PA_status_SSB)$year <- ac(an(dimnames(PA_status_SSB)$year) - 1)
      ### remove last year for F
      ### (no data, and to adapt dimension to PA_status_SSB)
      PA_status_F <- PA_status_F[, -dims(PA_status_F)$year]
    }
    ### combine SSB and F evaluation
    index(idx$PA_status)[, dimnames(PA_status_F)$year] <- 
      (PA_status_F & PA_status_SSB)
  }
  
  ### component r: index trend
  if (isTRUE(comp_r)) {
    r_res <- est_r(idx = idxB, ay = ay,
                   idxB_lag = idxB_lag, idxB_range_1 = idxB_range_1, 
                   idxB_range_2 = idxB_range_2)
  } else {
    r_res <- 1
  }
  tracking["comp_r", ac(ay)] <- r_res
  
  ### component f: length data
  if (isTRUE(comp_f)) {
    f_res <- est_f(idx = idxL, ay = ay,
                   Lref = Lref, idxL_range = idxL_range, idxL_lag = idxL_lag)
  } else {
    f_res <- 1
  }
  tracking["comp_f", ac(ay)] <- f_res
  
  ### component b: biomass safeguard
  if (isTRUE(comp_b)) {
    b_res <- est_b(idx = idxB, ay = ay,
                   I_trigger = I_trigger, idxB_lag = idxB_lag, 
                   idxB_range_3 = idxB_range_3)
  } else {
    b_res <- 1
  }
  
  ### PA buffer
  if (isTRUE(pa_buffer)) {
    b_res <- est_pa(idx = index(idx$PA_status), ay = ay, 
                    tracking = tracking, idxB_lag = idxB_lag,
                    pa_size = pa_size, pa_duration = pa_duration)
  }
  tracking["comp_b", ac(ay)] <- b_res
  
  ### component i: index value
  if (isTRUE(comp_i)) {
    i_res <- est_i(idx = idxB, ay = ay,
                   idxB_lag = idxB_lag, idxB_range_3 = idxB_range_3)
  } else {
    i_res <- 1
  }
  tracking["comp_i", ac(ay)] <- i_res
  
  ### current catch
  if (isTRUE(comp_c)) {
    c_res <- est_c(ay = ay, catch = catch(stk), catch_lag = catch_lag, 
                   catch_range = catch_range)
  } else {
    c_res <- 1
  }
  tracking["comp_c", ac(ay)] <- c_res
  
  ### component m: multiplier
  if (!isFALSE(comp_m)) {
    m_res <- comp_m
    ### subset to iteration when simultion is split into blocks
    if (isTRUE(length(comp_m) > dims(stk)$iter)) {
      m_res <- comp_m[as.numeric(dimnames(stk)$iter)]
    }
  } else {
    m_res <- 1
  }
  tracking["multiplier", ac(ay)] <- m_res
  
  ### component hr: harvest rate (catch/idx)
  if (!isFALSE(comp_hr)) {
    hr_res <- comp_hr
    ### subset to iteration when simultion is split into blocks
    if (isTRUE(length(comp_hr) > dims(stk)$iter)) {
      hr_res <- comp_hr[as.numeric(dimnames(stk)$iter)]
    }
  } else {
    hr_res <- 1
  }
  tracking["comp_hr", ac(ay)] <- hr_res
  
  return(list(stk = stk, tracking = tracking))
  
}

### biomass index trend
est_r <- function(idx, ay,
                  idxB_lag, idxB_range_1, idxB_range_2,
                  ...) {
  
  ### index ratio
  yrs_a <- seq(to = c(ay - idxB_lag), length.out = idxB_range_1)
  yrs_b <- seq(to = min(yrs_a) - 1, length.out = idxB_range_2)
  idx_a <- yearMeans(idx[, ac(yrs_a)])
  idx_b <- yearMeans(idx[, ac(yrs_b)])
  idx_ratio <- c(idx_a / idx_b)
  
  return(idx_ratio)
  
}

### length data
est_f <- function(idx, ay, 
                  Lref, idxL_range, idxL_lag,
                  ...) {
  
  ### if fewer iterations provided expand
  if (isTRUE(length(Lref) < dims(idx)$iter)) {
    Lref <- rep(Lref, dims(idx)$iter)
    ### if more iterations provided, subset
  } else if (isTRUE(length(Lref) > dims(idx)$iter)) {
    Lref <- Lref[an(dimnames(idx)$iter)]
  }
  
  ### get mean length in catch
  idx_yrs <- seq(to = ay - idxL_range, length.out = idxL_lag)
  idx_mean <- yearMeans(idx[, ac(idx_yrs)])
  ### length relative to reference
  idx_ratio <- c(idx_mean / Lref)
  ### avoid negative values
  idx_ratio <- ifelse(idx_ratio > 0, idx_ratio, 0)
  ### avoid NAs, happens if catch = 0
  idx_ratio <- ifelse(is.na(idx_ratio), 1, idx_ratio)
  return(idx_ratio)
}

### biomass index trend
est_b <- function(idx, ay, 
                  I_trigger, idxB_lag, idxB_range_3,
                  ...) {
  
  ### if fewer iterations provided expand
  if (isTRUE(length(I_trigger) < dims(idx)$iter)) {
    I_trigger <- rep(I_trigger, dims(idx)$iter)
    ### if more iterations provided, subset
  } else if (isTRUE(length(I_trigger) > dims(idx)$iter)) {
    I_trigger <- I_trigger[an(dimnames(idx)$iter)]
  }
  
  ### calculate index mean
  idx_yrs <- seq(to = ay - idxB_lag, length.out = idxB_range_3)
  idx_mean <- yearMeans(idx[, ac(idx_yrs)])
  ### ratio
  idx_ratio <- c(idx_mean / I_trigger)
  ### b is 1 or smaller
  idx_ratio <- ifelse(idx_ratio < 1, idx_ratio, 1)
  
  return(idx_ratio)
  
}

### biomass index trend
est_pa <- function(idx, ay, tracking, pa_size, pa_duration, idxB_lag,
                   ...) {
  
  ### find last year in which buffer was applied
  last <- apply(tracking["comp_b",,, drop = FALSE], 6, FUN = function(x) {#browser()
    ### positions (years) where buffer was applied
    yr <- dimnames(x)$year[which(x < 1)]
    ### return -Inf if buffer was never applied
    ifelse(length(yr) > 0, as.numeric(yr), -Inf)
  })
  ### find iterations to check 
  pos_check <- which(last <= (ay - pa_duration))
  ### find negative stock status (SSB<0.5Bmsy or F>Fmsy)
  pos_negative <- which(idx[, ac(ay - idxB_lag)] == 0)
  ### apply only if buffer applications need to be checked and status is negative
  pos_apply <- intersect(pos_check, pos_negative)
  
  return(ifelse(seq(dims(last)$iter) %in% pos_apply, pa_size, 1))
  
}

### index value
est_i <- function(idx, ay,
                  idxB_lag, idxB_range_3,
                  ...) {
  
  ### index ratio
  yrs_r <- seq(to = c(ay - idxB_lag), length.out = idxB_range_3)
  idx_i <- yearMeans(idx[, ac(yrs_r)])
  
  return(idx_i)
  
}

### recent catch
est_c <- function(catch, ay,
                  catch_lag, catch_range,
                  ...) {
  
  catch_yrs <- seq(to = ay - catch_lag, length.out = catch_range)
  catch_current <- yearMeans(catch[, ac(catch_yrs)])
  return(catch_current)
  
}

### ------------------------------------------------------------------------ ###
### phcr ####
### ------------------------------------------------------------------------ ###
### parametrization of HCR

phcr_comps <- function(tracking, args, 
                       exp_r = 1, exp_f = 1, exp_b = 1,
                       ...){
  
  ay <- args$ay
  
  hcrpars <- tracking[c("comp_r", "comp_f", "comp_b", "comp_i", 
                        "comp_hr", "comp_c", "multiplier",
                        "exp_r", "exp_f", "exp_b"), ac(ay)]
  hcrpars["exp_r", ] <- exp_r
  hcrpars["exp_f", ] <- exp_f
  hcrpars["exp_b", ] <- exp_b
  
  if (exp_r != 1) tracking["exp_r", ] <- exp_r
  if (exp_f != 1) tracking["exp_f", ] <- exp_f
  if (exp_b != 1) tracking["exp_b", ] <- exp_b
  
  ### return results
  return(list(tracking = tracking, hcrpars = hcrpars))
  
}

### ------------------------------------------------------------------------ ###
### hcr ####
### ------------------------------------------------------------------------ ###
### apply catch rule

hcr_comps <- function(hcrpars, args, tracking, interval = 2, 
                      ...) {
  
  ay <- args$ay ### current year
  iy <- args$iy ### first simulation year
  
  ### check if new advice requested
  if ((ay - iy) %% interval == 0) {
    
    ### calculate advice
    advice <- hcrpars["comp_c", ] *
      (hcrpars["comp_r", ]^hcrpars["exp_r", ]) *
      (hcrpars["comp_f", ]^hcrpars["exp_f", ]) *
      (hcrpars["comp_b", ]^hcrpars["exp_b", ]) *
      hcrpars["comp_i"] *
      hcrpars["comp_hr"] *
      hcrpars["multiplier", ] 
    #advice <- apply(X = hcrpars, MARGIN = 6, prod, na.rm = TRUE)
    
  } else {
    
    ### use last year's advice
    advice <- tracking["metric.hcr", ac(ay - 1)]
    
  }
  
  ctrl <- getCtrl(values = c(advice), quantity = "catch", years = ay + 1, 
                  it = dim(advice)[6])
  
  return(list(ctrl = ctrl, tracking = tracking))
  
}

### ------------------------------------------------------------------------ ###
### implementation ####
### ------------------------------------------------------------------------ ###
### no need to convert, already catch in tonnes
### apply TAC constraint, if required

is_comps <- function(ctrl, args, tracking, interval = 2, 
                     upper_constraint = Inf, lower_constraint = 0, 
                     cap_below_b = TRUE, ...) {
  
  ay <- args$ay ### current year
  iy <- args$iy ### first simulation year
  
  advice <- ctrl@trgtArray[ac(ay + args$management_lag), "val", ]
  
  ### check if new advice requested
  if ((ay - iy) %% interval == 0) {
    
    ### apply TAC constraint, if requested
    if (!is.infinite(upper_constraint) | lower_constraint != 0) {
      
      ### get last advice
      if (isTRUE(ay == iy)) {
        ### use OM value in first year of projection
        adv_last <- tracking["C.om", ac(iy)]
      } else {
        adv_last <- tracking["metric.is", ac(ay - 1)]
      }
      ### ratio of new advice/last advice
      adv_ratio <- advice/adv_last
      
      ### upper constraint
      if (!is.infinite(upper_constraint)) {
        ### find positions
        pos_upper <- which(adv_ratio > upper_constraint)
        ### turn of constraint when index below Itrigger?
        if (isFALSE(cap_below_b)) {
          pos_upper <- setdiff(pos_upper, 
                               which(c(tracking[, ac(ay)]["comp_b", ]) < 1))
        }
        ### limit advice
        if (length(pos_upper) > 0) {
          advice[pos_upper] <- adv_last[,,,,, pos_upper] * upper_constraint
        }
        ### lower constraint
      }
      if (lower_constraint != 0) {
        ### find positions
        pos_lower <- which(adv_ratio < lower_constraint)
        ### turn of constraint when index below Itrigger?
        if (isFALSE(cap_below_b)) {
          pos_lower <- setdiff(pos_lower, 
                               which(c(tracking[, ac(ay)]["comp_b", ]) < 1))
        }
        ### limit advice
        if (length(pos_lower) > 0) {
          advice[pos_lower] <- adv_last[,,,,, pos_lower] * lower_constraint
        }
      }
    }
    
    ### otherwise do nothing here and recycle last year's advice
  } else {
    
    advice <- tracking["metric.is", ac(ay - 1)]
    
  }
  ctrl@trgtArray[ac(ay + args$management_lag),"val",] <- advice
  
  return(list(ctrl = ctrl, tracking = tracking))
  
}

### ------------------------------------------------------------------------ ###
### implementation error ####
### ------------------------------------------------------------------------ ###

iem_comps <- function(ctrl, args, tracking, 
                      iem_dev = FALSE, use_dev, ...) {
  
  ay <- args$ay
  
  ### only do something if requested
  if (isTRUE(use_dev)) {
    
    ### get advice
    advice <- ctrl@trgtArray[ac(ay + args$management_lag), "val", ]
    ### get deviation
    dev <- c(iem_dev[, ac(ay)])
    ### implement deviation
    advice <- advice * dev
    ### insert into ctrl object
    ctrl@trgtArray[ac(ay + args$management_lag),"val",] <- advice
    
  }
  
  return(list(ctrl = ctrl, tracking = tracking))
  
}

### ------------------------------------------------------------------------ ###
### projection ####
### ------------------------------------------------------------------------ ###
fwd_attr <- function(stk, ctrl,
                     sr, ### stock recruitment model
                     sr.residuals, ### recruitment residuals
                     sr.residuals.mult = TRUE, ### are res multiplicative?
                     maxF = 5, ### maximum allowed Fbar
                     dupl_trgt = FALSE,
                     proc_res = NULL, ### process error noise,
                     migration = NULL, ### FLQuant with migration factor(s)
                     disc_survival = 0, ### discard survival proportion
                     dd_M = FALSE, ### density-dependent M
                     dd_M_relation = NULL, ### density-dependent M
                     dd_M_fun = calculate_ddM,
                     ...) {
  
  ### avoid the issue that the catch is higher than the targeted catch
  ### can happen due to bug in FLash if >1 iteration provided
  ### sometimes, FLash struggles to get estimates and then uses F estimate from
  ### previous iteration
  ### workaround: target same value several times and force FLash to try again
  if (isTRUE(dupl_trgt)) {
    
    ### duplicate target
    ctrl@target <- rbind(ctrl@target, ctrl@target, ctrl@target)
    ### replace catch in second row with landings
    ctrl@target$quantity[1] <- "landings"
    ctrl@target$quantity[3] <- "catch"
    
    ### extract target values
    val_temp <- ctrl@trgtArray[, "val", ]
    
    ### extend trgtArray
    ### extract dim and dimnames
    dim_temp <- dim(ctrl@trgtArray)
    dimnames_temp <- dimnames(ctrl@trgtArray)
    ### duplicate years
    dim_temp[["year"]] <- dim_temp[["year"]] * 3
    dimnames_temp$year <- rep(dimnames_temp$year, 3)
    
    ### create new empty array
    trgtArray <- array(data = NA, dim = dim_temp, dimnames = dimnames_temp)
    
    ### fill with values
    ### first as target
    trgtArray[1, "val", ] <- val_temp
    ### then again, but as max
    trgtArray[2, "max", ] <- val_temp
    ### min F
    trgtArray[3, "max", ] <- val_temp
    
    ### insert into ctrl object
    ctrl@trgtArray <- trgtArray
  }
  
  ### calculate density-dependent natural mortality if required
  if (isTRUE(dd_M)) {
    
    ### overwrite M in the target year before projecting forward
    m(stk)[, ac(ctrl@target[, "year"])] <- 
      dd_M_fun(stk, ctrl@target[, "year"], relation = dd_M_relation)
    
  }
  
  ### migration
  if (!is.null(migration)) {
    ### find year to adapt for migration
    yr_target <- ctrl@target[, "year"]
    yr_migr <- yr_target - 1 ### adapt year before target
    migr_factor <- FLQuant(NA, dimnames = list(age = dimnames(stk)$age, 
                                               year = yr_migr, 
                                               iter = dimnames(stk)$iter))
    migr_factor[] <- sr@ssb[, ac(yr_target)]
    ### update stock numbers
    stock.n_bckp <- stock.n(stk)[, ac(yr_migr)]
    stock.n(stk)[, ac(yr_migr)] <- stock.n(stk)[, ac(yr_migr)] * migr_factor
  }
  
  ### project forward with FLash::fwd
  if (!isTRUE(disc_survival > 0)) {
    
    stk[] <- fwd(object = stk, control = ctrl, sr = sr, 
                 sr.residuals = sr.residuals, 
                 sr.residuals.mult = sr.residuals.mult,
                 maxF = maxF)
  
  } else {
    ### account for discard survival
    stk_bckp <- stk
    stk_fc1 <- stk ### backup stock
    ### first projection with all discards
    stk_fc1[] <- fwd(object = stk_fc1, control = ctrl, sr = sr, 
                     sr.residuals = sr.residuals, 
                     sr.residuals.mult = sr.residuals.mult,
                     maxF = maxF)
    yr_fc <- ctrl@target$year
    landings <- landings(stk_fc1[, ac(yr_fc)])
    ### get discards after accounting for survival
    discards <- quantSums(discards.n(stk_fc1)[, ac(yr_fc)] * 
                            (1 - disc_survival) * 
                            discards.wt(stk_fc1)[, ac(yr_fc)])
    ### total catch (after accounting for discards survival)
    catch <- landings + discards
    ### ctrl object with new target
    ctrl_dead <- ctrl
    ctrl_dead@trgtArray[, "val", ] <- c(catch)
    ### prepare second forecast accounting for discard survival
    stk_fc2 <- stk
    ### update discard rate
    discards.n(stk_fc2)[, ac(yr_fc)] <- discards.n(stk_fc2)[, ac(yr_fc)] *
      (1 - disc_survival)
    ### re-standardise landings rate
    landings.n(stk_fc2)[, ac(yr_fc)] <- landings.n(stk_fc2)[, ac(yr_fc)] / 
      (landings.n(stk_fc2)[, ac(yr_fc)] + discards.n(stk_fc2)[, ac(yr_fc)])
    ### update catch weight
    catch.wt(stk_fc2)[, ac(yr_fc)] <- 
      ((landings.wt(stk_fc2)[, ac(yr_fc)] * landings.n(stk_fc2)[, ac(yr_fc)]) +
      (discards.wt(stk_fc2)[, ac(yr_fc)] * discards.n(stk_fc2)[, ac(yr_fc)])) /
      (landings.n(stk_fc2)[, ac(yr_fc)] + discards.n(stk_fc2)[, ac(yr_fc)])
    ### run second forecast
    stk_fc2[] <- fwd(object = stk_fc2, control = ctrl_dead, sr = sr, 
                     sr.residuals = sr.residuals, 
                     sr.residuals.mult = sr.residuals.mult,
                     maxF = maxF)
    ### insert projected values into OM stk
    ### catch includes all discards
    landings.n(stk)[, ac(yr_fc)] <- landings.n(stk_fc1)[, ac(yr_fc)]
    discards.n(stk)[, ac(yr_fc)] <- discards.n(stk_fc1)[, ac(yr_fc)]
    catch.n(stk)[, ac(yr_fc)] <- catch.n(stk_fc1)[, ac(yr_fc)]
    landings(stk)[, ac(yr_fc)] <- landings(stk_fc1)[, ac(yr_fc)]
    discards(stk)[, ac(yr_fc)] <- discards(stk_fc1)[, ac(yr_fc)]
    catch(stk)[, ac(yr_fc)] <- catch(stk_fc1)[, ac(yr_fc)]
    ### stock and harvest include only dead discards
    stock.n(stk)[, ac(yr_fc)] <- stock.n(stk_fc2)[, ac(yr_fc)]
    stock(stk)[, ac(yr_fc)] <- stock(stk_fc2)[, ac(yr_fc)]
    harvest(stk)[, ac(yr_fc)] <- harvest(stk_fc2)[, ac(yr_fc)]
    
  }
  
  ### if migration occurs, remove numbers for year before forecast
  ### (but they are kept for the forecast year)
  if (!is.null(migration)) {
    stock.n(stk)[, ac(yr_migr)] <- stock.n_bckp
  }
  
  ### add process error noise, if supplied
  if (!is.null(proc_res)) {
    
    ### projected years
    yrs_new <- seq(from = ctrl@target[, "year"], to = range(stk)[["maxyear"]])
    
    ### workaround to get residuals
    ### they are saved in the "fitted" slot of sr...
    if (!isTRUE(proc_res == "fitted")) {
      
      stop("survival process error inacessible")
      
    } else {
      
      ### implement process error
      stock.n(stk)[, ac(yrs_new)] <- stock.n(stk)[, ac(yrs_new)] *
        fitted(sr)[, ac(yrs_new)]
      ### update stock biomass
      stock(stk) <- computeStock(stk)
      
    }
    
  }
  
  ### return stock
  return(list(object = stk))
  
}


### ------------------------------------------------------------------------ ###
### iter subset  ####
### ------------------------------------------------------------------------ ###

iter_attr <- function(object, iters, subset_attributes = TRUE) {
  
  ### subset object to iter
  res <- FLCore::iter(object, iters)
  
  if (isTRUE(subset_attributes)) {
    
    ### get default attributes of object class
    attr_def <- names(attributes(new(Class = class(object))))
    
    ### get additional attributes
    attr_new <- setdiff(names(attributes(object)), attr_def)
    
    ### subset attributes
    for (attr_i in attr_new) {
      attr(res, attr_i) <- FLCore::iter(attr(res, attr_i), iters)
    }
    
  }
  
  return(res)
  
}

### ------------------------------------------------------------------------ ###
### estimtate steepness based on l50/linf ratio ####
### according to Wiff et al. 2018
### ------------------------------------------------------------------------ ###
h_Wiff <- function(l50, linf) {
  l50linf <- l50/linf
  ### linear model
  lin <- 2.706 - 3.698*l50linf
  ### logit
  h <- (0.2 + exp(lin)) / (1 + exp(lin))
  return(h)
}

### ------------------------------------------------------------------------ ###
### mean length in catch ####
### ------------------------------------------------------------------------ ###
lmean <- function(stk, params) {
  
  ### calculate length from age with a & b
  weights <- c(catch.wt(stk)[, 1,,,, 1])
  lengths <- (weights / c(params["a"]))^(1 / c(params["b"]))
  catch.n <- catch.n(stk)
  dimnames(catch.n)$age <- lengths
  ### subset to lengths > Lc
  catch.n <- catch.n[lengths > c(params["Lc"]),]
  
  ### calculate mean length
  lmean <- apply(X = catch.n, MARGIN = c(2, 6), FUN = function(x) {
    ### calculate
    res <- weighted.mean(x = an(dimnames(x)$age), 
                         w = ifelse(is.na(x), 0, x), na.rm = TRUE)
    ### check if result obtained
    ### if all catch at all lengths = 0, return 0 as mean length
    # if (is.nan(res)) {
    #   if (all(ifelse(is.na(x), 0, x) == 0)) {
    #     res[] <- 0
    #   }
    # }
    return(res)
  })
  return(lmean)
}

### ------------------------------------------------------------------------ ###
### length at first capture ####
### ------------------------------------------------------------------------ ###
calc_lc <- function(stk, a, b) {
  ### find position in age vector
  Ac <- apply(catch.n(stk), MARGIN = c(2, 6), function(x) {
    head(which(x >= (max(x, na.rm = TRUE)/2)), 1)
  })
  Ac <- an(median(Ac))
  ### calculate lengths
  weights <- c(catch.wt(stk)[, 1,,,, 1])
  lengths <- (weights / a)^(1 / b)
  ### length at Ac
  Lc <- floor(lengths[Ac]*10)/10
  return(Lc)
}

### ------------------------------------------------------------------------ ###
### inter-annual variability ####
### ------------------------------------------------------------------------ ###
#' calculate inter-annual variability of FLQuant
#'
#' This function calculates survey indices from the numbers at age of an 
#' FLStock object
#'
#' @param object Object of class \linkS4class{FLQuant} with values.
#' @param period Select every n-th year, e.g. biennial (optional).
#' @param from,to Optional year range for analysis.
#' @param summary_per_iter Function for summarising per iter. Defaults to mean.
#' @param summary Function for summarising over iter. Defaults to mean.
#' @return An object of class \code{FLQuant} with inter-annual variability.
#'
#' @export
#' 
setGeneric("iav", function(object, period, from, to, summary_per_iter, 
                           summary_year, summary_all) {
  standardGeneric("iav")
})

### object = FLQuant
#' @rdname iav
setMethod(f = "iav",
          signature = signature(object = "FLQuant"),
          definition = function(object, 
                                period, ### periodicity, e.g. use every 2nd value 
                                from, to,### year range
                                summary_per_iter, ### summarise values per iteration
                                summary_year,
                                summary_all) {
            
            ### subset years
            if (!missing(from)) object <- FLCore::window(object, start = from)
            if (!missing(to)) object <- FLCore::window(object, end = from)
            
            ### get years in object
            yrs <- dimnames(object)$year
            
            ### select every n-th value, if requested
            if (!missing(period)) {
              yrs <- yrs[seq(from = 1, to = length(yrs), by = period)]
            }
            
            ### reference years
            yrs_ref <- yrs[-length(yrs)]
            ### years to compare
            yrs_comp <- yrs[-1]
            
            ### calculate variation (absolute values, ignore pos/neg)
            res <- abs(1 - object[, yrs_comp] / object[, yrs_ref])
            
            ### replace Inf with NA (compared to 0 catch)
            res <- ifelse(is.finite(res), res, NA)
            
            ### summarise per iteration
            if (!missing(summary_per_iter)) {
              res <- apply(res, 6, summary_per_iter, na.rm = TRUE)
            }
            
            ### summarise per year
            if (!missing(summary_year)) {
              res <- apply(res, 1:5, summary_year, na.rm = TRUE)
            }
            
            ### summarise over everything
            if (!missing(summary_all)) {
              
              res <- summary_all(c(res), na.rm = TRUE)
              
            }
            
            return(res)
            
          })


### ------------------------------------------------------------------------ ###
### "correct" collapses ####
### ------------------------------------------------------------------------ ###

collapse_correction <- function(stk, quants = c("catch", "ssb", "fbar", "rec"),
                                threshold = 1, yrs = NULL) {
  
  if (is.null(yrs)) yrs <- dimnames(stk)$year
  
  names(quants) <- quants
  qnt_list <- lapply(quants, function(x) get(x)(stk))
  qnt_list <- lapply(qnt_list, function(x) x[, ac(yrs)])
  
  n_yrs <- dim(qnt_list[[1]])[2]
  n_its <- dim(qnt_list[[1]])[6]
  
  ### find collapses
  cd <- sapply(seq(n_its), function(x) {
    min_yr <- min(which(qnt_list$ssb[,,,,, x] < 1))
    if (is.finite(min_yr)) {
      all_yrs <- min_yr:n_yrs
    } else {
      all_yrs <- NA
    }
    all_yrs + (x - 1)*n_yrs
  })
  cd <- unlist(cd)
  cd <- cd[which(!is.na(cd))]
  ### remove values
  qnt_list <- lapply(qnt_list, function(x) {
    x@.Data[cd] <- 0
    return(x)
  })
  return(qnt_list)
}




