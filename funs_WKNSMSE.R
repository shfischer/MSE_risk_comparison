
### ------------------------------------------------------------------------ ###
### stock assessment: SAM wrapper ####
### ------------------------------------------------------------------------ ###
### wrapper for calling SAM
### this makes used of the SAM wrapper in FLfse,
### which in turn calls SAM from the package stockassessment
SAM_wrapper <- function(stk, idx, tracking,
                        args, ### contains ay (assessment year)
                        forecast = FALSE,
                        fwd_trgt = "fsq", ### what to target in forecast
                        fwd_yrs = 1, ### number of years to add
                        fwd_yrs_average = -3:0, ### years used for averages
                        fwd_yrs_rec_start = NULL, ### recruitment 
                        fwd_yrs_sel = NULL, ### selectivity
                        fwd_yrs_lf_remove = -2:-1,
                        fwd_splitLD = TRUE,
                        fwd_rec, fwd_rec_yrs, fwd_yrs_bio, ### for FLasher
                        fwd_yrs_mat, fwd_yrs_m,            ### forecast
                        fwd_disc_zero = TRUE, ### for FLasher only
                        parallel = FALSE,
                        conf = NULL, ### SAM configuration
                        par_ini = NULL, ### initial parameters
                        track_ini = FALSE, ### store ini for next year
                        ...){
  
  ### get additional arguments
  args <- c(args, list(...))
  
  ### get current (assessment) year
  ay <- args$ay
  
  ### check if initial parameter values for SAM exist from last year's fit
  ### and reuse if possible
  ### this overrides generic initial parameters, if they exist in par_ini
  ### (they are only used in first simulation year)
  if (isTRUE(track_ini) & !is.null(attr(tracking[[1]]@units, "par_ini"))) {
    
    par_ini <- attr(tracking[[1]]@units, "par_ini")
    
  }
  
  ### fit SAM to provided data
  fit <- FLR_SAM(stk = stk, idx = idx, conf = conf, par_ini = par_ini,
                 DoParallel = parallel, ...)
  
  ### store parameter values and provide them as initial values next year
  ### store in tracking object, this is the only object which does not get 
  ### overriden next year
  ### weird workaround: save as attribute of "unit" slot of tracking,
  ###                   otherwise the attributes will be lost later...
  if (isTRUE(track_ini)) {
    
    attr(tracking[[1]]@units, "par_ini") <- getpars(fit)
    
  }
  
  ### convert into FLStock
  stk0 <- SAM2FLStock(object = fit, stk = stk) 
  
  ### perform forecast to get SSB ay+1
  if (isTRUE(forecast)) {
    
    ### check how to do forecast
    ### currently, can only do F status quo
    if (!all(fwd_trgt %in% c("fsq","TAC"))) {
      stop("only fsq and TAC supported in forecast")
    }
    
    ### years for average values
    ave.years <- range(stk0)[["maxyear"]] + fwd_yrs_average
    ### years for sampling of recruitment years
    rec.years <- seq(from = fwd_yrs_rec_start, to = range(stk0)[["maxyear"]] - 1)
    ### years where selectivity is not used for mean in forecast
    overwriteSelYears <- NULL
    if (!is.null(fwd_yrs_sel)) 
      overwriteSelYears <- range(stk0)[["maxyear"]] + fwd_yrs_sel
    
    lst_yr <- range(stk0)[["maxyear"]]
    
    ### extend stk0
    stk0 <- window(stk0, end = range(stk0)[["maxyear"]] + fwd_yrs)
    
    ### modify fwd_yrs in case last data year does not include catch
    if (all(is.na(catch(stk0)[, ac(lst_yr)]))) fwd_yrs <- fwd_yrs + 1
    
    ### forecast years
    yrs <- seq(to = dims(stk0)$maxyear, length.out = fwd_yrs)
    
    ### coerce fit into list if only 1 iteration
    if (is(fit, "sam")) {
      fit <- list(fit)
      class(fit) <- "sam_list"
    }
    
    ### template for forecast targets
    fscale <- ifelse(fwd_trgt == "fsq", 1, NA)
    catchval <- ifelse(fwd_trgt == "TAC", -1, NA)
    ### recycle target if neccessary
    if (fwd_yrs > length(fwd_trgt)) {
      fscale <- c(fscale, tail(fscale, 1))
      catchval <- c(catchval, tail(catchval, 1))
    }
    ### get recent TAC
    if (args$iy == ay) {
      ### in first year of simulation, use value from OM saved earlier in ay
      TAC_last <- tracking[[1]]["C.om", ac(ay)]
    } else {
      ### in following years, use TAC advised the year before
      TAC_last <- tracking[[1]]["isys", ac(ay - 1)]
    }
    
    ### do forecast for all iterations
    fc <- foreach(fit_i = fit, iter_i = seq_along(fit),
                  .errorhandling = "pass") %do% {
                    
      ### overwrite landing fraction with last year, if requested
      if (!is.null(fwd_yrs_lf_remove)) {
        ### index for years to remove/overwrite
        idx_remove <- nrow(fit_i$data$landFrac) + args$fwd_yrs_lf_remove
        ### overwrite
        fit_i$data$landFrac[idx_remove, ] <- fit_i$data$landFrac[rep(nrow(fit_i$data$landFrac), length(idx_remove)), ]
      }
      
      ### define forecast targets for current iteration
      fscale_i <- fscale
      ### load TAC as catch target
      catchval_i <- ifelse(catchval == -1, c(TAC_last[,,,,, iter_i]), catchval)
      
      ### arguments for forecast
      fc_args <- list(fit = fit_i, fscale = fscale_i, catchval = catchval_i,
                      ave.years = ave.years, rec.years = rec.years,
                      overwriteSelYears = overwriteSelYears, 
                      splitLD = fwd_splitLD)
      ### for compatibility of stockassessment's commit a882a11 and later:
      if ("savesim" %in% names(formals(base::args(stockassessment::forecast)))) {
        fc_args$savesim <- TRUE
      }
      ### run forecast
      fc_i <- do.call(stockassessment::forecast, fc_args)
      
      ### get numbers at age for all forecast years
      numbers <- lapply(seq(fwd_yrs), function(x) {
        ### index for numbers at age
        idx <- seq(length(fit_i$conf$keyLogFsta[1, ]))
        ### get simulated numbers
        n <- exp(fc_i[[x]]$sim[, idx])
        ### median
        apply(n, 2, median)
      })
      numbers <- do.call(cbind, numbers)
      
      return(list(stock.n = numbers))
      
    }
    ### if forecast failed for a iteration, the list element will for this
    ### iteration will be an error message
    
    ### get numbers
    fc_stock.n <- lapply(fc, "[[", "stock.n")
    
    ### insert stock numbers
    for (iter_i in seq(dim(stk0)[6])) {
      ### do not insert numbers if forecast failed
      if (!is.numeric(fc_stock.n[[iter_i]])) next()
      stock.n(stk0)[, ac(yrs),,,, iter_i] <- fc_stock.n[[iter_i]]
    }
    
    ### extend stock characteristics required for calculation of SSB, 
    ### weights, etc.
    
    ### find years to fill (do not fill years, if there is already data inside)
    yrs_fill <- setdiff(yrs, lst_yr)
    
    stock.wt(stk0)[, ac(yrs_fill)] <- yearMeans(stock.wt(stk0)[, ac(ave.years)])
    m(stk0)[, ac(yrs_fill)] <- yearMeans(m(stk0)[, ac(ave.years)])
    mat(stk0)[, ac(yrs_fill)] <- yearMeans(mat(stk0)[, ac(ave.years)])
    
    harvest.spwn(stk0)[, ac(yrs_fill)] <- yearMeans(harvest.spwn(stk0)[, 
                                                                       ac(ave.years)])
    m.spwn(stk0)[, ac(yrs_fill)] <- yearMeans(m.spwn(stk0)[, ac(ave.years)])
    harvest(stk0)[, ac(yrs_fill)] <- yearMeans(harvest(stk0)[, ac(ave.years)])
    
    ### PLEASE NOTE:
    ### SSB value slightly different from SSB value generated from SAM:
    ### SAM calculates SSB per simulation and gives median
    ### here: calculate median of numbers at age and calculate SSB from
    ###       median numbers
    
  } else if (identical(forecast, "FLasher")) {
    ### deterministic forecast with FLasher
    yr_max <- range(stk)[["maxyear"]]
    yrs_fill <- seq(yr_max, to = yr_max + fwd_yrs)
    stk0 <- stf(stk0, nyears = fwd_yrs)
    ### fill with data
    stock.wt(stk0)[, ac(yrs_fill)] <- 
      yearMeans(stock.wt(stk0)[, ac(yr_max + fwd_yrs_bio)])
    landings.wt(stk0)[, ac(yrs_fill)] <- 
      yearMeans(landings.wt(stk0)[, ac(yr_max + fwd_yrs_bio)])
    discards.wt(stk0)[, ac(yrs_fill)] <- 
      yearMeans(discards.wt(stk0)[, ac(yr_max + fwd_yrs_bio)])
    harvest(stk0)[, ac(yrs_fill)] <- 
      yearMeans(harvest(stk0)[, ac(yr_max + fwd_yrs_sel)])
    mat(stk0)[, ac(yrs_fill)] <- 
      yearMeans(mat(stk0)[, ac(yr_max + fwd_yrs_mat)])
    m(stk0)[, ac(yrs_fill)] <- 
      yearMeans(m(stk0)[, ac(yr_max + fwd_yrs_m)])
    if (isTRUE(fwd_disc_zero)) discards.n(stk0)[] <- 0
    ### recruitment model - constant value
    sr_stf <- FLSR(model = "geomean", params = FLPar(1))
    residuals(sr_stf) <- FLQuant(1, dimnames = list(year = yrs_fill,
                                                    iter = seq(dim(stk0)[6])))
    if (identical(tail(fwd_rec, 1), "geomean_weighted")) {
      ### go through iterations and calculate recruitment
      R_insert <- foreach(fit_i = fit, .combine = c) %do% {
        . <- try({
          ### get index for recruitment values and subset to year to be used
          idx_R <- which(names(fit_i$sdrep$value) == "logR")
          idx_R_length <- length(idx_R)
          idx_R <- idx_R[idx_R_length + fwd_rec_yrs]
          ### R: weighted geometric mean
          R_fwd <- exp(weighted.mean(x = fit_i$sdrep$value[idx_R],
                                     w = (1/fit_i$sdrep$sd[idx_R])^2))
        })
        if (is(., "try-error")) return(0)
        return(R_fwd)
      }
      residuals(sr_stf)[] <- R_insert
    }
    if (identical(fwd_rec[1], "estimate")) {
      residuals(sr_stf)[, ac(yrs_fill[1])] <- rec(stk0)[, ac(yrs_fill[1])]
    }
    ### short-term forecast
    ftrgt <- fbar(stk0)[, ac(yrs_fill)]
    ctrl_stf <- fwdControl(ftrgt, quant = "f")
    stk0 <- fwd(stk0, control = ctrl_stf, sr = sr_stf, 
                deviances = residuals(sr_stf))
  }
  
  ### save convergence for all iterations
  tracking[[1]]["conv.est", ac(ay)] <- sapply(fit, function(x) {
    if (isTRUE(is(x, "sam"))) {
      return(x$opt$convergence)
    } else {
      return(2)
    }
  })
  ### add perceived F and SSB
  ### done in mp()
  #tracking["F.est", ac(ay)] <- fbar(stk0)[, ac(ay - 1)]
  #tracking["B.est", ac(ay)] <- tail(ssb(stk0))
  
  ### save model fit (list) as attribute in stk0
  attr(stk0, "fit") <- fit
  
  ### return assessed stock, tracking & model fit 
  ### (model fit required later for TAC calculation)
  return(list(stk = stk0, tracking = tracking))
  
}

### ------------------------------------------------------------------------ ###
### phcr: parameterize HCR ####
### ------------------------------------------------------------------------ ###
phcr_WKNSMSE <- function(Btrigger = NULL, Ftrgt = NULL, Bpa = NULL, Fpa = NULL,
                         Blim = NULL, tracking, ...) {
  
  ### coerce existing values into FLPar
  hcrpars <- list("Btrigger" = Btrigger, "Ftrgt" = Ftrgt, "Bpa" = Bpa, 
                  "Fpa" = Fpa, "Blim" = Blim)
  hcrpars <- hcrpars[!sapply(hcrpars, is.null)]
  hcrpars <- do.call(FLPar, hcrpars)
  
  ### if more iterations provided than neccessary, subset
  if (dims(hcrpars)$iter > dims(tracking[[1]])$iter) {
    hcrpars <- hcrpars[, dimnames(tracking[[1]])$iter]
  }
  
  ### return as list
  ### keep tracking unchanged
  return(list(hcrpars = hcrpars, tracking = tracking))
  
}

### ------------------------------------------------------------------------ ###
### hcr: ICES advice rule and modificatinons ####
### ------------------------------------------------------------------------ ###

### target Ftrgt
### A: if SSB < Btrigger, target reduced: Ftrgt * SSB/Btrigger
### B & C: different behaviour if SSB < Blim
hcr_WKNSME <- function(stk, args, hcrpars, tracking, 
                       option = "A", ### WKNSMSE options
                       ...) {
  
  ### target F
  ### reduce F if SSB below Btrigger
  
  ### get current (assessment) year
  ay <- args$ay
  ### number of iterations
  it <- dim(stk)[6]
  
  ### get reference values and extend iteration dimension
  Ftrgt <- propagate(FLPar(hcrpars["Ftrgt"]), it)
  Btrigger <- propagate(FLPar(hcrpars["Btrigger"]), it)
  
  ### get Blim
  ### if Blim not specified, set to zero
  if ("Blim" %in% dimnames(hcrpars)$params) {
    
    Blim <- propagate(FLPar(hcrpars["Blim"]), it)
    
  } else {
    
    Blim <- propagate(FLPar(0), it)
    
  } 
  
  ### SSB status (last year only) relative to Btrigger
  status_Btrigger <- tail(ssb(stk), 1) / Btrigger 
  ### SSB status (last year only) relative to Blim
  status_Blim <- tail(ssb(stk), 1) / Blim
  
  ### positions (iterations) where SSB is below Btrigger
  pos_Btrigger <- which(status_Btrigger < 1)
  ### below Blim
  pos_Blim <- which(status_Blim < 1)
  
  ### default ICES HCR (option A):
  ### if SSB<Btrigger => F = Ftrgt * SSB/Btrigger
  if (option == "A") {
    
    mult <- ifelse(status_Btrigger < 1, status_Btrigger, 1)
    
  } else if (option == "B") {
    ### option B:
    ### if SSB<Btrigger => F = Ftrgt * SSB/Btrigger
    ### if SSB<Blim => F = Ftrgt * 0.25
    ### SSB < Btrigger
    mult <- ifelse(status_Btrigger < 1, status_Btrigger, 1)
    ### set to 0.25 if SSB < Blim
    if (length(pos_Blim) > 0) {
      mult[,,,,, pos_Blim] <- 0.25
    }
    
  } else if (option == "C") {
    ### option C:
    ### if SSB<Btrigger => F = Ftrgt * SSB/Btrigger
    ### if SSB<Blim => F = max(Ftrgt * 0.25, Ftrgt * SSB/Btrigger)
    ### SSB < Btrigger
    mult <- ifelse(status_Btrigger < 1, status_Btrigger, 1)
    ### if SSB < Blim, use maximum of SSB/Btrigger or 0.25
    ### i.e. limit ratio to 0.25
    if (length(pos_Blim) > 0) {
      mult[,,,,, pos_Blim] <- c(ifelse(mult[,,,,, pos_Blim] < 0.25, 0.25, 
                                       mult[,,,,, pos_Blim]))
    }
    
  } else {
    
    ### unknown options
    mult <- 1
    
  }
  
  ### new target
  Ftrgt <- Ftrgt * mult
  
  ### create ctrl object
  advice <- FLQuant(c(Ftrgt), 
                    dimnames = list(year = ay + 1,
                                    iter = seq(it)))
  ctrl <- fwdControl(target = advice, quant = "f")
  
  ### save in tracking
  ### done later
  # tracking["advice", ac(ay)] <- ctrl@trgtArray[ac(ay + 1), "val", ]
  
  return(list(ctrl = ctrl, tracking = tracking))
  
}



### ------------------------------------------------------------------------ ###
### management implementation: short term forecast with SAM ####
### ------------------------------------------------------------------------ ###
### short term forecast with SAM
### including TAC constraint
is_WKNSMSE <- function(stk, tracking, ctrl,
                       args, ### contains ay (assessment year)
                       TAC_constraint = FALSE,
                       upper = Inf, lower = -Inf, Btrigger_cond = FALSE,
                       ### short term forecast
                       forecast = TRUE, ### SAM forecast
                       fwd_trgt = c("fsq", "hcr"), ### target in forecast
                       fwd_yrs = 2, ### number of years to add
                       fwd_yrs_average = -3:0, ### years used for averages
                       fwd_yrs_rec_start = NULL, ### recruitment 
                       fwd_yrs_sel = NULL, ### selectivity
                       fwd_yrs_lf_remove = -2:-1,
                       fwd_splitLD = TRUE, 
                       ### banking and borrowing
                       BB = FALSE, ### banking and borrowing
                       ### check stock status before applying BB
                       BB_check_hcr = FALSE, ### check status before forecast
                       BB_check_fc = FALSE, ### check status after forecast
                       BB_rho, ### definition of BB
                       ### reference points
                       hcrpars = list(),
                       ...) {
  
  ### get current (assessment) year
  ay <- args$ay
  ### number of iterations
  it <- dim(stk)[6]
  
  ### retrieve SAM model fit (list)
  fit <- attr(stk, "fit")
  
  ### check class of model fit(s)
  if (!class(fit) %in% c("sam", "sam_list")) 
    stop("attr(stk0, \"fit\") has to be class sam or sam_list")
  
  ### if single fit, turn into list
  if (is(fit, "sam")) fit <- list(fit)
  
  ### if conditional banking & borrowing applied, extend forecast for one more
  ### year to check SSB in year after advice year
  ### for this additional forecast assume Fsq as target
  ### i.e. target F from HCR twice, in analogy to intermediate year assumption
  if (isTRUE(BB) & isTRUE(BB_check_fc)) {
    
    ### duplicate last target value
    fwd_trgt <- c(fwd_trgt, tail(fwd_trgt, 1))
    
  }
  
  ### get recent TAC
  if (args$iy == ay) {
    ### in first year of simulation, use value from OM saved earlier in ay
    TAC_last <- tracking[[1]]["C.om", ac(ay)]
  } else {
    ### in following years, use TAC advised the year before
    TAC_last <- tracking[[1]]["isys", ac(ay - 1)]
  }
  
  ### forecast with SAM
  if (isTRUE(forecast)) {
    ### go through all model fits
    fc <- foreach(fit_i = fit, iter_i = seq_along(fit), 
                  .errorhandling = "pass") %do% {
                    
      ### overwrite landing fraction with last year, if requested
      if (!is.null(fwd_yrs_lf_remove)) {
        ### index for years to remove/overwrite
        idx_remove <- nrow(fit_i$data$landFrac) + fwd_yrs_lf_remove
        ### overwrite
        fit_i$data$landFrac[idx_remove, ] <- fit_i$data$landFrac[rep(nrow(fit_i$data$landFrac), length(idx_remove)), ]
      }
      
      ### check how to do forecast
      ### can handle F status quo, F target from ctrl object and TAC
      ### scaled F
      fscale <- ifelse(fwd_trgt == "fsq", 1, NA)
      ### target F values
      fval <- ifelse(fwd_trgt == "hcr", ctrl@iters[, "value", iter_i], NA)
      ### target catch values
      catchval <- ifelse(fwd_trgt == "TAC", c(TAC_last[,,,,, iter_i]), NA)
      
      ### years for average values
      ave.years <- max(fit_i$data$years) + fwd_yrs_average
      ### years for sampling of recruitment years
      if (is.null(fwd_yrs_rec_start)) {
        rec.years <- fit_i$data$years ### use all years, if not defined
      } else {
        rec.years <- seq(from = fwd_yrs_rec_start, max(fit_i$data$years))
      }
      
      ### years where selectivity is not used for mean in forecast
      overwriteSelYears <- NULL
      if (!is.null(fwd_yrs_sel)) 
        overwriteSelYears <- max(fit_i$data$years) + fwd_yrs_sel
      
      ### arguments for forecast
      fc_args <- list(fit = fit_i, fscale = fscale, fval = fval, 
                      catchval = catchval,
                      ave.years = ave.years, rec.years = rec.years,
                      overwriteSelYears = overwriteSelYears, 
                      splitLD = fwd_splitLD)
      ### for compatibility of stockassessment's commit a882a11 and later:
      if ("savesim" %in% names(formals(base::args(stockassessment::forecast)))) {
        fc_args$savesim <- TRUE
      }
      ### run forecast
      fc_i <- do.call(stockassessment::forecast, fc_args)
      
      ### return forecast table
      return(attr(fc_i, "tab"))
      
    }
    ### if forecast fails, error message returned
    ### replace error message with NA
    ### extract catch target
    catch_target <- sapply(fc, function(x) {
      if (is(x, "error")) {
        return(NA)
      } else {
        return(x[ac(ay + 1), "catch:median"])
      }
    })
  ### forecast with FLasher
  } else if (identical(forecast, "FLasher")) {
    ### define recruitment model
    sr_stf <- FLSR(model = "geomean", params = FLPar(1))
    ### keep recruitment values defined in estimator
    residuals(sr_stf) <- rec(stk)
    ### forecast observed stock
    stk <- fwd(stk, control = ctrl, sr = sr_stf, 
               deviances = residuals(sr_stf))
    ### extract catch
    catch_target <- c(catch(stk)[, ac(ay + 1)])
  }
  
  ### get reference points
  hcrpars <- hcrpars[!sapply(hcrpars, is.null)]
  hcrpars <- do.call(FLPar, hcrpars)
  ### if more iterations provided than neccessary, subset
  if (dims(hcrpars)$iter > dims(tracking[[1]])$iter) {
    hcrpars <- hcrpars[, dimnames(tracking[[1]])$iter]
  } else if (isTRUE(dim(hcrpars)[2] < it)) {
    hcrpars <- propagate(hcrpars, it)
  }
  
  ### ---------------------------------------------------------------------- ###
  ### TAC constraint ####
  if (isTRUE(TAC_constraint)) {
    
    ### target year
    yr_target <- ctrl@target$year
    
    ### get previous target from tracking
    if (args$iy == ay) {
      catch_prev <- TAC_last
    } else {
      catch_prev <- tracking[[1]]["isys", ac(yr_target - 1), drop = TRUE]
    }
    
    ### change in advice, % of last advice
    change <- (catch_target / catch_prev) * 100
    
    ### limit changes
    changes_new <- change
    ### upper limit
    changes_new <- ifelse(changes_new > upper, upper, changes_new) 
    ### lower limit
    changes_new <- ifelse(changes_new < lower, lower, changes_new) 
    
    ### find positions which exceed limits
    pos <- which(changes_new != change)
    
    ### conditional constraint based on SSB>=Btrigger?
    if (isTRUE(Btrigger_cond)) {
      
      ### get SSB in TAC year from forecast
      SSB_TACyr <- sapply(fc, function(x) {
        if (is(x, "error")) {
          return(NA)
        } else {
          return(x[ac(ay + 1), "ssb:median"])
        }
      })
      
      ### iterations where SSB is at or above Btrigger at start of TAC year
      pos_Btrigger <- which(SSB_TACyr >= c(hcrpars["Btrigger"]))
      ### only apply TAC constraint if both
      ### - TAC change exceeds limit
      ### - stock at or above Btrigger
      pos <- intersect(pos, pos_Btrigger)
      
    }
    
    ### modify advice
    catch_target[pos] <- catch_prev[pos] * changes_new[pos]/100
    
  }
  
  ### ---------------------------------------------------------------------- ###
  ### banking and borrowing ####
  
  if (isTRUE(BB)) {
    
    ### get current rho
    BB_rho_i <- tail(rep(BB_rho, length.out = (ay - args$y0)), 1)
    
    ### get catch borrowed last year
    BB_return <- tracking[[1]]["BB_borrow", ac(ay - 1)]
    ### assume nothing borrowed if NA
    BB_return <- ifelse(!is.na(BB_return), BB_return, 0)
    
    ### get catch banked last year
    BB_bank_use <- tracking[[1]]["BB_bank", ac(ay - 1)]
    BB_bank_use <- ifelse(!is.na(BB_bank_use), BB_bank_use, 0)
    
    ### bank for next year
    if (BB_rho_i < 0) {
      BB_bank <- catch_target * abs(BB_rho_i)
    } else {
      BB_bank <- rep(0, it)
    }
    
    ### borrow from next year
    if (BB_rho_i > 0) {
      BB_borrow <- catch_target * abs(BB_rho_i)
    } else {
      BB_borrow <- rep(0, it)
    }
    
    ### conditional banking and borrowing?
    ### first condition: for HCR option A (stability option D2)
    ### apply BB only if HCR option A1 (not A2) is applied
    ### i.e. stop BB if SSB is below Btrigger in TAC year
    
    ### find iterations where SSB is below Btriggerin TAC year
    if (isTRUE(BB_check_hcr)) {
      pos_hcr <- which(c(tail(ssb(stk), 1)) < c(hcrpars["Btrigger"]))
    } else {
      pos_hcr <- integer(0)
    }
    
    ### second condition: for HCR options B & C (stability option E2)
    ### stop BB if either
    ### - if SSB is below Bpa AND F above Fpa in TAC year
    ### - if SSB is below Bpa in TAC year and year after
    ### if TAC restricted by TAC constraint, additional forecast required
    ### to estimate stock status when only TAC is fished
    pos_fc <- integer(0)
    if (any(c(BB_bank, BB_borrow) > 0) & isTRUE(BB_check_fc)) {
      
      ### if TAC constraint activated, do a forecast
      ### check if TAC constraint used
      if (isTRUE(TAC_constraint)) {
        ### check where TAC constraint is implemented
        pos_constr <- which(changes_new != change)
        if (isTRUE(length(pos_constr) > 0)) {
          
          ### go through model fits of requested iterations
          fc_new <- foreach(fit_i = fit[pos_constr], 
                            iter_i = seq_along(fit)[pos_constr],
                            .errorhandling = "pass") %do% {
                              
            ### overwrite landing fraction with last year, if requested
            if (!is.null(fwd_yrs_lf_remove)) {
              ### index for years to remove/overwrite
              idx_remove <- nrow(fit_i$data$landFrac) + fwd_yrs_lf_remove
              ### overwrite
              fit_i$data$landFrac[idx_remove, ] <- fit_i$data$landFrac[rep(nrow(fit_i$data$landFrac), length(idx_remove)), ]
            }
            
            ### check how to do forecast
            ### target recently advised TAC
            catchval <- ifelse(fwd_trgt == "TAC", c(TAC_last[,,,,, iter_i]), NA)
            ### scaled F
            fscale <- ifelse(fwd_trgt == "fsq", 1, NA)
            ### target new TAC after implementation of TAC constraint TAC year
            catchval[which(fwd_trgt == "hcr")[1]] <- catch_target[iter_i]
            ### target fsq (F that corresponds to TAC) in following years
            fscale[which(fwd_trgt == "hcr")[-1]] <- 1
            
            ### years for average values
            ave.years <- max(fit_i$data$years) + fwd_yrs_average
            ### years for sampling of recruitment years
            if (is.null(fwd_yrs_rec_start)) {
              rec.years <- fit_i$data$years ### use all years, if not defined
            } else {
              rec.years <- seq(from = fwd_yrs_rec_start, max(fit_i$data$years))
            }
            
            ### years where selectivity is not used for mean in forecast
            overwriteSelYears <- NULL
            if (!is.null(fwd_yrs_sel)) 
              overwriteSelYears <- max(fit_i$data$years) + fwd_yrs_sel
            
            ### arguments for forecast
            fc_args <- list(fit = fit_i, fscale = fscale, catchval = catchval,
                            ave.years = ave.years, rec.years = rec.years,
                            overwriteSelYears = overwriteSelYears, 
                            splitLD = fwd_splitLD)
            ### for compatibility of stockassessment's commit a882a11 and later:
            if ("savesim" %in% names(formals(base::args(stockassessment::forecast)))) {
              fc_args$savesim <- TRUE
            }
            ### run forecast
            fc_i <- do.call(stockassessment::forecast, fc_args)
            
            ### return forecast table
            return(attr(fc_i, "tab"))
            
          }
          
          ### overwrite updated forecasts
          fc[pos_constr] <- fc_new
          
          ### get SSB in TAC year
          SSB_TACyr <- sapply(fc, function(x) {
            if (is(x, "error")) {
              return(NA)
            } else {
              return(x[ac(ay + 1), "ssb:median"])
            }
          })
          ### get SSB in year after TAC year
          SSB_TACyr1 <- sapply(fc, function(x) {
            if (is(x, "error")) {
              return(NA)
            } else {
              return(x[ac(ay + 2), "ssb:median"])
            }
          })
          ### get F in TAC year
          F_TACyr <- sapply(fc, function(x) {
            if (is(x, "error")) {
              return(NA)
            } else {
              return(x[ac(ay + 1), "fbar:median"])
            }
          })
          
          ### check if SSB in TAC year below Bpa AND F above Fpa
          pos_fc1 <- which(SSB_TACyr < c(hcrpars["Bpa"]) & 
                             F_TACyr > c(hcrpars["Fpa"]))
          ### check if SSB below Bpa in TAC year and year after
          pos_fc2 <- which(SSB_TACyr < c(hcrpars["Bpa"]) & 
                             SSB_TACyr1 < c(hcrpars["Bpa"]))
          ### combine both conditions
          pos_fc <- union(pos_fc1, pos_fc2)
          
        }
        
      }
      
    }
    
    ### get position where BB is stopped
    pos_stop <- union(pos_hcr, pos_fc)
    ### stop banking and borrowing
    ### (but still pay back/use from last year)
    BB_bank[pos_stop] <- 0
    BB_borrow[pos_stop] <- 0
    
    ### correct target catch later in iem module
    # catch_target <- catch_target - c(BB_return) +
    #   c(BB_bank_use) - c(BB_bank) + c(BB_borrow)
    
  } else {
    ### if B&B not applied, store 0 
    BB_return <- BB_bank_use <- BB_bank <- BB_borrow <- 0
    
  }
  
  ### save B&B transfers in  tracking
  tracking[[1]]["BB_return", ac(ay)] <- BB_return
  tracking[[1]]["BB_bank_use", ac(ay)] <- BB_bank_use
  tracking[[1]]["BB_bank", ac(ay)] <- BB_bank
  tracking[[1]]["BB_borrow", ac(ay)] <- BB_borrow
  
  ### create ctrl object
  advice <- FLQuant(c(catch_target), 
                    dimnames = list(year = ctrl@target$year,
                                    iter = seq(it)))
  ctrl <- fwdControl(target = advice, quant = "catch")
  
  ### return catch target and tracking
  return(list(ctrl = ctrl, tracking = tracking))
  
}

### ------------------------------------------------------------------------ ###
### implementation error: banking and borrowing ####
### ------------------------------------------------------------------------ ###
### so far only banking and borrowing (B&B) implemented

iem_WKNSMSE <- function(tracking, ctrl,
                        args, ### contains ay (assessment year)
                        BB = FALSE, ### apply banking and borrowing
                        ...) {
  
  if (isTRUE(BB)) {
    
    ### get current (assessment) year
    ay <- args$ay
    
    ### retrieve banking and borrowing values
    BB_return <- tracking[[1]]["BB_return", ac(ay)]
    BB_bank_use <- tracking[[1]]["BB_bank_use", ac(ay)]
    BB_bank <- tracking[[1]]["BB_bank", ac(ay)]
    BB_borrow <- tracking[[1]]["BB_borrow", ac(ay)]
    
    ### replace NAs with 0
    BB_return <- ifelse(!is.na(BB_return), BB_return, 0)
    BB_bank_use <- ifelse(!is.na(BB_bank_use), BB_bank_use, 0)
    BB_bank <- ifelse(!is.na(BB_bank), BB_bank, 0)
    BB_borrow <- ifelse(!is.na(BB_borrow), BB_borrow, 0)
    
    ### get catch target(s)
    catch_target <- ctrl@iters[, "value", ]
    
    ### implement B&B
    catch_target <- catch_target - c(BB_return) + c(BB_bank_use) - 
      c(BB_bank) + c(BB_borrow)
    
    ### save in ctrl object
    ctrl@iters[, "value", ] <- catch_target
    
  }
  
  ### return catch target and tracking
  return(list(ctrl = ctrl, tracking = tracking))
  
}

### ------------------------------------------------------------------------ ###
### density-dependent natural mortality ####
### ------------------------------------------------------------------------ ###
### calculate density-dependent M from M2 relationships and stock.n 

calculate_ddM <- function(stk,
                          yr,
                          relation,
                          migration = TRUE) {
  
  ### read in M2 relationships
  M2 <- relation
  M2$pM2 <- M2$Nprey <- NA
  M2 <- array(rep(unlist(M2), dim(stk)[6]), 
              dim = c(nrow(M2), ncol(M2), dim(stk)[6]))
  dimnames(M2)[[2]] <- c("age", "predator", "intercept", "logbPred", "logbPrey",
                         "nPred", "nPrey", "pM2") 
  
  M_new <- m(stk)[, ac(yr)] %=% 0
  
  for (iter_i in seq(dim(stk)[6])) {
    
    ### extract number of predator and prey cod-at-age from stock.n in the specified year (000s)
    M2[, "nPrey", iter_i] <- stock.n(stk)[M2[, "age", iter_i], ac(yr),,,, iter_i]
    M2[is.na(M2[, "nPred", iter_i]), "nPred", iter_i] <- 
      stock.n(stk)[M2[!is.na(M2[, "predator", iter_i]), "predator", iter_i], 
                   ac(yr),,,, iter_i]
    
    ### calculate partial M2s from predator and prey abundances
    M2[, "pM2", iter_i] <- M2[, "intercept", iter_i] + 
      M2[, "logbPrey", iter_i] * log(M2[, "nPrey", iter_i])
    M2[!is.na(M2[, "logbPred", iter_i]), "pM2", iter_i] <- 
      M2[!is.na(M2[, "logbPred", iter_i]), "pM2", iter_i] + 
      M2[!is.na(M2[, "logbPred", iter_i]), "logbPred", iter_i] * 
      log(M2[!is.na(M2[, "logbPred", iter_i]), "nPred", iter_i])
    M2[, "pM2", iter_i] <- exp(M2[, "pM2", iter_i])
    
    ### sum M2s for each prey age class
    M2age <- aggregate(M2[, "pM2", iter_i], 
                       by = list(age = M2[, "age", iter_i]), sum)
    
    ### overwrite M in the specified year
    M_new[M2age$age,,,,, iter_i] <- M2age$x
    
  }
  ### and add 0.2 for non-predation mortality
  M_new <- M_new + 0.2
  
  ### "migration" correction
  if (isTRUE(migration)) {
    ### 15% migration for ages 3+
    M_new[dimnames(M_new)$age[dimnames(M_new)$age >= 3]] <- 
      M_new[dimnames(M_new)$age[dimnames(M_new)$age >= 3]] - log(1 - 0.15)
  }
  
  return(M_new)
  
}
