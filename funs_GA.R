### ------------------------------------------------------------------------ ###
### objective function for multi species run ####
### ------------------------------------------------------------------------ ###
mp_fitness <- function(params, inp_file, path, check_file = FALSE,
                   scenario, MP, 
                   return_res = FALSE,
                   save_MP = FALSE, ### save MP results
                   collapse_correction = TRUE,
                   obj_fun = "ICES", ### objective function (or elements)
                   obj_yrs = "all", ### years to use in objective function
                   stat_yrs = "all", ### years for summary statistics
                   pen_neg = FALSE,
                   pen_max = 1,
                   pen_infl = 0.06,
                   pen_steep = 0.5e+3,
                   ...) {
  
  ### housekeeping
  invisible(gc())
  if (exists("res_mp")) {
    rm(res_mp)
    invisible(gc())
  }
  if (getDoParWorkers() > 1)
    . <- foreach(i = 1:getDoParWorkers()) %dopar% {invisible(gc())}
  
  ### rounding of arguments
  if (identical(MP, "rfb")) {
    params[1:4] <- round(params[1:4])
    params[5:7] <- round(params[5:7], 1)
    params[8] <- round(params[8])
    params[9] <- round(params[9], 2)
    params[10:11] <- round(params[10:11], 2)
    if (is.nan(params[10])) params[10] <- Inf ### fix NaN for upper_constraint
  } else if (identical(MP, "hr")) {
    ### idxB_lag, idxB_range_3, interval [years]
    params[c(1, 2, 5)] <- round(params[c(1, 2, 5)])
    ### exp_b, comp_b_multiplier
    params[c(3, 4)] <- round(params[c(3, 4)], 1)
    ### multiplier, upper_constraint, lower_constraint
    params[c(6, 7, 8)] <- round(params[c(6, 7, 8)], 2)
    ### fix NaN for upper_constraint
    if (is.nan(params[7])) params[7] <- Inf
  }
  
  ### check for files?
  run_mp <- TRUE
  ### current run
  run_i <- paste0(params, collapse = "_")
  if (isTRUE(check_file)) {
    ### check if path exists
    if (!dir.exists(path)) dir.create(path, recursive = TRUE)
    ### check if run already exists
    if (isTRUE(file.exists(paste0(path, run_i, ".rds")))) {
      ### load stats
      stats <- readRDS(paste0(path, run_i, ".rds"))
      ### set flag for running MP
      run_mp <- FALSE
      ### use different period to calculate stats?
      if (!any(stat_yrs %in% c("all", "multiple"))) {
        if (!any(grepl(x = rownames(stats), pattern = stat_yrs))) run_mp <- TRUE
      }
    }
  }
  
  if (isTRUE(run_mp)) {
    
    ### load input file from disk
    input <- readRDS(inp_file)
    
    ### insert arguments into input object for mp
    if (identical(MP, "rfb")) {
      input$ctrl$est@args$idxB_lag          <- params[1]
      input$ctrl$est@args$idxB_range_1      <- params[2]
      input$ctrl$est@args$idxB_range_2      <- params[3]
      input$ctrl$est@args$catch_range       <- params[4]
      input$ctrl$est@args$comp_m            <- params[9]
      input$ctrl$phcr@args$exp_r            <- params[5]
      input$ctrl$phcr@args$exp_f            <- params[6]
      input$ctrl$phcr@args$exp_b            <- params[7]
      input$ctrl$hcr@args$interval          <- params[8]
      input$ctrl$isys@args$interval         <- params[8]
      input$ctrl$isys@args$upper_constraint <- params[10]
      input$ctrl$isys@args$lower_constraint <- params[11]
    } else if (identical(MP, "hr")) {
      ### biomass index 
      input$ctrl$est@args$idxB_lag <- params[1]
      input$ctrl$est@args$idxB_range_3 <- params[2]
      ### biomass safeguard
      input$ctrl$phcr@args$exp_b <- params[3]
      ### change Itrigger? (default: Itrigger=1.4*Iloss)
      if (isFALSE(params[4] == 1.4)) {
        input$ctrl$est@args$I_trigger <- 
          input$ctrl$est@args$I_trigger/1.4*params[4]
      }
      ### multiplier
      input$ctrl$est@args$comp_m <- params[6]
      ### catch interval (default: 1)
      if (is.numeric(params[5])) {
        input$ctrl$hcr@args$interval <- params[5]
        input$ctrl$isys@args$interval <- params[5]
      }
      ### catch constraint
      input$ctrl$isys@args$upper_constraint <- params[7]
      input$ctrl$isys@args$lower_constraint <- params[8]
    }
    
    ### run MP
    res_mp <- do.call(mp, input)
    
    if (isTRUE(return_res)) {
      return(res_mp)
    }
    if (isTRUE(save_MP)) {
      saveRDS(res_mp, paste0(path, "mp_", run_i, ".rds"))
    }
    
    ### calculate stats
    stats <- mp_stats(input = input, res_mp = res_mp, stat_yrs = stat_yrs,
                      collapse_correction = collapse_correction)
    
    ### save result in file
    if (isTRUE(check_file)) {
      saveRDS(stats, paste0(path, run_i, ".rds"))
    }
    
  }
  
  ### prepare stats for objective function
  stat_names <- c("risk_Blim", "risk_Blim_max", "risk_Bmsy", 
                  "risk_halfBmsy", "risk_collapse", "SSB", "Fbar", 
                  "Catch", "SSB_rel", "Fbar_rel", "Catch_rel", "ICV")
  if (identical(obj_yrs, "all")) {
    stats_obj <- stats[stat_names]
  } else {
    stats_obj <- stats[paste0(obj_yrs, "_", stat_names)]
    names(stats_obj) <- gsub(x = names(stats_obj), 
                             pattern = paste0(obj_yrs, "_"), replacement = "")
  }
  ### initialise objective function
  obj <- 0
  ### MSY objectives: target MSY reference values
  if (isTRUE("SSB" %in% obj_fun)) 
    obj <- obj - sum(abs(unlist(stats_obj$SSB_rel) - 1))
  if (isTRUE("catch" %in% obj_fun)) 
    obj <- obj - sum(abs(unlist(stats_obj$Catch_rel) - 1))
  if (isTRUE("F" %in% obj_fun)) 
    obj <- obj - sum(abs(unlist(stats_obj$Fbar_rel) - 1))
  ### reduce risk & ICV
  if (isTRUE("risk" %in% obj_fun)) 
    obj <- obj - sum(unlist(stats_obj$risk_Blim))
  if (isTRUE("ICV" %in% obj_fun)) 
    obj <- obj - sum(unlist(stats_obj$ICV))
  ### MSY target but replace risk with PA objective
  if (isTRUE(obj_fun == "MSYPA")) {
    obj <- obj - sum(abs(unlist(stats_obj$SSB_rel) - 1)) -
      sum(abs(unlist(stats_obj$Catch_rel) - 1)) -
      sum(unlist(stats_obj$ICV)) -
      sum(penalty(x = unlist(stats_obj$risk_Blim), 
                  negative = pen_neg, max = pen_max, 
                  inflection = pen_infl, 
                  steepness = pen_steep))
      ### max penalty: 5
      ### for pollack zero catch has fitness of -4.7
  }
  ### ICES MSY approach, maximise catch while keeping risk <= 0.05
  if (isTRUE(obj_fun == "ICES")) {
    obj <- obj + stats_obj$Catch_rel -
      sum(penalty(x = stats_obj$risk_Blim_max, 
                  negative = pen_neg, max = pen_max, 
                  inflection = pen_infl, 
                  steepness = pen_steep))
  }
  
  
  ### housekeeping
  rm(res_mp, input)
  invisible(gc())
  if (getDoParWorkers() > 1)
    . <- foreach(i = 1:getDoParWorkers()) %dopar% {invisible(gc())}
  
  ### return objective function (fitness) value
  return(obj)
  
}

### ------------------------------------------------------------------------ ###
### stats from MSE run(s) ####
### ------------------------------------------------------------------------ ###

### function for calculating stats
mp_stats <- function(input, res_mp, stat_yrs = "multiple", 
                     collapse_correction = TRUE, start_yr = input$args$iy) {
  
  ### stock metrics
  SSBs <- FLCore::window(ssb(res_mp@om@stock), start = start_yr + 1)
  Fs <- FLCore::window(fbar(res_mp@om@stock), start = start_yr + 1)
  Cs <- FLCore::window(catch(res_mp@om@stock), start = start_yr + 1)
  yrs <- dim(SSBs)[2]
  its <- dim(SSBs)[6]
  ### collapse correction
  if (isTRUE(collapse_correction)) {
    ### find collapses
    cd <- sapply(seq(its), function(x) {
      min_yr <- min(which(SSBs[,,,,, x] < 1))
      if (is.finite(min_yr)) {
        all_yrs <- min_yr:yrs
      } else {
        all_yrs <- NA
      }
      all_yrs + (x - 1)*yrs
    })
    cd <- unlist(cd)
    cd <- cd[which(!is.na(cd))]
    ### remove values
    SSBs@.Data[cd] <- 0
    Cs@.Data[cd] <- 0
    Fs@.Data[cd] <- 0
  }
  ### extend Catch to include ICV calculation from last historical year
  Cs_long <- FLCore::window(Cs, start = start_yr)
  Cs_long[, ac(start_yr)] <- catch(res_mp@om@stock)[, ac(start_yr)]
  ### refpts
  Bmsy <- c(input$refpts["Bmsy"])
  Fmsy <- c(input$refpts["Fmsy"])
  Cmsy <- c(input$refpts["Cmsy"])
  Blim <- c(input$refpts["Blim"])
  ### TAC interval
  if (!is.null(input$ctrl$hcr@args$interval)) {
    TAC_intvl <- input$ctrl$hcr@args$interval
  } else {
    TAC_intvl <- 1
  }
  
  ### some stats
  stats_list <- function(SSBs, Cs, Fs, Cs_long, Blim, Bmsy, Fmsy, Cmsy,
                         TAC_intvl) {
    list(
      risk_Blim = mean(c((SSBs/Blim) < 1), na.rm = TRUE),
      risk_Blim_max = max(apply((SSBs/Blim) < 1, 2, mean, na.rm = TRUE), 
                          na.rm = TRUE),
      risk_Bmsy = mean(c((SSBs/Bmsy) < 1), na.rm = TRUE),
      risk_halfBmsy = mean(c((SSBs/(Bmsy/2)) < 1), na.rm = TRUE),
      risk_collapse = mean(c(SSBs < 1), na.rm = TRUE),
      SSB = median(c(SSBs), na.rm = TRUE), Fbar = median(c(Fs), na.rm = TRUE),
      Catch = median(c(Cs), na.rm = TRUE),
      SSB_rel = median(c(SSBs/Bmsy), na.rm = TRUE),
      Fbar_rel = median(c(Fs/Fmsy), na.rm = TRUE),
      Catch_rel = median(c(Cs/Cmsy), na.rm = TRUE),
      ICV = iav(Cs_long, from = start_yr, period = TAC_intvl,
                summary_all = median)
    )
  }
  ### stats for full period
  stats <- stats_list(SSBs = SSBs, Cs = Cs, Fs = Fs, 
                        Cs_long = Cs_long, 
                        Blim = Blim, Bmsy = Bmsy, Fmsy = Fmsy, Cmsy = Cmsy,
                        TAC_intvl = TAC_intvl)
  ### additional time periods?
  if (!identical(stat_yrs, "all")) {
    ### list of possible years
    yrs_labels <- c("1:5", "6:10", "11:20", "1:10", "1:20",
                    "1:30", "1:40", "1:50", "1:100", "40:50", "51:100", 
                    "91:100")
    names(yrs_labels) <- yrs_labels
    yrs_vals <- lapply(yrs_labels, function(x) eval(parse(text = x)))
    ### find available years and remove impossible years
    yrs_avail <- seq(dim(SSBs)[2])
    pos_keep <- which(sapply(yrs_vals, function(x) all(x %in% yrs_avail)))
    yrs_vals <- yrs_vals[pos_keep]
    yrs_labels <- yrs_labels[pos_keep]
    ### calculate stats for these years
    stats_add <- lapply(yrs_vals, function(x) {
      ### define years for summary statistics
      yrs_tmp <- x
      yrs_tmpp1 <- seq(from = min(as.numeric(yrs_tmp)),
                       to = max(as.numeric(yrs_tmp) + 1))
      stats_tmp <- c(stats_list(SSBs = SSBs[, yrs_tmp], 
                                Cs = Cs[, yrs_tmp],
                                Fs = Fs[, yrs_tmp], 
                                Cs_long = Cs_long[, yrs_tmpp1],
                                Blim = Blim, Bmsy = Bmsy, Fmsy = Fmsy, 
                                Cmsy = Cmsy, TAC_intvl = TAC_intvl))
      names(stats_tmp) <- paste0(paste0(head(yrs_tmp, 1), ":",
                                        tail(yrs_tmp, 1)),
                                 "_", names(stats_tmp))
      return(stats_tmp)
    })
    names(stats_add) <- NULL
    stats <- c(stats, unlist(stats_add))
  }
  
  return(stats)
  
}

### ------------------------------------------------------------------------ ###
### penalty function ####
### ------------------------------------------------------------------------ ###

penalty <- function(x, negative = FALSE, max = 1,
                    inflection = 0.06, steepness = 0.5e+3) {
  y <- max / (1 + exp(-(x - inflection)*steepness))
  if (isTRUE(negative)) y <- -y
  return(y)
}
