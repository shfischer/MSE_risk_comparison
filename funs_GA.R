### ------------------------------------------------------------------------ ###
### objective function for multi species run ####
### ------------------------------------------------------------------------ ###
mp_fitness <- function(params, inp_file, path, check_file = FALSE,
                   scenario, 
                   return_res = FALSE,
                   collapse_correction = TRUE,
                   obj_SSB = TRUE, obj_C = TRUE, obj_F = FALSE,
                   obj_risk = TRUE, obj_ICV = TRUE, obj_ICES_PA = FALSE,
                   obj_ICES_PA2 = FALSE, obj_ICES_MSYPA = FALSE,
                   stat_yrs = "all",
                   risk_threshold = 0.05,
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
  params[1:4] <- round(params[1:4])
  params[5:7] <- round(params[5:7], 1)
  params[8] <- round(params[8])
  params[9] <- round(params[9], 2)
  params[10:11] <- round(params[10:11], 2)
  ### fix NaN for upper_constraint
  if (is.nan(params[10])) params[10] <- Inf
  
  ### check for files?
  if (isTRUE(check_file)) {
    ### current run
    run_i <- paste0(params, collapse = "_")
    ### get current stock(s)
    stock_i <- strsplit(x = tail(strsplit(x = path, split = "/")[[1]], 1), 
                        split = "_")[[1]]
    base_path <- paste0(paste0(head(strsplit(x = path, split = "/")[[1]], -1), 
                               collapse = "/"), "/")
    ### check if path exists
    if (!dir.exists(path)) dir.create(path, recursive = TRUE)
    ### check if run already exists
    if (isTRUE(file.exists(paste0(path, run_i, ".rds")))) {
      ### load stats
      stats <- readRDS(paste0(path, run_i, ".rds"))
      ### set flag for running MP
      run_mp <- FALSE
      ### use different period to calculate stats?
      if (!any(stat_yrs %in% c("all", "more"))) {
        if (!any(grepl(x = rownames(stats), pattern = stat_yrs))) run_mp <- TRUE
      }
    } else {
      ### check if run exist in larger group
      dir_i <- paste0(stock_i, collapse = "_")
      dirs_i <- setdiff(x = dir(path = base_path, pattern = dir_i),
                        y = dir_i)
      if (isTRUE(length(dirs_i) > 0)) {
        dirs_i <- dirs_i[which(sapply(dirs_i, function(x) {
          tmp <- strsplit(x = x, split = "_")[[1]]
          ifelse(isFALSE(dir_i %in% tmp), FALSE, TRUE)
        }))]
        files_tmp <- lapply(dirs_i, function(x) {
          #browser()
          path_tmp <- paste0(base_path, x, "/", run_i, ".rds")
          if (isTRUE(file.exists(path = path_tmp))) {
            return(path_tmp)
          } else {
            return(NA)
          }
        })
        files_tmp[is.na(files_tmp)] <- NULL
        if (isTRUE(length(files_tmp) > 0)) {
          ### load stats from larger group
          stats <- readRDS(files_tmp[[1]])
          ### subset to current group
          stats <- stats[, stock_i]
          ### do not run MP
          run_mp <- FALSE
        } else {
          run_mp <- TRUE
        }
      } else {
        run_mp <- TRUE
      }
    }
  } else {
    run_mp <- TRUE
  }
  
  if (isTRUE(run_mp)) {
    
    ### load input file from disk
    input <- readRDS(inp_file)
    
    ### insert arguments into input object for mp
    input <- lapply(input, function(x) {
      x$ctrl$est@args$idxB_lag     <- params[1]
      x$ctrl$est@args$idxB_range_1 <- params[2]
      x$ctrl$est@args$idxB_range_2 <- params[3]
      x$ctrl$est@args$catch_range  <- params[4]
      x$ctrl$est@args$comp_m <- params[9]
      x$ctrl$phcr@args$exp_r <- params[5]
      x$ctrl$phcr@args$exp_f <- params[6]
      x$ctrl$phcr@args$exp_b <- params[7]
      x$ctrl$hcr@args$interval <- params[8]
      x$ctrl$isys@args$interval <- params[8]
      x$ctrl$isys@args$upper_constraint <- params[10]
      x$ctrl$isys@args$lower_constraint <- params[11]
      
      return(x)
    })
    
    ### if group of stocks, check if results for individual stocks exist
    group <- ifelse(isTRUE(length(input) > 1) & isTRUE(check_file), TRUE, FALSE)
    if (group) {
      ### get paths
      group_stocks <- names(input)
      path_base <- gsub(x = path, 
                        pattern = paste0(paste0(group_stocks, collapse = "_"), 
                                         "/"),
                        replacement = "")
      path_stocks <- paste0(path_base, group_stocks, "/")
      ### check for files
      run_exists <- file.exists(paste0(path_stocks, run_i, ".rds"))
      group <- ifelse(any(run_exists), TRUE, FALSE)
      
      ### do some results exist?
      if (group) {
        ### load results
        files_exist <- paste0(path_stocks, run_i, ".rds")[run_exists]
        stats_group <- lapply(files_exist, readRDS)
        names(stats_group) <- group_stocks[run_exists]
        ### get stocks which require simulation
        run_stocks <- group_stocks[!run_exists]
        ### subset input
        input <- input[run_stocks]
        
      }
      
    }
    
    ### run MP for each list element
    res_mp <- lapply(input, function(x) {
      if (getDoParWorkers() > 1)
        . <- foreach(i = 1:getDoParWorkers()) %dopar% {invisible(gc())}
      do.call(mp, x)
    })
    
    if (isTRUE(return_res)) {
      return(res_mp)
    }
    
    ### calculate stats
    stat_yrs_calc <- "more" ### always calculate all periods
    stats <- mp_stats(input = input, res_mp = res_mp, stat_yrs = stat_yrs_calc,
                      collapse_correction = collapse_correction)
    
    ### add existing results for stock groups
    if (group) {
      
      ### split old stats into list
      if (isTRUE(length(stats) > 0)) {
        stats <- asplit(stats, MARGIN = 2)
      }
      ### stats_group is already a list
      ### combine new and existing stats
      stats <- c(stats_group, stats)
      ### sort and coerce into matrix
      stats <- stats[group_stocks]
      stats <- do.call(cbind, stats)
      
    }
    
    ### save result in file
    if (isTRUE(check_file)) {
      saveRDS(stats, paste0(path, run_i, ".rds"))
    }
    
  }
  
  ### prepare stats for objective function
  if (identical(stat_yrs, "all") | identical(stat_yrs, "more")) {
    SSB_rel <- stats["SSB_rel", ]
    Catch_rel <- stats["Catch_rel", ]
    Fbar_rel <- stats["Fbar_rel", ]
    risk_Blim <- stats["risk_Blim", ]
    ICV <- stats["ICV", ]
  } else if (stat_yrs %in% c("first10", "41to50", "last10", "firsthalf",
                             "lastfhalf", "11to50")) {
    SSB_rel <- stats[paste0("SSB_rel_", stat_yrs), ]
    Catch_rel <- stats[paste0("Catch_rel_", stat_yrs), ]
    Fbar_rel <- stats[paste0("Fbar_rel_", stat_yrs), ]
    risk_Blim <- stats[paste0("risk_Blim_", stat_yrs), ]
    ICV <- stats[paste0("ICV_", stat_yrs), ]
    
  } else if (identical(stat_yrs, "last10")) {
    SSB_rel <- stats["SSB_rel_last10", ]
    Catch_rel <- stats["Catch_rel_last10", ]
    Fbar_rel <- stats["Fbar_rel_last10", ]
    risk_Blim <- stats["risk_Blim_last10", ]
    ICV <- stats["ICV_last10", ]
  }
  ### objective function
  obj <- 0
  ### MSY objectives: target MSY reference values
  if (isTRUE(obj_SSB)) obj <- obj - sum(abs(unlist(SSB_rel) - 1))
  if (isTRUE(obj_C)) obj <- obj - sum(abs(unlist(Catch_rel) - 1))
  if (isTRUE(obj_F)) obj <- obj - sum(abs(unlist(Fbar_rel) - 1))
  ### reduce risk & ICV
  if (isTRUE(obj_risk)) obj <- obj - sum(unlist(risk_Blim))
  if (isTRUE(obj_ICV)) obj <- obj - sum(unlist(ICV))
  ### ICES approach: maximise catch while keeping risk <5%
  if (isTRUE(obj_ICES_PA)) {
    obj <- obj + sum(unlist(Catch_rel))
    ### penalise risk above 5%
    obj <- obj - sum(ifelse(test = unlist(risk_Blim) <= 0.05,
                            yes = 0,
                            no = 10)) 
  }
  if (isTRUE(obj_ICES_PA2)) {
    obj <- obj + sum(unlist(Catch_rel))
    ### penalise risk above 5% - gradual
    obj <- obj - sum(penalty(x = unlist(risk_Blim), 
                             negative = FALSE, max = 1, inflection = 0.06, 
                             steepness = 0.5e+3))
  }
  ### MSY target but replace risk with PA objective
  if (isTRUE(obj_ICES_MSYPA)) {
    obj <- obj - sum(abs(unlist(SSB_rel) - 1)) -
      sum(abs(unlist(Catch_rel) - 1)) -
      sum(unlist(ICV)) -
      sum(penalty(x = unlist(risk_Blim), 
                             negative = FALSE, max = 5, 
                             inflection = risk_threshold + 0.01, 
                             steepness = 0.5e+3))
      ### max penalty: 5
      ### for pollack zero catch has fitness of -4.7
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
  SSBs <- FLCore::window(ssb(res_mp@stock), start = start_yr + 1)
  Fs <- FLCore::window(fbar(res_mp@stock), start = start_yr + 1)
  Cs <- FLCore::window(catch(res_mp@stock), start = start_yr + 1)
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
  Cs_long[, ac(start_yr)] <- catch(res_mp@stock)[, ac(start_yr)]
  ### refpts
  Bmsy <- c(input$refpts["Bmsy"])
  Fmsy <- c(input$refpts["Fmsy"])
  Cmsy <- c(input$refpts["Cmsy"])
  Blim <- c(input$refpts["Blim"])
  ### TAC interval
  TAC_intvl <- input$ctrl$hcr@args$interval
  
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
