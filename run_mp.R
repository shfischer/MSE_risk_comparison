### ------------------------------------------------------------------------ ###
### run MSE ####
### ------------------------------------------------------------------------ ###

### ------------------------------------------------------------------------ ###
### arguments ####
### ------------------------------------------------------------------------ ###

args <- commandArgs(TRUE)
print("arguments passed on to this script:")
print(args)

### evaluate arguments, if they are passed to R:
if (length(args) > 0) {
  
  ### extract arguments
  for (i in seq_along(args)) eval(parse(text = args[[i]]))
  ### set default arguments
  ### parallelization
  if (!exists("use_MPI")) use_MPI <- FALSE
  if (!exists("n_blocks")) n_blocks <- 1
  if (!exists("n_workers")) n_workers <- 0
  ### scenario definition
  if (!exists("n_iter")) n_iter <- 500
  if (!exists("n_yrs")) n_yrs <- 50
  if (!exists("fhist")) fhist <- "one-way"
  if (!exists("catch_rule")) catch_rule <- "catch_rule"
  if (!exists("comp_r")) comp_r <- TRUE
  if (!exists("comp_f")) comp_f <- TRUE
  if (!exists("comp_b")) comp_b <- TRUE
  if (!exists("scenario")) scenario <- "PA"
  if (!exists("cap_below_b")) cap_below_b <- TRUE
  ### GA search
  if (!exists("ga_search")) ga_search <- TRUE
  if (isTRUE(ga_search)) {
    if (!exists("popSize")) stop("popSize missing")
    if (!exists("maxiter")) stop("maxiter missing")
    if (!exists("stock_id")) stop("stock_id missing")
    if (!exists("run")) run <- maxiter
    if (!exists("collate")) collate <- TRUE
    ### objective function elements
    if (!exists("obj_SSB")) obj_SSB <- TRUE
    if (!exists("obj_F")) obj_F <- FALSE
    if (!exists("obj_C")) obj_C <- TRUE
    if (!exists("obj_risk")) obj_risk <- TRUE
    if (!exists("obj_ICV")) obj_ICV <- TRUE
    if (!exists("obj_ICES_PA")) obj_ICES_PA <- FALSE
    if (!exists("obj_ICES_PA2")) obj_ICES_PA2 <- FALSE
    if (!exists("obj_ICES_MSYPA")) obj_ICES_MSYPA <- FALSE
    if (!exists("risk_threshold")) risk_threshold <- 0.05
    ### GA
    if (!exists("add_suggestions")) add_suggestions <- TRUE
    if (!exists("stat_yrs")) stat_yrs <- "all"
  }
  
} else {
  
  stop("no argument passed to R")
  
}

### ------------------------------------------------------------------------ ###
### set up environment ####
### ------------------------------------------------------------------------ ###

### load packages
### GA fork from GitHub remotes::install_github("shfischer/GA")
### use mse fork from shfischer/mse, branch mseDL2.0 
### remotes::install_github("shfischer/mse", ref = "mseDL2.0)
req_pckgs <- c("FLCore", "FLash", "mse", "GA", "doParallel", "doRNG", "FLBRP")
for (i in req_pckgs) library(package = i, character.only = TRUE)

### load additional functions
source("funs.R")
source("funs_GA.R")

### ------------------------------------------------------------------------ ###
### setup parallel environment ####
### ------------------------------------------------------------------------ ###

### hybrid MPI
if (isTRUE(use_MPI)) {
  ### 1st: doMPI cluster with 1 worker per node
  message("starting doMPI")
  library(doMPI)
  cl1 <- startMPIcluster()
  message("startMPIcluster() succeeded")
  print(cl1)
  registerDoMPI(cl1)
  cl_length_1 <- cl1$workerCount
  cl_length_1
  
  ### 2nd: doParallel workers inside doMPI workers
  . <- foreach(i = seq(cl_length_1)) %dopar% {
    ### load packages and functions into MPI workers
    for (i in req_pckgs) library(package = i, character.only = TRUE,
                                 warn.conflicts = FALSE, verbose = FALSE,
                                 quietly = TRUE)
  }
  message("MPI package loading succeeded")
  . <- foreach(i = seq(cl_length_1)) %dopar% {
    source("funs.R", echo = FALSE)
    source("funs_GA.R", echo = FALSE)
  }
  message("MPI script loading succeeded")
  ### start doParallel inside MPI processes
  if (isTRUE(n_workers > 1)) {
    . <- foreach(i = seq(cl_length_1)) %dopar% {
      cl2 <- makeCluster(n_workers)
      registerDoParallel(cl2)
      cl_length_2 <- length(cl2)
      ### load packages and functions into parallel workers
      . <- foreach(i = seq(cl_length_2)) %dopar% {
        for (i in req_pckgs) library(package = i, character.only = TRUE,
                                     warn.conflicts = FALSE, verbose = FALSE,
                                     quietly = TRUE)
        source("funs.R", echo = FALSE)
        source("funs_GA.R", echo = FALSE)
      }
    }
  }
  message("setting up doParallel inside MPI succeeded")
} else {
  if (isTRUE(n_workers > 1)) {
    ### start doParallel cluster
    cl1 <- makeCluster(n_workers)
    registerDoParallel(cl1)
    cl_length_1 <- length(cl1)
    ### load packages and functions into parallel workers
    . <- foreach(i = seq(cl_length_1)) %dopar% {
      for (i in req_pckgs) library(package = i, character.only = TRUE,
                                   warn.conflicts = FALSE, verbose = FALSE,
                                   quietly = TRUE)
      source("funs.R", echo = FALSE)
      source("funs_GA.R", echo = FALSE)
    }
  } else {
    cl1 <- FALSE
  }
}

### ------------------------------------------------------------------------ ###
### load data ####
### ------------------------------------------------------------------------ ###

stocks <- read.csv("input/stocks.csv", stringsAsFactors = FALSE)
stock <- stocks$stock[stock_id]
names(stock) <- stock
input <- lapply(stock, function(x) {
  readRDS(paste0("input/", n_iter, "_", n_yrs, "/OM_2_mp_input/", fhist, "/", x,
                 ".rds"))
})


### ------------------------------------------------------------------------ ###
### specify scenario ####
### ------------------------------------------------------------------------ ###

### default catch rule
input <- lapply(input, function(x) {
  ### OEM: activate uncertainty
  x$oem@args$idx_dev <- TRUE
  x$oem@args$ssb_idx <- FALSE
  x$oem@args$tsb_idx <- FALSE
  x$oem@args$lngth <- TRUE
  x$oem@args$lngth_dev <- TRUE
  ### IEM: do not activate uncertainty
  x$iem@args$use_dev <- FALSE
  ### catch rule components
  x$ctrl$est@args$comp_r <- comp_r
  x$ctrl$est@args$comp_f <- comp_f
  x$ctrl$est@args$comp_b <- comp_b
  ### catch lag fixed
  x$ctrl$est@args$catch_lag <- 1
  ### turn of uncertainty cap when index below Itrigger?
  if (isFALSE(cap_below_b)) {
    x$ctrl$isys@args$cap_below_b <- cap_below_b
    #x$ctrl$isys@method <- is_comps
  }
  
  return(x)
})

### default ICES rule: 2 over 3
if (isTRUE(catch_rule == "2over3")) {
  input <- lapply(input, function(x) {
    ### OEM: turn of length index
    x$oem@args$lngth <- FALSE
    x$oem@args$lngth_dev <- FALSE
    ### add PA buffer stock status and deviation
    x$oem@args$PA_status <- TRUE
    x$oem@args$PA_status_dev <- TRUE
    ### catch rule components: turn of f & b, 2 over 3 rule
    x$ctrl$est@args$comp_f <- FALSE
    x$ctrl$est@args$comp_b <- FALSE
    x$ctrl$est@args$idxB_lag <- 1
    x$ctrl$est@args$idxB_range_1 <- 2
    x$ctrl$est@args$idxB_range_2 <- 3
    ### PA buffer
    x$ctrl$est@args$pa_buffer <- TRUE
    ###
    x$ctrl$phcr@args$exp_r <- 1
    x$ctrl$phcr@args$exp_f <- 0
    x$ctrl$phcr@args$exp_b <- 1 ### PA buffer
    ### biennial
    #x$ctrl$hcr@args$interval <- 2
    ### uncertainty cap
    x$ctrl$isys@args$upper_constraint <- 1.2
    x$ctrl$isys@args$lower_constraint <- 0.8
    return(x)
  })
  # input$pol$oem@method <- wklife_3.2.1_obs
  # input$pol$ctrl$est@method <- wklife_3.2.1_est
  # input$pol$ctrl$isys@method <- is_r
  #debugonce(goFishDL)
  #res <- do.call(mp, c(input$pol, cut_hist = FALSE))
  
} else if (isTRUE(catch_rule == "hr")) {
  input <- lapply(input, function(x) {
    ### OEM
    x$oem@method <- obs_generic
    x$oem@args$lngth <- FALSE
    x$oem@args$lngth_dev <- FALSE
    x$oem@args$ssb_idx <- FALSE
    x$oem@args$tsb_idx <- TRUE
    ### est
    x$ctrl$est@method <- est_hr
    x$ctrl$est@args$idxB_lag <- 1
    x$ctrl$est@args$idxB_range <- 1
    ### phcr
    x$ctrl$phcr@method <- phcr_hr
    x$ctrl$phcr@args$rate <- hr_rate
    ### hcr
    x$ctrl$hcr@method <- hcr_hr
    x$ctrl$hcr@args$interval <- 1
    ### is
    x$ctrl$isys@method <- is_r
    x$ctrl$isys@args$interval <- 1
    x$ctrl$isys@args$upper_constraint <- Inf
    x$ctrl$isys@args$lower_constraint <- 0
    return(x)
  })
}

### within scenario parallelisation?
if (isTRUE(n_workers > 1) & isTRUE(n_blocks > 1)) {
  ### use Iloss
  input <- lapply(input, function(x) {
    x$args$nblocks <- n_blocks
    return(x)
  })
}

### ------------------------------------------------------------------------ ###
### GA set-up ####
### ------------------------------------------------------------------------ ###
if (isTRUE(catch_rule == "catch_rule") & isTRUE(ga_search)) {
  
  ### GA arguments
  ga_names <- c("lag_idx", "range_idx_1", "range_idx_2", "range_catch",
                "exp_r", "exp_f", "exp_b", "interval", "multiplier",
                "upper_constraint", "lower_constraint")
  ga_suggestions <- rbind(c(0, 1, 1, 1, 1, 1, 1, 1, 1, Inf, 0), ### most current data
                          c(0, 1, 1, 1, 1, 1, 1, 2, 1, Inf, 0),
                          c(0, 1, 1, 1, 0, 0, 0, 2, 1, Inf, 0), ### constant catch
                          ### default, annual/biennial, turning off elements
                          expand.grid(1, 2, 3, 1, 0:1, 0:1, 0:1, 1:2, 0:1, 
                                      Inf, 0),
                          ### default uncertainty cap
                          c(1, 2, 3, 1, 1, 1, 1, 2, 1, 1.2, 0.8))
  ga_default <- c(1, 2, 3, 1, 1, 1, 1, 2, 1, Inf, 0)
  ga_lower <- c(0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 0)
  ga_upper <- c(1, 5, 5, 1, 2, 2, 2, 5, 2, 5, 1)
  ### turn of parameters not requested, i.e. limit to default value
  pos_default <- which(!unlist(mget(ga_names, ifnotfound = FALSE)))
  ga_lower[pos_default] <- ga_default[pos_default]
  ga_upper[pos_default] <- ga_default[pos_default]
  ### fix parameters?
  pos_fixed <- which(sapply(mget(ga_names, ifnotfound = FALSE), is.numeric))
  par_fixed <- names(pos_fixed)
  val_fixed <- as.vector(unlist(mget(ga_names, ifnotfound = FALSE)[pos_fixed]))
  ga_lower[pos_fixed] <- val_fixed
  ga_upper[pos_fixed] <- val_fixed
  ### remove not requested parameters from suggestions
  ga_suggestions[, pos_default] <- rep(ga_default[pos_default], 
                                       each = nrow(ga_suggestions))
  ga_suggestions[, pos_fixed] <- rep(val_fixed, 
                                     each = nrow(ga_suggestions))
  ga_suggestions <- unique(ga_suggestions)
  names(ga_suggestions) <- ga_names
  
  ### multiplier only: run all possible values
  if (exists("multiplier")) {
    if (isTRUE(multiplier) & 
        !any(sapply(mget(setdiff(ga_names, "multiplier"), ifnotfound = FALSE),
                    isTRUE))) {
      m_vals <- seq(from = ga_lower[9], to = ga_upper[9], by = 0.01)
      ga_suggestions[1, ] <- ga_lower
      ga_suggestions <- ga_suggestions[rep(1, length(m_vals)), ]
      ga_suggestions$multiplier <- m_vals
      ### adapt GA dimensions
      maxiter <- run <- 1
      popSize <- length(m_vals)
    }
  }
  
  ### ---------------------------------------------------------------------- ###
  ### paths ####
  ### ---------------------------------------------------------------------- ###
  
  ### output path
  ### set name depending on which GA parameters are used
  scn_pars <- ga_names[setdiff(seq_along(ga_names), pos_default)]
  scn_pars[which(scn_pars %in% par_fixed)] <- paste0(
    scn_pars[which(scn_pars %in% par_fixed)], val_fixed)
  scn_pars_c <- paste0(scn_pars, collapse = "-")
  #scenario <- "trial"
  path_out <- paste0("output/", n_iter, "_", n_yrs, "/", scenario, "/",
                     fhist, "/",
                     paste0(stock, collapse = "_"), "/")
  dir.create(path_out, recursive = TRUE)
  
  ### objective function elements
  obj_fun <- c("SSB", "F", "C", "risk", "ICV", "ICES_PA", "ICES_PA2",
               "ICES_MSYPA")
  obj_fun_use <- mget(x = paste0("obj_", obj_fun), 
                      ifnotfound = FALSE)
  for (i in seq_along(obj_fun)) {
    assign(x = paste0("obj_", obj_fun[i]), obj_fun_use[[i]])
  }
  obj_desc <- obj_fun[unlist(obj_fun_use)]
  obj_desc <- paste0("obj_", paste0(obj_desc, collapse = "_"), collapse = "")
  
  ### store input data in temp file
  inp_file <- tempfile()
  saveRDS(object = input, file = inp_file, compress = FALSE)
  rm(input)
  gc()
  
  ### ------------------------------------------------------------------------ ###
  ### check if previous solutions can be used as suggestions ####
  ### ------------------------------------------------------------------------ ###
  
  ### years for summary statistics
  file_ext <- ifelse(stat_yrs == "all", "_res", 
                     paste0("_res_", stat_yrs))
  ### suffix if different risk limit used
  file_ext <- ifelse(isTRUE(!identical(risk_threshold, 0.05) & 
                              isTRUE(obj_ICES_MSYPA)), 
                     paste0(file_ext, "_", risk_threshold), 
                     file_ext)
  file_ext <- paste0(file_ext, ".rds")
  
  if (isTRUE(add_suggestions)) {
    ### find files
    avail <- list.files(path_out, pattern = paste0("--", obj_desc, file_ext))
    avail <- gsub(x = avail, pattern = paste0("--", obj_desc, file_ext),
                  replacement = "")
    avail <- strsplit(x = avail, split = "-")
    ### need to have fewer parameters
    avail <- avail[which(sapply(avail, length) < length(scn_pars))]
    ### skip parameters not used
    if (isTRUE(length(avail) > 0)) {
      avail <- avail[which(sapply(avail, function(x) all(x %in% scn_pars)))]
      ### if some parameters fixed, remove suggestions without them
      avail <- avail[which(sapply(avail, function(x) 
        all(paste0(par_fixed, val_fixed) %in% x)))]
      if (isTRUE(length(avail) > 0)) {
        ### load results
        res_add <- lapply(avail, function(x) {
          tmp <- readRDS(file = 
                           paste0(path_out, paste0(x, collapse = "-"), "--", obj_desc, "_res",
                                  ifelse(identical(stat_yrs, "all"), "", paste0("_", stat_yrs)), 
                                  ".rds"))
          tmp <- tmp@solution[1, ]
          if (is.na(tmp[which("upper_constraint" == names(tmp))])) {
            tmp[which("upper_constraint" == names(tmp))] <- Inf
          }
          return(tmp)
        })
        res_add <- do.call(rbind, res_add)
        if (isTRUE(nrow(res_add) > 1)) {
          res_add <- data.frame(res_add, stringsAsFactors = FALSE)
        } else {
          res_add <- data.frame(res_add, stringsAsFactors = FALSE)
        }
        cat("adding GA suggestions:\n")
        print(res_add)
        ### add to GA suggestions
        ga_suggestions <- rbind(ga_suggestions, res_add)
        ga_suggestions <- unique(ga_suggestions)
      }
    }
  }
  
  ### ---------------------------------------------------------------------- ###
  ### run MSE with GA ####
  ### ---------------------------------------------------------------------- ###
  
  ### set random seed for reproducibility
  registerDoRNG(123)
  set.seed(1)
  
  ### run GA
  system.time({
    res <- ga(type = "real-valued", fitness = mp_fitness, inp_file = inp_file,
              obj_SSB = obj_SSB, obj_F = obj_F, obj_C = obj_C, 
              obj_risk = obj_risk, obj_ICV = obj_ICV, obj_ICES_PA = obj_ICES_PA,
              obj_ICES_PA2 = obj_ICES_PA2, obj_ICES_MSYPA = obj_ICES_MSYPA,
              stat_yrs = stat_yrs, risk_threshold = risk_threshold,
              path = path_out, check_file = TRUE,
              scenario = scenario,
              suggestions = ga_suggestions, lower = ga_lower, upper = ga_upper,
              names = ga_names,
              maxiter = maxiter, popSize = popSize, run = run,
              monitor = TRUE, keepBest = TRUE, parallel = cl1, seed = 1)
  })
  
  # debugonce(mse_r)
  # mse_r(c(0, 1, 1, 0, 1, 1), input = input, path = path_out, check_file = TRUE,
  #       scenario = "SSB_idx_r")
  
  
  ### save result
  saveRDS(object = res, file = paste0(path_out, scn_pars_c, 
                                      "--", obj_desc, file_ext))
  
  ### ---------------------------------------------------------------------- ###
  ### collate runs ####
  ### ---------------------------------------------------------------------- ###
  
  if (isTRUE(collate)) {
    files <- list.files(path = path_out, pattern = "[0-9]*[0-9].rds",
                        full.names = FALSE)
    files <- files[grep(x = files, pattern = "--", invert = TRUE)]
    names(files) <- sapply(files, function(x) {
      sub(x = x, pattern = ".rds", replacement = "", fixed = TRUE)
    })
    scns <- lapply(files, function(x) {
      pars <- an(strsplit(sub(x = x, pattern = ".rds", replacement = "", fixed = TRUE), 
                          split = "_")[[1]])
      names(pars) <- ga_names
      ### only keep scenarios where requested parameters are changed
      if (!all(ga_default[pos_default] == pars[pos_default])) return(NULL)
      if (!all(val_fixed == pars[pos_fixed])) return(NULL)
      stats <- readRDS(paste0(path_out, x))
      list(pars = pars, stats = stats)
    })
    scns[sapply(scns, is.null)] <- NULL
    #scns <- scns[order(sapply(scns, "[[", "obj"), decreasing = TRUE)]
    saveRDS(scns, 
            file = paste0(path_out, scn_pars_c, "--", obj_desc, "_runs",
                          ifelse(identical(stat_yrs, "last10"), "_last10", ""), 
                          ".rds"))
  }
  
  ### other catch rules
} else {
  
  ### output path
  path_out <- paste0("output/", n_iter, "_", n_yrs, "/", catch_rule, "/",
                     fhist, "/")
  dir.create(path_out, recursive = TRUE)
  
  ### run MSE
  ### run MP for each list element
  res_mp <- lapply(input, function(x) {
    do.call(mp, x)
  })
  file_name <- paste0(stock, collapse = "_")
  if (isTRUE(catch_rule == "hr")) {
    idx_quant <- i
    file_name <- paste0("int-", input[[1]]$ctrl$hcr@args$interval, "_",
                        "mult-", input[[1]]$ctrl$phcr@args$rate, "_",
                        file_name)
  }
  saveRDS(res_mp, paste0(path_out, file_name, "_mp.rds"))
  
  ### stats
  stats <- mp_stats(input = input, res_mp = res_mp, stat_yrs = stat_yrs)
  saveRDS(stats, paste0(path_out, file_name, "_stats", 
                        ifelse(stat_yrs == "last10", "_last10", ""), ".rds"))
  
}

### ------------------------------------------------------------------------ ###
### quit ####
### ------------------------------------------------------------------------ ###

quit(save = "no")

