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
  if (!exists("scenario")) scenario <- "multiplier"
  if (!exists("MP")) MP <- "rfb"
  ### rfb rule
  if (isTRUE(MP == "rfb")) {
    if (!exists("comp_r")) comp_r <- TRUE
    if (!exists("comp_f")) comp_f <- TRUE
    if (!exists("comp_b")) comp_b <- TRUE
    if (!exists("cap_below_b")) cap_below_b <- FALSE
    #if (!exists("upper_constraint")) upper_constraint <- 1.2
    #if (!exists("lower_constraint")) lower_constraint <- 0.7
  }
  ### OM
  if (!exists("stock_id")) stock_id <- "ple.27.7e"
  if (!exists("OM")) OM <- "baseline"
  ### GA search
  if (!exists("ga_search")) ga_search <- TRUE
  if (isTRUE(ga_search)) {
    if (!exists("popSize")) stop("popSize missing")
    if (!exists("maxiter")) stop("maxiter missing")
    if (!exists("run")) run <- maxiter
    if (!exists("collate")) collate <- TRUE
    ### objective function elements
    if (!exists("obj_fun")) obj_fun <- "ICES"
    if (!exists("obj_yrs")) obj_yrs <- "all"
    if (!exists("risk_threshold")) risk_threshold <- 0.05
    ### GA
    if (!exists("add_suggestions")) add_suggestions <- TRUE
    if (!exists("stat_yrs")) stat_yrs <- "multiple"
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
req_pckgs <- c("FLCore", "FLash", "FLBRP", "mse", "FLfse",
               "GA", "doParallel", "doRNG",
               "tidyr", "dplyr", "stockassessment")
for (i in req_pckgs) library(package = i, character.only = TRUE)

### load additional functions
source("funs.R")
source("funs_GA.R")
source("funs_WKNSMSE.R")

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
    source("funs.R")
    source("funs_GA.R")
    source("funs_WKNSMSE.R")
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
        source("funs_WKNSMSE.R", echo = FALSE)
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
      source("funs_WKNSMSE.R", echo = FALSE)
    }
  } else {
    cl1 <- FALSE
  }
}

### ------------------------------------------------------------------------ ###
### load OM ####
### ------------------------------------------------------------------------ ###

input <- readRDS(paste0("input/", stock_id, "/", OM, "/", n_iter, "_", n_yrs,
                        "/input_", MP, ".rds"))

### ------------------------------------------------------------------------ ###
### specify scenario ####
### ------------------------------------------------------------------------ ###

### modify MP parameters if requested
### rfb rule
if (isTRUE(MP == "rfb")) {
  ### catch rule components
  input$ctrl$est@args$comp_r <- comp_r
  input$ctrl$est@args$comp_f <- comp_f
  input$ctrl$est@args$comp_b <- comp_b
  ### turn of uncertainty cap when index below Itrigger?
  input$ctrl$isys@args$cap_below_b <- cap_below_b
  
### 2 over 3 rule
} else if (isTRUE(MP == "2over3")) {
  ### nothing to do
} else if (isTRUE(MP == "hr")) {
  ### to do
}

### within scenario parallelisation?
if (isTRUE(n_workers > 1) & isTRUE(n_blocks > 1)) {
  input$args$nblocks <- n_blocks
}

### ------------------------------------------------------------------------ ###
### GA set-up ####
### ------------------------------------------------------------------------ ###
if (isTRUE(MP == "rfb") & isTRUE(ga_search)) {
  
  ### GA arguments
  ga_names <- c("lag_idx", "range_idx_1", "range_idx_2", "range_catch",
                "exp_r", "exp_f", "exp_b", "interval", "multiplier",
                "upper_constraint", "lower_constraint")
  ga_suggestions <- rbind(
    c(0, 1, 1, 1, 1, 1, 1, 1, 1, 1.2, 0.7), ### most current data
    c(0, 1, 1, 1, 1, 1, 1, 2, 1, 1.2, 0.7),
    c(0, 1, 1, 1, 0, 0, 0, 2, 1, Inf, 0), ### constant catch
    c(0, 1, 1, 1, 0, 0, 0, 2, 1, 1.2, 0.7), ### constant catch (dummy cap)
    ### default, annual/biennial, turning off elements
    expand.grid(1, 2, 3, 1, 0:1, 0:1, 0:1, 1:2, 1, 1.2, 0.7),
    c(1, 2, 3, 1, 1, 1, 1, 2, 1, Inf, 0), ### default without cap
    c(1, 2, 3, 1, 1, 1, 1, 2, 0, Inf, 0), ### zero catch
    c(1, 2, 3, 1, 1, 1, 1, 2, 0, 1.2, 0.7) ### zero catch but capped
  )
  ga_default <- c(1, 2, 3, 1, 1, 1, 1, 2, 1, 1.2, 0.7)
  ga_lower <- c(0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 0)
  ga_upper <- c(1, 5, 5, 1, 2, 2, 2, 5, 2, 5, 1)
  ### turn of parameters not requested, i.e. limit to default value
  pos_default <- which(sapply(mget(ga_names, ifnotfound = FALSE), isFALSE))
  ga_lower[pos_default] <- ga_default[pos_default]
  ga_upper[pos_default] <- ga_default[pos_default]
  ### remove not requested parameters from suggestions
  ga_suggestions[, pos_default] <- rep(ga_default[pos_default], 
                                       each = nrow(ga_suggestions))
  ga_suggestions <- unique(ga_suggestions)
  ### fix parameters?
  pos_fixed <- which(sapply(mget(ga_names, ifnotfound = FALSE), is.numeric))
  if (isTRUE(length(pos_fixed) > 0)) {
    par_fixed <- names(pos_fixed)
    val_fixed <- mget(ga_names, ifnotfound = FALSE)[pos_fixed]
    par_fixed_single <- names(which(sapply(val_fixed, length) == 1))
    val_fixed_single <- val_fixed[par_fixed_single]
    if (isTRUE(length(par_fixed_single) > 0)) {
      pos_par_fixed_single <- match(par_fixed_single, ga_names)
      for (pos in seq_along(par_fixed_single)) {
        ga_suggestions[, pos_par_fixed_single[pos]] <- val_fixed_single[[pos]]
      }
      ga_suggestions <- unique(ga_suggestions)
    }
    par_fixed_multiple <- names(which(sapply(val_fixed, length) > 1))
    if (isTRUE(length(par_fixed_multiple) > 0)) {
      val_fixed_multiple <- val_fixed[par_fixed_multiple]
      pos_par_fixed_multiple <- match(par_fixed_multiple, ga_names)
      ga_suggestions[, pos_par_fixed_multiple] <- NA
      ga_suggestions <- unique(ga_suggestions)
      ga_suggestions <- ga_suggestions[rep(1, max(sapply(val_fixed_multiple, 
                                                         length))), ]
      for (pos in seq_along(par_fixed_multiple)) {
        ga_suggestions[, pos_par_fixed_multiple[pos]] <-
          val_fixed[[par_fixed_multiple[pos]]]
      }
      ga_suggestions <- unique(ga_suggestions)
    }
  } else {
    par_fixed_single <- par_fixed_multiple <- NULL
    val_fixed_single <- val_fixed_multiple <- NULL
  }
  ### finalise
  ga_suggestions <- unique(ga_suggestions)
  names(ga_suggestions) <- ga_names
  
  ### if only one parameter modified & fixed, 
  ### run all supplied values instead of GA search
  if (exists("par_fixed_multiple")) {
    if (isTRUE(length(par_fixed_multiple) > 0 )) {
      ### adapt GA dimensions
      maxiter <- run <- 1
      popSize <- nrow(ga_suggestions)
    }
  }
  
  ### ---------------------------------------------------------------------- ###
  ### paths ####
  ### ---------------------------------------------------------------------- ###
  
  ### output path
  ### set name depending on which GA parameters are used
  scn_pars <- ga_names[setdiff(seq_along(ga_names), pos_default)]
  if (isTRUE(length(pos_fixed) > 0)) {
    if (isTRUE(length(val_fixed_single) > 0)) {
      scn_pars[which(scn_pars %in% par_fixed_single)] <- paste0(
        scn_pars[which(scn_pars %in% par_fixed_single)], val_fixed_single)
    }
  }
  scn_pars_c <- paste0(scn_pars, collapse = "-")
  
  path_out <- paste0("output/", stock_id, "/", OM, "/", n_iter, "_", n_yrs, "/",
                     scenario, "/", MP, "/")
  dir.create(path_out, recursive = TRUE)
  
  ### objective function elements
  # obj
  obj_desc <- paste0("obj_", paste0(obj_fun, collapse = "_"), collapse = "")
  if (isTRUE(obj_yrs == "all")) obj_yrs <- paste0("1:", n_yrs)
  
  ### store input data in temp file
  inp_file <- tempfile()
  saveRDS(object = input, file = inp_file, compress = FALSE)
  rm(input)
  gc()
  
  ### ------------------------------------------------------------------------ ###
  ### check if previous solutions can be used as suggestions ####
  ### ------------------------------------------------------------------------ ###
  
  ### years for summary statistics
  file_ext <- ifelse(obj_yrs == "all", "_res",
                     paste0("_res_", 
                            gsub(x = obj_yrs, pattern = ":", replacement = "-")))
  ### suffix if different risk limit used
  file_ext <- ifelse(isTRUE(!identical(risk_threshold, 0.05) &
                              any(c("MSYPA", "ICES") %in% obj_fun)),
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
      if (isTRUE(length(par_fixed) > 0)) {
        if (isTRUE(length(par_fixed_single > 0))) {
          avail <- avail[which(sapply(avail, function(x)
            all(paste0(par_fixed_single, val_fixed_single) %in% x)))]
        }
      }
      if (isTRUE(length(avail) > 0)) {
        ### load results
        res_add <- lapply(avail, function(x) {
          tmp <- readRDS(file =
            paste0(path_out, paste0(x, collapse = "-"), "--", obj_desc, 
                   file_ext))
          tmp <- tmp@solution[1, ]
          if (is.na(tmp[which("upper_constraint" == names(tmp))])) {
            tmp[which("upper_constraint" == names(tmp))] <- Inf
          }
          return(tmp)
        })
        res_add <- do.call(rbind, res_add)
        res_add <- data.frame(res_add, stringsAsFactors = FALSE)
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
              obj_fun = obj_fun, obj_yrs = obj_yrs, stat_yrs = stat_yrs, 
              risk_threshold = risk_threshold,
              path = path_out, check_file = TRUE,
              scenario = scenario,
              suggestions = ga_suggestions, lower = ga_lower, upper = ga_upper,
              names = ga_names,
              maxiter = maxiter, popSize = popSize, run = run,
              monitor = TRUE, keepBest = TRUE, parallel = cl1, seed = 1)
  })
  
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
    scns <- lapply(files, function(x) {#browser()
      pars <- an(strsplit(sub(x = x, pattern = ".rds", replacement = "", 
                              fixed = TRUE),
                          split = "_")[[1]])
      names(pars) <- ga_names
      ### only keep scenarios where requested parameters are changed
      if (!all(ga_default[pos_default] == pars[pos_default])) return(NULL)
      if (isTRUE(length(pos_fixed) > 0)) {
        if (isTRUE(length(val_fixed_single) > 0)) {
          if (!all(val_fixed_single == pars[pos_par_fixed_single])) return(NULL)
        }
      }
      stats <- readRDS(paste0(path_out, x))
      list(pars = pars, stats = stats)
    })
    scns[sapply(scns, is.null)] <- NULL
    saveRDS(scns, 
            file = paste0(path_out, scn_pars_c, "--", obj_desc, "_runs", 
                          ".rds"))
  }
  
  ### other MPs
} else {
  
  ### output path
  path_out <- paste0("output/", stock_id, "/", OM, "/", n_iter, "_", n_yrs, "/",
                     scenario, "/", MP, "/")
  dir.create(path_out, recursive = TRUE)
  
  ### run MSE
  res_mp <- do.call(mp, input)
  
  file_name <- "mp"
  if (isTRUE(MP == "hr")) {
    idx_quant <- i
    file_name <- paste0("int-", input$ctrl$hcr@args$interval, "_",
                        "mult-", input$ctrl$phcr@args$rate, "_",
                        file_name)
  }
  saveRDS(res_mp, paste0(path_out, file_name, ".rds"))
  
  ### stats
  stats <- mp_stats(input = input, res_mp = res_mp, stat_yrs = stat_yrs)
  saveRDS(stats, paste0(path_out, "stats", ".rds"))
  
}

### ------------------------------------------------------------------------ ###
### quit ####
### ------------------------------------------------------------------------ ###

quit(save = "no")

