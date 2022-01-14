### ------------------------------------------------------------------------ ###
### script for running MSY estimations ####
### ------------------------------------------------------------------------ ###

### ------------------------------------------------------------------------ ###
### arguments ####
### ------------------------------------------------------------------------ ###
args <- commandArgs(TRUE)
print("arguments passed on to this script:")
print(args)

### extract arguments
for (i in seq_along(args)) eval(parse(text = args[[i]]))
### set default arguments
if (!exists("n_workers")) n_workers <- 10
if (!exists("stock_id")) stock_id <- "ple.27.7e"
if (!exists("OM")) OM <- "baseline"
#OM <- c("baseline", "M_low", "M_high", "M_Gislason", "rec_no_AC")
if (!exists("yr_start")) yr_start <- 2021
if (!exists("n_iter")) n_iter <- 1000
if (!exists("vals_ini")) vals_ini <- seq(0, 1, 0.1)
if (!exists("lower")) lower <- 0
if (!exists("upper")) upper <- 0.3
if (!exists("tol")) tol <- 0.001
if (!exists("plot")) plot <- TRUE
if (!exists("x_label")) x_label <- "F (ages 3-6)"
if (!exists("plot_res")) plot_res <- TRUE
if (!exists("save_res")) save_res <- TRUE


### ------------------------------------------------------------------------ ###
### prepare R session ####
### ------------------------------------------------------------------------ ###
req_pckgs <- c("FLCore", "FLasher", "FLBRP", "mse", "FLfse", "FLXSA",
               "doParallel", 
               "tidyr", "dplyr", "ggplot2")
for (i in req_pckgs) library(package = i, character.only = TRUE)

### load additional functions
req_scripts <- c("funs.R", "funs_GA.R", "funs_WKNSMSE.R", "funs_OM.R")
for (i in req_scripts) source(i)

### parallelisation
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
    for (i in req_scripts) source(i)
  }
}

### ------------------------------------------------------------------------ ###
### MSY reference values ####
### ------------------------------------------------------------------------ ###

### loop through OM scenarios
for (i in OM) {
  
  print(i)
  
  if (isTRUE(stock_id %in% c("ple.27.7e", "cod.27.47d20", "her.27.3a47d"))) {
    res <- est_MSY(stock_id = stock_id, OM = i,
                   yr_start = yr_start, n_blocks = n_workers, n_iter = n_iter,
                   vals_ini = vals_ini,
                   lower = lower, upper = upper, tol = tol,
                   plot = plot_res, x_label = x_label,
                   save = save_res)
    res$result[which.max(res$result$catch), ]
  #        Ftrgt   catch      ssb      tsb      rec
  # 15 0.1638845 1702.92 9536.053 11136.77 6542.729
  }
}

