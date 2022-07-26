### ------------------------------------------------------------------------ ###
### script for running MSE locally (not on HPC) ####
### ------------------------------------------------------------------------ ###
suppressMessages(library(FLCore))
suppressMessages(library(FLasher))
suppressMessages(library(FLBRP))
suppressMessages(library(mse))
suppressMessages(library(FLfse))
suppressMessages(library(FLXSA))
suppressMessages(library(GA))
suppressMessages(library(doParallel))
suppressMessages(library(doRNG))
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library(stockassessment))
suppressMessages(library(doParallel))

req_pckgs <- c("FLCore", "FLasher", "FLBRP", "mse", "FLfse", "FLXSA",
               "GA", "doParallel", "doRNG",
               "tidyr", "dplyr", "stockassessment")
for (i in req_pckgs) 
  suppressMessages(library(package = i, character.only = TRUE))
source("funs.R")
source("funs_GA.R")
source("funs_WKNSMSE.R")
source("funs_OM.R")

cl1 <- FALSE

### ------------------------------------------------------------------------ ###
### harvest rate & rfb rule: all stocks & OMs - default & optimised ####
### ------------------------------------------------------------------------ ###
stocks <- c("ple.27.7e", "cod.27.47d20", "her.27.3a47d")
OMs_cod <- c("baseline", "rec_higher", "M_dd", "M_no_migration", "rec_failure")
OMs_ple <- c("baseline", "M_low", "M_high", "M_Gislason", "no_discards",
             "rec_no_AC", "rec_failure")
OMs_her <- c("baseline", "rec_higher", "rec_failure",
             "M_high", "M_low")

. <- foreach(stock = stocks) %:%
  foreach(OM = switch(stock,
    "ple.27.7e" = OMs_ple,
    "cod.27.47d20" = OMs_cod,
    "her.27.3a47d" = OMs_her
  )) %:%
  foreach(MP = c("rfb", "hr"))  %:%
  foreach(optimised = c("default", "multiplier", "all")) %do% {
    #browser()
    cat(paste0("stock=", stock, " - OM=", OM, " - MP=", MP, " - optimised=",
               optimised, "\n"))
    suppressWarnings(
      rm(lag_idx, range_idx_1, range_idx_2, range_catch, exp_r, exp_f, exp_b, 
       interval, multiplier, upper_constraint, lower_constraint,
       idxB_lag, idxB_range_3, comp_b_multiplier))
    ### change if recruitment failure included
    scenario <- ifelse(OM == "rec_failure", "rec_failure", "")
    rec_failure <- ifelse(OM == "rec_failure", TRUE, FALSE)
    if (isTRUE(rec_failure)) rec_failure <- 2021:2025
    OM <- ifelse(OM == "rec_failure", "baseline", OM)
    ### default multiplier
    if (identical(optimised, "default")) {
      if (identical(MP, "hr")) {
        multiplier <- 1
      } else if (identical(MP, "rfb")) {
        multiplier <- switch(stock,
                             "ple.27.7e" = 0.95,
                             "cod.27.47d20" = 0.95,
                             "her.27.3a47d" = 0.90)
      } else {
        stop()
      }
    } else if (identical(optimised, "multiplier")) {
      if (identical(stock, "cod.27.47d20")) {
        multiplier <- switch(MP, "rfb" = 1.73, "hr" = 0.83)
      } else if (identical(stock, "ple.27.7e")) {
        multiplier <- switch(MP, "rfb" = 1.16, "hr" = 2.46)
      } else if (identical(stock, "her.27.3a47d")) {
        multiplier <- switch(MP, "rfb" = 0.93, "hr" = 1.55)
      }
    } else if (identical(optimised, "all")) {
      if (identical(stock, "cod.27.47d20")) {
        if (identical(MP, "rfb")) {
          lag_idx <- 0.4453777
          range_idx_1 <- 3.912436
          range_idx_2 <- 3.415026
          range_catch <- 1
          exp_r <- 0.1190475
          exp_f <- 1.310675
          exp_b <- 0.3605948
          interval <- 4.180232
          multiplier <- 1.055067
        } else if (identical(MP, "hr")) {
          idxB_lag <- 0.3107843
          idxB_range_3 <- 1.497951
          exp_b <- 1
          comp_b_multiplier <- 1.046093
          interval <- 3.79667
          multiplier <- 1.66751
        }
      } else if (identical(stock, "ple.27.7e")) {
        if (identical(MP, "rfb")) {
          lag_idx <- 0.3889051
          range_idx_1 <- 4.742625
          range_idx_2 <- 4.48785
          range_catch <- 1
          exp_r <- 1.720633
          exp_f <- 1.712309
          exp_b <- 1.924345
          interval <- 2.205795
          multiplier <- 1.648888
        } else if (identical(MP, "hr")) {
          idxB_lag <- 0.5264083
          idxB_range_3 <- 2.318901
          exp_b <- 1
          comp_b_multiplier <- 0.8147290
          interval <- 1.826095
          multiplier <- 2.559064
        }
      } else if (identical(stock, "her.27.3a47d")) {
        if (identical(MP, "rfb")) {
          lag_idx <- 0.2586439
          range_idx_1 <- 2.026223
          range_idx_2 <- 2.784513
          range_catch <- 1
          exp_r <- 1.205622
          exp_f <- 1.464026
          exp_b <- 1.367814
          interval <- 2.867753
          multiplier <- 0.9355739
        } else if (identical(MP, "hr")) {
          idxB_lag <- 0.3070643
          idxB_range_3 <- 2.295637
          exp_b <- 1
          comp_b_multiplier <- 1.035724
          interval <- 1.388251
          multiplier <- 1.692025
        }
      }
    }
    ### define local arguments
    args_local <- c("ga_search=TRUE","n_blocks=1", "n_workers=0", 
                    paste0("scenario=\'", scenario, "\'"),
                    paste0("MP=\'", MP, "\'"),
                    "n_yrs=20", "check_file=FALSE",
                    paste0("stock_id=\'", stock, "\'"),
                    paste0("OM=\'", OM, "\'"),
                    "save_MP=TRUE", "popSize=1", "maxiter=1",
                    paste0("multiplier=", multiplier, ""),
                    "add_suggestions=FALSE", "collate=FALSE")
    if (identical(MP, "rfb") & identical(optimised, "all")) 
      args_local <- append(args_local, 
                           c(paste0("lag_idx=", lag_idx),
                            paste0("range_idx_1=", range_idx_1),
                            paste0("range_idx_2=", range_idx_2),
                            paste0("range_catch=", range_catch),
                            paste0("exp_r=", exp_r),
                            paste0("exp_f=", exp_f),
                            paste0("exp_b=", exp_b),
                            paste0("interval=", interval)))
    if (identical(MP, "hr") & identical(optimised, "all"))
      args_local <- append(args_local, 
                           c(paste0("idxB_lag=", idxB_lag),
                            paste0("idxB_range_3=", idxB_range_3),
                            paste0("exp_b=", exp_b),
                            paste0("comp_b_multiplier=", comp_b_multiplier),
                            paste0("interval=", interval)))
    source("MP_run.R")
}
