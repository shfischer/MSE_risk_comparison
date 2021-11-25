### ------------------------------------------------------------------------ ###
### script for running MSE locally (not on HPC) ####
### ------------------------------------------------------------------------ ###
suppressMessages(library(FLCore))
suppressMessages(library(FLash))
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

req_pckgs <- c("FLCore", "FLash", "FLBRP", "mse", "FLfse", "FLXSA",
               "GA", "doParallel", "doRNG",
               "tidyr", "dplyr", "stockassessment")
for (i in req_pckgs) 
  suppressMessages(library(package = i, character.only = TRUE))
source("funs.R")
source("funs_GA.R")
source("funs_WKNSMSE.R")
source("funs_OM.R")
library(doParallel)
cl <- makeCluster(10)
registerDoParallel(cl)
cl_length <- length(cl)
. <- foreach(i = seq(cl_length)) %dopar% {
  for (i in req_pckgs) library(package = i, character.only = TRUE,
                               warn.conflicts = FALSE, verbose = FALSE,
                               quietly = TRUE)
  source("funs.R", echo = FALSE)
  source("funs_GA.R", echo = FALSE)
  source("funs_WKNSMSE.R", echo = FALSE)
}

### ------------------------------------------------------------------------ ###
### ple.27.7e: rfb rule - optimised parameterisations ####
### ------------------------------------------------------------------------ ###

### PA - multiplier = 0.95
args_local <- c("n_blocks=10", "n_workers=0", "scenario=''", "MP='rfb'",
                "n_yrs=20", "check_file=FALSE",
                "stock_id='ple.27.7e'", "OM='baseline'", "save_MP=TRUE",
                "popSize=1", "maxiter=1",
                "multiplier=0.95", "add_suggestions=FALSE", "collate=FALSE")
source("MP_run.R")
### optimised - multiplier = 1.16
args_local <- c("n_blocks=10", "n_workers=0", "scenario=''", "MP='rfb'",
                "n_yrs=20", "check_file=FALSE",
                "stock_id='ple.27.7e'", "OM='baseline'", "save_MP=TRUE",
                "popSize=1", "maxiter=1",
                "multiplier=1.16", "add_suggestions=FALSE", "collate=FALSE")
source("MP_run.R")
### optimised - all parameters
args_local <- c("n_blocks=10", "n_workers=0", "scenario=''", "MP='rfb'",
                "n_yrs=20", "check_file=FALSE",
                "stock_id='ple.27.7e'", "OM='baseline'", "save_MP=TRUE",
                "popSize=1", "maxiter=1",
                "add_suggestions=FALSE", "collate=FALSE",
                "lag_idx=1", "range_idx_1=2.593756", "range_idx_2=2.491983",
                "range_catch=1", "exp_r=1.002068", "exp_f=0.3327127", 
                "exp_b=1.084681", "interval=2.075026", "multiplier=1.095075",
                "upper_constraint=1.2", "lower_constraint=0.7")
source("MP_run.R")





### ------------------------------------------------------------------------ ###
### ple.27.7e: rfb rule & alternative OMs ####
### ------------------------------------------------------------------------ ###

### with multiplier = 1.16 (optimised)

### M high
rm(lag_idx, range_idx_1, range_idx_2, range_catch, exp_r, exp_f, exp_b, 
   interval, multiplier, upper_constraint, lower_constraint)
args_local <- c("n_blocks=10", "n_workers=0", "scenario=''", "MP='rfb'",
                "n_yrs=20", "check_file=FALSE",
                "stock_id='ple.27.7e'", "OM='M_high'", "save_MP=TRUE",
                "popSize=1", "maxiter=1",
                "multiplier=1.16", "add_suggestions=FALSE", "collate=FALSE")
source("MP_run.R")
### M low
rm(lag_idx, range_idx_1, range_idx_2, range_catch, exp_r, exp_f, exp_b, 
   interval, multiplier, upper_constraint, lower_constraint)
args_local <- c("n_blocks=10", "n_workers=0", "scenario=''", "MP='rfb'",
                "n_yrs=20", "check_file=FALSE",
                "stock_id='ple.27.7e'", "OM='M_low'", "save_MP=TRUE",
                "popSize=1", "maxiter=1",
                "multiplier=1.16", "add_suggestions=FALSE", "collate=FALSE")
source("MP_run.R")
### M Gislason
rm(lag_idx, range_idx_1, range_idx_2, range_catch, exp_r, exp_f, exp_b, 
   interval, multiplier, upper_constraint, lower_constraint)
args_local <- c("n_blocks=10", "n_workers=0", "scenario=''", "MP='rfb'",
                "n_yrs=20", "check_file=FALSE",
                "stock_id='ple.27.7e'", "OM='M_Gislason'", "save_MP=TRUE",
                "popSize=1", "maxiter=1",
                "multiplier=1.16", "add_suggestions=FALSE", "collate=FALSE")
source("MP_run.R")
### baseline without recruitment auto-correlation
rm(lag_idx, range_idx_1, range_idx_2, range_catch, exp_r, exp_f, exp_b, 
   interval, multiplier, upper_constraint, lower_constraint)
args_local <- c("n_blocks=10", "n_workers=0", "scenario=''", "MP='rfb'",
                "n_yrs=20", "check_file=FALSE",
                "stock_id='ple.27.7e'", "OM='rec_no_AC'", "save_MP=TRUE",
                "popSize=1", "maxiter=1",
                "multiplier=1.16", "add_suggestions=FALSE", "collate=FALSE")
source("MP_run.R")
### no discards
rm(lag_idx, range_idx_1, range_idx_2, range_catch, exp_r, exp_f, exp_b, 
   interval, multiplier, upper_constraint, lower_constraint)
args_local <- c("n_blocks=10", "n_workers=0", "scenario=''", "MP='rfb'",
                "n_yrs=20", "check_file=FALSE",
                "stock_id='ple.27.7e'", "OM='no_discards'", "save_MP=TRUE",
                "popSize=1", "maxiter=1",
                "multiplier=1.16", "add_suggestions=FALSE", "collate=FALSE")
source("MP_run.R")
### recruitment failure
rm(lag_idx, range_idx_1, range_idx_2, range_catch, exp_r, exp_f, exp_b, 
   interval, multiplier, upper_constraint, lower_constraint)
args_local <- c("n_blocks=10", "n_workers=0", "scenario='rec_failure'", 
                "MP='rfb'",
                "n_yrs=20", "check_file=FALSE",
                "stock_id='ple.27.7e'", "OM='baseline'", "save_MP=TRUE",
                "popSize=1", "maxiter=1",
                "multiplier=1.16", "add_suggestions=FALSE", "collate=FALSE",
                "rec_failure=2021:2025")
source("MP_run.R")


### with multiplier = 0.95 (generic PA)

### M high
rm(lag_idx, range_idx_1, range_idx_2, range_catch, exp_r, exp_f, exp_b, 
   interval, multiplier, upper_constraint, lower_constraint)
args_local <- c("n_blocks=10", "n_workers=0", "scenario=''", "MP='rfb'",
                "n_yrs=20", "check_file=FALSE",
                "stock_id='ple.27.7e'", "OM='M_high'", "save_MP=TRUE",
                "popSize=1", "maxiter=1",
                "multiplier=0.95", "add_suggestions=FALSE", "collate=FALSE")
source("MP_run.R")
### M low
rm(lag_idx, range_idx_1, range_idx_2, range_catch, exp_r, exp_f, exp_b, 
   interval, multiplier, upper_constraint, lower_constraint)
args_local <- c("n_blocks=10", "n_workers=0", "scenario=''", "MP='rfb'",
                "n_yrs=20", "check_file=FALSE",
                "stock_id='ple.27.7e'", "OM='M_low'", "save_MP=TRUE",
                "popSize=1", "maxiter=1",
                "multiplier=0.95", "add_suggestions=FALSE", "collate=FALSE")
source("MP_run.R")
### M Gislason
rm(lag_idx, range_idx_1, range_idx_2, range_catch, exp_r, exp_f, exp_b, 
   interval, multiplier, upper_constraint, lower_constraint)
args_local <- c("n_blocks=10", "n_workers=0", "scenario=''", "MP='rfb'",
                "n_yrs=20", "check_file=FALSE",
                "stock_id='ple.27.7e'", "OM='M_Gislason'", "save_MP=TRUE",
                "popSize=1", "maxiter=1",
                "multiplier=0.95", "add_suggestions=FALSE", "collate=FALSE")
source("MP_run.R")
### baseline without recruitment auto-correlation
rm(lag_idx, range_idx_1, range_idx_2, range_catch, exp_r, exp_f, exp_b, 
   interval, multiplier, upper_constraint, lower_constraint)
args_local <- c("n_blocks=10", "n_workers=0", "scenario=''", "MP='rfb'",
                "n_yrs=20", "check_file=FALSE",
                "stock_id='ple.27.7e'", "OM='rec_no_AC'", "save_MP=TRUE",
                "popSize=1", "maxiter=1",
                "multiplier=0.95", "add_suggestions=FALSE", "collate=FALSE")
source("MP_run.R")
### no discards
rm(lag_idx, range_idx_1, range_idx_2, range_catch, exp_r, exp_f, exp_b, 
   interval, multiplier, upper_constraint, lower_constraint)
args_local <- c("n_blocks=10", "n_workers=0", "scenario=''", "MP='rfb'",
                "n_yrs=20", "check_file=FALSE",
                "stock_id='ple.27.7e'", "OM='no_discards'", "save_MP=TRUE",
                "popSize=1", "maxiter=1",
                "multiplier=0.95", "add_suggestions=FALSE", "collate=FALSE")
source("MP_run.R")
### recruitment failure
rm(lag_idx, range_idx_1, range_idx_2, range_catch, exp_r, exp_f, exp_b, 
   interval, multiplier, upper_constraint, lower_constraint)
args_local <- c("n_blocks=10", "n_workers=0", "scenario='rec_failure'", 
                "MP='rfb'",
                "n_yrs=20", "check_file=FALSE",
                "stock_id='ple.27.7e'", "OM='baseline'", "save_MP=TRUE",
                "popSize=1", "maxiter=1",
                "multiplier=0.95", "add_suggestions=FALSE", "collate=FALSE",
                "rec_failure=2021:2025")
source("MP_run.R")

### all parameters
OMs <- c("M_low", "M_high", "M_Gislason", "rec_no_AC", 
         "no_discards", "rec_failure")
for (OM in OMs) {
   rm(lag_idx, range_idx_1, range_idx_2, range_catch, exp_r, exp_f, exp_b, 
      interval, multiplier, upper_constraint, lower_constraint)
   rec_failure <- ifelse(isTRUE(OM == "rec_failure"), "2021:2025", "FALSE")
   OM <- ifelse(isTRUE(OM == "rec_failure"), "baseline", OM)
   scenario <- ifelse(isTRUE(OM == "rec_failure"), "'rec_failure'", "''")
   args_local <- c("n_blocks=10", "n_workers=0", 
                   paste0("scenario=", scenario), "MP='rfb'",
                   "n_yrs=20", "check_file=FALSE",
                   "stock_id='ple.27.7e'", paste0("OM='", OM, "'"),
                   "save_MP=TRUE",
                   "popSize=1", "maxiter=1",
                   "add_suggestions=FALSE", "collate=FALSE",
                   "lag_idx=1", "range_idx_1=2.593756", "range_idx_2=2.491983",
                   "range_catch=1", "exp_r=1.002068", "exp_f=0.3327127", 
                   "exp_b=1.084681", "interval=2.075026", "multiplier=1.095075",
                   "upper_constraint=1.2", "lower_constraint=0.7",
                   paste0("rec_failure=", rec_failure))
   source("MP_run.R")
}


### ------------------------------------------------------------------------ ###
### ple.27.7e: 2 over 3 rule ####
### ------------------------------------------------------------------------ ###
### baseline OM
args_local <- c("n_blocks=10", "n_workers=0", "scenario=''", 
                "MP='2over3'",
                "n_yrs=20", "check_file=FALSE",
                "stock_id='ple.27.7e'", "OM='baseline'", "save_MP=TRUE",
                "ga_search=FALSE")
source("MP_run.R")

### alternative OMs
OMs <- c("baseline", "M_low", "M_high", "M_Gislason", "rec_no_AC", 
         "no_discards", "rec_failure")
for (OM in OMs) {
   rm(lag_idx, range_idx_1, range_idx_2, range_catch, exp_r, exp_f, exp_b, 
      interval, multiplier, upper_constraint, lower_constraint)
   rec_failure <- ifelse(isTRUE(OM == "rec_failure"), "2021:2025", "FALSE")
   OM <- ifelse(isTRUE(OM == "rec_failure"), "baseline", OM)
   scenario <- ifelse(isTRUE(OM == "rec_failure"), "'rec_failure'", "''")
   args_local <- c("n_blocks=10", "n_workers=0", 
                   paste0("scenario=", scenario), 
                   "MP='2over3'",
                   "n_yrs=20", "check_file=FALSE",
                   "stock_id='ple.27.7e'", 
                   paste0("OM='", OM, "'"),
                   paste0("rec_failure=", rec_failure), 
                   "save_MP=TRUE",
                   "ga_search=FALSE")
   source("MP_run.R")
}

### ------------------------------------------------------------------------ ###
### cod.27.47d20: rfb rule - optimised parameterisations ####
### ------------------------------------------------------------------------ ###

### PA - multiplier = 0.95
args_local <- c("n_blocks=10", "n_workers=0", "scenario=''", "MP='rfb'",
                "n_yrs=20", "check_file=FALSE",
                "stock_id='cod.27.47d20'", "OM='baseline'", "save_MP=TRUE",
                "popSize=1", "maxiter=1",
                "multiplier=0.95", "add_suggestions=FALSE", "collate=FALSE")
source("MP_run.R")
### optimised - multiplier = 1.63
args_local <- c("n_blocks=10", "n_workers=0", "scenario=''", "MP='rfb'",
                "n_yrs=20", "check_file=FALSE",
                "stock_id='cod.27.47d20'", "OM='baseline'", "save_MP=TRUE",
                "popSize=1", "maxiter=1",
                "multiplier=1.63", "add_suggestions=FALSE", "collate=FALSE")
source("MP_run.R")
### optimised - all parameters
args_local <- c("n_blocks=10", "n_workers=0", "scenario=''", "MP='rfb'",
                "n_yrs=20", "check_file=FALSE",
                "stock_id='cod.27.47d20'", "OM='baseline'", "save_MP=TRUE",
                "popSize=1", "maxiter=1",
                "add_suggestions=FALSE", "collate=FALSE",
                "lag_idx=0.7046149", "range_idx_1=2.635445", 
                "range_idx_2=3.240444",
                "range_catch=1", "exp_r=0.1971575", "exp_f=0.6958902", 
                "exp_b=1.129193", "interval=4.086768", "multiplier=1.176044",
                "upper_constraint=1.2", "lower_constraint=0.7")
source("MP_run.R")

rm(lag_idx, range_idx_1, range_idx_2, range_catch, exp_r, exp_f, exp_b, 
   interval, multiplier, upper_constraint, lower_constraint)

### ------------------------------------------------------------------------ ###
### cod.27.47d20: rfb rule - alternative OMs ####
### ------------------------------------------------------------------------ ###

### PA - multiplier = 0.95
### OM: rec_higher
rm(lag_idx, range_idx_1, range_idx_2, range_catch, exp_r, exp_f, exp_b, 
   interval, multiplier, upper_constraint, lower_constraint, rec_failure)
args_local <- c("n_blocks=10", "n_workers=0", "scenario=''", "MP='rfb'",
                "n_yrs=20", "check_file=FALSE",
                "stock_id='cod.27.47d20'", "OM='rec_higher'", "save_MP=TRUE",
                "popSize=1", "maxiter=1",
                "multiplier=0.95", "add_suggestions=FALSE", "collate=FALSE")
source("MP_run.R")
### OM: M_dd
rm(lag_idx, range_idx_1, range_idx_2, range_catch, exp_r, exp_f, exp_b, 
   interval, multiplier, upper_constraint, lower_constraint, rec_failure)
args_local <- c("n_blocks=10", "n_workers=0", "scenario=''", "MP='rfb'",
                "n_yrs=20", "check_file=FALSE",
                "stock_id='cod.27.47d20'", "OM='M_dd'", "save_MP=TRUE",
                "popSize=1", "maxiter=1",
                "multiplier=0.95", "add_suggestions=FALSE", "collate=FALSE")
source("MP_run.R")
### OM: M_no_migration
rm(lag_idx, range_idx_1, range_idx_2, range_catch, exp_r, exp_f, exp_b, 
   interval, multiplier, upper_constraint, lower_constraint, rec_failure)
args_local <- c("n_blocks=10", "n_workers=0", "scenario=''", "MP='rfb'",
                "n_yrs=20", "check_file=FALSE",
                "stock_id='cod.27.47d20'", "OM='M_no_migration'", "save_MP=TRUE",
                "popSize=1", "maxiter=1",
                "multiplier=0.95", "add_suggestions=FALSE", "collate=FALSE")
source("MP_run.R")
### OM: rec_failure
rm(lag_idx, range_idx_1, range_idx_2, range_catch, exp_r, exp_f, exp_b, 
   interval, multiplier, upper_constraint, lower_constraint, rec_failure)
args_local <- c("n_blocks=10", "n_workers=0", "scenario='rec_failure'", 
                "MP='rfb'", "rec_failure=2021:2025",
                "n_yrs=20", "check_file=FALSE",
                "stock_id='cod.27.47d20'", "OM='baseline'", "save_MP=TRUE",
                "popSize=1", "maxiter=1",
                "multiplier=0.95", "add_suggestions=FALSE", "collate=FALSE")
source("MP_run.R")

### optimised multiplier = 1.63
### OM: rec_higher
rm(lag_idx, range_idx_1, range_idx_2, range_catch, exp_r, exp_f, exp_b, 
   interval, multiplier, upper_constraint, lower_constraint, rec_failure)
args_local <- c("n_blocks=10", "n_workers=0", "scenario=''", "MP='rfb'",
                "n_yrs=20", "check_file=FALSE",
                "stock_id='cod.27.47d20'", "OM='rec_higher'", "save_MP=TRUE",
                "popSize=1", "maxiter=1",
                "multiplier=1.63", "add_suggestions=FALSE", "collate=FALSE")
source("MP_run.R")
### OM: M_dd
rm(lag_idx, range_idx_1, range_idx_2, range_catch, exp_r, exp_f, exp_b, 
   interval, multiplier, upper_constraint, lower_constraint, rec_failure)
args_local <- c("n_blocks=10", "n_workers=0", "scenario=''", "MP='rfb'",
                "n_yrs=20", "check_file=FALSE",
                "stock_id='cod.27.47d20'", "OM='M_dd'", "save_MP=TRUE",
                "popSize=1", "maxiter=1",
                "multiplier=1.63", "add_suggestions=FALSE", "collate=FALSE")
source("MP_run.R")
### OM: M_no_migration
rm(lag_idx, range_idx_1, range_idx_2, range_catch, exp_r, exp_f, exp_b, 
   interval, multiplier, upper_constraint, lower_constraint, rec_failure)
args_local <- c("n_blocks=10", "n_workers=0", "scenario=''", "MP='rfb'",
                "n_yrs=20", "check_file=FALSE",
                "stock_id='cod.27.47d20'", "OM='M_no_migration'", "save_MP=TRUE",
                "popSize=1", "maxiter=1",
                "multiplier=1.63", "add_suggestions=FALSE", "collate=FALSE")
source("MP_run.R")
### OM: rec_failure
rm(lag_idx, range_idx_1, range_idx_2, range_catch, exp_r, exp_f, exp_b, 
   interval, multiplier, upper_constraint, lower_constraint, rec_failure)
args_local <- c("n_blocks=10", "n_workers=0", "scenario='rec_failure'", 
                "MP='rfb'", "rec_failure=2021:2025",
                "n_yrs=20", "check_file=FALSE",
                "stock_id='cod.27.47d20'", "OM='baseline'", "save_MP=TRUE",
                "popSize=1", "maxiter=1",
                "multiplier=1.63", "add_suggestions=FALSE", "collate=FALSE")
source("MP_run.R")

### all parameters optimised
### OM: rec_higher
rm(lag_idx, range_idx_1, range_idx_2, range_catch, exp_r, exp_f, exp_b, 
   interval, multiplier, upper_constraint, lower_constraint, rec_failure)
args_local <- c("n_blocks=10", "n_workers=0", "scenario=''", "MP='rfb'",
                "n_yrs=20", "check_file=FALSE",
                "stock_id='cod.27.47d20'", "OM='rec_higher'", "save_MP=TRUE",
                "popSize=1", "maxiter=1",
                "add_suggestions=FALSE", "collate=FALSE",
                "lag_idx=0.7046149", "range_idx_1=2.635445", 
                "range_idx_2=3.240444",
                "range_catch=1", "exp_r=0.1971575", "exp_f=0.6958902", 
                "exp_b=1.129193", "interval=4.086768", "multiplier=1.176044",
                "upper_constraint=1.2", "lower_constraint=0.7")
source("MP_run.R")
### OM: M_dd
rm(lag_idx, range_idx_1, range_idx_2, range_catch, exp_r, exp_f, exp_b, 
   interval, multiplier, upper_constraint, lower_constraint, rec_failure)
args_local <- c("n_blocks=10", "n_workers=0", "scenario=''", "MP='rfb'",
                "n_yrs=20", "check_file=FALSE",
                "stock_id='cod.27.47d20'", "OM='M_dd'", "save_MP=TRUE",
                "popSize=1", "maxiter=1",
                "add_suggestions=FALSE", "collate=FALSE",
                "lag_idx=0.7046149", "range_idx_1=2.635445", 
                "range_idx_2=3.240444",
                "range_catch=1", "exp_r=0.1971575", "exp_f=0.6958902", 
                "exp_b=1.129193", "interval=4.086768", "multiplier=1.176044",
                "upper_constraint=1.2", "lower_constraint=0.7")
source("MP_run.R")
### OM: M_no_migration
rm(lag_idx, range_idx_1, range_idx_2, range_catch, exp_r, exp_f, exp_b, 
   interval, multiplier, upper_constraint, lower_constraint, rec_failure)
args_local <- c("n_blocks=10", "n_workers=0", "scenario=''", "MP='rfb'",
                "n_yrs=20", "check_file=FALSE",
                "stock_id='cod.27.47d20'", "OM='M_no_migration'", "save_MP=TRUE",
                "popSize=1", "maxiter=1",
                "add_suggestions=FALSE", "collate=FALSE",
                "lag_idx=0.7046149", "range_idx_1=2.635445", 
                "range_idx_2=3.240444",
                "range_catch=1", "exp_r=0.1971575", "exp_f=0.6958902", 
                "exp_b=1.129193", "interval=4.086768", "multiplier=1.176044",
                "upper_constraint=1.2", "lower_constraint=0.7")
source("MP_run.R")
### OM: rec_failure
rm(lag_idx, range_idx_1, range_idx_2, range_catch, exp_r, exp_f, exp_b, 
   interval, multiplier, upper_constraint, lower_constraint, rec_failure)
args_local <- c("n_blocks=10", "n_workers=0", "scenario='rec_failure'", 
                "MP='rfb'", "rec_failure=2021:2025",
                "n_yrs=20", "check_file=FALSE",
                "stock_id='cod.27.47d20'", "OM='baseline'", "save_MP=TRUE",
                "popSize=1", "maxiter=1",
                "add_suggestions=FALSE", "collate=FALSE",
                "lag_idx=0.7046149", "range_idx_1=2.635445", 
                "range_idx_2=3.240444",
                "range_catch=1", "exp_r=0.1971575", "exp_f=0.6958902", 
                "exp_b=1.129193", "interval=4.086768", "multiplier=1.176044",
                "upper_constraint=1.2", "lower_constraint=0.7")
source("MP_run.R")










### ------------------------------------------------------------------------ ###
### tmp - cod - SAM tests ####
### ------------------------------------------------------------------------ ###

if (FALSE) {
  args_local <- c("n_iter=1000", "n_blocks=1000", "n_workers=10", "scenario=''", 
                  "MP='ICES_SAM'",
                  "n_yrs=20", "check_file=FALSE", "ga_search=FALSE",
                  "stock_id='cod.27.47d20'", "OM='baseline'", "save_MP=TRUE",
                  "stat_yrs='multiple'", "collate=FALSE")
  source("MP_run.R")
}
