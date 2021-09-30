### ------------------------------------------------------------------------ ###
### script for running MSE locally (not on HPC) ####
### ------------------------------------------------------------------------ ###

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
### rfb rule - optimised parameterisations ####
### ------------------------------------------------------------------------ ###

### PA - multiplier = 0.95
args_local <- c("n_blocks=10", "n_workers=0", "scenario=''", "MP='rfb'",
                "n_yrs=100",
                "stock_id='ple.27.7e'", "OM='baseline'", "save_MP=TRUE",
                "popSize=1", "maxiter=1",
                "multiplier=0.95", "add_suggestions=FALSE", "collate=FALSE")
source("MP_run.R")
### optimised - multiplier = 1.16
args_local <- c("n_blocks=10", "n_workers=0", "scenario=''", "MP='rfb'",
                "n_yrs=100",
                "stock_id='ple.27.7e'", "OM='baseline'", "save_MP=TRUE",
                "popSize=1", "maxiter=1",
                "multiplier=1.16", "add_suggestions=FALSE", "collate=FALSE")
source("MP_run.R")
### optimised - all parameters
args_local <- c("n_blocks=10", "n_workers=0", "scenario=''", "MP='rfb'",
                "n_yrs=100",
                "stock_id='ple.27.7e'", "OM='baseline'", "save_MP=TRUE",
                "popSize=1", "maxiter=1",
                "add_suggestions=FALSE", "collate=FALSE",
                "lag_idx=1", "range_idx_1=2.995726", "range_idx_2=2.140256",
                "range_catch=1", "exp_r=0.9387879", "exp_f=0.8256765", 
                "exp_b=1.190055", "interval=2.466707", "multiplier=1.077025",
                "upper_constraint=1.2", "lower_constraint=0.7")
source("MP_run.R")





### ------------------------------------------------------------------------ ###
### rfb rule & alternative OMs ####
### ------------------------------------------------------------------------ ###

### with multiplier = 1.16 (optimised)

### baseline OM
args_local <- c("n_blocks=10", "n_workers=0", "scenario=''", "MP='rfb'",
                "n_yrs=20", "check_file=FALSE",
                "stock_id='ple.27.7e'", "OM='baseline'", "save_MP=TRUE",
                "popSize=1", "maxiter=1",
                "multiplier=1.16", "add_suggestions=FALSE", "collate=FALSE")
source("MP_run.R")
### M high
args_local <- c("n_blocks=10", "n_workers=0", "scenario=''", "MP='rfb'",
                "n_yrs=20", "check_file=FALSE",
                "stock_id='ple.27.7e'", "OM='M_high'", "save_MP=TRUE",
                "popSize=1", "maxiter=1",
                "multiplier=1.16", "add_suggestions=FALSE", "collate=FALSE")
source("MP_run.R")
### M low
args_local <- c("n_blocks=10", "n_workers=0", "scenario=''", "MP='rfb'",
                "n_yrs=20", "check_file=FALSE",
                "stock_id='ple.27.7e'", "OM='M_low'", "save_MP=TRUE",
                "popSize=1", "maxiter=1",
                "multiplier=1.16", "add_suggestions=FALSE", "collate=FALSE")
source("MP_run.R")
### M Gislason
args_local <- c("n_blocks=10", "n_workers=0", "scenario=''", "MP='rfb'",
                "n_yrs=20", "check_file=FALSE",
                "stock_id='ple.27.7e'", "OM='M_Gislason'", "save_MP=TRUE",
                "popSize=1", "maxiter=1",
                "multiplier=1.16", "add_suggestions=FALSE", "collate=FALSE")
source("MP_run.R")
### baseline without recruitment auto-correlation
args_local <- c("n_blocks=10", "n_workers=0", "scenario=''", "MP='rfb'",
                "n_yrs=20", "check_file=FALSE",
                "stock_id='ple.27.7e'", "OM='rec_no_AC'", "save_MP=TRUE",
                "popSize=1", "maxiter=1",
                "multiplier=1.16", "add_suggestions=FALSE", "collate=FALSE")
source("MP_run.R")



### with multiplier = 0.95 (generic PA)

### baseline OM
args_local <- c("n_blocks=10", "n_workers=0", "scenario=''", "MP='rfb'",
                "n_yrs=20", "check_file=FALSE",
                "stock_id='ple.27.7e'", "OM='baseline'", "save_MP=TRUE",
                "popSize=1", "maxiter=1",
                "multiplier=0.95", "add_suggestions=FALSE", "collate=FALSE")
source("MP_run.R")
### M high
args_local <- c("n_blocks=10", "n_workers=0", "scenario=''", "MP='rfb'",
                "n_yrs=20", "check_file=FALSE",
                "stock_id='ple.27.7e'", "OM='M_high'", "save_MP=TRUE",
                "popSize=1", "maxiter=1",
                "multiplier=0.95", "add_suggestions=FALSE", "collate=FALSE")
source("MP_run.R")
### M low
args_local <- c("n_blocks=10", "n_workers=0", "scenario=''", "MP='rfb'",
                "n_yrs=20", "check_file=FALSE",
                "stock_id='ple.27.7e'", "OM='M_low'", "save_MP=TRUE",
                "popSize=1", "maxiter=1",
                "multiplier=0.95", "add_suggestions=FALSE", "collate=FALSE")
source("MP_run.R")
### M Gislason
args_local <- c("n_blocks=10", "n_workers=0", "scenario=''", "MP='rfb'",
                "n_yrs=20", "check_file=FALSE",
                "stock_id='ple.27.7e'", "OM='M_Gislason'", "save_MP=TRUE",
                "popSize=1", "maxiter=1",
                "multiplier=0.95", "add_suggestions=FALSE", "collate=FALSE")
source("MP_run.R")
### baseline without recruitment auto-correlation
args_local <- c("n_blocks=10", "n_workers=0", "scenario=''", "MP='rfb'",
                "n_yrs=20", "check_file=FALSE",
                "stock_id='ple.27.7e'", "OM='rec_no_AC'", "save_MP=TRUE",
                "popSize=1", "maxiter=1",
                "multiplier=0.95", "add_suggestions=FALSE", "collate=FALSE")
source("MP_run.R")




