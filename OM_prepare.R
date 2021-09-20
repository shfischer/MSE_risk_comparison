### ------------------------------------------------------------------------ ###
### create input for mp() ####
### ------------------------------------------------------------------------ ###
### this function loads the elements required for running mse::mp(),
### adapts them if necessary (dimensions), sets the MP and 
### creates the input object for mp()

input_mp <- function(stock_id = "ple.27.7e", OM = "baseline", n_iter = 1000,
                     n_yrs = 100, yr_start = 2021, iy = yr_start - 1,
                     n_blocks = 1, seed = 1, cut_hist = TRUE, MP = "rfb") {
  
  ### path to input objects
  path_input <- paste0("input/", stock_id, "/", OM, "/1000_100/")
  
  ### load objects
  ### use full dimensions (100 years, 1000 iterations) - reduced later
  stk_fwd <- readRDS(paste0(path_input, "stk.rds"))
  sr <- readRDS(paste0(path_input, "sr.rds"))
  idx <- readRDS(paste0(path_input, "idx.rds"))
  idx_dev <- readRDS(paste0(path_input, "idx_dev.rds"))
  catch_res <- readRDS(paste0(path_input, "catch_res.rds"))
  proc_res <- readRDS(paste0(path_input, "proc_res.rds"))
  stk_oem <- readRDS(paste0(path_input, "stk_oem.rds"))
  SAM_pars_ini <- readRDS(paste0(path_input, "SAM_initial.rds"))
  SAM_conf <- readRDS(paste0(path_input, "SAM_conf.rds"))
  ALKs <- readRDS(paste0(path_input, "ALKs.rds"))
  refpts_mse <- readRDS(paste0(path_input, "refpts_mse.rds"))
  
  ### find last year
  yr_end <- yr_start + n_yrs - 1
  ### reduce dimensions, if requested
  if (isTRUE(n_yrs < 100)) {
    ### OM stock
    stk_fwd <- window(stk_fwd, end = yr_end)
    ### OEM stock
    stk_oem <- window(stk_oem, end = yr_end)
    ### stock-recruitment model - need to go through elements individually
    rec(sr) <- window(rec(sr), end = yr_end)
    ssb(sr) <- window(ssb(sr), end = yr_end)
    residuals(sr) <- window(residuals(sr), end = yr_end)
    fitted(sr) <- window(fitted(sr), end = yr_end)
    range(sr)[["maxyear"]] <- yr_end
    ### index
    idx <- window(idx, end = yr_end)
    ### residuals
    idx_dev <- window(idx_dev, end = yr_end)
    catch_res <- window(catch_res, end = yr_end)
    proc_res <- window(proc_res, end = yr_end)
  }
  if (isTRUE(n_iter < 1000)) {
    ### OM stock
    stk_fwd <- FLCore::iter(stk_fwd, seq(n_iter))
    ### OEM stock
    stk_oem <- FLCore::iter(stk_oem, seq(n_iter))
    ### stock-recruitment model 
    sr <- FLCore::iter(sr, seq(n_iter))
    ### index
    idx <- FLCore::iter(idx, seq(n_iter))
    ### residuals
    idx_dev <- FLCore::iter(idx_dev, seq(n_iter))
    catch_res <- FLCore::iter(catch_res, seq(n_iter))
    proc_res <- FLCore::iter(proc_res, seq(n_iter))
  }
  
  ### ---------------------------------------------------------------------- ###
  ### generic arguments ####
  args <- list(fy = yr_end, ### final simulation year
               y0 = dims(stk_fwd)$minyear, ### first data year
               iy = iy, ### first simulation (intermediate) year
               nsqy = 3, ### not used, but has to provided
               nblocks = n_blocks, ### block for parallel processing
               seed = seed ### random number seed before starting MSE
  )
  
  ### ---------------------------------------------------------------------- ###
  ### reference values ####
  refpts_mse <- FLPar(c(refpts_mse), params = row.names(refpts_mse), units = "",
                      iter = n_iter)
  
  ### ---------------------------------------------------------------------- ###
  ### Operating model (OM) ####
  om <- FLom(stock = stk_fwd, ### stock 
             sr = sr, ### stock recruitment and precompiled residuals
             projection = mseCtrl(method = fwd_attr, 
                                  args = list(maxF = 5,
                                              ### process noise on stock.n
                                              proc_res = "fitted",
                                              dupl_trgt = FALSE
                                  ))
  )
  
  ### ---------------------------------------------------------------------- ###
  ### Observation (error) model OEM ####
  
  ### default oem for rfb rule
  oem <- FLoem(method = obs_generic,
               observations = list(
                 stk = stk_oem, 
                 idx = idx), 
               deviances = list(
                 stk = FLQuants(catch.dev = catch_res), 
                 idx = idx_dev),
               args = list(use_catch_residuals = TRUE, 
                           use_idx_residuals = TRUE,
                           use_stk_oem = TRUE,
                           use_wt = TRUE,
                           alks = as.data.frame(ALKs),
                           use_age_idcs = c("Q1SWBeam", "FSP-7e"),
                           biomass_index = "FSP-7e",
                           length_idx = TRUE,
                           Lc = unique(c(refpts_mse["Lc"])),
                           lngth_samples = 2000))
  
  ### 2 over 3 rule: biomass index and PA status relative to OM + error
  if (isTRUE(MP == "2over3")) {
    oem@args$length_idx <- FALSE ### no length index
    oem@args$PA_status <- TRUE 
    oem@args$PA_status_dev <- TRUE
    oem@args$PA_Bmsy <- unique(c(refpts_mse["Bmsy"])) ### real MSY from OM
    oem@args$PA_Fmsy <- unique(c(refpts_mse["Fmsy"]))
    
  ### 2 over 3 based on XSA 
  } else if (isTRUE(MP == "2over3_XSA")) {
    oem@args$length_idx <- FALSE ### no length index
    oem@args$PA_status <- TRUE 
    oem@args$PA_status_dev <- TRUE
  } else if (isTRUE(MP == "ICES_SAM")) {
    oem@observations$idx <- oem@observations$idx[1:2] ### age indices only
    oem@deviances$idx <- oem@deviances$idx[1:2]
    oem@args <- list(cut_idx = TRUE, idx_timing = c(-1, -1),
                     catch_timing = -1, use_catch_residuals = TRUE,
                     use_idx_residuals = TRUE, use_stk_oem = TRUE)
  } else if (isTRUE(MP == "constF")) {
    oem <- FLoem(observations = list(stk = FLQuant(0),
                                     idx = FLQuant()),
                 deviances = list(stk = FLQuant(0),
                                  idx = FLQuant()))
  }
  
  ### ---------------------------------------------------------------------- ###
  ### Management procedure MP ####
  
  if (isTRUE(MP == "rfb")) {
    ### I_trigger = 1.4 * I_loss
    I_trigger = apply(quantSums(index(idx$`FSP-7e`) * idx_dev$`FSP-7e` * 
                                  catch.wt(idx$`FSP-7e`)),
                      6, min, na.rm = TRUE) * 1.4
    ctrl <- mpCtrl(list(
      est = mseCtrl(method = est_comps,
                    args = list(comp_r = TRUE, comp_f = TRUE, comp_b = TRUE,
                                comp_c = TRUE, comp_m = 1,
                                idxB_lag = 1, idxB_range_1 = 2, idxB_range_2 = 3,
                                idxB_range_3 = 1,
                                catch_lag = 0, ### 0 to mimic advice
                                catch_range = 1,
                                Lref = unique(c(refpts_mse["Lref"])), 
                                I_trigger = c(I_trigger),
                                idxL_lag = 1, idxL_range = 1)),
      phcr = mseCtrl(method = phcr_comps,
                     args = list(exp_r = 1, exp_f = 1, exp_b = 1)),
      hcr = mseCtrl(method = hcr_comps,
                    args = list(interval = 2)),
      isys = mseCtrl(method = is_comps,
                     args = list(interval = 2, 
                                 upper_constraint = 1.2, lower_constraint = 0.7, 
                                 cap_below_b = FALSE))
    ))
  } else if (isTRUE(MP == "2over3")) {
    ctrl <- mpCtrl(list(
      est = mseCtrl(method = est_comps,
                    args = list(comp_r = TRUE, comp_f = FALSE, comp_b = FALSE,
                                comp_c = TRUE, comp_m = 1, 
                                idxB_lag = 1, idxB_range_1 = 2, idxB_range_2 = 3,
                                catch_lag = 0, ### 0 to mimic advice
                                catch_range = 1, 
                                pa_buffer = TRUE, 
                                pa_size = 0.8, pa_duration = 3)),
      phcr = mseCtrl(method = phcr_comps),
      hcr = mseCtrl(method = hcr_comps,
                    args = list(interval = 2)),
      isys = mseCtrl(method = is_comps,
                     args = list(interval = 2, 
                                 upper_constraint = 1.2, lower_constraint = 0.8, 
                                 cap_below_b = TRUE))))
  } else if (isTRUE(MP == "2over3_XSA")) {
    FLXSA_control <- FLXSA.control(fse = 1.0, rage = 0, qage = 6, 
                                   shk.n = FALSE, 
                                   shk.ages = 3, shk.yrs = 3, 
                                   min.nse = 0.3, 
                                   tspower = 0, tsrange = 100,
                                   maxit = 100)
    ctrl <- mpCtrl(list(
      est = mseCtrl(method = est_comps,
                    args = list(comp_r = TRUE, comp_f = FALSE, comp_b = FALSE,
                                comp_c = TRUE, comp_m = 1, pa_buffer = TRUE,
                                idxB_lag = 1, 
                                idxB_range_1 = 2, idxB_range_2 = 3,
                                catch_lag = 0, ### 0 to mimic advice
                                catch_range = 1,
                                FLXSA = TRUE,
                                FLXSA_control = FLXSA_control,
                                FLXSA.control = NULL,
                                FLXSA_landings = TRUE, ### landings only?
                                FLXSA_idcs = c("Q1SWBeam", "FSP-7e"),
                                FLXSA_stf = TRUE,
        FLXSA_Btrigger = unique(c(refpts_mse["ICES_Btrigger"])), ### from WGCSE
        FLXSA_Ftrigger = unique(c(refpts_mse["ICES_Fmsy"]))
                                )),
      phcr = mseCtrl(method = phcr_comps),
      hcr = mseCtrl(method = hcr_comps,
                    args = list(interval = 1)), ### annual
      isys = mseCtrl(method = is_comps,
                     args = list(interval = 2, ### annual
                                 upper_constraint = 1.2, lower_constraint = 0.8, 
                                 cap_below_b = TRUE))))
  } else if (isTRUE(MP == "ICES_SAM")) {
    ### some specifications for short term forecast with SAM
    SAM_stf_def <- list(fwd_yrs_rec_start = 1980,
                        fwd_splitLD = TRUE,
                        fwd_yrs_average = -4:0,
                        fwd_yrs_sel = -4:0)
    ctrl <- mpCtrl(list(
      est = mseCtrl(method = SAM_wrapper,
                    args = c(### short term forecast specifications
                      forecast = TRUE, 
                      fwd_trgt = list(c("fsq", "fsq", "fsq")), fwd_yrs = 2, 
                      SAM_stf_def,
                      newtonsteps = 0, rel.tol = 0.001,
                      par_ini = list(SAM_pars_ini),
                      track_ini = TRUE, 
                      conf = list(SAM_conf)
                    )),
      phcr = mseCtrl(method = phcr_WKNSMSE,
                     args = list(
                       Btrigger = unique(c(refpts_mse["EqSim_Btrigger"])), 
                       Ftrgt = unique(c(refpts_mse["EqSim_Fmsy"])), 
                       Blim = unique(c(refpts_mse["EqSim_Blim"])))),
      hcr = mseCtrl(method = hcr_WKNSME, args = list(option = "A")),
      isys = mseCtrl(method = is_WKNSMSE, 
                     args = c(hcrpars = list(
                       Btrigger = unique(c(refpts_mse["EqSim_Btrigger"])), 
                       Ftrgt = unique(c(refpts_mse["EqSim_Fmsy"])), 
                       Blim = unique(c(refpts_mse["EqSim_Blim"]))),
                              fwd_trgt = list(c("fsq", "fsq", "hcr")), 
                              fwd_yrs = 3, SAM_stf_def
                     ))))
  } else if (isTRUE(MP == "constF")) {
    ctrl <- mpCtrl(list(hcr = mseCtrl(method = fixedF.hcr,
                                      args = list(ftrg = 0))))
  }
  
  ### ---------------------------------------------------------------------- ###
  ### additional tracking metrics ####
  tracking <- c("comp_c", "comp_i", "comp_r", "comp_f", "comp_b",
                "multiplier", "exp_r", "exp_f", "exp_b")
  if (isTRUE(MP == "ICES_SAM")) {
    tracking <- c("BB_return", "BB_bank_use", "BB_bank", "BB_borrow")
  }
  
  ### ---------------------------------------------------------------------- ###
  ### create input object ####
  
  input <- list(om = om, oem = oem, ctrl = ctrl,
                args = args, tracking = tracking, refpts = refpts_mse,
                cut_hist = cut_hist)
  
  ### within scenario parallelisation?
  if (isTRUE(n_blocks > 1)) {
    input$args$nblocks <- n_blocks
  }
  
  return(input)
  
}
