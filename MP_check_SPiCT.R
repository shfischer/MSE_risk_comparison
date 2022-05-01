### ------------------------------------------------------------------------ ###
### explore SPiCT for plaice, cod, and herring ####
### ------------------------------------------------------------------------ ###

library(spict)
library(FLCore)

### ------------------------------------------------------------------------ ###
### cod ####
### ------------------------------------------------------------------------ ###

stk <- readRDS("input/cod.27.47d20/preparation/stk.rds")
idx <- readRDS("input/cod.27.47d20/preparation/idx.rds")
### survey weights
idx_Q1_wts <- readRDS("input/cod.27.47d20/preparation/idx_weights_Q1.rds")
idx_Q3_wts <- readRDS("input/cod.27.47d20/preparation/idx_weights_Q3.rds")
### add to survey
catch.wt(idx$IBTS_Q1_gam) <- idx_Q1_wts
catch.wt(idx$IBTS_Q3_gam) <- idx_Q3_wts

input <- list(obsC = c(catch(stk)),
              timeC = as.numeric(dimnames(stk)$year),
              obsI = c(quantSums(idx$IBTS_Q3_gam@index * 
                                   idx$IBTS_Q3_gam@catch.wt)),
              timeI = as.numeric(dimnames(idx$IBTS_Q3_gam)$year) +
                mean(idx$IBTS_Q3_gam@range[c("startf", "endf")]))
input <- list(obsC = c(catch(stk)),
              timeC = as.numeric(dimnames(stk)$year),
              obsI = list(Q1 = c(quantSums(idx$IBTS_Q1_gam@index * 
                                             idx$IBTS_Q1_gam@catch.wt)),
                          Q3 = c(quantSums(idx$IBTS_Q3_gam@index * 
                                             idx$IBTS_Q3_gam@catch.wt))),
              timeI = list(Q1 = as.numeric(dimnames(idx$IBTS_Q1_gam)$year) +
                             mean(idx$IBTS_Q1_gam@range[c("startf", "endf")]),
                           Q3 = as.numeric(dimnames(idx$IBTS_Q3_gam)$year) +
                             mean(idx$IBTS_Q3_gam@range[c("startf", "endf")])))

input <- check.inp(input)

fit <- fit.spict(input)

plot(fit)
fit_res <- calc.osa.resid(fit)
plot(fit_res)

pdf(file = "output/plots/SPiCT/SPiCT_cod.pdf", 
    width = 16, height = 12)
plot(fit_res)
dev.off()

### ------------------------------------------------------------------------ ###
### plaice ####
### ------------------------------------------------------------------------ ###

stk <- readRDS("input/ple.27.7e/preparation/model_input_stk_d.RDS")
idx <- readRDS("input/ple.27.7e/preparation/model_input_idx.RDS")
### add weights to indices
catch.wt(idx$`FSP-7e`) <- catch.wt(stk)[ac(2:8), ac(2003:2020)]
catch.wt(idx$Q1SWBeam) <- catch.wt(stk)[ac(2:9), ac(2006:2020)]

input <- list(obsC = c(catch(stk)),
              timeC = as.numeric(dimnames(stk)$year),
              obsI = list(Q1 = c(quantSums(idx$`FSP-7e`@index * 
                                             idx$`FSP-7e`@catch.wt)),
                          Q3 = c(quantSums(idx$Q1SWBeam@index * 
                                             idx$Q1SWBeam@catch.wt))),
              timeI = list(Q1 = as.numeric(dimnames(idx$`FSP-7e`)$year) +
                             mean(idx$`FSP-7e`@range[c("startf", "endf")]),
                           Q3 = as.numeric(dimnames(idx$Q1SWBeam)$year) +
                             mean(idx$Q1SWBeam@range[c("startf", "endf")])))

input <- check.inp(input)

fit <- fit.spict(input)
plot(fit)
fit_res <- calc.osa.resid(fit)
plot(fit_res)

pdf(file = "output/plots/SPiCT/SPiCT_ple.pdf", 
    width = 16, height = 12)
plot(fit_res)
dev.off()

### ------------------------------------------------------------------------ ###
### herring ####
### ------------------------------------------------------------------------ ###


stk <- readRDS("input/her.27.3a47d/preparation/stk.rds")
idx <- readRDS("input/her.27.3a47d/preparation/idx_wts.rds")

input <- list(obsC = c(catch(stk)),
              timeC = as.numeric(dimnames(stk)$year),
              obsI = c(quantSums(idx$HERAS@index * 
                                   idx$HERAS@catch.wt)),
              timeI = as.numeric(dimnames(idx$HERAS)$year) +
                mean(idx$HERAS@range[c("startf", "endf")]))

input <- check.inp(input)

fit <- fit.spict(input)
plot(fit)
fit_res <- calc.osa.resid(fit)
plot(fit_res)

pdf(file = "output/plots/SPiCT/SPiCT_her.pdf", 
    width = 16, height = 12)
plot(fit_res)
dev.off()

plotspict.diagnostic(fit_res)

### check guidelines
fit$opt$convergence # ok
all(is.finite(fit$sd)) # ok
plotspict.diagnostic(fit_res) # failed normality
fit <- retro(fit)
plotspict.retro(fit) # retro pattern ok
calc.bmsyk(fit) # production curve ok
calc.om(fit) # uncertainty magnitude ok
fit
fit <- check.ini(fit)
fit$check.ini$resmat # issues:
# does not converge in 5/10 runs with different initial parameters
# K doubled in one run

system.time({fit2 <- fit.spict(input)})
