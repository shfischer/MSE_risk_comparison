
library(FLCore)
library(ggplotFL)
library(FLfse)
library(stockassessment)

### ------------------------------------------------------------------------ ###
### load input data ####
### ------------------------------------------------------------------------ ###

stk <- readRDS("input/model_input_stk_d.RDS")
idx <- readRDS("input/model_input_idx.RDS")

stk_XSA <- readRDS("input/stock_d.rds")

### ------------------------------------------------------------------------ ###
### discards trials ####
### ------------------------------------------------------------------------ ###

### half of discards
stk_d50 <- stk
discards.n(stk_d50) <- discards.n(stk_d50)/2
discards(stk_d50) <- computeDiscards(stk_d50)
catch(stk_d50) <- computeCatch(stk_d50, slot = "all")

plot(FLStocks(half = stk_d50, full = stk))
fit_d50 <- FLR_SAM(stk = stk_d50, idx = idx)
retro_d50 <- retro(fit_d50, year = 10)
plot(retro_d50)
retro_d50 <- retro(fit_d50, year = 5)
mohn(retro_d50)

### no discards before 2012
stk_d2012 <- stk
discards.n(stk_d2012)[, ac(1980:2011)] <- 0
discards.wt(stk_d2012)[, ac(1980:2011)] <- 0
discards(stk_d2012) <- computeDiscards(stk_d2012)
catch(stk_d2012) <- computeCatch(stk_d2012, slot = "all")
plot(FLStocks("2012" = stk_d2012, full = stk))
fit_d2012 <- FLR_SAM(stk = stk_d2012, idx = idx)
retro_d2012 <- retro(fit_d2012, year = 10)
plot(retro_d2012)
retro_d2012 <- retro(fit_d2012, year = 5)
mohn(retro_d2012)
# R(age 2)         SSB   Fbar(3-6) 
# -0.18063963  0.05647796 -0.09316801 

### no discards - retro worse & large patterns when going back >5 years
stk_d0 <- stk
discards.n(stk_d0) <- 0
discards.wt(stk_d0) <- 0
discards(stk_d0) <- computeDiscards(stk_d0)
catch(stk_d0) <- computeCatch(stk_d0, slot = "all")
plot(FLStocks("0" = stk_d0, full = stk))
fit_d0 <- FLR_SAM(stk = stk_d0, idx = idx)
retro_d0 <- retro(fit_d0, year = 10)
plot(retro_d0)
retro_d0 <- retro(fit_d0, year = 5)
mohn(retro_d0)
# R(age 2)        SSB  Fbar(3-6) 
# -0.0259818  0.1446794 -0.1680583

### catch at age 2 = 0
stk_d_a2_0 <- stk
catch.n(stk_d_a2_0)[1, ] <- NA
fit_d_a2_0 <- FLR_SAM(stk = stk_d_a2_0, idx = idx)
retro_d_a2_0 <- retro(fit_d_a2_0, year = 10)
plot(retro_d_a2_0)
retro_d_a2_0 <- retro(fit_d_a2_0, year = 5)
mohn(retro_d_a2_0)
# R(age 2)         SSB   Fbar(3-6) 
# -0.12957232  0.09250615 -0.12563298

### discard survival rates
library(doParallel)
cl <- makeCluster(5)
registerDoParallel(cl)

surv <- seq(0, 1, 0.05)
names(surv) <- surv
stks <- lapply(surv, function(x) {
  stk_i <- stk
  discards.n(stk_i) <- discards.n(stk_i) * (1 - x)
  discards(stk_i) <- computeDiscards(stk_i)
  catch(stk_i) <- computeCatch(stk_i, slot = "all")
  return(stk_i)
})
fits <- foreach(x = stks, 
                .packages = c("stockassessment", "FLfse")) %dopar% {
  FLR_SAM(x, idx)
}
ll <- sapply(fits, function(x) c(x$opt$objective))
# aics <- sapply(fits, AIC) # same pattern
plot(x = surv, y = ll, type = "l")
which.min(ll)
# min at surv=0.29 ~ 0.3
rets <- foreach(x = fits,
                .packages = c("stockassessment")) %dopar% {
  retro(x, year = 5)
}
mohns <- lapply(rets, mohn)
mohns

### lower plusgroup
### need to reduce PG to 7 to see improvement...
stk_pg <- setPlusGroup(stk, 7)
idx_pg <- idx
idx_pg$`FSP-7e` <- idx_pg$`FSP-7e`[ac(2:6), ]
idx_pg$Q1SWBeam <- idx_pg$Q1SWBeam[ac(2:6), ]
fit_pg <- FLR_SAM(stk = stk_pg, idx = idx_pg)
fit_pg
plot(fit_pg)
ret_pg <- retro(fit_pg, year = 5)
plot(ret_pg)
mohn(ret_pg)

### change idx ages
### removing young/old ages from surveys makes fit worse/increases retro pattern
idx_ages <- idx
idx_ages$`FSP-7e` <- idx_ages$`FSP-7e`[ac(3:8), ]
idx_ages$Q1SWBeam <- idx_ages$Q1SWBeam[ac(3:7), ]
fit_ages <- FLR_SAM(stk = stk, idx = idx_ages)
fit_ages
plot(fit_ages)
ret_ages <- retro(fit_ages, year = 5)
plot(ret_ages)
mohn(ret_ages)

### change SAM configuration
### keyVarObs - coupling of observation variance parameters - bit better
conf_i <- fit$conf
conf_i$keyVarObs[2, 1:7] <- c(1, 1, rep(2, 3), 3, 3)
conf_i$keyVarObs[3, 1:8] <- c(4, 4,  rep(5, 4), 6, 6)
conf_i$keyVarObs
fit_conf <- FLR_SAM(stk = stk, idx = idx, conf = conf_i, conf_full = TRUE)
fit_conf
plot(fit_conf)
ret_conf <- retro(fit_conf, year = 5)
plot(ret_conf)
mohn(ret_conf)

### keyLogFsta - coupling of F states - tiny improvement
conf_i <- fit$conf
conf_i$keyLogFsta[1, 1:9] <- c(0:5, 6, 6, 6)
conf_i$keyLogFsta
fit_conf <- FLR_SAM(stk = stk, idx = idx, conf = conf_i, conf_full = TRUE)
fit_conf
plot(fit_conf)
ret_conf <- retro(fit_conf, year = 5)
plot(ret_conf)
mohn(ret_conf)

### keyLogFpar - coupling of survey catchability - minor
conf_i <- fit$conf
conf_i$keyLogFpar[2, 1:7] <- c(0, 0, 1:3, 4, 4)
conf_i$keyLogFpar[3, 1:8] <- c(5, 5, 6:9, 10, 10)
conf_i$keyLogFpar
fit_conf <- FLR_SAM(stk = stk, idx = idx, conf = conf_i, conf_full = TRUE)
fit_conf
plot(fit_conf)
ret_conf <- retro(fit_conf, year = 5)
plot(ret_conf)
mohn(ret_conf)

### keyVarF - coupling of F random walk - minor
conf_i <- fit$conf
conf_i$keyVarF[1, ] <- c(0, 0, 1, 1, 1, 1, 1, 2, 2)
conf_i$keyVarF
fit_conf <- FLR_SAM(stk = stk, idx = idx, conf = conf_i, conf_full = TRUE)
fit_conf
plot(fit_conf)
ret_conf <- retro(fit_conf, year = 5)
plot(ret_conf)
mohn(ret_conf)

### remove surveys
### remove Q1SWBeam - minor impact
### remove FSP - large retro
lo <- leaveout(fit)
plot(lo)
fit_i <- FLR_SAM(stk = stk, idx = idx$`FSP-7e`)
fit_i
plot(fit_i)
ret_i <- retro(fit_i, year = 5)
plot(ret_i)
mohn(ret_i)
fit_i <- FLR_SAM(stk = stk, idx = idx$Q1SWBeam)
fit_i
plot(fit_i)
ret_i <- retro(fit_i, year = 5)
plot(ret_i)
mohn(ret_i)

### ------------------------------------------------------------------------ ###
### fit SAM with default configuration ####
### ------------------------------------------------------------------------ ###

fit <- FLR_SAM(stk = stk, idx = idx)
plot(fit)

plot(FLStocks(XSA = stk_XSA, SAM = SAM2FLStock(fit)))

### ------------------------------------------------------------------------ ###
### adapt SAM configuration ####
### ------------------------------------------------------------------------ ###

### get default configuration
conf <- fit$conf

### couple survey catchabilities
### adopt XSA settings: couple after age 6
names(idx)
dims(idx$`FSP-7e`)$age
dimnames(idx$`FSP-7e`)$age
conf$keyLogFpar[2, 1:7] <- c(0:4, 5, 5)
dims(idx$Q1SWBeam)$age
dimnames(idx$Q1SWBeam)$age
conf$keyLogFpar[3, 1:8] <- c(6:10, 11, 11, 11)

### ------------------------------------------------------------------------ ###
### fit SAM ####
### ------------------------------------------------------------------------ ###

fit <- FLR_SAM(stk = stk, idx = idx, conf = conf)
plot(fit)

plot(FLStocks(XSA = stk_XSA, SAM = SAM2FLStock(fit)))

### ------------------------------------------------------------------------ ###
### diagnostics ####
### ------------------------------------------------------------------------ ###

fit
summary(fit)

### ------------------------------------------------------------------------ ###
### model fit ####

plot(fit)
catchplot(fit)

png(filename = "output/SAM/fit.png", width = 25, height = 20, units = "cm", 
    res = 600, type = "cairo")
#windows()
par(mfrow = c(2, 2))
catchplot(fit)
recplot(fit)
fbarplot(fit)
ssbplot(fit)
dev.off()

### ------------------------------------------------------------------------ ###
### residuals ####
res <- residuals(fit)
plot(res)
png(filename = "output/SAM/residuals.png", width = 13, height = 15, units = "cm", 
    res = 600, type = "cairo")
plot(res)
dev.off()


### ------------------------------------------------------------------------ ###
### process residuals ####
resp <- procres(fit)
plot(resp)
png(filename = "output/SAM/process_residuals.png", width = 20, height = 15, 
    units = "cm", res = 600, type = "cairo")
plot(resp)
dev.off()

### ------------------------------------------------------------------------ ###
### between-age correlation by fleet ####
corplot(fit)
corplot(res)

### ------------------------------------------------------------------------ ###
### retro ####
retro <- retro(fit, year = 5)
plot(retro)
mohn(retro)
#    R(age 2)         SSB   Fbar(3-6) 
# -0.14044561  0.09763852 -0.13301581 
png(filename = "output/SAM/retro.png", width = 13, height = 20, units = "cm", 
    res = 600, type = "cairo")
plot(retro)
dev.off()

### ------------------------------------------------------------------------ ###
### survey leave out ####
lo <- leaveout(fit)
plot(lo)
png(filename = "output/SAM/leaveout.png", width = 13, height = 20, units = "cm", 
    res = 600, type = "cairo")
plot(lo)
dev.off()

### ------------------------------------------------------------------------ ###
### simulate data from fitted model and re-estimate from each run
sim <- simstudy(fit, nsim = 100)
plot(sim)

### ------------------------------------------------------------------------ ###
### jitter: start from random inital values
set.seed(1)
jit <- jit(fit, nojit = 100)
jit
plot(jit)

### ------------------------------------------------------------------------ ###
### YPR
ypr(fit)

### ------------------------------------------------------------------------ ###
### save fit ####
### ------------------------------------------------------------------------ ###

saveRDS(fit, "output/fit.rds")

### ------------------------------------------------------------------------ ###
### stock-recruit curve ####
### ------------------------------------------------------------------------ ###
library(ggrepel)

srplot(fit)

### get SR pairs with recruitment lag (rec age 2)
sr <- as.FLSR(SAM2FLStock(fit), model = "segreg")
sr <- fmle(sr)
plot(sr)
rec <- data.frame(year = as.numeric(dimnames(rec(sr))$year),
                  SSB = c(ssb(sr)),
                  rec = c(rec(sr)))
rec %>%
  ggplot(aes(x = SSB, y = rec)) +
  geom_path(colour = "grey", size = 0.2) + 
  geom_text_repel(aes(label = year), size = 2, segment.size = 0.3) +
  geom_point(size = 0.3) +
  ylim(0, NA) + xlim(0, NA) +
  labs(x = "SSB [tonnes]", y = "Recruitment (age 2) [thousands]") +
  theme_bw()
ggsave(filename = "output/SAM/sr.png", 
       width = 12, height = 8, units = "cm", dpi = 600, type = "cairo")

### ------------------------------------------------------------------------ ###
### compare to XSA results ####
### ------------------------------------------------------------------------ ###

### load XSA results
stk_d <- readRDS("input/stock_d.rds") ### total catch
stk_l <- readRDS("input/stock.rds") ### landings only

### list with stocks
stks <- FLStocks("XSA (landings)" = stk_l,
                 "XSA (catch)" = stk_d,
                 "SAM (catch)" = SAM2FLStock(fit))

### extract data for plotting
stks_df <- lapply(names(stks), function(x) {#browser()
  stk_i <- stks[[x]]
  df <- FLQuants(ssb = ssb(stk_i),
                 rec = rec(stk_i),
                 fbar = fbar(stk_i),
                 catch = catch(stk_i))
  df <- as.data.frame(df)
  df$ass <- x
  return(df)
})
stks_df <- do.call(rbind, stks_df)

stks_df %>%
  mutate(qname = factor(qname,
                           levels = c("catch", "rec", "fbar", "ssb"),
                           labels = c("Catch [tonnes]",
                                      "Recruitment (age 2)\n[thousands]",
                                      "F (ages 3-6)",
                                      "SSB [tonnes]")),
         ass = factor(ass, 
                      levels = c("XSA (landings)", "XSA (catch)", 
                                 "SAM (catch)"),
                      labels = c("XSA\n(landings)", "XSA\n(catch)", 
                                 "SAM\n(catch)"))) %>%
  ggplot(aes(x = year, y = data, colour = ass)) +
  geom_line() +
  scale_colour_discrete("") +
  facet_wrap(~ qname, scales = "free_y", strip.position = "left") +
  ylim(0, NA) +
  theme_bw() +
  theme(strip.placement = "outside",
        strip.background.y = element_blank(),
        axis.title.y = element_blank())
ggsave(filename = "output/SAM/SAM_XSA_comparison.png", 
       width = 15, height = 10, units = "cm", dpi = 600, type = "cairo")


### ------------------------------------------------------------------------ ###
### final fit ####
### ------------------------------------------------------------------------ ###

stk <- readRDS("input/model_input_stk_d.RDS")
idx <- readRDS("input/model_input_idx.RDS")

### use configuration similar to accepted XSA assessment
conf <- list(keyLogFpar = 
               matrix(data = c(rep(-1, 9),
                               0:5, 5, -1, -1,
                               6:11, 11, 11, -1),
                      ncol = 9, nrow = 3, byrow = TRUE))

fit <- FLR_SAM(stk, idx, conf = conf)

saveRDS(fit, file = "input/fit.rds")


