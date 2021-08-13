### ------------------------------------------------------------------------ ###
### run EqSim to estimate data-rich ICES reference points ####
### ------------------------------------------------------------------------ ###

library(stockassessment)
library(FLfse)
library(ggplotFL)
library(msy)
library(ggplot2)
library(tidyr)
library(dplyr)

### ------------------------------------------------------------------------ ###
### load SAM model fit ####
### ------------------------------------------------------------------------ ###
fit <- readRDS("input/fit.rds")
plot(fit)
srplot(fit)

stk <- SAM2FLStock(fit)


### ------------------------------------------------------------------------ ###
### basis ####
### ------------------------------------------------------------------------ ###
### follow ICES guidelines
### ICES. 2017. ICES fisheries management reference points for category 1 and 
###   2 stocks. In ICES Advice Technical Guidelines. ICES Advice 2017, pp. 1–19.


### ------------------------------------------------------------------------ ###
### Biomass limit reference points - Blim and Bpa ####
### ------------------------------------------------------------------------ ###
### Stock type: Type 5: Stocks showing no evidence of impaired recruitment or 
###   with no clear relation between stock and recruitment (no apparent S–R 
###   signal).
### -> Blim = Bloss

Blim <- min(ssb(stk))
# dimnames(stk)$year[which.min(ssb(stk))]
### 2109.681 tonnes in 2007
### rounded: 2100 

### set Bpa to Blim * 1.4 (default)
Bpa <- Blim * 1.4
### 2953.553 tonnes
### rounded: 2950

### ------------------------------------------------------------------------ ###
### F limit reference points - Flim and Fpa ####
### ------------------------------------------------------------------------ ###

### use preferred option (option a) from ICES guidelines:
### "The preferred method is simulating stock with segmented regression S–R 
### relationship with the point of inflection at Blim, thus determining the 
### F = Flim that, at equilibrium, gives a 50% probability of SSB > Blim. 
### Note: this simulation should be conducted based on a fixed F (i.e. without 
### inclusion of a Btrigger) and without inclusion of assessment/advice errors. 
### (In the EqSim software, this means Btrigger, Fcv, and Fphishould all be set 
### to zero.)"

### create stock recruitment model
set.seed(123)
FIT <- eqsr_fit(stk = stk, nsamp = 1000, 
                models = c("Ricker", "Segreg", "Bevholt"), 
                id.sr = "recruitment model exploration")
eqsr_plot(FIT)
png("input/eqsim/eqsim_rec_models_exploration.png", width = 20, height = 15, 
    res = 300, units = "cm", type = "cairo")
eqsr_plot(FIT)
dev.off()
### segmented regression only
set.seed(123)
FIT_sr <- eqsr_fit(stk = stk, nsamp = 1000, 
                models = c("Segreg"), 
                id.sr = "segmented regression")
eqsr_plot(FIT_sr)
png("input/eqsim/eqsim_rec_models_segreg.png", width = 20, height = 15, 
    res = 300, units = "cm", type = "cairo")
eqsr_plot(FIT_sr)
dev.off()
### segmented regression with breakpoint fixed to Blim
segreg_fixed  <- function(ab, ssb) {
  log(ifelse(ssb >= Blim, ab$a * Blim, ab$a * ssb))
}
set.seed(123)
FIT_fixed <- eqsr_fit(stk = stk, nsamp = 1000, 
                 models = c("segreg_fixed"), 
                 id.sr = "Blim fixed segmented regression")
eqsr_plot(FIT_fixed)
png("input/eqsim/eqsim_rec_models_fixed_Blim.png", width = 20, height = 15, 
    res = 300, units = "cm", type = "cairo")
eqsr_plot(FIT_fixed)
dev.off()

### project forward with various constant Fs without ICES advice rule
### no uncertainty
### use fixed segmented regression (Blim as breakpoint)
set.seed(123)
SIM_det <- eqsim_run(fit = FIT_sr, 
                     bio.years = c(2016, 2020), 
                     sel.years = c(2016, 2020),
                     Fscan = seq(0, 1.5, length.out = 151), 
                     Fcv = 0, Fphi = 0, Btrigger = 0,
                     Blim = Blim, Bpa = Bpa)
### find F at Blim
df <- SIM_det$rbp %>%
  filter(variable %in% c("Spawning stock biomass")) %>%
  select(Ftarget, p50) %>%
  transmute(F = Ftarget, SSB = p50)
b.lm <- loess(F ~ SSB, data = df, span = 0.2)
Flim <- predict(b.lm, Blim)
### 0.5487305

p <- df %>% ggplot(aes(x = F, y = SSB)) +
  geom_point(size = 0.4) + geom_line(size = 0.2) +
  theme_bw() +
  geom_hline(yintercept = Blim) +
  geom_vline(xintercept = Flim)
p
ggsave(file = "input/eqsim/eqsim_ssb_vs_F_Flim.png", 
       width = 13.4, height = 8, units = "cm", dpi = 300, plot = p)

### calculate Fpa as Flim * 1.4 (default)
Fpa <- Flim / 1.4
### 0.3919504

### ------------------------------------------------------------------------ ###
### Fmsy ####
### ------------------------------------------------------------------------ ###

### set MSY Btrigger to Bpa
MSYBtrigger <- Bpa

### step 1: without Btrigger
### including default uncertainty
### recruitment: use mix of models (Ricker, Segreg, BevHolt)
set.seed(123)
SIM_Fmsy <- eqsim_run(fit = FIT, 
                      bio.years = c(2016, 2020), 
                      sel.years = c(2016, 2020),
                      Fscan = seq(0, 1.5, length.out = 151), 
                      Fcv = 0.212, Fphi = 0.423, Btrigger = 0,
                      Blim = Blim, Bpa = Bpa)
### some plots of results
eqsim_plot(SIM_Fmsy)
png("input/eqsim/eqsim_Fmsy.png", width = 20, height = 15, res = 300, 
    units = "cm", type = "cairo")
eqsim_plot(SIM_Fmsy)
dev.off()
eqsim_plot_range(SIM_Fmsy)
png("input/eqsim/eqsim_Fmsy_range.png", width = 20, height = 15, res = 300, 
    units = "cm", type = "cairo")
eqsim_plot_range(SIM_Fmsy)
dev.off()

### extract Fmsy
#Fmsy <- SIM_Fmsy$Refs["catF","medianMSY"]
Fmsy <- SIM_Fmsy$refs_interval$FmsyMedianC
### 0.241206
### compare with Fpa
if (Fmsy > Fpa) Fmsy <- Fpa
### Fmsy < Fpa
Fmsy_lower <- SIM_Fmsy$refs_interval["FmsylowerMedianC"]
### 0.1658291
Fmsy_upper <- SIM_Fmsy$refs_interval["FmsyupperMedianC"]
### 0.3316583

### check risk with MSYBtrigger=Bpa
set.seed(123)
SIM_Fmsy_risk <- eqsim_run(fit = FIT, 
                           bio.years = c(2016, 2020), 
                           sel.years = c(2016, 2020),
                           Fscan = seq(0, 1.5, length.out = 151), 
                           Fcv = 0.212, Fphi = 0.423, Btrigger = Bpa,
                           Blim = Blim, Bpa = Bpa)
eqsim_plot(SIM_Fmsy_risk)
png("model/eqsim_Fmsy_with_Bpa.png", width = 20, height = 15, res = 300, 
    units = "cm", type = "cairo")
eqsim_plot(SIM_Fmsy_risk)
dev.off()
F0.5 <- SIM_Fmsy_risk$refs_interval$F5percRiskBlim
### 0.4245526
### Fsmy_upper < F0.5 (with Btrigger), so no need to cap values

### ------------------------------------------------------------------------ ###
### MSY Btrigger ####
### ------------------------------------------------------------------------ ###

df <- SIM_det$rbp %>%
  filter(variable %in% c("Spawning stock biomass")) %>%
  select(Ftarget, p05) %>%
  transmute(F = Ftarget, SSB = p05)
b.lm <- loess(SSB ~ F, data = df, span = 0.2)
MSYBtrigger <- predict(b.lm, Fmsy)
### 3894.557
### check risk
set.seed(123)
SIM_Fmsy_risk_Btrigger <- eqsim_run(fit = FIT, 
                                    bio.years = c(2016, 2020), 
                                    sel.years = c(2016, 2020),
                                    Fscan = seq(0, 1.5, length.out = 151), 
                                    Fcv = 0.212, Fphi = 0.423, 
                                    Btrigger = MSYBtrigger,
                                    Blim = Blim, Bpa = Bpa)
eqsim_plot(SIM_Fmsy_risk_Btrigger)
png("input/eqsim/eqsim_Fmsy_with_Btrigger.png", width = 20, height = 15, res = 300, 
    units = "cm", type = "cairo")
eqsim_plot(SIM_Fmsy_risk_Btrigger)
dev.off()
SIM_Fmsy_risk_Btrigger$refs_interval
### would be precautionary
### but keep MSYBtrigger at Bpa because F in recent years > Fmsy
MSYBtrigger <- Bpa
