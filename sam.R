
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



