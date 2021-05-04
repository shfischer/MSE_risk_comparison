### ------------------------------------------------------------------------ ###
### analyse MSE results ####
### ------------------------------------------------------------------------ ###

library(mse)
library(tidyr)
library(dplyr)

source("funs.R")

### ------------------------------------------------------------------------ ###
### load data ####
### ------------------------------------------------------------------------ ###



res_FLXSA <- readRDS("output/res_500_FLXSA.rds")
res_2over3 <- readRDS("output/res_500_2over3.rds")
res_rfb <- readRDS("output/res_500_rfb.rds")

res_list <- list("2 over 3 XSA" = res_FLXSA, 
                 "2 over 3" = res_2over3, 
                 "rfb" = res_rfb)

### ------------------------------------------------------------------------ ###
### stock plots ####
### ------------------------------------------------------------------------ ###


### extract SSB, catch, fbar and recruitment
### "correct": once SSB < 1, stays at 0
res_corrected <- lapply(res_list, function(x) {#browser()
  collapse_correction(x@stock, threshold = 1)
})
res <- lapply(names(res_list), function(x) {
  ### get percentiles
  tmp <- lapply(res_corrected[[x]], 
                quantile, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
  tmp <- as.data.frame(as(tmp, "FLQuants"))
  tmp$rule <- x
  return(tmp)
})
res <- do.call(rbind, res)
head(res)

res %>%
  mutate(data = ifelse(qname == "fbar" & data > 1, 1, data)) %>%
  pivot_wider(values_from = data, names_from = iter) %>%
  mutate(rule = factor(rule, 
                       levels = c("2 over 3 XSA", "2 over 3", "rfb"))) %>%
  mutate(qname = factor(qname,
                        levels = c("rec", "ssb", "catch", "fbar"),
                        labels = c("Recruitment\n[thousands]\n(age 2)",
                                   "SSB [tonnes]",
                                   "Catch [tonnes]",
                                   "F (ages 3-6)"
                                   ))) %>%
  ggplot(aes(x = year, y = `50%`)) +
  geom_ribbon(aes(ymin = `5%`, ymax = `95%`), 
              fill = "grey", alpha = 0.5) +
  geom_ribbon(aes(ymin = `25%`, ymax = `75%`), 
              fill = "grey", alpha = 0.99) +
  geom_line(size = 0.3) +
  geom_vline(xintercept = 2020.5, size = 0.3) + 
  facet_grid(qname ~ rule, scales = "free", switch = "y") +
  labs(y = "") +
  theme_bw() +
  theme(strip.placement = "outside",
        strip.background.y = element_blank(),
        axis.title.y = element_blank())
ggsave(filename = "output/comparison.png", 
       width = 17, height = 12, units = "cm", dpi = 600, type = "cairo")



### ------------------------------------------------------------------------ ###
### risk ####
### ------------------------------------------------------------------------ ###
library(stockassessment)
fit <- readRDS("output/fit.rds")
Blim <- min(ssbtable(fit)[,"Estimate"])
Blim

risks <- lapply(names(res_corrected), function(x) {#browser()
  ssb <- window(res_corrected[[x]]$ssb, start = 2021)
  annual <- c(apply(ssb < Blim, 2, mean))
  total <- sapply(seq(dims(ssb)$year), function(y) {
    mean(c(ssb[, seq(y)]) < Blim, na.rm = TRUE)
  })
  data.frame(year = as.numeric(dimnames(ssb)$year),
             annual = annual, total = total,
             rule = x)
})
risks <- do.call(rbind, risks)

risks %>% 
  pivot_longer(c(annual, total)) %>%
  mutate(rule = factor(rule, 
                          levels = c("2 over 3 XSA", "2 over 3", "rfb")),
         name = factor(name, levels = c("annual", "total"))) %>%
  ggplot(aes(x = year, y = value, colour = rule)) +
  geom_line() +
  geom_hline(yintercept = 0.05, colour = "black", linetype = "dashed", 
             size = 0.3) +
  scale_color_discrete("") +
  facet_wrap(~ name, ncol = 1) +
  labs(y = expression(B[lim]~risk)) + 
  scale_x_continuous(breaks = c(2020, 2040, 2060, 2080, 2100, 2120)) +
  theme_bw()
ggsave(filename = "output/comparison_risk.png", 
       width = 17, height = 12, units = "cm", dpi = 600, type = "cairo")

