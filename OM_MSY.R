### ------------------------------------------------------------------------ ###
### estimate MSY ####
### ------------------------------------------------------------------------ ###
library(mse)
source("funs.R")
library(doParallel)
library(tidyr)
library(dplyr)
library(ggplot2)


### ------------------------------------------------------------------------ ###
### set up parallelisation ####
### ------------------------------------------------------------------------ ###
cl <- makeCluster(10)
registerDoParallel(cl)
cl_length <- length(cl)
. <- foreach(i = seq(cl_length)) %dopar% {
  library(mse, quietly = TRUE)
  source("funs.R")
}

### ------------------------------------------------------------------------ ###
### load input data ####
### ------------------------------------------------------------------------ ###

input_constF <- input_mp(stock_id = "ple.27.7e", OM = "baseline", 
                         yr_start = 2021, MP = "constF")

### include parallelisation
input_constF$args$nblocks <- 10

### ------------------------------------------------------------------------ ###
### function for running MP and returning catch ####
### ------------------------------------------------------------------------ ###
### return median catch of last 10 years from 100-year projection

res_trace <- list()
mp_catch <- function(input, Ftrgt, minimise = FALSE) {#browser()
  if (isTRUE(Ftrgt %in% sapply(res_trace, function(x) x$Ftrgt))) {
    res_i <- res_trace[[which(Ftrgt == sapply(res_trace, function(x) x$Ftrgt))[[1]]]]
    catch_i <- res_i$catch
    ssb_i <- res_i$ssb
    tsb_i <- res_i$ssb
    rec_i <- res_i$rec
  } else {
    input$ctrl$hcr@args$ftrg <- Ftrgt
    res_i <- do.call(mp, input)
    catch_i <- median(tail(catch(res_i@stock), 10), na.rm = TRUE)
    ssb_i <- median(tail(ssb(res_i@stock), 10), na.rm = TRUE)
    tsb_i <- median(tail(tsb(res_i@stock), 10), na.rm = TRUE)
    rec_i <- median(tail(rec(res_i@stock), 10), na.rm = TRUE)
  }
  cat(paste0("Ftrgt=", Ftrgt, "; C=", catch_i, "; SSB=", ssb_i, "; R=", rec_i,
             "\n"))
  res_trace <<- append(res_trace, 
                       list(list(Ftrgt = Ftrgt, catch = catch_i,
                                 ssb = ssb_i, tsb = tsb_i, rec = rec_i)))
  if (isTRUE(minimise)) catch_i <- -catch_i
  return(catch_i)
}
#mp_catch(input = input_constF, Ftrgt = 0.2)

### ------------------------------------------------------------------------ ###
### optimise ####
### ------------------------------------------------------------------------ ###

### first, check some values
vals <- seq(0, 1, 0.1)
names(vals) <- vals
res_trace <- list()
res_0.1 <- lapply(vals, mp_catch, input = input_constF)
res_0.1_trace <- res_trace
saveRDS(res_0.1_trace, "input/ple.27.7e/baseline/1000_100/MSY_0.1.rds")

# set.seed(1)
# res_MSY <- optim(par = 0.2, fn = mp_catch, input = input_constF,
#                  method = "L-BFGS-B", ### to allow for constraints
#                  lower = 0, upper = 1,
#                  control = list(fnscale = -1,
#                                 factr = 1,
#                                 maxit = 5,
#                                 trace = 6))
# bckp <- res_trace

### try with optimise - 1D golden-section search
# set.seed(1)
# optimise(f = dnorm, 
#                 interval = c(-6, 5), 
#                 lower = -6, upper = 5,
#                 maximum = TRUE,
#                 tol = 0.001)
set.seed(1)
res_trace <- list()
res_optimise_MSY <- optimise(f = mp_catch, input = input_constF,
                             interval = c(0, 0.3), 
                             lower = 0, upper = 0.3,
                             maximum = TRUE,
                             tol = 0.001)
res_trace_optimise <- res_trace
saveRDS(res_trace_optimise, "input/ple.27.7e/baseline/1000_100/MSY_search.rds")

# set.seed(1)
# res_MSY <- nlm(f = mp_catch, p = 0.3, input = input_constF,
#                print.level = 2,
#                ndigit = 2, 
#                iterlim = 10,
#                  lower = 0, upper = 1,
#                  control = list(fnscale = -1,
#                                 abstol = 1,
#                                 maxit = 5))
# 
# blubb <- mp_catch(input = input_constF, Ftrgt = 0.5)

### ------------------------------------------------------------------------ ###
### plot MSY estimation ####
### ------------------------------------------------------------------------ ###

res_trace <- c(res_0.1_trace, res_trace_optimise)
res_trace <- bind_rows(res_trace)

res_trace %>%
  select(Ftrgt, catch, ssb, rec) %>%
  pivot_longer(cols = c("catch", "ssb", "rec")) %>%
  mutate(value = value/1000,
         name = factor(name, levels = c("catch", "ssb", "rec"), 
                       labels = c("Catch [1000t]", "SSB [1000t]",
                                  "Recruitment [millions]"))) %>%
  ggplot(aes(x = Ftrgt, y = value)) +
  geom_point(size = 0.8) +
  stat_smooth(aes(alpha = "loess smoother"), size = 0.5,
              se = FALSE, span = 0.4, n = 100, show.legend = TRUE) + 
  scale_alpha_manual("", values = 1) +
  facet_wrap(~ name, scales = "free_y", strip.position = "left") +
  labs(x = "F (ages 3-6)", y = "") +
  ylim(c(0, NA)) +
  scale_x_continuous(breaks = seq(0, 1, 0.2)) +
  theme_bw() +
  theme(legend.position = c(0.85, 0.2),
        legend.background = element_blank(),
        legend.key = element_blank(),
        strip.placement = "outside",
        strip.background = element_blank(),
        axis.title.y = element_blank())
ggsave("input/ple.27.7e/baseline/1000_100/MSY_search.png", 
       width = 17, height = 6, units = "cm", dpi = 300, type = "cairo")
ggsave("input/ple.27.7e/baseline/1000_100/MSY_search.pdf", 
       width = 17, height = 6, units = "cm")

### ------------------------------------------------------------------------ ###
### run MSY projection ####
### ------------------------------------------------------------------------ ###
### Fmsy ~ 0.1638845
res_trace[which.max(res_trace$catch), ]
Fmsy <- res_trace$Ftrgt[which.max(res_trace$catch)]

input_MSY <- input_constF
input_MSY$ctrl$hcr@args$ftrg <- Fmsy
input_MSY$cut_hist <- FALSE
res_mp_MSY <- do.call(mp, input_MSY)
saveRDS(res_mp_MSY, "input/ple.27.7e/baseline/1000_100/MSY_mp.rds")

median(tail(catch(res_mp_MSY@stock), 10), na.rm = TRUE)
# 1702.92
median(tail(ssb(res_mp_MSY@stock), 10), na.rm = TRUE)
# 9536.053
### check risk
mean(tail(ssb(res_mp_MSY@stock), 10) < 2110)
# 0 - no issue

### final MSY reference points:
### Fmsy = 0.164
### Bmsy = 9536 (t)
### Cmsy = 1703 (t)

