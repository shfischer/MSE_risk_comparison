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

input_constF <- readRDS("input/ple.27.7e/baseline/1000_100/input_constF.rds")

### include parallelisation
input_constF$args$nblocks <- 10

### ------------------------------------------------------------------------ ###
### function for running MP and returning catch ####
### ------------------------------------------------------------------------ ###
### return median catch of last 10 years from 100-year projection

res_trace <- list()
mp_catch <- function(input, Ftrgt, minimise = FALSE) {
  input$ctrl$hcr@args$ftrg <- Ftrgt
  res_i <- do.call(mp, input)
  catch_i <- median(tail(catch(res_i@stock), 10), na.rm = TRUE)
  cat(paste0("Ftrgt=", Ftrgt, " -- catch=", catch_i, "\n"))
  res_trace <<- append(res_trace, list(list(Ftrgt = Ftrgt, catch = catch_i)))
  if (isTRUE(minimise)) catch_i <- -catch_i
  return(catch_i)
}

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
set.seed(1)
res_trace <- list()
res_optimise_MSY <- optimise(f = mp_catch, input = input_constF,
                             interval = c(0, 1), 
                             lower = 0, upper = 1,
                             maximum = TRUE,
                             tol = 0.01)
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
  ggplot(aes(x = Ftrgt, y = catch)) +
  geom_point(size = 0.8) +
  stat_smooth(aes(alpha = "loess smoother"), size = 0.5,
              se = FALSE, n = 100, span = 0.2, show.legend = TRUE) + 
  scale_alpha_manual("", values = 1) +
  labs(x = "F (ages 3-6)", y = "catch [tonnes]") +
  theme_bw() +
  theme(legend.position = c(0.75, 0.85),
        legend.background = element_blank())
ggsave("input/ple.27.7e/baseline/1000_100/MSY_search.png", 
       width = 8.5, height = 6, units = "cm", dpi = 300, type = "cairo")
ggsave("input/ple.27.7e/baseline/1000_100/MSY_search.pdf", 
       width = 8.5, height = 6, units = "cm")

### ------------------------------------------------------------------------ ###
### run MSY projection ####
### ------------------------------------------------------------------------ ###
### Fmsy = 0.18

input_MSY <- input_constF
input_MSY$ctrl$hcr@args$ftrg <- 0.18
input_MSY$cut_hist <- FALSE
res_mp_MSY <- do.call(mp, input_MSY)
saveRDS(res_mp_MSY, "input/ple.27.7e/baseline/1000_100/MSY_mp.rds")

median(tail(catch(res_mp_MSY@stock), 10), na.rm = TRUE)
# 1666.445
median(tail(ssb(res_mp_MSY@stock), 10), na.rm = TRUE)
# 8543.38
### check risk
mean(tail(ssb(res_mp_MSY@stock), 10) < 2110)
# 0 - no issue

### final MSY reference points:
### Fmsy = 0.18
### Bmsy = 8543 (t)
### Cmsy = 1666 (t)

