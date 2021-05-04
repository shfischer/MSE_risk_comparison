library(ggplot2)
library(FLCore)
library(FLAssess)
library(FLXSA)
library(FLash)
library(FLfse)
library(ggplotFL)
library(stockassessment)
library(foreach)
library(dplyr)
library(tidyr)
library(mse)

source("funs.R")


rule <- "rfb"
n <- 500
m <- 0.8
# input_FLXSA <- readRDS(paste0("input/input_", n, "_FLXSA.rds"))
input <- readRDS(paste0("input/input_", n, "_", rule,  ".rds"))
input$ctrl$est@args$comp_m <- m
set.seed(1)
res <- do.call(mp, input)
saveRDS(res, file = paste0("output/res_", n, "_", rule, "_", m, ".rds"))

