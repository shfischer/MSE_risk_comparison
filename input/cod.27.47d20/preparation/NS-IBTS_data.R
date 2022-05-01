### ------------------------------------------------------------------------ ###
### IBTS weights and ALKs ####
### ------------------------------------------------------------------------ ###


library(tidyverse)

### ------------------------------------------------------------------------ ###
### weights at age ####
### ------------------------------------------------------------------------ ###


smalk <- read.csv("input/cod.27.47d20/preparation/DATRAS/NS-IBTS SMALK/SMALK    _2021-11-22 13_53_29.csv")

IBTS_weights <- smalk %>%
  filter(!is.na(Age) & !is.na(IndWgt)) %>%
  group_by(Survey, Year, Quarter, Age) %>%
  #summarise(weight = mean(IndWgt))
  summarise(weight = weighted.mean(x = IndWgt, w = CANoAtLngt, na.rm = TRUE))

  
# blubb <- smalk %>%
#     filter(!is.na(Age) & !is.na(IndWgt)) %>%
#     group_by(Survey, Year, Quarter, Age) %>%
#   filter(Year == 1991 & Quarter == 1 & Age == 1)
# weighted.mean(x = blubb$IndWgt, w = blubb$CANoAtLngt)

### load indices
idcs <- readRDS("input/cod.27.47d20/preparation/idx.rds")
range(idcs$IBTS_Q1_gam)
range(idcs$IBTS_Q3_gam)

### extract Q1 weights
IBTS_weights_Q1 <- IBTS_weights %>%
  ungroup() %>%
  filter(Quarter == 1 &
           Year %in% 1983:2021 &
           Age %in% 1:5) %>%
  select(Year, Age, weight) %>%
  rename(Year = Year, age = Age, data = weight) %>%
  mutate(data = data/1000) # in kg
### convert into FLQuant
IBTS_weights_Q1 <- window(as.FLQuant(IBTS_weights_Q1), start = 1983)
### weights missing for 1983-1996 & 1998-1999
### fill with averages
IBTS_weights_Q1[, ac(c(1983:1996, 1998:1999))] <- yearMeans(IBTS_weights_Q1)
plot(IBTS_weights_Q1)


### extract Q3 weights
IBTS_weights_Q3 <- IBTS_weights %>%
  ungroup() %>%
  filter(Quarter == 3 &
           Year %in% 1992:2020 &
           Age %in% 1:4) %>%
  select(Year, Age, weight) %>%
  rename(Year = Year, age = Age, data = weight) %>%
  mutate(data = data/1000) # in kg
### convert into FLQuant
IBTS_weights_Q3 <- window(start = 1992, as.FLQuant(IBTS_weights_Q3))
### weights missing for 1992-1998 & 2000
### fill with averages
IBTS_weights_Q3[, ac(c(1992:1998, 2000))] <- yearMeans(IBTS_weights_Q3)
plot(IBTS_weights_Q3)

### save weights
saveRDS(IBTS_weights_Q1, 
        file = "input/cod.27.47d20/preparation/idx_weights_Q1.rds")
saveRDS(IBTS_weights_Q3, 
        file = "input/cod.27.47d20/preparation/idx_weights_Q3.rds")

### ------------------------------------------------------------------------ ###
### ALKs ####
### ------------------------------------------------------------------------ ###

IBTS_alk <- read.csv("input/cod.27.47d20/preparation/DATRAS/NS-IBTS ALK/ALK_2021-11-22 14_25_14.csv")


IBTS_alk <- IBTS_alk %>%
  pivot_longer(Age_0:Age_10, names_prefix = "Age_", names_to = "age",
               values_to = "count") %>%
  mutate(age = as.numeric(age)) %>% 
  mutate(LngtClass = LngtClass/10) %>% # in cm
  filter(!is.na(count)) 

IBTS_alk %>%
  filter(Quarter == 1) %>%
  filter(Year >= 1990) %>%
  ggplot(aes(x = age, y = LngtClass, size = count)) +
  geom_point(shape = 1) +
  facet_wrap(~ Year)
IBTS_alk %>%
  filter(Quarter == 3) %>%
  filter(Year >= 1990) %>%
  ggplot(aes(x = age, y = LngtClass, size = count)) +
  geom_point(shape = 1) +
  facet_wrap(~ Year)

### combine Q1 & Q3 into single ALK - all ages
IBTS_alk_full <- IBTS_alk %>%
  #filter(Quarter == 3 & age == 3 & Year == 2020) %>%
  select(Year, Quarter, LngtClass, age, count) %>%
  rename(year = Year, length = LngtClass) %>%
  group_by(Quarter, year, age) %>%
  mutate(freq = count/sum(count)) %>% ### freq per Quarter
  ungroup(Quarter) %>%
  select(-Quarter, -count) %>% 
  group_by(year, age, length) %>%
  summarise(freq = sum(freq)) %>% ### combine Q1 & Q3
  ungroup(length) %>%
  mutate(freq = freq/sum(freq)) ### re-standardise

IBTS_alk_full %>%
  filter(year >= 1990) %>%
  ggplot(aes(x = age, y = length, size = freq)) +
  geom_point(shape = 1) +
  facet_wrap(~ year)

IBTS_alk_full %>%
  ungroup() %>% group_by(year, age) %>%
  summarise(sum = sum(freq)) %>% summary()

### combine Q1 & Q3 into single ALK - OM ages (1-6)
#IBTS_alk <- 
IBTS_alk <- IBTS_alk %>%
  #filter(Quarter == 3 & age == 3 & Year == 2020) %>%
  filter(age >= 1) %>%
  mutate(age = ifelse(age > 6, 6, age)) %>% ### age 6 is plusgroup
  mutate(age = ifelse(age < 1, 1, age)) %>%
  select(Year, Quarter, LngtClass, age, count) %>%
  rename(year = Year, length = LngtClass) %>%
  group_by(Quarter, year, age, length) %>%
  summarise(count = sum(count)) %>% ### combine age 6+
  ungroup() %>%
  group_by(Quarter, year, age) %>%
  mutate(freq = count/sum(count)) %>% ### standardise by age - per Quarter
  ungroup(Quarter) %>%
  select(-Quarter, -count) %>% 
  group_by(year, age, length) %>%
  summarise(freq = sum(freq)) %>% ### combine Q1 & Q3
  ungroup(length) %>%
  mutate(freq = freq/sum(freq)) %>%### re-standardise
  ungroup()
IBTS_alk %>%
  filter(year >= 1990) %>%
  ggplot(aes(x = age, y = length, size = freq)) +
  geom_point(shape = 1, show.legend = TRUE) +
  scale_size(breaks = c(0.25, 0.5, 0.75, 1)) +
  facet_wrap(~ year) +
  xlim(c(0, 10)) +
  theme(legend.position = "right")

### save ALK for MSE
saveRDS(IBTS_alk, file = "input/cod.27.47d20/preparation/ALK_MSE.rds")


### ------------------------------------------------------------------------ ###
### von Bertalanffy growth parameters ####
### ------------------------------------------------------------------------ ###

### von Bertalanffy growth function
vB_length <- function(L_inf, k, t, t_0 = 0){
  L_inf * (1 - exp(-k * (t - t_0)))
}

### fit to data from last 5 years
### (exclude 2021 because Q3 is missing)
data_vB <- IBTS_alk_full %>%
  filter(year %in% 2016:2020) %>%
  select(-year) %>%
  ungroup() %>%
  group_by(age, length) %>%
  summarise(freq = sum(freq)) %>%
  ungroup(length) %>%
  mutate(freq = freq/sum(freq)) %>%
  group_by(age, length) %>%
  summarise(freq = sum(freq))
# data_vB %>%
#   summarise(freq = sum(freq))


### fit model
### for each data year
vB_pars <- nls(length ~ vB_length(L_inf, k, t = age, t_0),
      data = as.data.frame(data_vB), weights = data_vB$freq,
      start = list(L_inf = 80, k = 0.1, t_0 = 0))
summary(vB_pars)
vB_pars
# Nonlinear regression model
# model: length ~ vB_length(L_inf, k, t = age, t_0)
# data: as.data.frame(data_vB)
# L_inf        k      t_0 
# 116.9165   0.1967  -0.2965 
# weighted residual sum-of-squares: 984.5
# 
# Number of iterations to convergence: 6 
# Achieved convergence tolerance: 7.316e-06

### plot data and fit
data_vB %>%
  ggplot(aes(x = age, y = length, size = freq)) +
  geom_point(shape = 21, stroke = 0.1) +
  geom_function(fun = vB_length, 
                args = list(k = summary(vB_pars)$parameters["k", "Estimate"],
                            t_0 = summary(vB_pars)$parameters["t_0", "Estimate"],
                            L_inf = summary(vB_pars)$parameters["L_inf", "Estimate"]),
                colour = "blue", show.legend = FALSE) +
  scale_size("Frequency", range = c(0.1, 5)) + 
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100, 125), limits = c(0, 125)) +
  scale_x_continuous(breaks = seq(0, 10, 2), limits = c(-0.5, 10.5)) +
  labs(x = "Age [years]", y = "Length [cm]") + 
  theme_bw(base_size = 8)
ggsave(filename = "input/cod.27.47d20/preparation/DATRAS/cod_alk_vB.png", 
       width = 8.5, height = 5, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "input/cod.27.47d20/preparation/DATRAS/cod_alk_vB.pdf",
       width = 8.5, height = 5, units = "cm")


### ------------------------------------------------------------------------ ###
### estimate Lc ####
### ------------------------------------------------------------------------ ###

### load input data
stk <- readRDS("input/cod.27.47d20/preparation/stk.rds")

data_length <- full_join(
  ### use observed catch
  x = as.data.frame(catch.n(stk)) %>%
    select(year, age, data) %>%
    rename("caa" = "data") %>%
    mutate(caa = caa + .Machine$double.eps), ### avoid 0s
  y = IBTS_alk %>%
    select(year, age, length, freq), 
  by = c("year", "age")) %>%
  ### calculate numbers at length
  mutate(cal = caa * freq) %>%
  filter(!is.na(cal)) %>%
  group_by(year, length) %>%
  summarise(cal = sum(cal))
data_length %>%
  ggplot(aes(x = length, y = cal)) +
  geom_col() +
  facet_wrap(~ year, scales = "free_y") +
  theme_bw(base_size = 8)

### find Lc by year
Lc_values <- data_length %>%
  group_by(year) %>%
  summarise(Lc = min(length[cal >= max(cal)/2])) %>%
  ungroup()
Lc_values %>% print(n = Inf)
### use last 10 years => 20 cm
mean(tail(Lc_values$Lc, 5)) # 19.2
mean(tail(Lc_values$Lc, 10)) # 19.8
mean(tail(Lc_values$Lc, 20)) # 19.9
mean(Lc_values$Lc) # 17.9


### mean catch length above Lc
data_length %>%
  filter(length > 20) %>%
  group_by(year) %>%
  summarise(length = weighted.mean(x = length, w = cal)) %>%
  print(n = Inf)
