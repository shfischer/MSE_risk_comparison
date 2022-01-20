### ------------------------------------------------------------------------ ###
### Preparation of length data ####
### ------------------------------------------------------------------------ ###
library(icesDatras)
library(tidyr)
library(dplyr)
library(ggplot2)
library(FLCore)

### ------------------------------------------------------------------------ ###
### Q1SWBeam ####
### ------------------------------------------------------------------------ ###
### use data age-length data from Q1SWBeam to supplement landings ALKs


if (FALSE) {
  ### get data from DATRAS
  Q1 <- getCAdata(survey = "BTS", year = 2020, quarter = 1)
  Q1 <- lapply(2006:2020, getCAdata, survey = "BTS", quarter = 1)
  Q1 <- do.call(bind_rows, Q1)
  saveRDS(Q1, "input/ple.27.7e/preparation/length/DATRAS_Q1SWBeam.rds")
  Q1 <- readRDS("input/ple.27.7e/preparation/length/DATRAS_Q1SWBeam.rds")
  
  ### select plaice
  Q1_ple <- Q1 %>%
    filter(Survey == "BTS" &
             Ship == "7.4e+10" & ### Cefas Endeavour
             Quarter == 1 &
             Country == "GB" &
             SpecCode == 127143 & ### plaice
             !is.na(Age)
    )
  saveRDS(object = Q1_ple, 
          file = "input/ple.27.7e/preparation/length/Q1SWBeam_ple.rds")
}

Q1_ple <- readRDS("input/ple.27.7e/preparation/length/Q1SWBeam_ple.rds")
Q1_ple <- Q1_ple %>%
  select(Year, Age, LngtClass, CANoAtLngt) %>%
  group_by(Year, Age, LngtClass) %>%
  summarise(CANoAtLngt = sum(CANoAtLngt))

### use 2cm length bins - same for landings ALK
Q1_ple <- Q1_ple %>% 
  mutate(LngtClass = ifelse((LngtClass %% 2) == 0, LngtClass - 1, LngtClass))
table(Q1_ple$LngtClass)

### check number of fish ages per year
Q1_ple %>% 
  group_by(Year) %>%
  summarise(count = sum(CANoAtLngt))
#     Year count
#  1  2006   276
#  2  2007   269
#  3  2008   200
#  4  2009   329
#  5  2010   375
#  6  2011   541
#  7  2012   394
#  8  2013   522
#  9  2014   938
# 10  2015   812
# 11  2016   870
# 12  2017   902
# 13  2018   835
# 14  2019   669
# 15  2020   371

### standardise column names to facilitate merging later
Q1_ple <- Q1_ple %>% 
  rename(year = Year, age = Age, length = LngtClass, count = CANoAtLngt)

### ------------------------------------------------------------------------ ###
### ALKs from commercial UK landings ####
### ------------------------------------------------------------------------ ###
### retrieved from ICES WGCSE accessions

### load ALKS (historical and current)
alks <- read.csv("input/ple.27.7e/preparation/length/ALKs.csv", 
                 stringsAsFactors = FALSE, as.is = TRUE)
names(alks)[1] <- "year"
alks <- alks %>% gather(key = "age", value = "count", X1:X26) %>%
  mutate(age = as.numeric(gsub(x = age, pattern = "X", replacement = ""))) %>%
  filter(!is.na(count))

### 1cm length bins for 2017, remaining years with 2cm bins -> standardise
alks <- alks %>% 
  mutate(length = ifelse((length %% 2) == 0, length - 1, length))
table(alks$length)

### check fish aged per year
alks %>% 
  group_by(year) %>%
  summarise(count = sum(count))
# year count
# 2013  2149
# 2014  2616
# 2015  2124
# 2016  1998
# 2017  2525
# 2018  1675
# 2019  2331
# 2020  1992

### check fish aged per year and quarter
alks %>% 
  group_by(year, quarter) %>%
  summarise(count = sum(count)) %>%
  print(n = Inf)

### summarise counts (remove season, sex, etc)
alks <- alks %>%
  filter(data == "landings") %>% ### discards only available for 2017 - remove
  group_by(year, age, length) %>%
  summarise(count = sum(count)) %>%
  ungroup() %>%
  filter(count > 0)

### ------------------------------------------------------------------------ ###
### combine landings & Q1SWBeam ALKs ####
### ------------------------------------------------------------------------ ###

range(alks$year)
range(Q1_ple$year)

### standardise year range
Q1_ple <- Q1_ple %>%
  filter(year >= min(alks$year))

### combine 
ALKs <- bind_rows(alks %>% mutate(data = "landings"),
                  Q1_ple %>% mutate(data = "Q1SWBeam"))

### standardise length frequencies by age - 1st by data source
ALKs <- ALKs %>%
  group_by(data, year, age) %>%
  mutate(freq = count / sum(count))
### then - combine data source and standardise again
ALKs <- ALKs %>%
  ungroup() %>%
  group_by(year, age) %>%
  mutate(freq = freq/sum(freq))
### remove redundant columns
ALKs <- ALKs %>%
  ungroup() %>%
  group_by(year, age, length) %>%
  summarise(freq = sum(freq))

### check total freq per year & age
ALKs %>% 
  group_by(year, age) %>%
  summarise(total = sum(freq)) %>%
  summary()

### trim ages for OM 2-10
ALKs <- ALKs %>%
  mutate(age = ifelse(age < 2, 2, age)) %>%
  mutate(age = ifelse(age > 10, 10, age))
### standardise again
ALKs <- ALKs %>%
  ungroup() %>%
  group_by(year, age) %>%
  mutate(freq = freq/sum(freq))

### save ALK for MSE
saveRDS(ALKs, file = "input/length/ALK_MSE.rds")

### ------------------------------------------------------------------------ ###
### check reproduction of length frequencies for historical period ####
### ------------------------------------------------------------------------ ###

stk <- readRDS("input/ple.27.7e/preparation/model_input_stk_d.RDS")
cn <- as.data.frame(catch.n(stk)[, ac(2013:2020)])

### catch numbers at age
cn <- cn %>%
  select(year, age, iter, data) %>%
  rename("caa" = "data") %>%
  mutate(iter = as.numeric(as.character(iter)))
### match ALKs with data
cn <- left_join(cn, ALKs)
### calculate numbers at length
cn <- cn %>%
  mutate(cal = caa * freq)


p <- cn %>% group_by(year, length) %>%
  summarise(cal = sum(cal)) %>%
  ggplot(aes(x = length, y = cal)) +
  geom_line() +
  facet_wrap(~ year)
p

### mean catch length above Lc
Lc <- 26
cn %>%
  filter(length >= Lc) %>%
  #filter(!is.na(length) & !is.na(cal)) %>%
  group_by(year) %>%
  summarise(L = weighted.mean(x = length, w = cal))


### sample
cn %>%
  group_by(year, iter) %>%
  filter(length >= Lc) %>%
  summarise(data = mean(sample(x = length, prob = cal, 
                               size = 2000, replace = TRUE))) %>%
  arrange(as.numeric(as.character(iter)))

