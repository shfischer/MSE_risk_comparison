library(tidyverse)

### ------------------------------------------------------------------------ ###
### collate data ####
### ------------------------------------------------------------------------ ###

path <- "input/her.27.3a47d/preparation/HERAS/biotic/"
files <- list.files(path, pattern = "*.zip")
path_out <- "input/her.27.3a47d/preparation/HERAS/biotic_unzipped"
dir.create(path_out, recursive = TRUE)

### unzip archives
for (file in files) {
  #browser()
  unzip(zipfile = paste0(path, file), exdir = path_out)
}

### combine records
path <- "input/her.27.3a47d/preparation/HERAS/biotic_unzipped/"
files <- list.files(path_out, pattern = "*.csv")

bio_data <- lapply(files, function(file) {
  #browser()
  ### get file header
  header <- read.csv(file = paste0(path_out, "/", file), nrows = 1)
  ### find biological data
  data <- readLines(paste0(path_out, "/", file))
  pos_bio <- grep(x = data, pattern = "^Biology,*")[1]
  ### load biological data
  bio <- read.csv(file = paste0(path_out, "/", file), skip = pos_bio - 1)
  ret <- full_join(header, bio)
  ret$CruiseLocalID <- as.character(ret$CruiseLocalID)
  ret$BiologyIndividualMaturity <- as.character(ret$BiologyIndividualMaturity)
  return(ret)
})
bio_data <- bind_rows(bio_data)
bio_data <- bio_data %>%
  filter(BiologyStockCode == "her-47d3")

table(bio_data$HaulGear)
table(bio_data$CatchSpeciesCode)

summary(bio_data$BiologyLengthClass)
summary(bio_data$BiologyIndividualAge)

### keep only individuals with age data
bio_data <- bio_data %>%
  filter(!is.na(BiologyIndividualAge))
nrow(bio_data)

summary(bio_data)

### add year
table(lubridate::as_date(bio_data$CruiseStartDate))
bio_data <- bio_data %>% 
  mutate(year = as.numeric(str_sub(CruiseStartDate, start = 1, end = 4)))
table(bio_data$year)

### create ALK
alk_full <- bio_data %>%
  mutate(age = BiologyIndividualAge, length = (BiologyLengthClass)/10) %>%
  select(age, year, length) %>%
  mutate(length = floor(length*2)/2) ### round to 0.5cm
alk <- alk_full %>%
  ### remove age 0
  ### only data for 2017-2020 exist
  ### age 0 corresponds to small fish, max of 17cm, most <=10cm
  ### remove for more robust estimation of Lc
  filter(age >= 1) %>%
  mutate(age = ifelse(age > 8, 8, age)) %>% ### age 8 is plusgroup
  #mutate(age = ifelse(age < 1, 1, age)) %>%
  group_by(year, age, length) %>%
  count() %>%
  ungroup() %>%
  group_by(year, age) %>%
  mutate(freq = n/sum(n)) %>%
  select(-n) %>% 
  ungroup()
View(alk)

alk %>% 
  group_by(year) %>%
  count()

saveRDS(alk, file = "input/her.27.3a47d/preparation/ALK_MSE.rds")
alk <- readRDS("input/her.27.3a47d/preparation/ALK_MSE.rds")

### ------------------------------------------------------------------------ ###
### von Bertalanffy growth parameters ####
### ------------------------------------------------------------------------ ###

### von Bertalanffy growth function
vB_length <- function(L_inf, k, t, t_0 = 0){
  L_inf * (1 - exp(-k * (t - t_0)))
}

### fit to data from last 5 years
data_vB <- alk_full %>%
  filter(year %in% 2016:2020)

### fit model
vB_pars <- nls(length ~ vB_length(L_inf, k, t = age, t_0),
               data = as.data.frame(data_vB),
               start = list(L_inf = 31, k = 0.3, t_0 = 0))
# Nonlinear regression model
# model: length ~ vB_length(L_inf, k, t = age, t_0)
# data: as.data.frame(data_vB)
# L_inf       k     t_0 
# 30.6699  0.4880 -0.7057 
# residual sum-of-squares: 124520
# 
# Number of iterations to convergence: 4 
# Achieved convergence tolerance: 2.58e-06
summary(vB_pars)

### plot data and fit
alk_full %>%
  filter(year %in% 2016:2020) %>%
  group_by(age, length) %>%
  count() %>%
  ungroup() %>%
  group_by(age) %>%
  mutate(freq = n/sum(n)) %>%
  select(-n) %>% 
  ungroup() %>%
  ggplot(aes(x = age, y = length, size = freq)) +
  geom_point(shape = 21, stroke = 0.1) +
  geom_function(
    fun = vB_length, 
    args = list(k = summary(vB_pars)$parameters["k", "Estimate"],
                t_0 = summary(vB_pars)$parameters["t_0", "Estimate"],
                L_inf = summary(vB_pars)$parameters["L_inf", "Estimate"]),
    colour = "blue", show.legend = FALSE) +
  scale_size("Frequency", range = c(0.1, 5), breaks = c(0.1, 0.5)) + 
  scale_y_continuous(limits = c(0, 40)) +
  scale_x_continuous(limits = c(-1, 15)) +
  labs(x = "Age [years]", y = "Length [cm]") + 
  theme_bw(base_size = 8)
ggsave(filename = "input/her.27.3a47d/preparation/her_alk_vB.png", 
       width = 8.5, height = 5, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "input/her.27.3a47d/preparation/her_alk_vB.pdf",
       width = 8.5, height = 5, units = "cm")



### ------------------------------------------------------------------------ ###
### estimate Lc ####
### ------------------------------------------------------------------------ ###

### load input data
stk <- readRDS("input/her.27.3a47d/preparation/stk.rds")

data_length <- full_join(
  ### use observed catch
  x = as.data.frame(catch.n(stk)) %>%
    select(year, age, data) %>%
    rename("caa" = "data") %>%
    mutate(caa = caa + .Machine$double.eps), ### avoid 0s
  y = alk %>%
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
### recruitment peak at small lengths
full_join(
  ### use observed catch
  x = as.data.frame(catch.n(stk)) %>%
    select(year, age, data) %>%
    rename("caa" = "data") %>%
    mutate(caa = caa + .Machine$double.eps), ### avoid 0s
  y = alk %>%
    select(year, age, length, freq) %>%
    filter(age == 0), 
  by = c("year", "age")) %>%
  ### calculate numbers at length
  mutate(cal = caa * freq) %>%
  filter(!is.na(cal)) %>%
  group_by(year, length) %>%
  summarise(cal = sum(cal)) %>%
  ggplot(aes(x = length, y = cal)) +
  geom_col() +
  facet_wrap(~ year, scales = "free_y") +
  theme_bw(base_size = 8)
### recruits (age 0) lengths only available 2017-2020

### find Lc by year
Lc_values <- data_length %>%
  group_by(year) %>%
  summarise(Lc = min(length[cal >= max(cal)/2])) %>%
  ungroup()
Lc_values %>% print(n = Inf)

### remove recruits
### find Lc by year
Lc_values <- data_length %>%
  group_by(year) %>%
  summarise(Lc = min(length[cal >= max(cal)/2])) %>%
  ungroup()
Lc_values %>% print(n = Inf)
### use last 10 years => 25 cm
mean(tail(Lc_values$Lc, 5)) # 25.2
mean(tail(Lc_values$Lc, 10)) # 25.45
mean(Lc_values$Lc) # 25.13636
### -> use 25 cm

### ------------------------------------------------------------------------ ###
### historical mean length in catch ####
### ------------------------------------------------------------------------ ###

### mean catch length above Lc
data_length %>%
  filter(length > 20) %>%
  group_by(year) %>%
  summarise(length = weighted.mean(x = length, w = cal)) %>%
  print(n = Inf)

### combine all ALKs for years <2010 without ALK data
alk_combined <- alk_full %>%
  filter(age >= 1) %>% ### 0 not included in assessment
  mutate(age = ifelse(age > 8, 8, age)) %>% ### age 8 is plusgroup
  mutate(age = ifelse(age < 1, 1, age)) %>%
  group_by(year, age, length) %>%
  count() %>%
  ungroup() %>%
  group_by(year, age) %>%
  mutate(freq = n/sum(n)) %>%
  select(-n) %>% 
  ungroup() %>%
  group_by(age, length) %>%
  summarise(freq = sum(freq)) %>%
  ungroup() %>%
  group_by(age) %>%
  mutate(freq = freq/sum(freq)) %>%
  ungroup()

### check
alk_combined %>%
  group_by(age) %>%
  summarise(sum(freq))

### get catch length distribution
data_length_hist <- full_join(
  ### use observed catch
  x = as.data.frame(catch.n(stk)) %>%
    select(year, age, data) %>%
    rename("caa" = "data") %>%
    mutate(caa = caa + .Machine$double.eps), ### avoid 0s
  y = alk %>%
    select(age, length, freq), 
  by = c("age")) %>%
  ### calculate numbers at length
  mutate(cal = caa * freq) %>%
  filter(!is.na(cal)) %>%
  group_by(year, length) %>%
  summarise(cal = sum(cal))

### combine all data
data_length_full <- bind_rows(data_length_hist, data_length)
data_length_full %>%
  ggplot(aes(x = length, y = cal)) +
  geom_col() +
  facet_wrap(~ year, scales = "free_y") +
  theme_bw(base_size = 8)
### mean catch length above Lc
data_length_full %>%
  filter(length > 20) %>%
  group_by(year) %>%
  summarise(length = weighted.mean(x = length, w = cal)) %>%
  print(n = Inf)
### always above Lc (25cm) since 1997
# # A tibble: 72 x 2
# year length
# <dbl>  <dbl>
#   1  1947   28.6
# 2  1948   28.5
# 3  1949   28.4
# 4  1950   27.6
# 5  1951   27.3
# 6  1952   27.0
# 7  1953   26.6
# 8  1954   26.5
# 9  1955   25.6
# 10  1956   25.7
# 11  1957   26.0
# 12  1958   25.1
# 13  1959   25.1
# 14  1960   25.5
# 15  1961   26.1
# 16  1962   26.6
# 17  1963   24.7
# 18  1964   25.2
# 19  1965   25.9
# 20  1966   25.8
# 21  1967   26.0
# 22  1968   25.4
# 23  1969   24.4
# 24  1970   24.7
# 25  1971   23.9
# 26  1972   23.7
# 27  1973   24.3
# 28  1974   24.7
# 29  1975   23.9
# 30  1976   24.6
# 31  1977   25.3
# 32  1980   24.8
# 33  1981   24.1
# 34  1982   24.0
# 35  1983   23.8
# 36  1984   24.6
# 37  1985   24.7
# 38  1986   24.4
# 39  1987   24.1
# 40  1988   24.1
# 41  1989   25.2
# 42  1990   25.4
# 43  1991   25.6
# 44  1992   25.6
# 45  1993   25.1
# 46  1994   25.0
# 47  1995   25.0
# 48  1996   24.7
# 49  1997   25.6
# 50  1998   25.2
# 51  1999   25.9
# 52  2000   25.8
# 53  2001   25.9
# 54  2002   25.9
# 55  2003   26.1
# 56  2004   27.0
# 57  2005   27.4
# 58  2006   27.8
# 59  2007   27.9
# 60  2008   27.2
# 61  2009   26.9
# 62  2010   26.4
# 63  2011   26.7
# 64  2012   27.0
# 65  2013   27.6
# 66  2014   27.5
# 67  2015   27.7
# 68  2016   27.4
# 69  2017   27.3
# 70  2018   28.1
# 71  2019   28.3
# 72  2020   27.7

