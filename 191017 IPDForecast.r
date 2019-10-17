rm(list=ls())
dir <- "P:/Pneumo"
setwd(dir)
library(plyr)
library(tidyverse)
library(lubridate)
library(reshape2)
library(timetk)
library(sweep)
library(forecast)

## Read in IPD incidnce data
data1 <- read_csv("190412IPD.csv", col_names=TRUE)
data1$Age <- as.integer(data1$Age)


## Read in data for predictor variables
prop0to14 <- read_csv("prop0to14.csv", col_names=TRUE)   ## Source: EUROSTAT https://ec.europa.eu/eurostat/databrowser/view/tps00010/default/table?lang=en
prop65plus <- read_csv("prop65plus.csv", col_names=TRUE)   ## Source: EUROSTAT https://ec.europa.eu/eurostat/databrowser/view/tps00010/default/table?lang=en
prop80plus <- read_csv("prop80plus.csv", col_names=TRUE)   ## Source: EUROSTAT https://ec.europa.eu/eurostat/databrowser/view/tps00010/default/table?lang=en
popProjections <- read_csv("popProjections.csv", col_names=TRUE) ## Source: UN DESA World Population Prospects 2019 https://population.un.org/wpp/Download/Probabilistic/Population/
childcare0to2 <- read_csv("childcare0to2.csv", col_names=TRUE)   ## Source: EUROSTAT https://ec.europa.eu/eurostat/databrowser/view/tps00185/default/table?lang=en
PCV7Data <- read_csv("PCV7coverage.csv", col_names=TRUE) ## Read in WHO PCV7 coverage data
PCV13Data <- read_csv("PCV13coverage.csv", col_names=TRUE) ## Read in WHO PCV13 coverage data
PPV23Data <- read_csv("PPV23coverage.csv", col_names=TRUE) ## Read in PPV23 coverage data
PCV13AData <- read_csv("PCV13Acoverage.csv", col_names=TRUE) ## Read in coverage data for PCV13 in adults

EEA <- c("AT","BE","BG", "HR", "CY", "CZ", "DK", "EE", "FI", "FR", "DE", "EL", "HU", "IS", "IE", "IT", "LV", "LI", "LT", "LU", "MT", "NL", "NO", "PL", "PT", "RO", "SK", "SI", "ES", "SE", "UK") 
countryList <- c("Austria", "Belgium", "Bulgaria", "Croatia", "Cyprus", "Czechia", "Denmark", "Estonia", "Finland", "France", "Germany", "Greece", "Hungary", "Iceland", "Ireland", "Italy", "Latvia", "Liechtenstein", "Lithuania","Luxembourg", "Malta", "Netherlands", "Norway", "Poland", "Portugal", "Romania", "Slovakia", "Slovenia", "Spain", "Sweden", "UK" )
countries <- tibble(code=EEA, name=countryList)

## Reformat predictor variables

# Childcare	
childcare0to2 <- gather(childcare0to2, countries$name, key="country", value="childcare0to2") %>%
				 select(country, year, childcare0to2) %>%
				 mutate(date = as.Date(ISOdate(year, 12, 31))) # set date to 31st December of given year


# Population projections
# Interpolate for single year predictions

# Timeframe for projections
years <- seq(min(popProjections$year), max(popProjections$year))

# Generate rows with missing numbers
popProjections <- popProjections %>% mutate(year = length(years)) %>% 
								  group_by(country) %>% 
								  expand(year = years) %>%
								  left_join(popProjections)
								  
# Interpolate with splines


