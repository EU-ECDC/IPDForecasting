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

PCVdates <- read_csv("VaccineDates.csv", col_names=TRUE)   # Read in dates of PCV introduction

## Read in and reformat IPD incidence data
IPDData <- read_csv("190412IPD.csv", col_names=TRUE) %>%
		   mutate(
				  ## Age groups
				  Age = as.integer(Age),
				  AgeGroup = case_when(Age < 5 ~ 1,		# 0-4
							  Age >= 5  & Age < 15 ~ 2,	# 5-14
							  Age >= 15 & Age < 65 ~ 3,	# 15-64
							  Age >= 65 & Age < 85 ~ 4,	# 65-84
							  Age >= 85 ~ 5),
							  
				  ## Countries
				  ReportingCountry = case_when(GeoCode == "UKM" ~ "SC"), # Recode region 'UKM' as Scotland 
							  
				  ## Serotypes
				  Serotype = recode(Serotype, "NT" = "UNK",           		 
											  "NTYP" = "UNK",		   
											  "O" = "OTHER"),
				  GroupType = case_when(Serotype == "4"  ~ "PCV7",    # Group serotypes according to vaccine 
										Serotype == "6B" ~ "PCV7",
										Serotype == "9V" ~ "PCV7",
										Serotype == "14"  ~ "PCV7",
										Serotype == "18C" ~ "PCV7",
										Serotype ==  "19F" ~ "PCV7",
										Serotype == "23F" ~ "PCV7",
										Serotype == "1" ~ "PCV10",
										Serotype == "5" ~ "PCV10",
										Serotype == "7F" ~ "PCV10",
										Serotype == "3" ~ "PCV13",
										Serotype == "6A" ~ "PCV13",
										Serotype == "19A" ~ "PCV13",
										Serotype == "2"  ~ "PPV23",
										Serotype == "8" ~ "PPV23",
										Serotype == "9N" ~ "PPV23",
										Serotype == "10A" ~ "PPV23",
										Serotype == "11A" ~ "PPV23",
										Serotype == "12F" ~ "PPV23",
										Serotype == "15B"  ~ "PPV23",
										Serotype == "17F" ~ "PPV23",
										Serotype == "20" ~ "PPV23",
										Serotype == "22F" ~ "PPV23",		   
										Serotype == "OTHER" ~ "OTHER"), 
										
					## Specimen type				   
					Specimen = recode(Specimen, "JOINT"  = "OTHER",    # Group samples from joint, peritoneal and pleural fluid
                                     "PERIT" = "OTHER",
									 "PLEURAL" = "OTHER",
									 "O"  = "OTHER",
                                     "NULL" = "UNK"),
					## Dates
					Date = as.Date(DateUsedForStatisticsISO, "%Y-%m-%d"),  # Convert to date format
					Date = case_when(nchar(DateUsedForStatisticsISO)==7 ~ paste(DateUsedForStatisticsISO,"-01",sep=""),    # Where no day specified, set to 1st.
									 nchar(DateUsedForStatisticsISO)==4 ~ paste(DateUsedForStatisticsISO,"-01-01",sep=""), # Where no day or month specified, set to 1st of the 1st.
									 TRUE ~ as.character(Date)),
				    Date = ymd(Date),
				    Year = substr(Date,1,4)) %>%    # Define year for matching 
				   
			select(-DateUsedForStatisticsISO) %>%   # Remove original date variable
			
			## PCV programmes
			left_join(PCVdates) %>%										# Add dates of PCV introduction 
			 mutate(StartPCV1013 = as.Date(StartPCV1013, "%d/%m/%Y"),
				    StartPCV7 = as.Date(StartPCV7, "%d/%m/%Y"),
					
				    TimePeriod = case_when(Date < StartSpID1 ~ 1,				# Before the inclusion period of SpIDnet1
								Date > StartSpID1   & Date < StartPCV7 ~ 2,		# SpIDnet pre-PCV7 period
								Date > StartPCV7    & Date < StartPCV1013 ~ 3,	# PCV7 period
								Date > StartPCV1013 & Date < EndDateSpID1 ~ 4,	# PCV10/13 period
								Date > EndDateSpID1 ~ 5),
								
					PCV1013Yr  = if_else(TimePeriod == 4 | TimePeriod==5,
											pmax(ceiling((Date - StartPCV1013)/365.25),1),
											NA_real_))


## Include only confirmed cases
IPDData <- IPDData %>% filter(Classification == "CONF")  
				   

## Read in data for predictor variables
prop0to14 <- read_csv("prop0to14.csv", col_names=TRUE)   # Source: EUROSTAT https://ec.europa.eu/eurostat/databrowser/view/tps00010/default/table?lang=en
prop65plus <- read_csv("prop65plus.csv", col_names=TRUE) # Source: EUROSTAT https://ec.europa.eu/eurostat/databrowser/view/tps00010/default/table?lang=en
prop80plus <- read_csv("prop80plus.csv", col_names=TRUE) # Source: EUROSTAT https://ec.europa.eu/eurostat/databrowser/view/tps00010/default/table?lang=en
popProjections <- read_csv("popProjections.csv", col_names=TRUE) # Source: UN DESA World Population Prospects 2019 https://population.un.org/wpp/Download/Probabilistic/Population/
childcare0to2 <- read_csv("childcare0to2.csv", col_names=TRUE)   # Source: EUROSTAT https://ec.europa.eu/eurostat/databrowser/view/tps00185/default/table?lang=en
PCV7Data <- read_csv("PCV7coverage.csv", col_names=TRUE) # Read in WHO PCV7 coverage data
PCV13Data <- read_csv("PCV13coverage.csv", col_names=TRUE) # Read in WHO PCV13 coverage data
PPV23Data <- read_csv("PPV23coverage.csv", col_names=TRUE) # Read in PPV23 coverage data
PCV13AData <- read_csv("PCV13Acoverage.csv", col_names=TRUE) # Read in coverage data for PCV13 in adults

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


