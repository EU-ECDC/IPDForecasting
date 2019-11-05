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

myColours <- c("26 107 133", "241 214 118", "168 45 23")
ECDCcol <- sapply(strsplit(myColours, " "), function(x)
    rgb(x[1], x[2], x[3], maxColorValue=255))  # convert to hexadecimal

colours1 <- colorRampPalette(ECDCcol)(11)

## List countries
EEA <- c("AT","BE","BG", "HR", "CY", "CZ", "DK", "EE", "FI", "FR", "DE", "EL", "HU", "IS", "IE", "IT", "LV", "LI", "LT", "LU", "MT", "NL", "NO", "PL", "PT", "RO", "SK", "SI", "ES", "SE", "UK") 
countryList <- c("Austria", "Belgium", "Bulgaria", "Croatia", "Cyprus", "Czechia", "Denmark", "Estonia", "Finland", "France", "Germany", "Greece", "Hungary", "Iceland", "Ireland", "Italy", "Latvia", "Liechtenstein", "Lithuania","Luxembourg", "Malta", "Netherlands", "Norway", "Poland", "Portugal", "Romania", "Slovakia", "Slovenia", "Spain", "Sweden", "UK" )
countries <- tibble(code=EEA, name=countryList)
countriesScot <- add_row(countries, code = "SC", name = "Scotland")

########################		
## IPD incidence data ##
########################

PCVdates <- read_csv("VaccineDates.csv", col_names=TRUE)   # Read in dates of PCV introduction

IPDData <- read_csv("190412IPD.csv", col_names=TRUE) %>%
		   mutate(
				  ## Age groups. N.B. These must correspond with population groups (see above)
				  Age = as.integer(Age),
				  AgeGroup = case_when(Age < 5 ~ 1,		# 0-4
							  Age >= 5  & Age < 15 ~ 2,	# 5-14
							  Age >= 15 & Age < 50 ~ 3,	# 15-49
							  Age >= 50 & Age < 65 ~ 4,	# 50-64
							  Age >= 65 & Age < 75 ~ 5,	# 65-74
							  Age >= 75 & Age < 85 ~ 6,	# 75-84
							  Age >= 85 ~ 7) %>% 
							  factor(labels = c("0-4", "5-14", "15-49", "50-64", "65-74", "75-84", "85+")),
							  
				  ## Countries
				  ReportingCountry = case_when(GeoCode == "UKM" ~ "SC",
											    TRUE ~ as.character(ReportingCountry)), # Recode region 'UKM' as Scotland 
				 
				 ## Serotypes
				  Serotype = recode(Serotype, "NT" = "UNK",           		 
											  "NTYP" = "UNK",		   
											  "O" = "OTHER") %>%
									as_factor(),
											  
				  GroupType = case_when(Serotype == "4"  ~ "PCV-7",    # Group serotypes according to vaccine 
										Serotype == "6B" ~ "PCV-7",
										Serotype == "9V" ~ "PCV-7",
										Serotype == "14"  ~ "PCV-7",
										Serotype == "18C" ~ "PCV-7",
										Serotype ==  "19F" ~ "PCV-7",
										Serotype == "23F" ~ "PCV-7",
										Serotype == "1" ~ "PCV-10",
										Serotype == "5" ~ "PCV-10",
										Serotype == "7F" ~ "PCV-10",
										Serotype == "3" ~ "PCV-13",
										Serotype == "6A" ~ "PCV-13",
										Serotype == "19A" ~ "PCV-13",
										Serotype == "2"  ~ "PPV-23",
										Serotype == "8" ~ "PPV-23",
										Serotype == "9N" ~ "PPV-23",
										Serotype == "10A" ~ "PPV-23",
										Serotype == "11A" ~ "PPV-23",
										Serotype == "12F" ~ "PPV-23",
										Serotype == "15B"  ~ "PPV-23",
										Serotype == "17F" ~ "PPV-23",
										Serotype == "20" ~ "PPV-23",
										Serotype == "22F" ~ "PPV-23",
										Serotype == "OTHER" ~ "Other") %>%
										as_factor() %>%
										fct_explicit_na(na_level = "Not serotyped"), 
										
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
				    
					Year = substr(Date,1,4) %>% as.integer()) %>%    # Define year for matching 
				   
			select(-DateUsedForStatisticsISO)   # Remove original date variable
			
firstYear <- IPDData %>% select(Year) %>% min()

#####################
## Population data ##
#####################

tessyPop <- read_csv("TESSy population.csv", col_names=TRUE) %>% # Read in TESSy population data
		    filter(ReportingCountry %in% EEA & ReportYear >= firstYear) %>%
			mutate(AgeGroup = case_when(AgeGroupId == 1   | AgeGroupId == 2  ~ 1,      # 0-4
										AgeGroupId >= 3   & AgeGroupId <= 5  ~ 2,      # 5-14
										AgeGroupId >= 6   & AgeGroupId <= 12 ~ 3,      # 15-49
										AgeGroupId >= 13  & AgeGroupId <= 15 ~ 4,      # 50-64
										AgeGroupId >= 16  & AgeGroupId <= 17 ~ 5,  	   # 65-74
										AgeGroupId >= 18  & AgeGroupId <= 19 ~ 6,      # 75-84
										AgeGroupId >= 20  & AgeGroupId <= 23 ~ 7) %>% # 85+
										factor(labels = c("0-4", "5-14", "15-49", "50-64", "65-74", "75-84", "85+"))) %>%
			filter(!is.na(AgeGroup))  %>% # Original population database included other ways of grouping ages. Drop these values	
			rename(Year = ReportYear) %>% # Rename to match with IPDData
			group_by(Year, ReportingCountry, AgeGroup) %>% # Group population by age groups
			summarise(Population = sum(Population))
	
# tessyPop %>% group_by(AgeGroupId, AgeGroup) %>% summarise() %>% print(n=63) #Table of TESSy age categories

eurostatPop <- read_csv("EUROSTATpop.csv", col_names=TRUE)   # Read in EUROSTAT population data
scotPop <- read_csv("ScotPop.csv", col_names=TRUE)   # Source: https://www.nrscotland.gov.uk/statistics-and-data/statistics/statistics-by-theme/population/population-estimates/mid-year-population-estimates

###############
## Incidence ##
###############

IPDData <- IPDData %>% filter(Classification == "CONF") %>% # Include only confirmed cases
						left_join(tessyPop) %>%
						mutate(ReportingCountry = as_factor(ReportingCountry)) %>%
						arrange(ReportingCountry, Year, AgeGroup) %>%
						filter(SubsetEpi == 1) # Select subset for epi. to avoid double-counting lab.results

		
incData <- IPDData %>% 	filter(Age >= 50)  %>%
						group_by(ReportingCountry, Year, AgeGroup, GroupType, Population) %>%
						#group_by(Year, AgeGroup, GroupType, Population) %>%
						summarise(Total = n()) %>% # Number of cases per year by age group, type and country
						mutate(Incidence = (Total/Population)*100000) # Incidence per 100,000	


## Incidence by country
  ggplot(data = incData, aes(Year, Incidence)) +
  geom_line(lwd = 0.7) +
  labs(x = "Year", y = "Incidence per 100,000") +
  facet_grid(GroupType ~ AgeGroup, labeller = label_wrap_gen(width = 2, multi_line = TRUE) ,  scales = "free_y") + 
  scale_x_continuous(breaks = seq(2007, 2017, by = 2)) +
  theme(panel.grid.major = element_blank(), 
					  panel.grid.minor = element_blank(), 
					  panel.background = element_blank(),
					  panel.spacing.x = unit(1, "lines"),
					  axis.line = element_line(colour = "black"),
					  text = element_text(size=12),
					  axis.text.x = element_text(size = rel(0.8)),
					  axis.text.x.top = element_text(vjust = -0.5),
					  strip.text.x = element_text(size=12, colour="white"),
					  strip.text.y = element_text(size=12, colour="white"),
					  strip.background = element_rect(fill=ECDCcol[1])) +
					  guides(fill=guide_legend(title="Vaccine type")) +
  ggtitle(paste("Rate of IPD incidence in", countries$name[countryNo])) +
  theme(plot.title = element_text(hjust = 0.5))
  
############################
## Proportion of serotype ##
############################

# Proportion of outcome by serotype
propTypes <- incData %>% filter(GroupType != "Not serotyped") %>% # 
			   group_by(ReportingCountry, Year, AgeGroup, GroupType) %>%
               dplyr::summarise(n = sum(Total)) %>%
			   mutate(freq = n / sum(n))
			   
countryNo <- 9
propTypes1 <- propTypes %>% filter(ReportingCountry == countries$code[countryNo])

## Trends in proportion of types, by country
ggplot(data = propTypes1, aes(Year, freq)) +
  geom_line(lwd = 0.7) +
  labs(x = "Year", y = "Proportion of typed samples") +
  facet_grid(GroupType ~ AgeGroup, labeller = label_wrap_gen(width = 2, multi_line = TRUE)) + # ,  scales = "free_y") + 
  scale_x_continuous(breaks = seq(2007, 2017, by = 2)) +
  theme(panel.grid.major = element_blank(), 
					  panel.grid.minor = element_blank(), 
					  panel.background = element_blank(),
					  panel.spacing = unit(1, "lines"),
					  axis.line = element_line(colour = "black"),
					  text = element_text(size=12),
					  axis.text.x = element_text(size = rel(0.8)),
					  axis.text.x.top = element_text(vjust = -0.5),
					  strip.text.x = element_text(size=12, colour="white"),
					  strip.text.y = element_text(size=12, colour="white"),
					  strip.background = element_rect(fill=ECDCcol[1])) +
					  guides(fill=guide_legend(title="Vaccine type")) +
  ggtitle(paste("Proportion of typed samples (", countries$name[countryNo], ")")) +
  theme(plot.title = element_text(hjust = 0.5))
  

#########################
## Predictor variables ##
#########################

## Read in data on predictor variables
prop0to14 <- read_csv("prop0to14.csv", col_names=TRUE)   # Source: EUROSTAT https://ec.europa.eu/eurostat/databrowser/view/tps00010/default/table?lang=en
prop65plus <- read_csv("prop65plus.csv", col_names=TRUE) # Source: EUROSTAT https://ec.europa.eu/eurostat/databrowser/view/tps00010/default/table?lang=en
prop80plus <- read_csv("prop80plus.csv", col_names=TRUE) # Source: EUROSTAT https://ec.europa.eu/eurostat/databrowser/view/tps00010/default/table?lang=en
popProjections <- read_csv("popProjections.csv", col_names=TRUE) # Source: UN DESA World Population Prospects 2019 https://population.un.org/wpp/Download/Probabilistic/Population/
childcare0to2 <- read_csv("childcare0to2.csv", col_names=TRUE)   # Source: EUROSTAT https://ec.europa.eu/eurostat/databrowser/view/tps00185/default/table?lang=en
PCV7Data <- read_csv("PCV7coverage.csv", col_names=TRUE) # Read in WHO PCV7 coverage data
PCV13Data <- read_csv("PCV13coverage.csv", col_names=TRUE) # Read in WHO PCV13 coverage data
PPV23Data <- read_csv("PPV23coverage.csv", col_names=TRUE) # Read in PPV23 coverage data
PCV13AData <- read_csv("PCV13Acoverage.csv", col_names=TRUE) # Read in coverage data for PCV13 in adults

## Reformat predictor variables

# Childcare	
childcare0to2 <- gather(childcare0to2, countries$name, key="country", value="childcare0to2") %>%
				 select(country, year, childcare0to2) %>%
				 mutate(date = as.Date(ISOdate(year, 12, 31))) # set date to 31st December of given year


# Age group proportions
prop0to14 <- gather(prop0to14, countries$name, key="country", value="prop0to14") %>%
				 select(country, year, prop0to14)

prop65plus <- gather(prop65plus, countries$name, key="country", value="prop65plus") %>%
				 select(country, year, prop65plus)
				 
prop80plus <- gather(prop80plus, countries$name, key="country", value="prop80plus") %>%
				 select(country, year, prop80plus)
				 
popProportions <- left_join(prop0to14, prop65plus) %>%
					left_join(prop80plus)
				 
# Population projections
# First interpolate for single year predictions
# Timeframe for projections
years <- seq(2019, max(popProjections$year))

# Interpolate missing years in projections with splines
popProjections <- popProjections %>% mutate(year = length(years)) %>% 
								  group_by(country) %>% 
								  expand(year = years) %>%
								  left_join(popProjections) %>%
								  mutate(prop0to14_int = spline(x=year, y=prop0to14, xout=year)$y,
										 prop65plus_int = spline(x=year, y=prop65plus, xout=year)$y,
										 prop80plus_int = spline(x=year, y=prop80plus, xout=year)$y) %>%
								  ungroup()
								  
popProportions <- bind_rows(popProportions, popProjections) %>% 
					arrange(country,year) 
					
popProportions <- popProportions %>%
					mutate(prop0to14_int  = case_when(year <= 2018 ~ prop0to14,
														TRUE ~ prop0to14_int),
						   prop65plus_int = case_when(year <= 2018 ~ prop65plus,
														TRUE ~ prop65plus_int),
						   prop80plus_int = case_when(year <= 2018 ~ prop80plus,
														TRUE ~ prop80plus_int)) 

ggplot(data = popProportions, aes(x=year)) +
  geom_line(aes(y=prop0to14_int), col=ECDCcol[1])

#########################
## Plot incidence data ##
#########################

 ## Incidence by country
  ggplot(data = incData, aes(Year, Incidence)) +
  geom_bar(stat="identity", aes(fill = GroupType)) +
  labs(	x = "Year", y = "Incidence per 100,000") +
  facet_wrap(~ ReportingCountry) + #,  scales = "free_y") +
  theme(panel.grid.major = element_blank(), 
					  panel.grid.minor = element_blank(), 
					  panel.background = element_blank(),
					  axis.line = element_line(colour = "black"),
					  text = element_text(size=14),
					  axis.text.x = element_text(size = rel(0.75), angle = 60),
					  axis.text.x.top = element_text(vjust = -0.5)) +
					  guides(fill=guide_legend(title="Vaccine type")) +
  ggtitle(paste("Incidence of IPD in 50+ year olds")) +
  theme(plot.title = element_text(hjust = 0.5))
  