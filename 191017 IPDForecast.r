rm(list=ls())
dir <- "P:/Pneumo"
setwd(dir)
library(plyr)
library(tidyverse)
library(lubridate)
library(reshape2)
library(tsibble)
library(fable)
#library(feasts)
library(fabletools)
library(zoo)
library(EcdcColors)

myColours <- c("26 107 133", "241 214 118", "168 45 23")
ECDCcol <- sapply(strsplit(myColours, " "), function(x)
   rgb(x[1], x[2], x[3], maxColorValue=255))  # convert to hexadecimal

colours1 <- colorRampPalette(ECDCcol)(11)

## List countries
EEA <- c("AT","BE","BG", "HR", "CY", "CZ", "DK", "EE", "FI", "FR", "DE", "EL", "HU", "IS", "IE", "IT", "LV", "LI", "LT", "LU", "MT", "NL", "NO", "PL", "PT", "RO", "SK", "SI", "ES", "SE", "UK") 
countryList <- c("Austria", "Belgium", "Bulgaria", "Croatia", "Cyprus", "Czechia", "Denmark", "Estonia", "Finland", "France", "Germany", "Greece", "Hungary", "Iceland", "Ireland", "Italy", "Latvia", "Liechtenstein", "Lithuania","Luxembourg", "Malta", "Netherlands", "Norway", "Poland", "Portugal", "Romania", "Slovakia", "Slovenia", "Spain", "Sweden", "UK" )
countries <- tibble(GeoCode = EEA, country = countryList)
countriesScot <- add_row(countries, GeoCode = "SC", country = "Scotland")

########################		
## IPD incidence data ##
########################

PCVdates <- read_csv("VaccineDates.csv", col_names=TRUE)   # Read in dates of PCV introduction

IPDData <- read_csv("190412IPD.csv", col_names=TRUE) %>%
		   mutate(
				  ## Age groups. N.B. These must correspond with population groups 
				  age = as.integer(Age),
				  ageGroup = case_when(age < 5 ~ 1,		# 0-4
							  age >= 5  & age < 15 ~ 2,	# 5-14
							  age >= 15 & age < 50 ~ 3,	# 15-49
							  age >= 50 & age < 65 ~ 4,	# 50-64
							  age >= 65 & age < 75 ~ 5,	# 65-74
							  age >= 75 & age < 85 ~ 6,	# 75-84
							  age >= 85 ~ 7) %>% # 85+
							  factor(labels = c("0-4", "5-14", "15-49", "50-64", "65-74", "75-84", "85+")),
							  
				  ## Countries
				  GeoCode = case_when(GeoCode == "UKM" ~ "SC", TRUE ~ as.character(ReportingCountry)), # Recode region 'UKM' as Scotland 
				 
				 ## Serotypes
				  serotype = recode(Serotype, "NT" = "UNK",           		 
											  "NTYP" = "UNK",		   
											  "O" = "OTHER") %>%
									as_factor(),
											  
				  groupType = case_when(serotype == "4"  ~ "PCV-7",    # Group serotypes according to vaccine 
										serotype == "6B" ~ "PCV-7",
										serotype == "9V" ~ "PCV-7",
										serotype == "14"  ~ "PCV-7",
										serotype == "18C" ~ "PCV-7",
										serotype ==  "19F" ~ "PCV-7",
										serotype == "23F" ~ "PCV-7",
										serotype == "1" ~ "PCV-10",
										serotype == "5" ~ "PCV-10",
										serotype == "7F" ~ "PCV-10",
										serotype == "3" ~ "PCV-13",
										serotype == "6A" ~ "PCV-13",
										serotype == "19A" ~ "PCV-13",
										serotype == "2"  ~ "PPV-23",
										serotype == "8" ~ "PPV-23",
										serotype == "9N" ~ "PPV-23",
										serotype == "10A" ~ "PPV-23",
										serotype == "11A" ~ "PPV-23",
										serotype == "12F" ~ "PPV-23",
										serotype == "15B"  ~ "PPV-23",
										serotype == "17F" ~ "PPV-23",
										serotype == "20" ~ "PPV-23",
										serotype == "22F" ~ "PPV-23",
										serotype == "OTHER" ~ "Other") %>%
										as_factor() %>%
										fct_explicit_na(na_level = "Not serotyped"), 
										
					## Specimen type				   
					specimen = recode(Specimen, "JOINT"  = "OTHER",    # Group samples from joint, peritoneal and pleural fluid
                                     "PERIT" = "OTHER",
									 "PLEURAL" = "OTHER",
									 "O"  = "OTHER",
                                     "NULL" = "UNK"),
					## Dates
					date = as.Date(DateUsedForStatisticsISO, "%Y-%m-%d"),  # Convert to date format
					date = case_when(nchar(DateUsedForStatisticsISO)==7 ~ paste(DateUsedForStatisticsISO,"-01",sep=""),    # Where no day specified, set to 1st.
									 nchar(DateUsedForStatisticsISO)==4 ~ paste(DateUsedForStatisticsISO,"-01-01",sep=""), # Where no day or month specified, set to 1st of the 1st.
									 TRUE ~ as.character(date)),
				    date = ymd(date),
				    
					year = substr(date,1,4) %>% as.integer()) %>%    # Define year for matching 
				   
					select(-DateUsedForStatisticsISO) %>%  # Remove original date variable

				    left_join(countries) # country names
			
firstYear <- IPDData %>% select(year) %>% min()

#####################
## Population data ##
#####################

tessyPop <- read_csv("TESSy population.csv", col_names=TRUE) %>% # Read in TESSy population data
		    filter(ReportingCountry %in% EEA & ReportYear >= firstYear) %>%
			mutate(ageGroup = case_when(AgeGroupId == 1   | AgeGroupId == 2  ~ 1,      # 0-4
										AgeGroupId >= 3   & AgeGroupId <= 5  ~ 2,      # 5-14
										AgeGroupId >= 6   & AgeGroupId <= 12 ~ 3,      # 15-49
										AgeGroupId >= 13  & AgeGroupId <= 15 ~ 4,      # 50-64
										AgeGroupId >= 16  & AgeGroupId <= 17 ~ 5,  	   # 65-74
										AgeGroupId >= 18  & AgeGroupId <= 19 ~ 6,      # 75-84
										AgeGroupId >= 20  & AgeGroupId <= 23 ~ 7) %>% # 85+
										factor(labels = c("0-4", "5-14", "15-49", "50-64", "65-74", "75-84", "85+"))) %>%
			filter(!is.na(ageGroup))  %>% # Original population database included other ways of grouping ages. Drop these values	
			rename(year = ReportYear, GeoCode = ReportingCountry) %>% # Rename to match with IPDData
			left_join(countries) %>% # Add in country names
			group_by(year, country, ageGroup) %>% # Group population by age groups 
			summarise(population = sum(Population)) %>%
			ungroup() %>%
			as_tsibble(index = year, key = c(country, ageGroup)) 
	
eurostatPop <- read_csv("EUROSTATpop.csv", col_names=TRUE)   # Read in EUROSTAT population data
scotPop <- read_csv("ScotPop.csv", col_names=TRUE)   # Source: https://www.nrscotland.gov.uk/statistics-and-data/statistics/statistics-by-theme/population/population-estimates/mid-year-population-estimates

###############
## Incidence ##
###############
years1 <- seq(2002, 2017) 

IPDData <- IPDData %>% filter(Classification == "CONF" & SubsetEpi == 1) %>% # Include only confirmed cases and select subset for epi. to avoid double-counting lab.results
						left_join(tessyPop) %>%
						arrange(country, year, ageGroup) 
		
incData <- IPDData %>% 	filter(age >= 50)  %>%
						group_by(country, year, ageGroup, groupType, population) %>%
						summarise(total = n()) %>% # Number of cases per year by age group, type and country
						mutate(incidence = (total/population)*100000) %>% # Incidence per 100,000	
						ungroup() 
						
incData <- incData %>% 	mutate(year = length(years1)) %>% 
						group_by(country, ageGroup, groupType) %>% 
						expand(year = years1) %>%
						left_join(incData) %>%	
						ungroup() %>%
						select(-population) %>%
						as_tsibble(index = year, key = c(country, ageGroup, groupType))
					
countryNo <- 9

## Incidence by country
  ggplot(data = incData, aes(year, incidence)) +
  geom_line(lwd = 0.7) +
  labs(x = "Year", y = "Incidence per 100,000") +
  facet_grid(groupType ~ ageGroup, labeller = label_wrap_gen(width = 2, multi_line = TRUE) ,  scales = "free_y") + 
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
  ggtitle(paste("Rate of IPD incidence in", countries$country[countryNo])) +
  theme(plot.title = element_text(hjust = 0.5))
  
############################
## Proportion of serotype ##
############################

# Proportion of outcome by serotype
propTypes <- as_tibble(incData) %>% filter(groupType != "Not serotyped") %>% # 
			   group_by(country, year, ageGroup, groupType) %>%
               dplyr::summarise(n = sum(total)) %>%
			   mutate(freq = n / sum(n))
			   

## Trends in proportion of types, by country
ggplot(data = propTypes, aes(year, freq)) +
  geom_line(lwd = 0.7) +
  labs(x = "Year", y = "Proportion of typed samples") +
  facet_grid(groupType ~ ageGroup, labeller = label_wrap_gen(width = 2, multi_line = TRUE)) + # ,  scales = "free_y") + 
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

## Population proportions by age group. Source: EUROSTAT https://ec.europa.eu/eurostat/databrowser/view/tps00010/default/table?lang=en
prop0to14 <- read_csv("prop0to14.csv", col_names=TRUE)   # Proportion of population aged 0 to 14
prop65plus <- read_csv("prop65plus.csv", col_names=TRUE) # Proportion of population aged over 65
prop80plus <- read_csv("prop80plus.csv", col_names=TRUE) # Proportion of population aged over 80

prop0to14 <- gather(prop0to14, countries$name, key="country", value="prop0to14") %>%
				 select(country, year, prop0to14)

prop65plus <- gather(prop65plus, countries$name, key="country", value="prop65plus") %>%
				 select(country, year, prop65plus)
				 
prop80plus <- gather(prop80plus, countries$name, key="country", value="prop80plus") %>%
				 select(country, year, prop80plus)

# Historical population proportions by age group				 
histPopProportions <- left_join(prop0to14, prop65plus) %>%
					left_join(prop80plus)
				 
# Population projections
# First interpolate for single year predictions
# Timeframe for projections
popProjections <- read_csv("popProjections.csv", col_names=TRUE) # Source: UN DESA World Population Prospects 2019 https://population.un.org/wpp/Download/Probabilistic/Population/
years <- seq(2019, max(popProjections$year)) # Extrapolate from 2019 (current year)

# Interpolate missing years in projections with splines
popProjections <- popProjections %>% mutate(year = length(years)) %>% 
								  group_by(country) %>% 
								  expand(year = years) %>%
								  left_join(popProjections) %>%
								  mutate(prop0to14_L95_int = spline(x=year, y=prop0to14_L95, xout=year)$y,
										 prop0to14_int = spline(x=year, y=prop0to14, xout=year)$y,
										 prop0to14_U95_int = spline(x=year, y=prop0to14_U95, xout=year)$y,
										 prop65plus_L95_int = spline(x=year, y=prop65plus_L95, xout=year)$y,
										 prop65plus_int = spline(x=year, y=prop65plus, xout=year)$y,
										 prop65plus_U95_int = spline(x=year, y=prop65plus_U95, xout=year)$y,
										 prop80plus_L95_int = spline(x=year, y=prop80plus_L95, xout=year)$y,
										 prop80plus_int = spline(x=year, y=prop80plus, xout=year)$y,
										 prop80plus_U95_int = spline(x=year, y=prop80plus_U95, xout=year)$y) %>%
								  ungroup()
								  
popProportions <- bind_rows(histPopProportions, popProjections) %>% 
					arrange(country,year) 
					
popProportions <- popProportions %>%
					mutate(prop0to14_int  = case_when(year <= 2018 ~ prop0to14,
														TRUE ~ prop0to14_int),
						   prop65plus_int = case_when(year <= 2018 ~ prop65plus,
														TRUE ~ prop65plus_int),
						   prop80plus_int = case_when(year <= 2018 ~ prop80plus,
														TRUE ~ prop80plus_int)) %>%
					filter(country != "Liechtenstein")

# Plot proportion aged 0-14
ggplot(data = popProportions %>% filter(year <= 2040), aes(x=year)) +
  facet_wrap(~ country) + #, labeller = label_wrap_gen(width = 2, multi_line = TRUE)) + # ,  scales = "free_y") + 
  geom_line(aes(y=prop0to14_int), col=ECDCcol[1]) +
  geom_ribbon(aes(x=year,ymin=prop0to14_L95_int, ymax=prop0to14_U95_int), fill=ECDCcol[1], alpha=.5) +
   ylim(0,25)+
   theme(panel.grid.major = element_blank(), 
		panel.grid.minor = element_blank(), 
		panel.background = element_blank(),
		axis.line = element_line(colour = "black"),
		text = element_text(size=14),
		axis.text.x.top = element_text(vjust = -0.5)) +
		scale_x_continuous(breaks = seq(2005, 2040, by = 5)) +
		guides(fill=guide_legend(title="")) +
	labs(title = "Proportion of population aged 0 to 14", y = "", x = "")

# Plot proportion aged 65 plus
ggplot(data = popProportions %>% filter(year <= 2040), aes(x=year)) +
  facet_wrap(~ country) + #, labeller = label_wrap_gen(width = 2, multi_line = TRUE)) + # ,  scales = "free_y") + 
  geom_line(aes(y=prop65plus_int), col=ECDCcol[1]) +
  geom_ribbon(aes(x=year,ymin=prop65plus_L95_int, ymax=prop65plus_U95_int), fill=ECDCcol[1], alpha=.5) +
   ylim(0,40)+
   theme(panel.grid.major = element_blank(), 
		panel.grid.minor = element_blank(), 
		panel.background = element_blank(),
		axis.line = element_line(colour = "black"),
		text = element_text(size=14),
		axis.text.x.top = element_text(vjust = -0.5)) +
		scale_x_continuous(breaks = seq(2005, 2040, by = 5)) +
		guides(fill=guide_legend(title="")) +
	labs(title = "Proportion of population aged over 65", y = "", x = "")


# Plot proportion aged 80 plus
ggplot(data = popProportions %>% filter(year <= 2040), aes(x=year)) +
  facet_wrap(~ country) + #, labeller = label_wrap_gen(width = 2, multi_line = TRUE)) + # ,  scales = "free_y") + 
  geom_line(aes(y=prop80plus_int), col=ECDCcol[1]) +
  geom_ribbon(aes(x=year,ymin=prop80plus_L95_int, ymax=prop80plus_U95_int), fill=ECDCcol[1], alpha=.5) +
   ylim(0,15)+
   theme(panel.grid.major = element_blank(), 
		panel.grid.minor = element_blank(), 
		panel.background = element_blank(),
		axis.line = element_line(colour = "black"),
		text = element_text(size=14),
		axis.text.x.top = element_text(vjust = -0.5)) +
		scale_x_continuous(breaks = seq(2005, 2040, by = 5)) +
		guides(fill=guide_legend(title="")) +
	labs(title = "Proportion of population aged over 80", y = "", x = "")

# Proportion of children under three not in formal childcare
childcare0to2 <- read_csv("childcare0to2.csv", col_names=TRUE) %>% # Source: EUROSTAT https://ec.europa.eu/eurostat/databrowser/view/tps00185/default/table?lang=en
					gather(countries$country, key="country", value="childcare0to2") %>%
					select(country, year, childcare0to2) %>%
					filter(country != "Liechtenstein")
					#mutate(date = as.Date(ISOdate(year, 12, 31))) # set date to 31st December of given year

# Forecast childcare
childcare_ts <- as_tsibble(childcare0to2, index=year, key=country) %>%# define tsibble with country as key and year as index
			   rename(proportion = childcare0to2) %>%
			   filter(proportion != "NA") 
			   
childcare_ts_short <- childcare_ts #%>% filter(year < 2016) # Remove later data points for out-of-sample validation. Forecast for the number of years that have bene excluded and then use accuracy() to assess.
			   
childcare_fit <- childcare_ts_short %>% model(ets = ETS(box_cox(proportion, lambda = 0) ~ error("A") + trend(("Ad"), phi_range = c(0.88, 0.98))))# Fit damped model with additive error term

childcare_fitted <- childcare_fit %>%
						fitted() %>%
						as_tibble() %>%
						rename(childcare_fit = .fitted) %>%
						select(country, year, childcare_fit)				

childcare_fcst <- childcare_fit %>%	forecast(h = "23 years") %>%
					mutate(interval = hilo(.distribution, 95))

childcare_extend <- as_tibble(childcare_fcst) %>% 
					select(country, year, proportion, interval) %>%
					unnest(interval) %>%
					select(- .level) %>%
					mutate(childcare_L95 = pmax(0, .lower), childcare_U95 = pmin(100, .upper)) %>%
					rename(childcare_fit = proportion)
					 
childcare_fitted <- bind_rows(childcare_fitted, childcare_extend) %>% # All model values (fit and forecast)
						arrange(country, year)
						
childcare0to2 <- right_join(childcare0to2, childcare_fitted) 

# Diagnostics
childcare_fcst %>% accuracy(childcare_ts) %>% 
					print(n=30)

childcare_fcst %>% filter(country == "Sweden") %>% 
					autoplot(childcare_ts)

childcare_fit %>% filter(country == "Sweden") %>% 
					report()
							
ggplot(data = childcare0to2, aes(x=year)) +
  facet_wrap(~ country) + #, labeller = label_wrap_gen(width = 2, multi_line = TRUE)) + # ,  scales = "free_y") + 
  geom_point(aes(y=childcare0to2), col=ECDCcol[1]) +
  geom_line(aes(y=childcare_fit), col=ECDCcol[1]) +
  geom_ribbon(aes(x=year,ymin=childcare_L95, ymax=childcare_U95), fill=ECDCcol[1], alpha=.5) +
  geom_line(aes(y=childcare_fit), col="black") +
    ylim(0,100)+
   theme(panel.grid.major = element_blank(), 
		panel.grid.minor = element_blank(), 
		panel.background = element_blank(),
		axis.line = element_line(colour = "black"),
		text = element_text(size=14),
		axis.text.x.top = element_text(vjust = -0.5)) +
		scale_x_continuous(breaks = seq(2005, 2040, by = 10)) +
  		guides(fill=guide_legend(title="")) +
	labs(title = "Proportion of children under three not enrolled in formal childcare", y = "")

					
# PCV-7 uptake 
PCV7Data <- read_csv("PCV7coverage.csv", col_names=TRUE) %>% # Read in PCV7 uptake data
				melt() %>%
				as_tibble() %>%
				rename(year = variable, PCV7cov = value) %>%
				mutate(year = as.numeric(levels(year))[year], country = as_factor(country)) %>%
				arrange(country) 
				
extend <- seq(min(PCV7Data$year), 2040)	
		
# Forecast PCV-7 uptake				
PCV7uptake <- PCV7Data %>%
					mutate(year = length(extend)) %>% 
					group_by(country) %>% 
					expand(year = extend) %>%
					left_join(PCV7Data) %>% 
					mutate(PCV7_fit = case_when(year >=2018 ~ 0, TRUE ~ PCV7cov),
						   PCV7_L95 = case_when(year >=2018 ~ 0, TRUE ~ (PCV7cov - 5)), # Define 95% CI as +/- 5% points
						   PCV7_U95 = case_when(year >=2018 ~ 0, TRUE ~ pmin((PCV7cov + 5),100))) 
						   
# PCV-10 uptake data
PCV10Data <- read_csv("PCV10coverage.csv", col_names=TRUE) %>% # Read in PCV10 uptake data
				melt() %>%
				as_tibble() %>%
				rename(year = variable, PCV10cov = value) %>%
				mutate(year = as.numeric(levels(year))[year], country = as_factor(country)) %>%
				arrange(country) %>%
				filter(PCV10cov != "NA")


# Forecast PCV-10 uptake
PCV10uptake <- PCV10Data %>%
					mutate(year = length(extend)) %>% 
					group_by(country) %>% 
					expand(year = extend) %>%
					left_join(PCV10Data) %>% 
					mutate(PCV10_fit = rollmean(PCV10cov, 2, align = "right", fill = NA) %>%  # Average of previous 2 years' uptake
									  lag(1)) %>%
					fill(PCV10_fit, .direction = "down") %>%
					mutate(PCV10_fit = case_when(year <= 2017 ~ PCV10cov, TRUE ~ PCV10_fit),
						   PCV10_L95 = case_when(year <= 2017 ~ NA_real_, TRUE ~ (PCV10_fit - 5)), # Define 95% CI as +/- 5% points
						   PCV10_U95 = case_when(year <= 2017 ~ NA_real_, TRUE ~ pmin((PCV10_fit + 5),100))) 

# N.B. Need manual intervention for countries with less than two years equilibrium uptake data: Poland
PCV10uptake <- PCV10uptake %>% 
				mutate(PCV10_fit  = case_when(country == "Poland" & year >= 2018 ~ 94, TRUE ~ PCV10_fit),
						PCV10_L95 = case_when(country == "Poland" & year >= 2018 ~ 94 - 5, TRUE ~ PCV10_L95), # Define 95% CI as +/- 5% points
						PCV10_U95 = case_when(country == "Poland" & year >= 2018 ~ pmin(94 + 5, 100), TRUE ~ PCV10_U95)) %>%
				mutate(PCV10_fit  = case_when(country == "Belgium" & year >= 2018 ~ 94, TRUE ~ PCV10_fit),
						PCV10_L95 = case_when(country == "Belgium" & year >= 2018 ~ 94 - 5, TRUE ~ PCV10_L95), # Define 95% CI as +/- 5% points
						PCV10_U95 = case_when(country == "Belgium" & year >= 2018 ~ pmin(94 + 5, 100), TRUE ~ PCV10_U95))
		
# PCV-13 uptake data
PCV13Data <- read_csv("PCV13coverage.csv", col_names=TRUE) %>% # Read in PCV13 uptake data
				melt() %>%
				as_tibble() %>%
				rename(year = variable, PCV13cov = value) %>%
				mutate(year = as.numeric(levels(year))[year], country = as_factor(country)) %>%
				arrange(country) 

# Forecast PCV-13 uptake
PCV13uptake <- PCV13Data %>%
					mutate(year = length(extend)) %>% 
					group_by(country) %>% 
					expand(year = extend) %>%
					left_join(PCV13Data) %>% 
					mutate(PCV13_fit = rollmean(PCV13cov, 2, align = "right", fill = NA) %>%  # Average of previous 2 years uptake
									  lag(1)) %>%
					fill(PCV13_fit, .direction = "down") %>%
					mutate(PCV13_fit = case_when(year <= 2017 ~ PCV13cov, TRUE ~ PCV13_fit),
						   PCV13_L95 = case_when(year <= 2017 ~ NA_real_, TRUE ~ (PCV13_fit - 5)), # Define 95% CI as +/- 5% points
						   PCV13_U95 = case_when(year <= 2017 ~ NA_real_, TRUE ~ pmin((PCV13_fit + 5),100))) 

# N.B. Need manual intervention for countries with less than two years equilibrium uptake data: Romania, Portugal, Belgium, Spain
PCV13uptake <- PCV13uptake %>% 
				mutate(PCV13_fit  = case_when(country == "Romania" & year >= 2018 ~ 75, TRUE ~ PCV13_fit),
						PCV13_L95 = case_when(country == "Romania" & year >= 2018 ~ 75 - 5, TRUE ~ PCV13_L95), # Define 95% CI as +/- 5% points
						PCV13_U95 = case_when(country == "Romania" & year >= 2018 ~ pmin(75 + 5, 100), TRUE ~ PCV13_U95)) %>%
				mutate(PCV13_fit  = case_when(country == "Portugal" & year >= 2018 ~ 97, TRUE ~ PCV13_fit),
						PCV13_L95 = case_when(country == "Portugal" & year >= 2018 ~ 97 - 5, TRUE ~ PCV13_L95), # Define 95% CI as +/- 5% points
						PCV13_U95 = case_when(country == "Portugal" & year >= 2018 ~ pmin(97 + 5, 100), TRUE ~ PCV13_U95)) %>%
				mutate(PCV13_fit  = case_when(country == "Belgium" & (year == 2017 | year == 2018) ~ 0, TRUE ~ PCV13_fit),
						PCV13_L95 = case_when(country == "Belgium" & year >= 2018 ~ 0, TRUE ~ PCV13_L95), # Define 95% CI as +/- 5% points
						PCV13_U95 = case_when(country == "Belgium" & year >= 2018 ~ 5, TRUE ~ PCV13_U95)) %>%
				mutate(PCV13_fit  = case_when(country == "Spain" & year >= 2018 ~ 48, TRUE ~ PCV13_fit),
						PCV13_L95 = case_when(country == "Spain" & year >= 2018 ~ 48 - 5, TRUE ~ PCV13_L95), # Define 95% CI as +/- 5% points
						PCV13_U95 = case_when(country == "Spain" & year >= 2018 ~ pmin(48 + 5, 100), TRUE ~ PCV13_U95))

PCVuptake <- full_join(PCV7uptake, PCV10uptake, by=c("country", "year")) %>%
			 full_join(PCV13uptake, by=c("country", "year"))

# Belgium is switching back to PCV-13 from PCV-10 (2020)
PCVuptake <- PCVuptake %>% mutate(PCV13_fit = case_when(country == "Belgium" & year >= 2020 ~ PCV10_fit, TRUE ~ PCV13_fit),
								  PCV13_L95 = case_when(country == "Belgium" & year >= 2020 ~ PCV10_L95, TRUE ~ PCV13_L95), # Define 95% CI as +/- 5% points
								  PCV13_U95 = case_when(country == "Belgium" & year >= 2020 ~ PCV10_U95, TRUE ~ PCV13_U95)) %>%
							mutate(PCV10_fit = case_when(country == "Belgium" & year >= 2020 ~ 0, TRUE ~ PCV10_fit),
								   PCV10_L95 = case_when(country == "Belgium" & year >= 2020 ~ 0, TRUE ~ PCV10_L95), # Define 95% CI as +/- 5% points
								   PCV10_U95 = case_when(country == "Belgium" & year >= 2020 ~ 0, TRUE ~ PCV10_U95)) 
					  
# Plot projections of uptake of PCV vaccination					
ggplot(data = PCVuptake, aes(x=year)) +
  facet_wrap(~ country) + #, labeller = label_wrap_gen(width = 2, multi_line = TRUE)) + # ,  scales = "free_y") + 
  #geom_line(aes(y=PCV10_fit), col=ECDCcol[1]) +
  geom_point(aes(y=PCV10cov), col=ECDCcol[1]) +
  geom_ribbon(aes(x=year,ymin=PCV10_L95, ymax=PCV10_U95), fill=ECDCcol[1], alpha=.5) +
  geom_line(aes(y=PCV10_fit), col="black") +
    ylim(0,100)+
   theme(panel.grid.major = element_blank(), 
		panel.grid.minor = element_blank(), 
		panel.background = element_blank(),
		axis.line = element_line(colour = "black"),
		text = element_text(size=14),
		axis.text.x.top = element_text(vjust = -0.5)) +
		scale_x_continuous(breaks = seq(2005, 2040, by = 10)) +
  		guides(fill=guide_legend(title="")) +
	labs(title = "Uptake of PCV-10 vaccination in children", y = "")

												
ggplot(data = PCVuptake, aes(x=year)) +
  facet_wrap(~ country) + #, labeller = label_wrap_gen(width = 2, multi_line = TRUE)) + # ,  scales = "free_y") + 
  geom_point(aes(y=PCV13cov), col=ECDCcol[1]) +
  geom_line(aes(y=PCV13_fit), col=ECDCcol[1]) +
  geom_ribbon(aes(x=year,ymin=PCV13_L95, ymax=PCV13_U95), fill=ECDCcol[1], alpha=.5) +
  geom_line(aes(y=PCV13_fit), col="black") +
    ylim(0,100)+
   theme(panel.grid.major = element_blank(), 
		panel.grid.minor = element_blank(), 
		panel.background = element_blank(),
		axis.line = element_line(colour = "black"),
		text = element_text(size=14),
		axis.text.x.top = element_text(vjust = -0.5)) +
		scale_x_continuous(breaks = seq(2005, 2040, by = 10)) +
  		guides(fill=guide_legend(title="")) +
	labs(title = "Uptake of PCV-13 vaccination in children", y = "")

		
ggplot(data = PCV7Data, aes(x=year)) +
	facet_wrap(~country) +
	geom_line(aes(y=PCV7cov), col=ECDCcol[1]) +
	#ylim(0,50)+
	scale_color_manual(values = colours1, guide = FALSE) +
	scale_fill_manual(values = colours1, guide = FALSE) +
	theme(panel.grid.major = element_blank(), 
		panel.grid.minor = element_blank(), 
		panel.background = element_blank(),
		axis.line = element_line(colour = "black"),
		text = element_text(size=14),
		axis.text.x.top = element_text(vjust = -0.5)) +
		guides(fill=guide_legend(title=""))+
	labs(title = "PCV7 coverage in children", y = "")

										
					
# PPV coverage
PPV23Data <- read_csv("PPV23coverage.csv", col_names=TRUE) # Read in PPV23 coverage data

# Influenza vaccination coverage
# Sources: https://www.ecdc.europa.eu/en/publications-data/seasonal-influenza-vaccination-europe-vaccination-recommendations-and-coverage-2007-2015 
# https://www.ecdc.europa.eu/sites/portal/files/documents/Seasonal-influenza-antiviral-use-EU-EEA-Member-States-December-2018_0.pdf 
# N.B. UK is population weighted average of four regions
fluCov55A <- read_csv("fluCov55A.csv", col_names=TRUE) %>%
			gather(countries$country, key="country", value="fluCov55A") %>%
			select(country, year, fluCov55A) %>%	# influenza vaccination coverage in over 55s, administrative sample
			mutate(country = as_factor(country))
			
fluCov59A <- read_csv("fluCov59A.csv", col_names=TRUE) %>%
			gather(countries$country, key="country", value="fluCov59A") %>%
			select(country, year, fluCov59A) %>%	# influenza vaccination coverage in over 59s, administrative sample
			mutate(country = as_factor(country))
			
fluCov60A <- read_csv("fluCov60A.csv", col_names=TRUE) %>%
			gather(countries$country, key="country", value="fluCov60A") %>%
			select(country, year, fluCov60A) %>%	# influenza vaccination coverage in over 60s, administrative sample
			mutate(country = as_factor(country))
						
fluCov60S <- read_csv("fluCov60S.csv", col_names=TRUE)  %>%
			gather(countries$country, key="country", value="fluCov60S") %>%
			select(country, year, fluCov60S) %>%	# influenza vaccination coverage in over 60s, survey sample
			mutate(country = as_factor(country))
			
fluCov65A <- read_csv("fluCov65A.csv", col_names=TRUE) %>%
			gather(countries$country, key="country", value="fluCov65A") %>%
			select(country, year, fluCov65A) %>%	# influenza vaccination coverage in over 65s, administrative sample
			mutate(country = as_factor(country))
			
fluCov65S <- read_csv("fluCov65S.csv", col_names=TRUE) %>%
			gather(countries$country, key="country", value="fluCov65S") %>%
			select(country, year, fluCov65S) %>%	# influenza vaccination coverage in over 65s, survey sample
			mutate(country = as_factor(country))

fluCoverage <- left_join(fluCov55A, fluCov59A, by = c("country", "year")) %>%
			   left_join(fluCov60A, by = c("country", "year")) %>%
			   left_join(fluCov60S, by = c("country", "year")) %>%
			   left_join(fluCov65A, by = c("country", "year")) %>%
			   left_join(fluCov65S, by = c("country", "year")) %>%
			   filter(country != "Liechtenstein") %>%
			   rowwise() %>% 
			   mutate(sumCov59 = case_when(is.na(fluCov55A) & is.na(fluCov59A) ~ NA_real_, TRUE ~ sum(fluCov55A,fluCov59A, na.rm=TRUE))) %>%
			   mutate(sumCov60 = case_when(is.na(fluCov55A) & is.na(fluCov59A) & is.na(fluCov59A) ~ NA_real_, TRUE ~ sum(fluCov55A,fluCov59A,fluCov60A, na.rm=TRUE))) %>%
     		   mutate(sumCov65 = case_when(is.na(fluCov55A) & is.na(fluCov59A) & is.na(fluCov60A) & is.na(fluCov65A) ~ NA_real_, TRUE ~ sum(fluCov55A,fluCov59A,fluCov60A,fluCov65A, na.rm=TRUE))) 
			   
			
# Plot influenza vaccination coverage data	
ggplot(data = fluCoverage, aes(x=year)) +
	facet_wrap(~country) +
	geom_line(aes(y=fluCov55A), col=ECDCcol[1]) +
	ylim(0,100)+
	scale_color_manual(values = colours1, guide = FALSE) +
	scale_fill_manual(values = colours1, guide = FALSE) +
	theme(panel.grid.major = element_blank(), 
		panel.grid.minor = element_blank(), 
		panel.background = element_blank(),
		axis.line = element_line(colour = "black"),
		text = element_text(size=12),
		axis.text.x.top = element_text(vjust = -0.5)) +
		scale_x_continuous(breaks = seq(2006, 2018, by = 2)) +
		guides(fill=guide_legend(title=""))+
	labs(title = "Seasonal influenza vaccination coverage in 55-58 year olds", y = "")

ggplot(data = fluCoverage, aes(x=year)) +
	facet_wrap(~country) +
	geom_line(aes(y=sumCov59), col=ECDCcol[1]) +
	ylim(0,100)+
	scale_color_manual(values = colours1, guide = FALSE) +
	scale_fill_manual(values = colours1, guide = FALSE) +
	theme(panel.grid.major = element_blank(), 
		panel.grid.minor = element_blank(), 
		panel.background = element_blank(),
		axis.line = element_line(colour = "black"),
		text = element_text(size=12),
		axis.text.x.top = element_text(vjust = -0.5)) +
		scale_x_continuous(breaks = seq(2006, 2018, by = 2)) +
		guides(fill=guide_legend(title=""))+
	labs(title = "Seasonal influenza vaccination coverage in 59 year olds", y = "")

ggplot(data = fluCoverage, aes(x=year)) +
	facet_wrap(~country) +
	geom_line(aes(y=sumCov60), col=ECDCcol[1]) +
	ylim(0,100)+
	scale_color_manual(values = colours1, guide = FALSE) +
	scale_fill_manual(values = colours1, guide = FALSE) +
	theme(panel.grid.major = element_blank(), 
		panel.grid.minor = element_blank(), 
		panel.background = element_blank(),
		axis.line = element_line(colour = "black"),
		text = element_text(size=12),
		axis.text.x.top = element_text(vjust = -0.5)) +
		scale_x_continuous(breaks = seq(2006, 2018, by = 2)) +
		guides(fill=guide_legend(title=""))+
	labs(title = "Seasonal influenza vaccination coverage in 60-64 year olds", y = "")

ggplot(data = fluCoverage, aes(x=year)) +
	facet_wrap(~country) +
	geom_line(aes(y=sumCov65), col=ECDCcol[1]) +
	ylim(0,100)+
	scale_color_manual(values = colours1, guide = FALSE) +
	scale_fill_manual(values = colours1, guide = FALSE) +
	theme(panel.grid.major = element_blank(), 
		panel.grid.minor = element_blank(), 
		panel.background = element_blank(),
		axis.line = element_line(colour = "black"),
		text = element_text(size=12),
		axis.text.x.top = element_text(vjust = -0.5)) +
		scale_x_continuous(breaks = seq(2006, 2018, by = 2)) +
		guides(fill=guide_legend(title=""))+
	labs(title = "Seasonal influenza vaccination coverage in the over 65s", y = "")

# Summarise coverage of seasonal influenza vaccination by age group
summFluCov <- fluCoverage %>% select(country, year, fluCov55A, sumCov59, sumCov60, sumCov65) %>%
								rename(fluCov55to58 = fluCov55A,
									   fluCov59 = sumCov59,
									   fluCov60to64 = sumCov60,
									   fluCov65plus = sumCov65) 

# Interpolate missing values. This is largely manual and very messy!
library(forecast)									   
summFluCov$fluCov65plus_int <- summFluCov %>% select(fluCov65plus) %>%
											  na.interp()
detach("package:forecast", unload=TRUE)

# Manually curate imputed values
summFluCov <- summFluCov %>%  mutate(fluCov65plus_int = case_when(country == "Austria" ~ NA_real_, TRUE ~ fluCov65plus_int)) %>%
							  mutate(fluCov65plus_int = case_when(country == "Belgium" ~ NA_real_, TRUE ~ fluCov65plus_int)) %>%
							  mutate(fluCov65plus_int = case_when(country == "Bulgaria" ~ NA_real_, TRUE ~ fluCov65plus_int)) %>%
							  mutate(fluCov65plus_int = case_when(country == "Croatia" & (year <= 2012 | year > 2015) ~ NA_real_, TRUE ~ fluCov65plus_int)) %>%
							  mutate(fluCov65plus_int = case_when(country == "Cyprus" ~ NA_real_, TRUE ~ fluCov65plus_int)) %>%
							  mutate(fluCov65plus_int = case_when(country == "Czechia" & (year <= 2014 | year >= 2016) ~ NA_real_, TRUE ~ fluCov65plus_int)) %>%
							  mutate(fluCov65plus_int = case_when(country == "Denmark" & year >= 2018 ~ NA_real_, TRUE ~ fluCov65plus_int)) %>%
							  mutate(fluCov65plus_int = case_when(country == "Estonia" & year >= 2018 ~ NA_real_, TRUE ~ fluCov65plus_int)) %>%
							  mutate(fluCov65plus_int = case_when(country == "Finland" & year >= 2018 ~ NA_real_, TRUE ~ fluCov65plus_int)) %>%
							  mutate(fluCov65plus_int = case_when(country == "France" & year >= 2016 ~ NA_real_, TRUE ~ fluCov65plus_int)) %>%
							  mutate(fluCov65plus_int = case_when(country == "Germany" & (year <= 2013 | year >= 2018) ~ NA_real_, TRUE ~ fluCov65plus_int)) %>%
							  mutate(fluCov65plus_int = case_when(country == "Greece" ~ NA_real_, TRUE ~ fluCov65plus_int)) %>%
							  mutate(fluCov65plus_int = case_when(country == "Hungary" & year <= 2008 ~ NA_real_, TRUE ~ fluCov65plus_int)) %>%
							  mutate(fluCov65plus_int = case_when(country == "Italy" & year >= 2018 ~ NA_real_, TRUE ~ fluCov65plus_int)) %>%
							  mutate(fluCov65plus_int = case_when(country == "Latvia" & year <= 2008 ~ NA_real_, TRUE ~ fluCov65plus_int)) %>%
							  mutate(fluCov65plus_int = case_when(country == "Lithuania" & year >= 2018 ~ NA_real_, TRUE ~ fluCov65plus_int)) %>%
							  mutate(fluCov65plus_int = case_when(country == "Luxembourg" & (year <= 2008 | year > 2014) ~ NA_real_, TRUE ~ fluCov65plus_int)) %>%
							  mutate(fluCov65plus_int = case_when(country == "Malta" & year >= 2016 ~ NA_real_, TRUE ~ fluCov65plus_int)) %>%
							  mutate(fluCov65plus_int = case_when(country == "Netherlands" & year >= 2018 ~ NA_real_, TRUE ~ fluCov65plus_int)) %>%
							  mutate(fluCov65plus_int = case_when(country == "Norway" & year >= 2018 ~ NA_real_, TRUE ~ fluCov65plus_int)) %>%
							  mutate(fluCov65plus_int = case_when(country == "Poland" & year >= 2018 ~ NA_real_, TRUE ~ fluCov65plus_int)) %>%
							  mutate(fluCov65plus_int = case_when(country == "Portugal" ~ NA_real_, TRUE ~ fluCov65plus_int)) %>%
							  mutate(fluCov65plus_int = case_when(country == "Romania" & year >= 2016 ~ NA_real_, TRUE ~ fluCov65plus_int)) %>%
							  mutate(fluCov65plus_int = case_when(country == "Slovakia" & year >= 2016 ~ NA_real_, TRUE ~ fluCov65plus_int)) %>%
							  mutate(fluCov65plus_int = case_when(country == "Slovenia" & year >= 2016 ~ NA_real_, TRUE ~ fluCov65plus_int)) %>%
							  mutate(fluCov65plus_int = case_when(country == "Sweden" & year <= 2008 ~ NA_real_, TRUE ~ fluCov65plus_int)) %>%
							  mutate(fluCov65plus_int = case_when(country == "UK" & year <= 2008 ~ NA_real_, TRUE ~ fluCov65plus_int)) 
							  
# Forecast influenza vaccination coverage data
fluCoverage_ts <- summFluCov %>% filter(!is.na(fluCov65plus_int)) %>%
								 as_tsibble(index=year, key=country)  # define tsibble with country as key and year as index

# Fit to interpolated data					
fluCoverage_fit <- fluCoverage_ts %>% model(ets = ETS(box_cox(fluCov65plus_int, lambda = 0) ~ error("A") + trend(("Ad"), phi_range = c(0.88, 0.98))))# Fit damped model with additive error term

fluCoverage_fitted <- fluCoverage_fit %>%
						fitted() %>%
						as_tibble() %>%
						rename(flu65plus_fit = .fitted) %>%
						select(country, year, flu65plus_fit)				

fluCoverage_fcst <- fluCoverage_fit %>%	forecast(h = "23 years") %>%
					mutate(interval = hilo(.distribution, 95))

fluCoverage_extend <- as_tibble(fluCoverage_fcst) %>% 
						select(country, year, fluCov65plus_int, interval) %>%
						unnest(interval) %>%
						select(- .level) %>%
						mutate(flu65plus_L95 = pmax(0, .lower), flu65plus_U95 = pmin(100, .upper)) %>%
						rename(flu65plus_fit = fluCov65plus_int) %>%
						select(-.lower, -.upper)
					 
fluCoverage_fitted <- bind_rows(fluCoverage_fitted, fluCoverage_extend) %>% # All model values (fit and forecast)
						arrange(country, year)
						
summFluCov <- right_join(summFluCov, fluCoverage_fitted) 
					
ggplot(data = summFluCov, aes(x=year)) +
  facet_wrap(~ country) + #, labeller = label_wrap_gen(width = 2, multi_line = TRUE)) + # ,  scales = "free_y") + 
  geom_point(aes(y=fluCov65plus_int), col=ECDCcol[1]) +
  geom_line(aes(y=flu65plus_fit), col=ECDCcol[1]) +
  geom_ribbon(aes(x=year,ymin=flu65plus_L95, ymax=flu65plus_U95), fill=ECDCcol[1], alpha=.5) +
  geom_line(aes(y=flu65plus_fit), col="black") +
    ylim(0,100)+
   theme(panel.grid.major = element_blank(), 
		panel.grid.minor = element_blank(), 
		panel.background = element_blank(),
		axis.line = element_line(colour = "black"),
		text = element_text(size=14),
		axis.text.x.top = element_text(vjust = -0.5)) +
		scale_x_continuous(breaks = seq(2005, 2040, by = 10)) +
  		guides(fill=guide_legend(title="")) +
	labs(title = "Uptake of seasonal influenza vaccination in over 65s", y = "")

# Combine all predictor variables (socio-demographic and vaccination uptake)
predictors <- full_join(popProportions, childcare0to2, by=c("country", "year")) %>%
			  full_join(PCVuptake) %>%
			  full_join(summFluCov)
			  
predictors$date <- as.Date(ISOdate(predictors$year, 12, 31))  

# Values of predictor variables up until present time for fitting model
predictors_past <- as_tibble(predictors) %>%
					arrange(country, date) %>%
					select(year, country, prop0to14_int, prop65plus_int, prop80plus_int, childcare_fit, PCV7_fit, PCV10_fit, PCV13_fit, flu65plus_fit) %>%
					filter(2002 <= year & year <= 2017)

# List of age groups and vaccine groups
groups <- IPDData %>% select(country, ageGroup, groupType) %>%
					 unique() %>%
					 expand(country, ageGroup, groupType)
					 
# Values of predictor variables from now until 2040 for making forecasts (splite in half due to memory limitations)
predictors_future_1 <- as_tibble(predictors) %>%
					 arrange(country, date) %>%
					 select(year, country, prop0to14_int, prop65plus_int, prop80plus_int, childcare_fit, PCV7_fit, PCV10_fit, PCV13_fit, flu65plus_fit) %>%
					 filter(2017 < year & year <= 2040) %>%
					 right_join(groups) %>%	
					 select(country, ageGroup, groupType, year, everything()) %>%
					 filter(country %in% countryList[c(1:15)]) %>%
					 as_tsibble(key = c(country, ageGroup, groupType), index = year)
					 
predictors_future_2 <- as_tibble(predictors) %>%
					 arrange(country, date) %>%
					 select(year, country, prop0to14_int, prop65plus_int, prop80plus_int, childcare_fit, PCV7_fit, PCV10_fit, PCV13_fit, flu65plus_fit) %>%
					 filter(2017 < year & year <= 2040) %>%
					 right_join(groups) %>%	
					 select(country, ageGroup, groupType, year, everything()) %>%
					 filter(country %in% countryList[c(16:length(countryList))]) %>%
					 as_tsibble(key = c(country, ageGroup, groupType), index = year)
				
predictors_future <- rbind(predictors_future_1, predictors_future_2)

				
################# 
## Projections ##
#################
inc_predictors <- full_join(incData, predictors_past)

agg1 <- inc_predictors %>% filter(!is.na(incidence)) %>%
							#select(country, ageGroup, groupType, year, incidence, PCV10_fit)%>%
							aggregate_key(country * ageGroup * groupType,
											incidence = mean(na.omit(incidence)),
											prop0to14_int = mean(na.omit(prop0to14_int)),
											prop65plus_int = mean(na.omit(prop65plus_int)),
											prop80plus_int = mean(na.omit(prop80plus_int)),
											childcare_fit = mean(na.omit(childcare_fit)),
											PCV7_fit = mean(na.omit(PCV7_fit)),
											PCV10_fit = mean(na.omit(PCV10_fit)),
											PCV13_fit = mean(na.omit(PCV13_fit)),
											flu65plus_fit = mean(na.omit(flu65plus_fit)))
											
exoReg <- c("prop0to14_int", "prop65plus_int", "prop80plus_int", "childcare_fit", "PCV7_fit", "PCV10_fit", "PCV13_fit", "flu65plus_fit")
							
arima0 <- agg1 %>% model(lag0 = ARIMA(incidence ~ exoReg)) %>%
					reconcile() 

# Forecast without exogenous regressors					
#forecast0 <- arima0	%>%	forecast(h=10)

# Select a model to forecast
select_model <- arima0 %>% filter(country == "Austria", ageGroup == "65-74", groupType == "PCV-13")

# Plot residual errors
bind_rows(
  `Regression Errors` = residuals(select_model, type="regression"),
  `ARIMA Errors` = residuals(select_model, type="innovation"),
  .id = "type"
) %>%
  ggplot(aes(x = year, y = .resid)) +
  geom_line() +
  facet_grid(vars(type), scales = "free_y") +
  xlab("Year") + ylab(NULL)

# Filter for regressor values that correspond to the model to be forecasted
select_future <- predictors_future %>%
					aggregate_key(country * ageGroup * groupType,
											prop0to14_int = mean(na.omit(prop0to14_int)),
											prop65plus_int = mean(na.omit(prop65plus_int)),
											prop80plus_int = mean(na.omit(prop80plus_int)),
											childcare_fit = mean(na.omit(childcare_fit)),
											PCV7_fit = mean(na.omit(PCV7_fit)),
											PCV10_fit = mean(na.omit(PCV10_fit)),
											PCV13_fit = mean(na.omit(PCV13_fit)),
											flu65plus_fit = mean(na.omit(flu65plus_fit))) %>%
					#select(country, ageGroup, groupType, year, PCV7_fit, PCV10_fit, PCV13_fit)  %>%
					filter(country == "Austria", ageGroup == "65-74", groupType == "PCV-13")

# Filter for preceding incidence data
select_past <- incData %>% 
				filter(country == "Austria", ageGroup == "65-74", groupType == "PCV-13")

# Forecast using projections of exogenous regressors
forecast1 <- select_model %>% forecast(new_data = select_future) 

forecast1 %>% autoplot()
							


#########################
## Plot incidence data ##
#########################

 ## Incidence by country
  ggplot(data = incData %>% filter(country != "NA"), aes(x = year, y = incidence)) +
  geom_bar(stat="identity", aes(fill = factor(groupType, levels=c("Not serotyped", "Other", "PPV-23", "PCV-13", "PCV-10", "PCV-7")))) +
  scale_fill_manual(values = c("lightgrey", EcdcColors("qual", n=5))) +
  labs(	x = "Year", y = "Incidence per 100,000") +
  facet_wrap(~ country,  scales = "free_y") +
  theme(panel.grid.major = element_blank(), 
					  panel.grid.minor = element_blank(), 
					  panel.background = element_blank(),
					  axis.line = element_line(colour = "black"),
					  text = element_text(size=14),
					  axis.text.x = element_text(size = rel(0.9), angle = 90, vjust = -0.4)) +
					  scale_x_continuous(breaks = seq(2006.5, 2017.5, by = 1), labels =  as.character(2006:2017)) +
					  guides(fill=guide_legend(title="Vaccine type")) +
  ggtitle(paste("Incidence of IPD in 50+ year olds")) +
  theme(plot.title = element_text(hjust = 0.5))

# Proportion of cases by serotype...
propTypes <- IPDData %>% filter(country != "NA") %>%
						 group_by(country, year, groupType) %>%
						 summarise(n = n()) %>%
						 mutate(freq = n / sum(n))
# ... and of those serotyped						 
propTypes1 <- IPDData %>% filter(country != "NA", groupType != "Not serotyped") %>%
						 group_by(country, year, groupType) %>%
						 summarise(n = n()) %>%
						 mutate(freq = n / sum(n))
			  
## Plot proportion of cases by grouped serotype...
ggplot(data = propTypes, aes(x=year, y=freq)) +
  geom_bar(stat="identity", aes(fill = factor(groupType, levels=c("Not serotyped", "Other", "PPV-23", "PCV-13", "PCV-10", "PCV-7")))) +  
  scale_fill_manual(values = c("lightgrey", EcdcColors("qual", n=5))) +
  labs(	x = "Year", y = "Proportion of cases") +
  facet_wrap(~ country,  scales = "free_y") +
  theme(panel.grid.major = element_blank(), 
					  panel.grid.minor = element_blank(), 
					  panel.background = element_blank(),
					  axis.line = element_line(colour = "black"),
					  text = element_text(size=14),
					  axis.text.x = element_text(size = rel(0.9), angle = 90, vjust = -0.4)) +
					  scale_x_continuous(breaks = seq(2006.5, 2017.5, by = 1), labels =  as.character(2006:2017)) +
					  guides(fill=guide_legend(title="Vaccine type")) +
  ggtitle(paste("Confirmed cases of IPD in adults aged 50 or over")) +
  theme(plot.title = element_text(hjust = 0.5))

# ...and of those serotyped
ggplot(data = propTypes1, aes(x=year, y=freq)) +
  geom_bar(stat="identity", aes(fill = factor(groupType, levels=c("Other", "PPV-23", "PCV-13", "PCV-10", "PCV-7")))) +  
  scale_fill_manual(values = EcdcColors("qual", n=5)) +
  labs(	x = "Year", y = "Proportion of cases") +
  facet_wrap(~ country,  scales = "free_y") +
  theme(panel.grid.major = element_blank(), 
					  panel.grid.minor = element_blank(), 
					  panel.background = element_blank(),
					  axis.line = element_line(colour = "black"),
					  text = element_text(size=14),
					  axis.text.x = element_text(size = rel(0.9), angle = 90, vjust = 0.4)) +
					  scale_x_continuous(breaks = seq(2006.5, 2017.5, by = 1), labels =  as.character(2006:2017)) +
					  guides(fill=guide_legend(title="Vaccine type")) +
  ggtitle(paste("Confirmed cases of IPD in adults aged 50 or over")) +
  theme(plot.title = element_text(hjust = 0.5))
  
  