rm(list=ls())
dir <- "P:/Pneumo"
setwd(dir)
library(plyr)
library(tidyverse)
library(lubridate)
library(reshape2)
library(tsibble)
library(fable)
library(feasts)

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
				  age = as.integer(Age),
				  ageGroup = case_when(age < 5 ~ 1,		# 0-4
							  age >= 5  & age < 15 ~ 2,	# 5-14
							  age >= 15 & age < 50 ~ 3,	# 15-49
							  age >= 50 & age < 65 ~ 4,	# 50-64
							  age >= 65 & age < 75 ~ 5,	# 65-74
							  age >= 75 & age < 85 ~ 6,	# 75-84
							  age >= 85 ~ 7) %>% 
							  factor(labels = c("0-4", "5-14", "15-49", "50-64", "65-74", "75-84", "85+")),
							  
				  ## Countries
				  country = case_when(GeoCode == "UKM" ~ "SC",
											    TRUE ~ as.character(ReportingCountry)), # Recode region 'UKM' as Scotland 
				 
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
				   
			select(-DateUsedForStatisticsISO)   # Remove original date variable
			
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
			rename(year = ReportYear, country = ReportingCountry) %>% # Rename to match with IPDData
			group_by(year, country, ageGroup) %>% # Group population by age groups
			summarise(population = sum(Population))
	
# tessyPop %>% group_by(AgeGroupId, ageGroup) %>% summarise() %>% print(n=63) #Table of TESSy age categories

eurostatPop <- read_csv("EUROSTATpop.csv", col_names=TRUE)   # Read in EUROSTAT population data
scotPop <- read_csv("ScotPop.csv", col_names=TRUE)   # Source: https://www.nrscotland.gov.uk/statistics-and-data/statistics/statistics-by-theme/population/population-estimates/mid-year-population-estimates

###############
## Incidence ##
###############

IPDData <- IPDData %>% filter(Classification == "CONF") %>% # Include only confirmed cases
						left_join(tessyPop) %>%
						mutate(country = as_factor(country)) %>%
						arrange(country, year, ageGroup) %>%
						filter(SubsetEpi == 1) # Select subset for epi. to avoid double-counting lab.results

		
incData <- IPDData %>% 	filter(age >= 50)  %>%
						group_by(country, year, ageGroup, groupType, population) %>%
						#group_by(year, ageGroup, groupType, Population) %>%
						summarise(total = n()) %>% # Number of cases per year by age group, type and country
						mutate(incidence = (total/population)*100000)  # Incidence per 100,000	
						
countryNo <- 9

## Incidence by country
  ggplot(data = incData, aes(year, Incidence)) +
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
  ggtitle(paste("Rate of IPD incidence in", countries$name[countryNo])) +
  theme(plot.title = element_text(hjust = 0.5))
  
############################
## Proportion of serotype ##
############################

# Proportion of outcome by serotype
propTypes <- incData %>% filter(groupType != "Not serotyped") %>% # 
			   group_by(country, year, ageGroup, groupType) %>%
               dplyr::summarise(n = sum(Total)) %>%
			   mutate(freq = n / sum(n))
			   

propTypes1 <- propTypes %>% filter(country == countries$code[countryNo])

## Trends in proportion of types, by country
ggplot(data = propTypes1, aes(year, freq)) +
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
years <- seq(2019, max(popProjections$year)) # Extrapolate to 2019 (current year)

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
ggplot(data = popProportions, aes(x=year)) +
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
		scale_x_continuous(breaks = seq(2005, 2105, by = 20)) +
		guides(fill=guide_legend(title="")) +
	labs(title = "Proportion of population aged 0 to 14", y = "")

# Plot proportion aged 65 plus
ggplot(data = popProportions, aes(x=year)) +
  facet_wrap(~ country) + #, labeller = label_wrap_gen(width = 2, multi_line = TRUE)) + # ,  scales = "free_y") + 
  geom_line(aes(y=prop65plus_int), col=ECDCcol[1]) +
  geom_ribbon(aes(x=year,ymin=prop65plus_L95_int, ymax=prop65plus_U95_int), fill=ECDCcol[1], alpha=.5) +
   ylim(0,60)+
   theme(panel.grid.major = element_blank(), 
		panel.grid.minor = element_blank(), 
		panel.background = element_blank(),
		axis.line = element_line(colour = "black"),
		text = element_text(size=14),
		axis.text.x.top = element_text(vjust = -0.5)) +
		scale_x_continuous(breaks = seq(2005, 2105, by = 20)) +
		guides(fill=guide_legend(title="")) +
	labs(title = "Proportion of population aged over 65", y = "")


# Plot proportion aged 80 plus
ggplot(data = popProportions, aes(x=year)) +
  facet_wrap(~ country) + #, labeller = label_wrap_gen(width = 2, multi_line = TRUE)) + # ,  scales = "free_y") + 
  geom_line(aes(y=prop80plus_int), col=ECDCcol[1]) +
  geom_ribbon(aes(x=year,ymin=prop80plus_L95_int, ymax=prop80plus_U95_int), fill=ECDCcol[1], alpha=.5) +
   ylim(0,30)+
   theme(panel.grid.major = element_blank(), 
		panel.grid.minor = element_blank(), 
		panel.background = element_blank(),
		axis.line = element_line(colour = "black"),
		text = element_text(size=14),
		axis.text.x.top = element_text(vjust = -0.5)) +
		scale_x_continuous(breaks = seq(2005, 2105, by = 20)) +
		guides(fill=guide_legend(title="")) +
	labs(title = "Proportion of population aged over 80", y = "")

# Proportion of children under three not in formal childcare
childcare0to2 <- read_csv("childcare0to2.csv", col_names=TRUE) %>% # Source: EUROSTAT https://ec.europa.eu/eurostat/databrowser/view/tps00185/default/table?lang=en
					gather(countries$name, key="country", value="childcare0to2") %>%
					select(country, year, childcare0to2) %>%
					filter(country != "Liechtenstein")
					#mutate(date = as.Date(ISOdate(year, 12, 31))) # set date to 31st December of given year

# Forecast childcare
childcareTS <- as_tsibble(childcare0to2, index=year, key=country) %>%# define tsibble with country as key and year as index
			   rename(proportion = childcare0to2)
			   
fitChildcare <- childcareTS %>% model(ets = ETS(proportion))

forecastChildcare <- fitChildcare %>% forecast(h = "80 years")


ggplot(data = childcare0to2, aes(x=year)) +
  facet_wrap(~ country) + #, labeller = label_wrap_gen(width = 2, multi_line = TRUE)) + # ,  scales = "free_y") + 
  geom_line(aes(y=childcare0to2), col=ECDCcol[1]) +
  ylim(0,100)+
   theme(panel.grid.major = element_blank(), 
		panel.grid.minor = element_blank(), 
		panel.background = element_blank(),
		axis.line = element_line(colour = "black"),
		text = element_text(size=14),
		axis.text.x.top = element_text(vjust = -0.5)) +
		scale_x_continuous(breaks = seq(2006, 2018, by = 2)) +
  		guides(fill=guide_legend(title="")) +
	labs(title = "Proportion of children under three not enrolled in formal childcare", y = "")

					
# PCV coverage 
PCV7Data <- read_csv("PCV7coverage.csv", col_names=TRUE) %>% # Read in PCV7 coverage data
				melt() %>%
				as_tibble() %>%
				rename(year = variable, PCV7cov = value) %>%
				mutate(year = as.numeric(levels(year))[year], country = as_factor(country)) %>%
				arrange(country) 

PCV10Data <- read_csv("PCV10coverage.csv", col_names=TRUE) %>% # Read in PCV10 coverage data
				melt() %>%
				as_tibble() %>%
				rename(year = variable, PCV10cov = value) %>%
				mutate(year = as.numeric(levels(year))[year], country = as_factor(country)) %>%
				arrange(country) 		

PCV13Data <- read_csv("PCV13coverage.csv", col_names=TRUE) %>% # Read in PCV13 coverage data
				melt() %>%
				as_tibble() %>%
				rename(year = variable, PCV13cov = value) %>%
				mutate(year = as.numeric(levels(year))[year], country = as_factor(country)) %>%
				arrange(country) 
				
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

ggplot(data = PCV10Data, aes(x=year)) +
	facet_wrap(~country) +
	geom_line(aes(y=PCV10cov), col=ECDCcol[1]) +
	scale_color_manual(values = colours1, guide = FALSE) +
	scale_fill_manual(values = colours1, guide = FALSE) +
	theme(panel.grid.major = element_blank(), 
		panel.grid.minor = element_blank(), 
		panel.background = element_blank(),
		axis.line = element_line(colour = "black"),
		text = element_text(size=14),
		axis.text.x.top = element_text(vjust = -0.5)) +
		guides(fill=guide_legend(title=""))+
	labs(title = "PCV10 coverage in children", y = "")
					
ggplot(data = PCV13Data, aes(x=year)) +
	facet_wrap(~country) +
	geom_line(aes(y=PCV13cov), col=ECDCcol[1]) +
	scale_color_manual(values = colours1, guide = FALSE) +
	scale_fill_manual(values = colours1, guide = FALSE) +
	theme(panel.grid.major = element_blank(), 
		panel.grid.minor = element_blank(), 
		panel.background = element_blank(),
		axis.line = element_line(colour = "black"),
		text = element_text(size=14),
		axis.text.x.top = element_text(vjust = -0.5)) +
		guides(fill=guide_legend(title=""))+
	labs(title = "PCV13 coverage in children", y = "")
										
					
# PPV coverage
PPV23Data <- read_csv("PPV23coverage.csv", col_names=TRUE) # Read in PPV23 coverage data

# Influenza vaccination coverage
# Sources: https://www.ecdc.europa.eu/en/publications-data/seasonal-influenza-vaccination-europe-vaccination-recommendations-and-coverage-2007-2015 
# https://www.ecdc.europa.eu/sites/portal/files/documents/Seasonal-influenza-antiviral-use-EU-EEA-Member-States-December-2018_0.pdf 
# N.B. UK is population weighted average of four regions
fluCov55A <- read_csv("fluCov55A.csv", col_names=TRUE) %>%
			gather(countries$name, key="country", value="fluCov55A") %>%
			select(country, year, fluCov55A) %>%	# influenza vaccination coverage in over 55s, administrative sample
			mutate(country = as_factor(country))
			
fluCov59A <- read_csv("fluCov59A.csv", col_names=TRUE) %>%
			gather(countries$name, key="country", value="fluCov59A") %>%
			select(country, year, fluCov59A) %>%	# influenza vaccination coverage in over 59s, administrative sample
			mutate(country = as_factor(country))
			
fluCov60A <- read_csv("fluCov60A.csv", col_names=TRUE) %>%
			gather(countries$name, key="country", value="fluCov60A") %>%
			select(country, year, fluCov60A) %>%	# influenza vaccination coverage in over 60s, administrative sample
			mutate(country = as_factor(country))
						
fluCov60S <- read_csv("fluCov60S.csv", col_names=TRUE)  %>%
			gather(countries$name, key="country", value="fluCov60S") %>%
			select(country, year, fluCov60S) %>%	# influenza vaccination coverage in over 60s, survey sample
			mutate(country = as_factor(country))
			
fluCov65A <- read_csv("fluCov65A.csv", col_names=TRUE) %>%
			gather(countries$name, key="country", value="fluCov65A") %>%
			select(country, year, fluCov65A) %>%	# influenza vaccination coverage in over 65s, administrative sample
			mutate(country = as_factor(country))
			
fluCov65S <- read_csv("fluCov65S.csv", col_names=TRUE) %>%
			gather(countries$name, key="country", value="fluCov65S") %>%
			select(country, year, fluCov65S) %>%	# influenza vaccination coverage in over 65s, survey sample
			mutate(country = as_factor(country))

fluCoverage <- left_join(fluCov55A, fluCov59A, by = c("country", "year")) %>%
			   left_join(fluCov60A, by = c("country", "year")) %>%
			   left_join(fluCov60S, by = c("country", "year")) %>%
			   left_join(fluCov65A, by = c("country", "year")) %>%
			   left_join(fluCov65S, by = c("country", "year")) %>%
			   filter(country != "Liechtenstein") %>%
			   rowwise() %>% 
			   mutate(sumCov59 = sum(fluCov55A,fluCov59A, na.rm=TRUE)) %>%
			   mutate(sumCov60 = sum(fluCov55A,fluCov59A,fluCov60A, na.rm=TRUE)) %>%
     		   mutate(sumCov65 = sum(fluCov55A,fluCov59A,fluCov60A,fluCov65A, na.rm=TRUE)) 
			   
fluCoverage[fluCoverage == 0] <- NA
			
			
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

# Combine vaccination coverage data
vaccCoverage <- full_join(PCV7Data, PCV10Data, by=c("country", "year")) %>%
				full_join(PCV13Data) %>%
				full_join(summFluCov)

# Combine demographic and vaccination predictors
predictors <- full_join(popProportions, childcare0to2) %>%
				full_join(vaccCoverage) 
predictors$date <- as.Date(ISOdate(predictors$year, 12, 31))  
 
predictors <- as_tibble(predictors) %>%
				arrange(country, date) %>%
				select(date, year, country, prop0to14_int, prop65plus_int, prop80plus_int, childcare0to2, PCV7cov, PCV10cov, PCV13cov, fluCov55to58, fluCov59, fluCov60to64, fluCov65plus)
	
################# 
## Projections ##
#################

firstDate <- as.Date("1998-12-31")
lastDate <- as.Date("2017-12-31")

## Dynamic regression model
# IPD incidence time series (currently one age group, 
incData$date <- as.Date(ISOdate(incData$year, 12, 31))   
incData <- as_tibble(incData)

inc_ts <- incData %>% complete(date = seq(from = firstDate, to = lastDate, by = "1 year"), fill = list(value = NA)) %>%
                      as_tsibble(index = date, key = c(country, ageGroup, groupType))


#########################
## Plot incidence data ##
#########################

 ## Incidence by country
  ggplot(data = incData, aes(year, Incidence)) +
  geom_bar(stat="identity", aes(fill = groupType)) +
  labs(	x = "Year", y = "Incidence per 100,000") +
  facet_wrap(~ country) + #,  scales = "free_y") +
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
  
  
  