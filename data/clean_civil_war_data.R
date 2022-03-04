#############################
# Clean Neilsen et al. data
#############################


# libraries ---------------------------------------------------------------

library(tidyverse) # for grammar
library(missRanger)    # for imputing missing values


# data --------------------------------------------------------------------

raw_data <- 
  haven::read_dta("cwdata.dta")


# narrow down to variables of interest ------------------------------------

raw_data <-
  raw_data %>%
   # keep only recipient countries:
  filter(inmysample==1) %>%
  transmute(
    # the response variable: civil war coded 1 if >= 25 battle deaths in a year (PRIO)
    civil_war = prio,
    
    # the causal variable of interest: change in aid received per gdp
    aid_change = lagaidgdpchange,
    
    # counfounding variables
    x1 = PTSave_filled, # HR violations
    x2 = assassinbanks, # no. of assassinations
    x3 = riotsbanks, # no. of riots
    x4 = strikesbanks, # no. of strikes
    x5 = demonstrationsbanks, # no. of demonstrations
    x6 = infantmort, # infant mortaility
    x7 = nciv, # no. of contiguous neighbors experiencing civil war
    x8 = partautocracy, # partial autocracy
    x9 = partdemocracy, # partial democracy
    x10 = factionaldemoc, # factional democracy
    x11 = fulldemocracy, # full democracy
    x12 = ln_rgdpc, # log of gdp/capita
    x13 = ln_population, # log of population
    x14 = oil, # oil exportation
    x15 = instab, # political instability
    x16 = ethfrac, # ethnic fractionalization
    x17 = relfrac, # religious fractionalization
    x18 = as.numeric(year <= 1991), # cold war
    x19 = logmtn, # mountains
    x20 = ncontig, # noncontiguous
    year = year,
    country = countryname
  )

 # save vector of names
covnames <- c(
  x1 = "HR Violations",
  x2 = "Assassinations",
  x3 = "Riots",
  x4 = "General Strikes",
  x5 = "Demonstrations",
  x6 = "Infant Mortality",
  x7 = "Bad Neighbors",
  x8 = "Partial Autocracy",
  x9 = "Partial Democracy",
  x10 = "Factional Democracy",
  x11 = "Full Democracy",
  x12 = "ln GDP/capita",
  x13 = "ln Population",
  x14 = "Oil Exportation",
  x15 = "Political Instability",
  x16 = "Ethnic Fractionalization",
  x17 = "Religious Fractionalization",
  x18 = "Cold War",
  x19 = "ln Mountainous",
  x20 = "Contiguous"
)

# multiple impuation ------------------------------------------------------

# fast missing value imputation by chained random forests

imp_data <-
  missRanger(
    data = raw_data
  )


# save raw and imputed data -----------------------------------------------

save(raw_data, imp_data, covnames, file = "data/civil_war_data.R")
