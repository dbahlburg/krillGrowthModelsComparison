# Load packages
library(tidyverse)
library(terra)
library(lubridate)
# Load auxiliary functions
source('functions/FachEtAl2002AuxiliaryFunctions.R')

# Implementation of the krill growth mdoel presented by Fach et al. (2002) MEPS
# Full reference:
# Fach, Bettina & Hofmann, Eileen & Murphy, Eugene. (2002). Modeling studies of Antarctic krill 
# Euphausia superba survival during transport across the Scotia Sea. Marine Ecology-progress Series -
# 231. 187-203. doi: 10.3354/meps231187.
# This model is a modified verion of Hofmann and Lascara (2000). Therefore, the same auxiliary functions are used
# but some modifications are included in the model code:
# 1. There is now heterotrophic food source which can be exploited by krill >18mm
# 2. Growth rates are scaled by ambient water temperature, according to a Q10 functional relationship
# Similar to Hofmann and Lascara (2000), the model uses lookup tables for size-specific traits of the modelled individuals. It is 
# recommended to use the lookup table-version of the model for faster simulations.
#--------------------------------------------------------------------------------------------------------------------#
# DEFINE PARAMETER VALUES
# conversion rate chlorophyll to carbon
chlorophyllToCarbon <- 50

# assimilation efficiency of ingested carbon
assimilationRate <- 0.8

# vector containing the different size classes (used in the model to assign current size class of individual)
sizeClassesVec <- seq(1,90, by = 0.25)

# Temporal settings for simulation (tmax, dt, startDay)
# define tmax and resolution
dt <- 1
startDay <- 306
#--------------------------------------------------------------------------------------------------------------------#
# Import lookup tables
lookUpTableTimeVariables <- read.csv('inputData/originalResults/timeDependentVariables.csv')
lookUpTableSizeClassTraits <- read.csv('inputData/originalResults/sizeClassTraits.csv')

templateGrid <- subset(rast('inputData/climatologies/photoperiodClim365.tif'),1)
chlorophyllClimatology <- rast('inputData/climatologies/climaChlaRollingAv15.nc')
chlorophyllClimatology <- flip(chlorophyllClimatology, direction = 'vertical')
crs(chlorophyllClimatology) <- 'epsg:4258'
sstClimatology <- rast('inputData/climatologies/climatologySST1_365.nc')

# Initialize state variables
initMass <- massAtLength(inputLength = 26, massType = 'carbonMass')
initMassRast <- templateGrid
initMassRast[!is.na(initMassRast)] <- initMass

# initialize environmental variables
initClimatologyChla <- subset(chlorophyllClimatology, startDay)

initIceAlgae <- templateGrid
initIceAlgae[!is.na(initIceAlgae)] <- 0

timeRast <- templateGrid
timeRast[!is.na(timeRast)] <- 1

heterotrophicCarbonRast <- templateGrid
heterotrophicCarbonRast[!is.na(heterotrophicCarbonRast)] <- 0

initSST <- subset(sstClimatology, startDay) - 273.15

# initialize simulation environmental conditions
simulationInit <- c(initMassRast,
                    timeRast,
                    initClimatologyChla,
                    initSST,
                    heterotrophicCarbonRast,
                    initIceAlgae)

# raster storing the simulation results
resultsRast <- initMassRast

for(i in 1:165){
  
  # determine current day of year for subsetting climatologies at correct position
  dayOfYear <- dayOfYearFunc(t = i, startDay = startDay)
  
  # in first timestep, use initial conditions to calculate growth
  if(i == 1){
    
    # add change to state variable
    newMass <- subset(resultsRast,i) + lapp(simulationInit, rk4stepFach)
    
    # update chlorophyll a
    newChla <- subset(chlorophyllClimatology,dayOfYear)
    
    # update sea surface temperature
    newSST <- subset(sstClimatology, dayOfYear) - 273.15
    
    # update time rast
    timeRast[!is.na(timeRast)] <- i + 1
    
    # create updated environment/state grid for following timestep
    newEnv <- c(newMass,
                timeRast,
                newChla,
                newSST,
                heterotrophicCarbonRast,
                initIceAlgae)
  }
  else{
    
    # add change to state variable
    newMass <- subset(resultsRast,i) + lapp(newEnv, rk4stepFach)
    
    # update chlorophyll a
    newChla <- subset(chlorophyllClimatology,dayOfYear)
    crs(newChla) <- 'epsg:4258'
    
    # update sea surface temperature
    newSST <- subset(sstClimatology, dayOfYear) - 273.15
    
    # update time rast
    timeRast[!is.na(timeRast)] <- i + 1
    
    # create updated environment/state grid for following timestep
    newEnv <- c(newMass,
                timeRast,
                newChla,
                newSST,
                heterotrophicCarbonRast,
                initIceAlgae)
    
  }
  
  # add new values for state variables to result-spatRaster
  add(resultsRast) <- newMass
  
  # when body mass becomes negative flag cell with -999
  newMass[newMass < 0] <- -999
  print(i)
}

# At this point, results are stored in carbon mass
# Convert from mass to length
resultsLength <- rast()
for(k in 1:166){
  newResultsLength <- app(subset(resultsRast,k), lengthAtMass)
  add(resultsLength) <- newResultsLength
  print(k)
}

# save results
writeCDF(resultsLength, 
         varname = 'FachEtAl2002_length_mm',
         overwrite=TRUE,
         'simulationResults/FachEtAl2002SummerGrowthdoy306_doy105.nc')
