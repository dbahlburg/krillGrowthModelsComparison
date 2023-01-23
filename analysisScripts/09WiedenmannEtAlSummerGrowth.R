# Load packages
library(tidyverse)
library(terra)
library(lubridate)
# Running the krill growth model presented in Wiedenmann et al. (2008)
# Full reference:
# Wiedenmann J, Cresswell K, Mangel M (2008) Temperature-dependent growth of Antarctic krill: 
# predictions for a changing climate from a cohort model. Mar Ecol Prog Ser 358:191-202.
# https://doi.org/10.3354/meps07350 
#--------------------------------------------------------------------------------------------------------------------#
# Simulation set-up:
# 1. initial size: 26mm
# 2. first simulation day: November 1st/day 306
# 3. last simulation day: April 15th/day 105 (total simulation time: 165 days)
# 4. when a cell is covered with sea ice, growth rate is set to 0
# 5. as a consequence, growth is delayed in areas that are still ice-covered 
# load model functions
source('functions/dayOfYear.R')
source('functions/WiedenmannEtAl2008AuxiliaryFunctions.R')
source('functions/AtkinsonEtAl2006AuxiliaryFunctions.R')

# start day
startDay <- 306

# load environmental data
templateGrid <- subset(rast('inputData/climatologies/photoperiodClim365.tif'),1)
chlorophyllClimatology <- rast('inputData/climatologies/climaChlaRollingAv15.nc')
chlorophyllClimatology <- flip(chlorophyllClimatology, direction = 'vertical')
crs(chlorophyllClimatology) <- 'epsg:4258'
sstClimatology <- rast('inputData/climatologies/climatologySST1_365.nc')

# initialize chlorophyll a
initClimatologyChla <- subset(chlorophyllClimatology, startDay)

# initialize sea surface temperature
initSST <- subset(sstClimatology, startDay) - 273.15

# initialize state variables
initLength <- 26
initLengthRast <- templateGrid
initLengthRast[!is.na(initLengthRast)] <- initLength

# calculate IMP based on initial environmental conditions
initMoultDay <- app(initSST, function(temp){round(exp(3.5371 - 0.5358 * log(temp + 2)))}) 
timeRast <- templateGrid
timeRast[!is.na(timeRast)] <- 1

initOldMoultDay <- app(initSST, function(temp){round(exp(3.5371 - 0.5358 * log(temp + 2)))}) 

# initialize envHistory raster that track the average environmental condition experienced
# by each individual throughout an IMP
initTemperatureHistory <- initSST*(1/initMoultDay)
initChlorophyllHistory <- initClimatologyChla*(1/initMoultDay)
initChlorophyllHistory[is.na(initChlorophyllHistory)] <- 0

# initialize simulation environmental conditions
simulationInit <- c(initLengthRast,
                    initSST,
                    initClimatologyChla,
                    timeRast,
                    initMoultDay,
                    initOldMoultDay,
                    initTemperatureHistory,
                    initChlorophyllHistory)

# create rasters containing the simulation results for the four state variables
resultsLengthRast <- initLengthRast
resultsImpRast <- initMoultDay
resultsOldMoultDay <- initOldMoultDay
resultsTemperatureHistory <- initTemperatureHistory
resultsChlorophyllHistory <- initChlorophyllHistory

# Simulate the model
for(i in 1:165){
  
  # determine current day of year for subsetting climatologies at correct position
  dayOfYear <- dayOfYearFunc(t = i, startDay = startDay)
  
  # in first timestep, use initial conditions to calculate growth
  if(i == 1){
    
    # returns change of state variables
    krillUpdate <- lapp(simulationInit, WiedenmannEtAl2008Model)
    
    # add changes to state variables
    newLength <- subset(resultsLengthRast,i) + subset(krillUpdate,1)
    newMoultDay <- subset(resultsImpRast,i) + subset(krillUpdate,2)
    newOldMoultDay <- subset(krillUpdate,3)
    newTHistory <- subset(krillUpdate,4)
    newChlaHistory <- subset(krillUpdate,5)
    
    # update chlorophyll a
    newChla <- subset(chlorophyllClimatology, dayOfYear)
    
    # update sea surface temperature
    newSST <- subset(sstClimatology,dayOfYear) - 273.15
    
    # update time raster
    timeRast[!is.na(timeRast)] <- i + 1
    
    # create updated environment/state grid for following timestep
    newEnv <- c(newLength,
                newSST,
                newChla,
                timeRast,
                newMoultDay,
                newOldMoultDay,
                newTHistory,
                newChlaHistory)
    
  }
  else{
    
    # returns change of state variables
    krillUpdate <- lapp(newEnv, WiedenmannEtAl2008Model)
    
    # add changes to state variables
    newLength <- subset(resultsLengthRast,i) + subset(krillUpdate,1)
    newMoultDay <- subset(resultsImpRast,i) + subset(krillUpdate,2)
    newOldMoultDay <- subset(krillUpdate,3)
    newTHistory <- subset(krillUpdate,4)
    newChlaHistory <- subset(krillUpdate,5)
    
    # update chlorophyll a
    newChla <- subset(chlorophyllClimatology, dayOfYear)
    
    # update sea surface temperature
    newSST <- subset(sstClimatology,dayOfYear) - 273.15
    
    # update time raster
    timeRast[!is.na(timeRast)] <- i + 1
    
    # create updated environment/state grid for following timestep
    newEnv <- c(newLength,
                newSST,
                newChla,
                timeRast,
                newMoultDay,
                newOldMoultDay,
                newTHistory,
                newChlaHistory)
    
  }
  
  # add new values for state variables to result-spatRaster
  add(resultsLengthRast) <- newLength
  add(resultsImpRast) <- newMoultDay
  add(resultsOldMoultDay) <- newOldMoultDay
  add(resultsTemperatureHistory) <- newTHistory
  add(resultsChlorophyllHistory) <- newChlaHistory
  
  # when body size becomes negative flag cell with -999
  newLength[newLength < 0] <- -999
  print(i)
}

# save results
writeCDF(resultsLengthRast, 
         varname = 'Wiedenmann_length_mm',
         overwrite=TRUE,
         'simulationResults/WiedenmannSummerGrowthdoy306_doy105.nc')
