library(tidyverse)
library(terra)
library(lubridate)
# Running the krill growth model presented in Tarling et al. (2006)
# Full reference:
#--------------------------------------------------------------------------------------------------------------------#
# Simulation set-up:
# 1. initial size: 15mm
# 2. first simulation day: September 1st/day 244
# 3. last simulation day: April 15th/day 105 (total simulation time: 227 days)
# 4. when a cell is covered with sea ice, growth rate is set to 0
# 5. as a consequence, growth is delayed in areas that are still ice-covered in early September
# load model functions
source('functions/dayOfYear.R')
source('functions/TarlingEtAl2006AuxiliaryFunctions.R')
source('functions/AtkinsonEtAl2006AuxiliaryFunctions.R')

# load environmental data
templateGrid <- subset(rast('inputData/climatologies/photoperiodClim365.tif'),1)
chlorophyllClimatology <- rast('inputData/climatologies/climaChlaRollingAv15.nc')
chlorophyllClimatology <- flip(chlorophyllClimatology, direction = 'vertical')
crs(chlorophyllClimatology) <- 'epsg:4258'
sstClimatology <- rast('inputData/climatologies/climatologySST1_365.nc')

# start day
startDay <- 306

# initialize chlorophyll a and sea surface temperature
initClimatologyChla <- subset(chlorophyllClimatology, startDay)
initSST <- subset(sstClimatology, startDay) - 273.15

# initialize krill, populate the grid
initLength <- 26
initLengthRast <- templateGrid
initLengthRast[!is.na(initLengthRast)] <- initLength

# initialize the krill stage (juvenile for <35mm)
initStage <- 1
initStageRast <- templateGrid
initStageRast[!is.na(initStageRast)] <- initStage

# calculate first IMP based on initial environmental conditions
initMoultDay <- lapp(c(initLengthRast, initSST), TarlingEtAl2006IMPRounded) 
timeRast <- templateGrid
timeRast[!is.na(timeRast)] <- 1

initOldMoultDay <- lapp(c(initLengthRast, initSST), TarlingEtAl2006IMPRounded) 

# initiate spatraster storing the experienced temperature- and chlorophyll a history
# by each individual
initTemperatureHistory <- initSST*(1/initMoultDay)
initChlorophyllHistory <- initClimatologyChla*(1/initMoultDay)
initChlorophyllHistory[is.na(initChlorophyllHistory)] <- 0

# initialize simulation environmental conditions
simulationInit <- c(initLengthRast,
                    initStageRast,
                    initSST,
                    initClimatologyChla,
                    timeRast,
                    initMoultDay,
                    initOldMoultDay,
                    initTemperatureHistory,
                    initChlorophyllHistory)

# create rasters containing the simulation results for the four state variables
resultsLengthRast <- initLengthRast
resultsStageRast <- initStageRast
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
    krillUpdate <- lapp(simulationInit, TarlingEtAl2006Model)
    
    # add changes to state variables
    newLength <- subset(resultsLengthRast,i) + subset(krillUpdate,1)
    newStage <- subset(krillUpdate,2)
    newMoultDay <- subset(resultsImpRast,i) + subset(krillUpdate,3)
    newOldMoultDay <- subset(krillUpdate,4)
    newTHistory <- subset(krillUpdate,5)
    newChlaHistory <- subset(krillUpdate,6)
    
    # update chlorophyll a
    newChla <- subset(chlorophyllClimatology, dayOfYear)
    
    # update sea surface temperature
    newSST <- subset(sstClimatology,dayOfYear) - 273.15
    
    # update time raster
    timeRast[!is.na(timeRast)] <- i + 1
    
    # create updated environment/state grid for following timestep
    newEnv <- c(newLength,
                newStage,
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
    krillUpdate <- lapp(newEnv, TarlingEtAl2006Model)
    
    # add changes to state variables
    newLength <- subset(resultsLengthRast,i) + subset(krillUpdate,1)
    newStage <- subset(krillUpdate,2)
    newMoultDay <- subset(resultsImpRast,i) + subset(krillUpdate,3)
    newOldMoultDay <- subset(krillUpdate,4)
    newTHistory <- subset(krillUpdate,5)
    newChlaHistory <- subset(krillUpdate,6)
    
    # update chlorophyll a
    newChla <- subset(chlorophyllClimatology, dayOfYear)
    
    # update sea surface temperature
    newSST <- subset(sstClimatology,dayOfYear) - 273.15
    
    # update time raster
    timeRast[!is.na(timeRast)] <- i + 1
    
    # create updated environment/state grid for following timestep
    newEnv <- c(newLength,
                newStage,
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
  add(resultsStageRast) <- newStage
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
         varname = 'Tarling_length_mm',
         overwrite=TRUE,
         'simulationResults/TarlingSummerGrowthdoy306_doy105.nc')
