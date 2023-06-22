# Load packages
library(tidyverse)
library(terra)
library(lubridate)
# Running the krill growth model presented in Jager and Ravagnan (2015)
# Full reference:
# Jager T, Ravagnan E (2015) Parameterising a generic model for the dynamic energy budget of Antarctic 
# krill Euphausia superba. Mar Ecol Prog Ser 519:115-128. https://doi.org/10.3354/meps11098
#--------------------------------------------------------------------------------------------------------------------#
# load model functions
source('functions/JagerRavagnanAuxiliaryFunctions.R')

# start day
startDay <- 306

# load environmental data
templateGrid <- subset(rast('inputData/climatologies/photoperiodClim365.tif'),1)
sstClimatology <- rast('inputData/climatologies/climatologySST1_365.nc')
chlorophyllClimatology <- rast('inputData/climatologies/climaChlaRollingAv15.nc')
chlorophyllClimatology <- flip(chlorophyllClimatology, direction = 'vertical')
crs(chlorophyllClimatology) <- 'epsg:4258'

# initialize chlorophyll a
initClimatologyChla <- subset(chlorophyllClimatology, startDay)
initSST <- subset(sstClimatology, startDay) - 273.15

# initialize raster storing time information
initTimeRast <- templateGrid
initTimeRast[!is.na(initTimeRast)] <- 1

# initialize state variables
initLength <- 26
initLengthRast <- templateGrid
initLengthRast[!is.na(initLengthRast)] <- initLength

initReproBuffer <- templateGrid
initReproBuffer[!is.na(initReproBuffer)] <- 0

initEggBuffer <- templateGrid
initEggBuffer[!is.na(initEggBuffer)] <- 0

# initialize simulation environmental conditions
simulationInit <- c(initLengthRast,
                    initTimeRast,
                    initEggBuffer,
                    initReproBuffer,
                    initSST,
                    initClimatologyChla)

# create rasters containing the simulation results of the state variables
resultsLength <- initLengthRast
resultsEggBuffer <- initEggBuffer
resultsReproBuffer <- initReproBuffer
resultsEggs <- initEggBuffer

# Run simulation
for(i in 1:165){
  
  # determine current day of year for subsetting climatologies at correct position
  dayOfYear <- dayOfYearFunc(i, startDay = startDay)
  
  # in first timestep, use initial conditions to calculate growth
  if(i == 1){
    
    # returns change of state variables
    change <- lapp(simulationInit, jagerRavagnanGrowth)
    
    # add change to state variable
    newLength <- subset(resultsLength,i) + subset(change,1)
    newEggBuffer <- subset(resultsEggBuffer,i) + subset(change,2)
    newReproBuffer <- subset(resultsReproBuffer,i) + subset(change, 3)
    newReleasedEggs <- subset(change, 4)
    
    # update chlorophyll a
    newChla <- subset(chlorophyllClimatology,dayOfYear)
    
    # update time rast
    timeRast <- initTimeRast
    timeRast[!is.na(timeRast)] <- i + 1
    
    # update sea surface temperature
    newSST <- subset(sstClimatology,dayOfYear) - 273.15
    
    # create updated environment/state grid for following timestep
    newEnv <- c(newLength,
                timeRast,
                newEggBuffer,
                newReproBuffer,
                newSST,
                newChla)
  }
  else{
    
    # returns change of state variables
    change <- lapp(newEnv, jagerRavagnanGrowth)
    
    # add change to state variable
    newLength <- subset(resultsLength,i) + subset(change,1)
    newEggBuffer <- subset(resultsEggBuffer,i) + subset(change,2)
    newReproBuffer <- subset(resultsReproBuffer,i) + subset(change, 3)
    newReleasedEggs <- subset(resultsEggs,i) + subset(change, 4)
    
    # update chlorophyll a
    newChla <- subset(chlorophyllClimatology,dayOfYear)
    
    # update time rast
    timeRast <- initTimeRast
    timeRast[!is.na(timeRast)] <- i + 1
    
    # update sea surface temperature
    newSST <- subset(sstClimatology,dayOfYear) - 273.15
    
    # create updated environment/state grid for following timestep
    newEnv <- c(newLength,
                timeRast,
                newEggBuffer,
                newReproBuffer,
                newSST,
                newChla)
  }
  
  # add new values for state variables to result-spatRaster
  add(resultsLength) <- newLength
  add(resultsEggBuffer) <- newEggBuffer
  add(resultsReproBuffer) <- newReproBuffer
  add(resultsEggs) <- newReleasedEggs
  
  # when body length becomes negative flag cell with -999
  newLength[newLength < 0] <- -999
  print(i)
}

# save results
writeCDF(resultsLength, 
         varname = 'JagerRavagnan_length_mm',
         overwrite=TRUE,
         'simulationResults/JagerRavagnanSummerGrowthdoy306_doy105.nc')
