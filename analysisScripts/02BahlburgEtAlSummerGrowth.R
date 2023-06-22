# Load packages
library(tidyverse)
library(terra)
library(lubridate)
# Running the krill growth model presented in Bahlburg et al. (2021)
# Full reference:
# Bahlburg, D. , Meyer, B. and Berger, U. (2021): The impact of seasonal regulation of metabolism on the 
# life history of Antarctic krill , Ecological Modelling, 442 . doi: 10.1016/j.ecolmodel.2021.109427 
#--------------------------------------------------------------------------------------------------------------------#
# load model functions
source('functions/BahlburgEtAl2021AuxiliaryFunctions.R')

# start day
startDay <- 306

# load environmental data
templateGrid <- subset(rast('inputData/climatologies/photoperiodClim365.tif'),1)
chlorophyllClimatology <- rast('inputData/climatologies/climaChlaRollingAv15.nc')
chlorophyllClimatology <- flip(chlorophyllClimatology, direction = 'vertical')
crs(chlorophyllClimatology) <- 'epsg:4258'
sstClimatology <- rast('inputData/climatologies/climatologySST1_365.nc')
photoClimatologies <- rast('inputData/climatologies/photoperiodClim365.tif')

# extract initial chlorophyll a 
initClimatologyChla <- subset(chlorophyllClimatology, startDay)

# extract initial sea surface temperature (in Kelvin)
initSST <- subset(sstClimatology, startDay) - 273.15

# extract initial photoperiod
initPhoto <- subset(photoClimatologies, startDay)

# initialize state variables length, time, reproBuffer, eggBuffer
initLength <- 26
initLengthRast <- templateGrid
initLengthRast[!is.na(initLengthRast)] <- initLength

initTimeRast <- templateGrid
initTimeRast[!is.na(initTimeRast)] <- 1

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
                    initPhoto,
                    initClimatologyChla)

# create rasters containing the simulation results for the four state variables
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
    
    # returns change in length, egg buffer, repro buffer and released eggs
    change <- lapp(simulationInit, bahlburgEtAlGrowth)
    
    # add changes to state variables
    newLength <- subset(resultsLength,i) + subset(change,1)
    newEggBuffer <- subset(resultsEggBuffer,i) + subset(change,2)
    newReproBuffer <- subset(resultsReproBuffer,i) + subset(change, 3)
    newReleasedEggs <- subset(change, 4)
    
    # update chlorophyll a
    newChla <- subset(chlorophyllClimatology,dayOfYear)
    
    # update time raster
    timeRast <- initTimeRast
    timeRast[!is.na(timeRast)] <- i + 1
    
    # update photoperiod
    newPhoto <- subset(photoClimatologies, dayOfYear)
    
    # update sea surface temperature
    newSST <- subset(sstClimatology,dayOfYear) - 273.15
    
    # create updated environment/state grid for following timestep
    newEnv <- c(newLength,
                timeRast,
                newEggBuffer,
                newReproBuffer,
                newSST,
                newPhoto,
                newChla)
  }else{
    
    # returns change in length, egg buffer, repro buffer and released eggs
    change <- lapp(newEnv, bahlburgEtAlGrowth)
    
    # add changes to state variables
    newLength <- subset(resultsLength,i) + subset(change,1)
    newEggBuffer <- subset(resultsEggBuffer,i) + subset(change,2)
    newReproBuffer <- subset(resultsReproBuffer,i) + subset(change, 3)
    newReleasedEggs <- subset(resultsEggs,i) + subset(change, 4)
    
    # update chlorophyll a
    newChla <- subset(chlorophyllClimatology,dayOfYear)
    
    # update time raster
    timeRast <- initTimeRast
    timeRast[!is.na(timeRast)] <- i + 1
    
    # update photoperiod
    newPhoto <- subset(photoClimatologies, dayOfYear)
    
    # update sea surface temperature
    newSST <- subset(sstClimatology,dayOfYear) - 273.15
    
    # create updated environment/state grid for following timestep
    newEnv <- c(newLength,
                timeRast,
                newEggBuffer,
                newReproBuffer,
                newSST,
                newPhoto,
                newChla)
  }
  
  # add new values for state variables to result-spatRaster
  add(resultsLength) <- newLength
  add(resultsEggBuffer) <- newEggBuffer
  add(resultsReproBuffer) <- newReproBuffer
  add(resultsEggs) <- newReleasedEggs
  
  # when body size becomes negative flag cell with -999
  newLength[newLength < 0] <- -999
  print(i)
}

# save results
writeCDF(resultsLength, 
         varname = 'Bahlburg_length_mm',
         overwrite=TRUE,
         'simulationResults/BahlburgSummerGrowthdoy306_doy105.nc')
