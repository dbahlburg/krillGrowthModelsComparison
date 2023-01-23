# Load packages
library(tidyverse)
library(terra)
library(lubridate)
# Implementation of the EMC model presented by Constable and Kawaguchi (2017) in the ICES Journal of Marine Science
# Full reference:
# Andrew John Constable, So Kawaguchi, Modelling growth and reproduction of Antarctic krill, Euphausia superba, 
# based on temperature, food and resource allocation amongst life history functions, ICES Journal of Marine Science, 
# Volume 75, Issue 2, March-April 2018, Pages 738–750, https://doi.org/10.1093/icesjms/fsx190
#--------------------------------------------------------------------------------------------------------------------#
# load model functions
source('functions/dayOfYear.R')
source('functions/ConstableKawaguchi2017AuxiliaryFunctions.R')

# start day
startDay <- 306

# Required input variables:
# 1. poc (chlorophyllToCarbon = 50)
# 2. sst
# 3. photoperiod
# the growth function expects the following input arguments
# bodyLength, stage, bodyCReserves, gonadCReserves, gPhase, lastMoult, moultDay, impSinceSwitch, temperature, pocConc, photoperiod, time

# load environmental data
templateGrid <- subset(rast('inputData/climatologies/photoperiodClim365.tif'),1)

# chlorophyll needs to be converted into poc. 
# as used in the original manuscript, I use a conversion constant of 50
chlorophyllToCarbon <- 50
PocClimatology <- rast('inputData/climatologies/climaChlaRollingAv15.nc') * chlorophyllToCarbon
PocClimatology <- flip(PocClimatology, direction = 'vertical')
crs(PocClimatology) <- 'epsg:4258'
sstClimatology <- rast('inputData/climatologies/climatologySST1_365.nc')
photoClimatologies <- rast('inputData/climatologies/photoperiodClim365.tif')

# initialize environmental variables
initClimatologyPoc <- subset(PocClimatology, startDay)
initSST <- subset(sstClimatology, startDay) - 273.15
initPhoto <- subset(photoClimatologies,1)
initTimeRast <- templateGrid
initTimeRast[!is.na(initTimeRast)] <- 1
# ------------------------------------------------------------------------- #
# initialize state variables
# length
initLength <- 26
initLengthRast <- templateGrid
initLengthRast[!is.na(initLengthRast)] <- initLength 

# stage, initialized as immature females
initStageRast <- templateGrid
initStageRast[!is.na(initStageRast)] <- 2

# body carbon reserves for a healthy individual of initLength
initBodyCReserves <- templateGrid
initBodyCReserves[!is.na(initBodyCReserves)] <- dryMassToCarbon(dryMass_mg = dryMassAtLength(bodyLength = initLength, stage = 2), tissue = 'bodyReserves')

# gonad carbon reserves are zero
initGonadCReserves <- templateGrid
initGonadCReserves[!is.na(initGonadCReserves)] <- 0

# init gPhase is "immature"
initGPhase <- templateGrid
initGPhase[!is.na(initGPhase)] <- 0

# individuals have moulted at t0
initLastMoult <- templateGrid
initLastMoult[!is.na(initLastMoult)] <- 0

# intermoult periods determined from initial temperature distribution
initIMP <- app(initSST, determineIMP)

# imp-counter for gPhase-Function
initImpSinceSwitch <- templateGrid
initImpSinceSwitch[!is.na(initImpSinceSwitch)] <- 1

# no eggs produced at t0
initEggsProduced <- templateGrid
initEggsProduced[!is.na(initEggsProduced)] <- 0

# order: bodyLength, stage, bodyCReserves, gonadCReserves, gPhase, lastMoult, moultDay, impSinceSwitch, temperature, pocConc, photoperiod, time
simulationInit <- c(initLengthRast,
                    initStageRast,
                    initBodyCReserves,
                    initGonadCReserves,
                    initGPhase,
                    initLastMoult,
                    initIMP,
                    initImpSinceSwitch,
                    initSST,
                    initClimatologyPoc,
                    initPhoto,
                    initTimeRast)

# create rasters containing the simulation results for the state variables
resultsLength <- initLengthRast
resultsStage <- initStageRast
resultsBodyCReserves <- initBodyCReserves
resultsGonadCReserves <- initGonadCReserves
resultsGPhase <- initGPhase
resultsLastMoult <- initLastMoult
resultsIMP <- initIMP
resultsImpSinceSwitch <- initImpSinceSwitch

# Run simulation
for(i in 1:165){
  
  # determine current day of year for subsetting climatologies at correct position
  dayOfYear <- dayOfYearFunc(i, startDay = startDay)
  
  # in first timestep, use initial conditions to calculate growth
  if(i == 1){
    
    # order of returned objects: dL, newStage, dB, dG, newGPhase, newLastMoult, newMoult, newImpSinceSwitch, eggsProduced
    change <- lapp(simulationInit, constableKawaguchiGrowth)
    
    # add changes to state variables
    newLength <- subset(resultsLength,i) + subset(change,1)
    newStage <- subset(change,2)
    newBodyCReserves <- subset(resultsBodyCReserves,i) + subset(change,3)
    newGonadCReserves <- subset(resultsGonadCReserves,i) + subset(change,4)
    newGPhase <- subset(change,5)
    newLastMoult <- subset(change,6)
    newIMP <- subset(resultsIMP,i) + subset(change,7)
    newImpSinceSwitch <- subset(change,8)
    
    # update poc concentration
    newPoc <- subset(PocClimatology,dayOfYear)
    
    # update time rast
    timeRast <- initTimeRast
    timeRast[!is.na(timeRast)] <- i + 1
    
    # update photoperiod
    newPhoto <- subset(photoClimatologies, dayOfYear)
    
    # update sea surface temperature
    newSST <- subset(sstClimatology,dayOfYear) - 273.15
    
    # create new environment/state raster to calculate growth at next timestep
    newEnv <- c(newLength,
                newStage,
                newBodyCReserves,
                newGonadCReserves,
                newGPhase,
                newLastMoult,
                newIMP,
                newImpSinceSwitch,
                newSST,
                newPoc,
                newPhoto,
                timeRast)
  }
  else{
    
    # order of returned objects: dL, newStage, dB, dG, newGPhase, newLastMoult, newMoult, newImpSinceSwitch, eggsProduced
    change <- lapp(newEnv, constableKawaguchiGrowth)
    
    # add changes to state variables
    newLength <- subset(resultsLength,i) + subset(change,1)
    newStage <- subset(change,2)
    newBodyCReserves <- subset(resultsBodyCReserves,i) + subset(change,3)
    newGonadCReserves <- subset(resultsGonadCReserves,i) + subset(change,4)
    newGPhase <- subset(change,5)
    newLastMoult <- subset(change,6)
    newIMP <- subset(resultsIMP,i) + subset(change,7)
    newImpSinceSwitch <- subset(change,8)
    
    # update poc concentration
    newPoc <- subset(PocClimatology,dayOfYear)
    
    # update time rast
    timeRast <- initTimeRast
    timeRast[!is.na(timeRast)] <- i + 1
    
    # update photoperiod
    newPhoto <- subset(photoClimatologies, dayOfYear)
    
    # update sea surface temperature
    newSST <- subset(sstClimatology,dayOfYear) - 273.15
    
    # create new environment/state raster to calculate growth at next timestep
    newEnv <- c(newLength,
                newStage,
                newBodyCReserves,
                newGonadCReserves,
                newGPhase,
                newLastMoult,
                newIMP,
                newImpSinceSwitch,
                newSST,
                newPoc,
                newPhoto,
                timeRast)
  }
  
  # add new values for state variables to result-spatRaster
  add(resultsLength) <- newLength
  add(resultsStage) <- newStage
  add(resultsBodyCReserves) <- newBodyCReserves
  add(resultsGonadCReserves) <- newGonadCReserves
  add(resultsGPhase) <- newGPhase
  add(resultsLastMoult) <- newLastMoult
  add(resultsIMP) <- newIMP
  add(resultsImpSinceSwitch) <- newImpSinceSwitch
  
  # when body size becomes negative flag cell with -999
  newLength[newBodyCReserves < 0] <- -999
  print(i)
}

# save results
writeCDF(resultsLength, 
         varname = 'ConstableKawaguchi_length_mm',
         overwrite=TRUE,
         'simulationResults/ConstableKawaguchiSummerGrowthdoy306_doy105.nc')



