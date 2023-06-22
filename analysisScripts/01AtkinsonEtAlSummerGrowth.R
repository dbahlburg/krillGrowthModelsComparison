# Load packages
library(tidyverse)
library(terra)
library(sf)
library(lubridate)
# Running the krill growth model presented in Atkinson et al. (2006)
# Full reference:
# Atkinson, Angus; Shreeve, Rachael S.; Hirst, Andrew G.; Rothery, Peter; Tarling, Geraint A.;
# Pond, David W.; Korb, Rebecca E.; Murphy, Eugene J. ; Watkins, Jonathon L.. 2006 Natural 
# growth rates in Antarctic krill (Euphausia superba): II. Predictive models based on food, 
# temperature, body length, sex, and maturity stage. Limnology and Oceanography, 51 (2). 973-987. 
#--------------------------------------------------------------------------------------------------------------------#
# Simulation set-up:
# 1. initial size: 26mm
# 2. first simulation day: November 1st/day 306
# 3. last simulation day: April 15th/day 105 (total simulation time: 165 days)
# 4. when a cell is covered with sea ice, growth rate is set to 0 (indicated by chlorophyll = NA)
# 5. as a consequence, growth is delayed in areas that are still ice-covered
# load model functions
source('functions/AtkinsonEtAl2006AuxiliaryFunctions.R')
source('functions/dayOfYear.R')

# start Day:
startDay <- 306

# load environmental data
templateGrid <- subset(rast('inputData/climatologies/photoperiodClim365.tif'),1)
chlorophyllClimatology <- rast('inputData/climatologies/climaChlaRollingAv15.nc')
chlorophyllClimatology <- flip(chlorophyllClimatology, direction = 'vertical')
crs(chlorophyllClimatology) <- 'epsg:4258'
sstClimatology <- rast('inputData/climatologies/climatologySST1_365.nc')

# extract initial chlorophyll a 
initClimatologyChla <- subset(chlorophyllClimatology, startDay)

# extract initial sea surface temperature
initSST <- subset(sstClimatology, startDay) - 273.15

# initialize krill length, add krill individuals to grid containing the individuals,
# only sea cells are populated with krill
initLength <- 26
initLengthRast <- templateGrid
initLengthRast[!is.na(initLengthRast)] <- initLength

# initialize krill stage (juvenile when krill is <35mm)
initStage <- 1
initStageRast <- templateGrid
initStageRast[!is.na(initStageRast)] <- initStage

# initialize simulation environmental conditions
simulationInit <- c(initClimatologyChla,
                    initStageRast,
                    initSST,
                    initLengthRast)

# create rasters containing the simulation results for the two state variables
# length and stage
resultsRast <- initLengthRast
resultsStageRast <- initStageRast

for(i in 1:165){
  
  # determine current day of year for subsetting climatologies at correct position
  dayOfYear <- dayOfYearFunc(t = i, startDay = startDay)
  
  # in first timestep, use initial conditions to calculate growth
  if(i == 1){
    
    # returns change in length and stage
    krillUpdate <-  lapp(simulationInit, atkinsonFunc)
    
    # extract new length and stage from results
    newLength <- subset(resultsRast,i) + subset(krillUpdate, 1)
    newStage <- subset(krillUpdate, 2)
    
    # update chlorophyll a
    newChla <- subset(chlorophyllClimatology, dayOfYear)
    
    # update temperature
    newSST <- subset(sstClimatology, dayOfYear) - 273.15
    
    # create updated environment/state grid for following timestep
    newEnv <- c(newChla,
                newStage,
                newSST,
                newLength)
  }
  else{
    # returns change in length and stage
    krillUpdate <-  lapp(newEnv, atkinsonFunc)
    
    # extract new length and stage from results
    newLength <- subset(resultsRast,i) + subset(krillUpdate, 1)
    newStage <- subset(krillUpdate, 2)
    
    # update chlorophyll a
    newChla <- subset(chlorophyllClimatology, dayOfYear)
    
    # update sea surface temperature
    newSST <- subset(sstClimatology, dayOfYear) - 273.15
    
    # create updated environment/state grid for following timestep
    newEnv <- c(newChla,
                newStage,
                newSST,
                newLength)
  }
  # add new values for state variables to result-spatRaster
  add(resultsRast) <- newLength
  add(resultsStageRast) <- newStage
  
  # when body size becomes negative flag cell with -999
  newLength[newLength < 0] <- -999
  print(i)
}

# save the results
writeCDF(resultsRast, 
         varname = 'Atkinson_length_mm',
         overwrite=TRUE,
         'simulationResults/AtkinsonSummerGrowthdoy306_doy105.nc')
