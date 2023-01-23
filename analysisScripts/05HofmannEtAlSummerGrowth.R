# Load packages
library(tidyverse)
library(terra)
library(lubridate)
# Load auxiliary functions
source('functions/HofmannLascara2000AuxiliaryFunctions.R')
# Implementation of the krill growth mdoel presented by Hofmann and Lascara (2000) MEPS
# Full reference:
# Hofmann EE, Lascara CM (2000) Modeling the growth dynamics of Antarctic krill Euphausia 
# superba. Mar Ecol Prog Ser 194:219-231. https://doi.org/10.3354/meps194219 
# This is a model version which uses lookup tables for size-specific traits of the modelled individuals.
# In a different version, these traits are explicitly calculated and returned by functions within 
# the simulation. The latter one allows for higher precision at the cost of additional computation time.
# Using lookup tables reduces the simulation time from 420s to 240s. The results from 
# both simulations are virtually indistinguishable. Therefore, it is recommended to use
# the lookup table-version of the model for faster simulations.
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

# Initialize state variables (start length 26mm)
initMass <- massAtLength(inputLength = 26, massType = 'carbonMass')
initMassRast <- templateGrid
initMassRast[!is.na(initMassRast)] <- initMass

# intialize chlorophyll a
initClimatologyChla <- subset(chlorophyllClimatology, startDay)

# initialize ice algae (zero)
initIceAlgae <- templateGrid
initIceAlgae[!is.na(initIceAlgae)] <- 0

# initialize time rast
timeRast <- templateGrid
timeRast[!is.na(timeRast)] <- 1

# initialize simulation environmental conditions
simulationInit <- c(initMassRast,
                    timeRast,
                    initClimatologyChla,
                    initIceAlgae)

# create raster containing the simulation results
resultsRast <- initMassRast

for(i in 1:165){
  
  # determine current day of year for subsetting climatologies at correct position
  dayOfYear <- dayOfYearFunc(t = i, startDay = startDay)
  
  # in first timestep, use initial conditions to calculate growth
  if(i == 1){
    
    # add change to state variable
    newMass <- subset(resultsRast,i) + lapp(simulationInit, rk4stepHofmann)
    
    # update chlorophyll a
    newChla <- subset(chlorophyllClimatology,dayOfYear)
    
    # update time rast
    timeRast[!is.na(timeRast)] <- i + 1
    
    # create updated environment/state grid for following timestep
    newEnv <- c(newMass,
                timeRast,
                newChla,
                initIceAlgae)
  }
  else{
    
    # add change to state variable
    newMass <- subset(resultsRast,i) + lapp(newEnv, rk4stepHofmann)
    
    # update chlorophyll a
    newChla <- subset(chlorophyllClimatology,dayOfYear)
    
    # update time rast
    timeRast[!is.na(timeRast)] <- i + 1
    
    # create updated environment/state grid for following timestep
    newEnv <- c(newMass,
                timeRast,
                newChla,
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
         varname = 'HofmannLascara_length_mm',
         overwrite=TRUE,
         'simulationResults/HofmannLascaraSummerGrowthdoy306_doy105.nc')

