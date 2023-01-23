# Load packages
library(tidyverse)
library(terra)
library(lubridate)
# Running the krill growth model presented in Ryabov et al. (2017)
# Full reference:
# Ryabov, A., de Roos, A., Meyer, B. et al. Competition-induced starvation drives large-scale population cycles in 
# Antarctic krill. Nat Ecol Evol 1, 0177 (2017). https://doi.org/10.1038/s41559-017-0177
#--------------------------------------------------------------------------------------------------------------------#
# Simulation set-up:
# 1. initial size: 26mm
# 2. first simulation day: November 1st/day 306
# 3. last simulation day: April 15th/day 105 (total simulation time: 165 days)
# 4. when a cell is covered with sea ice, growth rate is set to 0
# 5. as a consequence, growth is delayed in areas that are still ice-covered
# load model functions
source('functions/RyabovEtAl2017AuxiliaryFunctions.R')

# start day
startDay <- 306

# load environmental data
templateGrid <- subset(rast('inputData/climatologies/photoperiodClim365.tif'),1)
chlorophyllClimatology <- rast('inputData/climatologies/climaChlaRollingAv15.nc')
chlorophyllClimatology <- flip(chlorophyllClimatology, direction = 'vertical')
crs(chlorophyllClimatology) <- 'epsg:4258'

# initialize chlorophyll a
initClimatologyChla <- subset(chlorophyllClimatology, startDay)

# initialize ice algae (zero)
initIceAlgae <- templateGrid
initIceAlgae[!is.na(initIceAlgae)] <- 0

# initialize time raster
timeRast <- templateGrid
timeRast[!is.na(timeRast)] <- 1

# initialize state variables (start length 26mm)
initLength <- 26
initLengthRast <- templateGrid
initLengthRast[!is.na(initLengthRast)] <- initLength

#initial age is 200 days
initAge <- 200
ageRast <- templateGrid
ageRast[!is.na(timeRast)] <- initAge

# initialize simulation environmental conditions
simulationInit <- c(initLengthRast,
                    timeRast,
                    initClimatologyChla,
                    initIceAlgae,
                    ageRast)

# create raster containing the simulation results for the state variable
resultsRast <- initLengthRast

for(i in 1:165){
  
  # determine current day of year for subsetting climatologies at correct position
  dayOfYear <- dayOfYearFunc(t = i, startDay = 306)
  
  # in first timestep, use initial conditions to calculate growth
  if(i == 1){
    
    # add change in length to previous state
    newLength <- subset(resultsRast,i) + lapp(simulationInit, rk4stepRyabov)
    
    # update chlorophyll a
    newChla <- subset(chlorophyllClimatology, dayOfYear)
    
    # update time- and age rast
    timeRast[!is.na(timeRast)] <- dayOfYear + 1
    ageRast[!is.na(ageRast)] <- initAge + i
    
    # create updated environment/state grid for following timestep
    newEnv <- c(newLength,
                timeRast,
                newChla,
                initIceAlgae,
                ageRast)
  }
  else{
    # add change in length to previous state
    newLength <- subset(resultsRast,i) + lapp(newEnv, rk4stepRyabov)
    
    # update chlorophyll a
    newChla <- subset(chlorophyllClimatology, dayOfYear)
    
    # update time- and age rast
    timeRast[!is.na(timeRast)] <- dayOfYear + 1
    ageRast[!is.na(ageRast)] <- initAge + i
    
    # create updated environment/state grid for following timestep
    newEnv <- c(newLength,
                timeRast,
                newChla,
                initIceAlgae,
                ageRast)
  }
  # add current time step to results raster
  add(resultsRast) <- newLength
  
  # when body size becomes negative flag cell with -999
  newLength[newLength < 0] <- -999
  print(i)
}

# save results
writeCDF(resultsRast, 
         varname = 'Ryabov_length_mm',
         overwrite=TRUE,
         'simulationResults/ryabovEtAlSummerGrowthdoy306_doy105.nc')

