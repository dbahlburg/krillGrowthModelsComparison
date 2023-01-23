# July 18th, 2022
# Dominik Bahlburg
# In this script, I calculate the annual dynamics of photoperiod for the Southern Ocean South of 52 deg S. 
# I use a 1/4 degree grid resolution similar to the PISCES-NEMO output. 
# -------------------------------------------------------------------------------------------------------- #
# load required libraries
library(tidyverse)
library(lubridate)
library(suncalc)
library(terra)
# -------------------------------------------------------------------------------------------------------- #
# load pisces-nemo file to have grid template
piscesGrid <- subset(rast('/Volumes/T7/PISCES-NEMO-Climatology/climaChla_36.nc'),1)
piscesGrid <- flip(piscesGrid, direction = 'vertical')
crs(piscesGrid) <- 'epsg:4258'
piscesGridTb <- as.data.frame(piscesGrid, xy = T) %>% 
  dplyr::select(lat = y, lon = x)

photoStack <- rast()

for (i in 1:365){
  
    dayOfYear <- i
    photoperiod <- piscesGridTb %>% 
      distinct(lat) %>% 
      mutate(date = as.Date(dayOfYear, origin = "2021-12-31"),
             lon = 0) %>% 
      getSunlightTimes(data = ., keep = c('sunrise','sunset')) %>% 
      mutate(photoperiod = as.numeric(difftime(sunset, sunrise, unit = 'hours')),
             photoperiod = ifelse(between(dayOfYear, 110, 250) & is.na(photoperiod), 0, 
                                  ifelse(!between(dayOfYear, 110, 250) & is.na(photoperiod), 24, photoperiod))) %>% 
      dplyr::select(lat, photoperiod)
    photoperiodTib <- piscesGridTb %>% 
      left_join(., photoperiod, by = 'lat') %>% 
      dplyr::select(lon, lat, photoperiod)
    photoperiodRas <- rast(photoperiodTib)
    crs(photoperiodRas) <- 'epsg:4258'
    photoperiodRas <- project(photoperiodRas, piscesGrid)
    add(photoStack) <- photoperiodRas
    
    if (i %% 10 == 0){
      print(i)
    }
    if (i %% 20 == 0){
      plot(photoperiodRas)
    }
    if (i %% 365 == 0){
      #writeCDF(photoStack, '~/Desktop/photoperiod.nc', overwrite = T)
      writeRaster(photoStack, '/Volumes/T7/photoperiod/photoperiod365.tif', overwrite = T)
    }
}





