library(terra)
library(tidyverse)
library(lubridate)

# Data Sources:
# Global Ocean OSTIA Sea Surface Temperature and Sea Ice Reprocessed
# SST_GLO_SST_L4_REP_OBSERVATIONS_010_011 
# https://resources.marine.copernicus.eu/product-detail/SST_GLO_SST_L4_REP_OBSERVATIONS_010_011/INFORMATION

# create overview tibble containing the date info, day of year and indizes
sstFiles <- tibble(filepath = list.files('/Volumes/T7/SST_Data', pattern = '.nc', full.names = T),
                         filename = list.files('/Volumes/T7/SST_Data', pattern = '.nc')) %>%
  mutate(strDates = str_split(filename, pattern = '_')) %>% 
  rowwise() %>% 
  mutate(start = gsub("[^0-9]","", unlist(strDates)[1]),
         end = gsub("[^0-9]","", unlist(strDates)[2]),
         start = as.Date(paste(substr(start, 1, 4),substr(start, 5, 6),substr(start, 7, 8), sep = '-')),
         end = as.Date(paste(substr(end, 1, 4),substr(end, 5, 6),substr(end, 7, 8), sep = '-'))) %>% 
  mutate(intervalLength = end - start + 1) %>% 
  ungroup() %>% 
  slice(rep(1:n(), times = intervalLength)) %>% 
  group_by(start) %>% 
  mutate(date = seq(min(start), max(end), by = '1 day'),
         monthDay = paste(month(date), day(date), sep = '-')) %>% 
  filter(monthDay != '2-29') %>% 
  mutate(doy = as.Date(paste('1900', month(date), day(date), sep = '-')),
         sstIndex = 1:n()) %>% 
  mutate(doy = yday(doy)) %>% 
  ungroup() %>% 
  dplyr::select(-strDates, -monthDay, -filename, -start, -end, -intervalLength)

# --------------------------------------------------------------------------------- #
# Analysis - calculate climatology of SST and Sea Ice for each day of the year
climatologySST <- rast()
climatologySdSST <- rast()
logbook <- tibble() 

chlaRast <- subset(rast('/Volumes/T7/OceanColourL4-Output/climatologies/climaChlaRollingAv15.nc'),1)
for (i in 1:365){
  
  # filter files that contain data for respective day of year
  currentFiles <- sstFiles %>% 
    filter(doy == i)
  
  sstStack <- rast()
  
  # import the files
  for (k in 1:nrow(currentFiles)){
    sst <- aggregate(subset(rast(currentFiles$filepath[k]), currentFiles$sstIndex[k]),5)
    add(sstStack) <- sst
  }
  
  meanSST <- app(sstStack, mean)
  sdSST <- app(sstStack, sd)
  
  add(climatologySST) <- meanSST
  add(climatologySdSST) <- sdSST
  
  #create "logbook"
  logbook <- logbook %>% 
    bind_rows(., tibble(doy = i, 
                        date = currentFiles$date[1],
                        nFiles = nrow(currentFiles)))
  
  print(i)
  
  if(i%%365 == 0){
    climatologySST <- resample(climatologySST,chlaRast)
    climatologySdSST <- resample(climatologySdSST,chlaRast)
    terra::writeCDF(climatologySST, 
                    varname = 'climatologySST',
                    overwrite=TRUE,
                    filename = paste('/Volumes/T7/SST_Data/climaSST1_',i,'.nc',sep = ''))
    terra::writeCDF(climatologySdSST, 
                    varname = 'climatologySdSST',
                    overwrite=TRUE,
                    filename = paste('/Volumes/T7/SST_Data/climaSdSST1_',i,'.nc',sep = ''))
  }
}





