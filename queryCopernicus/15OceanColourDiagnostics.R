library(tidyverse)
library(lubridate)
library(terra)
library(tidyterra)
library(sf)

# ======================================================================================================================== #
# 1. Downsample original files to 0.25x0.25 degree grid - original files were deleted in the meantime to save disk space
# Data source: https://resources.marine.copernicus.eu/product-detail/OCEANCOLOUR_GLO_BGC_L4_MY_009_104/INFORMATION
# ======================================================================================================================== #
# create overview tibble containing the date info, day of year and indizes
# oceanColourFiles <- tibble(filepath = list.files('/Volumes/T7/OceanColourL4-Output', pattern = '.nc', full.names = T),
#                            filename = list.files('/Volumes/T7/OceanColourL4-Output', pattern = '.nc')) %>%
#   mutate(strDates = str_split(filename, pattern = '_')) %>% 
#   rowwise() %>% 
#   mutate(start = gsub("[^0-9]","", unlist(strDates)[2]),
#          end = gsub("[^0-9]","", unlist(strDates)[3]),
#          start = as.Date(paste(substr(start, 1, 4),substr(start, 5, 6),substr(start, 7, 8), sep = '-')),
#          end = as.Date(paste(substr(end, 1, 4),substr(end, 5, 6),substr(end, 7, 8), sep = '-'))) %>% 
#   mutate(intervalLength = end - start + 1) %>% 
#   ungroup() %>% 
#   slice(rep(1:n(), times = intervalLength)) %>% 
#   group_by(start) %>% 
#   mutate(date = seq(min(start), max(end), by = '1 day'),
#          monthDay = paste(month(date), day(date), sep = '-')) %>% 
#   filter(monthDay != '2-29') %>% 
#   mutate(doy = as.Date(paste('1900', month(date), day(date), sep = '-')),
#          doy = yday(doy),
#          chlaIndex = 1:n()) %>% 
#   ungroup() %>% 
#   dplyr::select(-strDates, -monthDay, -filename, -start, -end, -intervalLength) 
# 
# 
# piscesNemoGrid <- subset(rast('/Volumes/T7/PISCES-NEMO-Output/pisces_nemo2000-03-01_2000-03-31.nc'),1)
# oceanColourStack <- rast()
# oldDate <- oceanColourFiles$date[1]
# for (i in 1:nrow(oceanColourFiles)){
#   
#   oceanColour <- subset(rast(oceanColourFiles$filepath[i]), oceanColourFiles$chlaIndex[i])
#   oceanColour <- project(oceanColour, piscesNemoGrid)
#   add(oceanColourStack) <- oceanColour
#   
#   if(i%%25 == 0){
#     print(i)
#   }
#   
#   if(oceanColourFiles$doy[i] == 365){
#     writeCDF(oceanColourStack, paste('/Volumes/T7/OceanColourL4-Output/OCL4Resampled_',oldDate,'_',oceanColourFiles$date[i],'.nc', sep = ''))
#     oceanColourStack <- rast()
#     oldDate <- oceanColourFiles$date[i+1]
#   }
# }
# ======================================================================================================================== #
# import the resampled data
oceanColourFiles <- tibble(filepath = list.files('/Volumes/T7/OceanColourL4-Output', pattern = '.nc', full.names = T),
                           filename = list.files('/Volumes/T7/OceanColourL4-Output', pattern = '.nc')) %>%
  mutate(strDates = str_split(filename, pattern = '_')) %>%
  rowwise() %>%
  mutate(start = gsub("[^0-9]","", unlist(strDates)[2]),
         end = gsub("[^0-9]","", unlist(strDates)[3]),
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
         doy = yday(doy),
         chlaIndex = 1:n()) %>%
  ungroup() %>%
  dplyr::select(-strDates, -monthDay, -filename, -start, -end, -intervalLength)

# Calculate climatologies for each doy
climatologyMeanChla <- rast()
climatologySdChla <- rast()
climatologySampleSizes <- rast()
logbook <- tibble() 
oldDoy <- 1

# function that counts non-NAs when averaging
countFunc <- function(x){
  sum(!is.na(x))
}

for (i in 1:365){
  
  # filter files that contain data for respective day of year
  currentFiles <- oceanColourFiles %>% 
    filter(doy == i)
  
  chlaStack <- rast()
  
  # import the files
  for (k in 1:nrow(currentFiles)){
    chla <- subset(rast(currentFiles$filepath[k]), currentFiles$chlaIndex[k])
    add(chlaStack) <- chla
  }
  
  meanChla <- app(chlaStack, mean, na.rm = T)
  sdChla <- app(chlaStack, sd, na.rm = T)
  nChla <- app(chlaStack, countFunc)
  
  add(climatologyMeanChla) <- meanChla
  add(climatologySdChla) <- sdChla
  add(climatologySampleSizes) <- nChla
  
  #create "logbook"
  logbook <- logbook %>% 
    bind_rows(., tibble(doy = i, 
                        date = currentFiles$date[1],
                        nFiles = nrow(currentFiles)))
  
  print(i)
  
  if(i==365){
    terra::writeCDF(climatologyMeanChla, 
                    varname = 'climatologyChla',
                    overwrite=TRUE,
                    filename = paste('/Volumes/T7/OceanColourL4-Output/climaChla_',oldDoy,'_',i,'.nc',sep = ''))
    terra::writeCDF(climatologySdChla, 
                    varname = 'climatologySdChla',
                    overwrite=TRUE,
                    filename = paste('/Volumes/T7/OceanColourL4-Output/climaSdChla_',oldDoy,'_',i,'.nc',sep = ''))
    terra::writeCDF(climatologySampleSizes, 
                    varname = 'climatologySdChla',
                    overwrite=TRUE,
                    filename = paste('/Volumes/T7/OceanColourL4-Output/climatologySampleSizes_',oldDoy,'_',i,'.nc',sep = ''))
    
    climatologyMeanChla <- rast()
    climatologySdChla <- rast()
    oldDoy <- i+1
  }
}

chlaClim <- rast('/Volumes/T7/OceanColourL4-Output/climatologies/climaChla_doy_1_365.nc')
myCrs <- '+proj=stere +lat_0=-90 +lon_0=-15 +k=1 +x_0=-15 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0'
crs(chlaClim) <- 'epsg:4258'
crs(climatologySampleSizes) <- 'epsg:4258'

# Find day of year where satellite coverage of the Southern Ocean vanishes
# It seems that the timepoint is hard to determine in spring due to the vast sea ice cover of the Southern Ocean.
# It is not clear whether missing data come from missing satellite observations or ice cover.
# 
# It seems like day 240 could be a decent approximation:
chlaClimSub <- project(subset(chlaClim, 240),myCrs,gdal=F)
climatologySampleSizesSub <- project(flip(subset(climatologySampleSizes,240),direction='vertical'), myCrs, gdal = F)

# Import Sea ice climatology to have side-by-side plots
seaIceClim <- rast('/Volumes/T7/climatologySSTSeaIce/climaSeaIce_252.nc')
seaIceClim240 <- subset(seaIceClim, 24)
seaIceClim240 <- project(seaIceClim240, chlaClimSub)
coastline <- read_sf('/Users/dominik/Downloads/add_coastline_medium_res_polygon_v7_4/add_coastline_medium_res_polygon_v7_4.shp')

p1 <- ggplot() +
  geom_spatraster(data = seaIceClim240) +
  geom_sf(data = coastline, colour = '#c4c4c4', fill = '#4a4a4a', size = 0.08) +
  scale_fill_gradient(low = '#4f4f4f', high = '#49b4e6', 
                      na.value = NA, limits = c(0,1), breaks = c(0,0.5,1), labels = c(0,0.5,1)) +
  labs(fill = 'ice cover', title = 'sea ice climatology Aug 28th (doy 240)')   +
  theme(panel.background = element_rect(fill = NA, colour = '#303030'),
        panel.grid.major = element_blank(),
        plot.background = element_rect(fill = '#2b2b2b', colour = NA),
        legend.position = 'bottom',
        plot.subtitle = element_text(size = 15, hjust = 0.5, colour = '#f5f5f5'),
        legend.title = element_text(vjust = 1, size = 20, colour = '#f5f5f5'),
        legend.background = element_rect(fill = '#2b2b2b'),
        legend.text = element_text(size = 20, colour = '#f5f5f5'),
        axis.text = element_text(size = 25, colour = '#f5f5f5'),
        plot.title = element_text(hjust = 0.5, size = 35, colour = '#f5f5f5'))
ggsave('~/Desktop/seaIce240.png', plot = p1, width = 10, height = 10)


p2 <- ggplot() +
  geom_spatraster(data = chlaClimSub) +
  geom_sf(data = coastline, colour = '#c4c4c4', fill = '#4a4a4a', size = 0.08) +
  scale_fill_gradientn(colours = c('#2f3851','#6c8e69', '#ebb54b','#ff7b00','#ff2a00'),
                       na.value = NA, limits = c(0,3), breaks = c(0,3),
                       labels = c(0,3), oob = scales::squish) +
  labs(fill = expression(chlorophyll~a~mg/m^3), 
       title = 'chlorophyl a climatology Aug 28th (doy 240)')   +
  theme(panel.background = element_rect(fill = NA, colour = '#303030'),
        panel.grid.major = element_blank(),
        plot.background = element_rect(fill = '#2b2b2b', colour = NA),
        legend.position = 'bottom',
        plot.subtitle = element_text(size = 15, hjust = 0.5, colour = '#f5f5f5'),
        legend.title = element_text(vjust = 1, size = 20, colour = '#f5f5f5'),
        legend.background = element_rect(fill = '#2b2b2b'),
        legend.text = element_text(size = 20, colour = '#f5f5f5'),
        axis.text = element_text(size = 25, colour = '#f5f5f5'),
        plot.title = element_text(hjust = 0.5, size = 35, colour = '#f5f5f5'))
ggsave('~/Desktop/chlaClimatology240.png', plot = p2, width = 10, height = 10)



p3 <- ggplot() +
  geom_spatraster(data = climatologySampleSizesSub) +
  geom_sf(data = coastline, colour = '#c4c4c4', fill = '#4a4a4a', size = 0.08) +
  scale_fill_gradient(low = '#f2973d', high = '#303030',
                      na.value = NA, limits = c(0,18), breaks = c(0,18),
                      labels = c(0,18), oob = scales::squish) +
  labs(fill = 'number of samples', 
       title = 'climatology sample sizes Aug 28th (doy 240)')   +
  theme(panel.background = element_rect(fill = NA, colour = '#303030'),
        panel.grid.major = element_blank(),
        plot.background = element_rect(fill = '#2b2b2b', colour = NA),
        legend.position = 'bottom',
        plot.subtitle = element_text(size = 15, hjust = 0.5, colour = '#f5f5f5'),
        legend.title = element_text(vjust = 1, size = 20, colour = '#f5f5f5'),
        legend.background = element_rect(fill = '#2b2b2b'),
        legend.text = element_text(size = 20, colour = '#f5f5f5'),
        axis.text = element_text(size = 25, colour = '#f5f5f5'),
        plot.title = element_text(hjust = 0.5, size = 35, colour = '#f5f5f5'))
ggsave('~/Desktop/sampleSizes240.png', plot = p3, width = 10, height = 10)

# ======================================================================================================================== #
# Now find last available day
chlaClimSub <- project(subset(chlaClim, 100),myCrs,gdal=F)
climatologySampleSizesSub <- project(flip(subset(climatologySampleSizes,100),direction='vertical'), myCrs, gdal = F)

# Import Sea ice climatology to have side-by-side plots
seaIceClim <- rast('/Volumes/T7/climatologySSTSeaIce/climaSeaIce_108.nc')
seaIceClim100 <- subset(seaIceClim, 28)
seaIceClim100 <- project(seaIceClim100, chlaClimSub)
coastline <- read_sf('/Users/dominik/Downloads/add_coastline_medium_res_polygon_v7_4/add_coastline_medium_res_polygon_v7_4.shp')

p1 <- ggplot() +
  geom_spatraster(data = seaIceClim100) +
  geom_sf(data = coastline, colour = '#c4c4c4', fill = '#4a4a4a', size = 0.08) +
  scale_fill_gradient(low = '#4f4f4f', high = '#49b4e6', 
                      na.value = NA, limits = c(0,1), breaks = c(0,0.5,1), labels = c(0,0.5,1)) +
  labs(fill = 'ice cover', title = 'sea ice climatology Apr 10th (doy 100)')   +
  theme(panel.background = element_rect(fill = NA, colour = '#303030'),
        panel.grid.major = element_blank(),
        plot.background = element_rect(fill = '#2b2b2b', colour = NA),
        legend.position = 'bottom',
        plot.subtitle = element_text(size = 15, hjust = 0.5, colour = '#f5f5f5'),
        legend.title = element_text(vjust = 1, size = 20, colour = '#f5f5f5'),
        legend.background = element_rect(fill = '#2b2b2b'),
        legend.text = element_text(size = 20, colour = '#f5f5f5'),
        axis.text = element_text(size = 25, colour = '#f5f5f5'),
        plot.title = element_text(hjust = 0.5, size = 35, colour = '#f5f5f5'))
ggsave('~/Desktop/seaIce100.png', plot = p1, width = 10, height = 10)


p2 <- ggplot() +
  geom_spatraster(data = chlaClimSub) +
  geom_sf(data = coastline, colour = '#c4c4c4', fill = '#4a4a4a', size = 0.08) +
  scale_fill_gradientn(colours = c('#2f3851','#6c8e69', '#ebb54b','#ff7b00','#ff2a00'),
                       na.value = NA, limits = c(0,5), breaks = c(0,5),
                       labels = c(0,5), oob = scales::squish) +
  labs(fill = expression(chlorophyll~a~mg/m^3), 
       title = 'chlorophyl a climatology Apr 10th (doy 100)')   +
  theme(panel.background = element_rect(fill = NA, colour = '#303030'),
        panel.grid.major = element_blank(),
        plot.background = element_rect(fill = '#2b2b2b', colour = NA),
        legend.position = 'bottom',
        plot.subtitle = element_text(size = 15, hjust = 0.5, colour = '#f5f5f5'),
        legend.title = element_text(vjust = 1, size = 20, colour = '#f5f5f5'),
        legend.background = element_rect(fill = '#2b2b2b'),
        legend.text = element_text(size = 20, colour = '#f5f5f5'),
        axis.text = element_text(size = 25, colour = '#f5f5f5'),
        plot.title = element_text(hjust = 0.5, size = 35, colour = '#f5f5f5'))
ggsave('~/Desktop/chlaClimatology100.png', plot = p2, width = 10, height = 10)



p3 <- ggplot() +
  geom_spatraster(data = climatologySampleSizesSub) +
  geom_sf(data = coastline, colour = '#c4c4c4', fill = '#4a4a4a', size = 0.08) +
  scale_fill_gradient(low = '#f2973d', high = '#303030',
                      na.value = NA, limits = c(0,18), breaks = c(0,18),
                      labels = c(0,18), oob = scales::squish) +
  labs(fill = 'number of samples', 
       title = 'climatology sample sizes Apr 10th (doy 100)')   +
  theme(panel.background = element_rect(fill = NA, colour = '#303030'),
        panel.grid.major = element_blank(),
        plot.background = element_rect(fill = '#2b2b2b', colour = NA),
        legend.position = 'bottom',
        plot.subtitle = element_text(size = 15, hjust = 0.5, colour = '#f5f5f5'),
        legend.title = element_text(vjust = 1, size = 20, colour = '#f5f5f5'),
        legend.background = element_rect(fill = '#2b2b2b'),
        legend.text = element_text(size = 20, colour = '#f5f5f5'),
        axis.text = element_text(size = 25, colour = '#f5f5f5'),
        plot.title = element_text(hjust = 0.5, size = 35, colour = '#f5f5f5'))
ggsave('~/Desktop/sampleSizes100.png', plot = p3, width = 10, height = 10)

# Does deriving climatologies in 1 week intervals significantly improve things?
oceanColourFiles <- oceanColourFiles %>% 
  mutate(week = week(date))

# Calculate climatologies for each doy
climatologyMeanChla <- rast()
climatologySdChla <- rast()
climatologySampleSizes <- rast()
logbook <- tibble() 
oldWeek <- 1

# function that counts non-NAs when averaging
countFunc <- function(x){
  sum(!is.na(x))
}

for (i in 1:max(oceanColourFiles$week)){
  
  # filter files that contain data for respective day of year
  currentFiles <- oceanColourFiles %>% 
    filter(week == i)
  
  chlaStack <- rast()
  
  # import the files
  for (k in 1:nrow(currentFiles)){
    chla <- subset(rast(currentFiles$filepath[k]), currentFiles$chlaIndex[k])
    add(chlaStack) <- chla
  }
  
  meanChla <- app(chlaStack, mean, na.rm = T)
  sdChla <- app(chlaStack, sd, na.rm = T)
  nChla <- app(chlaStack, countFunc)
  
  add(climatologyMeanChla) <- meanChla
  add(climatologySdChla) <- sdChla
  add(climatologySampleSizes) <- nChla
  
  #create "logbook"
  logbook <- logbook %>% 
    bind_rows(., tibble(doy = i, 
                        date = currentFiles$date[1],
                        nFiles = nrow(currentFiles)))
  
  print(i)
  
  if(i==max(oceanColourFiles$week)){
    terra::writeCDF(climatologyMeanChla, 
                    varname = 'climatologyChla',
                    overwrite=TRUE,
                    filename = paste('/Volumes/T7/OceanColourL4-Output/climaChla_w',oldWeek,'_w',i,'.nc',sep = ''))
    terra::writeCDF(climatologySdChla, 
                    varname = 'climatologySdChla',
                    overwrite=TRUE,
                    filename = paste('/Volumes/T7/OceanColourL4-Output/climaSdChla_w',oldWeek,'_w',i,'.nc',sep = ''))
    terra::writeCDF(climatologySampleSizes, 
                    varname = 'climatologySdChla',
                    overwrite=TRUE,
                    filename = paste('/Volumes/T7/OceanColourL4-Output/climatologySampleSizes_w',oldWeek,'_w',i,'.nc',sep = ''))
    
    climatologyMeanChla <- rast()
    climatologySdChla <- rast()
    oldWeek <- i+1
  }
}

# ======================================================================================================================== #
# Plot the results and compare with 
chlaClimaWeekly <- rast('/Volumes/T7/OceanColourL4-Output/climatologies/climaChla_w1_w53.nc')
crs(chlaClimaWeekly) <- 'epsg:4258'
chlaClimaWeeklySub <- project(subset(chlaClimaWeekly, 15),myCrs,gdal=F)

p2 <- ggplot() +
  geom_spatraster(data = chlaClimaWeeklySub) +
  geom_sf(data = coastline, colour = '#c4c4c4', fill = '#4a4a4a', size = 0.08) +
  scale_fill_gradientn(colours = c('#2f3851','#6c8e69', '#ebb54b','#ff7b00','#ff2a00'),
                       na.value = NA, limits = c(0,5), breaks = c(0,5),
                       labels = c(0,5), oob = scales::squish) +
  labs(fill = expression(chlorophyll~a~mg/m^3), 
       title = 'chlorophyl a climatology week 15')   +
  theme(panel.background = element_rect(fill = NA, colour = '#303030'),
        panel.grid.major = element_blank(),
        plot.background = element_rect(fill = '#2b2b2b', colour = NA),
        legend.position = 'bottom',
        plot.subtitle = element_text(size = 15, hjust = 0.5, colour = '#f5f5f5'),
        legend.title = element_text(vjust = 1, size = 20, colour = '#f5f5f5'),
        legend.background = element_rect(fill = '#2b2b2b'),
        legend.text = element_text(size = 20, colour = '#f5f5f5'),
        axis.text = element_text(size = 25, colour = '#f5f5f5'),
        plot.title = element_text(hjust = 0.5, size = 35, colour = '#f5f5f5'))
ggsave('~/Desktop/chlaClimatologyWeek15.png', plot = p2, width = 10, height = 10)

climatologiesVec <- c(chlaClimaWeeklySub, chlaClimSub)
diffRast <- app(climatologiesVec,countFunc)
diffRast[diffRast!=1] <- NA

p3 <- ggplot() +
  geom_spatraster(data = diffRast) +
  scale_fill_gradient(low = '#36ff9e', high = '#36ff9e', na.value = NA) +
  geom_sf(data = coastline, colour = '#c4c4c4', fill = '#4a4a4a', size = 0.08) +
  labs(title = 'areas that could be filled with data')   +
  theme(panel.background = element_rect(fill = NA, colour = '#303030'),
        panel.grid.major = element_blank(),
        plot.background = element_rect(fill = '#2b2b2b', colour = NA),
        legend.position = 'none',
        plot.subtitle = element_text(size = 15, hjust = 0.5, colour = '#f5f5f5'),
        axis.text = element_text(size = 25, colour = '#f5f5f5'),
        plot.title = element_text(hjust = 0.5, size = 35, colour = '#f5f5f5'))
ggsave('~/Desktop/diffPlotWeek15.png', plot = p3, width = 10, height = 10)

# ------------------------------------------------------------------------------------ #
# how does this improve over time
missingDataTib <- tibble(date = seq(as.Date('2021-01-01'), as.Date('2021-12-31'), by = '1 day'),
                         doy = yday(date),
                         week = week(date))
diffTib <- NULL
lowResCoast <- vect('inputData/lowResCoastlineAntarctica/GSHHS_l_L5.shp')

chlaClimaWeekly[lowResCoast] <- -999
aggregatedDataGaps <- tibble(week = 1:53,
                             missingDataAgg = unlist(global(chlaClimaWeekly, fun = 'isNA')))

chlaClimCopy <- chlaClim
chlaClimCopy[lowResCoast] <- -999
dailyDataGaps <- tibble(doy = 1:365,
                        missingDataDaily = unlist(global(chlaClimCopy, fun = 'isNA')))

chlaClimCopy <- app(chlaClimCopy,  function(x){is.na(x) + 2})
chlaClimCopy[lowResCoast] <- -999
chlaClimCopy[chlaClimCopy == -999] <- NA
missingDataTib <- missingDataTib %>% 
  left_join(., aggregatedDataGaps, by = 'week') %>% 
  left_join(., dailyDataGaps, by = 'doy') %>% 
  mutate(totalData = unlist(global(chlaClimCopy, fun = 'notNA')))

p4 <- missingDataTib %>% 
  mutate(withoutAgg = missingDataDaily/totalData * 100,
         withAgg = missingDataAgg/totalData * 100) %>% 
  dplyr::select(doy, withoutAgg, withAgg) %>% 
  gather(type, value, -doy) %>% 
  ggplot(.,aes(x = doy, y = value, colour = type)) +
  geom_line(size = 1.5) +
  scale_colour_manual(values = c('#9c2424','#333333'), labels = c('aggregated data','daily data')) +
  scale_y_continuous(limits = c(0,100)) +
  labs(x = 'day of year',
       colour = '',
       y = 'data gaps as % of total grid cells',
       title = 'effect of aggregating weekly to close data gaps') +
  theme(panel.background = element_rect(fill = NA, colour = '#303030'),
        panel.grid.major = element_blank(),
        legend.position = c(0.18,0.85),
        legend.background = element_rect(fill = NA, colour = NA),
        plot.subtitle = element_text(size = 20, hjust = 0.5, colour = '#545454'),
        legend.title = element_text(vjust = 1, size = 24, colour = '#545454'),
        legend.key = element_rect(fill = NA),
        legend.text = element_text(size = 24, colour = '#545454'),
        axis.text = element_text(size = 25, colour = '#545454'),
        legend.key.width = unit(2, units = 'cm'),
        axis.title = element_text(size = 25, colour = '#545454'),
        plot.title = element_text(hjust = 0.5, size = 35, colour = '#545454')) +
  guides(linetype = guide_legend(override.aes = list(size = 5)))
ggsave('~/Desktop/Improvements3.png', width = 12, height = 8, plot = p4)

# ------------------------------------------------------------------------------------ #
# repeat the analysis but use monthly bins
oceanColourFiles <- oceanColourFiles %>% 
  mutate(month = month(date))

# Calculate climatologies for each doy
climatologyMeanChla <- rast()
climatologySdChla <- rast()
climatologySampleSizes <- rast()
logbook <- tibble() 
oldMonth <- 1

for (i in 1:max(oceanColourFiles$month)){
  
  # filter files that contain data for respective day of year
  currentFiles <- oceanColourFiles %>% 
    filter(month == i)
  
  chlaStack <- rast()
  
  # import the files
  for (k in 1:nrow(currentFiles)){
    chla <- subset(rast(currentFiles$filepath[k]), currentFiles$chlaIndex[k])
    add(chlaStack) <- chla
  }
  
  meanChla <- app(chlaStack, mean, na.rm = T)
  sdChla <- app(chlaStack, sd, na.rm = T)
  nChla <- app(chlaStack, countFunc)
  
  add(climatologyMeanChla) <- meanChla
  add(climatologySdChla) <- sdChla
  add(climatologySampleSizes) <- nChla
  
  #create "logbook"
  logbook <- logbook %>% 
    bind_rows(., tibble(doy = i, 
                        date = currentFiles$date[1],
                        nFiles = nrow(currentFiles)))
  
  print(i)
  
  if(i==max(oceanColourFiles$month)){
    terra::writeCDF(climatologyMeanChla, 
                    varname = 'climatologyChla',
                    overwrite=TRUE,
                    filename = paste('/Volumes/T7/OceanColourL4-Output/climaChla_m',oldMonth,'_m',i,'.nc',sep = ''))
    terra::writeCDF(climatologySdChla, 
                    varname = 'climatologySdChla',
                    overwrite=TRUE,
                    filename = paste('/Volumes/T7/OceanColourL4-Output/climaSdChla_m',oldMonth,'_m',i,'.nc',sep = ''))
    terra::writeCDF(climatologySampleSizes, 
                    varname = 'climatologySdChla',
                    overwrite=TRUE,
                    filename = paste('/Volumes/T7/OceanColourL4-Output/climatologySampleSizes_m',oldMonth,'_m',i,'.nc',sep = ''))
  }
}


# ======================================================================================================================== #
# Plot the results and compare with 
chlaClimaMonthly <- rast('/Volumes/T7/OceanColourL4-Output/climatologies/climaChla_m1_m12.nc')
crs(chlaClimaMonthly) <- 'epsg:4258'
chlaClimaMonthlySub <- project(subset(chlaClimaMonthly, 4),myCrs,gdal=F)
crs(chlaClim) <- 'epsg:4258'
chlaClimSub <- project(subset(chlaClim, 100),myCrs,gdal=F)

p2 <- ggplot() +
  geom_spatraster(data = chlaClimaMonthlySub) +
  geom_sf(data = coastline, colour = '#c4c4c4', fill = '#4a4a4a', size = 0.08) +
  scale_fill_gradientn(colours = c('#2f3851','#6c8e69', '#ebb54b','#ff7b00','#ff2a00'),
                       na.value = NA, limits = c(0,5), breaks = c(0,5),
                       labels = c(0,5), oob = scales::squish) +
  labs(fill = expression(chlorophyll~a~mg/m^3), 
       title = 'chlorophyl a climatology month 4')   +
  theme(panel.background = element_rect(fill = NA, colour = '#303030'),
        panel.grid.major = element_blank(),
        plot.background = element_rect(fill = '#2b2b2b', colour = NA),
        legend.position = 'bottom',
        plot.subtitle = element_text(size = 15, hjust = 0.5, colour = '#f5f5f5'),
        legend.title = element_text(vjust = 1, size = 20, colour = '#f5f5f5'),
        legend.background = element_rect(fill = '#2b2b2b'),
        legend.text = element_text(size = 20, colour = '#f5f5f5'),
        axis.text = element_text(size = 25, colour = '#f5f5f5'),
        plot.title = element_text(hjust = 0.5, size = 35, colour = '#f5f5f5'))
ggsave('~/Desktop/chlaClimatologyMonth4.png', plot = p2, width = 10, height = 10)

climatologiesVec <- c(chlaClimaMonthlySub, chlaClimSub)
diffRast <- app(climatologiesVec,countFunc)
diffRast[diffRast!=1] <- NA

p3 <- ggplot() +
  geom_spatraster(data = diffRast) +
  scale_fill_gradient(low = '#36ff9e', high = '#36ff9e', na.value = NA) +
  geom_sf(data = coastline, colour = '#c4c4c4', fill = '#4a4a4a', size = 0.08) +
  labs(title = 'areas that could be filled with data')   +
  theme(panel.background = element_rect(fill = NA, colour = '#303030'),
        panel.grid.major = element_blank(),
        plot.background = element_rect(fill = '#2b2b2b', colour = NA),
        legend.position = 'none',
        plot.subtitle = element_text(size = 15, hjust = 0.5, colour = '#f5f5f5'),
        axis.text = element_text(size = 25, colour = '#f5f5f5'),
        plot.title = element_text(hjust = 0.5, size = 35, colour = '#f5f5f5'))
ggsave('~/Desktop/diffPlotMonth4.png', plot = p3, width = 10, height = 10)

# ------------------------------------------------------------------------------------ #
# how does this improve over time
missingDataTib <- tibble(date = seq(as.Date('2021-01-01'), as.Date('2021-12-31'), by = '1 day'),
                         doy = yday(date),
                         month = month(date))
diffTib <- NULL

chlaClimaMonthly[lowResCoast] <- -999
aggregatedDataGaps <- tibble(month = 1:12,
                             missingDataAgg = unlist(global(chlaClimaMonthly, fun = 'isNA')))

chlaClimCopy <- chlaClim
chlaClimCopy[lowResCoast] <- -999
dailyDataGaps <- tibble(doy = 1:365,
                        missingDataDaily = unlist(global(chlaClimCopy, fun = 'isNA')))

chlaClimCopy <- app(chlaClimCopy,  function(x){is.na(x) + 2})
chlaClimCopy[lowResCoast] <- -999
chlaClimCopy[chlaClimCopy == -999] <- NA
missingDataTib <- missingDataTib %>% 
  left_join(., aggregatedDataGaps, by = 'month') %>% 
  left_join(., dailyDataGaps, by = 'doy') %>% 
  mutate(totalData = unlist(global(chlaClimCopy, fun = 'notNA')))

p4 <- missingDataTib %>% 
  mutate(withoutAgg = missingDataDaily/totalData * 100,
         withAgg = missingDataAgg/totalData * 100) %>% 
  dplyr::select(doy, withoutAgg, withAgg) %>% 
  gather(type, value, -doy) %>% 
  ggplot(.,aes(x = doy, y = value, colour = type)) +
  geom_line(size = 1.5) +
  scale_colour_manual(values = c('#9c2424','#333333'), labels = c('aggregated data','daily data')) +
  scale_y_continuous(limits = c(0,100)) +
  labs(x = 'day of year',
       colour = '',
       y = 'data gaps as % of total grid cells',
       title = 'effect of aggregating monthly to close data gaps') +
  theme(panel.background = element_rect(fill = NA, colour = '#303030'),
        panel.grid.major = element_blank(),
        legend.position = c(0.18,0.85),
        legend.background = element_rect(fill = NA, colour = NA),
        plot.subtitle = element_text(size = 20, hjust = 0.5, colour = '#545454'),
        legend.title = element_text(vjust = 1, size = 24, colour = '#545454'),
        legend.key = element_rect(fill = NA),
        legend.text = element_text(size = 24, colour = '#545454'),
        axis.text = element_text(size = 25, colour = '#545454'),
        legend.key.width = unit(2, units = 'cm'),
        axis.title = element_text(size = 25, colour = '#545454'),
        plot.title = element_text(hjust = 0.5, size = 35, colour = '#545454')) +
  guides(linetype = guide_legend(override.aes = list(size = 5)))
ggsave('~/Desktop/Improvements4.png', width = 12, height = 8, plot = p4)

# ======================================================================================================================== #
# Fill data gaps with Pisces-Nemo output on selected days (1,91,182,273)
chlaClimDaily <- subset(chlaClim, c(1,91,182,273))
chlaClimDaily[lowResCoast] <- -999

# import pisces-nemo climatologies
piscesNemoFiles <- list.files('/Volumes/T7/PISCES-NEMO-Climatology', pattern='climaChla', full.names = T)
piscesNemoClim <- rast(piscesNemoFiles)
piscesNemoClim <- flip(piscesNemoClim, direction = 'vertical')
ext(piscesNemoClim) <- ext(chlaClimDaily)
piscesNemoSub <- subset(piscesNemoClim, c(1,91,182,273))

# fill gaps in daily climatologies
mergedData <- rast()
for(i in 1:4){
  oceanColour <- subset(chlaClimDaily, i)
  piscesNemo <- subset(piscesNemoSub, i)
  oceanColour[is.na(oceanColour)] <- piscesNemo[is.na(oceanColour)]
  oceanColour[oceanColour == -999] <- NA
  add(mergedData) <- oceanColour
}

#transform coordinate system
chlaClimDaily[chlaClimDaily==-999] <- NA
chlaClimDaily <- project(chlaClimDaily,myCrs)
mergedData <- project(mergedData,myCrs)

# Plot data
plotFunc <- function(inputData, fileName, plotTitle, zRange){
  p2 <- ggplot() +
    geom_spatraster(data = inputData) +
    geom_sf(data = coastline, colour = '#c4c4c4', fill = '#4a4a4a', size = 0.08) +
    scale_fill_gradientn(colours = c('#2f3851','#6c8e69', '#ebb54b','#ff7b00','#ff2a00'),
                         na.value = NA, limits = zRange, breaks = zRange, labels = zRange) +
    labs(fill = expression(chlorophyll~a~mg/m^3), 
         title = plotTitle)   +
    theme(panel.background = element_rect(fill = NA, colour = '#303030'),
          panel.grid.major = element_blank(),
          plot.background = element_rect(fill = '#2b2b2b', colour = NA),
          legend.position = 'bottom',
          plot.subtitle = element_text(size = 15, hjust = 0.5, colour = '#f5f5f5'),
          legend.title = element_text(vjust = 1, size = 20, colour = '#f5f5f5'),
          legend.background = element_rect(fill = '#2b2b2b'),
          legend.text = element_text(size = 20, colour = '#f5f5f5'),
          axis.text = element_text(size = 25, colour = '#f5f5f5'),
          plot.title = element_text(hjust = 0.5, size = 35, colour = '#f5f5f5'))
  ggsave(paste('~/Desktop/chlaClimatology',fileName,'.png', sep = ''), plot = p2, width = 10, height = 10)
}

plotFunc(inputData = subset(chlaClimDaily,1), fileName = 'NotFilled1', plotTitle = 'Ocean Colour climatology doy 1', zRange = c(0,4))
plotFunc(inputData = subset(chlaClimDaily,2), fileName = 'NotFilled91',plotTitle = 'Ocean Colour climatology doy 91', zRange = c(0,4))
plotFunc(inputData = subset(chlaClimDaily,3), fileName = 'NotFilled182',plotTitle = 'Ocean Colour climatology doy 182', zRange = c(0,4))
plotFunc(inputData = subset(chlaClimDaily,4), fileName = 'NotFilled273',plotTitle = 'Ocean Colour climatology doy 273', zRange = c(0,4))

plotFunc(inputData = subset(mergedData,1), fileName = 'Filled1', plotTitle = 'Merged data doy 1', zRange = c(0,4))
plotFunc(inputData = subset(mergedData,2), fileName = 'Filled91',plotTitle = 'Merged data doy 91', zRange = c(0,4))
plotFunc(inputData = subset(mergedData,3), fileName = 'Filled182',plotTitle = 'Merged data doy 182', zRange = c(0,4))
plotFunc(inputData = subset(mergedData,4), fileName = 'Filled273',plotTitle = 'Merged data doy 273', zRange = c(0,4))






