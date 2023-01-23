# Calculate mean Southern Polar Front location based on Freeman and Lovenduski (2014)
library(ncdf4)
library(tidyterra)
library(terra)
library(sf)
library(tidyverse)

spf <- nc_open('~/github/krillGrowthModels/inputData/SouthernPolarFront/FreemanLovenduski_Polar_Front_weekly.nc')
print(spf)

coastline <- read_sf('/Users/dominik/Downloads/add_coastline_medium_res_polygon_v7_4/add_coastline_medium_res_polygon_v7_4.shp')

spfLat <- ncvar_get(spf, 'PFw')
spfLons <- ncvar_get(spf, 'longitude')
spfLons <- ifelse(spfLons > 180, spfLons-360, spfLons)
spfTimes <- ncvar_get(spf, 'time')

spfTib <- expand_grid(week = spfTimes,
                      lon = spfLons) %>% 
  mutate(lat = c(spfLat),
         part = 1) %>% 
  select(week, part, lon, lat) %>% 
  arrange(week, lon)

NAWeeks <- unique(spfTib$week[which(is.na(rowSums(spfTib)))])
spfMat <- spfTib %>% 
  filter(!week %in% NAWeeks) %>% 
  as.matrix(.)

spfVec <- vect(spfMat, 'polygons', crs='EPSG:4326')
spfVec <- terra::project(spfVec, '+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs')

# calculate mean SPF-position
meanSPF <- spfTib %>% 
  group_by(lon) %>% 
  summarize(lat = mean(lat)) %>% 
  mutate(part = 1, ID = 1) %>% 
  select(ID, part, lon, lat) %>% 
  as.matrix(.)

meanSPFVec <- vect(meanSPF, 'polygons', crs='EPSG:4326')
meanSPFVec <- terra::project(meanSPFVec, '+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs')

spfPlot <- ggplot() +
  geom_spatvector(data = spfVec, size = 0.25, colour = '#31313110', fill = NA) +
  geom_spatvector(data = meanSPFVec, size = 1.25, colour = '#ff9200', fill = NA) +
  geom_sf(data = coastline, fill = NA) +
  labs(title = 'Location of the Polar Front \n after Freeman and Lovenduski (2016)') +
  theme(panel.background = element_rect(fill = NA, colour = '#313131'),
        panel.grid.major = element_line(colour = '#99999960', linetype = '22'),
        legend.position = 'bottom',
        axis.text = element_text(size = 25, colour = '#313131'),
        plot.title = element_text(hjust = 0.5, size = 35, colour = '#313131'))
ggsave('~/github/krillGrowthModels/plots/spfLocations.jpg', plot = spfPlot, width = 10, height = 10)

# export mean SFP
writeVector(meanSPFVec, filename = 'inputData/SouthernPolarFront/meanSouthernPolarFront2002_2014.shp', overwrite = T)



