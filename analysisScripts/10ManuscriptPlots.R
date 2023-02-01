# This script contain all (sub)plots that are part of the manuscript and Appendix,
# except for the DGR-plains in Figure 1 as well as the chlorophyll a sensitivity plots
# in the Appendix
# Dominik Bahlburg
# 23.01.2023
# --------------------------------------------------------------------------------------- #
# load libraries
library(terra)
library(tidyterra)
library(sf)
library(scico)
library(tidyverse)
# --------------------------------------------------------------------------------------- #
# load auxiliary functions
# source day of year function
dayOfYear <- source('functions/dayOfYear.R')
# --------------------------------------------------------------------------------------- #
# import environmental data, Southern Polar Front and simulation results
# import chlorophyll climatologies
chlorophyllClimatology <- rast('inputData/climatologies/climaChlaRollingAv15.nc')
chlorophyllClimatology <- flip(chlorophyllClimatology, direction = 'vertical')
crs(chlorophyllClimatology) <- 'epsg:4258'

# import sst climatologies
sstClimatology <- rast('inputData/climatologies/climatologySST1_365.nc')

# import Southern Polar Front, derived from Freeman and Lovenduski (2016)
meanSPF <- vect('inputData/SouthernPolarFront/meanSouthernPolarFront2002_2014.shp')

# list files with model run results (exclude Constable and Kawaguchi (2018))
summerGrowthFiles <- list.files('simulationResults', pattern = 'SummerGrowthdoy306_doy105.nc', recursive = T, full.names = T)
summerGrowthFiles <- summerGrowthFiles[-which(str_detect(summerGrowthFiles, 'ConstableKawaguchi'))]

# import the results of the different model runs and stack them
modelResults <- rast(summerGrowthFiles)

# create vector for properly labelling the plots
nameVec <- c('Atkinson et al. (2006)*','Bahlburg et al. (2021)',
             'Fach et al. (2002)', 'Hofmann & Lascara (2000)','Jager & Ravagnan (2015)','Ryabov et al. (2017)',
             'Tarling et al. (2006)*','Wiedenmann et al. (2008)*')

# number of models
nModels <- 8

# import coastline
coastline <- read_sf('inputData/add_coastline_medium_res_polygon_v7_4/add_coastline_medium_res_polygon_v7_4.shp')

# transform projection of simulation results to stereographic for visualization
myCrs <- '+proj=stere +lat_0=-90 +lon_0=-15 +k=1 +x_0=-15 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0'

# --------------------------------------------------------------------------------------- #
# Figure 2: model ensemble plot for South Orkney location and environmental dynamics
# extract growth trajectories for  -60.959390, -45.458629 
coords <- matrix(c(-45.458629, -60.959390), ncol = 2)
SOTraj <- terra::extract(modelResults, coords) %>% 
  gather(source, value) %>% 
  mutate(model = rep(nameVec, each = 166),
         time = rep(1:166, times = 8),
         date = rep(seq(as.Date('2010-11-02'), as.Date('2011-04-16'), '1 day'), times = 8))

# calculate mean growth trajectory
SOMean <- SOTraj %>% 
  group_by(date) %>% 
  summarize(meanLength = mean(value),
            sdLength = sd(value),
            label = 'meanLength') %>% 
  ungroup()

# tibble that is used to slightly offset the labels of some models to avoid overlap
modLabel <- SOTraj %>% 
  filter(date == max(date)) %>% 
  mutate(valueMod = ifelse(model %in% c('Ryabov et al. (2017)',
                                        'Hofmann & Lascara (2000)'),
                           value - 1.25,
                           ifelse(model %in% c('Wiedenmann et al. (2008)*',
                                               'Fach et al. (2002)'),
                                  value + 1.25,
                                  value)))

# create plot
p3 <- SOTraj %>% 
  ggplot(.) +
  geom_ribbon(data = SOMean, aes(x = date, ymin = meanLength-sdLength,
                                 ymax = meanLength+sdLength), fill = '#58a6b885') +
  geom_line(data = SOMean, aes(x = date, y = meanLength), size = 0.9) +
  geom_line(aes(x = date, y = value, group = model), size = 0.4, colour = '#6e6e6e',
            linetype = '42') +
  geom_point(data = modLabel, aes(x = date, y = value), size = 2.5, colour = '#6e6e6e', pch=21,
             fill = '#ffffff50', stroke = 0.75) +
  geom_text(data = modLabel, aes(label = model, x = date + 3, y = valueMod),
            colour = '#6e6e6e', hjust = 0, size = 6.5) +
  geom_text(data = filter(SOMean, date == max(date)),
            aes(x = date + 3,
                y = meanLength + 0.4), label = 'mean length', size = 6.5, hjust = 0) +
  scale_x_date(breaks = seq(as.Date('2010-11-01'),as.Date('2011-04-01'), by = '1 month'), date_labels = '%b',
               limits = c(as.Date('2010-11-01'),as.Date('2011-06-25'))) +
  labs(x = 'time', y = 'length in mm',
       subtitle = '60.96°S, 45.46°W: Growth trajectories') +
  theme(panel.background = element_rect(fill = NA, colour = '#303030'),
        panel.grid.major = element_blank(),
        legend.position = 'none',
        legend.background = element_rect(fill = NA, colour = NA),
        plot.subtitle = element_text(size = 33, hjust = 0.5, colour = '#545454'),
        legend.title = element_text(vjust = 1, size = 24, colour = '#545454'),
        legend.key = element_rect(fill = NA),
        legend.text = element_text(size = 30, colour = '#545454'),
        axis.text = element_text(size = 25, colour = '#545454'),
        legend.key.width = unit(2, units = 'cm'),
        axis.title.y = element_text(size = 25, colour = '#545454'),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 35, colour = '#545454'))
ggsave('plots/modEnsemble.pdf', plot = p3, width = 13, height = 8)


# add dynamics of temperature and chlorophyll a, extract at same location
chlaDyn <- terra::extract(chlorophyllClimatology, coords) 
tempDyn <- terra::extract(sstClimatology, coords)

chlaDynSummer <- unlist(chlaDyn[1, (dayOfYearFunc(1:165, startDay = 306))])
tempDynSummer <- unlist(tempDyn[1, (dayOfYearFunc(1:165, startDay = 306))]) - 273.15

envDynamics <- as.data.frame(cbind(chlaDynSummer, tempDynSummer)) %>% 
  mutate(date = seq(as.Date('2020-11-02'),as.Date('2021-04-15'), by = '1 day')) 

# create dual-axis-plot
dualPlot <- envDynamics %>% 
  ggplot(.,aes(x = date)) +
  geom_line(aes(y = chlaDynSummer), colour = '#4287ed', size = 1) +
  geom_line(aes(y = (tempDynSummer+3)/4), colour = '#ed4e42', size = 1) +
  scale_y_continuous(limits = c(0.2,1.2), breaks = c(0.25,0.75,1.25), sec.axis = sec_axis(~.*4-3, name = 'water temperature °C')) +
  scale_x_date(breaks = seq(as.Date('2020-11-01'),as.Date('2021-04-01'), by = '1 month'), date_labels = '%b',
               limits = c(as.Date('2020-11-01'),as.Date('2021-04-15'))) +
  labs(x = 'simulation day', y = expression(chlorophyll~a~mg~m^{-3}), title = 'environmental dynamics') +
  theme(panel.background = element_rect(fill = NA, colour = '#303030'),
        panel.grid.major = element_blank(),
        plot.title = element_text(size = 30, colour = '#545454', hjust = 0.5),
        axis.text = element_text(size = 30, colour = '#545454'),
        axis.text.y.left = element_text(size = 30, colour = '#4287ed'),
        axis.text.y.right = element_text(size = 30, colour = '#ed4e42'),
        axis.title.y.left = element_text(size = 30, colour = '#4287ed'),
        axis.title.y.right = element_text(size = 30, colour = '#ed4e42'),
        axis.title.x = element_blank(),
        plot.margin = margin(t = 1, r = 0, b = 1, l = 0, unit = "cm"))

ggsave('plots/envDynamics.pdf', plot = dualPlot, width = 12, height = 5)
# --------------------------------------------------------------------------------------- #
# Figure 3: change-in-length-maps
# calculate coefficient of variation sd/mean
modelResultsCV <- app(modelResultsSub, sd, na.rm = T)/app(modelResultsSub, mean, na.rm = T)
modelResultsStereoCV <- terra::project(modelResultsCV, myCrs)

# Calculating length differences
extractSimDay0 <- 1
extractSimDay <- 166
modelResultsDay0 <- subset(modelResults, seq(extractSimDay0, extractSimDay0 + (166 * (nModels-1)), length.out = nModels))
modelResultsDayMax <- subset(modelResults, seq(extractSimDay, extractSimDay + (166 * (nModels-1)), length.out = nModels))

modelResultsDiff <- modelResultsDayMax - modelResultsDay0
modelResultsDiff[modelResultsDiff < -26] <- NA

# change projection to stereographic
modelResultsDiffStereo <- terra::project(modelResultsDiff, myCrs)

# remove data North of the Southern Polar Front
modelResultsDiffStereo <- mask(modelResultsDiffStereo, meanSPF)

# create binning function for binning the results since I am later using a discrete
# colour scale
binningFunc <- function(x){
  binSize <- 5
  return(round(x/binSize) * binSize)
}

# Create model overview plots
for(i in 1:8){
  
  resultsLengthRastStereoNew <- subset(modelResultsDiffStereo,i)
  darkMode <- ggplot() +
    geom_spatraster(data = resultsLengthRastStereoNew) +
    geom_sf(data = coastline, colour = '#cccccc', fill = '#61616150', size = 0.5) +
    scale_fill_stepsn(colours = scico(n = 11, palette = 'vikO', begin = 10/45, end = 1),
                      na.value = '#303030',
                      limits = c(-15,30),
                      breaks = seq(-15,30,by = 5)) +
    labs(fill = expression(Delta~L),
         caption = nameVec[i]) +
    theme(panel.background = element_rect(fill = NA, colour = NA),
          panel.grid.major = element_blank(),
          plot.background = element_rect(fill = 'transparent', colour = NA),
          legend.position = 'none',
          axis.text = element_blank(),
          plot.caption = element_text(hjust = 0.5, size = 42, colour = '#f5f5f5')) +
    guides(fill = guide_colourbar(title.position = 'top', title.hjust = 0.5))
  
  ggsave(paste('plots/modDiff',i,'.pdf',sep = ''), plot = darkMode, width = 10, height = 10, bg = 'transparent') 
}

# create one more map to extract the legend
legendPlot <- ggplot() +
  geom_spatraster(data = resultsLengthRastStereoNew) +
  geom_sf(data = coastline, colour = '#cccccc', fill = '#61616150', size = 0.5) +
  scale_fill_stepsn(colours = scico(n = 11, palette = 'vikO', begin = 10/45, end = 1),
                    na.value = '#303030',
                    limits = c(-15,30),
                    breaks = seq(-15,30,by = 5),
                    labels = c('<-15','-10','-5','0','5','10','15','20','25','>30')) +
  labs(fill = expression(Delta~L~'in'~mm),#'body length in mm',
       caption = nameVec[i]) +
  theme(panel.background = element_rect(fill = NA, colour = NA),
        panel.grid.major = element_blank(),
        plot.background = element_rect(fill = NA, colour = NA),
        legend.position = 'top',
        plot.subtitle = element_text(size = 42, hjust = 0.5, colour = '#f5f5f5'),
        legend.title = element_text(vjust = 1, size = 30, colour = '#f5f5f5'),
        legend.background = element_rect(fill = NA),
        legend.key.width = unit(2.5,'cm'),
        legend.text = element_text(size = 21, colour = '#f5f5f5'),
        axis.text = element_blank(),
        plot.caption = element_text(hjust = 0.5, size = 40, colour = '#f5f5f5', vjust = 0)) +
  guides(fill = guide_colourbar(title.position = 'top', title.hjust = 0.5))
ggsave('plots/legendDiff.pdf', plot = legendPlot, width = 10, height = 10, bg = 'transparent') 

# create the plot of coefficient of variation
modelResultsStereoCV <- mask(modelResultsStereoCV, meanSPF)
p4 <- ggplot() +
  geom_spatraster(data = modelResultsStereoCV) +
  geom_sf(data = coastline, colour = '#3d3d3d', fill = '#3d3d3d', size = 0.1) +
  scale_fill_gradientn(#colours = c('#405a85', '#f7f7f7', '#c74a4a'),
    colours = c('#405a85', '#ffd470','#c74a4a'),
    na.value = NA,
    limits = c(0,0.75), breaks = c(0,0.75), labels = c('0','>0.75'), 
    oob = scales::squish) +
  labs(fill = 'Coefficient of Variation',
       caption = 'Coefficient of Variation')  +
  theme(panel.background = element_rect(fill = 'transparent', colour = NA),
        panel.grid.major = element_blank(),
        plot.background = element_rect(fill = 'transparent', colour = NA),
        legend.position = 'top',
        plot.subtitle = element_text(size = 42, hjust = 0.5, colour = '#f5f5f5'),
        legend.title = element_text(vjust = 1, size = 30, colour = '#f5f5f5'),
        legend.background = element_rect(fill = 'transparent'),
        legend.key.width = unit(2.5,'cm'),
        legend.text = element_text(size = 30, colour = '#f5f5f5'),
        #axis.text = element_text(size = 25, colour = '#f5f5f5'),
        axis.text = element_blank(),
        plot.caption = element_text(hjust = 0.5, size = 40, colour = '#f5f5f5')) +
  guides(fill = guide_colourbar(title.position = 'top', title.hjust = 0.5))
ggsave('plots/ModelIntercopmCVAp15.pdf', plot = p4, width = 10, height = 10, bg = 'transparent') 

# create plot for mean chlorophyll a climatology
# first calculate mean, change projection to stereographic
chlorophyllClimatology[is.na(chlorophyllClimatology)] <- 0
chlaMean <- app(subset(chlorophyllClimatology, dayOfYearFunc(1:165, startDay = 306)), mean)
chlaMean <- terra::project(chlaMean, myCrs)
chlaMean <- mask(chlaMean, meanSPF)

p2 <- ggplot() +
  geom_spatraster(data = chlaMean) +
  geom_sf(data = coastline, colour = NA, fill = '#616161') +
  scale_fill_gradientn(colours = c('#2f3851','#6c8e69', '#ebb54b','#ff7b00','#ff2a00'),
                       na.value = NA, limits = c(0,2), breaks = c(0,2),
                       labels = c(0,2), oob = scales::squish) +
  labs(fill = expression(Chlorophyll~a~mg~m^{-3}), 
       caption = 'Mean Chl a Climatology') +
  theme(panel.background = element_rect(fill = NA, colour = NA),
        panel.grid.major = element_blank(),
        plot.background = element_rect(fill = NA, colour = NA),
        legend.position = 'top',
        plot.subtitle = element_text(size = 42, hjust = 0.5, colour = '#f5f5f5'),
        legend.title = element_text(vjust = 1, size = 30, colour = '#f5f5f5'),
        legend.background = element_rect(fill = NA),
        legend.key.width = unit(2.5,'cm'),
        legend.text = element_text(size = 30, colour = '#f5f5f5'),
        #axis.text = element_text(size = 25, colour = '#f5f5f5'),
        axis.text = element_blank(),
        plot.caption = element_text(hjust = 0.5, size = 40, colour = '#f5f5f5', vjust = 0)) +
  guides(fill = guide_colourbar(title.position = 'top', title.hjust = 0.5))
ggsave('plots/oceanColourMean.pdf', plot = p2, width = 10, height = 10, bg = 'transparent') 

# now repeat for sea surface temperature
sstMean <- app(subset(sstClimatology, dayOfYearFunc(1:165, startDay = 306)), mean, na.rm = T) - 273.15
sstMean <- project(sstMean, myCrs)
sstMean <- mask(sstMean, meanSPF)

p2 <- ggplot() +
  geom_spatraster(data = sstMean) +
  geom_sf(data = coastline, colour = NA, fill = '#616161') +
  scale_fill_scico(palette = 'batlow', limits = c(-2.5, 5), breaks = c(-2.5, 5), labels = c('-2.5°C','5°C'),
                   na.value = NA) +
  labs(fill = 'Sea Surface Temperature', 
       caption = 'Mean SST Climatology') +
  theme(panel.background = element_rect(fill = NA, colour = NA),
        panel.grid.major = element_blank(),
        plot.background = element_rect(fill = NA, colour = NA),
        legend.position = 'top',
        plot.subtitle = element_text(size = 42, hjust = 0.5, colour = '#f5f5f5'),
        legend.title = element_text(vjust = 1, size = 30, colour = '#f5f5f5'),
        legend.background = element_rect(fill = NA),
        legend.key.width = unit(2.5,'cm'),
        legend.text = element_text(size = 30, colour = '#f5f5f5'),
        #axis.text = element_text(size = 25, colour = '#f5f5f5'),
        axis.text = element_blank(),
        plot.caption = element_text(hjust = 0.5, size = 40, colour = '#f5f5f5', vjust = 0)) +
  guides(fill = guide_colourbar(title.position = 'top', title.hjust = 0.5))
ggsave('plots/sstMean.pdf', plot = p2, width = 10, height = 10, bg = 'transparent') 

# --------------------------------------------------------------------------------------- #
# Create plots contained in Table 2
# assign seasonal cover to each model
seasonalCover <- tibble(Atkinson = list(c('spring','summer')),
                        Bahlburg = list(c('spring','summer','winter','autumn')),
                        Constable = list(c('spring','summer','winter','autumn')),
                        Hofmann = list(c('spring','summer','winter','autumn')),
                        Fach = list(c('spring','summer','winter','autumn')),
                        Jager = list(c('spring','summer','winter','autumn')),
                        Ryabov = list(c('spring','summer','winter','autumn')),
                        Tarling = list(c('spring','summer')),
                        Wiedenmann = list(c('spring','summer')))

# define rectangles for the different seasons
seasonsRect <- tibble(season = factor(c('winter','spring','summer','autumn', 'winter'),
                                      levels = c('spring','summer','autumn', 'winter')),
                      xmin = c(-0.125,0,0.25,0.5,0.75),
                      xmax = c(0, 0.25,0.5,0.75,0.875),
                      ymin = 0.85,
                      ymax = 1.4)

# season labels with position
seasonsLabel <- tibble(season = c('spring','summer','autumn', 'winter'),
                       x = c(0.125, 0.375, 0.625, 0.875),
                       y = 1.5)

# resolution, 
reso <- 0.0005
seasonsLabelPlot <- filter(seasonsLabel, season %in% unlist(seasonalCover[1,i]))
seasonsCovered <- filter(seasonsRect, season %in% unlist(seasonalCover[1,i])) 

# function that lets colour fade in/out on season edges
fadeIn <- 1/121 * 1:121 
fadeOut <- rev(fadeIn)
seasonsRectFade <- tibble(xmin = seq(min(seasonsCovered$xmin), 
                                     max(seasonsCovered$xmax - reso), by = reso),
                          xmax = seq(min(seasonsCovered$xmin) + reso, 
                                     max(seasonsCovered$xmax), by = reso)) %>% 
  mutate(ID = 1:n(),
         alpha = 1)
seasonsRectFade$alpha[1:121] <- fadeIn
seasonsRectFade$alpha[(nrow(seasonsRectFade) - 120):nrow(seasonsRectFade)] <- fadeOut

# create plot for summer-spring coverage
p1 <- ggplot() +
  geom_path(aes(x = c(-0.125,0.875), y = 1.4), size = 30, colour = '#f5f5f5') +
  geom_path(data = seasonsRectFade, aes(x = xmin, y = 1.4, alpha = alpha), size = 30, colour = '#d66b6d') +
  scale_y_continuous(limits = c(0,1.8)) +
  scale_x_continuous(limits = c(-0.125,0.875)) +
  geom_textline(data = seasonsLabel, aes(x = x, y = y, group = season, label = season),
                linewidth = 0, size = 22, vjust = 1.5, colour = '#2b2b2b') +
  coord_polar() +
  theme(panel.background = element_rect(fill = 'transparent', colour = NA),
        plot.background = element_rect(fill = 'transparent', colour = NA),
        legend.position = 'none',
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank())
ggsave('plots/seasonalCoverSummerSpring.pdf', plot = p1,
       width = 7, height = 7)

# create plot for full seasonal coverage
p1 <- ggplot() +
  geom_path(aes(x = c(-0.125,0.875), y = 1.4), size = 30, colour = '#d66b6d') +
  scale_y_continuous(limits = c(0,1.8)) +
  scale_x_continuous(limits = c(-0.125,0.875)) +
  geom_textline(data = seasonsLabel, aes(x = x, y = y, group = season, label = season),
                linewidth = 0, size = 22, vjust = 1.5, colour = '#2b2b2b') +
  coord_polar() +
  theme(panel.background = element_rect(fill = 'transparent', colour = NA),
        plot.background = element_rect(fill = 'transparent', colour = NA),
        legend.position = 'none',
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank())
ggsave('plots/seasonalCoverAllYear.pdf', plot = p1, 
       width = 7, height = 7)

# --------------------------------------------------------------------------------------- #
# Create chlorophyll a density-curve plots
# import original files
RyabovPocData <- read_csv("inputData/originalResults/RyabovPocData.csv") %>% 
  mutate(Y = 10^(Y)/50)

# built loess smoother to estimage environmental model
RyabovEnvMod <- RyabovPocData %>% 
  select(X,Y) %>% 
  nest(data = c(X, Y)) %>% 
  mutate(loessModel = purrr::map(data, loess,
                                 formula = Y ~ X, span = 0.03))

# predict chlorophyll a for each day of the 11 years as was done in the original paper
RyabovPocRecon <- tibble(day = 1:4015) %>%
  mutate(loessModel = RyabovEnvMod$loessModel,
         chlorophyllA = unlist(map2(loessModel, day, predict)))

# visualize density of chlorophyll a values
p1 <- RyabovPocRecon %>% 
  ggplot(.,aes(x = chlorophyllA)) +
  geom_density(colour = '#d66b6d', size = 6) +
  scale_x_continuous(limits = c(0,7), breaks = c(0,7), labels = c(0,7)) +
  scale_y_continuous(expand = c(0.1,0.1), breaks = c(0,1), labels = c(0,1)) +
  labs(x = expression(mg~m^{-3})) +
  theme(axis.text.y = element_text(size = 55, colour = '#141414', face = 'bold', vjust = 0),#element_blank(),
        axis.text.x = element_text(size = 55, colour = '#141414', face = 'bold', vjust = 0),
        axis.line.x = element_line(size = 3.5, colour = '#141414'),
        axis.line.y = element_line(size = 3.5, colour = '#141414'),
        axis.title.x = element_text(size = 55, colour = '#141414', face = 'bold', vjust = 5, hjust = 0.5),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'transparent', colour = NA),
        plot.background = element_rect(fill = 'transparent', colour = NA),
        panel.grid = element_blank())
ggsave('~/github/krillGrowthModels/plots/RyabovChla.pdf', plot = p1, 
       width = 7, height = 6)

# now for Atkinson et al. (2006)
atkinsonChla <- read_csv("inputData/originalResults/seaWiFS_Table1.csv")
p2 <- atkinsonChla %>% 
  ggplot(.,aes(x = chlaSeaWiFS)) +
  geom_density(colour = '#d66b6d', size = 6) +
  scale_x_continuous(limits = c(0,7), breaks = c(0,7), labels = c(0,7)) +
  scale_y_continuous(expand = c(0.1,0.1), breaks = c(0,1), labels = c(0,1), limits = c(0,1)) +
  labs(x = expression(mg~m^{-3})) +
  theme(axis.text.y = element_text(size = 55, colour = '#141414', face = 'bold', vjust = 0),#element_blank(),
        axis.text.x = element_text(size = 55, colour = '#141414', face = 'bold', vjust = 0),
        axis.line.x = element_line(size = 3.5, colour = '#141414'),
        axis.line.y = element_line(size = 3.5, colour = '#141414'),
        axis.title.x = element_text(size = 55, colour = '#141414', face = 'bold', vjust = 5, hjust = 0.5),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'transparent', colour = NA),
        plot.background = element_rect(fill = 'transparent', colour = NA),
        panel.grid = element_blank())
ggsave('~/github/krillGrowthModels/plots/AtkinsonChla.pdf', plot = p2, 
       width = 7, height = 6)

# now for Hofmann and Lascara (2000)
HofmannLascaraEnv <- read_delim("inputData/originalResults/HofmannLascaraEnv.csv", 
                                ";", escape_double = FALSE, trim_ws = TRUE)
p3 <- HofmannLascaraEnv %>% 
  filter(scenario == 'coastal') %>% 
  ggplot(.,aes(x = chlorophyllA)) +
  geom_density(colour = '#d66b6d', size = 6) +
  scale_x_continuous(limits = c(0,7), breaks = c(0,7), labels = c(0,7)) +
  scale_y_continuous(expand = c(0.1,0.1), breaks = c(0,1), labels = c(0,1), limits = c(0,1)) +
  labs(x = expression(mg~m^{-3})) +
  theme(axis.text.y = element_text(size = 55, colour = '#141414', face = 'bold', vjust = 0),#element_blank(),
        axis.text.x = element_text(size = 55, colour = '#141414', face = 'bold', vjust = 0),
        axis.line.x = element_line(size = 3.5, colour = '#141414'),
        axis.line.y = element_line(size = 3.5, colour = '#141414'),
        axis.title.x = element_text(size = 55, colour = '#141414', face = 'bold', vjust = 5, hjust = 0.5),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'transparent', colour = NA),
        plot.background = element_rect(fill = 'transparent', colour = NA),
        panel.grid = element_blank())
ggsave('~/github/krillGrowthModels/plots/HofmannChla.pdf', plot = p3, 
       width = 7, height = 6)

# now for Bahlburg et al. (2021)
bahlburgEnv <- read_csv("inputData/originalResults/bahlburgEtAlTestEnv.csv")
p4 <- bahlburgEnv %>% 
  ggplot(.,aes(x = foodConc))+
  geom_density(colour = '#d66b6d', size = 6) +
  scale_x_continuous(limits = c(0,7), breaks = c(0,7), labels = c(0,7)) +
  scale_y_continuous(expand = c(0.1,0.1), breaks = c(0,1), labels = c(0,1), limits = c(0,1)) +
  labs(x = expression(mg~m^{-3})) +
  theme(axis.text.y = element_text(size = 55, colour = '#141414', face = 'bold', vjust = 0),#element_blank(),
        axis.text.x = element_text(size = 55, colour = '#141414', face = 'bold', vjust = 0),
        axis.line.x = element_line(size = 3.5, colour = '#141414'),
        axis.line.y = element_line(size = 3.5, colour = '#141414'),
        axis.title.x = element_text(size = 55, colour = '#141414', face = 'bold', vjust = 5, hjust = 0.5),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'transparent', colour = NA),
        plot.background = element_rect(fill = 'transparent', colour = NA),
        panel.grid = element_blank())
ggsave('~/github/krillGrowthModels/plots/BahlburgChla.pdf', plot = p4, 
       width = 7, height = 6)

# now create blank plot for models that do not report chlorophyll a values for calibration
p5 <- bahlburgEnv %>% 
  ggplot(.,aes(x = foodConc))+
  geom_density(colour = '#d66b6d', size = 0) +
  scale_x_continuous(limits = c(0,7), breaks = c(0,7), labels = c(0,7)) +
  scale_y_continuous(expand = c(0.1,0.1), breaks = c(0,1), labels = c(0,1), limits = c(0,1)) +
  labs(x = expression(mg~m^{-3})) +
  theme(axis.text.y = element_text(size = 55, colour = '#141414', face = 'bold', vjust = 0),#element_blank(),
        axis.text.x = element_text(size = 55, colour = '#141414', face = 'bold', vjust = 0),
        axis.line.x = element_line(size = 3.5, colour = '#141414'),
        axis.line.y = element_line(size = 3.5, colour = '#141414'),
        axis.title.x = element_text(size = 55, colour = '#141414', face = 'bold', vjust = 5, hjust = 0.5),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'transparent', colour = NA),
        plot.background = element_rect(fill = 'transparent', colour = NA),
        panel.grid = element_blank())
ggsave('~/github/krillGrowthModels/plots/BlankChla.pdf', plot = p5, 
       width = 7, height = 6)
# ----------------------------------------------------------------------------------------------------------- #
# Appendix plots
# sex-specific parameterizations for Atkinson et al (2006) and Tarling et al. (2006)
# list files with model run results
tarlingSexSpecific <- list.files('simulationResults', pattern = 'TarlingSummerGrowthdoy306_doy105', recursive = T, full.names = T)
atkinsonSexSpecific <- list.files('simulationResults', pattern = 'AtkinsonSummerGrowthdoy306_doy105', recursive = T, full.names = T)

# import the results of the different model runs and stack them
modelResults <- rast(c(atkinsonSexSpecific, tarlingSexSpecific))
nModels <- 6

# Plotting length differences instead of total body length
extractSimDay0 <- 1
extractSimDay <- 166
modelResultsDay0 <- subset(modelResults, seq(extractSimDay0, extractSimDay0 + (166 * (nModels-1)), length.out = nModels))
modelResultsDayMax <- subset(modelResults, seq(extractSimDay, extractSimDay + (166 * (nModels-1)), length.out = nModels))

modelResultsDiff <- modelResultsDayMax - modelResultsDay0
modelResultsDiff[modelResultsDiff < -26] <- NA


nameVec <- c('Atkinson et al. (2006) \n all krill', 'Atkinson et al. (2006) \n female krill','Atkinson et al. (2006) \n male krill',
             'Tarling et al. (2006) \n all krill', 'Tarling et al. (2006) \n female krill','Tarling et al. (2006) \n male krill'
)
modelResultsDiffStereo <- terra::project(modelResultsDiff, myCrs)
# remove data North of the Southern Polar Front
modelResultsDiffStereo <- mask(modelResultsDiffStereo, meanSPF)

# Create model overview plots
for(i in 1:6){
  
  resultsLengthRastStereoNew <- subset(modelResultsDiffStereo,i)
  darkMode <- ggplot() +
    geom_spatraster(data = resultsLengthRastStereoNew) +
    geom_sf(data = coastline, colour = '#cccccc', fill = '#61616150', size = 0.5) +
    scale_fill_stepsn(colours = scico(n = 11, palette = 'vikO', begin = 10/45, end = 1),
                      na.value = '#303030',
                      limits = c(-15,30),
                      breaks = seq(-15,30,by = 5)) +
    labs(fill = expression(Delta~L),#'body length in mm',
         caption = nameVec[i]) +
    theme(panel.background = element_rect(fill = NA, colour = NA),
          panel.grid.major = element_blank(),
          plot.background = element_rect(fill = 'transparent', colour = NA),
          legend.position = 'none',
          #legend.position = 'bottom',
          axis.text = element_blank(),
          plot.caption = element_text(hjust = 0.5, size = 42, colour = '#f5f5f5')) +
    guides(fill = guide_colourbar(title.position = 'top', title.hjust = 0.5))
  
  ggsave(paste('plots/modSensitivityDiff',i,'.pdf',sep = ''), plot = darkMode, width = 10, height = 10, bg = 'transparent') 
}




