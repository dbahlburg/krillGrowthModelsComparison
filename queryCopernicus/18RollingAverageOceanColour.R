# In this script, I process the ocean colour data to reduce sea-ice induced data gaps.
# For that, I'll use a rolling average with a width of 7 days as is common in other remote-sensing products
# This means, that the climatology of e.g. day 20 is calculated as the mean climatology of days 17-23, day 21 by averaging 
# 18-24 [...]
# First load packages
library(terra)
source('functions/dayOfYear.R')

# load ocean colour climatologies
oceanColourDaily <- rast('/Volumes/T7/OceanColourL4-Output/climatologies/climaChla_doy_1_365.nc')

oceanColourClimat7d <- rast()
oceanColourClimat15d <- rast()
oceanColourClimat23d <- rast()
oceanColourClimat31d <- rast()
for(i in 1:365){
  
  averagingWindow7d <- dayOfYearFunc(t = (i-1-3):(i-1+3), startDay = 1)
  averagingWindow15d <- dayOfYearFunc(t = (i-1-7):(i-1+7), startDay = 1)
  averagingWindow23d <- dayOfYearFunc(t = (i-1-11):(i-1+11), startDay = 1)
  averagingWindow31d <- dayOfYearFunc(t = (i-1-15):(i-1+15), startDay = 1)
  
  subsetData7d <-  subset(oceanColourDaily, averagingWindow7d)
  subsetData15d <- subset(oceanColourDaily, averagingWindow15d)
  subsetData23d <- subset(oceanColourDaily, averagingWindow23d)
  subsetData31d <- subset(oceanColourDaily, averagingWindow31d)
  
  averagedClimatology7d <- app(subsetData7d, mean, na.rm = T)
  averagedClimatology15d <- app(subsetData15d, mean, na.rm = T)
  averagedClimatology23d <- app(subsetData23d, mean, na.rm = T)
  averagedClimatology31d <- app(subsetData31d, mean, na.rm = T)
  
  add(oceanColourClimat7d) <- averagedClimatology7d
  add(oceanColourClimat15d) <- averagedClimatology15d
  add(oceanColourClimat23d) <- averagedClimatology23d
  add(oceanColourClimat31d) <- averagedClimatology31d
  print(i)
}
terra::writeCDF(oceanColourClimat7d, 
                varname = 'climat7d_rollingMean',
                overwrite=TRUE,
                filename = '/Volumes/T7/OceanColourL4-Output/climatologies/climaRollingAv7.nc')
terra::writeCDF(oceanColourClimat15d, 
                varname = 'climat15d_rollingMean',
                overwrite=TRUE,
                filename = '/Volumes/T7/OceanColourL4-Output/climatologies/climaRollingAv15.nc')
terra::writeCDF(oceanColourClimat23d, 
                varname = 'climat23d_rollingMean',
                overwrite=TRUE,
                filename = '/Volumes/T7/OceanColourL4-Output/climatologies/climaRollingAv23.nc')
terra::writeCDF(oceanColourClimat31d, 
                varname = 'climat31d_rollingMean',
                overwrite=TRUE,
                filename = '/Volumes/T7/OceanColourL4-Output/climatologies/climaRollingAv31.nc')
# Compare data quality of original dataset and rolling average
# ------------------------------------------------------------------------------------ #
# how does this improve over time
missingDataTib <- tibble(date = seq(as.Date('2021-01-01'), as.Date('2021-12-31'), by = '1 day'),
                         doy = yday(date),
                         week = week(date))
lowResCoast <- vect('inputData/lowResCoastlineAntarctica/GSHHS_l_L5.shp')

oceanColourClimat7d[lowResCoast] <- -999
aggregatedDataGaps7 <- tibble(doy = 1:365,
                             missingData7d = unlist(global(oceanColourClimat7d, fun = 'isNA')))

oceanColourClimat15d[lowResCoast] <- -999
aggregatedDataGaps15 <- tibble(doy = 1:365,
                             missingData15d = unlist(global(oceanColourClimat15d, fun = 'isNA')))

oceanColourClimat23d[lowResCoast] <- -999
aggregatedDataGaps23 <- tibble(doy = 1:365,
                             missingData23d = unlist(global(oceanColourClimat23d, fun = 'isNA')))

oceanColourClimat31d[lowResCoast] <- -999
aggregatedDataGaps31 <- tibble(doy = 1:365,
                             missingData31d = unlist(global(oceanColourClimat31d, fun = 'isNA')))

chlaClimCopy <- oceanColourDaily
chlaClimCopy[lowResCoast] <- -999
dailyDataGaps <- tibble(doy = 1:365,
                        missingDataDaily = unlist(global(chlaClimCopy, fun = 'isNA')))

chlaClimCopy <- app(chlaClimCopy,  function(x){is.na(x) + 2})
chlaClimCopy[lowResCoast] <- -999
chlaClimCopy[chlaClimCopy == -999] <- NA
missingDataTib <- missingDataTib %>% 
  left_join(., aggregatedDataGaps7, by = 'doy') %>% 
  left_join(., aggregatedDataGaps15, by = 'doy') %>% 
  left_join(., aggregatedDataGaps23, by = 'doy') %>% 
  left_join(., aggregatedDataGaps31, by = 'doy') %>% 
  left_join(., dailyDataGaps, by = 'doy') %>% 
  mutate(totalData = unlist(global(chlaClimCopy, fun = 'notNA'))) %>% 
  dplyr::select(c(2,4:9)) %>% 
  gather(meanWindow, value, -doy, -totalData) %>% 
  mutate(dataGaps = value/totalData * 100,
         meanWindow = factor(meanWindow, levels = c("missingDataDaily","missingData7d","missingData15d","missingData23d","missingData31d"),
                             labels = c('1 day','7 days','15 days','23 days','31 days')))

colourFunc <- colorRampPalette(c("#4685e3", "#e83a3a"))
p4 <- missingDataTib  %>% 
  filter(meanWindow != '1 day') %>% 
  ggplot(.,aes(x = doy, y = dataGaps, colour = meanWindow)) +
  geom_line(size = 1, linetype = '22') +
  geom_line(data = filter(missingDataTib, meanWindow == '1 day'), size = 1, colour = '#383838') +
  scale_colour_manual(values = colourFunc(4)) +
  scale_linetype_manual(values = c('42','42','42','42','42')) +
  scale_y_continuous(limits = c(0,100)) +
  labs(x = 'day of year',
       colour = '',
       y = 'data gaps as % of ocean grid cells',
       title = 'data gaps in ocean colour: \n using rolling means vs. daily climatologies') +
  theme(panel.background = element_rect(fill = NA, colour = '#303030'),
        panel.grid.major = element_blank(),
        legend.position = 'bottom',
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


