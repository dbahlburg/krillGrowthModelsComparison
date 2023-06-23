# Create growth plains
# To shed more light into the question why different krill growth models predict drastically differing
# growth trajectories, we create "growth plains". These plains map predicted daily growth rate for krill
# individuals onto a grid of chlorophyll a concentrations over a range of temperatures.
# load required packages:
library(tidyverse) 
library(lubridate)
library(metR)
library(scico)
library(terra)
# ------------------------------------------------------------------------------------------------------- #
# source the different krill growth models
source('functions/AtkinsonEtAl2006AuxiliaryFunctions.R')
source('functions/BahlburgEtAl2021AuxiliaryFunctions.R')
source('functions/ConstableKawaguchi2017AuxiliaryFunctions.R')
source('functions/FachEtAl2002AuxiliaryFunctions.R')
source('functions/HofmannLascara2000AuxiliaryFunctions.R')
source('functions/JagerRavagnanAuxiliaryFunctions.R')
source('functions/RyabovEtAl2017AuxiliaryFunctions.R')
source('functions/TarlingEtAl2006AuxiliaryFunctions.R')
source('functions/WiedenmannEtAl2008AuxiliaryFunctions.R')
lookUpTableTimeVariables <- read.csv('inputData/originalResults/timeDependentVariables.csv')
lookUpTableSizeClassTraits <- read.csv('inputData/originalResults/sizeClassTraits.csv')
# ------------------------------------------------------------------------------------------------------- #
# Create a grid of temperature and chlorophyll a gradients
# Some models require POC concentrations as input - chlorophyll a will be converted to POC using 
# a conversion factor of 50 (Hofmann and Lascara 2000)
growthPlain <- expand_grid(chla = seq(0,5, length.out = 200),
                           temp = seq(-2,3,length.out = 200)) %>% 
  mutate(poc = chla * 50)
# ------------------------------------------------------------------------------------------------------- #
# Set some parameters... 
# Some models require more input variables than just temperature and chlorophyll a
# We set fixed values for photoperiod, day of year, body length, ice algae
# The values are chosen to resemble peak summer growth conditions
photoperiod <- 20
dayOfYear <- 15
inputDayOfYear <- 15
startDay <- 1
bodyLength <- 26
iceAlgae <- 0
hetC <- 0

# conversion rate chlorophyll to carbon
chlorophyllToCarbon <- 50
# assimilation efficiency of ingested carbon
assimilationRate <- 0.8
# vector containing the different size classes (used in the model to assign current size class of individual)
sizeClassesVec <- seq(1,90, by = 0.25)
# ------------------------------------------------------------------------------------------------------- #
# predict the growth rates
growthPlainDGR <- growthPlain %>% 
  mutate(#Fach et al. (2002) DGR
         fachChangeMass = rk4stepFach(state = massAtLength(inputLength = bodyLength, massType = 'carbonMass'),
                                      time = dayOfYear,
                                      inputChla = chla,
                                      inputTemp = temp,
                                      inputHetC = hetC,
                                      inputIceAlgae = iceAlgae),
         `Fach et al. (2002)` = lengthAtMass(inputMass = massAtLength(inputLength = bodyLength, massType = 'carbonMass') + fachChangeMass) -
           lengthAtMass(inputMass = massAtLength(inputLength = bodyLength, massType = 'carbonMass')),
         #Hofmann and Lascara (2000) DGR
         hofmannChangeMass = rk4stepHofmann(state = massAtLength(inputLength = bodyLength, massType = 'carbonMass'),
                                            time = dayOfYear,
                                            inputChla = chla,
                                            inputIceAlgae = iceAlgae),
         `Hofmann and Lascara (2000)` = lengthAtMass(inputMass = massAtLength(inputLength = bodyLength, massType = 'carbonMass') + hofmannChangeMass) -
           lengthAtMass(inputMass = massAtLength(inputLength = bodyLength, massType = 'carbonMass'))) 

growthPlainDGR <- growthPlainDGR %>% 
  rowwise() %>% 
  mutate(#Atkinson et al. (2006) DGR
          `Atkinson et al. (2006)*` = unlist(atkinsonFunc(inputFood = chla, inputTemperature = temp, inputStage = 1, inputLength = bodyLength))[1],
         #Bahlburg et al. (2021) DGR
         `Bahlburg et al. (2021)` = bahlburgEtAlGrowth(inputLength = bodyLength,
                                          time = dayOfYear,
                                          eggBuffer = 0,
                                          reproBuffer = 0,
                                          temperature = temp,
                                          photoperiod = photoperiod,
                                          chlorophyllA = chla)[1],
         #Jager and Ravagnan (2015) DGR
         `Jager and Ravagnan (2015)` = jagerRavagnanGrowth(inputLength = bodyLength,
                                        time = dayOfYear,
                                        eggBuffer = 0,
                                        reproBuffer = 0,
                                        temperature = temp,
                                        chlorophyllA = chla)[1],
         #Ryabov et al. (2017)
         `Ryabov et al. (2017)` = rk4stepRyabov(state = bodyLength,
                                   time = dayOfYear,
                                   inputChla = chla,
                                   inputIceAlgae = iceAlgae,
                                   inputAge = 400),
         #Tarling et al. (2006)
         `Tarling et al. (2006)*` = TarlingEtAl2006Model(inputLength = bodyLength,
                                inputTemperature = temp,
                                inputChla = chla,
                                inputStage = 1,
                                time = TarlingEtAl2006IMP(bodyLength = bodyLength,
                                                          temperature = temp,
                                                          stage = 1, roundValues = F),
                                oldMoultDay = moultDay,
                                temperatureHistory = temp,
                                chlorophyllHistory = chla,
                                moultDay = TarlingEtAl2006IMP(bodyLength = bodyLength,
                                                              temperature = temp,
                                                              stage = 1, roundValues = F))[1]/TarlingEtAl2006IMP(bodyLength = bodyLength,
                                                                                                temperature = temp,
                                                                                                stage = 1, roundValues = F),
         #Wiedenmann et al. (2008)
         `Wiedenmann et al. (2008)*` = WiedenmannEtAl2008Model(inputLength = bodyLength, 
                                                 inputTemperature = temp, 
                                                 inputChla = chla, 
                                                 time = exp(3.5371 - 0.5358 * log(temp + 2)),
                                                 oldMoultDay = moultDay,
                                                 temperatureHistory = temp,
                                                 chlorophyllHistory = chla,
                                                 moultDay = exp(3.5371 - 0.5358 * log(temp + 2)))[1]/exp(3.5371 - 0.5358 * log(temp + 2)))


# create long format for tibble
growthPlainDGR <- growthPlainDGR %>%
  select(-hofmannChangeMass, -fachChangeMass) %>% 
  gather(model, dgr, -chla, -temp, -poc)

# extract growth trajectories for -60.959390, -45.458629 (South Orkney Islands)
chlorophyllClimatology <- rast('inputData/climatologies/climaChlaRollingAv15.nc')
sstClimatology <- rast('inputData/climatologies/climatologySST1_365.nc')
chlorophyllClimatology <- flip(chlorophyllClimatology, direction = 'vertical')
crs(chlorophyllClimatology) <- 'epsg:4258'

coords <- data.frame(x = -45.458629, y =  -60.959390)
chlaDyn <- terra::extract(chlorophyllClimatology, coords) 
tempDyn <- terra::extract(sstClimatology, coords)

chlaDynSummer <- unlist(chlaDyn[1, (dayOfYearFunc(1:165, startDay = 306) + 1)])
tempDynSummer <- unlist(tempDyn[1, (dayOfYearFunc(1:165, startDay = 306) + 1)]) - 273.15

envDynamics <- as.data.frame(cbind(chlaDynSummer, tempDynSummer))
envDynamicsStSt <- envDynamics %>% 
  mutate(timeInd = 1:n()) %>% 
  filter(timeInd %in% c(min(timeInd), max(timeInd))) %>% 
  mutate(type = c('start','stop'))

# create plot
p1 <- growthPlainDGR %>% 
  mutate(model = factor(model, levels = levels(factor(growthPlainDGR$model))[c(1,7,8,6,2,5,4,3)])) %>% 
  ggplot(.) +
  geom_raster(aes(x = chla, y = temp, fill = dgr)) +
  geom_path(data = envDynamics, aes(x = chlaDynSummer, y = tempDynSummer), 
            size = 0.75,
            colour = '#49bee6') +
  geom_point(data = envDynamicsStSt, aes(x = chlaDynSummer, y = tempDynSummer, colour = type),
             size = 1.3,
             show.legend=FALSE) +
  geom_contour(aes(x = chla, y = temp, z = dgr), breaks = c(0, 0.1, 0.2, 0.3, 0.4), colour = '#f5f5f5',
               linetype = '22') +
  geom_text_contour(aes(x = chla, y = temp, z = dgr), colour = '#f5f5f5',
                    label.placer = label_placer_fraction(
                      frac = 0.75,
                      rot_adjuster = isoband::angle_halfcircle_bottom()
                    ),
                    nudge_y = 0.12,
                    nudge_x = 0.12,
                    skip = 0,
                    breaks = c(0, 0.1,0.2,0.3, 0.4)) +
  facet_wrap(~model, ncol = 2) +
  scale_colour_manual(values = c('#242424','#f5f5f5')) +
  scale_fill_scico(palette = 'romaO', 
                   na.value = NA, 
                   midpoint = 0,
                   limits = c(-0.2,0.4), 
                   breaks = c(-0.2,0,0.2,0.4),
                   labels = c(-0.2,0,0.2,0.4), 
                   oob = scales::squish) +
  labs(x = expression(chlorophyll~a~concentration~mg~m^{-3}),
       y = 'temperature °C',
       fill = expression(growth~mm~d^{-1})) +
  theme(panel.background = element_rect(fill = NA, colour = '#2b2b2b'),
        panel.grid = element_blank(),
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 18),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(size = 16),
        legend.position = 'bottom',
        legend.key.width=unit(1.15,"cm"),
        legend.key.height=unit(0.175,"cm"),
        legend.title = element_text(size = 18, hjust = 0.5),
        legend.text = element_text(size = 18)) +
  guides(fill = guide_colourbar(title.position = "top"))

# save plot
ggsave('plots/dgrPlains.pdf', plot = p1, width = 9, height = 12)


