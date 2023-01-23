# 16.8.2022
# Dominik Bahlburg
# Here, I will check all the different krill growth models for correct implementations.
# The procedure works as follows:
# 1. extract the environmental input data used in the original papers (using data extractor tools from ImageJ)
# 2. run my own implementation of the respective krill growth model using the original environmental data
# 3. compare my predicted growth trajectories with those reported in the original paper
# --------------------------------------------------------------------------------------------------------------- #
# load packages
library(tidyverse)
# --------------------------------------------------------------------------------------------------------------- #
# Hofmann and Lascara (2000)
# load original results extracted from Figure 7 and Figure 8
hofmannLascaraResultsFig7 <- read_delim("inputData/originalResults/hofmannLascara2000Fig7.csv", ";", escape_double = FALSE, trim_ws = TRUE)
hofmannLascaraResultsFig8 <- read_delim("inputData/originalResults/hofmannLascara2000Fig8.csv", ";", escape_double = FALSE, trim_ws = TRUE)

# load model and auxiliary functions
source('functions/hofmannLascara2000AuxiliaryFunctions.R')

# load lookup-tables to reduce simulation time (this is how the original model was run)
lookUpTableTimeVariables <- read.csv('inputData/originalResults/timeDependentVariables.csv') 
lookUpTableSizeClassTraits <- read.csv('inputData/originalResults/sizeClassTraits.csv')

# in the original paper, the model is run using the output of an environmental model of the Western Antarctic Peninsula region
# several scenarios were simulated, I test my implementation using the "coastal" scenario
# DEFINE PARAMETER VALUES
# conversion rate chlorophyll to carbon
chlorophyllToCarbon <- 50

# assimilation efficiency of ingested carbon
assimilationRate <- 0.8

# vector containing the different size classes (used in the model to assign current size class of individual)
sizeClassesVec <- seq(1,90, by = 0.25)

# Initialize state variables
dt <- 1
timeVec <- seq(1,500-32, by = dt)
index <- 2
nIndividuals <- 3
results <- matrix(ncol = nIndividuals, nrow = length(timeVec))
results[1,1:nIndividuals] <- massAtLength(inputLength = c(2,22,45), massType = 'carbonMass')
resultsWithIceAlgae <- matrix(ncol = nIndividuals, nrow = length(timeVec))
resultsWithIceAlgae[1,1:nIndividuals] <- massAtLength(inputLength = c(2,22,45), massType = 'carbonMass')
timeInd <- timeVec[1]
startDay <- 32

# execute the model
for(i in timeVec[-1]){
  doy <- dayOfYearFunc(t = i, startDay = startDay)
  results[index,] <- results[index-1,] + rk4stepHofmann(dt = dt, 
                                                 state = results[index-1,],
                                                 time = timeInd,
                                                 inputChla = lookUpTableTimeVariables$coastal[doy],
                                                 inputIceAlgae = 0)
  
  resultsWithIceAlgae[index,] <- resultsWithIceAlgae[index-1,] + rk4stepHofmann(dt = dt, 
                                                                         state = resultsWithIceAlgae[index-1,],
                                                                         time = timeInd,
                                                                         inputChla = lookUpTableTimeVariables$coastal[doy],
                                                                         inputIceAlgae = lookUpTableTimeVariables$iceAlgae[doy])
  timeInd <- timeInd + dt
  index <- index+1
  print(i)
}

# transform the results to dataframe, convert carbon mass to length
hofmannLascaraFig7 <- as.data.frame(results) %>% 
  mutate(doy = timeVec+startDay) %>% 
  gather(individual, cMass, -doy) %>% 
  mutate(length = lengthAtMass(inputMass = cMass, massType = 'carbonMass'))
hofmannLascaraFig8 <- as.data.frame(resultsWithIceAlgae) %>% 
  mutate(doy = timeVec+startDay) %>% 
  gather(individual, cMass, -doy) %>% 
  mutate(length = lengthAtMass(inputMass = cMass, massType = 'carbonMass'))

hofmannLascaraComparisonFig7 <- hofmannLascaraFig7 %>%
  ggplot(.,aes(x = doy, y = length)) +
  geom_point(data = hofmannLascaraResultsFig7, aes(x = doy, y = length), colour = '#d43535') +
  geom_line(aes(group = individual), size = 1) +
  scale_y_continuous(limits = c(0,70)) +
  scale_x_continuous(limits = c(1,500), breaks = c(1,seq(50,500,by = 50))) +
  labs(x = 'Time (Julian Day)',
       y = 'Krill size (mm)',
       title = 'Hofmann and Lascara (2000) \n quality control Fig 7') +
  theme(panel.background = element_rect(fill = NA, colour = '#3d3d3d', size = 1.5),
        legend.title = element_text(vjust = 1, size = 24, colour = '#3d3d3d'),
        legend.text = element_text(size = 24, colour = '#3d3d3d'),
        axis.text = element_text(size = 25, colour = '#3d3d3d'),
        axis.title = element_text(size = 25, colour = '#3d3d3d'),
        plot.title = element_text(hjust = 0.5, size = 30, colour = '#3d3d3d'))
ggsave('plots/hofmannLascara2000QC1.pdf', plot = hofmannLascaraComparisonFig7, width = 10, height = 7.5)

hofmannLascaraComparisonFig8 <- hofmannLascaraFig8 %>%
  ggplot(.,aes(x = doy, y = length)) +
  geom_point(data = hofmannLascaraResultsFig8, aes(x = doy, y = length), colour = '#d43535') +
  geom_line(aes(group = individual), size = 1) +
  scale_y_continuous(limits = c(0,70)) +
  scale_x_continuous(limits = c(1,500), breaks = c(1,seq(50,500,by = 50))) +
  labs(x = 'Time (Julian Day)',
       y = 'Krill size (mm)',
       title = 'Hofmann and Lascara (2000) \n quality control Fig 8') +
  theme(panel.background = element_rect(fill = NA, colour = '#3d3d3d', size = 1.5),
        legend.title = element_text(vjust = 1, size = 24, colour = '#3d3d3d'),
        legend.text = element_text(size = 24, colour = '#3d3d3d'),
        axis.text = element_text(size = 25, colour = '#3d3d3d'),
        axis.title = element_text(size = 25, colour = '#3d3d3d'),
        plot.title = element_text(hjust = 0.5, size = 30, colour = '#3d3d3d'))
ggsave('plots/hofmannLascara2000QC2.pdf', plot = hofmannLascaraComparisonFig8, width = 10, height = 7.5)

# --------------------------------------------------------------------------------------------------------------- #
# Ryabov et al. (2017)
# load original results, extracted from Alex' Matlab-Script
ryabovEtAlResults <- read_delim("inputData/originalResults/ryabovEtAl2017.csv", ";", escape_double = FALSE, trim_ws = TRUE)
# load model and auxiliary functions
source('functions/RyabovEtAl2017AuxiliaryFunctions.R')

# load environmental time series sued in Alex' Matlab-Script
ryabovEnv <- read_csv("inputData/originalResults/ryabovTestEnvironment.csv") %>% 
  mutate(iceAlgae = iceAlgae/50)

# Initialize state variables
dt <- 1
timeVec <- seq(1,365 * 6, by = dt)
index <- 2
nIndividuals <- 3
results <- matrix(ncol = nIndividuals, nrow = length(timeVec))
results[1,1:nIndividuals] <- c(10,30,50)
timeInd <- timeVec[1]
startDay <- 1
krillAge <- seq(50, 5 * 365, length.out = nIndividuals)

# execute the model
for(i in timeVec[-1]){
  doy <- dayOfYearFunc(t = i, startDay = startDay)
  results[index,] <- results[index-1,] + rk4stepRyabov(state = results[index-1,],
                                                 time = doy,
                                                 inputChla = ryabovEnv$chlA[doy],
                                                 inputIceAlgae = ryabovEnv$iceAlgae[doy],
                                                 inputAge = krillAge)
  timeInd <- timeInd + dt
  krillAge <- krillAge + dt
  index <- index+1
  print(i)
}

# transform the results to dataframe, convert carbon mass to length
ryabovEtAlSolved <- as.data.frame(results) %>% 
  mutate(age = (timeVec + startDay)/365) %>% 
  gather(individual, length, -age) 

ryabovEtAlComparison <- ryabovEtAlSolved %>%
  ggplot(.,aes(x = age, y = length)) +
  geom_point(data = ryabovEtAlResults, aes(x = time, y = length), colour = '#d43535') +
  geom_line(aes(group = individual), size = 1) +
  scale_y_continuous(limits = c(0,70)) +
  scale_x_continuous(limits = c(0,6), breaks = c(0:6)) +
  labs(x = 'Age (years)',
       y = 'Krill size (mm)',
       title = 'Ryabov et al. (2017) quality control') +
  theme(panel.background = element_rect(fill = NA, colour = '#3d3d3d', size = 1.5),
        legend.title = element_text(vjust = 1, size = 24, colour = '#3d3d3d'),
        legend.text = element_text(size = 24, colour = '#3d3d3d'),
        axis.text = element_text(size = 25, colour = '#3d3d3d'),
        axis.title = element_text(size = 25, colour = '#3d3d3d'),
        plot.title = element_text(hjust = 0.5, size = 30, colour = '#3d3d3d'))

ggsave('plots/ryabovEtAl2017QC.pdf', plot = ryabovEtAlComparison, width = 10, height = 7.5)
# --------------------------------------------------------------------------------------------------------------- #
# Atkinson et al. (2006)
# load original results, extracted from Figure 6
atkinsonEtAlResults <- read_delim("inputData/originalResults/atkinsonEtAl2006.csv", ";", escape_double = FALSE, trim_ws = TRUE) %>% 
  mutate(temp = paste(temp, '°C', sep = ' '))
# load model an auxiliary functions
source('functions/AtkinsonEtAl2006AuxiliaryFunctions.R')

# Calculate predicted daily growth rates at 0.5 and 4 deg C over a chla-gradient of 0-4.5 mg/m-3
# and 3 different body lengths
atkinsonControl <- expand_grid(chla = seq(0,4.5, by = 0.05),
                               temp = c(0.5, 4),
                               length = c(25,40,60)) %>% 
  rowwise() %>% 
  mutate(dgr = atkinsonFunc(inputFood = chla, inputTemperature = temp, inputStage = 3, inputLength = length)[1],
         temp = paste(temp, '°C', sep = ' '))

atkLabels <- atkinsonControl %>% 
  group_by(length) %>% 
  filter(chla == max(chla)) %>% 
  mutate(label = paste(length,'mm', sep = ' '))

atkinsonComparison <- atkinsonControl %>% 
  ggplot(.) +
  geom_point(data = atkinsonEtAlResults, aes(x = food, y = dgr), colour = '#d43535', size = 1.8) +
  geom_line(aes(x = chla, y = dgr, group = length), size = 1) +
  geom_text(data = atkLabels, aes(x = 5, y = dgr, label = label), size = 6) +
  facet_grid(~temp) +
  scale_y_continuous(limits = c(-0.3,0.5)) +
  scale_x_continuous(limits = c(0,5.5), breaks = c(0:5.5)) +
  labs(x = expression(chlorophyll~a~'('*mg~m^{-3}*')'),
       y = expression(daily~growth~rate~"("*mm~d^{-1}*")"),
       title = 'Atkinson et al. (2006) quality control') +
  theme(panel.background = element_rect(fill = NA, colour = '#3d3d3d', size = 1.5),
        legend.title = element_text(vjust = 1, size = 24, colour = '#3d3d3d'),
        legend.text = element_text(size = 24, colour = '#3d3d3d'),
        axis.text = element_text(size = 25, colour = '#3d3d3d'),
        axis.title = element_text(size = 25, colour = '#3d3d3d'),
        strip.text = element_text(size = 25, colour = '#3d3d3d'),
        strip.background = element_rect(fill = NA, colour = NA),
        plot.title = element_text(hjust = 0.5, size = 30, colour = '#3d3d3d'))
atkinsonComparison
ggsave('plots/atkinsonEtAl2006QC.pdf', plot = atkinsonComparison, width = 14, height = 7.5)
# --------------------------------------------------------------------------------------------------------------- #
# Tarling et al. (2006)
# load original results, extracted from Figure 7 and 8
tarlingEtAlFig7 <- read_delim("inputData/originalResults/tarlingEtAl2006Fig7.csv", ";", escape_double = FALSE, trim_ws = TRUE) %>% 
  mutate(temp = paste(temp, '°C', sep = ' '))
tarlingEtAlFig8 <- read_delim("inputData/originalResults/tarlingEtAl2006Fig8.csv", ";", escape_double = FALSE, trim_ws = TRUE) %>% 
  mutate(length = paste(length, 'mm', sep = ' '))
# load model an auxiliary functions
source('functions/TarlingEtAl2006AuxiliaryFunctions.R')

# Calculate predicted daily growth rates at 0.5 and 4 deg C over a chla-gradient of 0-4.5 mg/m-3
# and 3 different body lengths
tarlingFig7Ctrl <- expand_grid(length = seq(30,62,by = 0.05),
                               temp = c(0,2,4)) %>% 
  mutate(imp = TarlingEtAl2006IMP(bodyLength = length, stage = 3, temperature = temp),
         temp = paste(temp, '°C', sep = ' '))

tarlingFig7 <- tarlingFig7Ctrl %>% 
  ggplot(.,aes(x = length, y = imp)) +
  geom_point(data = tarlingEtAlFig7, aes(x = length, y = imp), colour = '#d43535', size = 1.8) +
  geom_line() +
  facet_grid(~temp)+
  labs(x = 'length (mm)',
       y = 'IMP (days)',
       title = 'Tarling et al. (2006) quality control I') +
  theme(panel.background = element_rect(fill = NA, colour = '#3d3d3d', size = 1.5),
        legend.title = element_text(vjust = 1, size = 24, colour = '#3d3d3d'),
        legend.text = element_text(size = 24, colour = '#3d3d3d'),
        axis.text = element_text(size = 25, colour = '#3d3d3d'),
        axis.title = element_text(size = 25, colour = '#3d3d3d'),
        strip.text = element_text(size = 25, colour = '#3d3d3d'),
        strip.background = element_rect(fill = NA, colour = NA),
        plot.title = element_text(hjust = 0.5, size = 30, colour = '#3d3d3d'))
ggsave('plots/tarlingEtAl2006QC1.pdf', plot = tarlingFig7, width = 14, height = 6)

# repeat for Figure 8
tarlingFig8Ctrl <- expand_grid(temp = seq(-2,4.8,by = 0.05),
                               length = c(50,60)) %>% 
  mutate(imp = TarlingEtAl2006IMP(bodyLength = length, stage = 3, temperature = temp),
         length = paste(length, 'mm', sep = ' '))

tarlingFig8 <- tarlingFig8Ctrl %>% 
  ggplot(.,aes(x = temp, y = imp)) +
  geom_point(data = tarlingEtAlFig8, aes(x = temp, y = imp), colour = '#d43535', size = 1.8) +
  geom_line() +
  facet_grid(~length)+
  labs(x = 'temperature [°C]',
       y = 'IMP (days)',
       title = 'Tarling et al. (2006) quality control II') +
  theme(panel.background = element_rect(fill = NA, colour = '#3d3d3d', size = 1.5),
        legend.title = element_text(vjust = 1, size = 24, colour = '#3d3d3d'),
        legend.text = element_text(size = 24, colour = '#3d3d3d'),
        axis.text = element_text(size = 25, colour = '#3d3d3d'),
        axis.title = element_text(size = 25, colour = '#3d3d3d'),
        strip.text = element_text(size = 25, colour = '#3d3d3d'),
        strip.background = element_rect(fill = NA, colour = NA),
        plot.title = element_text(hjust = 0.5, size = 30, colour = '#3d3d3d'))
ggsave('plots/tarlingEtAl2006QC2.pdf', plot = tarlingFig8, width = 14, height = 7.5)

# add direct comparison with IMP-function from Kawaguchi et al. (2006)
tarlingVsKawaguchi <- expand_grid(length = seq(30,62, by = 0.05),
                                  temp = seq(-1, 4, length.out = 350)) %>% 
  mutate(tarlingIMP = TarlingEtAl2006IMP(bodyLength = length, stage = 3, temperature = temp),
         kawaguchiIMP = exp(3.5371 - 0.5358 * log(temp + 2))) %>% 
  gather(model, impValue, -length, -temp)


tarlingKawaLines <- tarlingVsKawaguchi %>% 
  filter(length %in% c(30,40,50)) %>% 
  mutate(length = paste(length, ' mm'),
         model = ifelse(model == 'kawaguchiIMP', 'Kawaguchi et al. (2006)','Tarling et al. (2006)')) %>% 
  ggplot(.,aes(x = temp, y = impValue, colour = model)) +
  geom_line() +
  scale_colour_manual(values = c('#384863','#f58905')) +
  facet_grid(~length) +
  labs(x = 'temperature °C',
       y = 'IMP [d]',
       colour = '') +
  theme(panel.background = element_rect(fill = NA, colour = '#2b2b2b'),
        panel.grid = element_blank(),
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 18),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(size = 16),
        legend.position = 'bottom',
        legend.key.width=unit(2,"cm"),
        legend.key = element_rect(fill = NA),
        legend.key.height=unit(0.2,"cm"),
        legend.title = element_text(size = 18, hjust = 0.5),
        legend.text = element_text(size = 18)) +
  guides(fill = guide_colourbar(title.position = "top"))
ggsave('plots/tarlingKawaLines.pdf', plot = tarlingKawaLines, width = 12, height = 5)

# IMP-plains for Kawaguchi et al. (2006) and Tarling et al. (2006)
tarlingKawaPlains <- tarlingVsKawaguchi %>% 
  ggplot(.,aes(x = temp, y = length, fill = impValue)) +
  geom_raster() +
  geom_contour(aes(x = temp, y = length, z = impValue), breaks = seq(10,50, by = 10), colour = '#f5f5f5',
               linetype = '22') +
  geom_text_contour(aes(x = temp, y = length, z = impValue), colour = '#f5f5f5',
                    label.placer = label_placer_fraction(
                      frac = 0.75,
                      rot_adjuster = isoband::angle_halfcircle_bottom()
                    ),
                    nudge_y = 0.12,
                    nudge_x = 0.12,
                    skip = 0,
                    size = 6,
                    breaks = seq(10,50, by = 10)) +
  scale_fill_scico(palette = 'bamO') +
  facet_grid(~model) +
  labs(x = 'temperature °C',
       y = 'body length [mm]',
       fill = 'IMP [d]') +
  theme(panel.background = element_rect(fill = NA, colour = '#2b2b2b'),
        panel.grid = element_blank(),
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 18),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(size = 16),
        legend.position = 'bottom',
        legend.key.width=unit(3,"cm"),
        legend.key.height=unit(0.2,"cm"),
        legend.title = element_text(size = 18, hjust = 0.5),
        legend.text = element_text(size = 18)) +
  guides(fill = guide_colourbar(title.position = "top"))
ggsave('plots/tarlingKawaPlains.pdf', plot = tarlingKawaPlains, width = 12, height = 5.5)
# --------------------------------------------------------------------------------------------------------------- #
# Wiedenmann et al. (2008)
# load original results, extracted from Figure 7 and 8
kawaguchiEtAlFig5 <- read_csv("inputData/originalResults/kawaguchiEtAl2006Fig5.csv") 

wiedenmannCtrl <- tibble(temperature = seq(-1.8,5,by = 0.05)) %>% 
  mutate(imp = exp(3.5371 - 0.5358 * log(temperature + 2)))

wiedenmannCtrlPlot <- wiedenmannCtrl %>% 
  ggplot(.,aes(x = temperature, y = imp)) +
  geom_point(data = kawaguchiEtAlFig5, aes(x = temp, y = imp), colour = '#d43535') +
  geom_line() +
  labs(x = 'temperature (°C)',
       y = 'IMP (days)',
       title = 'Wiedenmann et al. (2008) quality control') +
  theme(panel.background = element_rect(fill = NA, colour = '#3d3d3d', size = 1.5),
        legend.title = element_text(vjust = 1, size = 24, colour = '#3d3d3d'),
        legend.text = element_text(size = 24, colour = '#3d3d3d'),
        axis.text = element_text(size = 25, colour = '#3d3d3d'),
        axis.title = element_text(size = 25, colour = '#3d3d3d'),
        strip.text = element_text(size = 25, colour = '#3d3d3d'),
        strip.background = element_rect(fill = NA, colour = NA),
        plot.title = element_text(hjust = 0.5, size = 30, colour = '#3d3d3d'))
ggsave('plots/wiedenmannEtAl2008QC.pdf', plot = wiedenmannCtrlPlot, width = 10, height = 7.5)

# --------------------------------------------------------------------------------------------------------------- #
# Jager and Ravagnan (2015)
source('functions/JagerRavagnanAuxiliaryFunctions.R')
dgrDat <- read.csv('inputData/originalResults/jagerRavagnandgrData.csv', sep = ';') %>% 
  rowwise() %>% 
  mutate(modelOriginal = dgr,
         modelRebuilt = jagerRavagnanGrowth(inputLength = length,
                                            time = 10,
                                            eggBuffer = 0,
                                            reproBuffer = 0,
                                            temperature = ifelse(type == 'dgr0deg', 0, 5),
                                            chlorophyllA = 1)[1])

jagerRavagnanCtrlPlot <- dgrDat %>% 
  ggplot(.,aes(x = length, shape = type)) +
  geom_point(aes(y = modelOriginal), colour = '#d43535', fill = '#d43535', size = 3) +
  geom_line(aes(y = modelRebuilt)) +
  scale_shape_manual(values = c(21,24), labels = c('0°C','5°C')) +
  labs(x = 'body length [mm]',
       y = 'daily growth rate [mm d-1]',
       colour = '') +
  theme(panel.background = element_rect(fill = NA, colour = '#3d3d3d', size = 1.5),
        legend.title = element_text(vjust = 1, size = 24, colour = '#3d3d3d'),
        legend.position = 'bottom',
        legend.text = element_text(size = 24, colour = '#3d3d3d'),
        axis.text = element_text(size = 25, colour = '#3d3d3d'),
        axis.title = element_text(size = 25, colour = '#3d3d3d'),
        strip.text = element_text(size = 25, colour = '#3d3d3d'),
        strip.background = element_rect(fill = NA, colour = NA),
        plot.title = element_text(hjust = 0.5, size = 30, colour = '#3d3d3d'))

ggsave('plots/jagerRavagnan2015QC2.pdf', plot = jagerRavagnanCtrlPlot, width = 10, height = 7.5)
# --------------------------------------------------------------------------------------------------------------- #
# Constable and Kawaguchi (2018)
# load original results, extracted from Figure 4 and 5
constableKawaguchiFig4 <- read_delim("inputData/originalResults/constableKawaguchi2017Fig4.csv", ";", escape_double = FALSE, trim_ws = TRUE) %>% 
  mutate(value = as.numeric(value))
constableKawaguchiFig5 <- read_delim("inputData/originalResults/constableKawaguchi2017Fig5.csv", ";", escape_double = FALSE, trim_ws = TRUE) 
# load model and auxiliary functions
source('functions/ConstableKawaguchi2017AuxiliaryFunctions.R')
# load original environmental data
constableKawaguchiEnv <- read_csv('inputData/originalresults/constableKawaguchiEnvData.csv') %>% 
  mutate(pocConcRed = pocConc/12,
         chlAProx = pocConc/50,
         chlAProxRed = pocConcRed/50)

# Settings of the simulation
dt <- 1
startDay = 1
timeVec <- seq(1,365*6, by = dt)
# ------------------------------------------------------------------------------------------------------------------ #
# initialize state variables
# length
initLength <- 20

# stage, initialized as immature females
initStage <- 2

# body carbon reserves for a healthy individual of initLength
initBodyCReserves <- dryMassToCarbon(dryMass_mg = dryMassAtLength(bodyLength = initLength, stage = 2), tissue = 'bodyReserves')

# gonad carbon reserves are zero
initGonadCReserves <- 0

# init gPhase is "immature"
initGPhase <- 0

# individuals have moulted at t0
initLastMoult <- 0

# intermoult periods determined from initial temperature distribution
initIMP <- determineIMP(temperature = constableKawaguchiEnv$temperature[startDay])

# imp-counter for gPhase-Function
initImpSinceSwitch <- 1

# no eggs produced at t0
initEggsProduced <- 0

results <- matrix(ncol = 9, nrow = length(timeVec))
results[1,] <- c(initLength, initStage, initBodyCReserves, initGonadCReserves, initGPhase, initLastMoult, initIMP, initImpSinceSwitch, initEggsProduced)

# Simulate the model
index <- 1
for(i in timeVec){
  doy <- dayOfYearFunc(t = i, startDay = startDay)
  
  change <- constableKawaguchiGrowth(bodyLength = results[index,1], stage = results[index,2], bodyCReserves = results[index,3], 
                                     gonadCReserves = results[index,4], gPhase = results[index,5], lastMoult = results[index,6], 
                                     moultDay = results[index,7], impSinceSwitch = results[index,8], temperature = constableKawaguchiEnv$temperature[doy],
                                     pocConc = constableKawaguchiEnv$pocConc[doy], photoperiod = constableKawaguchiEnv$dayLength[doy], time = i)
  index <- index + 1
  
  results[index,] <- c(results[index - 1,1] + change[1],
                       change[2],
                       results[index - 1,3] + change[3],
                       results[index - 1,4] + change[4],
                       change[5],
                       change[6],
                       results[index - 1,7] + change[7],
                       change[8],
                       change[9])
  print(i)
}

resultsDf <- as.data.frame(results) %>% 
  mutate(age = (timeVec/365) + 1) %>% 
  mutate(dryMass = carbonToDryMass(carbonMass_mg = V3, tissue = 'bodyReserves'),
         gonadMass = carbonToDryMass(carbonMass_mg = V4, tissue = 'gonadReserves'))
names(resultsDf) <- c('Length (mm)','stage','bodyCReserves','gonadCReserves','gPhase','lastMoult','moultDay','impSinceSwitch','Eggs','age','Dry Mass (mg)','Gonad Dry Mass (mg)')

constableKawaguchiCtrlFig4 <- resultsDf %>% 
  select(1,11,10,9) %>% #age, length, dryMass, eggs) %>% 
  gather(variable, value, -age) %>% 
  ggplot(.,aes(x = age, y =  value)) +
  geom_point(data = constableKawaguchiFig4, colour = '#d43535', size = 1.8) +
  geom_line() +
  #scale_x_continuous(limits = c(1,2)) +
  facet_wrap(~variable, scales = 'free_y') +
  labs(x = 'Age (Years)',
       y = '',
       title = 'constable and kawaguchi (2018) \n quality control Fig 4') +
  theme(panel.background = element_rect(fill = NA, colour = '#3d3d3d', size = 1.5),
        legend.title = element_text(vjust = 1, size = 24, colour = '#3d3d3d'),
        legend.text = element_text(size = 24, colour = '#3d3d3d'),
        axis.text = element_text(size = 25, colour = '#3d3d3d'),
        axis.title = element_text(size = 25, colour = '#3d3d3d'),
        strip.text = element_text(size = 25, colour = '#3d3d3d'),
        strip.background = element_rect(fill = NA, colour = NA),
        plot.title = element_text(hjust = 0.5, size = 30, colour = '#3d3d3d'))
ggsave('plots/constableKawaguchi2017QC1.pdf', plot = constableKawaguchiCtrlFig4, width = 14, height = 6)

constableKawaguchiCtrlFig5 <- resultsDf %>% 
  select(1,12,10) %>% #age, length, dryMass, eggs) %>% 
  filter(between(age,5,6)) %>%
  gather(variable, value, -age) %>% 
  mutate(propOfYear = age - 5) %>% 
  ggplot(.,aes(x = propOfYear, y =  value)) +
  geom_point(data = constableKawaguchiFig5, colour = '#d43535', size = 1.8) +
  geom_line() +
  facet_wrap(~variable, scales = 'free_y') +
  labs(x = 'proportion of year',
       y = '',
       title = 'constable and kawaguchi (2018) \n quality control Fig 5') +
  theme(panel.background = element_rect(fill = NA, colour = '#3d3d3d', size = 1.5),
        legend.title = element_text(vjust = 1, size = 24, colour = '#3d3d3d'),
        legend.text = element_text(size = 24, colour = '#3d3d3d'),
        axis.text = element_text(size = 25, colour = '#3d3d3d'),
        axis.title = element_text(size = 25, colour = '#3d3d3d'),
        strip.text = element_text(size = 25, colour = '#3d3d3d'),
        strip.background = element_rect(fill = NA, colour = NA),
        plot.title = element_text(hjust = 0.5, size = 30, colour = '#3d3d3d'))
ggsave('plots/constableKawaguchi2017QC2.pdf', plot = constableKawaguchiCtrlFig5, width = 14, height = 7.5)

# Follow up - plot manipulated environmental time series
pocDynamics <- constableKawaguchiEnv %>% 
  select(dayOfYear, pocConc,pocConcRed,chlAProx,chlAProxRed) %>% 
  gather(variable,value,-dayOfYear) %>% 
  mutate(class = ifelse(str_detect(variable, pattern = 'poc'), 'POC','Chlorophyll a')) %>% 
  filter(class == 'POC') %>% 
  ggplot(.,aes(x = dayOfYear, y = value, colour = variable)) +
  geom_line(size = 1.4) +
  scale_colour_manual(values = c('#333333','#346acf'), labels = c('original POC','POC/12')) +
  scale_y_continuous(breaks = c(0,100,200,300,400,500,600,700,800)) +
  labs(x = 'day of year',
       y = expression(POC~concentration~'('*mg~m^{-3}*')'),
       title = 'constable and kawaguchi (2018) \n POC dynamics',
       colour = 'legend') +
  theme(panel.background = element_rect(fill = NA, colour = '#3d3d3d', size = 1.5),
        legend.title = element_text(vjust = 1, size = 24, colour = '#3d3d3d'),
        legend.text = element_text(size = 24, colour = '#3d3d3d'),
        axis.text = element_text(size = 25, colour = '#3d3d3d'),
        axis.title = element_text(size = 25, colour = '#3d3d3d'),
        strip.text = element_text(size = 25, colour = '#3d3d3d'),
        strip.background = element_rect(fill = NA, colour = NA),
        legend.position = c(0.75,0.85),
        plot.title = element_text(hjust = 0.5, size = 30, colour = '#3d3d3d'))
ggsave('plots/constableKawaguchi2017QC3.pdf', plot = pocDynamics, width = 9, height = 7.5)


# Follow up - plot manipulated environmental time series
chlADynamics <- constableKawaguchiEnv  %>% 
  select(dayOfYear, pocConc,pocConcRed,chlAProx,chlAProxRed) %>% 
  gather(variable,value,-dayOfYear) %>% 
  mutate(class = ifelse(str_detect(variable, pattern = 'poc'), 'POC','Chlorophyll a')) %>% 
  filter(class == 'Chlorophyll a') %>% 
  ggplot(.,aes(x = dayOfYear, y = value, colour = variable)) +
  geom_line(size = 1.4) +
  scale_colour_manual(values = c('#333333','#346acf'), labels = c('original chlorophyll','reduced chlorophyll')) +
  geom_hline(yintercept = 1, size = 1.2, linetype = '22', colour = '#666666') +
  labs(x = 'day of year',
       y = expression(chlorophyll~a~'('*mg~m^{-3}*')'),
       title = 'constable and kawaguchi (2018) \n chlorophyll a dynamics',
       colour = 'legend') +
  theme(panel.background = element_rect(fill = NA, colour = '#3d3d3d', size = 1.5),
        legend.title = element_text(vjust = 1, size = 24, colour = '#3d3d3d'),
        legend.text = element_text(size = 24, colour = '#3d3d3d'),
        axis.text = element_text(size = 25, colour = '#3d3d3d'),
        axis.title = element_text(size = 25, colour = '#3d3d3d'),
        strip.text = element_text(size = 25, colour = '#3d3d3d'),
        strip.background = element_rect(fill = NA, colour = NA),
        legend.position = c(0.75,0.85),
        plot.title = element_text(hjust = 0.5, size = 30, colour = '#3d3d3d'))
ggsave('plots/constableKawaguchi2017QC4.pdf', plot = chlADynamics, width = 9, height = 7.5)

# repeat simulations with modified environment
results <- matrix(ncol = 9, nrow = length(timeVec))
results[1,] <- c(initLength, initStage, initBodyCReserves, initGonadCReserves, initGPhase, initLastMoult, initIMP, initImpSinceSwitch, initEggsProduced)

# Simulate the model
index <- 1
for(i in timeVec){
  doy <- dayOfYearFunc(t = i, startDay = startDay)
  
  change <- constableKawaguchiGrowth(bodyLength = results[index,1], stage = results[index,2], bodyCReserves = results[index,3], 
                                     gonadCReserves = results[index,4], gPhase = results[index,5], lastMoult = results[index,6], 
                                     moultDay = results[index,7], impSinceSwitch = results[index,8], temperature = constableKawaguchiEnv$temperature[doy],
                                     pocConc = constableKawaguchiEnv$pocConcRed[doy], photoperiod = constableKawaguchiEnv$dayLength[doy], time = i)
  index <- index + 1
  
  results[index,] <- c(results[index - 1,1] + change[1],
                       change[2],
                       results[index - 1,3] + change[3],
                       results[index - 1,4] + change[4],
                       change[5],
                       change[6],
                       results[index - 1,7] + change[7],
                       change[8],
                       change[9])
  print(i)
}

resultsDf <- as.data.frame(results) %>% 
  mutate(age = (timeVec/365) + 1) %>% 
  mutate(dryMass = carbonToDryMass(carbonMass_mg = V3, tissue = 'bodyReserves'),
         gonadMass = carbonToDryMass(carbonMass_mg = V4, tissue = 'gonadReserves'))
names(resultsDf) <- c('Length (mm)','stage','bodyCReserves','gonadCReserves','gPhase','lastMoult','moultDay','impSinceSwitch','Eggs','age','Dry Mass (mg)','Gonad Dry Mass (mg)')

constableKawaguchiCtrlFig4 <- resultsDf %>% 
  select(1,11,10,9) %>% #age, length, dryMass, eggs) %>% 
  gather(variable, value, -age) %>% 
  ggplot(.,aes(x = age, y =  value)) +
  geom_point(data = constableKawaguchiFig4, colour = '#d43535', size = 1.8) +
  geom_line() +
  #scale_x_continuous(limits = c(1,2)) +
  facet_wrap(~variable, scales = 'free_y') +
  labs(x = 'Age (Years)',
       y = '',
       title = 'constable and kawaguchi (2018) \n quality control Fig 4') +
  theme(panel.background = element_rect(fill = NA, colour = '#3d3d3d', size = 1.5),
        legend.title = element_text(vjust = 1, size = 24, colour = '#3d3d3d'),
        legend.text = element_text(size = 24, colour = '#3d3d3d'),
        axis.text = element_text(size = 25, colour = '#3d3d3d'),
        axis.title = element_text(size = 25, colour = '#3d3d3d'),
        strip.text = element_text(size = 25, colour = '#3d3d3d'),
        strip.background = element_rect(fill = NA, colour = NA),
        plot.title = element_text(hjust = 0.5, size = 30, colour = '#3d3d3d'))
ggsave('plots/constableKawaguchi2017QC5.pdf', plot = constableKawaguchiCtrlFig4, width = 14, height = 6)


constableKawaguchiCtrlFig5 <- resultsDf %>% 
  select(1,12,10) %>% #age, length, dryMass, eggs) %>% 
  filter(between(age,5,6)) %>%
  gather(variable, value, -age) %>% 
  mutate(propOfYear = age - 5) %>% 
  ggplot(.,aes(x = propOfYear, y =  value)) +
  geom_point(data = constableKawaguchiFig5, colour = '#d43535', size = 1.8) +
  geom_line() +
  facet_wrap(~variable, scales = 'free_y') +
  labs(x = 'proportion of year',
       y = '',
       title = 'constable and kawaguchi (2018) \n quality control Fig 5') +
  theme(panel.background = element_rect(fill = NA, colour = '#3d3d3d', size = 1.5),
        legend.title = element_text(vjust = 1, size = 24, colour = '#3d3d3d'),
        legend.text = element_text(size = 24, colour = '#3d3d3d'),
        axis.text = element_text(size = 25, colour = '#3d3d3d'),
        axis.title = element_text(size = 25, colour = '#3d3d3d'),
        strip.text = element_text(size = 25, colour = '#3d3d3d'),
        strip.background = element_rect(fill = NA, colour = NA),
        plot.title = element_text(hjust = 0.5, size = 30, colour = '#3d3d3d'))
ggsave('constableKawaguchi2017QC6.pdf', plot = constableKawaguchiCtrlFig5, width = 14, height = 7.5)
