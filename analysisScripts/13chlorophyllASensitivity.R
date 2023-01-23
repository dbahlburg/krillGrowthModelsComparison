# Calculate chlorophyll a concentrations required by each model to predict fixed growth rates
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
lookUpTableTimeVariables <- read.csv('inputData/HofmannLascara2000/timeDependentVariables.csv')
lookUpTableSizeClassTraits <- read.csv('inputData/HofmannLascara2000/sizeClassTraits.csv')
# ------------------------------------------------------------------------------------------------------- #
# Create a grid of temperature and chlorophyll a gradients
# Some models require POC concentrations as input - chlorophyll a will be converted to POC using 
# a conversion factor of 50 (Hofmann and Lascara 2000)
growthPlain <- expand_grid(chla = seq(0,50, length.out = 1000),
                           temp = 1) %>% 
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
  mutate(
    #Fach et al. (2002) DGR
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
    `Atkinson et al. (2006)` = atkinsonFunc(inputFood = chla, inputTemperature = temp, inputStage = 3, inputLength = bodyLength)[1],
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
    `Tarling et al. (2006)` = TarlingEtAl2006Model(inputLength = bodyLength, 
                                                   inputTemperature = temp, 
                                                   inputChla = chla, 
                                                   inputStage = 3,
                                                   time = TarlingEtAl2006IMPRounded(bodyLength = bodyLength,
                                                                                    temperature = temp),
                                                   oldMoultDay = moultDay,
                                                   temperatureHistory = temp,
                                                   chlorophyllHistory = chla,
                                                   moultDay = TarlingEtAl2006IMPRounded(bodyLength = bodyLength,
                                                                                        temperature = temp))[1]/TarlingEtAl2006IMPRounded(bodyLength = bodyLength,
                                                                                                                                          temperature = temp),
    #Wiedenmann et al. (2008)
    `Wiedenmann et al. (2008)` = WiedenmannEtAl2008Model(inputLength = bodyLength, 
                                                         inputTemperature = temp, 
                                                         inputChla = chla, 
                                                         time = round(exp(3.5371 - 0.5358 * log(temp + 2))),
                                                         oldMoultDay = moultDay,
                                                         temperatureHistory = temp,
                                                         chlorophyllHistory = chla,
                                                         moultDay = round(exp(3.5371 - 0.5358 * log(temp + 2))))[1]/round(exp(3.5371 - 0.5358 * log(temp + 2))))


growthPlainDGR <- growthPlainDGR %>%
  select(-hofmannChangeMass, -fachChangeMass) %>% 
  gather(model, dgr, -chla, -temp, -poc) 

# extract environmental conditions for each model where distance to the wanted DGR is minimal
chlaGrowth <- growthPlainDGR %>% 
  mutate(dgr01 = abs(dgr - 0.1),
         dgr015 = abs(dgr - 0.15),
         dgr02 = abs(dgr - 0.2)) %>% 
  group_by(model) %>% 
  filter(dgr01 == min(dgr01) | dgr015 == min(dgr015) | dgr02 == min(dgr02)) %>% 
  select(model, chla, dgr01, dgr015, dgr02) %>% 
  gather(growthRate, distance, -model, -chla) %>% 
  group_by(model, growthRate) %>% 
  filter(distance == min(distance) & growthRate != 'dgr02') %>% 
  mutate(growthRate = ifelse(growthRate == 'dgr01', 'dgr = 0.1mm/d',
                             'dgr = 0.15mm/d'),
         growthRate = factor(growthRate, levels = c('dgr = 0.1mm/d','dgr = 0.15mm/d')),
         model = factor(model))

sensitivityPlot <- chlaGrowth %>% 
  mutate(model = fct_relevel(model, levels(chlaGrowth$model)[c(6,2,5,4,3,7,9,1,8)])) %>% 
  ggplot(.,aes(x = chla, y = model)) +
  geom_col(fill = '#458255') +
  facet_wrap(~growthRate, scales = 'free_x') +
  labs(x = expression(chlorophyll~a~concentration~mg~m^{-3}),
       y = '') +
  theme(panel.background = element_rect(fill = NA, colour = '#2b2b2b'),
        panel.grid = element_blank(),
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 18),
        strip.background = element_rect(fill = NA),
        strip.text = element_text(size = 18))

ggsave('plots/chlaSensitivity.pdf', plot = sensitivityPlot, width = 13, height = 5)
