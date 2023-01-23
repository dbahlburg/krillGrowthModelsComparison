# Implementation of the krill growth model presented in Bahlburg et al. (2021), Ecological Modelling
# Full reference:
# Jager T, Ravagnan E (2015) Parameterising a generic model for the dynamic energy budget of Antarctic 
# krill Euphausia superba. Mar Ecol Prog Ser 519:115-128. https://doi.org/10.3354/meps11098
#--------------------------------------------------------------------------------------------------------------------#
# This is the auxiliary feeding function for the bioenergetic growth model of Antarctic Krill
# The function allows to switch between different options:
# JagerRavagnan: The default from the Jager & Ravagnan paper assuming optimum food concentrations
# HollingType2: type 2 functional response. Default parameterisation after Groeneveld et al. (2015)
# scaling in stage 1 krill is 0.28 according to the fB parameter in the JagerRavagnan-Model
feedingFunction <- function(chosenFunction, stage, feedingMode, foodConc = 10, k = 0.3){
  if(chosenFunction == 'JagerRavagnan'){
    ifelse(feedingMode > 0,
           0.28,
           1)
  }
  else if(chosenFunction == 'HollingType2'){
    ifelse(feedingMode > 0,
           0.28,
           foodConc / (foodConc + k) )
  }
}

#This is the temperature correction function based on Arrhenius temperature TA.
#Parameterisation follows estimates from experimental studies where oxygen consumption rates
#were measured as a function of temperature. The unit of the input temperature is Celsius which 
#is converted to Kelvin within the function.
temperatureCorrection <- function(selectedFunction = 'aTemperature', switched = 'on', 
                                  temperature = 0, TA = 8000, T1 = 273.15,
                                  eA = -0.662249,
                                  k = 8.617333262145*10^(-5)){
  if (switched == 'on' & selectedFunction == 'aTemperature') {
    return(exp(TA/T1 - TA/(temperature+273.15)))
  } else if (switched == 'on' & selectedFunction == 'aEnergy'){
    return(exp(27.660421) * exp(eA/(k * (temperature + 273.15)))/(exp(27.660421) * exp(eA/(k * 273.15))))
  } else if(switched == 'off'){
    return(1)
  }
}

# function that triggers reproduction (Spawning occurs every 90d in the spawning window, clutch size simply depends on assimilated energy)
triggerReproductionJager <- function(length, reproState, currentMonth, inputDayOfYear,
                                reproScenario = 'fixedSpawningInterval', weightEgg = 0.028, 
                                spawningWindow = c(10,11,12,1,2,3)){
  if (reproScenario == 'fixedClutchSize'){
    triggerEnergy <<- (length * 150.83 - 3027.244) * weightEgg
    return(reproState > triggerEnergy & currentMonth %in% spawningWindow)
  } else if (reproScenario == 'fixedSpawningInterval'){
    return(inputDayOfYear %in% c(334, 60))
  }
}

# Define a function that returns the dayOfYear at a given simulation timepoint
dayOfYearFunc <- function(t, startDay){
  return(ceiling((t+startDay) - ceiling(((t+startDay)-365)/365) * 365))
}

jagerRavagnanGrowth <- function(inputLength, time, eggBuffer, reproBuffer, temperature, chlorophyllA){
  
  # Define constants
  #-------------------------------------------------------------------------------#
  #Definition of parameters
  #-------------------------------------------------------------------------------#
  # ============================================================================= #
  # Fixed model parameters 
  # ----------------------------------------------------------------------------- #
  # maximum area-specific assimilation rate (scales assimilation of energy)
  # mgDW mm-2 d-1
  maxAssimilation <- 0.045
  # volume-specific somatic maintenance flux
  Jsm<-0.0032 #mgDW mm-3 d-1
  # yield of assimilates on food (carbon based)
  foodConversion<-0.8 # mgC mg-1C
  # yield of structure on assimilates
  growthConversion<-0.8 #mgDW mgDW-1
  # yield of assimilates on structure (starvation)
  structureConversion<-0.8 #mgDW mgDW-1
  # yield of storage buffer on assimilates
  storageConversion<-0.95 #mgDW mgDW-1
  # fraction allocation to soma
  kappaVal <- 0.8
  # physical length at start of juvenile stage
  lengthJuvenile<-11 #mm
  # physical length at start of adult stage
  lengthAdult<-35 #mm
  # temperature effect
  tempSwitch <- 'on'
  # Which scenario should trigger reproduction (Default is "fixedClutchSize")
  reproScenario <- 'fixedSpawningInterval'
  # =============================================================================
  # Conversion factors
  # -----------------------------------------------------------------------------
  # dry weight density
  dv<-0.22 #mgDW mm-3
  # shape correction coefficient (for Lw > 2mm)
  shapeCorrection<-0.2
  # =============================================================================
  #variables that need to be updated in each timestep:
  # update stage information
  stage <- ifelse(eggBuffer > 0,
                  1,
                  ifelse(inputLength > lengthAdult,
                         3,2))
  
  volumetricLength <- inputLength*shapeCorrection # convert actual length to volumetric length
  dayOfYear <- dayOfYearFunc(t = time, startDay = startDay)
  structuralBiomass <- (inputLength*shapeCorrection)^3 * dv
  # =============================================================================
  # calculate assimilated energy (carbon equivalent)
  assEnergy <- maxAssimilation * volumetricLength^2 * 
    feedingFunction(chosenFunction = 'HollingType2', #'JagerRavagnan'
                    stage = stage, 
                    feedingMode = eggBuffer,
                    foodConc = chlorophyllA, 
                    k = 1) *
    temperatureCorrection(switched = tempSwitch, 
                          temperature = temperature) #* foodConversion
  
  #minimum required energy (in carbon equivalent): somatic maintenance
  somMaintenance <- Jsm * volumetricLength^3 *
    temperatureCorrection(switched = tempSwitch, temperature = temperature) 
  
  #calculate potential energy deficit. Because the conversion efficiency of structure into energy is not perfect (0.8),
  #the amount of missing carbon is divided by structureConversion
  energyDeficit <- ifelse((assEnergy * kappaVal - somMaintenance) < 0, 
                          abs(assEnergy * kappaVal - somMaintenance), 0)/structureConversion
  
  #now calculate where the energy deficit is drawn from: either from the reproduction buffer or from the 
  #structural biomass in case the reproduction buffer is depleted
  reproBufferAss <- ifelse(energyDeficit > 0 & reproBuffer >= energyDeficit, 
                           energyDeficit, 
                           ifelse(energyDeficit > 0 & reproBuffer < energyDeficit, 
                                  reproBuffer, 
                                  0)
  )
  structBiomassAss <- ifelse(reproBuffer < energyDeficit, 
                             abs(reproBuffer - energyDeficit),
                             0)
  
  #calculate energy that will be allocated to somatic growth.
  #three scenarios. Each scenario is (de)activated by boolean switches
  #1) enough energy has been assimilated, a fraction covers somMaintenance, the remainder goes into growth (line 1)
  #2) assimilated food (fraction kappa) does not suffice to cover somatic maintenance: energy is drawn from the reproductive branch (growth = 0)
  #3) assimilated food + energy in reproduction buffer do not suffice to cover somatic maintenance: energy is drawn from the structural biomass
  growthEnergy <- ifelse(stage == 1, 
                         assEnergy - somMaintenance, 
                         ifelse(abs(reproBufferAss + structBiomassAss)>0, 0,
                                ((assEnergy * kappaVal) - somMaintenance) * growthConversion))
  
  #Energy that is allocated to the reproduction buffer.
  #There are three possible scenarios (corresponding with the three mentioned above):
  #1) enough energy has been assimilated to cover somatic maintenance with the kappa-fraction. The remainder goes into reproduction
  #2) kappa * assEnergy is not sufficient to cover somatic maintenance: The missing energy is re-allocated from the reproduction branch of energy allocation and subtracted from reproEnergy
  #3) assimilated energy is not sufficient to cover somatic maintenance: No energy goes into reproduction
  reproEnergy <- (1 - kappaVal) * assEnergy * growthConversion
  
  #update the state variables:
  #Structural biomass
  changeStructuralBiomass <- growthEnergy - structBiomassAss #krill$structuralBiomass[index-1] + 
  
  #Egg buffer (for embryos).
  #There are three phases of buffer depletion:
  #1) buffer can cover somatic maintenance costs: energy is drawn from buffer
  #2) buffer still has energy but not sufficient to cover somatic maintenance: remaining energy is utilized
  #3) buffer is depleted: it remains zero
  changeEggBuffer <- ifelse(eggBuffer - somMaintenance > 0,
                            -1 * (assEnergy/storageConversion),
                            0)
  
  #Reproduction Buffer
  changeReproBuffer <- ifelse(stage == 3 & triggerReproductionJager(inputDayOfYear = dayOfYear,
                                                               reproScenario = reproScenario) == FALSE, 
                              reproEnergy - reproBufferAss,
                              ifelse(stage == 3 & triggerReproductionJager(inputDayOfYear = dayOfYear,
                                                                      reproScenario = reproScenario) == TRUE,
                                     -reproBuffer,0))
  # released Eggs
  releasedEggs <- ifelse(stage == 3 & triggerReproductionJager(inputDayOfYear = dayOfYear,
                                                          reproScenario = reproScenario) == TRUE,
                         round(reproBuffer/0.028),
                         0)
  # change in length
  changeLength <- ((((structuralBiomass + changeStructuralBiomass)/dv)^(1/3))/shapeCorrection) -  ((((structuralBiomass)/dv)^(1/3))/shapeCorrection)
  
  # all processes are set to zero under sea ice cover
  changeLength <- ifelse(is.na(chlorophyllA), 0, changeLength)
  changeEggBuffer <- ifelse(is.na(chlorophyllA), 0, changeEggBuffer)
  changeReproBuffer <- ifelse(is.na(chlorophyllA), 0, changeReproBuffer)
  releasedEggs <- ifelse(is.na(chlorophyllA), 0, releasedEggs)
  return(c(changeLength, changeEggBuffer, changeReproBuffer, releasedEggs))
}

# Runge-Kutta 4th order solver
rk4step <- function(state, time, inputEggBuffer, inputReproBuffer, inputTemperature, inputChla){
  
  # inputLength, time, eggBuffer, reproBuffer, temperature, photoperiod, chlorophyllA
  dt <- 1
  k1 <- jagerRavagnanGrowth(inputLength = state, time = time, eggBuffer = inputEggBuffer, reproBuffer = inputReproBuffer, temperature = inputTemperature, chlorophyllA = inputChla)
  k2 <- jagerRavagnanGrowth(inputLength = state+0.5*k1[1]*dt, time = time + 0.5 * dt, eggBuffer = inputEggBuffer, reproBuffer = inputReproBuffer, temperature = inputTemperature, chlorophyllA = inputChla)
  k3 <- jagerRavagnanGrowth(inputLength = state+0.5*k2[1]*dt, time = time + 0.5 * dt, eggBuffer = inputEggBuffer, reproBuffer = inputReproBuffer, temperature = inputTemperature, chlorophyllA = inputChla)
  k4 <- jagerRavagnanGrowth(inputLength = state+k3[1]*dt, time = time + dt, eggBuffer = inputEggBuffer, reproBuffer = inputReproBuffer, temperature = inputTemperature, chlorophyllA = inputChla)
  
  return(dt * (k1 + 2*k2 + 2*k3 + k4) /6)
}

