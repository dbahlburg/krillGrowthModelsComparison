# Implementation of the krill growth model presented in Ryabov et al. (2017), Nature EcoEvo
# Full reference:
# Ryabov, A., de Roos, A., Meyer, B. et al. Competition-induced starvation drives large-scale population cycles in 
# Antarctic krill. Nat Ecol Evol 1, 0177 (2017). https://doi.org/10.1038/s41559-017-0177
#--------------------------------------------------------------------------------------------------------------------#
# function that defines stage
determineStage <- function(age){
  # ifelse(age < 30, return('E'),
  #        ifelse(between(age, 30, 365), return('L'),
  #               ifelse(between(age,365, 730),return('J'),
  #                      ifelse(age > 730, return('A'), return(NA))
  #                      )
  #               )
  #        )
  stages <- ifelse(age < 30, 'E',
                   ifelse(between(age, 30, 365), 'L',
                          ifelse(between(age,365, 730),'J',
                                 ifelse(age > 730, 'A', NA)
                          )
                   )
  )
  return(stages)
}

# function that converts dry mass to length
lengthAtMassRyabov <- function(inputMass){
  #mass type: dry mass
  return((inputMass > 0.0277)*exp(0.5 * (-1.8050 + sqrt(1.8050 * 1.8050  + 4 * 0.2380 * log(inputMass/0.0058)))/0.2380) + 2.2 * (inputMass <= 0.0277))
}

# function that converts length to dry mass
massAtLengthRyabov <- function(inputLength){
  #mass type: dry mass
  return((inputLength > 2.2) * 0.0058 * inputLength^(1.8050 + 0.2380 * log(inputLength)) + (inputLength <= 2.2) * 0.0277)
}

# function that scales metabolism according to day of year for juveniles and adults
metabolicRate <- function(stage, dayOfYear){
  #variation in activity for adults and juveniles
  activityMin <- 0.2   #of the max ingestion rate during the period of inactivity
  activityBestDay <- -30/365  #max activity day, 1st of December
  activityDuration <- 200  #duration of the high activity period, approx 200 days
  activityDurationConst <- (pi/2 * activityDuration/365)^2  
  
  metRate <- ifelse(stage %in% c('J','A'),
              activityMin + (1 - activityMin) * (exp((cos(2 * pi * (dayOfYear/365 - activityBestDay)) - 1)/activityDurationConst)),
              1)
  return(metRate)
}

# function that returns the fraction of body weight that is needed for daily maintenance
maintenanceCostsFrac <- function(stage){
  maintenanceCosts <- ifelse(stage == 'E', 0,
                             ifelse(stage == 'L', 0.01,
                                    ifelse(stage %in% c('J','A'), 0.003, NA
                                           )
                                    )
                             )
  return(maintenanceCosts)
}

ingestionFunc <- function(inputChla, inputIceAlgae, foodSource, stage){
  
  # defined constants
  I_max <- 0.04 #(mg Carbon /(mg dry weight * day))   %10% of body carbon per day  -> 10/100 * Carb2DryW= 0.04 (mg Carbon /(mg dry weight * day));
  H_S <- 0.5 #25 mg C/m^3 == 0.5 mg Chla/m^3 assuming carb:Chl = 50;
  
  
  # Calculate ingestion rates larvae
  iceAlgLarv <- inputIceAlgae * inputIceAlgae / (inputChla + inputIceAlgae)
  chlLarv     <- inputChla * inputChla / (inputChla + inputIceAlgae)
  
  # return results
  ingestions <- ifelse(stage == 'E', 0,
                       ifelse(stage == 'L' & foodSource == 'iceAlgae', 
                              I_max * iceAlgLarv/(chlLarv + iceAlgLarv + H_S),
                              ifelse(stage == 'L' & foodSource == 'waterColumn',
                                     I_max * chlLarv/(chlLarv + iceAlgLarv + H_S),
                                     ifelse(stage %in% c('J','A') & foodSource == 'iceAlgae',
                                            0,
                                            ifelse(stage %in% c('J','A') & foodSource == 'waterColumn',
                                                   I_max*inputChla/(inputChla + H_S),
                                                   NA
                                                   )
                                            )
                                     )
                              )
                       )
  return(ingestions)
}

# Define a function that returns the dayOfYear at a given simulation timepoint
dayOfYearFunc <- function(t, startDay){
  return(ceiling((t+startDay) - ceiling(((t+startDay)-365)/365) * 365))
}

# function that determines how much energy is allocated to growth 
growthAlloc <- function(inputLength){
  return(1/(1+exp(0.2 * (inputLength-43))))
}

# the growth model
ryabovEtAlGrowth <- function(y, time, chlA, iceAlgae, age, parameters = c('startDay' = 1, 'assimEff' = 1.04)){
  # parameter assimEff stands for assimilation efficiency: [mg DW/mg C consumed] assimilation efficiency, note they consume not only carbon
  
  # determine current day of year
  currentDOY <- time #dayOfYearFunc(t = time, startDay = parameters['startDay'])
  
  # determine stage
  krillStage <- determineStage(age = age)
  
  # get dry weight
  dryMass <- massAtLengthRyabov(inputLength = y)
  
  # energy assimilation
  assEnergy <- (parameters['assimEff'] * metabolicRate(stage = krillStage, dayOfYear = currentDOY) * 
                  (ingestionFunc(inputChla = chlA, inputIceAlgae = iceAlgae, foodSource = 'iceAlgae', stage = krillStage) + 
                     ingestionFunc(inputChla = chlA, inputIceAlgae = iceAlgae, foodSource = 'waterColumn', stage = krillStage)) - 
                  maintenanceCostsFrac(stage = krillStage)) * dryMass
  
  # determine mass-specific growth rate
  dydt_Mass <- ifelse(assEnergy < 0, assEnergy,
                      growthAlloc(inputLength = y) * assEnergy)
  
  # convert into length-specific growth rate
  dydt_Length <- lengthAtMassRyabov(inputMass = (dryMass + dydt_Mass)) - lengthAtMassRyabov(inputMass = dryMass)
  
  #return growth rate
  # When sea ice cover, return zero growth
  chlaNA <- which(is.na(chlA))
  dydt_Length[chlaNA] <- 0
  return(as.numeric(dydt_Length))
}

# Runge Kutta 4th-order solver of the model
rk4stepRyabov <- function(state, time, inputChla, inputIceAlgae, inputAge){
  
  dt <- 1
  k1 <- ryabovEtAlGrowth(y = state, time = time, chlA = inputChla, iceAlgae = inputIceAlgae, age = inputAge)
  k2 <- ryabovEtAlGrowth(y = state+0.5*k1*dt, time = time + 0.5 * dt, chlA = inputChla, iceAlgae = inputIceAlgae, age = inputAge)
  k3 <- ryabovEtAlGrowth(y = state+0.5*k2*dt, time = time + 0.5 * dt, chlA = inputChla, iceAlgae = inputIceAlgae, age = inputAge)
  k4 <- ryabovEtAlGrowth(y = state+k3*dt, time = time + dt, chlA = inputChla, iceAlgae = inputIceAlgae, age = inputAge)

  return(dt * (k1 + 2*k2 + 2*k3 + k4) /6)
}






