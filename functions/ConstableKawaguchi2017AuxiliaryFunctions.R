# Define a function that returns the dayOfYear at a given simulation timepoint
dayOfYearFunc <- function(t, startDay){
  return(ceiling((t+startDay) - ceiling(((t+startDay)-365)/365) * 365))
}

# Equations used in the following lines: eq. 4.5 (d_1 and d_2 corrected), eq. 4.6, eq. 4.7, eq. 4.8 (beta corrected according to Hofmann and Lascara)
# filtrationRate is converted from m3 to liters.
filtrationRateLiter <- function(dryWeight_mg){
  modelDataFiltration <-
    tibble(
      modDryMass = c(seq(0.01, 26, length.out = 30), seq(84, 100, length.out = 30)),
      filtrationRate = ifelse(
        between(modDryMass, 0, 26),
        0.00085 * modDryMass^0.825,
        0.00343 * modDryMass^0.514))
  filtrationLoess <- loess(filtrationRate ~ modDryMass, data = modelDataFiltration)
  
  conversionTable <- tibble(dryMass = dryWeight_mg) %>% 
    mutate(filtrationRateWaterColumn = ifelse(dryMass<26,  0.00085 * dryWeight_mg^0.825,
                                              ifelse(between(dryMass, 26, 84), as.numeric(predict(filtrationLoess, data.frame(modDryMass = dryMass))), 
                                                     ifelse(dryMass > 84, 0.00343 * dryWeight_mg^0.514, NA))))
  
  filtrationRateWaterColumn <- as.numeric(conversionTable$filtrationRateWaterColumn)
  return(filtrationRateWaterColumn)
}

# eq. 4.1
# Definition stages: 
# 1: male
# 2: immature female
# 3: mature female
dryMassAtLength <- function(bodyLength, stage){
  ifelse(stage == 1, return(1000 * (0.195 * (21.38e-6 * bodyLength^2.76) - 0.014)),
         ifelse(stage == 2, return(1000 * (0.251 * (30.2e-6 * bodyLength^2.62) - 0.013)),
                ifelse(stage == 3, return(1000 * (0.275 * (9.77e-6 * bodyLength^2.98) - 0.014)), return(NA))))
}

# feeding and metabolic activity from Hofmann and Lascara
detFeedingActivity <- function(dailyRatio){
  return(ifelse(dailyRatio <=10, 0.1 * dailyRatio, 1))
}

# eq. 4.1 inverted
# Definition stages: 
# 1: male
# 2: immature female
# 3: mature female
lengthAtDryMass <- function(dryMass_mg, stage){
  #((a * (dryMass_mg/1000 - d2))/d1)^(1/b)
  
  ifelse(stage == 1, return((((dryMass_mg/1000) + 0.014)/(21.38e-6 * 0.195))^(1/2.76)),
         ifelse(stage == 2, return(((dryMass_mg/1000 + 0.013)/(30.2e-6 * 0.251))^(1/2.62)),
                ifelse(stage == 3, return(((dryMass_mg/1000 + 0.014)/(9.77e-6 * 0.275))^(1/2.98)), return(NA))))
}

# Equation 2.4
carbonResp_perDay <- function(temperature, dryWeight_mg){
  
  return(1e-3 * 24 * 0.5357 * (10^(0.02438 * temperature - 0.1838) * dryWeight_mg^(-0.0109 * temperature + 0.8918)))
  #return(10^((-0.01089 * temperature + 0.8918) * log10(dryWeight_mg) + (0.02348 * temperature - 0.1838))* 24 * 0.5357)# * 1e-3 )
}

# Eq. 4.3 and 4.4
dryMassToCarbon <- function(dryMass_mg, tissue){
  ifelse(tissue == 'bodyReserves', return(0.366 * dryMass_mg^1.037),
         ifelse(tissue == 'gonadReserves', return(dryMass_mg * 0.551), NA))
}

# Eq. 4.3 inverted
carbonToDryMass <- function(carbonMass_mg, tissue){
  ifelse(tissue == 'bodyReserves', return((carbonMass_mg/0.366)^(1/1.037)),
         ifelse(tissue == 'gonadReserves', return(carbonMass_mg/0.551), NA))
}

# max gonadal mass (equation 2.7 derived from combining equations 4.2 and 4.4 and Morris 1988)
# WS4, WS1, according to my interpretation stand for the weight of the head+carapax and the weight of the total krill
# Therefore, max gonad mass is defined as the weight of the abdomen + telson in Morris 1988
# maxGonadMass <- function(bodyLength){
#   
#   gonadDryMass <- 
#   
# }

rk4step <- function(dt, state, time, inputFunc){
  k1 <- inputFunc(y = state, time = time, parameters = parameters)
  k2 <- inputFunc(y = state+0.5*k1*dt, time = time + 0.5 * dt, parameters = parameters)
  k3 <- inputFunc(y = state+0.5*k2*dt, time = time + 0.5 * dt, parameters = parameters)
  k4 <- inputFunc(y = state+k3*dt, time = time + dt, parameters = parameters)
  
  return(dt * (k1 + 2*k2 + 2*k3 + k4) /6)
}

# determine reproductive state
# one of: reproductive, regression, postregression
# rules: During reproductive period:
# 1. Spawning occurs at every second moult if M_Gonads > M_crit
# 2. If (1.) is not true, M_Gonads is transferred to next IMP
# 3. Regression phase starts when: 1. more than 4 moults since last spawning
#                                  2. at least two moults since last spawning and day length <5 hours
#                                  3. day length <2 hours
# 4. Regression phase ends after 5 moult cycles if initiated by 1. or 2. or after 3 moult cycles when initiated by 3.
# 5. Post-regression phase follows regression phase and lasts for 5 moult cycles, then reproductive period is initiated
# Definitions reproductive states:
# 0: immature
# 1: reproductive
# 2: regression1
# 3: regression2
# 4: postregression
determineGPhase <- function(gPhase, stage, impSinceSwitch, dayLength){
  
  newgPhaseRepro <- function(impSinceSwitch, dayLength){
    newGPhase <- ifelse(impSinceSwitch > 4, 2,
                        ifelse(impSinceSwitch >= 2 & dayLength < 5, 2,
                               ifelse(dayLength < 2, 3, 1)
                        )
    ) 
    return(newGPhase)
  }
  
  newgPhaseReg <- function(gPhase, impSinceSwitch){
    newGPhase <- ifelse(gPhase == 2 & impSinceSwitch >= 5, 4,
                        ifelse(gPhase == 3 & impSinceSwitch >= 3, 4, gPhase)
    )
    return(newGPhase)
  }
  
  newgPhasePostReg <- function(impSinceSwitch){
    newGPhase <- ifelse(impSinceSwitch >= 5, 1, 4)
    return(newGPhase)
  }
  
  newgPhase <- ifelse(gPhase == 1, newgPhaseRepro(impSinceSwitch = impSinceSwitch, dayLength = dayLength), 
                        ifelse(gPhase %in% c(2,3), newgPhaseReg(gPhase = gPhase, impSinceSwitch = impSinceSwitch),
                               ifelse(gPhase == 4, newgPhasePostReg(impSinceSwitch = impSinceSwitch), 0
                               )
                        )
  )
  
  return(newgPhase)
}

# IMP-model from Kawaguchi et al. (2006) - Equation 1.5
determineIMP <- function(temperature){
  return(round(exp(3.5371 - 0.5358 * log(temperature + 2))))
}


# Exoskeleton carbon mass (Equation 4.9)
cMassExo <- function(dryMass_mg){
  return(0.238 * 0.075 * dryMass_mg)
}


# function for sexual maturity
# stage: 0 male, 1 immature female, 2 mature female
# this function is added and not included in the original paper (stage switch between immature females and mature females
# is not described)
sexMaturity <- function(inputLength, stage){
  ifelse(stage == 1, return(1),
         ifelse(stage == 2 & inputLength < 35, return(2), return(3)))
}

# function returning number of eggs released during spawning
eggsProducedFunc <- function(gPhase, inputGonadCReserves, input_dGMax, impSinceSwitch){
  gonadDryMass <- carbonToDryMass(carbonMass_mg = inputGonadCReserves, tissue = 'gonadReserves')
  eggsProduced <- ifelse(gPhase == 1 & inputGonadCReserves > (0.1 * input_dGMax) & impSinceSwitch >=2 , 
                         round((gonadDryMass * 25.8)/0.6824), 0)
  return(eggsProduced)
}

# metabolic activity function after Hofmann and Lascara (2000)
detMetabolicActivity <- function(dayOfYear = NA, dayLength = NA){
    ifelse(between(dayOfYear, 1, 105), return(0),
           ifelse(between(dayOfYear, 105, 150), return(-0.5/45 * (dayOfYear-105)),
                  ifelse(between(dayOfYear, 150, 250), return(-0.5),
                         ifelse(between(dayOfYear, 250, 295), return(-0.5 + 0.5/45 * (dayOfYear-250)),
                                ifelse(between(dayOfYear, 295, 365), return(0), return(NA)
                                       )
                                )
                         )
                  )
           )
}


moultFunc <- function(bodyLength, moultDay, lastMoult, bodyCReserves, moultCCosts, stage, gPhase, gonadCReserves, dG_max, impSinceSwitch, photoperiod, temperature, pocConc){
  
  # determine maxiumum growth and shrinkage. Shrinkage has to be multiplied with (moultDay - lastMoult) since its 
  # unit is mm d-1
  maxGrowth <- 8.63 - 0.14 * bodyLength
  maxShrinkage <- 1e-3 - 1.418e-5 * bodyLength * (moultDay - lastMoult)
  
  # new length is determined by the body reserves at the day of moulting
  # 1. When body reserves are lower than the weight of a healthy individual at the given length, shrinkage occurs
  # 2. When the body reserves are higher than the weight of a healthy individual at the given length, growth occurs
  # In case 2. is true, growth is capped by the maximum growth increment. When the maximum growth is not reached, 
  # growth is scaled according to the predicted length at the given dry mass
  # In case 1. is true, shrinkage occurs based on the weight loss since the last moult. It is capped by the maximum
  # shrinkage. When the dry weight is less than that of a healthy krill with the new length, thin krill will be the outcome
  # 1. subtract the carbon costs of moulting from body reserves
  dL <- lengthAtDryMass(carbonToDryMass(carbonMass_mg = (bodyCReserves - moultCCosts), tissue = 'bodyReserves'), stage = 2) - bodyLength
  
  
  # Update the state variables:
  # 1. rate of change body length
  # 2. rate of change body c reserves
  # 3. new stage (1 - male, 2 - immature female, 3 - mature female)
  # 4. number of produced eggs
  # 5. rate of change gonad c reserves
  # 6. update impSinceSwitch
  # 7. update lastMoult (last timepoint for moulting - alias the current timepoint)
  # 8. determine newMoult (new timepoint for moulting)
  
  dL <- ifelse(dL > maxGrowth, maxGrowth,
               ifelse(dL < maxShrinkage, maxShrinkage, dL)
  )
  # body C-reserves are decreased by the C-costs of moulting
  dB <- -1 * moultCCosts
  
  # change stage- and reproduction specific dynamics
  # update reproductive state
  oldStage <- stage
  newStage <- sexMaturity(inputLength = bodyLength, stage = stage)
  
  # check reproduction
  eggsProduced <- eggsProducedFunc(gPhase = gPhase, inputGonadCReserves = gonadCReserves, 
                                   input_dGMax = dG_max, impSinceSwitch = impSinceSwitch)
  
  # change in gonad reserves is zero, when no spawning occurs. Otherwise, they do not change
  dG <- ifelse(eggsProduced > 0, -1 * gonadCReserves, 0)
  
  # if spawning occurs, impSinceSwitch is set to 1 
  newImpSinceSwitch <- ifelse(eggsProduced > 0, 1, impSinceSwitch)
  
  # determine gonadal phase (when individual just reached sexual maturity, it is set to 1, otherwise the rep-function is called)
  oldGPhase <- gPhase
  newGPhase <- ifelse(oldStage == 2 & newStage == 3, 1, determineGPhase(gPhase = gPhase, impSinceSwitch = impSinceSwitch, dayLength = photoperiod))
  
  # when stage switch occurred, impSinceSwitch is set to 1 (the logic is that the next moulting will be the first since the 
  # the stage switch)
  newImpSinceSwitch <- ifelse(eggsProduced > 0, newImpSinceSwitch,
                              ifelse(newGPhase != oldGPhase, 1, impSinceSwitch + 1))
  
  # determine new imp:
  newLastMoult <- moultDay
  newIMP <- determineIMP(temperature)
  
  # deal with na's
  dL <- ifelse(is.na(pocConc),0,dL)
  dB <- ifelse(is.na(pocConc),0,dB)
  dG <- ifelse(is.na(pocConc),0,dG)
  newStage <- ifelse(is.na(pocConc),stage,newStage)
  eggsProduced <- ifelse(is.na(pocConc),0,eggsProduced)
  newImpSinceSwitch <- ifelse(is.na(pocConc),impSinceSwitch,newImpSinceSwitch)
  newLastMoult <- ifelse(is.na(pocConc),lastMoult,newLastMoult)
  newIMP <- ifelse(is.na(pocConc), newIMP, newIMP)
  newGPhase <- ifelse(is.na(pocConc),gPhase,newGPhase)
  return(list(dL = dL, newStage = newStage, dB = dB, dG = dG, newGPhase = newGPhase, newLastMoult = newLastMoult, 
           newIMP = newIMP, newImpSinceSwitch = newImpSinceSwitch, eggsProduced = eggsProduced))
}


# Constable and Kawaguchi model.
# at each timestep the model evaluates
# 1. the energy dynamics (intake and demand, investment into growth and reproduction)
# 2. moult dynamics:
#   2.1 Does moulting take place?
#   2.2 What is the body condition of the krill - determines growth increment
#   2.3 Does reproduction take place?
#   2.4 Update reproductive status
#   2.5 Determine new moult day
# Definition of state variables:
# bodyLength: body length [mm]
# stage: development stage (1-3)
# bodyReserves: body reserves [mg]
# gonadReserves: gonad reserves [mg]
# gPhase: reproductive state (0-4)
# moultDay: day at which moulting occurs
# impSinceSwitch: how many imps have past since last time reproductive state switched or spawning occurred
constableKawaguchiGrowth <- function(bodyLength, stage, bodyCReserves, gonadCReserves, gPhase, lastMoult, moultDay, impSinceSwitch, temperature, pocConc, photoperiod, time){
  
  dayOfYear <- dayOfYearFunc(t = time, startDay = startDay)
  assimEfficiency <- 0.8
  actualBodyMass_mg <- carbonToDryMass(carbonMass_mg = bodyCReserves, tissue = 'bodyReserves') + 
    carbonToDryMass(carbonMass_mg = gonadCReserves, tissue = 'gonadReserves')
  # ---------------------------------------------------------------------------------------------------------------- #
  # DAILY BUSINESS OF THE KRILL
  # ---------------------------------------------------------------------------------------------------------------- #
  # 1. assimilated carbon in mg
  assimilatedCarbon <- filtrationRateLiter(dryWeight_mg = dryMassAtLength(bodyLength = bodyLength,
                                                                          stage = stage)) * pocConc * assimEfficiency
  
  # Bioenergetics of the model (Figure 2)
  # POC is ingested and assimilated at a certain rate which leads to a certain carbon budget
  # The budget is calculated as A' = (Body Mass + Gonad Mass - Healthy Body Mass - Moulting Costs) + (Assimilated Food - Metabolic costs)
  # When A' < 0, the deficit is drawn from the body reserves (negative growth)
  # When A' > 0, the energy can be invested into growth and reproduction.
  # When A' > 0, it is checked whether A' exceeds the maximum amount of reserves that can be allocated to growth and gonads (dB_max, dG_max)
  # If this is true, A' is set to dB_max + dG_max (and now referred to as A'')
  # If this is not true, A' remains as it is (and now referred to as A'')
  # ---------------------------------------------------------------------------------------------------------------- #
  # ## energy investment into growth ##
  # 1. From the available pool of energy A'', a certain fraction (1-P) is allocated to growth
  # 2. It is checked whether A'' * (1-P) > dB_max
  # 3. If this is true, A'' is set to dB_max and the excess reserves (dB_exc) are being allocated to reproduction
  # 4. If A'' * (1-P) < dB_max, it is checked whether there are excess gonad reserves that have been allocated to growth. 
  # 5. is A'' * (1-P) + dG_exc <= dB_max? (when excess gonad reserves are being added, is the energy still below dB_max)
  # 6. If it is true, A'' * (1-P) + dG_exc is invested into body reserves
  # 7. If dB_max was exceeded with dG_exc being added, energy investment into body reserves is set to 0
  # ---------------------------------------------------------------------------------------------------------------- #
  # ## energy investment into gonads ##
  # 1. From the available pool of energy A'', a certain fraction (P) is allocated to the gonads
  # 2. It is checked whether A'' * P > dG_max
  # 3. If this is true, A'' is set to dG_max and the excess reserves (dG_exc) are being allocated to body reserves
  # 4. If A'' * P < dG_max, it is checked whether there are excess body reserves (dB_exc) that have been allocated to the gonads 
  # 5. is A'' * P + dB_exc <= dG_max (when excess gonad reserves are being added, is the energy still below dG_max)
  # 6. If it is true, A'' * P + dB_exc is invested into body reserves
  # 7. If dG_max was exceeded with dB_exc being added, energy investment into gonads is set to 0
  # ---------------------------------------------------------------------------------------------------------------- #
  # healthy carbon mass is calculated as the dry mass of a stage2-individual with the given length to exclude gonad mass
  healthyCMass <- dryMassToCarbon(dryMass_mg = dryMassAtLength(bodyLength = bodyLength, stage = 2), tissue = 'bodyReserves')
  moultCCosts <- cMassExo(dryMass_mg = dryMassAtLength(bodyLength = bodyLength, stage = stage))
  
  dailyRation <- assimilatedCarbon/actualBodyMass_mg * 100
  dailyRation <- ifelse(dailyRation < 1.5, 0, dailyRation)
  CMetabolism <- carbonResp_perDay(temperature = temperature, dryWeight_mg = actualBodyMass_mg) * (1 + detMetabolicActivity(dayOfYear = dayOfYear))# + detFeedingActivity(dailyRatio = dailyRation))
  # CMetabolism <- lookUpTableSizeClassTraits$carbonRespiration[match(round(bodyLength*4)/4, sizeClassesVec)] *
  #   (1 + metabolicActivity + detFeedingActivity(dailyRatio = assimilatedCarbon/bodyCReserves * 100))
  # dB_max is defined as the body body reserves necessary to achieve maximum growth increment and ending up as a healthy krill post-moulting.
  # Daily carbon allocation to growth is evaluated against dB_max in such way that the already present body reserves + the allocated carbon
  # must not exceed the value of dB_max
  # current length + maximum growth increment
  maxGrowthIncrement <- 8.63 - 0.14 * bodyLength
  dB_max <- dryMassToCarbon(dryMass_mg = dryMassAtLength(bodyLength = (bodyLength + maxGrowthIncrement), stage = 2), tissue = 'bodyReserves') + moultCCosts
  P <- ifelse(gPhase == 1, 0.95,0)#0.95,0) #0.8, 0) not clear if allocation constant is 0.8 or 0.95 (both are mentioned in the manuscript)
  
  # dG_max is defined as the maximum gonadal carbon mass. Daily carbon allocation to gonads is evaluated against dG_max in such way that
  # the already present gonad reserves + allocated carbon must not exceed the value of dG_max
  dG_max <- dryMassToCarbon(dryMass_mg = dryMassAtLength(bodyLength = bodyLength, stage = stage), tissue = 'gonadReserves') -
    dryMassToCarbon(dryMass_mg = dryMassAtLength(bodyLength = bodyLength, stage = 2), tissue = 'gonadReserves')
  
  # calculate available energy budget in carbon units
  #A_1 <- (bodyCReserves + gonadCReserves - healthyCMass - moultCCosts) + (assimilatedCarbon - CMetabolism)
  A_1 <- assimilatedCarbon - CMetabolism
  A_2 <- ifelse(A_1 < (dB_max + dG_max), A_1, (dB_max + dG_max))
  # ---------------------------------------------------------------------------------------------------------------- #
  # energy investment into growth and reserves according to the rules above
  # in contrast to the manuscript, dB and dG here refer to the rate of change in body and gonad carbon reserves
  # they are set to zero when dB_max or dG_max is reached, respectively
  dB <- ifelse((bodyCReserves + (A_2 * (1-P))) <= dB_max, A_2 * (1-P), dB_max - bodyCReserves)
  dG <- ifelse((gonadCReserves + (A_2 * P)) <= dG_max, A_2 * P, dG_max - gonadCReserves)
  
  dB_exc <- ifelse(dB == 0, (bodyCReserves + (A_2 * (1-P))) - dB_max, 0)
  dG_exc <- ifelse(P>0 & dG == 0, (gonadCReserves + (A_2 * P)) - dG_max, 0)
  
  dB <- ifelse((bodyCReserves + (dB + dG_exc)) < dB_max, dB + dG_exc, dB_max - bodyCReserves)
  dG <- ifelse((gonadCReserves + (dG + dB_exc)) < dG_max, dG + dB_exc,  dG_max - gonadCReserves)
  
  # If there is energy deficit for body reserves and gonad reserves are >0, energy will be drawn from gonads
  A_1 <- ifelse(A_1 < 0 & gonadCReserves > abs(A_1), 0, 
                ifelse(A_1 < 0 & gonadCReserves < abs(A_1), gonadCReserves - A_1, A_1)
  )
  dG <- ifelse(A_1 < 0 & gonadCReserves > 0, ((gonadCReserves - abs(A_1))>0) * A_1 + ((gonadCReserves - abs(A_1))<0) * gonadCReserves, dG)
  
  # change in gonad, body reserves
  dG <- (dG<0) * 0 + (dG>0) * dG
  dB <- (A_1<0) * A_1 + (A_1>0) * dB
  
  # update state variables that only change during moulting
  dL <- ifelse(is.na(pocConc),0,0)
  dB <- ifelse(is.na(pocConc),0,dB)
  dG <- ifelse(is.na(pocConc),0,dG)
  newStage <- ifelse(is.na(pocConc),stage,stage)
  eggsProduced <- ifelse(is.na(pocConc),0,0)
  newImpSinceSwitch <- ifelse(is.na(pocConc),impSinceSwitch,impSinceSwitch)
  newLastMoult <- ifelse(is.na(pocConc),lastMoult,lastMoult)
  newIMP <- ifelse(is.na(pocConc), 1, 0)
  newGPhase <- ifelse(is.na(pocConc),gPhase,gPhase)
  # ---------------------------------------------------------------------------------------------------------------- #
  # PROCESSES AT MOULTING
  # ---------------------------------------------------------------------------------------------------------------- #
  # During moulting, 
  # determine growth increment based on reserves
  # if(moultDay == time){
  # 
  #   # determine maxiumum growth and shrinkage. Shrinkage has to be multiplied with (moultDay - lastMoult) since its 
  #   # unit is mm d-1
  #   maxGrowth <- 8.63 - 0.14 * bodyLength
  #   maxShrinkage <- 1e-3 - 1.418e-5 * bodyLength * (moultDay - lastMoult)
  #   
  #   # new length is determined by the body reserves at the day of moulting
  #   # 1. When body reserves are lower than the weight of a healthy individual at the given length, shrinkage occurs
  #   # 2. When the body reserves are higher than the weight of a healthy individual at the given length, growth occurs
  #   # In case 2. is true, growth is capped by the maximum growth increment. When the maximum growth is not reached, 
  #   # growth is scaled according to the predicted length at the given dry mass
  #   # In case 1. is true, shrinkage occurs based on the weight loss since the last moult. It is capped by the maximum
  #   # shrinkage. When the dry weight is less than that of a healthy krill with the new length, thin krill will be the outcome
  #   # 1. subtract the carbon costs of moulting from body reserves
  #   dL <- lengthAtDryMass(carbonToDryMass(carbonMass_mg = (bodyCReserves - moultCCosts), tissue = 'bodyReserves'), stage = 2) - bodyLength
  #   
  #   
  #   # Update the state variables:
  #   # 1. rate of change body length
  #   # 2. rate of change body c reserves
  #   # 3. new stage (1 - male, 2 - immature female, 3 - mature female)
  #   # 4. number of produced eggs
  #   # 5. rate of change gonad c reserves
  #   # 6. update impSinceSwitch
  #   # 7. update lastMoult (last timepoint for moulting - alias the current timepoint)
  #   # 8. determine newMoult (new timepoint for moulting)
  #   
  #   dL <- ifelse(dL > maxGrowth, maxGrowth,
  #                ifelse(dL < maxShrinkage, maxShrinkage, dL)
  #   )
  #   # body C-reserves are decreased by the C-costs of moulting
  #   dB <- -1 * moultCCosts
  #   
  #   # change stage- and reproduction specific dynamics
  #   # update reproductive state
  #   oldStage <- stage
  #   newStage <- sexMaturity(inputLength = bodyLength, stage = stage)
  #   
  #   # check reproduction
  #   eggsProduced <- eggsProducedFunc(gPhase = gPhase, inputGonadCReserves = gonadCReserves, 
  #                                    input_dGMax = dG_max, impSinceSwitch = impSinceSwitch)
  #   
  #   # change in gonad reserves is zero, when no spawning occurs. Otherwise, they do not change
  #   dG <- ifelse(eggsProduced > 0, -1 * gonadCReserves, 0)
  #   
  #   # if spawning occurs, impSinceSwitch is set to 1 
  #   newImpSinceSwitch <- ifelse(eggsProduced > 0, 1, impSinceSwitch)
  #   
  #   # determine gonadal phase (when individual just reached sexual maturity, it is set to 1, otherwise the rep-function is called)
  #   oldGPhase <- gPhase
  #   newGPhase <- ifelse(oldStage == 2 & newStage == 3, 1, determineGPhase(gPhase = gPhase, impSinceSwitch = impSinceSwitch, dayLength = photoperiod))
  #   
  #   # when stage switch occurred, impSinceSwitch is set to 1 (the logic is that the next moulting will be the first since the 
  #   # the stage switch)
  #   newImpSinceSwitch <- ifelse(eggsProduced > 0, newImpSinceSwitch,
  #                               ifelse(newGPhase != oldGPhase, 1, impSinceSwitch + 1))
  #   
  #   # determine new imp:
  #   newLastMoult <- moultDay
  #   newMoult <- determineIMP(temperature) + moultDay
  # }
  
  # return the results
  # moultFunc inputs: bodyLength, moultDay, lastMoult, bodyCReserves, moultCCosts, stage, gPhase, gonadCReserves, dG_max, impSinceSwitch, photoperiod, temperature
  changeWhenMoult <- moultFunc(bodyLength = bodyLength, moultDay = moultDay, lastMoult = lastMoult, bodyCReserves = bodyCReserves,
                             moultCCosts = moultCCosts, stage = stage, gPhase = gPhase, gonadCReserves = gonadCReserves,
                             dG_max = dG_max, impSinceSwitch = impSinceSwitch, photoperiod = photoperiod, temperature = temperature,
                             pocConc = pocConc)
  
  #return(c(dL, newStage, dB, dG, newGPhase, newLastMoult, newIMP, newImpSinceSwitch, eggsProduced))
  dL_final <- ifelse(time == moultDay, changeWhenMoult[[1]], dL)
  newStage_final <- ifelse(time == moultDay,changeWhenMoult[[2]], newStage)
  dB_final <- ifelse(time == moultDay,changeWhenMoult[[3]],dB)
  dG_final <- ifelse(time == moultDay,changeWhenMoult[[4]],dG)
  newGPhase_final <- ifelse(time == moultDay,changeWhenMoult[[5]],newGPhase)
  newLastMoult_final <- ifelse(time == moultDay,changeWhenMoult[[6]],lastMoult)
  newIMP_final <- ifelse(time == moultDay, changeWhenMoult[[7]], newIMP)
  newImpSinceSwitch_final <- ifelse(time == moultDay,changeWhenMoult[[8]],newImpSinceSwitch)
  eggsProduced_final <- ifelse(time == moultDay,changeWhenMoult[[9]],eggsProduced)
  
  return(c(dL_final, newStage_final, dB_final, dG_final, newGPhase_final, newLastMoult_final, 
           newIMP_final, newImpSinceSwitch_final, eggsProduced_final))
  #return(modelOutput)
  #return(c(dL, newStage, dB, dG, newGPhase, newLastMoult, newMoult, newImpSinceSwitch, eggsProduced))
}

