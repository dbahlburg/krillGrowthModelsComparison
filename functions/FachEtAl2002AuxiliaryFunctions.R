# feeding and metabolic activity (Figure 5)
detFeedingActivity <- function(dailyRatio){
  return(ifelse(dailyRatio <=10, 0.1 * dailyRatio, 1))
}

detMetabolicActivity <- function(dayOfYear = NA, dayLength = NA){
  if(!is.na(dayOfYear)){
    return(ifelse(between(dayOfYear, 1, 105), 0,
                  ifelse(between(dayOfYear, 105, 150), -0.5/45 * (dayOfYear-105),
                         ifelse(between(dayOfYear, 150, 250), -0.5,
                                ifelse(between(dayOfYear, 250, 295), -0.5 + 0.5/45 * (dayOfYear-250),
                                       ifelse(between(dayOfYear, 295, 365), 0, NA))))))
  }
}

# allometric functions to convert between body length, dry weight, wet weight and carbon mass
# (Table 2)
massAtLength <- function(inputLength, massType) {
  
  conversionTable <- tibble(length = inputLength) %>% 
    mutate(wetWeight = ifelse(between(length, 1, 5), 0.047 * length ^ 2.121,
                              ifelse(between(length, 5, 10), -6.75+3.30859*length-0.479004*length^2+0.0292664*length^3, 
                                     ifelse(between(length,10,40), 0.0072 * length ^ 3.021, ifelse(length> 40,  0.0016 * length ^ 3.423, NA)))))
  
  wetMass <- as.numeric(conversionTable$wetWeight)
  dryMass <- 0.216 * wetMass
  carbonMass <- 0.366 * dryMass^1.037
  ifelse(massType == 'wetMass', return(wetMass), 
         ifelse(massType == 'dryMass', return(dryMass),
                ifelse(massType == 'carbonMass', return(carbonMass), NA)))
}

lengthAtMass <- function(inputMass, massType = 'carbonMass'){
  loessModelData <- tibble(inputLength = seq(1,90,by = 0.05),
                           carbonMass = massAtLength(inputLength = inputLength, massType = 'carbonMass'),
                           dryMass = massAtLength(inputLength = inputLength, massType = 'dryMass'),
                           wetMass = massAtLength(inputLength = inputLength, massType = 'wetMass'))
  loessModel <- loess(inputLength ~ carbonMass, data = loessModelData, span = 0.05)
  ifelse(massType == 'carbonMass', 
         return(as.numeric(predict(loessModel, tibble(carbonMass = inputMass)))),NA)
  
}

# Filtration rates (Table 3)
metabolicRates <- function(inputDryMass, rateType){
  
  #define filtration Rates
  if (rateType %in% c('filtrationRateWaterColumn','filtrationRateSeaIce')){
    
    modelDataFiltration <-
      tibble(
        modDryMass = c(seq(0.01, 26, length.out = 30), seq(84, 100, length.out = 30)),
        filtrationRate = ifelse(
          between(modDryMass, 0, 26),
          0.00085 * modDryMass^0.825,
          0.00343 * modDryMass^0.514))
    filtrationLoess <- loess(filtrationRate ~ modDryMass, data = modelDataFiltration)
    
    conversionTable <- tibble(dryMass = inputDryMass) %>% 
      mutate(filtrationRateWaterColumn = ifelse(dryMass<26,  0.00085 * inputDryMass^0.825,
                                                ifelse(between(dryMass, 26, 84), as.numeric(predict(filtrationLoess, data.frame(modDryMass = dryMass))), 
                                                       ifelse(dryMass > 84, 0.00343 * inputDryMass^0.514, NA))))
    
    filtrationRateWaterColumn <- as.numeric(conversionTable$filtrationRateWaterColumn)
    filtrationRateSeaIce <- 0.05 * filtrationRateWaterColumn
    
    ifelse(rateType == 'filtrationRateWaterColumn', return(filtrationRateWaterColumn), 
           ifelse(rateType == 'filtrationRateSeaIce', return(filtrationRateSeaIce), NA))
  }
  
  # define respiration rates
  if(rateType %in% c('carbonRespirationMlD','oxygenRespirationUlH')){
    
    modelDataRespiration <- tibble(modDryMass = c(seq(0.01, 1.5, length.out = 30), seq(6.3, 7.8, length.out = 30)),
                                   respirationRate = ifelse(
                                     between(modDryMass, 0, 1.5),
                                     0.686*modDryMass^1.031,
                                     0.847 * modDryMass^0.85))
    loessRespiration <- loess(respirationRate ~ modDryMass, data = modelDataRespiration)
    
    conversionTable <- tibble(dryMass = inputDryMass) %>% 
      mutate(respirationRate = ifelse(dryMass<1.5,  0.686*dryMass^1.031,
                                      ifelse(between(dryMass, 1.5, 6.3), as.numeric(predict(loessRespiration, data.frame(modDryMass = dryMass))), 
                                             ifelse(dryMass > 6.3, 0.847 * dryMass^0.85, NA))))
    
    oxygenRate <- as.numeric(conversionTable$respirationRate)
    carbonRate <- oxygenRate * 1e-3 * 24 * 0.5357
    ifelse(rateType == 'oxygenRespirationUlH', return(oxygenRate),
           ifelse(rateType == 'carbonRespirationMlD', return(carbonRate), NA))
  }
}

# gamma (equation 4) which stands for the time spend on feeding
feedingTime <- function(inputLength){
  return(ifelse(inputLength < 12, 0.9, 0.75))
}

# Define a function that returns the dayOfYear at a given simulation timepoint
dayOfYearFunc <- function(t, startDay){
  return(ceiling((t+startDay) - ceiling(((t+startDay)-365)/365) * 365))
}

# growth model
fachEtAl2002Model <- function(y, time, chlA, temperature, heterotrophicCarbon, iceAlgae, parameters = c(chlorophyllToCarbon = 50, assimilationRate = 0.8)){
  
  # ------------------------------------------------------------------------------------------------------------------------------------------ #
  # calculate length and size class at current time step
  currentLength <- lengthAtMass(inputMass = y, massType = 'carbonMass')
  currentSizeClass <- match(round(currentLength*4)/4, sizeClassesVec)
  dayOfYear <- dayOfYearFunc(time, startDay = startDay)
  
  # calculate amount of assimilated carbon (heterotrophic food included for individuals>18mm)
  # temperature scaling of ingestion rates using a q10-function
  ingestionRate <- lookUpTableSizeClassTraits$filtrationRateWaterColumn[currentSizeClass] * 3.5^((temperature - 0)/10)
  ingestionRateSeaIce <- lookUpTableSizeClassTraits$filtrationRateSeaIce[currentSizeClass] * 3.5^((temperature - 0)/10)
  ingestedCarbonWater <- lookUpTableSizeClassTraits$feedingTime[currentSizeClass] * ingestionRate * parameters["chlorophyllToCarbon"] * chlA + 
    (currentLength > 18) * (heterotrophicCarbon * lookUpTableSizeClassTraits$feedingTime[currentSizeClass] * ingestionRate)
  ingestedCarboniceAlgae <- 0.05 * lookUpTableSizeClassTraits$feedingTime[currentSizeClass] * ingestionRateSeaIce * parameters["chlorophyllToCarbon"] * iceAlgae
  ingestedCarbon <- ingestedCarbonWater + ingestedCarboniceAlgae
  assimilatedCarbon <- parameters["assimilationRate"] * ingestedCarbon
  
  # calculate metabolic costs
  metabolicActivity <- lookUpTableTimeVariables$feedingActivity[dayOfYear]
  dailyRation <- assimilatedCarbon/lookUpTableSizeClassTraits$carbonWeight[currentSizeClass] * 100
  dailyRation <- ifelse(dailyRation < 1.5, 0, dailyRation)
  
  # temperature scaling of carbon respiration using a q10-function
  carbonRespiration <- lookUpTableSizeClassTraits$carbonRespiration[currentSizeClass] *
    (1 + metabolicActivity + detFeedingActivity(dailyRatio = dailyRation)) * 3.5^((temperature - 0)/10)
  
  # net production equals the difference between assimilated and metabolized carbon.
  # the net balance is added (or subtracted) to carbon weight. Growth increments are capped at 0.25mm per day
  # and 20% of daily ration to avoid unrealistically high growth rates
  netProduction <- (ifelse(dailyRation > 20, 0.2 * y, assimilatedCarbon) - carbonRespiration)
  netProduction <- ifelse((lengthAtMass(inputMass = (y + netProduction), massType = 'carbonMass') - currentLength) > 0.25,
                          massAtLength(inputLength = (currentLength+0.25), massType = 'carbonMass') -  massAtLength(inputLength = currentLength, massType = 'carbonMass'), netProduction)
  
  netProduction <- ifelse(is.na(chlA), 0, netProduction)
  return(netProduction) 
}

#Runge-Kutta-Solver
rk4stepFach <- function(state, time, inputChla, inputTemp, inputHetC, inputIceAlgae){
  
  dt <- 1
  k1 <- fachEtAl2002Model(y = state, time = time, chlA = inputChla, temperature = inputTemp, heterotrophicCarbon = inputHetC, iceAlgae = inputIceAlgae)
  k2 <- fachEtAl2002Model(y = state+0.5*k1*dt, time = time + 0.5 * dt, chlA = inputChla, temperature = inputTemp, heterotrophicCarbon = inputHetC, iceAlgae = inputIceAlgae)
  k3 <- fachEtAl2002Model(y = state+0.5*k2*dt, time = time + 0.5 * dt, chlA = inputChla, temperature = inputTemp, heterotrophicCarbon = inputHetC, iceAlgae = inputIceAlgae)
  k4 <- fachEtAl2002Model(y = state+k3*dt, time = time + dt, chlA = inputChla, temperature = inputTemp, heterotrophicCarbon = inputHetC, iceAlgae = inputIceAlgae)
  
  return(dt * (k1 + 2*k2 + 2*k3 + k4) /6)
}

