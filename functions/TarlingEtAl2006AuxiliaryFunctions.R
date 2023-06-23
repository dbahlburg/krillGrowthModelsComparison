TarlingEtAl2006IMP <- function(bodyLength, stage, temperature, d = 5, roundValues = T){
  
  newIMP <- 
    #juveniles
    (stage == 1) * (d/(exp(-0.392315 + 0.0021159 * bodyLength - 0.404726 * temperature + 0.0687522 * temperature^2)/
                       (1 + exp(-0.392315 + 0.0021159 * bodyLength - 0.404726 * temperature + 0.0687522 * temperature^2)))) +
    
    #male
    (stage == 2) * (d/(exp(1.52581 - 0.0529790 * bodyLength - 0.213042 * temperature + 0.0350464 * temperature^2)/
                         (1 + exp(1.52581 - 0.0529790 * bodyLength - 0.213042 * temperature + 0.0350464 * temperature^2)))) +
    
    #immature female
    (stage == 4) * (d/(exp(-1.55926 + 0.0093231 * bodyLength + 0.375765 * temperature - 0.0733018 * temperature^2)/
                        (1 + exp(-1.55926 + 0.0093231 * bodyLength + 0.375765 * temperature - 0.0733018 * temperature^2)))) +
    
    #female
    (stage %in% c(3,5)) * (d/(exp(2.00098 - 0.0566740 * bodyLength + 0.152815 * temperature - 0.0786357 * temperature^2)/
                         (1 + exp(2.00098 - 0.0566740 * bodyLength + 0.152815 * temperature - 0.0786357 * temperature^2))))
  newIMP <- (roundValues == T) * round(newIMP) + (roundValues == F) * newIMP
  newIMP <- ifelse(newIMP == 0, NA, newIMP)
  #newIMP <- ifelse(newIMP == 0, NA, round(newIMP))
  return(newIMP)
}

TarlingEtAl2006Model <- function(inputLength, inputStage, inputTemperature, inputChla, time, moultDay, oldMoultDay, temperatureHistory, chlorophyllHistory){
  
  growthIncrement <- ifelse(time == moultDay,
                            inputLength * 0.01 * atkinsonEtAlModels(growthType = 'growthIncrement',
                                                     stage = inputStage, #ifelse(inputLength < 35, 'juvenile', 'adultFemale'),
                                                     inputLength = inputLength,
                                                     food = chlorophyllHistory,
                                                     temperature = temperatureHistory),
                            0) 
  
  newMoultDay <- ifelse(is.na(inputChla), 1,
                        ifelse(time == moultDay,
                        TarlingEtAl2006IMP(bodyLength = inputLength, 
                                           stage =  inputStage, #'female'), 
                                           temperature = temperatureHistory),
                        0))
  
  changeOldMoultDay <- ifelse(time == moultDay, newMoultDay, oldMoultDay)
  
  temperatureHistory <- ifelse(time == moultDay, 1/changeOldMoultDay * inputTemperature,
                               ifelse(is.na(inputChla), temperatureHistory,
                               temperatureHistory + 1/changeOldMoultDay * inputTemperature))
  
  chlorophyllHistory <- ifelse(time == moultDay, ifelse(!is.na(inputChla), 1/changeOldMoultDay * inputChla, chlorophyllHistory),
                               ifelse(is.na(inputChla), chlorophyllHistory,
                               chlorophyllHistory + 1/changeOldMoultDay * inputChla))
  
  # under sea ice: growth is zero and moult day is shifted by 1 day
  growthIncrement <- ifelse((growthIncrement + inputLength) > 60, 0, growthIncrement)  
  chlaNAs <- which(is.na(inputChla))
  growthIncrement[chlaNAs] <- 0
  newStage <- ifelse((inputLength + growthIncrement) < 35, 1, 5)
  
  return(c(growthIncrement, newStage, newMoultDay, changeOldMoultDay, temperatureHistory, chlorophyllHistory))
}
