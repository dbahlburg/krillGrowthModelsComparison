WiedenmannEtAl2008Model <- function(inputLength, inputTemperature, inputChla, time, moultDay, oldMoultDay, temperatureHistory, chlorophyllHistory){
  
  # calculate growth increment based on Atkinson et al. (2006) model
  # environmental conditions are averaged over the last IMP
  growthIncrement <- ifelse(time == moultDay,
                            inputLength * 0.01 * atkinsonEtAlModels(growthType = 'growthIncrement',
                                                                    stage = 3, 
                                                                    inputLength = inputLength,
                                                                    food = chlorophyllHistory,
                                                                    temperature = temperatureHistory),
                            0) 
  
  
  # when current timepoint corresponds to day of moulting, new IMP is determined
  # based on the Kawaguchi et al. (2006) IMP-model, otherwise it is zero
  newMoultDay <- ifelse(time == moultDay,
                        round(exp(3.5371 - 0.5358 * log(inputTemperature + 2))),
                        0)
  
  changeOldMoultDay <- ifelse(time == moultDay, newMoultDay, oldMoultDay)
  
  # Update temperature history which contains the averaged temperatures experienced by the krill
  # during the current IMP
  temperatureHistory <- ifelse(time == moultDay, 1/changeOldMoultDay * inputTemperature,
                               ifelse(is.na(inputChla), temperatureHistory,
                                      temperatureHistory + 1/changeOldMoultDay * inputTemperature))
  
  # Update chlorophyll a history which contains the averaged chlorophyll a concentrations
  # experienced by the krill during the current IMP
  chlorophyllHistory <- ifelse(time == moultDay, ifelse(!is.na(inputChla), 1/changeOldMoultDay * inputChla, chlorophyllHistory),
                               ifelse(is.na(inputChla), chlorophyllHistory,
                                      chlorophyllHistory + 1/changeOldMoultDay * inputChla))
  
  # Calculate growth increment, capped at 60mm
  growthIncrement <- ifelse((growthIncrement + inputLength) > 60, 0, growthIncrement)
  
  # Under sea ice growth is zero and the imp is extended by 1 day
  chlaNAs <- which(is.na(inputChla))
  newMoultDay[chlaNAs] <- 1
  growthIncrement[chlaNAs] <- 0
  return(c(growthIncrement, newMoultDay, changeOldMoultDay, temperatureHistory, chlorophyllHistory))
}
