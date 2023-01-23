library(tidyverse)
# Auxiliary functions needed to run the larval individual based energy budget model of Atkinson et al. (2006)
# Full reference:
# Atkinson, Angus; Shreeve, Rachael S.; Hirst, Andrew G.; Rothery, Peter; Tarling, Geraint A.;
# Pond, David W.; Korb, Rebecca E.; Murphy, Eugene J. ; Watkins, Jonathon L.. 2006 Natural 
# growth rates in Antarctic krill (Euphausia superba): II. Predictive models based on food, 
# temperature, body length, sex, and maturity stage. Limnology and Oceanography, 51 (2). 973-987. 
# -------------------------------------------------------------------------------------------------------------------- #

# The paper provides 5 different parameterizations of equation 4 (listed in Table 5)
# Model 1 and 2 are for estimating growth increment whereas model 1 averages for all krill while 
# model 2 provides maturity-stage-specific parameterizatinos
# 
# Model 3-5 are for estimating daily growth rates whereas 
# model 3 averages for all krill
# model 4 provides maturity-stage-specific parameterizations
# model 5 averages for all krill but is fitted using incubation temperatures (from experiments)
atkinsonEtAlModels <- function(growthType, stage, inputLength, food, temperature){
  
  # DEFINE PARAMETER VALUES
  # equation 1
  dailyGrowthRate <- function(preMoultLength, growthIncrement, interMoultPeriod){
    return(preMoultLength * growthIncrement/(interMoultPeriod * 100))
  }
  
  # equation 2
  preMoultLength <- function(inputLength, growthIncrement){
    return(inputLength * 100/(100 + growthIncrement))
  }
  
  # equation 4
  growth <- function(c_a, c_b, inputLength, c_c, c_d, food, c_e, c_f, temperature, c_g){
    return(c_a + c_b * inputLength + c_c * inputLength^2 + (c_d * food/(c_e + food)) + c_f * temperature + c_g * temperature^2)
  }
  
  if(growthType == 'growthIncrement'){
    
   growthRate <-
      #allKrill
      (stage == 3) * growth(c_a = 6.6, c_b = -0.385, c_c = 0.00259, 
                            c_d = 17.53, c_e = 0.332, c_f = 0.595, 
                            c_g = -0.477, food = food, inputLength = inputLength, temperature = temperature) +
      #juveniles
      (stage == 1) * growth(c_a = 2.52, c_b = -0.184, c_c = 0.00107, 
                            c_d = 17.24, c_e = 0.323, c_f = 0.845, 
                            c_g = -0.548, food = food, inputLength = inputLength, temperature = temperature) +
      #male
      (stage == 2) * growth(c_a = 1.13, c_b = -0.184, c_c = 0.00107, 
                            c_d = 17.24, c_e = 0.323, c_f = 0.845, 
                            c_g = -0.548, food = food, inputLength = inputLength, temperature = temperature) +
      #immature female
      (stage == 4) * growth(c_a = 0.99, c_b = -0.184, c_c = 0.00107, 
                            c_d = 17.24, c_e = 0.323, c_f = 0.845, 
                            c_g = -0.548, food = food, inputLength = inputLength, temperature = temperature) +
      #female
      (stage == 5) * growth(c_a = -0.4, c_b = -0.184, c_c = 0.00107, 
                            c_d = 17.24, c_e = 0.323, c_f = 0.845, 
                            c_g = -0.548, food = food, inputLength = inputLength, temperature = temperature)
   growthRate <- ifelse(growthRate == 0, NA, growthRate)
   return(growthRate)

  }
  
  if(growthType == 'dailyGrowthRate'){
    
   growthRate <- 
      #allKrill
      (stage == 3) * growth(c_a = -0.066, c_b = 0.002, c_c = -0.000061, 
                            c_d = 0.385, c_e = 0.328, c_f = 0.0078, 
                            c_g = -0.0101, food = food, inputLength = inputLength, temperature = temperature) +
      #juveniles
      (stage == 1) * growth(c_a = -0.158, c_b = 0.00674, c_c = -0.000101, 
                            c_d = 0.377, c_e = 0.321, c_f = 0.013, 
                            c_g = -0.0115, food = food, inputLength = inputLength, temperature = temperature) +
      #male
      (stage == 2) * growth(c_a = -0.193, c_b = 0.00674, c_c = -0.000101, 
                            c_d = 0.377, c_e = 0.321, c_f = 0.013, 
                            c_g = -0.0115, food = food, inputLength = inputLength, temperature = temperature) +
      #immature female
      (stage == 4) * growth(c_a = -0.192, c_b = 0.00674, c_c = -0.000101, 
                            c_d = 0.377, c_e = 0.321, c_f = 0.013, 
                            c_g = -0.0115, food = food, inputLength = inputLength, temperature = temperature) +
      #female
      (stage == 5) * growth(c_a = -0.216, c_b = 0.00674, c_c = -0.000101, 
                            c_d = 0.377, c_e = 0.321, c_f = 0.013, 
                            c_g = -0.0115, food = food, inputLength = inputLength, temperature = temperature) +
      #allKrillIncubation
      (stage == 6) * growth(c_a = 0.057, c_b = -0.0019, c_c = -0.000036, 
                                  c_d = 0.345, c_e = 0.297, c_f = -0.011, 
                                  c_g = -0.00569, food = food, inputLength = inputLength, temperature = temperature)

    growthRate <- ifelse(growthRate == 0, NA, growthRate)
    return(growthRate)
    
  }
}

atkinsonFunc <- function(inputFood, inputStage, inputTemperature, inputLength){
  
  growth <- atkinsonEtAlModels(growthType = 'dailyGrowthRate', 
                     stage = inputStage, 
                     inputLength = inputLength,
                     food = inputFood,
                     temperature = inputTemperature)
  
  growth <- ifelse(is.na(inputFood), 0, growth)
  newStage <- ifelse((inputLength + growth) < 35, 1, 3)
  return(c(growth, newStage))
}
