# Define a function that returns the dayOfYear at a given simulation timepoint
dayOfYearFunc <- function(t, startDay){
  return(ceiling((t+startDay) - ceiling(((t+startDay)-365)/365) * 365))
}
