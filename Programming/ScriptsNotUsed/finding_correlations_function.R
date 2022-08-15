### Script for finding the correlations   ###
### and extrapolate the important columns ### 
### in a datasheet                        ### 

# loading libraries
# library(dplyr)

# beginning function
extrapolating_correlations <- function(givenData) {
  # creating a new, empty data frame for output data
  outputData <- matrix(0, nrow = length(givenData[1,]), ncol = 1)
  
  # finding length of the dataset
  len <- length(givenData[1,]) 
  
  # looping twice over the data
  for (i in 1:len) {
    for (j in 1:len) {
      
      # if its the same variables, don't continue
      if (i != j & all(givenData[,i] == 0) == FALSE | 
          all(givenData[,j] == 0) == FALSE) {
        
        # finding the correlation coefficient
        corrCoefficient = cor(givenData[,i], givenData[,j])
        print("Correlation coefficient")
        print(corrCoefficient)
        
        
        
        if (abs(corrCoefficient) >= 0.8) {
          # if there is a strong correlation, then include it in a new data frame
          # first, check if its already included
          if (any(names(outputData) == colnames(givenData[i])) == FALSE) {
            outputData$colnames(givenData[i]) <- givenData[,i]
          }
          
          if (any(names(outputData) == colnames(givenData[j])) == FALSE) {
            outputData$colnames(givenData[j]) <- givenData[,j]
          }
          
          print("Dataframe")
          head(outputData, n = 10)
        }
      }
    }
  }
  return(outputData)
}