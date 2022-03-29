### DATA FILTERING ###
### The goal of this function is to take in a file of multiple omics data ### 
### and then filter it by a given limit of coefficient ###

# Loading libraries
library(ggplot2)
library(gsubfn)

options(warn = 0) # set=2 to end loops when warnings occur

data_filtering <- function(data) {
  # Gonna write it out first, then make it into a function
  inData = data
  testID  <- inData$PIG_DATE
  OUA     <- inData$OUA
  
  inData = subset(inData, select = -c(PIG_DATE, OUA))
  
  # if the variance is equal to zero, then its a zero-column
  #which(unname(apply(inData, 2, var))==0)
  throwAway <- which(apply(inData, 2, var)==0)[]
  
  # dropping columns by index
  inData = inData[, -throwAway]
  rm(throwAway) # to save memory
  
  ## finding different interesting information on the set
  ## mean, var, STD, ratio of zeros, 
  meanData  <- unname(apply(inData, 2, mean))
  stdData   <- unname(apply(inData, 2, sd))
  varData   <- unname(apply(inData, 2, var))
  ratioData <- unname(apply(inData, 2, function(x){length(which(x==0))/length(x)}))
  CVData    <- unname(apply(inData, 2, function(x){sd(x)/mean(x)}))
                                        
  
  plot(sort(meanData), main = "Mean Data")
  plot(sort(stdData), main = "Standard Variations")
  plot(sort(varData), main = "Variance")
  plot(sort(ratioData), main = "Singleton ratios")
  plot(sort(CVData), main = "Relative standard deviation")
  
  ratData = list(ratioData)
  
  
  # fitting Michaelis-Menten function to data
  time_frame <- seq(from = 1, to = length(ratioData), by = 1)
  Vmax  <- max(ratioData)
  k     <- max(time_frame) / 2
  
  mmModel <- nls(sort(ratioData) ~ Vm*time_frame/(K+time_frame), start = list(Vm=Vmax, K=k))
  
  # parameters estimated including confidence interval
  coef(mmModel)
  confint(mmModel, level = 0.9)
  
  # updating to estimated parameters
  Vmax = coef(mmModel)[[1]]
  K    = coef(mmModel)[[2]]
  
  # visualizing the model in relation to ratio Data
  plot(sort(ratioData) ~ time_frame, col = "grey")
  lines(predict(mmModel) ~ time_frame, lwd = 3, col = "dark red")
  abline(h = 0.51111)
  
  # finding local maxima, potentiel knee points
  kneePoints_collected <- which(diff(sign(diff(ratioData)))==-2)+1
  
  # finding the linear fit from origo to asymptote
  a_linearFit <- max(ratioData) / max(time_frame)
  b_linearFit <- min(ratioData) - a_linearFit*min(time_frame)
  
  # using distance formular on all knee points to line
  # distance formula: d = abs(a*x+b-y)/sqrt(a^2+1)
  longest_distance  <- 0
  longest_x         <- 0
  
  for (i in kneePoints_collected) {
    temp_x = i
    temp_y = Vmax*temp_x / (K*temp_x)
    
    dist = abs(a_linearFit*temp_x + b_linearFit-temp_y)/sqrt(a_linearFit^2+1)
    print(dist)
    if (dist >= longest_distance) {longest_distance = dist; longest_x = temp_x}
  }
  
  # defining the cutoff threshold
  threshold_ratio <- ratioData[longest_x]
  
  # if the ratioData is above the found cutoff, it is to be removed
  throwAway <- which(ratioData <= threshold_ratio)
  
  # dropping columns by index
  inData = inData[, throwAway]
  rm(throwAway) # to save memory
  
  outData = inData
  
  rm(list=setdiff(ls(), c("outData", "testID", "OUA")))
  list(outData, testID, OUA)
}