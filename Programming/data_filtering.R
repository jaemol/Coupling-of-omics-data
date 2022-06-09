### DATA FILTERING ###
### The goal of this function is to take in a file of multiple omics data ### 
### and then filter it by a given limit of coefficient ###

# Loading libraries
#library(ggplot2)
library(gsubfn)

options(warn = 0) # set=2 to end loops when warnings occur

data_filtering <- function(data) {
  # Gonna write it out first, then make it into a function
  inData = data
  testID  <- inData$testID
  OUA     <- inData$OUA
  
  inData = subset(inData, select = -c(testID, OUA))
  
  # if the variance is equal to zero, then its a zero-column
  #which(unname(apply(inData, 2, var))==0)
  throwAway <- which(apply(inData, 2, var)==0)[]
  
  # dropping columns by index
  inData = inData[, -throwAway]
  rm(throwAway) # to save memory
  
  print("Generating metadata...")
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
  
  
  # fitting Michaelis-Menten function to data
  print("Fitting Michaelis-Menten function...")
  time_frame <- seq(from = 1, to = length(ratioData), by = 1)
  
  #mmModel <- nls(sort(ratioData) ~ Vm*time_frame/(K+time_frame), start = list(Vm=max(ratioData), K=max(time_frame) / 2))
  mmModel <- nls(sort(ratioData) ~ Vm/(1+exp(K+time_frame)), start = list(Vm=max(ratioData), K=max(time_frame) / 2))
  
  # parameters estimated including confidence interval
  coef(mmModel)
  confint(mmModel, level = 0.9)
  
  # defining to estimated parameters
  Vmax = coef(mmModel)[[1]]
  K    = coef(mmModel)[[2]]
  
  # visualizing the model in relation to ratio Data
  plot(sort(ratioData) ~ time_frame, col = "grey")
  lines(predict(mmModel) ~ time_frame, lwd = 3, col = "dark red")
  #abline(h = threshold_ratio)
  
  # finding local maxima, potentiel knee points
  kneePoints_collected <- which(diff(sign(diff(ratioData)))==-2)+1
  
  # finding the linear fit from origo to asymptote
  a_linearFit <- max(ratioData) / max(time_frame)
  b_linearFit <- min(ratioData) - a_linearFit*min(time_frame)
  
  max_x = which(sort(ratioData)>=max(ratioData)*0.97)[1]
  max_y = Vmax*max_x / (K+max_x)
  a_linearFit <- (max_y - min(ratioData)) / (max_x-min(time_frame))
  b_linearFit <- min(ratioData) - a_linearFit*min(time_frame)
  
  abline(a = b_linearFit, b=a_linearFit)
  
  # using distance formular on all knee points to line
  # distance formula: d = abs(a*x+b-y)/sqrt(a^2+1)
  longest_distance  <- 0
  longest_x         <- 0
  
  #i=kneePoints_collected[400]
  for (i in which(kneePoints_collected<=max_x)) {
    temp_x = i
    temp_y = Vmax*temp_x / (K+temp_x)
    
    #points(temp_x, temp_y)
    
    dist = abs(a_linearFit*temp_x + b_linearFit-temp_y)/sqrt(a_linearFit^2+1)
    
    if (dist >= longest_distance) {longest_distance = dist; longest_x = temp_x}
  }
  
  # defining the cutoff threshold
  threshold_ratio <- Vmax*longest_x / (K+longest_x)
  
  points(longest_x, threshold_ratio)
  abline(h = threshold_ratio)
  
  # if the ratioData is above the found cutoff, it is to be removed
  keepIn <- which(ratioData <= threshold_ratio)
  
  # dropping columns by index
  sprintf("Filtering out data with zero-ratio above %f...", threshold_ratio)
  inData = inData[, keepIn]
  rm(keepIn) # to save memory
  
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
  
  
  inData = cbind(testID, OUA, inData)
  
  outData = inData
  
  rm(list=setdiff(ls(), "outData"))
  print("Filtering done...")
  return(outData)
}