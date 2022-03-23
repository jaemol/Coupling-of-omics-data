### DATA FILTERING ###
### The goal of this function is to take in a file of multiple omics data ### 
### and then filter it by a given limit of coefficient ###

# Loading libraries
library(drc)
library(ggplot2)

options(warn = 0) # set=2 to end loops when warnings occur

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
## mean, var, STD, ratio of zeros
meanData  <- unname(apply(inData, 2, mean))
stdData   <- unname(apply(inData, 2, sd))
varData   <- unname(apply(inData, 2, var))
ratioData <- unname(apply(inData, 2, function(x){length(which(x==0))/length(x)}))
CV        <- unname(apply(inData, 2, function(x){sd(x)/mean(x)}))
                                      

plot(sort(meanData), main = "Mean Data")
plot(sort(stdData), main = "Standard Variations")
plot(sort(varData), main = "Variance")
plot(sort(ratioData), main = "Singleton ratios")
plot(sort(CV), main = "Relative standard deviation")


# fitting Michaelis-Menten function to data
Vmax  <- max(ratioData)
k     <- Vmax / 2

mmModel <- nls(data = sort(ratioData), start = list(Vm=Vmax, K=k))
#http://strata.uga.edu/8370/lecturenotes/nonlinearRegression.html
"""
par(mfrow=c(2,1), mar=c(2,2,0,0))
plot(sort(ratioData))
plot(diff(sort(ratioData))~sort(ratioData)[-1])
abline(h=0.85)

par(mfrow=c(1,1), mar=c(2,2,2,2))
"""




outData = inData





