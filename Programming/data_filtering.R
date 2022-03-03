### DATA FILTERING ###
### The goal of this function is to take in a file of multiple omics data ### 
### and then filter it by a given limit of coefficient ###

options(warn = 0) # set=2 to end loops when warnings occur

# Gonna write it out first, then make it into a function

dataLength = len_16s + len_qpcr
inData = data
testID  <- inData$PIG_DATE
OUA     <- inData$OUA

inData = subset(inData, select = -c(PIG_DATE, OUA))

## We start by running through and filtering out zero-columns
throwAway <- c()
#throwAway=seq(-1, NCOL(inData ))
for (i in 1:NCOL(inData)) {
  # if the variance is equal to zero, then its a zero-column
  if (var(inData[,i]) == 0) {
    throwAway = c(throwAway, i)
    #throwAway[i]=i
  }
}

# dropping columns by index
inData = inData[, -throwAway]


# checking how many of the respective omics attributes were dropped
#print("Number of metataxonomic attributes dropped") 
#length(throwAway[throwAway <= len_16s]) # 12196

#print("Number of genomic attributes dropped")
#length(throwAway[throwAway > len_16s])  # 31


## finding different interesting information on the set
## mean, var, STD, ratio of zeros
meanData  <- unname(apply(inData, 2, mean))
stdData   <- unname(apply(inData, 2, sd))
varData   <- unname(apply(inData, 2, var))
ratioData <- unname(apply(inData, 2, function(x){length(which(x==0))/length(x)}))
                                      

plot(meanData, main = "Mean Data")
plot(stdData, main = "Standard Variations")
plot(varData, main = "Variance")
plot(ratioData, main = "Singleton ratios")


# filtering out all of the singleton-zero ratios
singletonValue = max(ratioData)

throwAway <- c()

for (i in 1:NCOL(inData)) {
  # if the value of ratioData is equal to singletonValue, throw away
  if (ratioData[i] == singletonValue) {
    throwAway = c(throwAway, i)
    #throwAway[i]=i
  }
}

# dropping columns by index
inData = inData[, -throwAway]

# singleton-zero ratio filtering complete





