### DATA FILTERING ###
### The goal of this function is to take in a file of multiple omics data ### 
### and then filter it by a given limit of coefficient ###


# Gonna write it out first, then make it into a function

dataLength = len_16s + len_qpcr
inData = data


## We start by running through and filtering out zero-columns
throwAway <- c()
# disregarding the first two columns - ID and OUA
for (i in 3:length(inData[1,])) {
  # if the variance is equal to zero, then its a zero-column
  if (var(inData[,i]) == 0 | is.na(var(inData[,i]))) {
    throwAway = rbind(throwAway, i)
  }
}

# dropping columns by index
inData = inData[, -throwAway] # losing 12230 attributes

# checking how many of the respective omics attributes were dropped
print("Number of metataxonomic attributes dropped") 
length(throwAway[throwAway <= len_16s]) # 12199

print("Number of genomic attributes dropped")
length(throwAway[throwAway > len_16s])  # 31

## next, we filtrate based on the worth of the attributes, by finding a baseline
# coefficient for checking the worth of the attributes: 
## standard variation / mean
listCoef <- c()
# for loop running through the given data set, finding the mean coefficient values for particular set
for (i in 3:length(inData[1,])){
  tempCoef = (sd(inData[,i]) / mean(inData[,i]))
  listCoef = rbind(listCoef, tempCoef)
}
  
# this is the baseline, that we can now use to filter through data
baseCoef = mean(listCoef)


# filtering through data, keeping all attributes, that are above or equal to baseline coefficient
throwAway <- c()
for (i in 3: length(inData[1,])) {
  tempCoef = (sd(inData[,i]) / mean(inData[,i]))
  if (tempCoef < baseCoef) {
    throwAway = rbind(throwAway, i) 
  }
}

# excluding the attributes with too low interest
inData = inData[, -throwAway] # losing 631, now have 716 attributes!


