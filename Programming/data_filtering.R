### DATA FILTERING ###
### The goal of this function is to take in a file of multiple omics data ### 
### and then filter it by a given limit of coefficient ###


# Gonna write it out first, then make it into a function

dataLength = len_16s + len_qpcr
inData = data

# coefficient for checking the worth of the attributes: 
## variance / mean
sumCoef = 0
# for loop running through the given data set, finding the mean coefficient values for particular set
for (i in 3:length(inData[,1])){
  tempCoef = (sd(data[,i]) / mean(data[,i]))
  sumCoef = sumCoef + tempCoef
}
  
sumCoef = sumCoef / dataLength

# now we have the baseline, to filter the data set
# now we run through the data set again, and drop uninteresting attributes
for (j in 3:length(inData[,1])) {
  currCoef = (sd(data[,i]) / mean())
}
