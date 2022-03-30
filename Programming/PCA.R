### Principal component analysis ### 
# Data needs to transposed..!
# Using 'spectral decomposition', examining covariances / correlation between variables
# Using R-function 'princomp()'

# loading libraries
library(factoextra) # used to create ggplot2-based visualization

# excluding testID and OUA (antibiotics used or not)
inData = data
testID  <- inData$PIG_DATE
OUA     <- inData$OUA

inData = subset(inData, select = -c(testID, OUA))

# transposing data; samples on columns, features on rows
inData = t(inData) 

# making PCA model
pcaModel <- princomp(inData, cor = FALSE, scores = TRUE)
