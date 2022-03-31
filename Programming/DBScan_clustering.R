### DBScan density-based spatial clustering ###

# loading libraries
library(fpc)
library(dbscan)
library(factoextra)

# setting seed
set.seed(123)

# excluding testID and OUA (antibiotics used or not)
inData = data
testID  <- inData$testID
OUA     <- inData$OUA

inData = subset(inData, select = -c(testID, OUA))

# selecting arbitrary cluster number, just for visualization
kmeans_model <- kmeans(inData, 7, nstart = 25)

# generating cluster
fviz_cluster(kmeans_model, inData, frame = FALSE, geom = "point")


# making DBSCAN model


