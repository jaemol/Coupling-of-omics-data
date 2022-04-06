### DBScan density-based spatial clustering ###

# loading libraries
library(dbscan)
library(fpc)
library(factoextra)

# setting seed
set.seed(123)

# excluding testID and OUA (antibiotics used or not)
inData = data
testID  <- inData$testID
OUA     <- inData$OUA

inData = subset(inData, select = -c(testID, OUA))

# selecting arbitrary cluster number, just for visualization
kmeans_model <- kmeans(inData, 3, nstart = 10)

# generating cluster
fviz_cluster(kmeans_model, inData, frame = FALSE, geom = "point")


# making DBSCAN model

findEpsi=function(L, minRange=0,maxRange=3, steps=0.1, maxY=200, minP=10) {
  plot(0,0, col=0,xlim=c(minRange,maxRange), ylim=c(0,maxY))
  legend("topright", legend = c("nClust","nOutliers"), col=1:2, pch=16)
  for(i in seq(minRange,maxRange,steps)) {
    
    DB=dbscan::dbscan(L,i,minP)
    
    points(i,length(unique(DB$cluster)), pch=16, col=1)
    points(i,length(which(DB$cluster==0)), pch=16,col=2)
    
    print(i)
    print(length(unique(DB$cluster)))
    
  }
}

inDataScale=scale(inData, center = T, scale = T)
findEpsi(inDataScale, minRange = 0, maxRange = 65, steps = 1 ,maxY=900)

dbscan::dbscan(inDataScale,6,65)

