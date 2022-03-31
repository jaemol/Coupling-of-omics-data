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

findEpsi=function(L, minRange=0,maxRange=3, steps=0.1, maxY=200, minP=10) {
  plot(0,0, col=0,xlim=c(minRange,maxRange), ylim=c(0,maxY))
  legend("topright", legend = c("nClust","nOutliers"), col=1:2, pch=16)
  for(i in seq(minRange,maxRange,steps)) {
    
    DB=dbscan::dbscan(L,i,minP)
    
    points(i,length(unique(DB$cluster)), pch=16, col=1)
    points(i,length(which(DB$cluster==0)), pch=16,col=2)
    
  }
}

inDataScale=scale(inData, center = T, scale = T)
findEpsi(inDataScale, minRange = 0, maxRange = 100, steps = 10,maxY=30)

dbscan::dbscan(scale(inData),10000,10)

head(colnames(inData))

origNames=colnames(inData)
newNames=apply(stringr::str_split_fixed(string = origNames, pattern = "_",8)[,1:6],1, paste, collapse="_")

length(unique(newNames))
uniqNames=unique(newNames)

newDat=data.frame(dummy=1:NROW(inData))

j=uniqNames[1]
for(j in uniqNames) {
  
  jIndx=grep(j,origNames )
  if(length(jIndx)>1) {
    newDat=cbind(newDat,rowSums(inData[,jIndx]))
  } else {
    newDat=cbind(newDat,(inData[,jIndx]))
  }
}

newDat=newDat[,-1]
colnames(newDat)=uniqNames
