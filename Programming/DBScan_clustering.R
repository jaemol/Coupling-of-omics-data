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
#inData = t(inData)

colnames(inData) <- testID

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
    
    #print(i)
    #sprintf("clusters: %f",length(unique(DB$cluster)))
    
  }
  abline(v=seq(minRange,maxRange,steps))
  
}

inDataScale=data.frame(scale(inData, center = T, scale = T))
findEpsi(t(inDataScale), minRange = 30, maxRange = 32, steps = .1 ,maxY=100,minP = 3)

DB=dbscan::dbscan(t(inDataScale),31.3,3)

clusVars=data.frame(Clus=DB$cluster,Vars=colnames(inDataScale))

subset(clusVars, Clus==0)


plot(log1p(inData$DATA.Bacteria_Proteobacteria_Gammaproteobacteria_Enterobacterales_Enterobacteriaceae_Escherichia.Shigella),
     log1p(inData$Increp64))
