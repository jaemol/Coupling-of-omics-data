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