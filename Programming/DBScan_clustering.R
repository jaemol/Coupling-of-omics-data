### DBScan density-based spatial clustering ###

# loading libraries
library(dbscan)
library(fpc)
library(factoextra)

# loading function
source("Programming/extracting_data_KAT.R")
source("Programming/extracting_data_NATH.R")
source("Programming/data_filtering.R")
source("Programming/Functions.R")


# setting seed
set.seed(123)

# loading data
#chosenWeek <- "Week 04"
#inData <- extracting_data_KAT(whichWeek = chosenWeek)
inData <- extracting_data_NATH()

# filtering data
inData <- data_filtering(inData)

# excluding testID and OUA (antibiotics used or not)
#inData = data
testID  <- inData$testID
OUA     <- inData$OUA

inData = subset(inData, select = -c(testID, OUA))
#inData = t(inData)

#rownames(inData) <- testID


data = outData
phaobacBin <- data$array.phaebac.bin
inData = subset(data, select = -c(array.phaebac.bin))
#rownames(inData) <- commonIDs

# selecting arbitrary cluster number, just for visualization
#kmeans_model <- kmeans(inData, 3, nstart = 10)

# generating cluster
#fviz_cluster(kmeans_model, inData, frame = FALSE, geom = "point")


# making DBSCAN model
inDataScale=data.frame(scale(inData, center = T, scale = T))
findEpsi(t(inDataScale), minRange = 2, maxRange = 10, steps = .1 ,maxY=15,minP = 3)

### results
# week 02: epsilon = 9.1
# Week 04: epsilon = 8.7
###

DB=dbscan::dbscan(t(inDataScale),3.7,3); print(DB)

clusVars=data.frame(Clus=DB$cluster,Vars=colnames(inDataScale))

#sub <- subset(clusVars, Clus==5); print(sub)
subset(clusVars, Clus==1)

#plot(log1p(inData$DATA.Bacteria_Proteobacteria_Gammaproteobacteria_Enterobacterales_Enterobacteriaceae_Escherichia.Shigella),
 #    log1p(inData$Increp64))
#plot(log1p(inData$noquote(sub[1,2])), log1p(inData$noquote(sub[5,2])))

plot(inDataScale[DB$cluster %in% 1])
     