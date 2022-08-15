### Script to implement NetCoMi network analysis ### 
### Using Git-package NetCoMi ### 
### https://github.com/stefpeschel/NetCoMi ###

# loading NetCoMi library
library(NetCoMi)

# loading functions
source("Programming/extracting_data_RWA.R")
source("Programming/data_filtering.R")
source("Programming/Functions.R")

# loading data
chosenWeek      <- "Week 02"
chosenTaxonomy  <- "genus"
inData <- extracting_data_KAT(whichWeek = chosenWeek, whichTaxLevel = chosenTaxonomy, loadOrigData = TRUE)

# filtering data
inData <- data_filtering(inData, whichDataSet = "genom")

# resetting plot window
par(mfrow=c(1,1))

# excluding testID and OUA (antibiotics used or not)
testID  <- inData$testID
OUA     <- inData$OUA
data = subset(inData, select = -c(testID, OUA))

# to use for setting specific colors or shapes to the nodes
dataOneColumns <- grep(x = colnames(data), pattern = "DATA.*")
colVector <- 1:length(colnames(data))
shapeArray <- 1:length(colnames(data))
dataOrigin <- 1:length(colnames(data))
for (i in 1:length(colVector)) {
  if (i<=max(dataOneColumns)) {
    colVector[i] = "green"
    shapeArray[i] = "circle"
    dataOrigin[i] = "16s"
  } else {
    colVector[i] = "red"
    shapeArray[i] = "triangle"
    dataOrigin[i] = "qPCR"
  }
}

# simplifying the metataxonomic names to either genus or species
if (chosenTaxonomy=="genus"){
  colnames(data)=gsub("DATA.Bacteria_[A-Za-z]*[_.][A-Za-z]*[_.][A-Za-z]*[_.][A-Za-z]*[_.](.*)","\\1",colnames(data))
  colnames(data)=gsub("DATA.Archaea_[A-Za-z]*[_.][A-Za-z]*[_.][A-Za-z]*[_.][A-Za-z]*[_.](.*)","\\1",colnames(data))
} else {
  colnames(data)=gsub("DATA.Bacteria_[A-Za-z]*[_.][A-Za-z]*[_.][A-Za-z]*[_.][A-Za-z]*[_.](.*)[_.](.*)","\\1",colnames(data))
  colnames(data)=gsub("DATA.Archaea_[A-Za-z]*[_.][A-Za-z]*[_.][A-Za-z]*[_.][A-Za-z]*[_.](.*)[_.](.*)","\\1",colnames(data))
}


# network untreated, single network with spearman association
data_untreated <- data[which(OUA==1),]
net_single_untreated <- netConstruct(data_untreated,
                                   #filtTax = "highestFreq",
                                   #filtTaxPar = list(highestFreq = 100),
                                   #filtSamp = "totalReads",
                                   #filtSampPar = list(totalReads = 1000),
                                   #filtTax = "highestVar",
                                   #filtTaxPar = list(highestVar = 50),
                                   measure = "spearman",thresh = 0.5,
                                   measurePar = list(nlambda=10, 
                                                     rep.num=10),
                                   normMethod = "none", 
                                   zeroMethod = "none",
                                   sparsMethod = "threshold", 
                                   dissFunc = "signed",
                                   verbose = 3,
                                   seed = 123456)

props_single_untreated <- netAnalyze(net_single_untreated, 
                                   centrLCC = TRUE,
                                   clustMethod = "cluster_fast_greedy",
                                   hubPar = "eigenvector",
                                   weightDeg = FALSE, normDeg = FALSE)

plot(props_single_untreated,
     labelScale = F,
     cexLabels = 1.3,
     title1 = paste("Single network with Spearman, untreated", chosenWeek, chosenTaxonomy),
     showTitle = T,
     nodeSize = "mclr",
     shortenLabels = "none",
     cexNodes = 1.5,
     cexTitle = 2.3,
     nodeColor = "feature",
     featVecCol = dataOrigin,
     colVector = c("green", "red"))

# network treated, single network with spearman association
data_treated <- data[which(OUA==0),]
net_single_treated <- netConstruct(data_treated,
                                     #filtTax = "highestFreq",
                                     #filtTaxPar = list(highestFreq = 100),
                                     #filtSamp = "totalReads",
                                     #filtSampPar = list(totalReads = 1000),
                                     filtTax = "highestVar",
                                     filtTaxPar = list(highestVar = 100),
                                     measure = "spearman",thresh = 0.4,
                                     measurePar = list(nlambda=10, 
                                                       rep.num=10),
                                     normMethod = "none", 
                                     zeroMethod = "none",
                                     sparsMethod = "threshold", 
                                     dissFunc = "signed",
                                     verbose = 3, weighted = T,
                                     seed = 123456)

props_single_treated <- netAnalyze(net_single_treated, 
                                     centrLCC = TRUE,
                                     clustMethod = "cluster_fast_greedy",
                                     hubPar = "eigenvector",
                                     weightDeg = FALSE, normDeg = FALSE)

plot(props_single_treated,
     labelScale = F,
     cexLabels = 1.3,
     title1 = paste("Single network with Spearman, treated", chosenWeek, chosenTaxonomy),
     showTitle = T,
     nodeSize = "mclr",
     shortenLabels = "none",
     cexNodes = 1.5,
     cexTitle = 2.3,
     nodeColor = colVector)
     #colorVec = colVector)

summary(props_single_treated, numbNodes = 5L)


## Comparative network - treated vs untreated ## 
data_untreated <- data[which(OUA==1),]
data_treated <- data[which(OUA==0),]

# Network construction - TDAvsControl
net_untreated_treated <- netConstruct(data = data_untreated, 
                              data2 = data_treated,  
                              #filtTax = "highestVar",
                              #filtTaxPar = list(highestVar = 76),
                              #filtTax = "highestFreq",
                              #filtTaxPar = list(highestFreq = 50),
                              measure = "pear", thresh = 0.2,
                              measurePar = list(nlambda=10, 
                                                rep.num=10),
                              normMethod = "none", 
                              zeroMethod = "none",
                              sparsMethod = "threshold", 
                              dissFunc = "signed",
                              verbose = 3, weighted = T,
                              seed = 123456)

props_untreated_treated <- netAnalyze(net_untreated_treated, 
                              centrLCC = FALSE,
                              avDissIgnoreInf = TRUE,
                              sPathNorm = FALSE,
                              clustMethod = "cluster_fast_greedy",
                              #hubPar = c("degree", "between", "closeness"),
                              hubPar = "eigenvector",
                              hubQuant = 0.9,
                              lnormFit = TRUE,
                              normDeg = FALSE,
                              normBetw = FALSE,
                              normClose = FALSE,
                              normEigen = FALSE)

summary(props_untreated_treated)

plot(props_untreated_treated, 
     #layout = lay_fr,
     #layout = "triangle",
     rmSingles = "inboth",
     sameLayout = TRUE, 
     #nodeColor = colVector,
     #featVecShape = shapeArray,
     nodeSize = "mclr",
     labelScale = F,
     shortenLabels = "none",
     cexNodes = 1.5, 
     cexLabels = 1.3,
     cexHubLabels = 1,
     cexTitle = 3.7,
     groupNames = c("Untreated", "Treated"),
     hubBorderCol  = "gray40")

comp_untreated_treated <- netCompare(props_untreated_treated, permTest = FALSE, verbose = FALSE)

summary(comp_untreated_treated, 
        groupNames = c("Untreated", "Treated"),
        showCentr = c("degree", "between", "closeness"), 
        #showCentr = c("eigenvector"),
        numbNodes = 5)



#####
# using data_analyze
feat1 <- "lachnoclostridium"
feat2 <- "tetO.2"

data_analyze(data = inData, feature1 = feat1, feature2 = feat2)
