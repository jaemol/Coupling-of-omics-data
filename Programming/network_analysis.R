### Script to implement NetCoMi network analysis ### 
### Using Git-package NetCoMi ### 
### https://github.com/stefpeschel/NetCoMi ###

## installing package first time 

if(!requireNamespace("BiocManager", quietly = TRUE)){
  utils::install.packages("BiocManager")
}

BiocManager::install(pkgs = c("Biobase", "doSNOW", "fdrtool", "filematrix",
                              "foreach", "graphics", "grDevices", "gtools",
                              "huge", "igraph", "MASS", "Matrix", "phyloseq",
                              "pulsar", "qgraph", "Rdpack", "snow", "SPRING",
                              "stats", "utils", "vegan", "WGCNA"))

BiocManager::install("GO.db")


devtools::install_github("GraceYoon/SPRING")
devtools::install_github("zdk123/SpiecEasi")

devtools::install_github("stefpeschel/NetCoMi", 
                         dependencies = c("Depends", "Imports"),
                         repos = c("https://cloud.r-project.org/",
                                   BiocManager::repositories()))


# loading NetCoMi library
library(NetCoMi)

# loading functions
source("Programming/extracting_data_KAT.R")
source("Programming/data_filtering.R")

# loading data
chosenWeek <- "Week 04"
inData <- extracting_data_KAT(whichWeek = chosenWeek)

# filtering data
inData <- data_filtering(inData)

# excluding testID and OUA (antibiotics used or not)
#inData = data
testID  <- inData$testID
OUA     <- inData$OUA
data = subset(inData, select = -c(testID, OUA))

dataOneColumns <- grep(x = colnames(data), pattern = "DATA.*")


colVector <- 1:length(colnames(data))
for (i in 1:length(colVector)) {
  if (i<=max(dataOneColumns)) {colVector[i] = "green"} else {colVector[i] = "red"}
}



colnames(data)=gsub("DATA.Bacteria_[A-Za-z]*[_.][A-Za-z]*[_.][A-Za-z]*[_.][A-Za-z]*[_.](.*)","\\1",colnames(data))
colnames(data)=gsub("DATA.Archaea_[A-Za-z]*[_.][A-Za-z]*[_.][A-Za-z]*[_.][A-Za-z]*[_.](.*)","\\1",colnames(data))



# building single network with SPRING as association measure - Full dataset, both treated and untreated
net_single_fullSet <- netConstruct(data,
                           #filtTax = "highestFreq",
                           #filtTaxPar = list(highestFreq = 100),
                           #filtSamp = "totalReads",
                           #filtSampPar = list(totalReads = 1000),
                           measure = "spearman",thresh = 0.6,
                           measurePar = list(nlambda=10, 
                                             rep.num=10),
                           normMethod = "none", 
                           zeroMethod = "none",
                           sparsMethod = "threshold", 
                           dissFunc = "signed",
                           verbose = 3,
                           seed = 123456)

props_single_fullSet <- netAnalyze(net_single_fullSet, 
                           centrLCC = TRUE,
                           clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector",
                           weightDeg = FALSE, normDeg = FALSE)

plot(props_single_fullSet,
     labelScale = F,
     cexLabels = 1.3,
     title1 = paste("Single network with Spearman",chosenWeek),
     showTitle = T,
     cexTitle = 2.3,
     nodeColor = colVector)
     

#cor(inData$TolC1, inData$Escherichia.Shigella, method = "spear")


# network untreated, single network with spearman association
data_untreated <- data[which(OUA==0),]
net_single_untreated <- netConstruct(data_untreated,
                                   #filtTax = "highestFreq",
                                   #filtTaxPar = list(highestFreq = 100),
                                   #filtSamp = "totalReads",
                                   #filtSampPar = list(totalReads = 1000),
                                   measure = "spearman",thresh = 0.6,
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
     title1 = "Single network with Spearman, untreated",
     showTitle = T,
     cexTitle = 2.3,
     nodeColor = colVector)

# network treated, single network with spearman association
data_treated <- data[which(OUA==1),]
net_single_treated <- netConstruct(data_treated,
                                     #filtTax = "highestFreq",
                                     #filtTaxPar = list(highestFreq = 100),
                                     #filtSamp = "totalReads",
                                     #filtSampPar = list(totalReads = 1000),
                                     measure = "spearman",thresh = 0.5,
                                     measurePar = list(nlambda=10, 
                                                       rep.num=10),
                                     normMethod = "none", 
                                     zeroMethod = "none",
                                     sparsMethod = "threshold", 
                                     dissFunc = "signed",
                                     verbose = 3,
                                     seed = 123456)

props_single_treated <- netAnalyze(net_single_treated, 
                                     centrLCC = TRUE,
                                     clustMethod = "cluster_fast_greedy",
                                     hubPar = "eigenvector",
                                     weightDeg = FALSE, normDeg = FALSE)

plot(props_single_treated,
     labelScale = F,
     cexLabels = 1.3,
     title1 = "Single network with Spearman, treated",
     showTitle = T,
     cexTitle = 2.3,
     nodeColor = colVector)
