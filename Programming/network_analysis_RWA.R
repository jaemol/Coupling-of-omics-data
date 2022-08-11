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
#inData = data
testID  <- inData$testID
OUA     <- inData$OUA
data = subset(inData, select = -c(testID, OUA))

dataOneColumns <- grep(x = colnames(data), pattern = "DATA.*")
colVector <- 1:length(colnames(data))
for (i in 1:length(colVector)) {
  if (i<=max(dataOneColumns)) {colVector[i] = "green"} else {colVector[i] = "red"}
}


if (chosenTaxonomy=="genus"){
  colnames(data)=gsub("DATA.Bacteria_[A-Za-z]*[_.][A-Za-z]*[_.][A-Za-z]*[_.][A-Za-z]*[_.](.*)","\\1",colnames(data))
  colnames(data)=gsub("DATA.Archaea_[A-Za-z]*[_.][A-Za-z]*[_.][A-Za-z]*[_.][A-Za-z]*[_.](.*)","\\1",colnames(data))
} else {
  colnames(data)=gsub("DATA.Bacteria_[A-Za-z]*[_.][A-Za-z]*[_.][A-Za-z]*[_.][A-Za-z]*[_.](.*)[_.](.*)","\\1",colnames(data))
  colnames(data)=gsub("DATA.Archaea_[A-Za-z]*[_.][A-Za-z]*[_.][A-Za-z]*[_.][A-Za-z]*[_.](.*)[_.](.*)","\\1",colnames(data))
}


# building single network with SPRING as association measure - Full dataset, both treated and untreated
net_single_fullSet <- netConstruct((data),
                           filtTax = "highestVar",
                           filtTaxPar = list(highestVar = 150),
                           #filtSamp = "totalReads",
                           #filtSampPar = list(totalReads = 1000),
                           measure = "spearman",thresh = 0.5,
                           measurePar = list(nlambda=10, 
                                             rep.num=10),
                           normMethod = "none", 
                           zeroMethod = "none",
                           sparsMethod = "threshold", 
                           dissFunc = "signed",
                           verbose = 3,weighted = T,
                           seed = 123456)

props_single_fullSet <- netAnalyze(net_single_fullSet, 
                           centrLCC = TRUE,
                           clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector",
                           weightDeg = FALSE, normDeg = FALSE)

#?summary.microNetProps
summary(props_single_fullSet, numbNodes = 5L)


plot(props_single_fullSet,
     labelScale = F,
     cexLabels = 1.3,
     #title1 = paste("Single network with Spearman",chosenWeek, chosenTaxonomy),
     title1 = "Single network with Spearman",
     showTitle = T,
     cexTitle = 2.3)
     #nodeColor = colVector)
     
# adding legend
#legend("topright", cex = 0.5, title = "estimated association:",
 #      legend = c("Metataxonomic","Genomic"), lty = 1, lwd = 1, col = c("green","red"), 
  #     bty = "n", horiz = TRUE)

legend(x=0.85,y=0.9,legend=c("Positive","Negative"),
       cex=0.6,col=c("green","red"),pch=c(".","."),lwd = c(3,3))
legend(x=0.85,y=0.6,legend=c("Metataxonomic","Genomic"),
       cex=0.6,col=c("green","red"),pch=c(16,16),lwd = c(3,3))



###

cor(log1p(data$Jannaschia[-13]),log1p(data$`943.99194`[-13]),method = "pear" )^2
plot(log1p(data$Jannaschia[-13]),log1p(data$`943.99194`[-13]))
LM=lm(log1p(data$Jannaschia[-13])~log1p(data$`943.99194`[-13]))
summary(LM)
plot(LM)
###

#cor(inData$TolC1, inData$Escherichia.Shigella, method = "spear")


# network untreated, single network with spearman association
data_untreated <- data[which(OUA==0),]
net_single_untreated <- netConstruct(data_untreated,
                                   #filtTax = "highestFreq",
                                   #filtTaxPar = list(highestFreq = 100),
                                   #filtSamp = "totalReads",
                                   #filtSampPar = list(totalReads = 1000),
                                   filtTax = "highestVar",
                                   filtTaxPar = list(highestVar = 50),
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
     cexTitle = 2.3)
#nodeColor = colVector)


legend(x=0.85,y=0.9,legend=c("Metataxonomic","Genomic"),
       cex=0.6,col=c("green","red"),pch=c(16,16),lwd = c(3,3))


# network treated, single network with spearman association
data_treated <- data[which(OUA==1),]
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
     cexTitle = 2.3)
     #nodeColor = colVector)



legend(x=0.85,y=0.9,legend=c("Metataxonomic","Genomic"),
       cex=0.6,col=c("green","red"),pch=c(16,16),lwd = c(3,3))

summary(props_single_treated, numbNodes = 5L)

## Comparative network - treated vs untreated
data_untreated <- data[which(OUA==0),]
data_treated <- data[which(OUA==1),]

# Network construction - TDAvsControl
net_untreated_treated <- netConstruct(data = data_untreated, 
                              data2 = data_treated,  
                              filtTax = "highestVar",
                              filtTaxPar = list(highestVar = 76),
                              #filtTax = "highestFreq",
                              #filtTaxPar = list(highestFreq = 50),
                              measure = "pear", thresh = 0.7,
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

#summary(props_untreated_treated)

# for saving the plot as an image
# png(filename = "Treated_Vs_Untreated_thres7.png",
#     width = 4000, height = 3000,units = "px", pointsize = 12,
#     bg = "white", res = 300, family = "", restoreConsole = TRUE,
#     type = c("windows", "cairo", "cairo-png"),
#     symbolfamily = "default")

plot(props_untreated_treated, 
     #layout = lay_fr,
     #layout = "circle",
     sameLayout = TRUE, 
     nodeColor = "cluster",
     nodeSize = "mclr",
     labelScale = F,
     shortenLabels = "none",
     cexNodes = 1.5, 
     cexLabels = 1.3,
     cexHubLabels = 1,
     cexTitle = 3.7,
     groupNames = c("Untreated", "Treated"),
     hubBorderCol  = "gray40")
#dev.off() # shutting off image saving

# legend("bottomleft", title = "estimated association:", legend = c("+","-"), 
#        col = c("#009900","red"), inset = 0.02, cex = 2, lty = 1, lwd = 4, 
#        bty = "n", horiz = TRUE)


comp_untreated_treated <- netCompare(props_untreated_treated, permTest = FALSE, verbose = FALSE)

summary(comp_untreated_treated, 
        groupNames = c("Untreated", "Treated"),
        showCentr = c("degree", "between", "closeness"), 
        #showCentr = c("eigenvector"),
        numbNodes = 5)



#####
# trying data_analyze

#feat1 <- "Cohaesibacter"
#feat2 <- "611.1929419"
feat1 <- "lachnoclostridium"
feat2 <- "tetO.2"

data_analyze(data = inData, feature1 = feat1, feature2 = feat2)
