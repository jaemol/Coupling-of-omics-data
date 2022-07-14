### Script to implement NetCoMi network analysis ### 
### Using Git-package NetCoMi ### 
### https://github.com/stefpeschel/NetCoMi ###

## installing package first time 

# if(!requireNamespace("BiocManager", quietly = TRUE)){
#   utils::install.packages("BiocManager")
# }
# 
# BiocManager::install(pkgs = c("Biobase", "doSNOW", "fdrtool", "filematrix",
#                               "foreach", "graphics", "grDevices", "gtools",
#                               "huge", "igraph", "MASS", "Matrix", "phyloseq",
#                               "pulsar", "qgraph", "Rdpack", "snow", "SPRING",
#                               "stats", "utils", "vegan", "WGCNA"))
# 
# BiocManager::install("GO.db")
# 
# 
# devtools::install_github("GraceYoon/SPRING")
# devtools::install_github("zdk123/SpiecEasi")
# 
# devtools::install_github("stefpeschel/NetCoMi", 
#                          dependencies = c("Depends", "Imports"),
#                          repos = c("https://cloud.r-project.org/",
#                                    BiocManager::repositories()))
# 

# loading NetCoMi library
library(NetCoMi)

# loading functions
source("Programming/extracting_data_NATH.R")
source("Programming/data_filtering.R")

# loading data
chosenDataSet     = "metab"       # "metab" or "genom"
chosenTaxonomy    <- "species"    # "species" or "genus"   
chosenWeek        <- "null"       # "1", "4", or "10"
chosenCutoffMass  <- 200          # arbitrary value, removing based on column name
inData <- extracting_data_NATH(whichWeek=chosenWeek, whichTaxLevel=chosenTaxonomy, cutOffMetabMass=chosenCutoffMass)

# filtering data
inData <- data_filtering(data=inData, whichDataSet=chosenDataSet, whichWeek=chosenWeek)

# inData_week1 <- data_filtering(data = extracting_data_NATH(whichWeek = "1"), whichDataSet = "metab", whichWeek = 1)


### for choosing the different days, with the maximum filtering of the full data set

# can only select week number 1, 4, or 10
choiceOfWeekHere <- "null"
if (choiceOfWeekHere != "null") {
  data = inData[gsub(".+-(?=\\d+$)", "", rownames(inData), perl = TRUE)==choiceOfWeekHere,]  
} else {
  data = inData  
}

# or if wanted, make your own data here: 
# data_week1  = inData[gsub(".+-(?=\\d+$)", "", rownames(inData), perl = TRUE)=="1",]
# data_week4  = inData[gsub(".+-(?=\\d+$)", "", rownames(inData), perl = TRUE)=="4",]
# data_week10 = inData[gsub(".+-(?=\\d+$)", "", rownames(inData), perl = TRUE)=="10",]



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


############################ 
# building single network with spearman as association measure - Full dataset, both treated and untreated
net_single_fullSet <- netConstruct((data),
                                   #filtTax = "highestFreq",
                                   #filtTaxPar = list(highestFreq = 100),
                                   #filtSamp = "totalReads",
                                   #filtSampPar = list(totalReads = 1000),
                                   measure = "spearman",thresh = 0.65,
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
     shortenLabels = "none",
     #nodeFilter = "clustMin",
     #nodeFilter = "highestBetween",
     #nodeFilterPar = 50,
     cexLabels = 1.3,
     title1 = paste("Single network with Spearman\nWeek:", choiceOfWeekHere, "taxonomy:", chosenTaxonomy),
     #title1 = "Single network with Spearman",
     showTitle = T,
     cexTitle = 1.7)
#nodeColor = colVector)

plot.microNetProps


# adding legend
#legend("topright", cex = 0.5, title = "estimated association:",
#      legend = c("Metataxonomic","Genomic"), lty = 1, lwd = 1, col = c("green","red"), 
#     bty = "n", horiz = TRUE)

legend(x=0.85,y=0.9,legend=c("Positive","Negative"),
       cex=0.6,col=c("green","red"),pch=c(".","."),lwd = c(3,3))
legend(x=0.85,y=0.6,legend=c("Metataxonomic","Metabolomic"),
       cex=0.6,col=c("green","red"),pch=c(16,16),lwd = c(3,3))



###

cor(log1p(data$Jannaschia[-13]),log1p(data$`943.99194`[-13]),method = "pear" )^2
plot(log1p(data$Jannaschia[-13]),log1p(data$`943.99194`[-13]))
LM=lm(log1p(data$Jannaschia[-13])~log1p(data$`943.99194`[-13]))
summary(LM)
plot(LM)
###

#cor(inData$TolC1, inData$Escherichia.Shigella, method = "spear")



############################ 
# compare two networks differentiated upon presence of TDA or not
# splitting the data set of all weeks into two; TDA and noTDA
data_TDA    <- data[substr(rownames(data), 1, 1) == "D",]
data_noTDA  <- data[substr(rownames(data), 1, 1) == "P",]

# Network construction
net_TDA <- netConstruct(data = data_noTDA, 
                        data2 = data_TDA,  
                        filtTax = "highestVar",
                        filtTaxPar = list(highestVar = 50),
                        measure = "spearman", thresh = 0.65,
                        measurePar = list(nlambda=10, 
                                         rep.num=10),
                        normMethod = "none", 
                        zeroMethod = "none",
                        sparsMethod = "threshold", 
                        dissFunc = "signed",
                        verbose = 3, weighted = T,
                        seed = 123456)

props_TDA <- netAnalyze(net_TDA, 
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

summary(props_TDA)

# for saving the plot as an image
# png(filename = "TDA_Vs_NoTDA_FullSet.png",
#     width = 4000, height = 3000,units = "px", pointsize = 12,
#     bg = "white", res = 300, family = "", restoreConsole = TRUE,
#     type = c("windows", "cairo", "cairo-png"),
#     symbolfamily = "default")

plot(props_TDA, 
     sameLayout = TRUE, 
     nodeColor = "cluster",
     nodeSize = "mclr",
     labelScale = FALSE,
     shortenLabels = "none",
     cexNodes = 1.5, 
     cexLabels = 1.3,
     cexHubLabels = 1,
     cexTitle = 3.7,
     groupNames = c("No TDA", "TDA"),
     hubBorderCol  = "gray40")
# dev.off() # shutting off image saving

legend("bottomleft", title = "estimated association:", legend = c("+","-"), 
       col = c("#009900","red"), inset = 0.02, cex = 2, lty = 1, lwd = 4, 
       bty = "n", horiz = TRUE)


comp_TDA <- netCompare(props_TDA, permTest = FALSE, verbose = FALSE)

summary(comp_TDA, 
        groupNames = c("No TDA", "TDA"),
        showCentr = c("degree", "between", "closeness"), 
        #showCentr = c("eigenvector"),
        numbNodes = 5)

############################ 
# compare two networks differentiated upon presence of TDA or not
# splitting the data set of the chosen week into two; TDA and noTDA


# can only select week number 1, 4, or 10
choiceOfWeekHere <- "1"
if (choiceOfWeekHere != "null") {
  data_weekly = inData[gsub(".+-(?=\\d+$)", "", rownames(inData), perl = TRUE)==choiceOfWeekHere,]  
} else {
  data_weekly = inData  
}


if (chosenTaxonomy=="genus"){
  colnames(data_weekly)=gsub("DATA.Bacteria_[A-Za-z]*[_.][A-Za-z]*[_.][A-Za-z]*[_.][A-Za-z]*[_.](.*)","\\1",colnames(data_weekly))
  colnames(data_weekly)=gsub("DATA.Archaea_[A-Za-z]*[_.][A-Za-z]*[_.][A-Za-z]*[_.][A-Za-z]*[_.](.*)","\\1",colnames(data_weekly))
} else {
  colnames(data_weekly)=gsub("DATA.Bacteria_[A-Za-z]*[_.][A-Za-z]*[_.][A-Za-z]*[_.][A-Za-z]*[_.](.*)[_.](.*)","\\1",colnames(data_weekly))
  colnames(data_weekly)=gsub("DATA.Archaea_[A-Za-z]*[_.][A-Za-z]*[_.][A-Za-z]*[_.][A-Za-z]*[_.](.*)[_.](.*)","\\1",colnames(data_weekly))
}


data_weekly_TDA    <- data_weekly[substr(rownames(data_weekly), 1, 1) == "D",]
data_weekly_noTDA  <- data_weekly[substr(rownames(data_weekly), 1, 1) == "P",]

# Network construction
net_weekly_TDA <- netConstruct(data = data_weekly_noTDA, 
                        data2 = data_weekly_TDA,  
                        #filtTax = "highestVar",
                        #filtTaxPar = list(highestVar = 50),
                        measure = "spearman", thresh = 0.8,
                        measurePar = list(nlambda=10, 
                                          rep.num=10),
                        normMethod = "none", 
                        zeroMethod = "none",
                        sparsMethod = "threshold", 
                        dissFunc = "signed",
                        verbose = 3, weighted = T,
                        seed = 123456)

props_weekly_TDA <- netAnalyze(net_weekly_TDA, 
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

#summary(props_weekly_TDA)

# using png to save plot
# png(filename = paste("TDA_Vs_NoTDA_Week_",choiceOfWeekHere,".png"), 
#     width = 4000, height = 3000,units = "px", pointsize = 12,
#     bg = "white", res = 300, family = "", restoreConsole = TRUE,
#     type = c("windows", "cairo", "cairo-png"),
#     symbolfamily = "default")


plot(props_weekly_TDA, 
     sameLayout = TRUE, 
     nodeColor = "cluster",
     nodeSize = "mclr",
     labelScale = FALSE,
     shortenLabels = "none",
     cexNodes = 1.5, 
     cexLabels = 1.3,
     cexHubLabels = 1,
     # showTitle = TRUE,
     # title2 = paste("Week",choiceOfWeekHere),
     #title(main = paste("Week",choiceOfWeekHere)),
     #cexTitle = 3.7,
     cexTitle = 2.5,
     groupNames = c(paste("No TDA\nWeek:",choiceOfWeekHere), paste("TDA\nWeek:",choiceOfWeekHere)),
     hubBorderCol  = "gray40")

# using png() and dev.off to save the plot
# dev.off()


# legend("bottomleft", title = "estimated association:", legend = c("+","-"), 
#        col = c("#009900","red"), inset = 0.02, cex = 2, lty = 1, lwd = 4, 
#        bty = "n", horiz = TRUE)


comp_weekly_TDA <- netCompare(props_weekly_TDA, permTest = FALSE, verbose = FALSE)

summary(comp_weekly_TDA, 
        groupNames = c(paste("No TDA\nWeek:",choiceOfWeekHere),paste("TDA\nWeek:",choiceOfWeekHere)),
        showCentr = c("degree", "between", "closeness"), 
        #showCentr = c("eigenvector"),
        numbNodes = 5)
