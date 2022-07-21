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

# install.packages("devtools")
devtools::install_github("GraceYoon/SPRING")
devtools::install_github("zdk123/SpiecEasi")

devtools::install_github("stefpeschel/NetCoMi",
                         dependencies = c("Depends", "Imports"),
                         repos = c("https://cloud.r-project.org/",
                                   BiocManager::repositories()))


# loading NetCoMi library
library(NetCoMi)

# loading functions
source("Programming/extracting_data_NATH.R")
source("Programming/data_filtering.R")
source("Programming/data_analyze.R")

# loading data
chosenDataSet       = "metab"       # "metab" or "genom"
chosenTaxonomy      <- "species"    # "species" or "genus"   
chosenWeek          <- "null"       # "1", "4", or "10"
chosenCutoffMass    <- 200          # arbitrary value, removing based on column name
chosenNormalization <- "mad"        # can either be 'mad', 'median' or 'peak'
inData <- extracting_data_NATH(whichWeek=chosenWeek, whichTaxLevel=chosenTaxonomy, 
                               cutOffMetabMass=chosenCutoffMass, whichNormalization=chosenNormalization)

# filtering data
chosenCutoffFiltering <- 0.85
inData <- data_filtering(data=inData, whichDataSet=chosenDataSet, whichWeek=chosenWeek, cutOffOrAuto=chosenCutoffFiltering)

### for choosing the different days, with the maximum filtering of the full data set

# can only select week number 1, 4, 6, or 10
choiceOfWeekHere <- "null"
if (choiceOfWeekHere != "null") {
  data = inData[gsub(".+-(?=\\d+$)", "", rownames(inData), perl = TRUE)==choiceOfWeekHere,]  
} else {
  data = inData  
}

# or if wanted, make your own data here: 
# data_week1  = inData[gsub(".+-(?=\\d+$)", "", rownames(inData), perl = TRUE)=="1",]
# data_week4  = inData[gsub(".+-(?=\\d+$)", "", rownames(inData), perl = TRUE)=="4",]
# data_week6  = inData[gsub(".+-(?=\\d+$)", "", rownames(inData), perl = TRUE)=="6",]
# data_week10 = inData[gsub(".+-(?=\\d+$)", "", rownames(inData), perl = TRUE)=="10",]



#dataOneColumns <- grep(x = colnames(data), pattern = "DATA.*")
# colVector <- 1:length(colnames(data))
# for (i in 1:length(colVector)) {
#   if (i<=max(dataOneColumns)) {colVector[i] = "green"} else {colVector[i] = "red"}
# }


if (chosenTaxonomy=="genus"){
  colnames(data)=gsub("DATA.Bacteria_[A-Za-z]*[_.][A-Za-z]*[_.][A-Za-z]*[_.][A-Za-z]*[_.](.*)","\\1",colnames(data))
  colnames(data)=gsub("DATA.Archaea_[A-Za-z]*[_.][A-Za-z]*[_.][A-Za-z]*[_.][A-Za-z]*[_.](.*)","\\1",colnames(data))
} else {
  colnames(data)=gsub("DATA.Bacteria_[A-Za-z]*[_.][A-Za-z]*[_.][A-Za-z]*[_.][A-Za-z]*[_.](.*)[_.](.*)","\\1",colnames(data))
  colnames(data)=gsub("DATA.Archaea_[A-Za-z]*[_.][A-Za-z]*[_.][A-Za-z]*[_.][A-Za-z]*[_.](.*)[_.](.*)","\\1",colnames(data))
}


############################ 
# building single network with spearman as association measure - Full dataset, both treated and untreated
# net_single_fullSet <- netConstruct((data),
#                                    #filtTax = "highestFreq",
#                                    #filtTaxPar = list(highestFreq = 100),
#                                    #filtSamp = "totalReads",
#                                    #filtSampPar = list(totalReads = 1000),
#                                    measure = "spearman",thresh = 0.65,
#                                    measurePar = list(nlambda=10, 
#                                                      rep.num=10),
#                                    normMethod = "none", 
#                                    zeroMethod = "none",
#                                    sparsMethod = "threshold", 
#                                    dissFunc = "signed",
#                                    verbose = 3,weighted = T,
#                                    seed = 123456)
# 
# props_single_fullSet <- netAnalyze(net_single_fullSet, 
#                                    centrLCC = TRUE,
#                                    clustMethod = "cluster_fast_greedy",
#                                    hubPar = "eigenvector",
#                                    weightDeg = FALSE, normDeg = FALSE)
# 
# #?summary.microNetProps
# summary(props_single_fullSet, numbNodes = 5L)
# 
# 
# plot(props_single_fullSet,
#      labelScale = F,
#      shortenLabels = "none",
#      #nodeFilter = "clustMin",
#      #nodeFilter = "highestBetween",
#      #nodeFilterPar = 50,
#      cexLabels = 1.3,
#      title1 = paste("Single network with Spearman\nWeek:", choiceOfWeekHere, "taxonomy:", chosenTaxonomy),
#      #title1 = "Single network with Spearman",
#      showTitle = T,
#      cexTitle = 1.7)
# #nodeColor = colVector)
# 
# plot.microNetProps
# 
# 
# # adding legend
# #legend("topright", cex = 0.5, title = "estimated association:",
# #      legend = c("Metataxonomic","Genomic"), lty = 1, lwd = 1, col = c("green","red"), 
# #     bty = "n", horiz = TRUE)
# 
# legend(x=0.85,y=0.9,legend=c("Positive","Negative"),
#        cex=0.6,col=c("green","red"),pch=c(".","."),lwd = c(3,3))
# legend(x=0.85,y=0.6,legend=c("Metataxonomic","Metabolomic"),
#        cex=0.6,col=c("green","red"),pch=c(16,16),lwd = c(3,3))





############################ 
# compare two networks differentiated upon presence of TDA or not
# splitting the data set of all weeks into three; TDA, noTDA and control
data_TDA      <- data[substr(rownames(data), 1, 1) == "D",]
data_noTDA    <- data[substr(rownames(data), 1, 1) == "P",]
data_control  <- data[substr(rownames(data), 1, 1) == "C",]

# Network construction - noTDAvsControl
net_noTDA_control <- netConstruct(data = data_noTDA, 
                        data2 = data_control,  
                        filtTax = "highestVar",
                        filtTaxPar = list(highestVar = 50),
                        measure = "spearman", thresh = 0.4,
                        measurePar = list(nlambda=10, 
                                         rep.num=10),
                        normMethod = "none", 
                        zeroMethod = "none",
                        sparsMethod = "threshold", 
                        dissFunc = "signed",
                        verbose = 3, weighted = T,
                        seed = 123456)

props_noTDA_control <- netAnalyze(net_noTDA_control, 
                        centrLCC = FALSE,
                        avDissIgnoreInf = TRUE,
                        sPathNorm = FALSE,
                        clustMethod = "cluster_fast_greedy",
                        #hubPar = c("degree", "between", "closeness"),
                        hubPar = "eigenvector",
                        hubQuant = 01.9,
                        lnormFit = TRUE,
                        normDeg = FALSE,
                        normBetw = FALSE,
                        normClose = FALSE,
                        normEigen = FALSE)

summary(props_noTDA_control)

# for saving the plot as an image
# png(filename = "noTDA_Vs_Control_thres40_mad.png",
#     width = 4000, height = 3000,units = "px", pointsize = 12,
#     bg = "white", res = 300, family = "", restoreConsole = TRUE,
#     type = c("windows", "cairo", "cairo-png"),
#     symbolfamily = "default")

plot(props_noTDA_control, 
     sameLayout = TRUE, 
     nodeColor = "cluster",
     nodeSize = "mclr",
     labelScale = FALSE,
     shortenLabels = "none",
     cexNodes = 1.5, 
     cexLabels = 1.3,
     cexHubLabels = 1,
     cexTitle = 3.7,
     groupNames = c("No TDA", "Control"),
     hubBorderCol  = "gray40")
#dev.off() # shutting off image saving

# legend("bottomleft", title = "estimated association:", legend = c("+","-"), 
#        col = c("#009900","red"), inset = 0.02, cex = 2, lty = 1, lwd = 4, 
#        bty = "n", horiz = TRUE)


comp_noTDA_control <- netCompare(props_noTDA_control, permTest = FALSE, verbose = FALSE)

summary(comp_noTDA_control, 
        groupNames = c("No TDA", "TDA"),
        showCentr = c("degree", "between", "closeness"), 
        #showCentr = c("eigenvector"),
        numbNodes = 5)

# Network construction - TDAvsControl
net_TDA_control <- netConstruct(data = data_TDA, 
                        data2 = data_control,  
                        filtTax = "highestVar",
                        filtTaxPar = list(highestVar = 50),
                        measure = "spearman", thresh = 0.4,
                        measurePar = list(nlambda=10, 
                                          rep.num=10),
                        normMethod = "none", 
                        zeroMethod = "none",
                        sparsMethod = "threshold", 
                        dissFunc = "signed",
                        verbose = 3, weighted = T,
                        seed = 123456)

props_TDA_control <- netAnalyze(net_TDA_control, 
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

summary(props_TDA_control)

# for saving the plot as an image
# png(filename = "TDA_Vs_Control_thres40_mad.png",
#     width = 4000, height = 3000,units = "px", pointsize = 12,
#     bg = "white", res = 300, family = "", restoreConsole = TRUE,
#     type = c("windows", "cairo", "cairo-png"),
#     symbolfamily = "default")

plot(props_TDA_control, 
     sameLayout = TRUE, 
     nodeColor = "cluster",
     nodeSize = "mclr",
     labelScale = FALSE,
     shortenLabels = "none",
     cexNodes = 1.5, 
     cexLabels = 1.3,
     cexHubLabels = 1,
     cexTitle = 3.7,
     groupNames = c("TDA", "Control"),
     hubBorderCol  = "gray40")
#dev.off() # shutting off image saving

# legend("bottomleft", title = "estimated association:", legend = c("+","-"), 
#        col = c("#009900","red"), inset = 0.02, cex = 2, lty = 1, lwd = 4, 
#        bty = "n", horiz = TRUE)


comp_TDA_control <- netCompare(props_TDA_control, permTest = FALSE, verbose = FALSE)

summary(comp_TDA_control, 
        groupNames = c("TDA", "Control"),
        showCentr = c("degree", "between", "closeness"), 
        #showCentr = c("eigenvector"),
        numbNodes = 5)

# Network construction - TDAvsControl
net_TDA_noTDA <- netConstruct(data = data_TDA, 
                                data2 = data_noTDA,  
                                filtTax = "highestVar",
                                filtTaxPar = list(highestVar = 50),
                                measure = "spearman", thresh = 0.4,
                                measurePar = list(nlambda=10, 
                                                  rep.num=10),
                                normMethod = "none", 
                                zeroMethod = "none",
                                sparsMethod = "threshold", 
                                dissFunc = "signed",
                                verbose = 3, weighted = T,
                                seed = 123456)

props_TDA_noTDA <- netAnalyze(net_TDA_noTDA, 
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

summary(props_TDA_noTDA)

# for saving the plot as an image
# png(filename = "TDA_Vs_NoTDA_thres40_mad.png",
#     width = 4000, height = 3000,units = "px", pointsize = 12,
#     bg = "white", res = 300, family = "", restoreConsole = TRUE,
#     type = c("windows", "cairo", "cairo-png"),
#     symbolfamily = "default")

plot(props_TDA_noTDA, 
     sameLayout = TRUE, 
     nodeColor = "cluster",
     nodeSize = "mclr",
     labelScale = FALSE,
     shortenLabels = "none",
     cexNodes = 1.5, 
     cexLabels = 1.3,
     cexHubLabels = 1,
     cexTitle = 3.7,
     groupNames = c("TDA", "No TDA"),
     hubBorderCol  = "gray40")
#dev.off() # shutting off image saving

# legend("bottomleft", title = "estimated association:", legend = c("+","-"), 
#        col = c("#009900","red"), inset = 0.02, cex = 2, lty = 1, lwd = 4, 
#        bty = "n", horiz = TRUE)


comp_TDA_noTDA <- netCompare(props_TDA_noTDA, permTest = FALSE, verbose = FALSE)

summary(comp_TDA_noTDA, 
        groupNames = c("TDA", "No TDA"),
        showCentr = c("degree", "between", "closeness"), 
        #showCentr = c("eigenvector"),
        numbNodes = 5)



#####
# trying data_analyze
# henriciella + 770.4851152
feat1 <- "Jannaschia"
feat2 <- "645.2042601"
data_analyze(data = inData, feature1 = feat1, feature2 = feat2)
