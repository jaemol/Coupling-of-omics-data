### Script to implement NetCoMi network analysis ### 
### Using Git-package NetCoMi ### 
### https://github.com/stefpeschel/NetCoMi ###

# loading NetCoMi library
library(NetCoMi)

# loading functions
source("Programming/extracting_data_NATH.R")
source("Programming/data_filtering.R")
source("Programming/Functions.R")

# loading data
chosenDataSet       = "metab"       # "metab" or "genom"
chosenTaxonomy      <- "species"    # "species" or "genus"   
chosenWeek          <- "null"       # "1", "4", or "10"
chosenCutoffMass    <- 200          # arbitrary value, removing based on column name
chosenNormalization <- "peak"        # can either be 'mad', 'median', 'mean' or 'peak'
inData <- extracting_data_NATH(whichWeek=chosenWeek, whichTaxLevel=chosenTaxonomy, 
                               cutOffMetabMass=chosenCutoffMass, whichNormalization=chosenNormalization)

# filtering data
chosenCutoffFiltering <- 0.78
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

# run this for to simplify the naming of the 16s rRNA sequencing data 
if (chosenTaxonomy=="genus"){
  colnames(data)=gsub("DATA.Bacteria_[A-Za-z]*[_.][A-Za-z]*[_.][A-Za-z]*[_.][A-Za-z]*[_.](.*)","\\1",colnames(data))
  colnames(data)=gsub("DATA.Archaea_[A-Za-z]*[_.][A-Za-z]*[_.][A-Za-z]*[_.][A-Za-z]*[_.](.*)","\\1",colnames(data))
} else {
  colnames(data)=gsub("DATA.Bacteria_[A-Za-z]*[_.][A-Za-z]*[_.][A-Za-z]*[_.][A-Za-z]*[_.](.*)[_.](.*)","\\1",colnames(data))
  colnames(data)=gsub("DATA.Archaea_[A-Za-z]*[_.][A-Za-z]*[_.][A-Za-z]*[_.][A-Za-z]*[_.](.*)[_.](.*)","\\1",colnames(data))
}


############################ 
# compare two networks differentiated upon presence of TDA or not
# splitting the data set of all weeks into three; TDA, noTDA and control
data_TDA      <- data[substr(rownames(data), 1, 1) == "D",]
data_noTDA    <- data[substr(rownames(data), 1, 1) == "P",]
data_control  <- data[substr(rownames(data), 1, 1) == "C",]

chosenThreshold <- 0.45

# Network construction - noTDAvsControl
net_noTDA_control <- netConstruct(data = data_noTDA, 
                        data2 = data_control,  
                        filtTax = "highestVar",
                        filtTaxPar = list(highestVar = 76),
                        measure = "spearman", thresh = chosenThreshold,
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
                        hubQuant = 0.9,
                        lnormFit = TRUE,
                        normDeg = FALSE,
                        normBetw = FALSE,
                        normClose = FALSE,
                        normEigen = FALSE)

summary(props_noTDA_control)

# for saving the plot as an image
# png(filename = paste("NoTDA_Vs_Control_thres",chosenThreshold,chosenNormalization,".png"),
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
                        filtTaxPar = list(highestVar = 76),
                        measure = "spearman", thresh = chosenThreshold,
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
# png(filename = paste("TDA_Vs_Control_thres",chosenThreshold,chosenNormalization,".png"),
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
                                filtTaxPar = list(highestVar = 76),
                                measure = "spearman", thresh = chosenThreshold,
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
# png(filename = paste("TDA_Vs_NoTDA_thres",chosenThreshold,chosenNormalization,".png"),
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
feat1 <- "Cohaesibacter"
feat2 <- "257.1857486"
data_analyze(data = inData, feature1 = feat1, feature2 = feat2)















#############################
# making a frontpage-figure # 
#############################

net_single_front <- netConstruct(data_noTDA,
                                   #filtTax = "highestFreq",
                                   #filtTaxPar = list(highestFreq = 100),
                                   #filtSamp = "totalReads",
                                   #filtSampPar = list(totalReads = 1000),
                                   filtTax = "highestVar",
                                   filtTaxPar = list(highestVar = 50),
                                   measure = "spearman",thresh = chosenThreshold,
                                   measurePar = list(nlambda=10, 
                                                     rep.num=10),
                                   normMethod = "none", 
                                   zeroMethod = "none",
                                   sparsMethod = "threshold", 
                                   dissFunc = "signed",
                                   verbose = 3, weighted = T,
                                   seed = 123456)

props_single_front <- netAnalyze(net_single_front, 
                                   centrLCC = TRUE,
                                   clustMethod = "cluster_fast_greedy",
                                   hubPar = "eigenvector",
                                   weightDeg = FALSE, normDeg = FALSE)

plot(props_single_front,
     labelScale = F,
     cexLabels = 0,
     #title1 = paste("Single network with Spearman, treated", chosenWeek, chosenTaxonomy),
     #showTitle = T,
     nodeSize = "mclr",
     shortenLabels = "none",
     cexNodes = 1.5,
     #cexTitle = 2.3)
#nodeColor = colVector)
