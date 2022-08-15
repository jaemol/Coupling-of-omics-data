findEpsi=function(L, minRange=0,maxRange=3, steps=0.1, maxY=200, minP=10) {
  # Function, that helps determine the optimal epsilon for use in DBSCAN
  # optimal, meaning the one that will generate the most clusters,
  # with given data and minimal points in radius (epsilon)
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


getMetabDataNormEval=function(cutOffMetabMass=200){
  # Function, for extracting only metabolomic data, for evaluating the 
  # normalization done, to be used in PerMANOVA.
  # Output is unnormalized data frame and classes within data (groups) (tda, no tda, control + day of sampling)
  # loading in data
  library(gsubfn)
  df_metab_original <- read.csv("Data/allData_LCMS_metabolomics.csv", header = FALSE,
                                sep = ";", stringsAsFactors = FALSE, strip.white = TRUE, skip = 1)
  
  colnames(df_metab_original) = read.csv("Data/allData_LCMS_metabolomics.csv", header = F,
                                         sep = ";", stringsAsFactors = FALSE, strip.white = TRUE, nrows = 1)
  rownames(df_metab_original) = read.csv("Data/allData_LCMS_metabolomics.csv", header = F,
                                         sep = ";", stringsAsFactors = FALSE, strip.white = TRUE)[-1,1]
  
  # removing first column (masses)
  df_metab_original = df_metab_original[,-1]
  
  # transposing the data frame
  df_metab_tmp1 = as.data.frame(t(df_metab_original))
  
  # loading in metadata metabolomics sheet - changed data to .csv first, to make it work
  metadata_metabolomics  <- read.csv("Data/Metadata-metabolomics.csv", fill = TRUE, header = TRUE, sep = ";")
  
  # fetching first row
  rownam_samples_metab  <- rownames(df_metab_tmp1)
  
  # running through the names, matching them with phyloseq naming
  metab_new_names <- rownam_samples_metab
  groups    <- array(0, length(metab_new_names))
  days    <- array(0, length(metab_new_names))
  
  # getting sample names and groups
  for (i in 1:length(rownam_samples_metab)) {
    
    num_grep = as.numeric(unlist(regmatches(rownam_samples_metab[i], gregexpr("[[:digit:]]+", rownam_samples_metab[i]))))
    
    #if (length(num_grep) == 1 && num_grep > 400) { 
    # find which sample is talked about
    metaDataRow = which(metadata_metabolomics$Sample.no...MCCe. == num_grep)
    
    # ensuring no medium control (blank samples)
    #if (!(num_grep==449||num_grep==458)){
    
    # checking for TDA or control
    if(metadata_metabolomics$System[metaDataRow]=="TDA"){
      tdaBin="P"
      groups[i] <- "TDA"
    } else if (metadata_metabolomics$System[metaDataRow]=="NoTDA"){
      tdaBin="D"
      groups[i] <- "noTDA"
    } else {
      tdaBin="C"
      groups[i] <- "Control"
    }
    
    # finding biorep
    biorepSample  = metadata_metabolomics$Bio.Rep[metaDataRow]
    
    # finding time
    timeSample    = metadata_metabolomics$Time[metaDataRow] / 7
    days[i] <- gsub(" ","",paste("day",timeSample))
    
    # inserting into new name format
    metab_new_names[i] = paste(tdaBin,biorepSample,timeSample, sep = "-")
    
    #} 
  }
  
  # removing samples with NA
  whichNaNRemove <- which(stringr::str_split_fixed(metab_new_names, pattern = "-", 3)[,2] == "NA")
  metab_new_names <- metab_new_names[-whichNaNRemove]
  df_metab_tmp2   <- df_metab_tmp1[-whichNaNRemove,]
  groups    <- groups[-whichNaNRemove]
  days      <- days[-whichNaNRemove]
  
  # inserting the correct names
  rownames(df_metab_tmp2) = metab_new_names
  
  # removing all features with masses <= 200 m/Z (mass over charge)
  metabFeatToDrop <- which(as.numeric(colnames(df_metab_tmp2)) <= cutOffMetabMass)
  data_metab   <- subset(df_metab_tmp2, select = -c(metabFeatToDrop))
  
  list(data_metab, groups, days)
}

data_analyze <- function(data, feature1, feature2, sampleRemovalOutliers="null") {
  # finding the correct features from the data set
  feat1 = data[which(stringr::str_split_fixed(string = colnames(data), 
                                              pattern = "[_]",7)==feature1, arr.ind = TRUE)[1]]
  feat2 = data[which(stringr::str_split_fixed(string = colnames(data), 
                                              pattern = "[_]",7)==feature2, arr.ind = TRUE)[1]]
  
  if (sampleRemovalOutliers != "null" && length(sampleRemovalOutliers)>0) {
    feat1 = feat1[-sampleRemovalOutliers,]
    feat2 = feat2[-sampleRemovalOutliers,]
  } 
  
  # finding the results and printing them out
  print(paste("Generating results for features", feature1, "and", feature2))
  print(paste("The ordinal spearman correlation is:", cor(feat1, feat2, method = "spear")))
  print(paste("The log1p spearman correlation is:", cor(log1p(feat1), log1p(feat2), method = "spear")))
  print(paste("The ordinal pearson correlation is:", cor(feat1, feat2, method = "pear")))
  print(paste("The log1p pearson correlation is:", cor(log1p(feat1), log1p(feat2), method = "pear")))
  print("Plotting the correlation between the two features")
  plot(unlist(feat1) ~ unlist(feat2), ylab = feature1, xlab = feature2)
  print("Building a linear model between the two features")
  LM = lm(log1p(unlist(feat1))~log1p(unlist(feat2)))
  LM = lm(unlist(feat1)~unlist(feat2))
  print("Summary of the linear model")
  print(summary(LM))
  print("Plotting the related plots to the linear model")
  par(mfrow=c(2,2))
  plot(LM)
  par(mfrow=c(1,1))
}