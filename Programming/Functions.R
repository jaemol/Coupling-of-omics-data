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


getMetabDataNormEval=function(cutOffMetabMass=200, whichNormalization){
  # Function, for extracting only metabolomic data, for evaluating the 
  # normalization done, to be used in PerMANOVA.
  # Output is normalized data frame and field (groups) (tda, no tda, control)
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
  strata_field    <- array(0, length(metab_new_names))
  
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
      strata_field[i] <- 3
    } else if (metadata_metabolomics$System[metaDataRow]=="NoTDA"){
      tdaBin="D"
      strata_field[i] <- 2
    } else {
      tdaBin="C"
      strata_field[i] <- 1
    }
    
    # finding biorep
    biorepSample  = metadata_metabolomics$Bio.Rep[metaDataRow]
    
    # finding time
    timeSample    = metadata_metabolomics$Time[metaDataRow] / 7
    
    # inserting into new name format
    metab_new_names[i] = paste(tdaBin,biorepSample,timeSample, sep = "-")
    
    #} 
  }
  
  # removing samples with NA
  #whichNaNRemove  <- which(gsub("MCCe", x = metab_new_names, replacement = "", perl = TRUE)==TRUE)
  whichNaNRemove  <- c(1,10)
  metab_new_names <- metab_new_names[-whichNaNRemove]
  df_metab_tmp2   <- df_metab_tmp1[-whichNaNRemove,]
  strata_field    <- strata_field[-whichNaNRemove]
  
  # inserting the correct names
  rownames(df_metab_tmp2) = metab_new_names
  
  # removing all features with masses <= 200 m/Z (mass over charge)
  metabFeatToDrop <- which(as.numeric(colnames(df_metab_tmp2)) <= cutOffMetabMass)
  df_metab_tmp3   <- subset(df_metab_tmp2, select = -c(metabFeatToDrop))
  
  ## normalization of the metabolomic LCMS data
  if (whichNormalization == "peak") {
    # normalizing the metabolomic data, percentage-based according to max peak per feature
    ### OBS COMMENT: Maybe we need to normalize per median ###
    data_metab = as.data.frame(apply(df_metab_tmp3,MARGIN = 2, function(x){x/max(x)})) 
    
  } else if (whichNormalization == "median") {
    # normalizing per median
    data_metab = as.data.frame(apply(df_metab_tmp3, MARGIN = 2, function(x){x/(median(x)+1)}))
    
  } else if (whichNormalization == "mean") {
    # normalizing per mean
    data_metab = as.data.frame(apply(df_metab_tmp3, MARGIN = 2, function(x){x/(mean(x)+1)}))
    
  } else if (whichNormalization == "mad") {
    # normalizing per mad (median absolute deviation)
    data_metab = as.data.frame(apply(df_metab_tmp3, MARGIN = 2, function(x){x/(mad(x, center = median(x), na.rm = FALSE, constant = 1)+1)}))
    
  } else {
    print("No correct normalization were chosen, no normalization performed.\n Enter either: 'median', 'mad' or 'peak'")
  }
  
  # giving the proper names to the data set
  colnames(data_metab) = colnames(df_metab_tmp3)
  rownames(data_metab) = rownames(df_metab_tmp3)
  
  # rm(list=setdiff(ls(), c("complete_data", "strata_field")))
  list(data_metab, strata_field)
}