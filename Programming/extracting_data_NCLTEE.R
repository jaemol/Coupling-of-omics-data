### Script to extract data from Nathalie and Scott ###
### Metataxonomic and metabolomic data ###

# loading libraries
library(phyloseq)
library(stringr)
library(gsubfn)

# beginning function
extracting_data_NCLTEE <- function(whichWeek="null", whichTaxLevel="species", 
                                 cutOffMetabMass=200, whichNormalization) {
  ################## Metataxonomic data ##################
  # loading data
  load("Data/allData_16S_NCLTEE_Reduced.RData")
  data_phys_original <- ps.new
  
  # The actual data
  df.otu.metatax      <- as.data.frame(data_phys_original@otu_table)
  
  # The results / variables
  df.tax.metatax      <- as.data.frame(data_phys_original@tax_table)
  
  # The sample IDs
  sampleID.metatax    <- data_phys_original@sam_data$sample.name
  
  ## setting up dataframe similar to Katrines
  df.fulldata.metax <- as.data.frame(t(df.otu.metatax))
  rownames(df.fulldata.metax) = sampleID.metatax
  
  # making string array with colnames
  colnames_array_metatax <- str_c("DATA.",str_replace_na(df.tax.metatax$domain, replacement="NA"),"_",
                                  str_replace_na(df.tax.metatax$phylum,  replacement="NA"),"_",
                                  str_replace_na(df.tax.metatax$class,   replacement="NA"),"_",
                                  str_replace_na(df.tax.metatax$order,   replacement="NA"),"_",
                                  str_replace_na(df.tax.metatax$family,  replacement="NA"),"_",
                                  str_replace_na(df.tax.metatax$genus,   replacement="NA"),"_",
                                  str_replace_na(df.tax.metatax$species, replacement="NA"))
                            
  colnames(df.fulldata.metax) = colnames_array_metatax
  
  # saving information on phaeobacter presence
  array.phaebac.bin <- data_phys_original@sam_data$phaeobacter
  
  # moving up to genus level or staying at species
  if (whichTaxLevel=="genus") {
    # moving up in taxonomy for the 16s data, going from species to genus
    origNames <- colnames(df.fulldata.metax)
    newNames <- apply(str_split_fixed(string = origNames, pattern = "[_]",7)[,1:6],1, paste, collapse="_")
    
    length(unique(newNames))
    uniqNames=unique(newNames)
    
    newDat=data.frame(dummy=1:NROW(df.fulldata.metax))
    
    #j=uniqNames[1]
    for(j in uniqNames) {
      
      jIndx=grep(j,origNames )
      if(length(jIndx)>1) {
        newDat=cbind(newDat,rowSums(df.fulldata.metax[,jIndx]))
      } else {
        newDat=cbind(newDat,(df.fulldata.metax[,jIndx]))
      }
    }
    
    newDat=newDat[,-1]
    colnames(newDat)=uniqNames
    
    # now making that new data into 16s data set
    df.fulldata.metax = newDat
  }
  df.full.metax <- df.fulldata.metax
  
  ################## Metabolomic data ##################
  # loading in data
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
  
  # loading in metadata metabolomics sheet
  metadata_metabolomics  <- read.csv("Data/Metadata-metabolomics.csv", fill = TRUE, header = TRUE, sep = ";")
  
  # fetching first row
  rownam_samples_metab  <- rownames(df_metab_tmp1)
  
  # running through the names, matching them with phyloseq naming
  metab_new_names <- rownam_samples_metab
  
  for (i in 1:length(rownam_samples_metab)) {
    
    num_grep = as.numeric(unlist(regmatches(rownam_samples_metab[i], gregexpr("[[:digit:]]+", rownam_samples_metab[i]))))

    # find which sample is relevant
    metaDataRow = which(metadata_metabolomics$Sample.no...MCCe. == num_grep)
    
    # checking for TDA, no TDA or control
    if(metadata_metabolomics$System[metaDataRow]=="TDA"){
      tdaBin="P"
    } else if (metadata_metabolomics$System[metaDataRow]=="NoTDA"){
      tdaBin="D"
    } else {tdaBin="C"}
    
    # finding biorep
    biorepSample  = metadata_metabolomics$Bio.Rep[metaDataRow]
    
    # finding time
    timeSample    = metadata_metabolomics$Time[metaDataRow] / 7
    
    # inserting into new name format
    metab_new_names[i] = paste(tdaBin,biorepSample,timeSample, sep = "-")
      
    #} 
  }
  
  # removing samples with NA
  whichNaNRemove <- which(stringr::str_split_fixed(metab_new_names, pattern = "-", 3)[,2] == "NA")
  metab_new_names <- metab_new_names[-whichNaNRemove]
  df_metab_tmp2   <- df_metab_tmp1[-whichNaNRemove,]
  
  
  
  # inserting the correct names
  rownames(df_metab_tmp2) = metab_new_names
  
  # finding the common test IDs, to make a full dataset
  commonIDs <- intersect(sampleID.metatax, metab_new_names)
  
  # removing all features with masses <= 200 m/Z (mass over charge)
  metabFeatToDrop <- which(as.numeric(colnames(df_metab_tmp2)) <= cutOffMetabMass)
  df_metab_tmp3   <- subset(df_metab_tmp2, select = -c(metabFeatToDrop))
  
  data_metab = df_metab_tmp3
  data_metax = df.full.metax
    
  # only keeping the relevant testIDs
  df_metab = data_metab[metab_new_names %in% commonIDs,]
  df_metax = data_metax[sampleID.metatax %in% commonIDs,]
  
  ## normalization of the metabolomic LC-MS data ##
  if (whichNormalization == "peak") {
    # normalizing the metabolomic data, percentage-based according to max peak per feature
    df_metab = as.data.frame(apply(df_metab_tmp3,MARGIN = 2, function(x){x/max(x)})) 
    
  } else if (whichNormalization == "median") {
    # normalizing per median
    df_metab = as.data.frame(apply(df_metab_tmp3, MARGIN = 2, function(x){x/(median(x)+1)}))
    
  } else if (whichNormalization == "mean") {
    # normalizing per mean
    df_metab = as.data.frame(apply(df_metab_tmp3, MARGIN = 2, function(x){x/(mean(x))}))
    
  } else if (whichNormalization == "mad") {
    # normalizing per mad (median absolute deviation)
    df_metab = as.data.frame(apply(df_metab_tmp3, MARGIN = 2, function(x){x/(mad(x, center = median(x), na.rm = FALSE, constant = 1)+1)}))
    
  } else {
    print("No correct normalization were chosen, no normalization performed.\n Enter either: 'median', 'mad' or 'peak'")
  }
  
  
  # giving the proper names to the data set
  colnames(df_metab) = colnames(df_metab)
  rownames(df_metab) = rownames(df_metab)
  
  # sorting for combining
  df_metab = df_metab[sort(commonIDs, decreasing = FALSE),]
  df_metax = df_metax[sort(commonIDs, decreasing = FALSE),]
   
  complete_data <- cbind(df_metax, df_metab, deparse.level = 1)
  
  # choosing a specific week of data if wanted
  if (whichWeek != "null") {
    # can only select week number 1, 4, or 10
    complete_data = complete_data[gsub(".+-(?=\\d+$)", "", rownames(complete_data), perl = TRUE)==whichWeek,]
  }
  
  return(complete_data)
}









