### Script to extract data from Nathalie and Scott ###
### Metataxonomic and metabolomic data ###

# loading libraries
library(phyloseq)
library(stringr)
library(gsubfn)

# loading data
#data_phys_original  <- readRDS("Data/allDataMetataxonomicNCLTEE.rds")
load("Data/ps.asv.reduced.wTree.RData")
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
whichTaxLevel <- "species"

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


df.full.metax <- cbind(array.phaebac.bin,df.fulldata.metax)


#### Metabolomic data

# char_data <- read.csv("Data/metabolomic_day7,28,70.csv", stringsAsFactors = F)
# num_data <- data.frame(data.matrix(char_data))
# numeric_columns <- sapply(num_data,function(x){mean(as.numeric(is.na(x)))<0.5})
# final_data <- data.frame(num_data[,numeric_columns], char_data[,!numeric_columns])


# loading in data
df_metab_original <- read.csv("Data/metabolomic_day7,28,70.csv", header = TRUE, 
                              sep = ";", stringsAsFactors = FALSE, na.strings = "NA", strip.white = TRUE)

df_metab_numOnly = read.csv("Data/metabolomic_day7,28,70_NumOnly.csv", header = FALSE,
                            sep = ";", stringsAsFactors = FALSE, strip.white = TRUE)

# naming columns and rows
colnames(df_metab_numOnly) = colnames(df_metab_original)[-1]
rownames(df_metab_numOnly) = df_metab_original[,1][-1]

# transposing the data frame
df_metab_tmp1 = as.data.frame(t(df_metab_numOnly))

# loading in metadata metabolomics sheet - changed data to .csv first, to make it work
metadata_metabolomics  <- read.csv("Data/Metadata-metabolomics.csv", fill = TRUE, header = TRUE, sep = ";")

# fetching two first rows / names
rownam_samples_metab  <- rownames(df_metab_tmp1)
addit_info_samples    <- df_metab_original[,1]

# running through the names, matching them with phyloseq naming
metab_new_names <- rownam_samples_metab

for (i in 1:length(rownam_samples_metab)) {
  
  num_grep = as.numeric(unlist(regmatches(rownam_samples_metab[i], gregexpr("[[:digit:]]+", rownam_samples_metab[i]))))
    
  if (length(num_grep) == 1 && num_grep > 400) { 
    # find which sample is talked about
    metaDataRow = which(metadata_metabolomics$ï..Sample.no...MCCe. == num_grep)
    
    # checking for TDA
    if(metadata_metabolomics$System[metaDataRow]=="TDA"){tdaBin="P"}else{tdaBin="D"}
    
    # finding biorep
    biorepSample  = metadata_metabolomics$Bio.Rep[metaDataRow]
    
    # finding time
    timeSample    = metadata_metabolomics$Time[metaDataRow] / 7
    
    # inserting into new name format
    metab_new_names[i] = paste(tdaBin,biorepSample,timeSample, sep = "-")
    
  } 
}

# inserting the correct names
rownames(df_metab_tmp1) = metab_new_names

# removing the additional column with information (1)
df_metab_tmp2 = df_metab_tmp1[,-c(1)]

# finding the common test IDs, to make a full dataset
commonIDs <- intersect(sampleID.metatax, metab_new_names)

#data_metab = cbind(metab_new_names,df_metab_noblanks_tmp3)
#data_metax = cbind(sampleID.metatax,df.full.metax)

data_metab = df_metab_tmp2
data_metax = df.full.metax
  
# only keeping the relevant testIDs
df_metab = data_metab[metab_new_names %in% commonIDs,]
df_metax = data_metax[sampleID.metatax %in% commonIDs,]

# normalizing the metabolomic data, percentage-based according to max peak overall
maxPeak <- max(df_metab)
df_metab_tmp3 = as.data.frame(lapply(df_metab, function(x){x/maxPeak}))
colnames(df_metab_tmp3) = colnames(df_metab)
rownames(df_metab_tmp3) = rownames(df_metab)


df_metab = df_metab_tmp3[sort(commonIDs, decreasing = FALSE),]
df_metax = df_metax[sort(commonIDs, decreasing = FALSE),]
 
complete_data <- cbind(df_metax, df_metab, deparse.level = 1)


# jogging around, finding max metabolite peak
#maxValuesMetab <- unlist(apply(df_metab, 2, max))
#max(maxValuesMetab)










