### Script to extract data from Nathalie and Scott ###
### Metataxonomic and metabolomic data ###

# loading libraries
library(phyloseq)
library(stringr)
library(gsubfn)

# loading data
data_phys_original  <- readRDS("Data/allDataMetataxonomicNCLTEE.rds")

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



