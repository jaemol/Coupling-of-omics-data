### Script to extract data from Nathalie and Scott ###
### Metataxonomic and metabolomic data ###

library(phyloseq)
library(stringr)

data_phys_original  <- readRDS("Data/allDataMetataxonomicNCLTEE.rds")

# The actual data
df.otu.metatax      <- as.data.frame(data_phys_original@otu_table)

# The results / variables
df.tax.metatax      <- as.data.frame(data_phys_original@tax_table)

# The sample IDs
sampleID.metatax    <- data_phys_original@sam_data$sample.name

## setting up dataframe similar to Katrines
df.full.metax <- as.data.frame(t(df.otu.metatax))
rownames(df.full.metax) = sampleID.metatax

# making string array with colnames
colnames_array_metatax <- str_c("DATA.",str_replace_na(df.tax.metatax$domain, replacement="NA"),"_",
                                str_replace_na(df.tax.metatax$phylum,  replacement="NA"),"_",
                                str_replace_na(df.tax.metatax$class,   replacement="NA"),"_",
                                str_replace_na(df.tax.metatax$order,   replacement="NA"),"_",
                                str_replace_na(df.tax.metatax$family,  replacement="NA"),"_",
                                str_replace_na(df.tax.metatax$genus,   replacement="NA"),".",
                                str_replace_na(df.tax.metatax$species, replacement="NA"))
                          
colnames(df.full.metax) = colnames_array_metatax
