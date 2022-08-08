###### normalization evaluation script ####
# This script is for importing 
# metabolomic (and metataxonomic?)
# data, and running perMANOVA to
# benchmark them and choosing the best 
###########################################

# loading libraries
library(vegan)
library(gsubfn)

# loading functions
source("Programming/extracting_data_NATH.R")
source("Programming/data_filtering.R")
source("Programming/data_analyze.R")
source("Programming/Functions.R")


### 
# the data has nestedness (split in groups), so a field is needed
# the field will be given values for the samples, based upon their group, so:
# 1 = control
# 2 = no TDA present
# 3 = TDA present
# and so, this need to be incorporated via 'STRATA' in the PerMANOVA

# getting the data
list[data_metab, strata_field]  <- getMetabDataNormEval(cutOffMetabMass = 200)

# implementing the different normalizations
data_peak   <- as.data.frame(apply(data_metab,MARGIN = 2, function(x){x/max(x)})); rownames(data_peak)=rownames(data_metab); colnames(data_peak)=colnames(data_metab)
data_median <- as.data.frame(apply(data_metab,MARGIN = 2, function(x){x/(median(x)+1)})); rownames(data_peak)=rownames(data_metab); colnames(data_peak)=colnames(data_metab)
data_mean   <- as.data.frame(apply(data_metab,MARGIN = 2, function(x){x/(mean(x)+1)})); rownames(data_peak)=rownames(data_metab); colnames(data_peak)=colnames(data_metab)
data_mad    <- as.data.frame(apply(data_metab,MARGIN = 2, function(x){x/(mad(x, center = median(x), na.rm = FALSE, constant = 1)+1)})); rownames(data_peak)=rownames(data_metab); colnames(data_peak)=colnames(data_metab)

# list[data_peak, strata_field]   <- getMetabDataNormEval(cutOffMetabMass = 200, whichNormalization = "peak")
# list[data_median, strata_field] <- getMetabDataNormEval(cutOffMetabMass = 200, whichNormalization = "median")
# list[data_mad, strata_field]    <- getMetabDataNormEval(cutOffMetabMass = 200, whichNormalization = "mad")
# list[data_mean, strata_field]    <- getMetabDataNormEval(cutOffMetabMass = 200, whichNormalization = "mean")

# adonis(Y ~ NO3, data=dat, strata=dat$field, perm=999)
adonis2(data_peak ~ strata_field, strata = NULL, permutations = 999, by = NULL, method = "euclidean")
adonis2(data_median ~ strata_field, strata = NULL, permutations = 999, by = NULL, method = "euclidean")
adonis2(data_mad ~ strata_field, strata = NULL, permutations = 999, by = NULL, method = "euclidean")
adonis2(data_mean ~ strata_field, strata = NULL, permutations = 999, by = NULL, method = "euclidean")


