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
# 3 = TDA presen
# and so, this need to be incorporated via 'STRATA' in the PerMANOVA

# getting the data
list[data_peak, strata_field]   <- getMetabDataNormEval(cutOffMetabMass = 200, whichNormalization = "peak")
list[data_median, strata_field] <- getMetabDataNormEval(cutOffMetabMass = 200, whichNormalization = "median")
list[data_mad, strata_field]    <- getMetabDataNormEval(cutOffMetabMass = 200, whichNormalization = "mad")


