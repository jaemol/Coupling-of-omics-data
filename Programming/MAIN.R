### This is the main programming script for handling the data ### 
rm(list = ls())

# loading libraries
#library(dplyr)
library(gsubfn)

# loading in the sources of other functions

## Labtop
#source("C:/Users/Jmoll/Documents/GitHub/Coupling-of-omics-data/Programming/extracting_data_KAT.R")
#source("C:/Users/Jmoll/Documents/GitHub/Coupling-of-omics-data/Programming/data_filtering.R")

## Stationary
source("Programming/extracting_data_KAT.R")
source("Programming/data_filtering.R")

# loading in the data
list[data, len_16s, len_qpcr] <- extracting_data_KAT()

list[data, testID, OUA]       <- data_filtering(data)

