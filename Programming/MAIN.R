### This is the main programming script for handling the data ### 
rm(list = ls())

# loading libraries
#library(dplyr)
library(gsubfn)

# loading in the sources of other functions

source("C:/Users/Jmoll/Documents/GitHub/Coupling-of-omics-data/Programming/extracting_data_KAT.R")
#source("C:/Users/Jmoll/Documents/GitHub/Coupling-of-omics-data/Programming/data_filtering.R")


# loading in the data
list[data, len_16s, len_qpcr] <- extracting_data_KAT()







# Great, now trying to cluster
# Finding distance matrix
#dist_mat_16s <- dist(clust_16s, method = 'euclidean')

# Fitting Hierarchical clustering Model
#set.seed(240) # Setting seed
#hierar_cl <- hclust(dist_mat_16s, method = "average") 

# Plotting dendrogram
#plot(hierar_cl)

# For fun, making a scatterplot
#scatterplotMatrix(data_week02_16s[20:30])

