### Script to implement NetCoMi network analysis ### 
### Using Git-package NetCoMi ### 
### https://github.com/stefpeschel/NetCoMi ###

# loading library
library(NetCoMi)

#install.packages("devtools")

devtools::install_github("stefpeschel/NetCoMi", dependencies = TRUE,
                         repos = c("https://cloud.r-project.org/",
                                   BiocManager::repositories()))


# manual installation
install.packages("stats") # containing pearson and spearman
install.packages("WGCNA") # containing bicor
install.packages("vegan") # containing vegdist
install.packages("LaplaceDemon") # containing KLD (divergences)
