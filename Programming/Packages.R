### Packages used in 'Coupling of omics data' ###
library(devtools)
library(BiocManager)
library(NetCoMi)
library(phyloseq)
library(stringr)
library(gsubfn)
library(vegan)
library(stringr)


## Installing NetCoMi for the first time 
if(!requireNamespace("BiocManager", quietly = TRUE)){
  utils::install.packages("BiocManager")
}

BiocManager::install(pkgs = c("Biobase", "doSNOW", "fdrtool", "filematrix",
                              "foreach", "graphics", "grDevices", "gtools",
                              "huge", "igraph", "MASS", "Matrix", "phyloseq",
                              "pulsar", "qgraph", "Rdpack", "snow", "SPRING",
                              "stats", "utils", "vegan", "WGCNA"))

BiocManager::install("GO.db")


devtools::install_github("GraceYoon/SPRING")
devtools::install_github("zdk123/SpiecEasi")

devtools::install_github("stefpeschel/NetCoMi", 
                         dependencies = c("Depends", "Imports"),
                         repos = c("https://cloud.r-project.org/",
                                   BiocManager::repositories()))