### Principal component analysis ### 
# Data needs to transposed..!
# Using 'spectral decomposition', examining covariances / correlation between variables
# Using R-function 'princomp()'

# loading libraries
library(pbkrtest)
library(ggplot2)
library(ggpubr)
library(factoextra) # used to create ggplot2-based visualization
library(data.table) # used to keep structure of data frame when transposing

# excluding testID and OUA (antibiotics used or not)
# loading functions
source("Programming/extracting_data_KAT.R")
source("Programming/data_filtering.R")

# loading data
chosenWeek      <- "Week 02"
chosenTaxonomy  <- "species"
inData <- extracting_data_KAT(whichWeek = chosenWeek, whichTaxLevel = chosenTaxonomy)

# filtering data
inData <- data_filtering(inData)

# excluding testID and OUA (antibiotics used or not)
#inData = data
testID  <- inData$testID
OUA     <- inData$OUA
data = subset(inData, select = -c(testID, OUA))

dataOneColumns <- grep(x = colnames(data), pattern = "DATA.*")


colVector <- 1:length(colnames(data))
for (i in 1:length(colVector)) {
  if (i<=max(dataOneColumns)) {colVector[i] = "green"} else {colVector[i] = "red"}
}


if (chosenTaxonomy=="genus"){
  colnames(data)=gsub("DATA.Bacteria_[A-Za-z]*[_.][A-Za-z]*[_.][A-Za-z]*[_.][A-Za-z]*[_.](.*)","\\1",colnames(data))
  colnames(data)=gsub("DATA.Archaea_[A-Za-z]*[_.][A-Za-z]*[_.][A-Za-z]*[_.][A-Za-z]*[_.](.*)","\\1",colnames(data))
} else {
  colnames(data)=gsub("DATA.Bacteria_[A-Za-z]*[_.][A-Za-z]*[_.][A-Za-z]*[_.][A-Za-z]*[_.](.*)[_.](.*)","\\1",colnames(data))
  colnames(data)=gsub("DATA.Archaea_[A-Za-z]*[_.][A-Za-z]*[_.][A-Za-z]*[_.][A-Za-z]*[_.](.*)[_.](.*)","\\1",colnames(data))
}

# transposing data; samples on columns, features on rows
#inData = transpose(inData) 

# adding the names to the transposed data
#rownames(inData) <- colnames(data[-(1:2)])
#colnames(inData) <- testID
rownames(inData) <- testID

# making PCA model
pcaModel <- prcomp(data, center = TRUE, scale = TRUE)

# making scree plot, to visualize eigenvalues
fviz_eig(pcaModel)

# making a graph of the variables
# visualizing the correlation between the variables
# positive correlated variables point to the same side of the plot
# negative correlated variables point to opposite sides of the plot
fviz_pca_var(pcaModel,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = FALSE     # Avoid text overlapping
             #,select.var = list(name = rownames(inData))
)

plot(pcaModel$x[,1], pcaModel$x[,2])
plot(pcaModel$x[,1], pcaModel$x[,3])



####################################
pca <- prcomp((data), scale=TRUE) 

## plot pc1 and pc2
plot(pca$x[,1], pca$x[,2])
plot(pca$x[,1], pca$x[,3])

## make a scree plot
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)

barplot(pca.var.per, main="Scree Plot", xlab="Principal Component", ylab="Percent Variation")
abline(h=1/ncol(inData)*100, col='red')
# other way! - making scree plot, to visualize eigenvalues
fviz_eig(pcaModel)

## now make a fancy looking plot that shows the PCs and the variation:
library(ggplot2)

pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,1],
                       Y=pca$x[,2])
pca.data

ggplot(data=pca.data, aes(x=X, y=Y, label=Sample)) +
  geom_text() +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
  theme_bw() +
  ggtitle("My PCA Graph")

## get the name of the top 10 measurements (genes) that contribute
## most to pc1.
loading_scores <- pca$rotation[,1]
gene_scores <- abs(loading_scores) ## get the magnitudes
gene_score_ranked <- sort(gene_scores, decreasing=TRUE)
top_10_genes <- names(gene_score_ranked[1:10])

top_10_genes ## show the names of the top 10 genes

pca$rotation[top_10_genes,1] ## show the scores (and +/- sign)
