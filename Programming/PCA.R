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
inData = data
testID  <- inData$testID
OUA     <- inData$OUA

inData = subset(inData, select = -c(testID, OUA))

# transposing data; samples on columns, features on rows
inData = transpose(inData) 

# adding the names to the transposed data
rownames(inData) <- colnames(data[-(1:2)])

# making PCA model
pcaModel <- prcomp(inData, scale = TRUE)

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

