data_analyze <- function(data, feature1, feature2, sampleRemovalOutliers="null") {
  # finding the correct features from the data set
  feat1 = data[which(stringr::str_split_fixed(string = colnames(data), 
                                         pattern = "[_]",7)==feature1, arr.ind = TRUE)[1]]
  feat2 = data[which(stringr::str_split_fixed(string = colnames(data), 
                                         pattern = "[_]",7)==feature2, arr.ind = TRUE)[1]]
  
  if (sampleRemovalOutliers != "null" && length(sampleRemovalOutliers)>0) {
    feat1 = feat1[-sampleRemovalOutliers,]
    feat2 = feat2[-sampleRemovalOutliers,]
  } 
  
  # finding the results and printing them out
  paste("Generating results for features", feature1, "and", feature2)
  paste("The ordinal spearman correlation is:", cor(feat1, feat2, method = "spear"))
  paste("The log1p spearman correlation is:", cor(log1p(feat1), log1p(feat2), method = "spear"))
  paste("The ordinal pearson correlation is:", cor(feat1, feat2, method = "pear"))
  paste("The log1p pearson correlation is:", cor(log1p(feat1), log1p(feat2), method = "pear"))
  print("Plotting the correlation between the two features\n")
  print("Building a linear model between the two features")
  LM = lm(log1p(feat1)~log1p(feat2))
  summary(LM)
  print("Plotting the related plots to the linear model")
  plot(LM)
}