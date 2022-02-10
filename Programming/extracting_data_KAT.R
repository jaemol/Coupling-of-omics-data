### Script to extract data from testdatasheet from Katrine ###
### Only want the data from e.g. week 2 and the appropriate columns ###

# loading libraries
#library(dplyr)
library(gsubfn)

# beginnning function # 
extracting_data_KAT <- function() {
  # loading in the two datasheets
  data_16s  <- read.table("G:/Mit drev/Specialeprojekt/Testdata/allData_16S_cleaned.txt")
  data_qpcr <- read.table("G:/Mit drev/Specialeprojekt/Testdata/allData_qPCR_cleaned.txt")
  
  # Finding just the results from week two
  data_week02_16s   <- data_16s[data_16s$SAMPLEWEEK == 'Week 02',]
  data_week02_qpcr  <- data_qpcr[data_qpcr$SAMPLEWEEK == 'Week 02',]
  
  # next, include only testID, OUA and the results of the analyses
  pure_data_16s   <- subset(data_week02_16s, select = -c(PIG, DATE, WEANING_TIME, GROUP, Newlytreated, 
                                                   WEEK, SAMPLEWEEK, Florkem, Metacam, Zactran, 
                                                   Antibiotic, Treatment_group, Treatment_date, CorrectedGroup))
  pure_data_qpcr  <- subset(data_week02_qpcr, select = -c(PIG, DATE, WEANING_TIME, GROUP, Newlytreated, 
                                                      WEEK, SAMPLEWEEK, Florkem, Metacam, Zactran, 
                                                      Antibiotic, Treatment_group, Treatment_date, CorrectedGroup, sample))
  
  # extract two arrays of testIDs, and find which are common between them both
  testID_16s  <- pure_data_16s$PIG_DATE
  testID_qpcr <- pure_data_qpcr$PIG_DATE
  
  common_IDs  <- intersect(testID_16s, testID_qpcr)
  
  # only keep relevant testIDs
  pure_data_16s   = pure_data_16s[pure_data_16s$PIG_DATE %in% common_IDs, ]
  pure_data_qpcr  = pure_data_qpcr[pure_data_qpcr$PIG_DATE %in% common_IDs, ]
  
  # sort the datasets based testID, append whilst removing 1 edition of testID and OUA
  pure_data_16s   = pure_data_16s[order(pure_data_16s$PIG_DATE), ]
  pure_data_qpcr  = pure_data_qpcr[order(pure_data_qpcr$PIG_DATE), ]
  pure_data_qpcr  = subset(pure_data_qpcr, select = -c(PIG_DATE, OUA))
  
  # now, append the datasets, to get one set
  complete_data <- cbind(pure_data_16s, pure_data_qpcr, deparse.level = 1)
  
  # keeping the length of 16s and qpcr data - JUST the results, no ID or OUA
  len_16s   <- length(pure_data_16s[1,]) - 2
  len_qpcr  <- length(pure_data_qpcr[1,])
  
  # dropping the variables, except the resulting data set
  rm(list=setdiff(ls(), c("complete_data", "len_16s", "len_qpcr")))
  
  list(complete_data, len_16s, len_qpcr)
}
