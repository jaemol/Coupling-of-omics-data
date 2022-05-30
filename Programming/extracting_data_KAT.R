### Script to extract data from testdatasheet from Katrine ###
### Only want the data from e.g. week 2 and the appropriate columns ###

# loading libraries
#library(dplyr)
library(gsubfn)

# beginning function # 
extracting_data_KAT <- function(whichWeek="null", loadOrigData=FALSE) {
  # loading in the two datasheets
  print("Loading in data...")
  if (loadOrigData == TRUE) {
    data_16s_original  <- read.table("Testdata/allData_16S_cleaned.txt")
    data_qpcr_original <- read.table("Testdata/allData_qPCR_cleaned.txt") 
  }
  
  if (whichWeek != "null"){
    # Finding just the results from week two
    data_16s   <- data_16s_original[data_16s_original$SAMPLEWEEK == whichWeek,]
    data_qpcr  <- data_qpcr_original[data_qpcr_original$SAMPLEWEEK == whichWeek,]
  } else {
    data_16s   <- data_16s_original
    data_qpcr  <- data_qpcr_original
  }
  
  print("Removing unimportant attributes...")
  # next, include only testID (PIG_DATE), OUA and the results of the analyses
  pure_data_16s   <- subset(data_16s, select = -c(PIG, DATE, WEANING_TIME, GROUP, Newlytreated, 
                                                   WEEK, SAMPLEWEEK, Florkem, Metacam, Zactran, 
                                                   Antibiotic, Treatment_group, Treatment_date, CorrectedGroup,
                                                   PIG.1, DATE.1, OUA.1, PIG_DATE.1, OriginalNAME))
  pure_data_qpcr  <- subset(data_qpcr, select = -c(PIG, DATE, WEANING_TIME, GROUP, Newlytreated, 
                                                      WEEK, SAMPLEWEEK, Florkem, Metacam, Zactran, 
                                                      Antibiotic, Treatment_group, Treatment_date, CorrectedGroup, sample))
  
  print("Only keeping common testIDs from both datasets...")
  # extract two arrays of testIDs, and find which are common between them both
  testID_16s  <- pure_data_16s$PIG_DATE
  testID_qpcr <- pure_data_qpcr$PIG_DATE
  
  common_IDs  <- intersect(testID_16s, testID_qpcr)
  
  # only keep relevant testIDs
  pure_data_16s   = pure_data_16s[data_16s$PIG_DATE %in% common_IDs, ]
  pure_data_qpcr  = pure_data_qpcr[data_qpcr$PIG_DATE %in% common_IDs, ]
  
  # sort the datasets based testID, append whilst removing 1 edition of testID and OUA
  pure_data_16s   = pure_data_16s[order(pure_data_16s$PIG_DATE), ]
  pure_data_qpcr  = pure_data_qpcr[order(pure_data_qpcr$PIG_DATE), ]
  pure_data_qpcr  = subset(pure_data_qpcr, select = -c(PIG_DATE, OUA))
  
  # extracting testID and OUA from 16s data
  testID  <- pure_data_16s$PIG_DATE
  OUA     <- pure_data_16s$OUA
  
  pure_data_16s = subset(pure_data_16s, select = -c(PIG_DATE, OUA))
  
  # moving up in taxonomy for the 16s data, going from species to genus
  origNames <- colnames(pure_data_16s)
  newNames <- apply(stringr::str_split_fixed(string = origNames, pattern = "_",8)[,1:6],1, paste, collapse="_")
  
  length(unique(newNames))
  uniqNames=unique(newNames)
  
  newDat=data.frame(dummy=1:NROW(pure_data_16s))
  
  #j=uniqNames[1]
  for(j in uniqNames) {
    
    jIndx=grep(j,origNames )
    if(length(jIndx)>1) {
      newDat=cbind(newDat,rowSums(pure_data_16s[,jIndx]))
    } else {
      newDat=cbind(newDat,(pure_data_16s[,jIndx]))
    }
  }
  
  newDat=newDat[,-1]
  colnames(newDat)=uniqNames
  
  # now making that new data into 16s data set
  pure_data_16s = newDat
  
  # now, append the datasets, to get one set, along with testID and OUA
  complete_data <- cbind(testID, OUA, pure_data_16s, pure_data_qpcr, deparse.level = 1)
  
  # keeping the length of 16s and qpcr data - JUST the results, no ID or OUA
  len_16s   <- length(data_16s[1,]) - 2
  len_qpcr  <- length(data_qpcr[1,])
  
  # dropping the variables, except the resulting data set and the respective length of the datasets
  rm(list=setdiff(ls(), c("complete_data", "len_16s", "len_qpcr")))
  
  print("Data loading finished...")
  #list(complete_data, len_16s, len_qpcr)
  return(complete_data)
}
