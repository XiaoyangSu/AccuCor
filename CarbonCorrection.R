###Please note that this code only deals with carbon labeling experiment.
###Please use other version of correction codes for other experiments.
#Please make sure you have the packages below installed.###
#Use install.packages("") to install packages###

library(accucor)
library(xlsx)

###Please make sure these parameters are accurate.#
C13Purity <- 0.99
Resolution <- 140000 #For Exactive, the Resolution is 100000, defined at Mw 200#
ResDefAt <- 200
ReportPoolSize <- TRUE

#Please make sure the paths below are correct.
#Make sure to use / in the paths
InputFile <- "test/JCGC_test_2.xlsx"
InputSheetName <- "JCGC_test"

OutputFile <- "test/JCGC_test_2_corrected.xlsx"

InputFile <- "inst/extdata/C_Sample_Input_Simple.xlsx"
InputSheetName <- "13C"
OutputDataFrames <- carbon_correction(InputFile, InputSheetName,
                                      Resolution=Resolution, ResDefAt=ResDefAt, ReportPoolSize = ReportPoolSize)


# Write to output file using list of dataframes
# https://github.com/ropensci/writexl/issues/3
# Get original data as read in by readxl (need to modify io function), store as original sheet name
write.xlsx2(OutputDataFrames$Corrected, file=OutputFile, sheetName = "Corrected", row.names=FALSE, append=TRUE)
write.xlsx2(OutputDataFrames$Normalized, file=OutputFile, sheetName = "Normalized", row.names=FALSE, append=TRUE)
write.xlsx2(OutputDataFrames$PoolAfterDF, file=OutputFile, sheetName = "Pool Size", row.names=FALSE, append=TRUE)
