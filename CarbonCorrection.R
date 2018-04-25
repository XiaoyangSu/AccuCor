###Please note that this code only deals with carbon labeling experiment.
###Please use other version of correction codes for other experiments.
#Please make sure you have the packages below installed.###
#Use install.packages("") to install packages###

library(accucor)

###Please make sure these parameters are accurate.#
C13Purity <- 0.99
Resolution <- 140000 #For Exactive, the Resolution is 100000, defined at Mw 200#
ResDefAt <- 200
ReportPoolSize <- TRUE

#Please make sure the paths below are correct.
#Make sure to use / in the paths
InputFile <- "test/JCGC_test_2.xlsx"
InputSheetName <- "JCGC_test"

OutputFile <- "JCGC_test_2_corrected.xlsx"

InputFile <- "inst/extdata/C_Sample_Input_Simple.xlsx"
InputSheetName <- "13C"
OutputDataFrames <- carbon_correction(InputFile, InputSheetName, output_path = OutputFile,
                                      Resolution = Resolution, ResDefAt = ResDefAt,
                                      C13Purity = C13Purity, ReportPoolSize = ReportPoolSize)
# Output is written to InputFile_corrected.xlsx
# Set output_path to override location, set output_path to FALSE to avoid writing
