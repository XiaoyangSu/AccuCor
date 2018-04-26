###Please note that this code only deals with carbon labeling experiment.
###Please use other version of correction codes for other experiments.

# Use install.packages("accucor") to install accucor
library(accucor)

# Please make sure these parameters are accurate.
Resolution <- 140000  # For Exactive, the Resolution is 100000, defined at Mw 200#
ResDefAt <- 200
InputFile <- "inst/extdata/C_Sample_Input_Simple.xlsx"

# Output is written to InputFile_corrected.xlsx
OutputDataFrames <- natural_abundance_correction(path = InputFile, Resolution = Resolution, ResDefAt = ResDefAt)


#####
# Optional parameters
# Default is first sheet, specify name or number
InputSheetName <- "13C"
# Set output_path to override location, set output_path to FALSE to avoid writing
OutputFile <- "JCGC_test_2_corrected.xlsx"
# C13Purity defaults to 0.99
C13Purity <- 0.99
# Default is to ReportPoolSize (NOT USED?)
ReportPoolSize <- TRUE

OutputDataFrames <- carbon_correction(InputFile, sheet = InputSheetName, output_path = OutputFile,
                                      Resolution = Resolution, ResDefAt = ResDefAt,
                                      C13Purity = C13Purity, ReportPoolSize = ReportPoolSize)

# Test Deuterium
deuterium_input_file = "inst/extdata/D_Sample_Input_Simple.xlsx"
deuterium_corrected <- natural_abundance_correction(path = deuterium_input_file,
                                                    Resolution = Resolution, ResDefAt = ResDefAt)

# Test Nitrogen
nitrogen_input_file = "inst/extdata/N_Sample_Input_Simple.xlsx"
nitrogen_corrected <- natural_abundance_correction(path = nitrogen_input_file,
                                                    Resolution = Resolution, ResDefAt = ResDefAt)
