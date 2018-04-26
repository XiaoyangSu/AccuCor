###Please note that this code only deals with carbon labeling experiment.
###Please use other version of correction codes for other experiments.

# Use install.packages("accucor") to install accucor
library(accucor)

# Please make sure these parameters are accurate.
Resolution <- 140000  # For Exactive, the Resolution is 100000, defined at Mw 200#
ResDefAt <- 200
carbon_input_file <- system.file("extdata", "C_Sample_Input_Simple.xlsx", package = "accucor")

# Output is written to InputFile_corrected.xlsx
carbon_corrected <- natural_abundance_correction(path = carbon_input_file, Resolution = Resolution, ResDefAt = ResDefAt)


#####
# Optional parameters
# Default is first sheet, specify name or number
jcgc_input_file < system.file("extdata", "JCGC_test_2.xlsx", package = "accucor")
sheet <- "13C"
# Set output_path to override location, set output_path to FALSE to avoid writing
jcgc_output_file <- "JCGC_test_2_corrected.xlsx"
# C13Purity defaults to 0.99
c13_purity <- 0.99
# Default is to ReportPoolSize (NOT USED?)
ReportPoolSize <- TRUE

jcgc_corrected <- natural_abundance_correction(jcgc_input_file, sheet = sheet, output_path = jcgc_output_file,
                                                 Resolution = Resolution, ResDefAt = ResDefAt,
                                                 purity = c13_purity, ReportPoolSize = ReportPoolSize)

# Test Deuterium
deuterium_input_file = system.file("extdata", "D_Sample_Input_Simple.xlsx", package = "accucor")
deuterium_corrected <- natural_abundance_correction(path = deuterium_input_file,
                                                    Resolution = Resolution, ResDefAt = ResDefAt)

# Test Nitrogen
nitrogen_input_file = system.file("extdata", "N_Sample_Input_Simple.xlsx", package = "accucor")
nitrogen_corrected <- natural_abundance_correction(path = nitrogen_input_file,
                                                   Resolution = Resolution, ResDefAt = ResDefAt)
