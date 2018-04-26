###Please note that this code only deals with carbon labeling experiment.
###Please use other version of correction codes for other experiments.

# Use install.packages("accucor") to install accucor
library(accucor)

# Please make sure these parameters are accurate.
resolution <- 140000  # For Exactive, the Resolution is 100000, defined at Mw 200#
resolution_defined_at <- 200
carbon_input_file <- system.file("extdata", "C_Sample_Input_Simple.xlsx", package = "accucor")

# Output is written to InputFile_corrected.xlsx
carbon_corrected <- natural_abundance_correction(
  path = carbon_input_file,
  resolution = resolution, resolution_defined_at = resolution_defined_at)


#####
# Optional parameters
# Default is first sheet, specify name or number
jcgc_input_file <- system.file("extdata", "JCGC_test_2.xlsx", package = "accucor")
sheet <- "JCGC_test"
# Set output_path to override location, set output_path to FALSE to avoid writing
jcgc_output_base <- "JCGC_test_2_output.xlsx"
# C13Purity defaults to 0.99
c13_purity <- 0.99
# Default is to ReportPoolSize (NOT USED?)
report_pool_size <- TRUE

jcgc_corrected <- natural_abundance_correction(
  jcgc_input_file, sheet = sheet, output_base = jcgc_output_base,
  resolution = resolution, resolution_defined_at = resolution_defined_at,
  purity = c13_purity, report_pool_size = report_pool_size)

# Test Deuterium
deuterium_input_file <- system.file(
  "extdata", "D_Sample_Input_Simple.xlsx", package = "accucor")
deuterium_corrected <- natural_abundance_correction(
  path = deuterium_input_file, resolution = resolution,
  resolution_defined_at = resolution_defined_at)

# Test Nitrogen
nitrogen_input_file <- system.file(
  "extdata", "N_Sample_Input_Simple.xlsx", package = "accucor")
nitrogen_corrected <- natural_abundance_correction(
  path = nitrogen_input_file,
  resolution = resolution, resolution_defined_at = resolution_defined_at)
