context("Natural abundance correction")
library(accucor)

test_that("Carbon correction is accurate", {
  resolution <- 100000
  resolution_defined_at <- 200
  carbon_input_file <- system.file("extdata", "C_Sample_Input_Simple.xlsx",
                                   package = "accucor")

  carbon_corrected <- natural_abundance_correction(
    path = carbon_input_file,
    output_base = FALSE,
    resolution = resolution, resolution_defined_at = resolution_defined_at)

  read_expected <- function(file, sheet) {
    expected <- readxl::read_excel(path = file, sheet = sheet)
    expected <- dplyr::mutate_at(expected,
                                 dplyr::vars(dplyr::ends_with("_Label")),
                                 as.integer)
  }
  carbon_expected <- list(
    "Original" = read_expected(
      system.file("extdata", "C_Sample_Input_Simple.xlsx", package = "accucor"),
      sheet = 1),
    "Corrected" = read_expected(
      system.file("extdata", "C_Sample_Input_Simple_corrected.xlsx", package = "accucor"),
      sheet = "Corrected"),
    "Normalized" = read_expected(
      system.file("extdata", "C_Sample_Input_Simple_corrected.xlsx", package = "accucor"),
      sheet = "Normalized"),
    "PoolAfterDF" = read_expected(
      system.file("extdata", "C_Sample_Input_Simple_corrected.xlsx", package = "accucor"),
      sheet = "PoolAfterDF")
  )

  # Must convert to dataframe due to https://github.com/tidyverse/dplyr/issues/2751
  expect_equal(as.data.frame(carbon_corrected$Original),
               as.data.frame(carbon_expected$Original))
  expect_equal(as.data.frame(carbon_corrected$Corrected),
               as.data.frame(carbon_expected$Corrected))
  expect_equal(as.data.frame(carbon_corrected$Normalized),
               as.data.frame(carbon_expected$Normalized))
  expect_equal(as.data.frame(carbon_corrected$PoolAfterDF),
               as.data.frame(carbon_expected$PoolAfterDF))
})


