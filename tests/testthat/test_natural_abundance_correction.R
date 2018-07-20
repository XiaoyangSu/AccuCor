context("Natural abundance correction")
library(accucor)

test_that("Carbon correction (Excel, simple format)", {
  resolution <- 100000
  input_file <- system.file("extdata", "C_Sample_Input_Simple.xlsx",
                                   package = "accucor")

  corrected <- natural_abundance_correction(
    path = input_file,
    output_base = FALSE,
    resolution = resolution)

  read_expected <- function(file, sheet) {
    expected <- readxl::read_excel(path = file, sheet = sheet)
    expected <- dplyr::mutate_at(expected,
                                 dplyr::vars(dplyr::ends_with("_Label")),
                                 as.integer)
  }
  expected_output <- list(
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
  expect_equal(as.data.frame(corrected$Original),
               as.data.frame(expected_output$Original))
  expect_equal(as.data.frame(corrected$Corrected),
               as.data.frame(expected_output$Corrected))
  expect_equal(as.data.frame(corrected$Normalized),
               as.data.frame(expected_output$Normalized))
  expect_equal(as.data.frame(corrected$PoolAfterDF),
               as.data.frame(expected_output$PoolAfterDF))
})


test_that("PoolBeforeDF parameter", {
  resolution <- 100000
  input_file <- system.file("extdata", "C_Sample_Input_Simple.xlsx",
                            package = "accucor")

  corrected <- natural_abundance_correction(
    path = input_file,
    output_base = FALSE,
    report_pool_size_before_df = TRUE,
    resolution = resolution)

  read_expected <- function(file, sheet) {
    expected <- readxl::read_excel(path = file, sheet = sheet)
    expected <- dplyr::mutate_at(expected,
                                 dplyr::vars(dplyr::ends_with("_Label")),
                                 as.integer)
  }
  expected_output <- list(
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
      sheet = "PoolAfterDF"),
    "PoolBeforeDF" = read_expected(
      system.file("extdata", "C_Sample_Input_Simple_corrected.xlsx", package = "accucor"),
      sheet = "PoolBeforeDF")
  )

  # Must convert to dataframe due to https://github.com/tidyverse/dplyr/issues/2751
  expect_equal(as.data.frame(corrected$Original),
               as.data.frame(expected_output$Original))
  expect_equal(as.data.frame(corrected$Corrected),
               as.data.frame(expected_output$Corrected))
  expect_equal(as.data.frame(corrected$Normalized),
               as.data.frame(expected_output$Normalized))
  expect_equal(as.data.frame(corrected$PoolAfterDF),
               as.data.frame(expected_output$PoolAfterDF))
  expect_equal(as.data.frame(corrected$PoolBeforeDF),
               as.data.frame(expected_output$PoolBeforeDF))
})

test_that("Carbon correction (csv, simple format)", {
  resolution <- 100000
  resolution_defined_at <- 200
  input_file <- system.file("extdata", "C_Sample_Input_Simple.csv",
                            package = "accucor")

  corrected <- natural_abundance_correction(
    path = input_file,
    output_base = FALSE,
    resolution = resolution, resolution_defined_at = resolution_defined_at)

  read_expected <- function(file, sheet) {
    expected <- readxl::read_excel(path = file, sheet = sheet)
    expected <- dplyr::mutate_at(expected,
                                 dplyr::vars(dplyr::ends_with("_Label")),
                                 as.integer)
  }
  expected_output <- list(
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
  expect_equal(as.data.frame(corrected$Original),
               as.data.frame(expected_output$Original))
  expect_equal(as.data.frame(corrected$Corrected),
               as.data.frame(expected_output$Corrected))
  expect_equal(as.data.frame(corrected$Normalized),
               as.data.frame(expected_output$Normalized))
  expect_equal(as.data.frame(corrected$PoolAfterDF),
               as.data.frame(expected_output$PoolAfterDF))
})


test_that("Carbon correction (Excel, Classic MAVEN copy/paste)", {
  resolution <- 100000
  resolution_defined_at <- 200
  input_file <- system.file("extdata", "C_Sample_Input.xlsx",
                            package = "accucor")
  knowns_file <- system.file("extdata", "KNOWNS.csv", package = "accucor")
  corrected <- natural_abundance_correction(
    path = input_file,
    compound_database = knowns_file,
    output_base = FALSE,
    resolution = resolution, resolution_defined_at = resolution_defined_at)

  read_expected <- function(file, sheet) {
    expected <- readxl::read_excel(path = file, sheet = sheet)
    expected <- dplyr::mutate_at(expected,
                                 dplyr::vars(dplyr::ends_with("_Label")),
                                 as.integer)
  }
  expected_output <- list(
    "Original" = read_expected(
      system.file("extdata", "C_Sample_Input_Simple_corrected.xlsx", package = "accucor"),
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
  expect_equal(as.data.frame(corrected$Original),
               as.data.frame(expected_output$Original))
  expect_equal(as.data.frame(corrected$Corrected),
               as.data.frame(expected_output$Corrected))
  expect_equal(as.data.frame(corrected$Normalized),
               as.data.frame(expected_output$Normalized))
  expect_equal(as.data.frame(corrected$PoolAfterDF),
               as.data.frame(expected_output$PoolAfterDF))
})


test_that("Deuterium correction (Excel, simple format)", {
  resolution <- 100000
  resolution_defined_at <- 200
  input_file <- system.file("extdata", "D_Sample_Input_Simple.xlsx",
                                   package = "accucor")

  corrected <- natural_abundance_correction(
    path = input_file,
    output_base = FALSE,
    resolution = resolution, resolution_defined_at = resolution_defined_at)

  read_expected <- function(file, sheet) {
    expected <- readxl::read_excel(path = file, sheet = sheet)
    expected <- dplyr::mutate_at(expected,
                                 dplyr::vars(dplyr::ends_with("_Label")),
                                 as.integer)
  }
  expected_output <- list(
    "Original" = read_expected(
      system.file("extdata", "D_Sample_Input_Simple.xlsx", package = "accucor"),
      sheet = 1),
    "Corrected" = read_expected(
      system.file("extdata", "D_Sample_Input_Simple_corrected.xlsx", package = "accucor"),
      sheet = "Corrected"),
    "Normalized" = read_expected(
      system.file("extdata", "D_Sample_Input_Simple_corrected.xlsx", package = "accucor"),
      sheet = "Normalized"),
    "PoolAfterDF" = read_expected(
      system.file("extdata", "D_Sample_Input_Simple_corrected.xlsx", package = "accucor"),
      sheet = "PoolAfterDF")
  )

  # Must convert to dataframe due to https://github.com/tidyverse/dplyr/issues/2751
  expect_equal(as.data.frame(corrected$Original),
               as.data.frame(expected_output$Original))
  expect_equal(as.data.frame(corrected$Corrected),
               as.data.frame(expected_output$Corrected))
  expect_equal(as.data.frame(corrected$Normalized),
               as.data.frame(expected_output$Normalized))
  expect_equal(as.data.frame(corrected$PoolAfterDF),
               as.data.frame(expected_output$PoolAfterDF))
})


test_that("Deuterium correction (Excel, Classic Maven Cut/Paste)", {
  resolution <- 100000
  resolution_defined_at <- 200
  input_file <- system.file("extdata", "D_Sample_Input.xlsx",
                            package = "accucor")
  knowns_file <- system.file("extdata", "KNOWNS.csv", package = "accucor")

  corrected <- natural_abundance_correction(
    path = input_file,
    compound_database = knowns_file,
    output_base = FALSE,
    resolution = resolution, resolution_defined_at = resolution_defined_at)

  read_expected <- function(file, sheet) {
    expected <- readxl::read_excel(path = file, sheet = sheet)
    expected <- dplyr::mutate_at(expected,
                                 dplyr::vars(dplyr::ends_with("_Label")),
                                 as.integer)
  }
  expected_output <- list(
    "Original" = read_expected(
      system.file("extdata", "D_Sample_Input_Simple.xlsx", package = "accucor"),
      sheet = 1),
    "Corrected" = read_expected(
      system.file("extdata", "D_Sample_Input_Simple_corrected.xlsx", package = "accucor"),
      sheet = "Corrected"),
    "Normalized" = read_expected(
      system.file("extdata", "D_Sample_Input_Simple_corrected.xlsx", package = "accucor"),
      sheet = "Normalized"),
    "PoolAfterDF" = read_expected(
      system.file("extdata", "D_Sample_Input_Simple_corrected.xlsx", package = "accucor"),
      sheet = "PoolAfterDF")
  )

  # Must convert to dataframe due to https://github.com/tidyverse/dplyr/issues/2751
  expect_equal(as.data.frame(corrected$Original),
               as.data.frame(expected_output$Original))
  expect_equal(as.data.frame(corrected$Corrected),
               as.data.frame(expected_output$Corrected))
  expect_equal(as.data.frame(corrected$Normalized),
               as.data.frame(expected_output$Normalized))
  expect_equal(as.data.frame(corrected$PoolAfterDF),
               as.data.frame(expected_output$PoolAfterDF))
})


test_that("Nitrogen correction (Excel, simple format)", {
  resolution <- 140000
  resolution_defined_at <- 200
  input_file <- system.file("extdata", "N_Sample_Input_Simple.xlsx",
                            package = "accucor")

  corrected <- natural_abundance_correction(
    path = input_file,
    output_base = FALSE,
    resolution = resolution, resolution_defined_at = resolution_defined_at)

  read_expected <- function(file, sheet) {
    expected <- readxl::read_excel(path = file, sheet = sheet)
    expected <- dplyr::mutate_at(expected,
                                 dplyr::vars(dplyr::ends_with("_Label")),
                                 as.integer)
  }
  expected_output <- list(
    "Original" = read_expected(
      system.file("extdata", "N_Sample_Input_Simple.xlsx", package = "accucor"),
      sheet = 1),
    "Corrected" = read_expected(
      system.file("extdata", "N_Sample_Input_Simple_corrected.xlsx", package = "accucor"),
      sheet = "Corrected"),
    "Normalized" = read_expected(
      system.file("extdata", "N_Sample_Input_Simple_corrected.xlsx", package = "accucor"),
      sheet = "Normalized"),
    "PoolAfterDF" = read_expected(
      system.file("extdata", "N_Sample_Input_Simple_corrected.xlsx", package = "accucor"),
      sheet = "PoolAfterDF")
  )

  # Must convert to dataframe due to https://github.com/tidyverse/dplyr/issues/2751
  expect_equal(as.data.frame(corrected$Original),
               as.data.frame(expected_output$Original))
  expect_equal(as.data.frame(corrected$Corrected),
               as.data.frame(expected_output$Corrected))
  expect_equal(as.data.frame(corrected$Normalized),
               as.data.frame(expected_output$Normalized))
  expect_equal(as.data.frame(corrected$PoolAfterDF),
               as.data.frame(expected_output$PoolAfterDF))
})


test_that("Nitrogen correction (Excel, Classic Maven Cut/Paste)", {
  resolution <- 140000
  resolution_defined_at <- 200
  input_file <- system.file("extdata", "N_Sample_Input.xlsx",
                            package = "accucor")
  knowns_file <- system.file("extdata", "KNOWNS.csv", package = "accucor")

  corrected <- natural_abundance_correction(
    path = input_file,
    compound_database = knowns_file,
    output_base = FALSE,
    resolution = resolution, resolution_defined_at = resolution_defined_at)

  read_expected <- function(file, sheet) {
    expected <- readxl::read_excel(path = file, sheet = sheet)
    expected <- dplyr::mutate_at(expected,
                                 dplyr::vars(dplyr::ends_with("_Label")),
                                 as.integer)
  }
  expected_output <- list(
    "Original" = read_expected(
      system.file("extdata", "N_Sample_Input_Simple.xlsx", package = "accucor"),
      sheet = 1),
    "Corrected" = read_expected(
      system.file("extdata", "N_Sample_Input_Simple_corrected.xlsx", package = "accucor"),
      sheet = "Corrected"),
    "Normalized" = read_expected(
      system.file("extdata", "N_Sample_Input_Simple_corrected.xlsx", package = "accucor"),
      sheet = "Normalized"),
    "PoolAfterDF" = read_expected(
      system.file("extdata", "N_Sample_Input_Simple_corrected.xlsx", package = "accucor"),
      sheet = "PoolAfterDF")
  )

  # Must convert to dataframe due to https://github.com/tidyverse/dplyr/issues/2751
  expect_equal(as.data.frame(corrected$Original),
               as.data.frame(expected_output$Original))
  expect_equal(as.data.frame(corrected$Corrected),
               as.data.frame(expected_output$Corrected))
  expect_equal(as.data.frame(corrected$Normalized),
               as.data.frame(expected_output$Normalized))
  expect_equal(as.data.frame(corrected$PoolAfterDF),
               as.data.frame(expected_output$PoolAfterDF))
})


test_that("Carbon correction (csv, El-MAVEN export (with set names))", {
  resolution <- 140000
  resolution_defined_at <- 200
  input_file <- system.file("extdata", "elmaven_export.csv",
                            package = "accucor")

  # This file includes the empty line from exported "set names" which is also
  # incorrectly formatted (one too few columns)
  # Warnings from parsing this input are expected
  corrected <- expect_warning(natural_abundance_correction(
    path = input_file,
    output_base = FALSE,
    resolution = resolution, resolution_defined_at = resolution_defined_at))

  read_expected <- function(file, sheet) {
    expected <- readxl::read_excel(path = file, sheet = sheet)
    expected <- dplyr::mutate_at(expected,
                                 dplyr::vars(dplyr::ends_with("_Label")),
                                 as.integer)
  }
  expected_output <- list(
    "Original" = read_expected(
      system.file("extdata", "elmaven_export_corrected.xlsx", package = "accucor"),
      sheet = 1),
    "Corrected" = read_expected(
      system.file("extdata", "elmaven_export_corrected.xlsx", package = "accucor"),
      sheet = "Corrected"),
    "Normalized" = read_expected(
      system.file("extdata", "elmaven_export_corrected.xlsx", package = "accucor"),
      sheet = "Normalized"),
    "PoolAfterDF" = read_expected(
      system.file("extdata", "elmaven_export_corrected.xlsx", package = "accucor"),
      sheet = "PoolAfterDF")
  )

  # Must convert to dataframe due to https://github.com/tidyverse/dplyr/issues/2751
  # Blank label column parses as different type
  expected_output$Original$label <-
    sapply(expected_output$Original$label, as.character)
  expect_equal(as.data.frame(corrected$Original),
               as.data.frame(expected_output$Original))
  expect_equal(as.data.frame(corrected$Corrected),
               as.data.frame(expected_output$Corrected))
  expect_equal(as.data.frame(corrected$Normalized),
               as.data.frame(expected_output$Normalized))
  expect_equal(as.data.frame(corrected$PoolAfterDF),
               as.data.frame(expected_output$PoolAfterDF))
})


test_that("Carbon correction (Excel, El-MAVEN export (with set names))", {
  resolution <- 140000
  resolution_defined_at <- 200
  input_file <- system.file("extdata", "elmaven_export.xlsx",
                            package = "accucor")

  corrected <- natural_abundance_correction(
    path = input_file,
    output_base = FALSE,
    resolution = resolution, resolution_defined_at = resolution_defined_at)

  read_expected <- function(file, sheet) {
    expected <- readxl::read_excel(path = file, sheet = sheet)
    expected <- dplyr::mutate_at(expected,
                                 dplyr::vars(dplyr::ends_with("_Label")),
                                 as.integer)
  }
  expected_output <- list(
    "Original" = read_expected(
      system.file("extdata", "elmaven_export_corrected.xlsx", package = "accucor"),
      sheet = 1),
    "Corrected" = read_expected(
      system.file("extdata", "elmaven_export_corrected.xlsx", package = "accucor"),
      sheet = "Corrected"),
    "Normalized" = read_expected(
      system.file("extdata", "elmaven_export_corrected.xlsx", package = "accucor"),
      sheet = "Normalized"),
    "PoolAfterDF" = read_expected(
      system.file("extdata", "elmaven_export_corrected.xlsx", package = "accucor"),
      sheet = "PoolAfterDF")
  )

  # Must convert to dataframe due to https://github.com/tidyverse/dplyr/issues/2751
  expect_equal(as.data.frame(corrected$Original),
               as.data.frame(expected_output$Original))
  expect_equal(as.data.frame(corrected$Corrected),
               as.data.frame(expected_output$Corrected))
  expect_equal(as.data.frame(corrected$Normalized),
               as.data.frame(expected_output$Normalized))
  expect_equal(as.data.frame(corrected$PoolAfterDF),
               as.data.frame(expected_output$PoolAfterDF))
})


test_that("Carbon correction (csv, El-MAVEN export (w/o names))", {
  resolution <- 140000
  input_file <- system.file("extdata", "elmaven_d2_export.csv",
                            package = "accucor")

  corrected <- natural_abundance_correction(
    path = input_file,
    resolution = resolution,
    output_base = FALSE)

  read_expected <- function(file, sheet) {
    expected <- readxl::read_excel(path = file, sheet = sheet)
    expected <- dplyr::mutate_at(expected,
                                 dplyr::vars(dplyr::ends_with("_Label")),
                                 as.integer)
  }
  expected_output <- list(
    "Original" = read_expected(
      system.file("extdata", "elmaven_d2_export_corrected.xlsx", package = "accucor"),
      sheet = 1),
    "Corrected" = read_expected(
      system.file("extdata", "elmaven_d2_export_corrected.xlsx", package = "accucor"),
      sheet = "Corrected"),
    "Normalized" = read_expected(
      system.file("extdata", "elmaven_d2_export_corrected.xlsx", package = "accucor"),
      sheet = "Normalized"),
    "PoolAfterDF" = read_expected(
      system.file("extdata", "elmaven_d2_export_corrected.xlsx", package = "accucor"),
      sheet = "PoolAfterDF")
  )

  # Must convert to dataframe due to https://github.com/tidyverse/dplyr/issues/2751
  # Label column parses as different type
  expected_output$Original$label <-
    sapply(expected_output$Original$label, as.character)
  expect_equal(as.data.frame(corrected$Original),
               as.data.frame(expected_output$Original))
  expect_equal(as.data.frame(corrected$Corrected),
               as.data.frame(expected_output$Corrected))
  expect_equal(as.data.frame(corrected$Normalized),
               as.data.frame(expected_output$Normalized))
  expect_equal(as.data.frame(corrected$PoolAfterDF),
               as.data.frame(expected_output$PoolAfterDF))
})


test_that("Carbon correction (csv, El-MAVEN export, multiple groups per compound)", {
  resolution <- 140000
  input_file <- system.file("extdata", "alanine_three_peak_groups.csv",
                            package = "accucor")

  corrected <- natural_abundance_correction(
    path = input_file,
    resolution = resolution,
    output_base = FALSE)

  read_expected <- function(file, sheet) {
    expected <- readxl::read_excel(path = file, sheet = sheet)
    expected <- dplyr::mutate_at(expected,
                                 dplyr::vars(dplyr::ends_with("_Label")),
                                 as.integer)
  }
  expected_output <- list(
    "Original" = read_expected(
      system.file("extdata", "alanine_three_peak_groups_corrected.xlsx", package = "accucor"),
      sheet = 1),
    "Corrected" = read_expected(
      system.file("extdata", "alanine_three_peak_groups_corrected.xlsx", package = "accucor"),
      sheet = "Corrected"),
    "Normalized" = read_expected(
      system.file("extdata", "alanine_three_peak_groups_corrected.xlsx", package = "accucor"),
      sheet = "Normalized"),
    "PoolAfterDF" = read_expected(
      system.file("extdata", "alanine_three_peak_groups_corrected.xlsx", package = "accucor"),
      sheet = "PoolAfterDF")
  )

  # Must convert to dataframe due to https://github.com/tidyverse/dplyr/issues/2751
  # Label column parses as different type
  expected_output$Original$label <-
    sapply(expected_output$Original$label, as.character)
  expect_equal(as.data.frame(corrected$Original),
               as.data.frame(expected_output$Original))
  expect_equal(as.data.frame(corrected$Corrected),
               as.data.frame(expected_output$Corrected))
  expect_equal(as.data.frame(corrected$Normalized),
               as.data.frame(expected_output$Normalized))
  expect_equal(as.data.frame(corrected$PoolAfterDF),
               as.data.frame(expected_output$PoolAfterDF))
})
