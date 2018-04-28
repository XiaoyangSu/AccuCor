context("Natural abundance correction")
library(accucor)

test_that("Carbon correction (Excel, El-MAVEN export)", {
  resolution <- 100000
  resolution_defined_at <- 200
  input_file <- system.file("extdata", "C_Sample_Input_Simple.xlsx",
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


test_that("Deuterium correction (Excel, El-MAVEN export)", {
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


test_that("Nitrogen correction (Excel, El-MAVEN Export)", {
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
