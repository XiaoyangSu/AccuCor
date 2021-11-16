#' Natural Abundance carbon isotope correction for one metabolite
#'
#' @param formula String representing molecular formula
#' @param datamatrix Matrix of abundances for each sample for each isotope
#' @param label vector of integer labels
#' @param Resolution For Exactive, the Resolution is 100000, defined at Mw 200
#' @param ResDefAt Resolution defined at (in Mw), e.g. 200 Mw
#' @param purity Carbon 13 purity, default: 0.99
#' @param ReportPoolSize default: TRUE
#' @importFrom rlang .data
#' @importFrom dplyr "%>%"
#' @return Named list of matrices: 'Corrected', 'Normalized',
#'      'PoolBeforeDF', and 'PoolAfterDF'.
#' @export
#' @examples
#' \dontrun{
#' carbon_isotope_correction(formula = "C6H13O9P",
#'                           datamatrix = DataMatrix,
#'                           label = c(0, 1, 2, 3, 4, 5),
#'                           Resolution = 100000)
#' }
carbon_isotope_correction <- function(formula,
                                      datamatrix,
                                      label,
                                      Resolution,
                                      ResDefAt = 200,
                                      purity = 0.99,
                                      ReportPoolSize = TRUE) {
  CarbonNaturalAbundace <- c(0.9893, 0.0107)
  HydrogenNaturalAbundace <- c(0.999885, 0.000115)
  NitrogenNaturalAbundace <- c(0.99636, 0.00364)
  OxygenNaturalAbundace <- c(0.99757, 0.00038, 0.00205)
  SulfurNaturalAbundace <- c(0.95, 0.0075, 0.0425)
  SiliconNaturalAbundace <- c(0.92223, 0.04685, 0.03092)
  ChlorineNaturalAbundance <- c(0.7576, 0.2424)
  BromineNaturalAbundance <- c(0.5069, 0.4931)

  AtomNumber <- rep(0, 9)
  names(AtomNumber) <- c("C", "H", "N", "O", "P", "S", "Si", "Cl", "Br")
  # Required to ensure the "thermo" object is created and defaults are used
  suppressMessages(CHNOSZ::reset())
  AtomicComposition <- CHNOSZ::makeup(formula)

  for (i in names(AtomicComposition)) {
    AtomNumber[i] <- AtomicComposition[i]
  }

  MolecularWeight <- sum(AtomNumber * c(12, 1, 14, 16, 31, 32, 28, 35.5, 80))
  Mass.Limit <- 1.66 * MolecularWeight^1.5 / Resolution / sqrt(ResDefAt)

  ExpMatrix <- matrix(0, ncol = ncol(datamatrix), nrow = AtomNumber["C"] + 1)
  CorrectedMatrix <- matrix(
    0,
    ncol = ncol(datamatrix),
    nrow = AtomNumber["C"] + 1
  )
  if (AtomNumber["C"] < max(label)) {
    stop(
      paste(
        formula,
        ":the number of labeling exceeded the number of carbon atoms,",
        "check the input file."
      )
    )
  }
  else {
    for (i in seq_along(label)) {
      ExpMatrix[label[i] + 1, ] <- datamatrix[i, ]
    }
  }


  PurityMatrix <- diag(AtomNumber["C"] + 1)
  CarbonMatrix <- diag(AtomNumber["C"] + 1)
  NonTracerMatrix <- matrix(
    0,
    ncol = (AtomNumber["C"] + 1),
    nrow = (AtomNumber["C"] + 1)
  )

  for (i in 1:(AtomNumber["C"] + 1)) {
    PurityMatrix[i, ] <- sapply(
      0:(AtomNumber["C"]),
      function(x) {
        stats::dbinom(
          x - i + 1,
          x,
          (1 - purity)
        )
      }
    )
  }

  for (i in 1:(AtomNumber["C"] + 1)) {
    CarbonMatrix[, i] <- sapply(
      0:AtomNumber["C"],
      function(x) {
        stats::dbinom(
          x - i + 1,
          AtomNumber["C"] - i + 1,
          CarbonNaturalAbundace[2]
        )
      }
    )
  }

  Isotope.Combinations <- expand.grid(
    H2 = c(0:AtomNumber["H"]),
    N15 = c(0:AtomNumber["N"]),
    O17 = c(0:AtomNumber["O"]),
    O18 = c(0:AtomNumber["O"]),
    S33 = c(0:AtomNumber["S"]),
    S34 = c(0:AtomNumber["S"]),
    Si29 = c(0:AtomNumber["Si"]),
    Si30 = c(0:AtomNumber["Si"]),
    Cl37 = c(0:AtomNumber["Cl"]),
    Br81 = c(0:AtomNumber["Br"])
  )

  Isotope.Combinations <- Isotope.Combinations %>%
    dplyr::mutate(
      MassSum = .data$H2 +
        .data$N15 +
        .data$O17 +
        .data$O18 * 2 +
        .data$S33 +
        .data$S34 * 2 +
        .data$Si29 +
        .data$Si30 * 2 +
        .data$Cl37 * 2 +
        .data$Br81 * 2
    ) %>%
    dplyr::filter(
      (.data$O17 + .data$O18) <= AtomNumber["O"] &
        (.data$S33 + .data$S34) <= AtomNumber["S"] &
        (.data$Si29 + .data$Si30) <= AtomNumber["Si"] &
        .data$MassSum <= AtomNumber["C"]
    ) %>%
    dplyr::mutate(
      MassDiff = (2.0141 - 1.00783) * .data$H2 +
        (15.00011 - 14.00307) * .data$N15 +
        (16.99913 - 15.99491) * .data$O17 +
        (17.99916 - 15.99491) * .data$O18 +
        (32.97146 - 31.97207) * .data$S33 +
        (33.96787 - 31.97207) * .data$S34 +
        (28.97649 - 27.97693) * .data$Si29 +
        (29.97377 - 27.97693) * .data$Si30 +
        (36.96590 - 34.96885) * .data$Cl37 +
        (80.91629 - 78.91833) * .data$Br81 -
        (13.00335 - 12) * .data$MassSum
    ) %>%
    dplyr::filter(abs(.data$MassDiff) < Mass.Limit)

  for (i in seq_len(nrow(Isotope.Combinations))) {
    p <- stats::dbinom(
      Isotope.Combinations[i, 1],
      AtomNumber["H"],
      HydrogenNaturalAbundace[2]
    ) *
      stats::dbinom(
        Isotope.Combinations[i, 2],
        AtomNumber["N"],
        NitrogenNaturalAbundace[2]
      ) *
      stats::dmultinom(
        unlist(c(
          AtomNumber["O"] - sum(Isotope.Combinations[i, 3:4]),
          Isotope.Combinations[i, 3:4]
        )),
        AtomNumber["O"],
        OxygenNaturalAbundace
      ) *
      stats::dmultinom(
        unlist(c(
          AtomNumber["S"] - sum(Isotope.Combinations[i, 5:6]),
          Isotope.Combinations[i, 5:6]
        )),
        AtomNumber["S"],
        SulfurNaturalAbundace
      ) *
      stats::dmultinom(
        unlist(c(
          AtomNumber["Si"] - sum(Isotope.Combinations[i, 7:8]),
          Isotope.Combinations[i, 7:8]
        )),
        AtomNumber["Si"],
        SiliconNaturalAbundace
      ) *
      stats::dbinom(
        Isotope.Combinations[i, 9],
        AtomNumber["Cl"],
        ChlorineNaturalAbundance[2]
      ) *
      stats::dbinom(
        Isotope.Combinations[i, 10],
        AtomNumber["Br"],
        BromineNaturalAbundance[2]
      )

    MSum <- Isotope.Combinations[i, 11]
    for (j in 1:(AtomNumber["C"] + 1 - MSum)) {
      NonTracerMatrix[MSum + j, j] <- NonTracerMatrix[MSum + j, j] + p
    }
  }

  for (i in seq_len(ncol(datamatrix))) {
    CorrectedMatrix[, i] <- stats::coef(nnls::nnls(
      NonTracerMatrix %*% CarbonMatrix %*% PurityMatrix, ExpMatrix[, i]
    ))
  }

  return(CorrectedMatrix)
}

#' Natural Abundance deuterium isotope correction for one metabolite
#'
#' @param formula String representing molecular formula
#' @param datamatrix Matrix of abundances for each sample for each isotope
#' @param label vector of integer labels
#' @param Resolution For Exactive, the Resolution is 100000, defined at Mw 200
#' @param ResDefAt Resolution defined at (in Mw), e.g. 200 Mw
#' @param purity Deuterium purity, default: 0.99
#' @param ReportPoolSize default: TRUE
#' @importFrom rlang .data
#' @return Named list of matrices: 'Corrected', 'Normalized',
#'      'PoolBeforeDF', and 'PoolAfterDF'.
#' @export
#' @examples
#' \dontrun{
#' deuterium_isotope_correction(formula = "C6H13O9P",
#'                              datamatrix = DataMatrix,
#'                              label = c(0, 1),
#'                              Resolution = 100000)
#' }
deuterium_isotope_correction <- function(formula,
                                         datamatrix,
                                         label,
                                         Resolution,
                                         ResDefAt = 200,
                                         purity = 0.99,
                                         ReportPoolSize = TRUE) {
  CarbonNaturalAbundace <- c(0.9893, 0.0107)
  HydrogenNaturalAbundace <- c(0.999885, 0.000115)
  NitrogenNaturalAbundace <- c(0.99636, 0.00364)
  OxygenNaturalAbundace <- c(0.99757, 0.00038, 0.00205)
  SulfurNaturalAbundace <- c(0.95, 0.0075, 0.0425)
  SiliconNaturalAbundace <- c(0.92223, 0.04685, 0.03092)
  ChlorineNaturalAbundance <- c(0.7576, 0.2424)
  BromineNaturalAbundance <- c(0.5069, 0.4931)

  AtomNumber <- rep(0, 9)
  names(AtomNumber) <- c("C", "H", "N", "O", "P", "S", "Si", "Cl", "Br")
  # Required to ensure the "thermo" object is created and defaults are used
  suppressMessages(CHNOSZ::reset())
  AtomicComposition <- CHNOSZ::makeup(formula)

  for (i in names(AtomicComposition)) {
    AtomNumber[i] <- AtomicComposition[i]
  }

  MolecularWeight <- sum(AtomNumber * c(12, 1, 14, 16, 31, 32, 28, 35.5, 80))
  Mass.Limit <- 1.66 * MolecularWeight^1.5 / Resolution / sqrt(ResDefAt)

  ExpMatrix <- matrix(0, ncol = ncol(datamatrix), nrow = AtomNumber["H"] + 1)
  CorrectedMatrix <- matrix(
    0,
    ncol = ncol(datamatrix),
    nrow = AtomNumber["H"] + 1
  )
  if (AtomNumber["H"] < max(label)) {
    stop(
      paste(
        formula,
        ":the number of labeling exceeded the number of hydrogen atoms,",
        "check the input file."
      )
    )
  }
  else {
    for (i in seq_along(label)) {
      ExpMatrix[label[i] + 1, ] <- datamatrix[i, ]
    }
  }


  PurityMatrix <- diag(AtomNumber["H"] + 1)
  HydrogenMatrix <- diag(AtomNumber["H"] + 1)

  for (i in 1:(AtomNumber["H"] + 1)) {
    PurityMatrix[i, ] <- sapply(
      0:AtomNumber["H"],
      function(x) stats::dbinom(x - i + 1, x, (1 - purity))
    )
  }

  for (i in 1:(AtomNumber["H"] + 1)) {
    HydrogenMatrix[, i] <- sapply(
      0:AtomNumber["H"],
      function(x) {
        stats::dbinom(
          x - i + 1,
          AtomNumber["H"] - i + 1,
          HydrogenNaturalAbundace[2]
        )
      }
    )
  }

  NonTracerMatrix <- matrix(
    0,
    ncol = (AtomNumber["H"] + 1),
    nrow = (AtomNumber["H"] + 1)
  )

  Isotope.Combinations <- expand.grid(
    C13 = c(0:AtomNumber["C"]),
    N15 = c(0:AtomNumber["N"]),
    O17 = c(0:AtomNumber["O"]),
    O18 = c(0:AtomNumber["O"]),
    S33 = c(0:AtomNumber["S"]),
    S34 = c(0:AtomNumber["S"]),
    Si29 = c(0:AtomNumber["Si"]),
    Si30 = c(0:AtomNumber["Si"]),
    Cl37 = c(0:AtomNumber["Cl"]),
    Br81 = c(0:AtomNumber["Br"])
  )

  Isotope.Combinations <- Isotope.Combinations %>%
    dplyr::mutate(
      MassSum = .data$C13 +
        .data$N15 +
        .data$O17 +
        .data$O18 * 2 +
        .data$S33 +
        .data$S34 * 2 +
        .data$Si29 +
        .data$Si30 * 2 +
        .data$Cl37 * 2 +
        .data$Br81 * 2
    ) %>%
    dplyr::filter(
      (.data$O17 + .data$O18) <= AtomNumber["O"] &
        (.data$S33 + .data$S34) <= AtomNumber["S"] &
        (.data$Si29 + .data$Si30) <= AtomNumber["Si"] &
        .data$MassSum <= AtomNumber["H"]
    ) %>%
    dplyr::mutate(
      MassDiff =
        (13.00335 - 12) * .data$C13 +
          (15.00011 - 14.00307) * .data$N15 +
          (16.99913 - 15.99491) * .data$O17 +
          (17.99916 - 15.99491) * .data$O18 +
          (32.97146 - 31.97207) * .data$S33 +
          (33.96787 - 31.97207) * .data$S34 +
          (28.97649 - 27.97693) * .data$Si29 +
          (29.97377 - 27.97693) * .data$Si30 +
          (36.96590 - 34.96885) * .data$Cl37 +
          (80.91629 - 78.91833) * .data$Br81 -
          (2.0141 - 1.00783) * .data$MassSum
    ) %>%
    dplyr::filter(abs(.data$MassDiff) < Mass.Limit)

  for (i in seq_len(nrow(Isotope.Combinations))) {
    p <- stats::dbinom(
      Isotope.Combinations[i, 1],
      AtomNumber["C"],
      CarbonNaturalAbundace[2]
    ) *
      stats::dbinom(
        Isotope.Combinations[i, 2],
        AtomNumber["N"],
        NitrogenNaturalAbundace[2]
      ) *
      stats::dmultinom(
        unlist(c(
          AtomNumber["O"] - sum(Isotope.Combinations[i, 3:4]),
          Isotope.Combinations[i, 3:4]
        )),
        AtomNumber["O"],
        OxygenNaturalAbundace
      ) *
      stats::dmultinom(
        unlist(c(
          AtomNumber["S"] - sum(Isotope.Combinations[i, 5:6]),
          Isotope.Combinations[i, 5:6]
        )),
        AtomNumber["S"],
        SulfurNaturalAbundace
      ) *
      stats::dmultinom(
        unlist(c(
          AtomNumber["Si"] - sum(Isotope.Combinations[i, 7:8]),
          Isotope.Combinations[i, 7:8]
        )),
        AtomNumber["Si"],
        SiliconNaturalAbundace
      ) *
      stats::dbinom(
        Isotope.Combinations[i, 9],
        AtomNumber["Cl"],
        ChlorineNaturalAbundance[2]
      ) *
      stats::dbinom(
        Isotope.Combinations[i, 10],
        AtomNumber["Br"],
        BromineNaturalAbundance[2]
      )

    MSum <- Isotope.Combinations[i, 11]
    for (j in 1:(AtomNumber["H"] + 1 - MSum)) {
      NonTracerMatrix[MSum + j, j] <- NonTracerMatrix[MSum + j, j] + p
    }
  }

  for (i in seq_len(ncol(datamatrix))) {
    CorrectedMatrix[, i] <- stats::coef(nnls::nnls(
      NonTracerMatrix %*% HydrogenMatrix %*% PurityMatrix,
      ExpMatrix[, i]
    ))
  }

  return(CorrectedMatrix)
}


#' Natural Abundance deuterium isotope correction for one metabolite
#'
#' @param formula String representing molecular formula
#' @param datamatrix Matrix of abundances for each sample for each isotope
#' @param label vector of integer labels
#' @param Resolution For Exactive, the Resolution is 100000, defined at Mw 200
#' @param ResDefAt Resolution defined at (in Mw), e.g. 200 Mw
#' @param purity Nitrogen purity, default: 0.99
#' @param ReportPoolSize default: TRUE
#' @importFrom rlang .data
#' @return Named list of matrices: 'Corrected', 'Normalized',
#'      'PoolBeforeDF', and 'PoolAfterDF'.
#' @export
#' @examples
#' \dontrun{
#' nitrogen_isotope_correction(formula = "C23H38N7O17P3S",
#'                              datamatrix = DataMatrix,
#'                              label = c(0, 1, 2, 3, 4, 5, 6, 7),
#'                              Resolution = 140000)
#' }
nitrogen_isotope_correction <- function(formula,
                                        datamatrix,
                                        label,
                                        Resolution,
                                        ResDefAt = 200,
                                        purity = 0.99,
                                        ReportPoolSize = TRUE) {
  CarbonNaturalAbundace <- c(0.9893, 0.0107)
  HydrogenNaturalAbundace <- c(0.999885, 0.000115)
  NitrogenNaturalAbundace <- c(0.99636, 0.00364)
  OxygenNaturalAbundace <- c(0.99757, 0.00038, 0.00205)
  SulfurNaturalAbundace <- c(0.95, 0.0075, 0.0425)
  SiliconNaturalAbundace <- c(0.92223, 0.04685, 0.03092)
  ChlorineNaturalAbundance <- c(0.7576, 0.2424)
  BromineNaturalAbundance <- c(0.5069, 0.4931)

  AtomNumber <- rep(0, 9)
  names(AtomNumber) <- c("C", "H", "N", "O", "P", "S", "Si", "Cl", "Br")
  # Required to ensure the "thermo" object is created and defaults are used
  suppressMessages(CHNOSZ::reset())
  AtomicComposition <- CHNOSZ::makeup(formula)

  for (i in names(AtomicComposition)) {
    AtomNumber[i] <- AtomicComposition[i]
  }

  MolecularWeight <- sum(AtomNumber * c(12, 1, 14, 16, 31, 32, 28, 35.5, 80))
  Mass.Limit <- 1.66 * MolecularWeight^1.5 / Resolution / sqrt(ResDefAt)

  ExpMatrix <- matrix(0, ncol = ncol(datamatrix), nrow = AtomNumber["N"] + 1)
  CorrectedMatrix <- matrix(
    0,
    ncol = ncol(datamatrix),
    nrow = AtomNumber["N"] + 1
  )
  if (AtomNumber["N"] < max(label)) {
    stop(
      paste(
        formula,
        ":the number of labeling exceeded the number of nitrogen atoms,",
        "check the input file."
      )
    )
  }
  else {
    for (i in seq_along(label)) {
      ExpMatrix[label[i] + 1, ] <- datamatrix[i, ]
    }
  }


  PurityMatrix <- diag(AtomNumber["N"] + 1)
  NitrogenMatrix <- diag(AtomNumber["N"] + 1)


  for (i in 1:(AtomNumber["N"] + 1)) {
    PurityMatrix[i, ] <- sapply(
      0:(AtomNumber["N"]),
      function(x) stats::dbinom(x - i + 1, x, (1 - purity))
    )
  }

  for (i in 1:(AtomNumber["N"] + 1)) {
    NitrogenMatrix[, i] <- sapply(
      0:AtomNumber["N"],
      function(x) {
        stats::dbinom(
          x - i + 1,
          AtomNumber["N"] - i + 1,
          NitrogenNaturalAbundace[2]
        )
      }
    )
  }

  NonTracerMatrix <- matrix(
    0,
    ncol = (AtomNumber["N"] + 1),
    nrow = (AtomNumber["N"] + 1)
  )

  Isotope.Combinations <- expand.grid(
    C13 = c(0:AtomNumber["C"]),
    H2 = c(0:AtomNumber["H"]),
    O17 = c(0:AtomNumber["O"]),
    O18 = c(0:AtomNumber["O"]),
    S33 = c(0:AtomNumber["S"]),
    S34 = c(0:AtomNumber["S"]),
    Si29 = c(0:AtomNumber["Si"]),
    Si30 = c(0:AtomNumber["Si"]),
    Cl37 = c(0:AtomNumber["Cl"]),
    Br81 = c(0:AtomNumber["Br"])
  )

  Isotope.Combinations <- Isotope.Combinations %>%
    dplyr::mutate(
      MassSum =
        .data$C13 +
          .data$H2 +
          .data$O17 +
          .data$O18 * 2 +
          .data$S33 +
          .data$S34 * 2 +
          .data$Si29 +
          .data$Si30 * 2 +
          .data$Cl37 * 2 +
          .data$Br81 * 2
    ) %>%
    dplyr::filter((
      .data$O17 + .data$O18) <= AtomNumber["O"] &
      (.data$S33 + .data$S34) <= AtomNumber["S"] &
      (.data$Si29 + .data$Si30) <= AtomNumber["Si"] &
      .data$MassSum <= AtomNumber["N"]) %>%
    dplyr::mutate(
      MassDiff =
        (13.00335 - 12) * .data$C13 +
          (2.0141 - 1.00783) * .data$H2 +
          (16.99913 - 15.99491) * .data$O17 +
          (17.99916 - 15.99491) * .data$O18 +
          (32.97146 - 31.97207) * .data$S33 +
          (33.96787 - 31.97207) * .data$S34 +
          (28.97649 - 27.97693) * .data$Si29 +
          (29.97377 - 27.97693) * .data$Si30 +
          (36.96590 - 34.96885) * .data$Cl37 +
          (80.91629 - 78.91833) * .data$Br81 -
          (15.00011 - 14.00307) * .data$MassSum
    ) %>%
    dplyr::filter(abs(.data$MassDiff) < Mass.Limit)

  for (i in seq_len(nrow(Isotope.Combinations))) {
    p <- stats::dbinom(
      Isotope.Combinations[i, 1],
      AtomNumber["C"],
      CarbonNaturalAbundace[2]
    ) *
      stats::dbinom(
        Isotope.Combinations[i, 2],
        AtomNumber["H"],
        HydrogenNaturalAbundace[2]
      ) *
      stats::dmultinom(
        unlist(c(
          AtomNumber["O"] - sum(Isotope.Combinations[i, 3:4]),
          Isotope.Combinations[i, 3:4]
        )),
        AtomNumber["O"],
        OxygenNaturalAbundace
      ) *
      stats::dmultinom(
        unlist(c(
          AtomNumber["S"] - sum(Isotope.Combinations[i, 5:6]),
          Isotope.Combinations[i, 5:6]
        )),
        AtomNumber["S"],
        SulfurNaturalAbundace
      ) *
      stats::dmultinom(
        unlist(c(
          AtomNumber["Si"] - sum(Isotope.Combinations[i, 7:8]),
          Isotope.Combinations[i, 7:8]
        )),
        AtomNumber["Si"],
        SiliconNaturalAbundace
      ) *
      stats::dbinom(
        Isotope.Combinations[i, 9],
        AtomNumber["Cl"],
        ChlorineNaturalAbundance[2]
      ) *
      stats::dbinom(
        Isotope.Combinations[i, 10],
        AtomNumber["Br"],
        BromineNaturalAbundance[2]
      )

    MSum <- Isotope.Combinations[i, 11]
    for (j in 1:(AtomNumber["N"] + 1 - MSum)) {
      NonTracerMatrix[MSum + j, j] <- NonTracerMatrix[MSum + j, j] + p
    }
  }

  for (i in seq_len(ncol(datamatrix))) {
    CorrectedMatrix[, i] <- stats::coef(nnls::nnls(
      NonTracerMatrix %*% NitrogenMatrix %*% PurityMatrix, ExpMatrix[, i]
    ))
  }

  return(CorrectedMatrix)
}


#' Natural Abundance correction for mass spectrometry data
#'
#' \code{natural_abundance_correction} returns the corrected and normalized
#' intensities of isotopically labeled mass spectrometry data. It was designed
#' to work with input data from
#' \href{https://elucidatainc.github.io/ElMaven/}{El-MAVEN} and
#' \href{http://maven.princeton.edu}{MAVEN} software.
#'
#' C13, H2, and N15 isotopes are supported. The isotopes are detected from the
#' \code{isotopeLabel} column of the input file. The expected label text is
#' \code{C13-label-#}. \code{D-label-#}. or \code{N15-label-#}. Parent
#' (unlabeled) compounds are specified by \code{C12 PARENT}.
#'
#' @param data Path to input data file (xlsx, xls, csv, txt, or tsv) OR
#'   dataframe. If dataframe is specified, specify output_base to output files
#'   automatically written.
#' @param sheet Name of sheet in xlsx file with columns 'compound',
#'   'formula', 'isotopelabel', and one column per sample. Defaults to the
#'   first sheet.
#' @param compound_database Path to compound database in csv format. Only used
#'   for classic MAVEN style input when formula is not specified.
#' @param resolution For Exactive, the resolution is 100000, defined at Mw 200
#' @param resolution_defined_at Mw at which the resolution is defined, default
#'   200 Mw
#' @param purity Isotope purity, default: Carbon 0.99; Deuterium 0.98;
#'   Nitrogen 0.99
#' @param output_base Path to basename of output file, default is the basename
#'   of the input path. `_corrected` will be appended. If `FALSE` then no
#'   output file is written.
#' @param output_filetype Filetype of the output file, one of: 'xls', xlsx',
#'   'csv', or 'tsv'. The default is 'xlsx'.
#' @param columns_to_skip Specify column heading to skip. All other columns not
#'   named 'compound', 'formula', and 'isotopelabel' will be assumed to be
#'   sample names.
#' @param report_pool_size_before_df Report PoolSizeBeforeDF, default = FALSE
#' @param path Deprecated. Specify path to input data file (alias for `data`).
#' @importFrom rlang .data
#' @return Named list of matrices: 'Corrected', 'Normalized',
#'   'PoolBeforeDF', and 'PoolAfterDF'.
#' @export
#' @examples
#' \dontrun{
#' natural_abundance_correction("inst/extdata/C_Sample_Input_Simple.xlsx",
#'   Resolution = 100000, ResDefAt = 200
#' )
#' }
natural_abundance_correction <- function(data,
                                         sheet = NULL,
                                         compound_database = NULL,
                                         output_base = NULL,
                                         output_filetype = "xlsx",
                                         columns_to_skip = NULL,
                                         resolution,
                                         resolution_defined_at = 200,
                                         purity = NULL,
                                         report_pool_size_before_df = FALSE,
                                         path = NULL) {
  default_purity <- list("C" = 0.99, "D" = 0.98, "N" = 0.99)

  # Specify `path` as named variable, backwards compatibility
  if (missing(data)) {
    if (!missing(path) || path != "") {
      data <- path
    }
  }

  if (missing(data)) {
    stop("Must specify 'data' (data frame or path to input file)")
  }

  # Determine if data is data.frame or path
  if (is.data.frame(data)) {
    # Use data frame directly
    # Created "cleaned" data frame with normalized labels
    cleaned_dataframe <- clean_data_frame(
      data,
      columns_to_skip = columns_to_skip
    )

    # Determine isotope
    isotope <- determine_isotope(cleaned_dataframe$isotope_label)$isotope

    input_data <- (list(
      original = tibble::as_tibble(data),
      cleaned = cleaned_dataframe,
      isotope = isotope
    ))
    # Only output files if output_base is specified
    if (is.null(output_base)) {
      output_base <- FALSE
    }
  } else if (is.character(data) & length(data) == 1) {
    # Data is a string
    if (!file.exists(data)) {
      stop(sprintf("Unable to find file '%s'", path))
    }
    if (is.null(output_filetype)) {
      output_filetype <- tools::file_ext(path)
    }
    if (!(output_filetype %in% c("xls", "xlsx", "csv", "tsv"))) {
      stop(paste("Unsupported output_filetype: '", output_filetype, "'",
        sep = ""
      ))
    }
    # Parse input file
    input_data <- read_elmaven(
      path = path, sheet = sheet,
      compound_database = compound_database,
      columns_to_skip = columns_to_skip
    )
  } else {
    stop(sprintf("Argument 'data' must be either a data.frame or a string"))
  }

  if (missing(resolution)) {
    stop("Must specify 'resolution'")
  }
  if (!is.numeric(resolution)) {
    stop("'resolution' must be an integer")
  }
  if (as.numeric(resolution) %% 1 != 0) {
    stop("'resolution' must be an integer")
  }

  if (!is.numeric(resolution_defined_at)) {
    stop("'resolution_defined_at' must be an integer")
  }
  if (resolution_defined_at %% 1 != 0) {
    stop("'resolution_defined_at' must be an integer")
  }
  if (!is.null(purity)) {
    if (!is.numeric(resolution_defined_at)) {
      stop("'purity' must be a number")
    } else if ((purity < 0) | (purity > 1)) {
      stop("'purity' must be between 0 and 1")
    }
  }

  # Determine sample columns
  sample_col_names <- names(input_data$cleaned)[
    which(!(tolower(names(input_data$cleaned))
    %in% tolower(
        c("compound", "formula", "isotope_label", "label_index", "metaGroupId")
      )))
  ]

  if (!(input_data$isotope %in% names(default_purity))) {
    stop(paste(
      "Unsupported isotope '",
      input_data$isotope,
      "' detected",
      sep = ""
    ))
  }
  if (is.null(purity)) {
    purity <- default_purity[[input_data$isotope]]
  }

  # Setup empty matrices for output
  # TODO Refactor this to preallocate or better yet, use sapply
  OutputMatrix <- matrix(NA, nrow = 0, ncol = length(sample_col_names))
  OutputPercentageMatrix <- matrix(
    NA,
    nrow = 0,
    ncol = length(sample_col_names)
  )
  OutputPoolBefore <- matrix(NA, nrow = 0, ncol = length(sample_col_names))
  OutputPoolAfter <- matrix(NA, nrow = 0, ncol = length(sample_col_names))
  OutputCompound <- NULL
  OutputLabel <- NULL
  OutputPoolCompound <- NULL

  for (i in unique(input_data$cleaned$metaGroupId)) {
    CurrentMetabolite <- dplyr::filter(
      input_data$cleaned,
      .data$metaGroupId == i
    )
    CurrentCompoundName <- CurrentMetabolite$compound[1]
    Formula <- as.character(CurrentMetabolite$formula[1])
    if (length(Formula) == 0 || is.na(Formula)) {
      print(paste("The formula of", i, "is unknown", sep = " "))
      break
    }
    DataMatrix <- data.matrix(dplyr::select(
      CurrentMetabolite, -.data$compound, -.data$formula,
      -.data$isotope_label, -.data$label_index, -.data$metaGroupId
    ))
    DataMatrix[is.na(DataMatrix)] <- 0
    if (input_data$isotope == "C") {
      Corrected <- carbon_isotope_correction(
        Formula,
        DataMatrix,
        CurrentMetabolite$label_index,
        Resolution = resolution,
        ResDefAt = resolution_defined_at,
        purity = purity,
        ReportPoolSize = report_pool_size_before_df
      )
    } else if (input_data$isotope == "D") {
      Corrected <- deuterium_isotope_correction(
        Formula,
        DataMatrix,
        CurrentMetabolite$label_index,
        Resolution = resolution,
        ResDefAt = resolution_defined_at,
        purity = purity,
        ReportPoolSize = report_pool_size_before_df
      )
    } else if (input_data$isotope == "N") {
      Corrected <- nitrogen_isotope_correction(
        Formula,
        DataMatrix,
        CurrentMetabolite$label_index,
        Resolution = resolution,
        ResDefAt = resolution_defined_at,
        purity = purity,
        ReportPoolSize = report_pool_size_before_df
      )
    } else {
      stop(paste(
        "Unsupported isotope '",
        input_data$isotope,
        "' detected",
        sep = ""
      ))
    }
    CorrectedPercentage <- scale(
      Corrected,
      scale = colSums(Corrected),
      center = FALSE
    )
    OutputMatrix <- rbind(OutputMatrix, Corrected)
    OutputPercentageMatrix <- rbind(
      OutputPercentageMatrix,
      CorrectedPercentage
    )
    OutputPoolBefore <- rbind(OutputPoolBefore, colSums(DataMatrix))
    OutputPoolAfter <- rbind(OutputPoolAfter, colSums(Corrected))
    OutputCompound <- append(
      OutputCompound,
      rep(
        CurrentCompoundName,
        nrow(Corrected)
      )
    )
    OutputLabel <- append(OutputLabel, c(0:(nrow(Corrected) - 1)))
    OutputPoolCompound <- append(OutputPoolCompound, CurrentCompoundName)
  }

  compound_label_tbl <- dplyr::tibble(OutputCompound, OutputLabel)
  colnames(OutputMatrix) <- sample_col_names
  OutputDF <- dplyr::bind_cols(
    compound_label_tbl,
    dplyr::as_tibble(OutputMatrix,
      .name_repair = "minimal"
    )
  )
  colnames(OutputPercentageMatrix) <- sample_col_names
  OutputPercentageDF <- dplyr::bind_cols(
    compound_label_tbl,
    dplyr::as_tibble(OutputPercentageMatrix,
      .name_repair = "minimal"
    )
  )
  colnames(OutputPoolBefore) <- sample_col_names
  OutputPoolBeforeDF <- dplyr::bind_cols(
    dplyr::tibble(OutputPoolCompound),
    dplyr::as_tibble(OutputPoolBefore,
      .name_repair = "minimal"
    )
  )
  colnames(OutputPoolAfter) <- sample_col_names
  OutputPoolAfterDF <- dplyr::bind_cols(
    dplyr::tibble(OutputPoolCompound),
    dplyr::as_tibble(OutputPoolAfter,
      .name_repair = "minimal"
    )
  )
  names(OutputDF) <- c(
    "Compound",
    paste(input_data$isotope, "Label", sep = "_"),
    sample_col_names
  )
  names(OutputPercentageDF) <- names(OutputDF)
  names(OutputPoolBeforeDF) <- c("Compound", sample_col_names)
  names(OutputPoolAfterDF) <- c("Compound", sample_col_names)

  OutputDataFrames <- list(
    "Original" = input_data$original,
    "Corrected" = OutputDF,
    "Normalized" = OutputPercentageDF,
    "PoolAfterDF" = OutputPoolAfterDF
  )
  if (report_pool_size_before_df) {
    OutputDataFrames$PoolBeforeDF <- OutputPoolBeforeDF
  }

  if (!identical(FALSE, output_base)) {
    if (is.null(output_base)) {
      output_base <- path
    }
    write_output(OutputDataFrames, output_base, filetype = output_filetype)
  }

  return(OutputDataFrames)
}
