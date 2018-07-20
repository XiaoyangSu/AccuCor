#' Natural Abundance carbon isotope correction for one metabolite
#'
#' @param formula String representing molecular formula
#' @param datamatrix Matrix of abundnaces for each sample for each isotope
#' @param label vector of integer labels
#' @param Resolution For Exactive, the Resolution is 100000, defined at Mw 200
#' @param ResDefAt Resolution defined at (in Mw), e.g. 200 Mw
#' @param purity Carbon 13 purity, default: 0.99
#' @param ReportPoolSize default: TRUE
#' @importFrom rlang .data
#' @return Named list of matrices: 'Corrected', 'Normalized',
#'      'PoolBeforeDF', and 'PoolAfterDF'.
carbon_isotope_correction <- function(formula, datamatrix, label, Resolution,
                                      ResDefAt, purity=0.99, ReportPoolSize=TRUE) {

  CarbonNaturalAbundace <- c(0.9893, 0.0107)
  HydrogenNaturalAbundace <- c(0.999885, 0.000115)
  NitrogenNaturalAbundace <- c(0.99636, 0.00364)
  OxygenNaturalAbundace <- c(0.99757, 0.00038, 0.00205)
  SulfurNaturalAbundace <- c(0.9493, 0.00762, 0.0429)

  AtomNumber <- rep(0,6)
  names(AtomNumber) <- c("C","H","N","O","P","S")
  MassDifference <- abs(c((2.0141-1.00783),(15.00011-14.00307),(16.99913-15.99491),
                          (17.99916-15.99491),(32.97146-31.97207),(33.96787-31.97207))-
                          ((13.00335-12)*c(1,1,1,2,1,2)))
  names(MassDifference) <- c("H2","N15","O17","O18","S33","S34")
  CorrectionLimit <- rep(0,6)
  names(CorrectionLimit) <- c("H2","N15","O17","O18","S33","S34")

  if.not.null <- function(x) if(!is.null(x)) x else 0
  AtomNumber["C"] <- if.not.null(unlist(gsubfn::strapply(formula, "(C)(\\d*)", ~ as.numeric(if (nchar(..2)) ..2 else 1))))
  AtomNumber["H"] <- if.not.null(unlist(gsubfn::strapply(formula, "(H)(\\d*)", ~ as.numeric(if (nchar(..2)) ..2 else 1))))
  AtomNumber["N"] <- if.not.null(unlist(gsubfn::strapply(formula, "(N)(\\d*)", ~ as.numeric(if (nchar(..2)) ..2 else 1))))
  AtomNumber["O"] <- if.not.null(unlist(gsubfn::strapply(formula, "(O)(\\d*)", ~ as.numeric(if (nchar(..2)) ..2 else 1))))
  AtomNumber["P"] <- if.not.null(unlist(gsubfn::strapply(formula, "(P)(\\d*)", ~ as.numeric(if (nchar(..2)) ..2 else 1))))
  AtomNumber["S"] <- if.not.null(unlist(gsubfn::strapply(formula, "(S)(\\d*)", ~ as.numeric(if (nchar(..2)) ..2 else 1))))

  MolecularWeight <- sum(AtomNumber*c(12,1,14,16,31,32))

  CorrectionLimit <- floor(MolecularWeight^(3/2)*1.66/(Resolution*sqrt(ResDefAt))/MassDifference)

  ExpMatrix <- matrix(0, ncol=ncol(datamatrix), nrow=AtomNumber["C"]+1)
  CorrectedMatrix <- matrix(0, ncol=ncol(datamatrix), nrow=AtomNumber["C"]+1)
  for (i in 1:length(label)) {
    ExpMatrix[label[i]+1,] <- datamatrix[i,]
  }


  PurityMatrix <- diag(AtomNumber["C"]+1)
  CarbonMatrix <- diag(AtomNumber["C"]+1)
  NitrogenMatrix <- diag(AtomNumber["C"]+1)
  HydrogenMatrix <- diag(AtomNumber["C"]+1)
  OxygenMatrix <- matrix(0,ncol=(AtomNumber["C"]+1),nrow=(AtomNumber["C"]+1))
  SulfurMatrix <- matrix(0,ncol=(AtomNumber["C"]+1),nrow=(AtomNumber["C"]+1))


  for(i in 1:(AtomNumber["C"]+1)){
    PurityMatrix[i,] <- sapply(0:(AtomNumber["C"]), function(x) stats::dbinom(x-i+1, x , (1-purity)))
  }

  for(i in 1:(AtomNumber["C"]+1)){
    CarbonMatrix[,i] <- sapply(0:AtomNumber["C"], function(x) stats::dbinom(x-i+1, AtomNumber["C"]-i+1 , CarbonNaturalAbundace[2]))
  }

  for(j in 0:min(AtomNumber["N"], CorrectionLimit["N15"], AtomNumber["C"]))
    for(i in 1:(AtomNumber["C"]-j+1)){
      NitrogenMatrix[(i+j),i] <- stats::dbinom(j, AtomNumber["N"], NitrogenNaturalAbundace[2])
    }

  for(j in 0:min(AtomNumber["H"], CorrectionLimit["H2"], AtomNumber["C"]))
    for(i in 1:(AtomNumber["C"]-j+1)){
      HydrogenMatrix[(i+j),i] <- stats::dbinom(j, AtomNumber["H"], HydrogenNaturalAbundace[2])
    }

  for(i in 0:min(AtomNumber["O"],CorrectionLimit["O17"])) {
    for(j in 0:min(AtomNumber["O"],CorrectionLimit["O18"])){
      k<-(i+j*2)
      if ((i+j)>AtomNumber["O"]|k>AtomNumber["C"]) {
        break
      }
      else {
        for (m in 1:(AtomNumber["C"]-k+1)) {
          OxygenMatrix[(m+k),m] <- OxygenMatrix[(m+k),m] + stats::dmultinom(c((AtomNumber["O"]-i-j),i,j), AtomNumber["O"], OxygenNaturalAbundace)
        }
      }
    }
  }

  for(i in 0:min(AtomNumber["S"],CorrectionLimit["S33"])) {
    for(j in 0:min(AtomNumber["S"],CorrectionLimit["S34"])){
      k<-(i+j*2)
      if ((i+j)>AtomNumber["S"]|k>AtomNumber["C"]) {
        break
      }
      else {
        for (m in 1:(AtomNumber["C"]-k+1)) {
          SulfurMatrix[(m+k),m] <- SulfurMatrix[(m+k),m] + stats::dmultinom(c((AtomNumber["S"]-i-j),i,j), AtomNumber["S"], SulfurNaturalAbundace)
        }
      }
    }
  }

  for(i in 1:ncol(datamatrix)) {
    CorrectedMatrix[,i] <- stats::coef(nnls::nnls(SulfurMatrix %*% OxygenMatrix %*% NitrogenMatrix %*%
                                           HydrogenMatrix %*% CarbonMatrix %*% PurityMatrix, ExpMatrix[,i]))
  }

  return(CorrectedMatrix)

}

#' Natural Abundance deuterium isotope correction for one metabolite
#'
#' @param formula String representing molecular formula
#' @param datamatrix Matrix of abundnaces for each sample for each isotope
#' @param label vector of integer labels
#' @param Resolution For Exactive, the Resolution is 100000, defined at Mw 200
#' @param ResDefAt Resolution defined at (in Mw), e.g. 200 Mw
#' @param purity Deuterium purity, default: 0.99
#' @param ReportPoolSize default: TRUE
#' @importFrom rlang .data
#' @return Named list of matrices: 'Corrected', 'Normalized',
#'      'PoolBeforeDF', and 'PoolAfterDF'.
deuterium_isotope_correction <- function(formula, datamatrix, label, Resolution, ResDefAt,
                                         purity=0.99, ReportPoolSize=TRUE) {

  CarbonNaturalAbundace <- c(0.9893, 0.0107)
  HydrogenNaturalAbundace <- c(0.999885, 0.000115)
  NitrogenNaturalAbundace <- c(0.99636, 0.00364)
  OxygenNaturalAbundace <- c(0.99757, 0.00038, 0.00205)
  SulfurNaturalAbundace <- c(0.9493, 0.00762, 0.0429)

  AtomNumber <- rep(0,6)
  names(AtomNumber) <- c("C","H","N","O","P","S")
  MassDifference <- abs(c((13.00335-12),(15.00011-14.00307),(16.99913-15.99491),
                          (17.99916-15.99491),(32.97146-31.97207),(33.96787-31.97207))-
                          ((2.0141-1.00783)*c(1,1,1,2,1,2)))
  names(MassDifference) <- c("C13","N15","O17","O18","S33","S34")
  CorrectionLimit <- rep(0,6)
  names(CorrectionLimit) <- c("C13","N15","O17","O18","S33","S34")

  if.not.null <- function(x) if(!is.null(x)) x else 0
  AtomNumber["C"] <- if.not.null(unlist(gsubfn::strapply(formula, "(C)(\\d*)", ~ as.numeric(if (nchar(..2)) ..2 else 1))))
  AtomNumber["H"] <- if.not.null(unlist(gsubfn::strapply(formula, "(H)(\\d*)", ~ as.numeric(if (nchar(..2)) ..2 else 1))))
  AtomNumber["N"] <- if.not.null(unlist(gsubfn::strapply(formula, "(N)(\\d*)", ~ as.numeric(if (nchar(..2)) ..2 else 1))))
  AtomNumber["O"] <- if.not.null(unlist(gsubfn::strapply(formula, "(O)(\\d*)", ~ as.numeric(if (nchar(..2)) ..2 else 1))))
  AtomNumber["P"] <- if.not.null(unlist(gsubfn::strapply(formula, "(P)(\\d*)", ~ as.numeric(if (nchar(..2)) ..2 else 1))))
  AtomNumber["S"] <- if.not.null(unlist(gsubfn::strapply(formula, "(S)(\\d*)", ~ as.numeric(if (nchar(..2)) ..2 else 1))))

  MolecularWeight <- sum(AtomNumber*c(12,1,14,16,31,32))

  CorrectionLimit <- floor(MolecularWeight^(3/2)*1.66/(Resolution*sqrt(ResDefAt))/MassDifference)

  ExpMatrix <- matrix(0, ncol=ncol(datamatrix), nrow=AtomNumber["H"]+1)
  CorrectedMatrix <- matrix(0, ncol=ncol(datamatrix), nrow=AtomNumber["H"]+1)
  for (i in 1:length(label)) {
    ExpMatrix[label[i]+1,] <- datamatrix[i,]
  }


  PurityMatrix <- diag(AtomNumber["H"]+1)
  CarbonMatrix <- diag(AtomNumber["H"]+1)
  NitrogenMatrix <- diag(AtomNumber["H"]+1)
  HydrogenMatrix <- diag(AtomNumber["H"]+1)
  OxygenMatrix <- matrix(0,ncol=(AtomNumber["H"]+1),nrow=(AtomNumber["H"]+1))
  SulfurMatrix <- matrix(0,ncol=(AtomNumber["H"]+1),nrow=(AtomNumber["H"]+1))


  for(i in 1:(AtomNumber["H"]+1)){
    PurityMatrix[i,] <- sapply(0:AtomNumber["H"], function(x) stats::dbinom(x-i+1, x , (1-purity)))
  }

  for(j in 0:min(AtomNumber["C"], CorrectionLimit["C13"], AtomNumber["H"]))
    for(i in 1:(AtomNumber["H"]-j+1)){
      CarbonMatrix[(i+j),i] <- stats::dbinom(j, AtomNumber["C"], CarbonNaturalAbundace[2])
    }

  for(i in 1:(AtomNumber["H"]+1)){
    HydrogenMatrix[,i] <- sapply(0:AtomNumber["H"], function(x) stats::dbinom(x-i+1, AtomNumber["H"]-i+1 , HydrogenNaturalAbundace[2]))
  }

  for(j in 0:min(AtomNumber["N"], CorrectionLimit["N15"], AtomNumber["H"]))
    for(i in 1:(AtomNumber["H"]-j+1)){
      NitrogenMatrix[(i+j),i] <- stats::dbinom(j, AtomNumber["N"], NitrogenNaturalAbundace[2])
    }

  for(i in 0:min(AtomNumber["O"],CorrectionLimit["O17"])) {
    for(j in 0:min(AtomNumber["O"],CorrectionLimit["O18"])){
      k<-(i+j*2)
      if ((i+j)>AtomNumber["O"]|k>AtomNumber["H"]) {
        break
      }
      else {
        for (m in 1:(AtomNumber["H"]-k+1)) {
          OxygenMatrix[(m+k),m] <- OxygenMatrix[(m+k),m] + stats::dmultinom(c((AtomNumber["O"]-i-j),i,j), AtomNumber["O"], OxygenNaturalAbundace)
        }
      }
    }
  }

  for(i in 0:min(AtomNumber["S"],CorrectionLimit["S33"])) {
    for(j in 0:min(AtomNumber["S"],CorrectionLimit["S34"])){
      k<-(i+j*2)
      if ((i+j)>AtomNumber["S"]|k>AtomNumber["H"]) {
        break
      }
      else {
        for (m in 1:(AtomNumber["H"]-k+1)) {
          SulfurMatrix[(m+k),m] <- SulfurMatrix[(m+k),m] + stats::dmultinom(c((AtomNumber["S"]-i-j),i,j), AtomNumber["S"], SulfurNaturalAbundace)
        }
      }
    }
  }

  for(i in 1:ncol(datamatrix)) {
    CorrectedMatrix[,i] <- stats::coef(nnls::nnls(SulfurMatrix %*% CarbonMatrix %*% OxygenMatrix %*%
                                       NitrogenMatrix %*% HydrogenMatrix %*% PurityMatrix, ExpMatrix[,i]))
  }

  return(CorrectedMatrix)

}


#' Natural Abundance deuterium isotope correction for one metabolite
#'
#' @param formula String representing molecular formula
#' @param datamatrix Matrix of abundnaces for each sample for each isotope
#' @param label vector of integer labels
#' @param Resolution For Exactive, the Resolution is 100000, defined at Mw 200
#' @param ResDefAt Resolution defined at (in Mw), e.g. 200 Mw
#' @param purity Nitrogen purity, default: 0.99
#' @param ReportPoolSize default: TRUE
#' @importFrom rlang .data
#' @return Named list of matrices: 'Corrected', 'Normalized',
#'      'PoolBeforeDF', and 'PoolAfterDF'.
nitrogen_isotope_correction <- function(formula, datamatrix, label, Resolution, ResDefAt,
                                        purity=0.99, ReportPoolSize=TRUE) {

  CarbonNaturalAbundace <- c(0.9893, 0.0107)
  HydrogenNaturalAbundace <- c(0.999885, 0.000115)
  NitrogenNaturalAbundace <- c(0.99636, 0.00364)
  OxygenNaturalAbundace <- c(0.99757, 0.00038, 0.00205)
  SulfurNaturalAbundace <- c(0.9493, 0.00762, 0.0429)

  AtomNumber <- rep(0,6)
  names(AtomNumber) <- c("C","H","N","O","P","S")
  MassDifference <- abs(c((13.00335-12),(2.0141-1.00783),(16.99913-15.99491),
                          (17.99916-15.99491),(32.97146-31.97207),(33.96787-31.97207))-
                          ((15.00011-14.00307)*c(1,1,1,2,1,2)))
  names(MassDifference) <- c("C13","H2","O17","O18","S33","S34")
  CorrectionLimit <- rep(0,6)
  names(CorrectionLimit) <- c("C13","H2","O17","O18","S33","S34")

  if.not.null <- function(x) if(!is.null(x)) x else 0
  AtomNumber["C"] <- if.not.null(unlist(gsubfn::strapply(formula, "(C)(\\d*)", ~ as.numeric(if (nchar(..2)) ..2 else 1))))
  AtomNumber["H"] <- if.not.null(unlist(gsubfn::strapply(formula, "(H)(\\d*)", ~ as.numeric(if (nchar(..2)) ..2 else 1))))
  AtomNumber["N"] <- if.not.null(unlist(gsubfn::strapply(formula, "(N)(\\d*)", ~ as.numeric(if (nchar(..2)) ..2 else 1))))
  AtomNumber["O"] <- if.not.null(unlist(gsubfn::strapply(formula, "(O)(\\d*)", ~ as.numeric(if (nchar(..2)) ..2 else 1))))
  AtomNumber["P"] <- if.not.null(unlist(gsubfn::strapply(formula, "(P)(\\d*)", ~ as.numeric(if (nchar(..2)) ..2 else 1))))
  AtomNumber["S"] <- if.not.null(unlist(gsubfn::strapply(formula, "(S)(\\d*)", ~ as.numeric(if (nchar(..2)) ..2 else 1))))

  MolecularWeight <- sum(AtomNumber*c(12,1,14,16,31,32))

  CorrectionLimit <- floor(MolecularWeight^(3/2)*1.66/(Resolution*sqrt(ResDefAt))/MassDifference)

  ExpMatrix <- matrix(0, ncol=ncol(datamatrix), nrow=AtomNumber["N"]+1)
  CorrectedMatrix <- matrix(0, ncol=ncol(datamatrix), nrow=AtomNumber["N"]+1)
  for (i in 1:length(label)) {
    ExpMatrix[label[i]+1,] <- datamatrix[i,]
  }


  PurityMatrix <- diag(AtomNumber["N"]+1)
  CarbonMatrix <- diag(AtomNumber["N"]+1)
  NitrogenMatrix <- diag(AtomNumber["N"]+1)
  HydrogenMatrix <- diag(AtomNumber["N"]+1)
  OxygenMatrix <- matrix(0,ncol=(AtomNumber["N"]+1),nrow=(AtomNumber["N"]+1))
  SulfurMatrix <- matrix(0,ncol=(AtomNumber["N"]+1),nrow=(AtomNumber["N"]+1))


  for(i in 1:(AtomNumber["N"]+1)){
    PurityMatrix[i,] <- sapply(0:(AtomNumber["N"]), function(x) stats::dbinom(x-i+1, x , (1-purity)))
  }

  for(i in 1:(AtomNumber["N"]+1)){
    NitrogenMatrix[,i] <- sapply(0:AtomNumber["N"], function(x) stats::dbinom(x-i+1, AtomNumber["N"]-i+1 , NitrogenNaturalAbundace[2]))
  }

  for(j in 0:min(AtomNumber["C"], CorrectionLimit["C13"], AtomNumber["N"]))
    for(i in 1:(AtomNumber["N"]-j+1)){
      CarbonMatrix[(i+j),i] <- stats::dbinom(j, AtomNumber["C"], CarbonNaturalAbundace[2])
    }

  for(j in 0:min(AtomNumber["H"], CorrectionLimit["H2"], AtomNumber["N"]))
    for(i in 1:(AtomNumber["N"]-j+1)){
      HydrogenMatrix[(i+j),i] <- stats::dbinom(j, AtomNumber["H"], HydrogenNaturalAbundace[2])
    }

  for(i in 0:min(AtomNumber["O"],CorrectionLimit["O17"])) {
    for(j in 0:min(AtomNumber["O"],CorrectionLimit["O18"])){
      k<-(i+j*2)
      if ((i+j)>AtomNumber["O"]|k>AtomNumber["N"]) {
        break
      }
      else {
        for (m in 1:(AtomNumber["N"]-k+1)) {
          OxygenMatrix[(m+k),m] <- OxygenMatrix[(m+k),m] + stats::dmultinom(c((AtomNumber["O"]-i-j),i,j), AtomNumber["O"], OxygenNaturalAbundace)
        }
      }
    }
  }

  for(i in 0:min(AtomNumber["S"],CorrectionLimit["S33"])) {
    for(j in 0:min(AtomNumber["S"],CorrectionLimit["S34"])){
      k<-(i+j*2)
      if ((i+j)>AtomNumber["S"]|k>AtomNumber["N"]) {
        break
      }
      else {
        for (m in 1:(AtomNumber["N"]-k+1)) {
          SulfurMatrix[(m+k),m] <- SulfurMatrix[(m+k),m] + stats::dmultinom(c((AtomNumber["S"]-i-j),i,j), AtomNumber["S"], SulfurNaturalAbundace)
        }
      }
    }
  }

  for(i in 1:ncol(datamatrix)) {
    CorrectedMatrix[,i] <- stats::coef(nnls::nnls(SulfurMatrix %*% OxygenMatrix %*% CarbonMatrix %*%
                                       HydrogenMatrix %*% NitrogenMatrix %*% PurityMatrix, ExpMatrix[,i]))
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
#' @param path Path to xlsx file.
#' @param sheet Name of sheet in xlsx file with columns 'compound',
#'   'formula', 'isotopelabel', and one column per sample.
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
#' @importFrom rlang .data
#' @return Named list of matrices: 'Corrected', 'Normalized',
#'   'PoolBeforeDF', and 'PoolAfterDF'.
#' @export
#' @examples
#' \dontrun{
#' natural_abundance_correction("inst/extdata/C_Sample_Input_Simple.xlsx",
#'                              Resolution=100000, ResDefAt=200)
#' }
natural_abundance_correction <- function(path,
                                         sheet = NULL,
                                         compound_database = NULL,
                                         output_base = NULL,
                                         output_filetype = 'xlsx',
                                         columns_to_skip = NULL,
                                         resolution,
                                         resolution_defined_at = 200,
                                         purity = NULL,
                                         report_pool_size_before_df = FALSE) {

  default_purity <- list("C" = 0.99, "D" = 0.98, "N" = 0.99)

  if (missing(path) || path == "") {
    stop("Must specify 'path' to input file")
  }

  if (!file.exists(path)) {
    stop(sprintf("Unable to find file '%s'", path))
  }

  if (missing(resolution)) {
    stop("Must specify 'resolution'")
  }
  if (!is.numeric(resolution)) {
    stop("'resolution' must be an integer")
  }
  if (as.numeric(resolution)%%1!=0) {
    stop("'resolution' must be an integer")
  }

  if (!is.numeric(resolution_defined_at)) {
    stop("'resolution_defined_at' must be an integer")
  }
  if (resolution_defined_at%%1!=0) {
    stop("'resolution_defined_at' must be an integer")
  }
  if (!is.null(purity)) {
    if (!is.numeric(resolution_defined_at)) {
      stop("'purity' must be a number")
    } else if ((purity < 0) | (purity > 1)) {
      stop("'purity' must be between 0 and 1")
    }
  }

  if (is.null(output_filetype)) {
    output_filetype = tools::file_ext(path)
  }
  if (! (output_filetype %in% c('xls', 'xlsx', 'csv', 'tsv'))) {
    stop(paste("Unsupported output_filetype: '", output_filetype, "'",
               sep = ""))
  }

  input_data <- read_elmaven(path = path, sheet = sheet,
                             compound_database = compound_database,
                             columns_to_skip = columns_to_skip)
  sample_col_names <- names(input_data$cleaned)[
    which( !(tolower(names(input_data$cleaned))
             %in% tolower(
               c("compound", "formula", "isotope_label", "label_index", "metaGroupId"))))]

  if ( !(input_data$isotope %in% names(default_purity)) ) {
    stop(paste("Unsupported isotope '", input_data$isotope, "' detected", sep = ""))
  }
  if (is.null(purity)) {
    purity = default_purity[[input_data$isotope]]
  }

  # Setup empty matrices for output
  # TODO Refactor this to preallocate or better yet, use sapply
  OutputMatrix <- matrix(NA, nrow = 0, ncol = length(sample_col_names))
  OutputPercentageMatrix <- matrix(NA, nrow = 0, ncol = length(sample_col_names))
  OutputPoolBefore <- matrix(NA, nrow = 0, ncol = length(sample_col_names))
  OutputPoolAfter <- matrix(NA, nrow = 0, ncol = length(sample_col_names))
  OutputCompound <- NULL
  OutputLabel <- NULL
  OutputPoolCompound <- NULL

  for (i in unique(input_data$cleaned$metaGroupId)) {
    CurrentMetabolite <- dplyr::filter(input_data$cleaned, .data$metaGroupId==i)
    CurrentCompoundName <- CurrentMetabolite$compound[1]
    Formula=as.character(CurrentMetabolite$formula[1])
    if(length(Formula)==0 || is.na(Formula)) {
      print(paste("The formula of",i,"is unknown",sep=" "))
      break
    }
    DataMatrix <- data.matrix(dplyr::select(CurrentMetabolite, -.data$compound, -.data$formula,
                                            -.data$isotope_label, -.data$label_index, -.data$metaGroupId)
    )
    DataMatrix[is.na(DataMatrix)] <- 0
    if (input_data$isotope == "C") {
      Corrected <- carbon_isotope_correction(Formula, DataMatrix, CurrentMetabolite$label_index,
                                             Resolution = resolution, ResDefAt = resolution_defined_at,
                                             purity = purity, ReportPoolSize = report_pool_size_before_df)
    } else if (input_data$isotope == "D") {
      Corrected <- deuterium_isotope_correction(Formula, DataMatrix, CurrentMetabolite$label_index,
                                             Resolution = resolution, ResDefAt = resolution_defined_at,
                                             purity = purity, ReportPoolSize = report_pool_size_before_df)
    } else if (input_data$isotope == "N") {
      Corrected <- nitrogen_isotope_correction(Formula, DataMatrix, CurrentMetabolite$label_index,
                                                Resolution = resolution, ResDefAt = resolution_defined_at,
                                                purity = purity, ReportPoolSize = report_pool_size_before_df)
    } else {
      stop(paste("Unsupported isotope '", input_data$isotope, "' detected", sep = ""))
    }
    CorrectedPercentage <- scale(Corrected,scale=colSums(Corrected),center=FALSE)
    OutputMatrix <- rbind(OutputMatrix, Corrected)
    OutputPercentageMatrix <- rbind(OutputPercentageMatrix, CorrectedPercentage)
    OutputPoolBefore <- rbind(OutputPoolBefore, colSums(DataMatrix))
    OutputPoolAfter <- rbind(OutputPoolAfter, colSums(Corrected))
    OutputCompound <- append(OutputCompound, rep(CurrentCompoundName, nrow(Corrected)))
    OutputLabel <- append(OutputLabel, c(0:(nrow(Corrected)-1)))
    OutputPoolCompound <- append(OutputPoolCompound, CurrentCompoundName)
  }

  compound_label_tbl <- dplyr::tibble(OutputCompound, OutputLabel)
  OutputDF <- dplyr::bind_cols(compound_label_tbl,
                               dplyr::as_tibble(OutputMatrix))
  OutputPercentageDF <- dplyr::bind_cols(compound_label_tbl,
                                         dplyr::as_tibble(OutputPercentageMatrix))
  OutputPoolBeforeDF <- dplyr::bind_cols(dplyr::tibble(OutputPoolCompound),
                                         dplyr::as_tibble(OutputPoolBefore))
  OutputPoolAfterDF <- dplyr::bind_cols(dplyr::tibble(OutputPoolCompound),
                                        dplyr::as_tibble(OutputPoolAfter))
  names(OutputDF) <- c("Compound",
                       paste(input_data$isotope,"Label", sep="_"),
                       sample_col_names)
  names(OutputPercentageDF) <- names(OutputDF)
  names(OutputPoolBeforeDF) <- c("Compound", sample_col_names)
  names(OutputPoolAfterDF) <- c("Compound", sample_col_names)

  OutputDataFrames <- list("Original" = input_data$original,
                           "Corrected" = OutputDF,
                           "Normalized" = OutputPercentageDF,
                           "PoolAfterDF" = OutputPoolAfterDF)
  if(report_pool_size_before_df) {
    OutputDataFrames$PoolBeforeDF = OutputPoolBeforeDF
  }

  if(!identical(FALSE, output_base)) {
    if(is.null(output_base)) {
      output_base = path
    }
    write_output(OutputDataFrames, output_base, filetype = output_filetype)
  }

  return(OutputDataFrames)
}
