#' Natural Abundance isotope correction for one metabolite
#'
#' @param formula String representing molecular formula
#' @param datamatrix Matrix of abundnaces for each sample for each isotope
#' @param label vector of integer labels
#' @param Resolution For Exactive, the Resolution is 100000, defined at Mw 200
#' @param ResDefAt Resolution defined at (in Mw), e.g. 200 Mw
#' @param C13Purity Carbon 13 purity, default: 0.99
#' @param ReportPoolSize default: TRUE
#' @importFrom rlang .data
#' @return Named list of matrices: 'Corrected', 'Normalized',
#'      'PoolBeforeDF', and 'PoolAfterDF'.
#' @examples
#' \dontrun{
#' carbon_correction("ExcelFile", Resolution=100000, ResDefAt=200)
#' }
IsotopeCorrection <- function(formula, datamatrix, label, Resolution, ResDefAt,
                              C13Purity=0.99, ReportPoolSize=TRUE) {

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
    PurityMatrix[i,] <- sapply(0:(AtomNumber["C"]), function(x) stats::dbinom(x-i+1, x , (1-C13Purity)))
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

#' Natural Abundance correction for Carbon labeled samples
#'
#' @param path Path to xlsx file.
#' @param sheet Name of sheet in xlsx file with columns 'compound',
#'      'formula', 'isotopelabel', and one column per sample.
#' @param output_path Path to output file, default is input path with
#'      `_corrected` appended. If `FALSE` then no output file is written.
#' @param ColumnsToSkip Specify column heading to skip. All other columns not
#'      named 'compound', 'formula', and 'isotopelabel' will be assumed to be
#'      sample names.
#' @param Resolution For Exactive, the Resolution is 100000, defined at Mw 200
#' @param ResDefAt Resolution defined at (in Mw), e.g. 200 Mw
#' @param C13Purity Carbon 13 purity, default: 0.99
#' @param ReportPoolSize default: TRUE
#' @importFrom rlang .data
#' @return Named list of matrices: 'Corrected', 'Normalized',
#'      'PoolBeforeDF', and 'PoolAfterDF'.
#' @export
#' @examples
#' \dontrun{
#' carbon_correction("ExcelFile", Resolution=100000, ResDefAt=200)
#' }
carbon_correction <- function(path, sheet = NULL, output_path = NULL, ColumnsToSkip = NULL,
                              Resolution, ResDefAt, C13Purity = 0.99, ReportPoolSize = TRUE) {

  input_data <- read_elmaven_xlsx(path = path, sheet = sheet)
  sample_col_names <- names(input_data$cleaned)[which(!(tolower(names(input_data$cleaned)) %in%
                                               tolower(c("compound", "formula", "isotope_label", "label_index"))))]

  # Setup empty matrices for output
  # TODO Refactor this to preallocate or better yet, use sapply
  OutputMatrix <- matrix(NA, nrow = 0, ncol = length(sample_col_names))
  OutputPercentageMatrix <- matrix(NA, nrow = 0, ncol = length(sample_col_names))
  OutputPoolBefore <- matrix(NA, nrow = 0, ncol = length(sample_col_names))
  OutputPoolAfter <- matrix(NA, nrow = 0, ncol = length(sample_col_names))
  OutputCompound <- NULL
  OutputLabel <- NULL
  OutputPoolCompound <- NULL

  for (i in unique(input_data$cleaned$compound)) {
    CurrentMetabolite <- dplyr::filter(input_data$cleaned, .data$compound==i)
    Formula=as.character(CurrentMetabolite$formula[1])
    if(length(Formula)==0 || is.na(Formula)) {
      print(paste("The formula of",i,"is unknown",sep=" "))
      break
    }
    DataMatrix <- data.matrix(dplyr::select(CurrentMetabolite, -.data$compound, -.data$formula,
                                            -.data$isotope_label, -.data$label_index)
    )
    DataMatrix[is.na(DataMatrix)] <- 0
    Corrected <- IsotopeCorrection(Formula, DataMatrix, CurrentMetabolite$label_index,
                                   Resolution = Resolution, ResDefAt = ResDefAt,
                                   C13Purity = C13Purity, ReportPoolSize = ReportPoolSize)
    CorrectedPercentage <- scale(Corrected,scale=colSums(Corrected),center=FALSE)
    OutputMatrix <- rbind(OutputMatrix, Corrected)
    OutputPercentageMatrix <- rbind(OutputPercentageMatrix, CorrectedPercentage)
    OutputPoolBefore <- rbind(OutputPoolBefore, colSums(DataMatrix))
    OutputPoolAfter <- rbind(OutputPoolAfter, colSums(Corrected))
    OutputCompound <- append(OutputCompound, rep(i,nrow(Corrected)))
    OutputLabel <- append(OutputLabel, c(0:(nrow(Corrected)-1)))
    OutputPoolCompound <- append(OutputPoolCompound, i)
  }

  OutputDF <- data.frame(OutputCompound, OutputLabel, OutputMatrix)
  OutputPercentageDF <- data.frame(OutputCompound, OutputLabel, OutputPercentageMatrix)
  OutputPoolBeforeDF <- data.frame(OutputPoolCompound, OutputPoolBefore)
  OutputPoolAfterDF <- data.frame(OutputPoolCompound, OutputPoolAfter)
  names(OutputDF) <- c("Compound", "Label", sample_col_names)
  names(OutputPercentageDF) <- names(OutputDF)
  names(OutputPoolBeforeDF) <- c("Compound", sample_col_names)
  names(OutputPoolAfterDF) <- c("Compound", sample_col_names)

  OutputDataFrames <- list("Original" = input_data$original,
                           "Corrected" = OutputDF,
                           "Normalized" = OutputPercentageDF,
                           "PoolBeforeDF" = OutputPoolBeforeDF,
                           "PoolAfterDF" = OutputPoolAfterDF)

  if(!identical(FALSE, output_path)) {
    if(is.null(output_path)) {
      output_path = paste(tools::file_path_sans_ext(path), "_corrected.", tools::file_ext(path), sep="")
    }
    writexl::write_xlsx(OutputDataFrames, output_path)
  }

  return(OutputDataFrames)
}
