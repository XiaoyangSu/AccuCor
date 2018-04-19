###Please note that this code only deals with carbon labeling experiment.
###Please use other version of correction codes for other experiments.
#Please make sure you have the packages below installed.###
#Use install.packages("") to install packages###

library(gsubfn)
library(nnls)
library(dplyr)
library(gdata)
library(xlsx)
library(stringr)

CarbonNaturalAbundace <- c(0.9893, 0.0107)
HydrogenNaturalAbundace <- c(0.999885, 0.000115)
NitrogenNaturalAbundace <- c(0.99636, 0.00364)
OxygenNaturalAbundace <- c(0.99757, 0.00038, 0.00205)
SulfurNaturalAbundace <- c(0.9493, 0.00762, 0.0429)

###Please make sure these parameters are accurate.#
C13Purity <- 0.99
Resolution <- 140000
ResDefAt <- 200
ReportPoolSize <- TRUE
#For Exactive, the Resolution is 100000, defined at Mw 200#

#Please make sure the paths below are correct.
#Make sure to use / in the paths

InputFile <- "test/JCGC_test_2.xlsx"
InputSheetName <- "JCGC_test"

InputDF <- read.xlsx(InputFile, InputSheetName, header = TRUE, check.names=FALSE, stringsAsFactors=FALSE)


IsotopeCorrection <- function(formula, datamatrix, label) {
  
  AtomNumber <- rep(0,6)
  names(AtomNumber) <- c("C","H","N","O","P","S")
  MassDifference <- abs(c((2.0141-1.00783),(15.00011-14.00307),(16.99913-15.99491),
                          (17.99916-15.99491),(32.97146-31.97207),(33.96787-31.97207))-
                          ((13.00335-12)*c(1,1,1,2,1,2)))
  names(MassDifference) <- c("H2","N15","O17","O18","S33","S34")
  CorrectionLimit <- rep(0,6)
  names(CorrectionLimit) <- c("H2","N15","O17","O18","S33","S34")
  
  if.not.null <- function(x) if(!is.null(x)) x else 0
  AtomNumber["C"] <- if.not.null(unlist(strapply(formula, "(C)(\\d*)", ~ as.numeric(if (nchar(..2)) ..2 else 1))))
  AtomNumber["H"] <- if.not.null(unlist(strapply(formula, "(H)(\\d*)", ~ as.numeric(if (nchar(..2)) ..2 else 1))))
  AtomNumber["N"] <- if.not.null(unlist(strapply(formula, "(N)(\\d*)", ~ as.numeric(if (nchar(..2)) ..2 else 1))))
  AtomNumber["O"] <- if.not.null(unlist(strapply(formula, "(O)(\\d*)", ~ as.numeric(if (nchar(..2)) ..2 else 1))))
  AtomNumber["P"] <- if.not.null(unlist(strapply(formula, "(P)(\\d*)", ~ as.numeric(if (nchar(..2)) ..2 else 1))))
  AtomNumber["S"] <- if.not.null(unlist(strapply(formula, "(S)(\\d*)", ~ as.numeric(if (nchar(..2)) ..2 else 1))))
  
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
    PurityMatrix[i,] <- sapply(0:(AtomNumber["C"]), function(x) dbinom(x-i+1, x , (1-C13Purity)))
  }
  
  for(i in 1:(AtomNumber["C"]+1)){
    CarbonMatrix[,i] <- sapply(0:AtomNumber["C"], function(x) dbinom(x-i+1, AtomNumber["C"]-i+1 , CarbonNaturalAbundace[2]))
  }
  
  for(j in 0:min(AtomNumber["N"], CorrectionLimit["N15"], AtomNumber["C"]))
    for(i in 1:(AtomNumber["C"]-j+1)){
      NitrogenMatrix[(i+j),i] <- dbinom(j, AtomNumber["N"], NitrogenNaturalAbundace[2])
    }
  
  for(j in 0:min(AtomNumber["H"], CorrectionLimit["H2"], AtomNumber["C"]))
    for(i in 1:(AtomNumber["C"]-j+1)){
      HydrogenMatrix[(i+j),i] <- dbinom(j, AtomNumber["H"], HydrogenNaturalAbundace[2])
    }
  
  for(i in 0:min(AtomNumber["O"],CorrectionLimit["O17"])) {
    for(j in 0:min(AtomNumber["O"],CorrectionLimit["O18"])){
      k<-(i+j*2)
      if ((i+j)>AtomNumber["O"]|k>AtomNumber["C"]) {
        break
      }
      else {
        for (m in 1:(AtomNumber["C"]-k+1)) {
          OxygenMatrix[(m+k),m] <- OxygenMatrix[(m+k),m] + dmultinom(c((AtomNumber["O"]-i-j),i,j), AtomNumber["O"], OxygenNaturalAbundace)
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
          SulfurMatrix[(m+k),m] <- SulfurMatrix[(m+k),m] + dmultinom(c((AtomNumber["S"]-i-j),i,j), AtomNumber["S"], SulfurNaturalAbundace)
        }
      }
    }
  }
  
  for(i in 1:ncol(datamatrix)) {
    CorrectedMatrix[,i] <- coef(nnls(SulfurMatrix %*% OxygenMatrix %*% NitrogenMatrix %*%
                                       HydrogenMatrix %*% CarbonMatrix %*% PurityMatrix, ExpMatrix[,i]))
  }
  
  return(CorrectedMatrix)
  
}


CarbonCorrection <- function(InputDF, ColumnsToSkip=NULL) {
  # Remove columns that are not needed
  if(is.null(ColumnsToSkip)) {
    ColumnsToSkip = c(
      "label", "metaGroupId", "groupId", "goodPeakCount", "medMz", "medRt", 
      "maxQuality", "isotopeLabel", "compoundId", "expectedRtDiff", "ppmDiff",
      "parent", "junk"
    )
  }
  keep_col_nums <- which(!(tolower(names(InputDF)) %in% tolower(ColumnsToSkip)))
  formula_col_num <- which(tolower(names(InputDF)) == 'formula')
  formula_col_name <- names(InputDF)[formula_col_num]
  compound_col_num <- which(tolower(names(InputDF)) == 'compound')
  compound_col_name <- names(InputDF)[compound_col_num]
  sample_col_names <- names(InputDF)[which(!(tolower(names(InputDF)) %in% tolower(c(ColumnsToSkip, "compound", "formula"))))] 
  # Generate label column as incremental index
  isotope_label_col_num <- which(tolower(names(InputDF)) == 'isotopelabel')
  isotope_label_col_name <- names(InputDF)[isotope_label_col_num]
  
  InputDF <- InputDF %>% 
    mutate(label_idx = str_replace(UQ(as.name(isotope_label_col_name)), "C13-label-", "")) %>% 
    mutate(label_idx = str_replace(label_idx, "C12 PARENT", "0")) %>%
    mutate(label_idx = as.numeric(label_idx))
  
  for (i in unique(InputDF[,compound_col_num])) {
    CurrentMetabolite <- filter(InputDF, UQ(as.name(compound_col_name))==i)
    Formula=CurrentMetabolite[1, formula_col_num]
    if(length(Formula)==0 || is.na(Formula)) {
      print(paste("The formula of",i,"is unknown",sep=" "))
      break
    }
    DataMatrix <- data.matrix(CurrentMetabolite %>%
                                select(keep_col_nums) %>%
                                select(-UQ(as.name(compound_col_name)),
                                       -UQ(as.name(formula_col_name)))
    )
    DataMatrix[is.na(DataMatrix)] <- 0
    Corrected <- IsotopeCorrection(Formula, DataMatrix, CurrentMetabolite$label_idx)
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
  
  return(list("Corrected" = OutputDF,
              "Normalized" = OutputPercentageDF,
              "PoolBeforeDF" = OutputPoolBeforeDF,
              "PoolAfterDF" = OutputPoolAfterDF))
}

OutputFile <- "test/JCGC_test_2_corrected.xlsx"
OutputDataFrames <- CarbonCorrection(InputDF)


write.xlsx2(OutputDataFrames$Corrected, file=OutputFile, sheetName = "Corrected", row.names=FALSE, append=TRUE)
write.xlsx2(OutputDataFrames$Normalized, file=OutputFile, sheetName = "Normalized", row.names=FALSE, append=TRUE)
write.xlsx2(OutputDataFrames$PoolAfterDF, file=OutputFile, sheetName = "Pool Size", row.names=FALSE, append=TRUE)
