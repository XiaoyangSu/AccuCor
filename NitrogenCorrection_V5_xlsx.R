###Please note that this code only deals with nitrogen labeling experiment.
###Please use other version of correction codes for other experiments.
#Please make sure you have the packages below installed.### 
#Use install.packages("") to install packages###

require(gsubfn)
require(nnls)
require(dplyr)
require(gdata)
require(xlsx)

CarbonNaturalAbundace <- c(0.9893, 0.0107)
HydrogenNaturalAbundace <- c(0.999885, 0.000115)
NitrogenNaturalAbundace <- c(0.99636, 0.00364)
OxygenNaturalAbundace <- c(0.99757, 0.00038, 0.00205)
SulfurNaturalAbundace <- c(0.9493, 0.00762, 0.0429)

###Please make sure these parameters are accurate.#
N15Purity <- 0.99
Resolution <- 140000
ResDefAt <- 200
ReportPoolSize <- TRUE
#For Exactive, the Resolution is 100000, defined at Mw 200#

#Please make sure the paths below are correct. 
#Make sure to use / in the paths
#The input data file name should NOT contain the ".csv"

InputFile <- "D:/Sorcerers/Princeton/R/Correction Release/N_Sample_Input.xlsx"
InputSheetName <- "N_Input"
MetaboliteList <- read.csv("D:/Sorcerers/Princeton/R/Correction Release/KNOWNS.csv", header = TRUE, check.names = FALSE)

InputDF <- read.xlsx(InputFile, InputSheetName, header = FALSE, check.names=FALSE, stringsAsFactors=FALSE)
InputDF[,1] <- as.character(InputDF[,1])
OutputMatrix <- matrix(0, ncol=(ncol(InputDF)-1), nrow=0)
OutputPercentageMatrix <- matrix(0, ncol=(ncol(InputDF)-1), nrow=0)
OutputPoolBefore <- matrix(0, ncol=(ncol(InputDF)-1), nrow=0)
OutputPoolAfter <- matrix(0, ncol=(ncol(InputDF)-1), nrow=0)
OutputCompound <- NULL
OutputLabel <- NULL
OutputPoolCompound <- NULL
names(InputDF) <- sapply(InputDF[1,], as.character)
names(InputDF)[1] <- "Compound"
InputDF$Label <- rep(NA,nrow(InputDF))

if.not.null <- function(x) if(!is.null(x)) x else 0

for (i in 1:nrow(InputDF)) {
  if(startsWith(InputDF[i,1], "C12")) {
    InputDF$Label[i]=0
    InputDF[i,1] <- InputDF[i-1,1]
  }
  else if(startsWith(InputDF[i,1], "N15-label")) {
    InputDF$Label[i] <- unlist(strapply(InputDF[i,1], "(-)(\\d*)", ~ as.numeric(if (nchar(..2)) ..2)))
    InputDF[i,1] <- InputDF[i-1,1]  }
}

IsotopeCorrection <- function(formula, datamatrix, label) {
  
  AtomNumber <- rep(0,6)
  names(AtomNumber) <- c("C","H","N","O","P","S")
  MassDifference <- abs(c((13.00335-12),(2.0141-1.00783),(16.99913-15.99491),
                      (17.99916-15.99491),(32.97146-31.97207),(33.96787-31.97207))-
                      ((15.00011-14.00307)*c(1,1,1,2,1,2)))
  names(MassDifference) <- c("C13","H2","O17","O18","S33","S34")
  CorrectionLimit <- rep(0,6)
  names(CorrectionLimit) <- c("C13","H2","O17","O18","S33","S34")
  
  AtomNumber["C"] <- if.not.null(unlist(strapply(formula, "(C)(\\d*)", ~ as.numeric(if (nchar(..2)) ..2 else 1))))
  AtomNumber["H"] <- if.not.null(unlist(strapply(formula, "(H)(\\d*)", ~ as.numeric(if (nchar(..2)) ..2 else 1))))
  AtomNumber["N"] <- if.not.null(unlist(strapply(formula, "(N)(\\d*)", ~ as.numeric(if (nchar(..2)) ..2 else 1))))
  AtomNumber["O"] <- if.not.null(unlist(strapply(formula, "(O)(\\d*)", ~ as.numeric(if (nchar(..2)) ..2 else 1))))
  AtomNumber["P"] <- if.not.null(unlist(strapply(formula, "(P)(\\d*)", ~ as.numeric(if (nchar(..2)) ..2 else 1))))
  AtomNumber["S"] <- if.not.null(unlist(strapply(formula, "(S)(\\d*)", ~ as.numeric(if (nchar(..2)) ..2 else 1))))
  
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
    PurityMatrix[i,] <- sapply(0:(AtomNumber["N"]), function(x) dbinom(x-i+1, x , (1-N15Purity)))
  }
  
  for(i in 1:(AtomNumber["N"]+1)){
    NitrogenMatrix[,i] <- sapply(0:AtomNumber["N"], function(x) dbinom(x-i+1, AtomNumber["N"]-i+1 , NitrogenNaturalAbundace[2]))
  }
  
  for(j in 0:min(AtomNumber["C"], CorrectionLimit["C13"], AtomNumber["N"]))
    for(i in 1:(AtomNumber["N"]-j+1)){
      CarbonMatrix[(i+j),i] <- dbinom(j, AtomNumber["C"], CarbonNaturalAbundace[2])
    }
  
  for(j in 0:min(AtomNumber["H"], CorrectionLimit["H2"], AtomNumber["N"]))
    for(i in 1:(AtomNumber["N"]-j+1)){
      HydrogenMatrix[(i+j),i] <- dbinom(j, AtomNumber["H"], HydrogenNaturalAbundace[2])
    }
  
  for(i in 0:min(AtomNumber["O"],CorrectionLimit["O17"])) {
    for(j in 0:min(AtomNumber["O"],CorrectionLimit["O18"])){
      k<-(i+j*2)
      if ((i+j)>AtomNumber["O"]|k>AtomNumber["N"]) {
        break
      }
      else {
        for (m in 1:(AtomNumber["N"]-k+1)) {
          OxygenMatrix[(m+k),m] <- OxygenMatrix[(m+k),m] + dmultinom(c((AtomNumber["O"]-i-j),i,j), AtomNumber["O"], OxygenNaturalAbundace)      
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
          SulfurMatrix[(m+k),m] <- SulfurMatrix[(m+k),m] + dmultinom(c((AtomNumber["S"]-i-j),i,j), AtomNumber["S"], SulfurNaturalAbundace)      
        }
      }
    }
  }
  
  for(i in 1:ncol(datamatrix)) {
    CorrectedMatrix[,i] <- coef(nnls(SulfurMatrix %*% OxygenMatrix %*% CarbonMatrix %*% 
        HydrogenMatrix %*% NitrogenMatrix %*% PurityMatrix, ExpMatrix[,i]))
  }
  
  return(CorrectedMatrix)
  
}

for (i in unique(InputDF$Compound)) {
  Formula=as.character(MetaboliteList$formula[MetaboliteList$compound==i])
  if(length(Formula)==0) {
    print(paste("The formula of",i,"is unknown",sep=" "))
    break
  }
  CurrentMetabolite <- filter(InputDF, Compound==i)
  DataMatrix <- data.matrix(CurrentMetabolite[2:nrow(CurrentMetabolite),2:(ncol(CurrentMetabolite)-1)])
  DataMatrix[is.na(DataMatrix)] <- 0
  Corrected <- IsotopeCorrection(Formula, DataMatrix, CurrentMetabolite$Label[2:nrow(CurrentMetabolite)])
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
names(OutputDF) <- c("Compound", "Label", names(InputDF)[2:(length(names(InputDF))-1)])
names(OutputPercentageDF) <- names(OutputDF)
names(OutputPoolBeforeDF) <- c("Compound", names(InputDF)[2:(length(names(InputDF))-1)])
names(OutputPoolAfterDF) <- c("Compound", names(InputDF)[2:(length(names(InputDF))-1)])

write.xlsx2(OutputDF, file=InputFile, sheetName = "Corrected", row.names=FALSE, append=TRUE)
write.xlsx2(OutputPercentageDF, file=InputFile, sheetName = "Normalized", row.names=FALSE, append=TRUE)
write.xlsx2(OutputPoolAfterDF, file=InputFile, sheetName = "Pool Size", row.names=FALSE, append=TRUE)