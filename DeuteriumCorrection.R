###Please note that this code only deals with carbon labeling experiment.
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
H2Purity <- 0.99
Resolution <- 140000
ResDefAt <- 200
ReportPoolSize <- TRUE
#For Exactive, the Resolution is 100000, defined at Mw 200#

#Please make sure the paths below are correct. 
#Make sure to use / in the paths

InputFile <- "C:/Users/lc8/Desktop/0414RtNADPH.xlsx"
InputSheetName <- "Sheet1"

InputDF <- read.xlsx(InputFile, InputSheetName, header = TRUE, check.names=FALSE, stringsAsFactors=FALSE)
#InputDF <- InputTable %>% filter(isotopeLabel=="C12 PARENT") %>% 
#  mutate(CompoundLabel=paste(compound," (",round(medRt,2),"min)",sep=""))
#InputDF <- InputDF[,c(ncol(InputDF,))]

InputDF[,1] <- as.character(InputDF[,1])
OutputMatrix <- matrix(0, ncol=(ncol(InputDF)-14), nrow=0)
OutputPercentageMatrix <- matrix(0, ncol=(ncol(InputDF)-14), nrow=0)
OutputPoolBefore <- matrix(0, ncol=(ncol(InputDF)-14), nrow=0)
OutputPoolAfter <- matrix(0, ncol=(ncol(InputDF)-14), nrow=0)
OutputCompound <- NULL
OutputLabel <- NULL
OutputPoolCompound <- NULL
names(InputDF)[9] <- "Compound"
InputDF$Label <- rep(NA,nrow(InputDF))

if.not.null <- function(x) if(!is.null(x)) x else 0

for (i in 1:nrow(InputDF)) {
  if(startsWith(InputDF[i,8], "C12")) {
    InputDF$Label[i]=0
    InputDF$Compound[i]=paste(InputDF[i,10]," (",round(InputDF[i,6],2),"min)",sep="")
  }
  else if(startsWith(InputDF[i,8], "D2-label-")) {
    InputDF$Label[i] <- unlist(strapply(InputDF[i,8], "(-)(\\d*)", ~ as.numeric(if (nchar(..2)) ..2)))
    InputDF$Compound[i] <- InputDF$Compound[i-1]  }
}

IsotopeCorrection <- function(formula, datamatrix, label) {
  
  AtomNumber <- rep(0,6)
  names(AtomNumber) <- c("C","H","N","O","P","S")
  MassDifference <- abs(c((13.00335-12),(15.00011-14.00307),(16.99913-15.99491),
                          (17.99916-15.99491),(32.97146-31.97207),(33.96787-31.97207))-
                          ((2.0141-1.00783)*c(1,1,1,2,1,2)))
  names(MassDifference) <- c("C13","N15","O17","O18","S33","S34")
  CorrectionLimit <- rep(0,6)
  names(CorrectionLimit) <- c("C13","N15","O17","O18","S33","S34")
  
  AtomNumber["C"] <- if.not.null(unlist(strapply(formula, "(C)(\\d*)", ~ as.numeric(if (nchar(..2)) ..2 else 1))))
  AtomNumber["H"] <- if.not.null(unlist(strapply(formula, "(H)(\\d*)", ~ as.numeric(if (nchar(..2)) ..2 else 1))))
  AtomNumber["N"] <- if.not.null(unlist(strapply(formula, "(N)(\\d*)", ~ as.numeric(if (nchar(..2)) ..2 else 1))))
  AtomNumber["O"] <- if.not.null(unlist(strapply(formula, "(O)(\\d*)", ~ as.numeric(if (nchar(..2)) ..2 else 1))))
  AtomNumber["P"] <- if.not.null(unlist(strapply(formula, "(P)(\\d*)", ~ as.numeric(if (nchar(..2)) ..2 else 1))))
  AtomNumber["S"] <- if.not.null(unlist(strapply(formula, "(S)(\\d*)", ~ as.numeric(if (nchar(..2)) ..2 else 1))))
  
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
    PurityMatrix[i,] <- sapply(0:AtomNumber["H"], function(x) dbinom(x-i+1, x , (1-H2Purity)))
  }
  
  for(j in 0:min(AtomNumber["C"], CorrectionLimit["C13"], AtomNumber["H"]))
    for(i in 1:(AtomNumber["H"]-j+1)){
      CarbonMatrix[(i+j),i] <- dbinom(j, AtomNumber["C"], CarbonNaturalAbundace[2])
    }
  
  for(i in 1:(AtomNumber["H"]+1)){
    HydrogenMatrix[,i] <- sapply(0:AtomNumber["H"], function(x) dbinom(x-i+1, AtomNumber["H"]-i+1 , HydrogenNaturalAbundace[2]))
  }
  
  for(j in 0:min(AtomNumber["N"], CorrectionLimit["N15"], AtomNumber["H"]))
    for(i in 1:(AtomNumber["H"]-j+1)){
      NitrogenMatrix[(i+j),i] <- dbinom(j, AtomNumber["N"], NitrogenNaturalAbundace[2])
    }
  
  for(i in 0:min(AtomNumber["O"],CorrectionLimit["O17"])) {
    for(j in 0:min(AtomNumber["O"],CorrectionLimit["O18"])){
      k<-(i+j*2)
      if ((i+j)>AtomNumber["O"]|k>AtomNumber["H"]) {
        break
      }
      else {
        for (m in 1:(AtomNumber["H"]-k+1)) {
          OxygenMatrix[(m+k),m] <- OxygenMatrix[(m+k),m] + dmultinom(c((AtomNumber["O"]-i-j),i,j), AtomNumber["O"], OxygenNaturalAbundace)      
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
          SulfurMatrix[(m+k),m] <- SulfurMatrix[(m+k),m] + dmultinom(c((AtomNumber["S"]-i-j),i,j), AtomNumber["S"], SulfurNaturalAbundace)      
        }
      }
    }
  }
  
  for(i in 1:ncol(datamatrix)) {
    CorrectedMatrix[,i] <- coef(nnls(SulfurMatrix %*% CarbonMatrix %*% OxygenMatrix %*% 
                                       NitrogenMatrix %*% HydrogenMatrix %*% PurityMatrix, ExpMatrix[,i]))
  }
  
  return(CorrectedMatrix)
  
}
i=unique(InputDF$Compound)[2]
for (i in unique(InputDF$Compound)) {
  CurrentMetabolite <- filter(InputDF, Compound==i)
  Formula=CurrentMetabolite[1,11]
  if(length(Formula)==0) {
    print(paste("The formula of",i,"is unknown",sep=" "))
    break
  }

  DataMatrix <- data.matrix(CurrentMetabolite[,15:(ncol(CurrentMetabolite)-1)])
  DataMatrix[is.na(DataMatrix)] <- 0
  Corrected <- IsotopeCorrection(Formula, DataMatrix, CurrentMetabolite$Label)
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
names(OutputDF) <- c("Compound", "Label", names(InputDF)[15:(length(names(InputDF))-1)])
names(OutputPercentageDF) <- names(OutputDF)
names(OutputPoolBeforeDF) <- c("Compound", names(InputDF)[15:(length(names(InputDF))-1)])
names(OutputPoolAfterDF) <- c("Compound", names(InputDF)[15:(length(names(InputDF))-1)])

write.xlsx2(OutputDF, file=InputFile, sheetName = "Corrected", row.names=FALSE, append=TRUE)
write.xlsx2(OutputPercentageDF, file=InputFile, sheetName = "Normalized", row.names=FALSE, append=TRUE)
write.xlsx2(OutputPoolAfterDF, file=InputFile, sheetName = "Pool Size", row.names=FALSE, append=TRUE)
