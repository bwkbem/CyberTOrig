

Lambda <- 11;

BioSamples <- 4;

MinProbes <- 2;
MinArrays <- 8;

#PARAMETERS FOR LOWESS
LowessfParam <- 0.4;

#PARAMETERS FOR HUBER/MAD ESTIMATOR
HuberConf <- 1.345; #(ie 95% Confidence Limit)
HuberTol <- 1e-6;


            ################# INPUT / OUTPUT  #################

PlotResults <- TRUE;

InputDirectory <- "/home/kirkb/Playarea/Bayesian Clustering/RockWeb/Linux/Data/Lowess/";

ProbesetFiles <- c("S730.LowPset", "F730.LowPset",
                   "S806.LowPset", "F806.LowPset",
                   "S812.LowPset", "F812.LowPset",
                   "S924.LowPset", "F924.LowPset");
HeaderLength <- 9;  #DO NOT EDIT

GenomeDirectory <- "/home/kirkb/Playarea/Bayesian Clustering/Distribution Depend/Data/Input/GenomeCrawlerSuite/"

GenomeFile <- "SPyM1Table.txt";
GenomeIDPos <- 6;   #DO NOT EDIT
GenomeHeaderLength <- 0;   #DO NOT EDIT



OutputDirectory <- "/home/kirkb/Playarea/Bayesian Clustering/Presentations/Thinking/temp/";
OutputFile <- "G8bA8P2tStat";

#This software is a rerelease of CyberTOR_WY.R module of Microarray Expression
#Potential Suite with alterations in the input and ouput format in order to
#improve transition between modules of eMESS.

#CyberT algorithm from:
#CyberT for use with spotted microaray data based on the derivation in
#"Improved Statistical Inference from DNA Microarray Data Using Analysis of
#Variance and Bayesian Statistical Framework" A. D. Long et al, J. Bio. Chem.,
#276(22), pp 19937-19944, 2001.
#Resampling-based FDR control of Benjamini-Hochberg (third alternative, p371)
#for use with microarray data based on the derivation in "Identifying
#differentially expressed genes using false discovery rate controlling
#procedures", A. Reiner et al, Bioinformatics, 19(3), pp 368-375, 2003.


#No part of this text or software may be used or reproduced in any matter
#whatsoever without written permission, except in the case of brief quotations
#embodied in critical articles of or for reviews.  For information address
#Brian W. Kirk (bwkbem@yahoo.com).
#Software Copyright (C) October 2002  Brian W. Kirk
#Form/Text Copyright (C) October 2002  Brian W. Kirk
#Released for use with R v1.5.1 (C)

#Copyright (C) October 2002  Brian W. Kirk
#Released for R v1.5.1 (C)

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  

###############################################################################
#                                FUNCTIONS                                    #
###############################################################################

library(MASS);

#FUNCTION THAT CALCULATES FOR A VECTOR OF SIGNALS A ROBUST SUMMARY USING THE
#HUBER M-ESTIMATOR AND THE MAD ESTIMATOR FOR SCALE.  RESULTS ARE COMBINED INTO
#A DATA.FRAME, WITH INCOMPATIBLE SETS (IE WITH ONLY 1 ENTRY) AMENDED TO THE END
#OF THE DATA.FRAME.

  fhuber.summarize <- function(SignalVec=vec, Filter=filter,
                               NameIdentifier=name, HuberConf=k, HuberTol=tol)
  {
   #Prime Header of Output Table
   huber.report.template <- c("ProbeID", "mu", "s", "probes");
   
   #Determine the Number of Probes Present for each Probeset 
   probes <- NameIdentifier[Filter];
   probesetf <- factor(probes);
   probesets <- levels(probesetf);

   probe.number <- as.vector(tapply(probes, probesetf, length));
   
   #probe.number.table <- data.frame(cbind(probesets, probe.number))

   #Select Porbesets With More Than 1 Probe For Huber Analysis
   huber.ok <- probesets[which(as.numeric(
                                     as.vector(probe.number)) > 1)];

   #Apply Estimators To Qualified Probesets
   if (length(huber.ok) > 0) {
     for (i in (1:length(huber.ok))) {

       #Test That More Than Half of the Probes Are Not Identical
       test.probe <- SignalVec[as.logical(
                        match(NameIdentifier, huber.ok[i], nomatch=0)*Filter)];
     
       test.probef <- factor(test.probe);
       test.probes <- levels(test.probef);
       level.lengths <- tapply(test.probe, test.probef, length);

       continue <- TRUE;
       for (j in 1:length(test.probes)) {
         if ((length(test.probe) - level.lengths[j]) < 0.5*length(test.probe)){
           continue <- FALSE
         }
       }

       #Process Accordingly Based on the Previous Test
       if (continue) {
         xtemp.matrix <-  as.matrix(huber(test.probe, k=HuberConf,
                                          tol=HuberTol));
         huber.report.template <- cbind(huber.report.template,
                        c(as.character(huber.ok[i]),
                          as.numeric(xtemp.matrix[1,1]),
                          as.numeric(xtemp.matrix[2,1]),
                          as.numeric(as.vector(
                              probe.number[as.logical(match(
                              probesets, huber.ok[i], nomatch=0))]))));
       } else {
         mu <- median(test.probe);
         scale <- sd(test.probe);

         huber.report.template <- cbind(huber.report.template,
                        c(as.character(huber.ok[i]),
                          mu,
                          scale,
                          as.numeric(as.vector(
                              probe.number[as.logical(match(
                              probesets, huber.ok[i], nomatch=0))]))));
   
       } 
     }
   }
   
   #Organize Porbesets That Did Not Qualify For Estimator Analysis
   huber.lost <- probesets[which(as.numeric(
                                    as.vector(probe.number)) <= 1)];
   if (length(huber.lost) > 0) {
     for (i in (1:length(huber.lost))) {
       huber.report.template <- cbind(huber.report.template,
                            c(as.character(huber.lost[i]),
                              as.numeric(SignalVec
                                          [as.logical(match(NameIdentifier,
                                           huber.lost[i], nomatch=0)
                                           * Filter)]), NA, 1))
     }
   }

   #Combine Both Probeset Types Into A Single Report
   huber.report.template2 <- t(huber.report.template);
   huber.report.template3 <-
                        huber.report.template2[2:(length(probesets)+1),];
   huber.sort.filter <- order(as.character(huber.report.template3[,1]));
   
   return(data.frame(cbind(
                  ProbeID=
                    as.character(huber.report.template3[huber.sort.filter, 1]),
                  mu=
                    as.numeric(huber.report.template3[huber.sort.filter, 2]),
                  s=
                    as.numeric(huber.report.template3[huber.sort.filter, 3]),
                  probes=
                    as.numeric(
                    as.vector(huber.report.template3[huber.sort.filter, 4])))))
  }


###############################################################################
#                           CONSTANT DECLARATION                              #
###############################################################################

#CONSTANTS AND VARIABLES OF TEMPORARY USE BEGIN WITH "x, y, or z" AND ALL
#USER DEFINED FUNCTIONS BEGIN WITH "f" or "g" IN ORDER TO FACILITATE REMOVABLE
#FROM MEMORY.  WHEN USING VARIABLES WITH THESE RESERVED INTIAL LETTERS,
#NEAREST CLEANUP SHOULD BE IDENTIFIED TO INSURE PROPER PERMANENCY

                 ############################################
                 #######  DEFINE PERMANENT CONSTANTS  #######
                 ############################################

#Define Output Dirctory with Prefix
OutputPrefix <- paste(OutputDirectory, OutputFile, sep="");

#DEFINE GRAPHICS WINDOW SIZE
WinWidth <- 8.5;           
WinHeight <- 11;


             #####################################################
             ##########  DEFINE TEMPORARY CONSTANTS  #############
             #####################################################


#Column position in LNProbesetDatawBGSubtract.RData
xIDPos <- 1;
MuFoldPos <- 2;
xMuAmpPos <- 3;
xProbeNumberPos <- 4;


###############################################################################
#                             VARIABLE DECLARATION                            #
###############################################################################

#CONSTANTS AND VARIABLES OF TEMPORARY USE BEGIN WITH "x", "y", or "z" AND ALL
#USER DEFINED FUNCTIONS BEGIN WITH "f" or "g" IN ORDER TO FACILITATE REMOVABLE
#FROM MEMORY.  WHEN USING VARIABLES WITH THESE RESERVED INTIAL LETTERS,
#NEAREST CLEANUP SHOULD BE IDENTIFIED TO INSURE PROPER PERMANCY



#Load LNProbesetData for each array into a list of data.frames
TotalArrayNumber <- length(ProbesetFiles);
ArrayList <- vector("list", TotalArrayNumber);
m <- 0;
for (Array in ProbesetFiles) {
  m <- m + 1;
  zNormalizedData <- paste(InputDirectory, Array, sep="");
  ArrayList[[m]] <- data.frame(read.delim(zNormalizedData, skip=HeaderLength));
}

#Concatenate respective vectors from each element of the ArrayList
IDVec <- as.character(as.vector(
                             ArrayList[[1]][,xIDPos]));
LogFoldVec <- as.numeric(as.vector(
                          ArrayList[[1]][,MuFoldPos]));
zRMSIntensityVec <- as.numeric(as.vector(
                               ArrayList[[1]][,xMuAmpPos]));
zProbeNumberVec  <- as.numeric(as.vector(
                               ArrayList[[1]][,xProbeNumberPos]));

ArrayID <- rep(1, length(IDVec));

for (n in 2:TotalArrayNumber) {
  
  tempIDVec <- as.character(as.vector(
                                 ArrayList[[n]][,xIDPos]));
  
  IDVec <- c(IDVec, tempIDVec);
  
  LogFoldVec <- c(LogFoldVec,
                   as.numeric(as.vector(
                              ArrayList[[n]][,MuFoldPos])));
  
  zRMSIntensityVec <- c(zRMSIntensityVec,
                        as.numeric(as.vector(
                                   ArrayList[[n]][,xMuAmpPos])));
  
  zProbeNumberVec  <- c(zProbeNumberVec,
                       as.numeric(as.vector(
                                  ArrayList[[n]][,xProbeNumberPos])));
  
  ArrayID <- c(ArrayID, rep(n, length(tempIDVec)));
}

#Filter for data that passes the minimum number of probes required
ProbeFilter <- as.logical(zProbeNumberVec >= MinProbes);

IDVec <- IDVec[ProbeFilter];
LogFoldVec <- LogFoldVec[ProbeFilter];
zRMSIntensityVec <- zRMSIntensityVec[ProbeFilter];
zProbeNumberVec <- zProbeNumberVec[ProbeFilter];
ArrayID <- ArrayID[ProbeFilter];
  

#Input Ordered Genes from Annotated Genome
GenomeDataIn <- paste(GenomeDirectory, GenomeFile, sep="");
xGenomeData <- data.frame(read.delim(GenomeDataIn,
                                       skip=GenomeHeaderLength));
GenomeGenes <- as.character(as.vector(xGenomeData[,GenomeIDPos]));

rm(list=ls(pat="^x"));
###############################################################################
#                               MAIN BODY                                     #
###############################################################################

                 #############################################
                 ####     Apply huberM/MAD Estimators     #### 
                 #############################################

DummyFilter <- rep(TRUE, length(IDVec));

zLogFoldHuberReport <- fhuber.summarize(LogFoldVec, DummyFilter,
                                IDVec, HuberConf, HuberTol);

zRMSHuberReport <- fhuber.summarize(zRMSIntensityVec, DummyFilter,
                                IDVec, HuberConf, HuberTol);


#Match Huber Reports for Intensity data and LogFold Change Data
#(ie if gene is present in only one, then remove it)
zRMSHuberSetFilter <- !as.logical(
                             match(zRMSHuberReport[,3], c(NA, 0), nomatch=0) |
                             as.numeric(as.vector(zRMSHuberReport[,4])) <
                                        MinArrays);

zLogFoldHuberSetFilter <- !as.logical(
                         match(zLogFoldHuberReport[,3], c(NA, 0), nomatch=0) |
                         as.numeric(as.vector(zLogFoldHuberReport[,4])) <
                                    MinArrays);

zTotHuberSetFilter <- as.logical(zRMSHuberSetFilter * zLogFoldHuberSetFilter);

HuberGenes <- as.character(zLogFoldHuberReport[zTotHuberSetFilter, 1]);
            
Intensity <- as.numeric(as.vector(zRMSHuberReport[zTotHuberSetFilter, 2]));

LogFoldChange <- as.numeric(as.vector(
                                  zLogFoldHuberReport[zTotHuberSetFilter, 2]));

LogFoldStndDev <- as.numeric(as.vector(
                                  zLogFoldHuberReport[zTotHuberSetFilter, 3]));

ArrayNumber <- as.numeric(as.vector(zRMSHuberReport[zTotHuberSetFilter, 4]));

StndErr <- LogFoldStndDev/sqrt(BioSamples);

### Remove non Gene Probes and unrepresented genes in annotation
GenomeFilter <- as.logical(as.vector(
                             match(HuberGenes, GenomeGenes, nomatch=0)));

HuberGenes <- HuberGenes[GenomeFilter];
            
Intensity <- Intensity[GenomeFilter];

LogFoldChange <- LogFoldChange[GenomeFilter];

LogFoldStndDev <- LogFoldStndDev[GenomeFilter];

ArrayNumber <- ArrayNumber[GenomeFilter];

StndErr <- StndErr[GenomeFilter];


rm(list=ls(pat="^z"));
rm(list=ls(pat="^f"));
rm(list=ls(pat="^g"));
      
                 ##############################################
                 ####    CALCULATE ADJ STANDRD DEVIATION   #### 
                 ##############################################

# LOWESS fit data and organize for extraction of LOWESS values
xBGSigmaLowess <- lowess(Intensity, LogFoldStndDev, f=LowessfParam);

BGStDev <- 0;
xLowessFilter <- order(Intensity);
for (j in 1:length(xLowessFilter)) {
  tempBGStDev <- xBGSigmaLowess$y[as.logical(
                                      match(xLowessFilter, j, nomatch=0))]
  BGStDev <- c(BGStDev, tempBGStDev);
}
                                
BGStDev <- BGStDev[2:(length(BGStDev))]

#Calculate Adjusted Stndard Err
#nu is defined to result in a constant lambda for all genes
nu <- 0;
xAdjustedVar <- 0;
for (k in 1:length(HuberGenes)) {
  if (BioSamples >= Lambda) {
    tempnu <- 0;
    } else {
      tempnu <- Lambda - BioSamples;
    }
  xSigmaSquared <- (tempnu * BGStDev[k]^2 +
       (BioSamples-1) * LogFoldStndDev[k]^2)/
           (tempnu + BioSamples-2);
    
  xAdjustedVar <- c(xAdjustedVar, xSigmaSquared);
  nu <- c(nu, tempnu);
}
      
xAdjustedVar <- xAdjustedVar[2:(length(xAdjustedVar))];
nu <- nu[2:length(nu)];
      
AdjustedStndDev <- sqrt(xAdjustedVar);
AdjustedStndErr <- sqrt(xAdjustedVar/BioSamples);

rm(list=ls(pat="^temp"));
      
              #### Plot Fits Before and After Adjustment  ####
if (PlotResults) {
  
  #Plot to file
  postscript(file=paste(OutputPrefix, "CyberT", "L",
             as.character(Lambda), ".eps", sep=""),  horizontal=FALSE,
             onefile=TRUE);

  chh <-par()$cxy[2];
  chw <-par()$cxy[2];
  par(mar=c(4,5,3,2), omi=c(0.5,0.5,0.5,0.5));

  par(mfrow=c(2,1));

  y0 <- as.numeric(as.vector(LogFoldStndDev));
  ymax <- max(y0);
  x0 <- as.numeric(as.vector(Intensity));
 
  plot(y0 ~ x0, type="p",  main="StndDev", 
     xlab="A", ylab = "StndDev \nLogFoldProbeSets",
     cex= 1.0, xlim=c(0, 15),
     ylim=c(0, ymax),
     pch=16);
  lines(xBGSigmaLowess$x, xBGSigmaLowess$y, col=2, lwd=3);

  y0 <- AdjustedStndDev;
  x0 <- as.numeric(as.vector(Intensity));

  #Best fit after adjustment
  xCorrectedLowess <- lowess(Intensity, AdjustedStndDev, f=LowessfParam);
  
  plot(y0 ~ x0, type="p",  main="Adj Stnd Dev", 
     xlab="A", ylab = "AdjStndDev \nLogFoldProbeSets",
     cex= 1.0, xlim=c(0, 15),
     ylim=c(0, ymax),
     pch=16);
  lines(xBGSigmaLowess$x, xBGSigmaLowess$y, col=2, lwd=3);
  lines(xCorrectedLowess$x, xCorrectedLowess$y, col=3, lwd=3);
 
  dev.off()
}
rm(list=ls(pat="^x"));
rm(list=ls(pat="^y"));


           ###################################################
           ####    PERFORM t-test AND CALCULATE pVALUES   #### 
           ###################################################

#Calculate t-Test statistic using a pairwise t-Test with log fold data
tStatNull <- abs(LogFoldChange/StndErr);
tStat <- abs(LogFoldChange/AdjustedStndErr);

#Calculate pValue with Bonferroni Correction for comparison
tStatFilterNull <- order(tStatNull, decreasing=TRUE);
RankedtStatNull <- tStatNull[tStatFilterNull];

#Calculate pValue with Bonferroni Correction for comparison
tStatFilter <- order(tStat, decreasing=TRUE);
RankedtStat <- tStat[tStatFilter];


xNumOfHyp <- length(RankedtStat);
Rank <- 1:xNumOfHyp;

FreqpValue <- Rank/xNumOfHyp;

DegreesOfFreedom <- Lambda - 2;
DegreesOfFreedomNull <- BioSamples - 1;

pValue <- 2 * (1 - pt(RankedtStat, DegreesOfFreedom));
pValueNull <- 2 * (1 - pt(RankedtStatNull, DegreesOfFreedomNull));


BHpValue <- pValue * xNumOfHyp / Rank;
BHpValue[which(BHpValue > 1)] <- 1;
for (i in Rank) {
  BHpValue[i] <- min(BHpValue[i:xNumOfHyp]);   
}

PKu <- 1-(2 * (1 - pt(sum(RankedtStat), DegreesOfFreedom)));
Pg <- 1 - pValue;
PhenompValue <- PKu - Pg;

PKuNull <- 1-(2 * (1 - pt(sum(RankedtStatNull), DegreesOfFreedomNull)));
PgNull <- 1 - pValueNull;
PhenompValueNull <- PKuNull - PgNull;


BonpValue <- xNumOfHyp * pValue;
BonpValue[which(BonpValue > 1)] <- 1;

SidakpValue <- 1 - Pg^xNumOfHyp;

HolmpValue <- (xNumOfHyp - Rank + 1) * pValue;
HolmpValue[which(HolmpValue > 1)] <- 1;
for (i in xNumOfHyp:1) {
  HolmpValue[i] <- max(HolmpValue[1:i]);   
}

HochpValue <- (xNumOfHyp - Rank + 1) * pValue;
HochpValue[which(HochpValue > 1)] <- 1;
for (i in Rank) {
  HochpValue[i] <- min(HochpValue[i:xNumOfHyp]);   
}

TRankedtStat <-  RankedtStat;
TBHpValue <- BHpValue;
TFreqpValue <- FreqpValue;
TPhenompValue <- PhenompValue;
TSidakpValue <- SidakpValue;


rm(ArrayList)
rm(list=ls(pat="^x"));
rm(list=ls(pat="^temp"));
rm(list=ls(pat="^Per"));
rm(list=ls(pat="^Next"));



postscript(file=paste("/home/kirkb/Playarea/Bayesian Clustering/Presentations/Thinking/temp/logpvslogTMicrodof3L14.eps", sep=""),  horizontal=TRUE, onefile=TRUE);

 
x11(width=11, height=8);


y <- log10(TSidakpValue);
x <- log10(TRankedtStat);

plot(y~x, type="p",
     xlab="", ylab="", xlim=c(min(x), max(x)),
                       ylim=c(min(log10(TPhenompValue )),0),
     cex=2, pch=1, col="red");


y <- log10(TFreqpValue);
x <- log10(TRankedtStat);

points(x, y, cex=2, pch=1, col="blue");


y <- log10(TBHpValue);
x <- log10(TRankedtStat);

points(x, y, cex=2, pch=1, col="green");

y <- log10(PhenompValueNull);
x <- log10(RankedtStatNull);

points(x, y, cex=2, pch=1, col="black");


legend(-3.3, -4, c("Phenom", "BH", "Freq", "Sidak"), cex=2, pch=1, col=c("black", "green", "blue", "red"));

dev.off();


postscript(file=paste("/home/kirkb/Playarea/Bayesian Clustering/Presentations/Thinking/temp/pvsTMicrodof3L14.eps", sep=""),  horizontal=TRUE, onefile=TRUE);


#x11(width=11, height=8);

Column <- 12;

y <- TSidakpValue ;
x <- TRankedtStat ;

plot(y~x, type="p",
     xlab="", ylab="", xlim=c(min(x), max(x)),
                       ylim=c(min(TPhenompValue),
                              max(TPhenompValue)),
     cex=2, pch=1, col="red");


y <- TFreqpValue;
x <- TRankedtStat;

points(x, y, cex=2, pch=1, col="blue");


y <- TBHpValue;
x <- TRankedtStat;

points(x, y, cex=2, pch=1, col="green");

y <- TPhenompValue;
x <- TRankedtStat;

points(x, y, cex=2, pch=1, col="black");


legend(7, 1, c("Phenom", "BH", "Freq", "Sidak"), cex=2, pch=1, col=c("black", "green", "blue", "red"));

dev.off();


postscript(file=paste("/home/kirkb/Playarea/Bayesian Clustering/Presentations/Thinking/temp/pvsTMicrodof3", as.character(Lambda), ".eps", sep=""),  horizontal=TRUE, onefile=TRUE);


#x11(width=11, height=8);


y <- SidakpValue;
x <- RankedtStat;

plot(y~x, type="p",
     xlab="", ylab="", xlim=c(min(x), max(x)),
                       ylim=c(min(PhenompValue), max(y)),
     cex=2, pch=1, col="red");


y <- FreqpValue;
x <- RankedtStat;

points(x, y, cex=2, pch=1, col="blue");


y <- BHpValue;
x <- RankedtStat;

points(x, y, cex=2, pch=1, col="green");

y <- PhenompValue;
x <- RankedtStat;

points(x, y, cex=2, pch=1, col="black");


legend(12.6, 1, c("Phenom", "BH", "Freq", "Sidak"), cex=2, pch=1, col=c("black", "green", "blue", "red"));

dev.off()




 
                       ############################
                       #####  OUTPUT RESULTS  #####
                       ############################
      


#Combine vectors into a table for output
CyberTBHTable <- cbind(Rank, ID=HuberGenes[tStatFilter], pValue, PhenompValue,
                       FreqpValue,
                       EstpValue, BHpValueT, WYpValue, BHpValue, BonpValue,
                       SidakpValue, HolmpValue, HochpValue,
                       tStat=tStat[tStatFilter],
                       LogFoldChange=LogFoldChange[tStatFilter],
                       AdjStndErr=AdjustedStndErr[tStatFilter],
                       ArrayNumber=ArrayNumber[tStatFilter],
                       RMSIntensity=Intensity[tStatFilter], BStat, 
                       BGStDev=BGStDev[tStatFilter],
                       StDev=LogFoldStndDev[tStatFilter],
                       StndErr=StndErr[tStatFilter]);

FilterHeader1 <- paste("Input Directory = ", InputDirectory, sep="");
FilterHeader2 <- c("NormFiles = ", paste(as.character(ProbesetFiles),
                                             sep=""));
FilterHeader3 <- date();
FilterHeader4 <- paste("Lambda = ", as.character(Lambda), sep="");
FilterHeader5 <- paste("MinArrays = ", as.character(MinArrays), sep="");
FilterHeader6 <- paste("MinProbes = ", as.character(MinProbes), sep="");   
FilterHeader7 <- paste("LowessfParam = ", as.character(LowessfParam),
                         sep="");          
FilterHeader8 <- c(paste("HuberConf = ", as.character(HuberConf), sep=""),
                   paste("HuberTol = ", as.character(HuberTol), sep=""));
FilterHeader9 <- c("Algorithm = Benjamini-Hochberg",
                   paste("Permutation Number = ",
                       as.character(TotalPermutations), sep=""));
FilterHeaderTable <-  c("Rank", "ID", "pValue", "PhenompValue",
                        "FreqpValue", "EstpValue", "BHpValueT", "WYpValue",
                        "BHpValue", "BonpValue", "SidakpValue", "HolmpValue",
                        "HochpValue", "tStat", "LogFoldChng", "AdjStndErr",
                        "ArrayNumber", "RMSIntens", "BStat", "BGStDev",
                        "StDev", "StndErr");

cat(FilterHeader1, file=paste(OutputPrefix, "L",
                     as.character(Lambda), ".CyTBH", sep=""), sep="\n");
cat(FilterHeader2, file=paste(OutputPrefix, "L",
                     as.character(Lambda), ".CyTBH", sep=""), sep="\t",
                     append=TRUE);
cat(c(""), file=paste(OutputPrefix, "L", as.character(Lambda),
                     ".CyTBH", sep=""), sep="\n", append=TRUE);
cat(FilterHeader3, file=paste(OutputPrefix, "L",
                     as.character(Lambda), ".CyTBH", sep=""), sep="\t",
                     append=TRUE);
cat(c(""), file=paste(OutputPrefix, "L", as.character(Lambda),
                     ".CyTBH", sep=""), sep="\n", append=TRUE);
cat(FilterHeader4, file=paste(OutputPrefix, "L",
                     as.character(Lambda), ".CyTBH", sep=""), sep="\t",
                     append=TRUE);
cat(c(""), file=paste(OutputPrefix, "L", as.character(Lambda),
                     ".CyTBH", sep=""), sep="\n", append=TRUE);
cat(FilterHeader5, file=paste(OutputPrefix, "L",
                     as.character(Lambda), ".CyTBH", sep=""), sep="\t",
                     append=TRUE);
cat(c(""), file=paste(OutputPrefix, "L", as.character(Lambda),
                     ".CyTBH", sep=""), sep="\n", append=TRUE);
cat(FilterHeader6, file=paste(OutputPrefix, "L",
                     as.character(Lambda), ".CyTBH", sep=""), sep="\t",
                     append=TRUE);
cat(c(""), file=paste(OutputPrefix, "L", as.character(Lambda),
                     ".CyTBH", sep=""), sep="\n", append=TRUE);
cat(FilterHeader7, file=paste(OutputPrefix, "L",
                     as.character(Lambda), ".CyTBH", sep=""), sep="\t",
                     append=TRUE);
cat(c(""), file=paste(OutputPrefix, "L", as.character(Lambda),
                     ".CyTBH", sep=""), sep="\n", append=TRUE);
cat(FilterHeader8, file=paste(OutputPrefix, "L",
                     as.character(Lambda), ".CyTBH", sep=""), sep="\t",
                     append=TRUE);
cat(c(""), file=paste(OutputPrefix, "L", as.character(Lambda),
                     ".CyTBH", sep=""), sep="\n", append=TRUE);
cat(FilterHeader9, file=paste(OutputPrefix, "L",
                     as.character(Lambda), ".CyTBH", sep=""), sep="\t",
                     append=TRUE);
cat(c(""), file=paste(OutputPrefix, "L", as.character(Lambda),
                     ".CyTBH", sep=""), sep="\n", append=TRUE);
cat(FilterHeaderTable, file=paste(OutputPrefix, "L",
                         as.character(Lambda), ".CyTBH", sep=""), sep="\t",
                         append=TRUE);
cat(c(""), file=paste(OutputPrefix, "L", as.character(Lambda),
             ".CyTBH", sep=""), sep="\n", append=TRUE);

write.table(CyberTBHTable, file=paste(OutputPrefix, "L",
            as.character(Lambda), ".CyTBH", sep=""), quote=FALSE, sep="\t",
            col.names=FALSE, row.names=FALSE, append=TRUE);



##################
##  Single Genes
##################
RankedGeneNames <- HuberGenes[tStatFilter];
RankedArrayNumber <-  ArrayNumber[tStatFilter];

GenetStat <- rep(0, length(GenomeGenes));
GeneArrayNumber <- rep(0, length(GenomeGenes));
GeneBHpValue <- rep(1, length(GenomeGenes));
GenepValue <- rep(1, length(GenomeGenes));
GenepValueBon <- rep(1, length(GenomeGenes));
GenePos <- 1:length(GenomeGenes);
GeneNames <- GenomeGenes;


DegreesOfFreedom <- Lambda - 1;
#DegreesOfFreedom <- GeneArrayNumber - 1;

pValue <- 2 * (1 - pt(RankedtStat, DegreesOfFreedom));
pValue[which(pValue > 1)] <- 1;

pValueBon <- length(tStatFilter)*(2 * (1 - pt(RankedtStat, DegreesOfFreedom)));
pValueBon[which(pValueBon > 1)] <- 1;

#DegreesOfFreedom <- RankedArrayNumber - 1;
#pValueN <- length(tStatFilter)*(2 * (1 - pt(RankedtStat, DegreesOfFreedom)));
#pValueN[which(pValueN > 1)] <- 1;

#### Output

for (i in 1:length(RankedGeneNames)) {
  GeneFilter <- as.logical(match(GeneNames, RankedGeneNames[i], nomatch=0));

  GeneBHpValue[GeneFilter] <- BHpValue[i];
  GenepValue[GeneFilter] <- pValue[i];
  GenepValueBon[GeneFilter] <- pValueBon[i];
  GenetStat[GeneFilter] <- RankedtStat[i];
}



SingleGeneTable <- cbind(GeneNames, GenetStat, GenepValue, GeneBHpValue,
                         GenepValueBon, GenePos);

FilterHeader1 <- paste("Input Directory = ", InputDirectory, sep="");
FilterHeader2 <- c("NormFiles = ", paste(as.character(ProbesetFiles),
                                             sep=""));
FilterHeader3 <- date();
FilterHeader4 <- paste("Lambda = ", as.character(Lambda), sep="");
FilterHeader5 <- paste("MinArrays = ", as.character(MinArrays), sep="");
FilterHeader6 <- paste("MinProbes = ", as.character(MinProbes), sep="");   
FilterHeader7 <- paste("LowessfParam = ", as.character(LowessfParam),
                         sep="");          

FilterHeader8 <- c(paste("HuberConf = ", as.character(HuberConf), sep=""),
                   paste("HuberTol = ", as.character(HuberTol), sep=""));
FilterHeader9 <- c("Algorithm = Benjamini-Hochberg",
                   paste("Permutation Number = ",
                       as.character(TotalPermutations), sep=""));
FilterHeaderTable <-  c("ID", "tStat", "pValue","BHpValue", "pValue Bon",  "Gene Pos");

cat(FilterHeader1, file=paste(OutputPrefix, "L",
                     as.character(Lambda), ".GeCyTBH", sep=""), sep="\n");
cat(FilterHeader2, file=paste(OutputPrefix, "L",
                     as.character(Lambda), ".GeCyTBH", sep=""), sep="\t",
                     append=TRUE);
cat(c(""), file=paste(OutputPrefix, "L", as.character(Lambda),
                     ".GeCyTBH", sep=""), sep="\n", append=TRUE);
cat(FilterHeader3, file=paste(OutputPrefix, "L",
                     as.character(Lambda), ".GeCyTBH", sep=""), sep="\t",
                     append=TRUE);
cat(c(""), file=paste(OutputPrefix, "L", as.character(Lambda),
                     ".GeCyTBH", sep=""), sep="\n", append=TRUE);
cat(FilterHeader4, file=paste(OutputPrefix, "L",
                     as.character(Lambda), ".GeCyTBH", sep=""), sep="\t",
                     append=TRUE);
cat(c(""), file=paste(OutputPrefix, "L", as.character(Lambda),
                     ".GeCyTBH", sep=""), sep="\n", append=TRUE);
cat(FilterHeader5, file=paste(OutputPrefix, "L",
                     as.character(Lambda), ".GeCyTBH", sep=""), sep="\t",
                     append=TRUE);
cat(c(""), file=paste(OutputPrefix, "L", as.character(Lambda),
                     ".GeCyTBH", sep=""), sep="\n", append=TRUE);
cat(FilterHeader6, file=paste(OutputPrefix, "L",
                     as.character(Lambda), ".GeCyTBH", sep=""), sep="\t",
                     append=TRUE);
cat(c(""), file=paste(OutputPrefix, "L", as.character(Lambda),
                     ".GeCyTBH", sep=""), sep="\n", append=TRUE);
cat(FilterHeader7, file=paste(OutputPrefix, "L",
                     as.character(Lambda), ".GeCyTBH", sep=""), sep="\t",
                     append=TRUE);
cat(c(""), file=paste(OutputPrefix, "L", as.character(Lambda),
                     ".GeCyTBH", sep=""), sep="\n", append=TRUE);
cat(FilterHeader8, file=paste(OutputPrefix, "L",
                     as.character(Lambda), ".GeCyTBH", sep=""), sep="\t",
                     append=TRUE);
cat(c(""), file=paste(OutputPrefix, "L", as.character(Lambda),
                     ".GeCyTBH", sep=""), sep="\n", append=TRUE);
cat(FilterHeader9, file=paste(OutputPrefix, "L",
                     as.character(Lambda), ".GeCyTBH", sep=""), sep="\t",
                     append=TRUE);
cat(c(""), file=paste(OutputPrefix, "L", as.character(Lambda),
                     ".GeCyTBH", sep=""), sep="\n", append=TRUE);
cat(FilterHeaderTable, file=paste(OutputPrefix, "L",
                         as.character(Lambda), ".GeCyTBH", sep=""),
                         sep="\t", append=TRUE);
cat(c(""), file=paste(OutputPrefix, "L", as.character(Lambda),
             ".GeCyTBH", sep=""), sep="\n", append=TRUE);

write.table(SingleGeneTable, file=paste(OutputPrefix, "L",
            as.character(Lambda), ".GeCyTBH", sep=""), quote=FALSE,
            sep="\t", col.names=FALSE, row.names=FALSE, append=TRUE);

rm(SingleGeneTable);

#rm(list=ls());

}
