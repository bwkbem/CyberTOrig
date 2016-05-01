

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

InputDirectory <- "/home/anon/Projects/Data/SPyM1/Lowess";

ProbesetFiles <- c("S730.LowPset", "F730.LowPset",
                   "S806.LowPset", "F806.LowPset",
                   "S812.LowPset", "F812.LowPset",
                   "S924.LowPset", "F924.LowPset");
HeaderLength <- 9;  #DO NOT EDIT

GenomeDirectory <- "/home/anon/Projects/Data/SPyM1"

GenomeFile <- "SPyM1Table.txt";
GenomeIDPos <- 6;   #DO NOT EDIT
GenomeHeaderLength <- 0;   #DO NOT EDIT



OutputDirectory <- "/home/anon/Projects/Data/SPyM1/Development";
OutputFile <- "G8bA8P2tStat";

FilterHeader4 <- paste("Lambda = ", as.character(Lambda), sep="");
FilterHeader5 <- paste("MinArrays = ", as.character(MinArrays), sep="");
FilterHeader6 <- paste("MinProbes = ", as.character(MinProbes), sep="");   
