args <- commandArgs(T)
x <- as.numeric(args[1])

source("BestP_8_17.R")
AUC_faster(x)