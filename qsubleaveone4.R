args <- commandArgs(T)
x <- as.numeric(args[1])

source("batches.R")
#batch_ASD_leaveone4(x)
batch_leaveone4(x)

