args <- commandArgs(T)
x <- as.numeric(args[1])

source("batches.R")
#batch_ASD_randset_1_faster(x)
batch_randset_1_faster(x)

