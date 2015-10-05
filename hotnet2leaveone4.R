args <- commandArgs(T)
x <- as.numeric(args[1])

source("batches.R")
hotnet2_leaveone4(x)
