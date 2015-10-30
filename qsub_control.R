args <- commandArgs(T)
x <- as.numeric(args[1])

source("batches.R")
#control_case_meta(x)
batch_control()
