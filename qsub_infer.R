args <- commandArgs(T)
x <- as.numeric(args[1])

source("batch_inference.R")
##Network_inference(x)
GMRF_inference(x)

