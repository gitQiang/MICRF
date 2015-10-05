#####################################################################################################
##DAWN package example code
#####################################################################################################

##load module assignment, the module is assigned by WGCNA. 
##first column is the gene name, second column is the module membership.
module=read.csv("module.csv")
 
##load gene expression data, p*n matrix, p is number of genes and n is number of samples. the row names 
##are the gene names.
load("geneexp.RData")
 

##load p-value
##first column is the gene name, second column is the gene p-value
pvalue_data=read.csv("pvalue.csv") 
  

source("DAWN_sourcefile.R")

################################################################################
##main function of DAWN, output a csv report file in the working directory
##
#function parameters listed below
#expression_data: gene expression data, p*n matrix, p is number of genes and n 
#is number of samples. the row names are the gene names.
#module_data: module assignment, the module is assigned by WGCNA. first column 
#is the gene name, second column is the module membership.
#pvalue_data: p-value file,first column is the gene name, second column is the 
#gene p-value
#cor_thres1: threshold of supernode, default 0.75
#cor_thres2: threshold of edge, default 0.7
#nASDthres: threshold of posterior probability for network ASD genes, default 0.5
#fdrthres: threshold of FDR for risk ASD genes, default 0.1
################################################################################
DAWN_package(expression_data=data,module_data=module, pvalue_data=pvalue_data)



###############################################################################
##Output file:
#Gene: Gene name 
#nASD: indicator of network ASD genes
#rASD: indicator of risk ASD genes
#original_pvalue: input p-value
#Module: module membership
#Class: supernode membership
#FDR: FDR value of the second stage
#Stage1_posterior: posterior probability of the stage I
#Stage2_posterior: posterior probability of the stage II
#risk_neighbor: number of risk neighbors in the network
#total_neighbor: number of total neighbors in the network
###############################################################################