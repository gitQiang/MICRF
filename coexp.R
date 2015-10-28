coexp <- function(){
    load("DAWNgeneexp.Rdata")
    


}


# check not expressed genes -----------------------------------------------
checkg <- function(){
    load("DAWNgeneexp.Rdata")
    genes <- rownames(data)
    source("enrichana_6_12.R")
    filenames <- c("result/Nat1377_meta/TADAdenovo_Meta_dmis.csv","result/Nat1377_meta/DAWN/DAWN_Meta_dmis.csv","result/Nat1377_meta/MAGI/RandomGeneList.3","result/Nat1377_meta/RWR/coexp1/hotnetresult1meta_dmis21.txt","result/Nat1377_meta/RWR/HPRD/hotnetresult1meta_dmis20.txt","result/Nat1377_meta/RWR/STRING/hotnetresult1meta_dmis31.txt","result/Nat1377_meta/RWR/iRef/hotnetresult1meta_dmis6.txt","result/Nat1377_meta/RWR/coexp/hotnetresult1meta_dmis7.txt","result/Nat1377_meta/MICRF/v4/coexp/CRFresult_0.2meta_dmisLBP_7.txt") 
    tmp <- read.csv(filenames[1])
    
    Tadag <- tmp[tmp[,"qvalue.dn"]<0.3,1] ## 400 genes
    
    length(setdiff(Tadag,genes)) # 74

}

# brain 3-6 ---------------------------------------------------------------
expressed <- function(){
    ## RNA-seq data
    # to data matrix 
    filename <- "data/brain_expression/expression_matrix.csv"
    genefile <- "data/brain_expression/rows_metadata.csv"
    samplefile <- "data/brain_expression/columns_metadata.csv"
    regions <- c("OFC","DFC","VFC","MFC","M1C","S1C","IPC","A1C","STC","ITC","V1C")
    periods <- paste(10:24," pcw",sep="")
    genesM <- as.matrix(read.csv(genefile))
    genenames <- genesM[,"gene_symbol"]
    sampleM <- as.matrix(read.csv(samplefile))
    subs1 <- sampleM[,"structure_acronym"] %in% regions
    subs2 <- sampleM[,"age"] %in% periods
    subs <- which(subs1 & subs2)
    samplenames <- gsub(" ","",sampleM[subs,2])
    
    datExpr <- read.csv(filename,header=FALSE,row.names=1)
    datExpr <- as.matrix(datExpr[,subs])
    rownames(datExpr) <- genenames
    colnames(datExpr) <- samplenames
    
    # delete no expression gene
    datExpr <- datExpr[rowSums(datExpr)>0,]
    
    # log2 transform
    datExpr1 <- log2(datExpr+1)
    #===============================================
    source("Network_analysis.R")
    genes <- rownames(datExpr1)
    genes <- mapping_to(genes)
    rownames(datExpr1) <- genes

    Tadag <- read.csv("result/Nat1377_meta/TADAdenovo_Meta_dmis.csv")[,1]
    datExpr2 <- datExpr1[genes %in% Tadag,]
    save(datExpr2,file="data/brain_expression/Brainexp")
    
    ##--------------------------------------------------------------
    ## microarray data
    # to data matrix 
    filename <- "data/brain_expression/gene_array_matrix_csv/expression_matrix.csv"
    genefile <- "data/brain_expression/gene_array_matrix_csv/rows_metadata.csv"
    samplefile <- "data/brain_expression/gene_array_matrix_csv/columns_metadata.csv"
    regions <- c("OFC","DFC","VFC","MFC","M1C","S1C","IPC","A1C","STC","ITC","V1C")
    periods <- paste(10:24," pcw",sep="")
    genesM <- as.matrix(read.csv(genefile))
    genenames <- genesM[,"gene_symbol"]
    sampleM <- as.matrix(read.csv(samplefile))
    subs1 <- sampleM[,"structure_acronym"] %in% regions
    subs2 <- sampleM[,"age"] %in% periods
    subs <- which(subs1 & subs2)
    samplenames <- gsub(" ","",sampleM[subs,2])
    
    datExpr <- read.csv(filename,header=FALSE,row.names=1)
    datExpr <- as.matrix(datExpr[,subs])
    rownames(datExpr) <- genenames
    colnames(datExpr) <- samplenames
    
    source("Network_analysis.R")
    genes <- rownames(datExpr)
    genes <- mapping_to(genes)
    rownames(datExpr) <- genes
    
    Tadag <- read.csv("result/Nat1377_meta/TADAdenovo_Meta_dmis.csv")[,1]
    sum(genes %in% Tadag)
    save(datExpr,file="data/brain_expression/Brainexp_micro")

}

# build coexp modules -----------------------------------------------------
modules <- function(){
    load("data/brain_expression/Brainexp_micro")
    module <- build_module(datExpr)
    save(module,file="data/module_coexp")
}

build_module <- function(datExpr){
    
    library(WGCNA)
    
#     distgene = adjacency(t(datExpr), type = "unsigned", power=1);
#     powers = c(c(1:10), seq(from = 12, to=20, by=2))
#     sft = pickSoftThreshold.fromSimilarity(distgene, powerVector = powers, verbose = 5)
#     softPower <- sft$fitIndices[sft$fitIndices[,2]==max(sft$fitIndices[,2]),1] #sft$fitIndices[,1],-sign(sft$fitIndices[,3]*sft$fitIndices[,2]), #sft$fitIndices[,5]
    
    # WGCNA to module relationships
    MEDissThres = 0.15
    minModuleSize = 200; # set the minimum module size relatively high:
    softPower = 6;
    adjacency = adjacency(t(datExpr), power = softPower);
    TOM = TOMsimilarity(adjacency);
    dissTOM = 1-TOM
    # Call the hierarchical clustering function
    geneTree = hclust(as.dist(dissTOM), method = "average");
    #geneTree = flashClust(as.dist(dissTOM), method = "average"); # Call the hierarchical clustering function
    dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize); # Module identification using dynamic tree cut
    dynamicColors = labels2colors(dynamicMods) # Convert numeric lables into colors
    merge = mergeCloseModules(t(datExpr), dynamicColors, cutHeight = MEDissThres, verbose = 3) # Call an automatic merging function
    moduleColors = merge$colors; # The merged module colors
    
    colorOrder = c("grey", standardColors());
    Labels = match(moduleColors, colorOrder)-1;
    module <- cbind(rownames(datExpr),Labels)
    module
}

build_net <- function(){
    load("data/brain_expression/Brainexp_micro")
    load("data/module_coexp")
    
    library(WGCNA)
    
    ## cor.test.pvalue ??
    Labels <- unique(module[,2])
    allnet<- c()
    for(i in Labels){
        modgenes <- module[module[,2]==i,1]
        modexp <- datExpr[match(modgenes,rownames(datExpr)),]
        corm=abs(cor(t(modexp),use='pair',nThreads=3))
        diag(corm) <- 0
        edges <- which(corm > 0.7,arr.ind=TRUE)
        net <- cbind(modgenes[edges[,1]],modgenes[edges[,2]])
        net <- cbind(net,corm[edges])
        allnet <- rbind(allnet,net)
        print(i)
    }
    write.table(allnet,file="data/network_inference/brainspan_net_cor.txt",quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
    
    ###===================
    source("ASD_data_set.R")
    Labels <- unique(module[,2])
    allnet<- c()    
    for(i in Labels){
        modgenes <- module[module[,2]==i,1]
        modexp <- datExpr[match(modgenes,rownames(datExpr)),]
        corm=abs(cor(t(modexp),use='pair',nThreads=3))
        diag(corm) <- 0
        
        corm1 <- select_para(corm,d=5)$corm1
        edges <- which(corm1==1,arr.ind=TRUE)
        net <- cbind(modgenes[edges[,1]],modgenes[edges[,2]])
        net <- cbind(net,corm[edges])
        allnet <- rbind(allnet,net)
        print(i)
    }
    write.table(allnet,file="data/network_inference/brainspan_net_top5.txt",quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")


}