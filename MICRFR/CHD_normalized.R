CHD_normalized <- function(){
    options(stringsAsFactors=FALSE)
    ## GTEx herat expression data: count table; two batches: Heart - Left Ventricle and Heart - Atrial Appendage;
    ## DESeq normalization
    ## edgeR remove batch effect
    dataFile <- c("data/expressiondata/GTEx.heart.countstable.txt")
    batchFile <- "data/expressiondata/renameFileWithHeartTissue1125.list"
    datExpr <- read.delim(dataFile,sep="\t",row.names=1)
    datExpr <- datExpr[rowSums(abs(datExpr))>0,]
    
    # get the batch information
    samplenames <- as.matrix(read.delim(batchFile,header=FALSE,sep="\t"))
    samplenames <- cbind(samplenames,paste(samplenames[,1],1,sep=""))
    samplenames <- cbind(samplenames,paste(samplenames[,1],1,sep="_"))
    # gene name mapping
    genes <- rownames(datExpr)
    datExpr1 <- datExpr[,samplenames[samplenames[,4]=="Heart - Left Ventricle",5]]
    datExpr2 <- datExpr[,samplenames[samplenames[,4]=="Heart - Atrial Appendage",5]]
  
    library(DESeq)
    condition <- rep(1,157)
    condition[115:157] <- 2
    cds = newCountDataSet( cbind(datExpr1,datExpr2), condition )
    cds = estimateSizeFactors( cds )
    exprM = counts( cds, normalized=TRUE )
    write.table(exprM[,1:114],file="data/expressiondata/GTExbatch1.txt",quote=FALSE,sep="\t")
    write.table(exprM[,115:157],file="data/expressiondata/GTExbatch2.txt",quote=FALSE,sep="\t")
    
    library(edgeR)
    datExpr <- removeBatchEffect(exprM,batch=condition)
    loc <- prcomp(t(datExpr)) # PCA plot
    plot(loc$x[,1],loc$x[,2],main="PCA plot",xlab="PC1",ylab="PC2")
    subs <- loc$x[,1] <1e7 & loc$x[,2]> -1e6
    datExpr <- datExpr[,subs]
    loc <- prcomp(t(datExpr)) # PCA plot
    plot(loc$x[,1],loc$x[,2],main="PCA plot",xlab="PC1",ylab="PC2")
    write.table(datExpr,file="data/expressiondata/GTExnormal.txt",quote=FALSE,sep="\t")
    
    # =========================================================================
    # PCGC fetal expression data: FPKM values; tow batches: 5pr and Nextera
    # log(FPKM+1) transform
    filename <- "data/expressiondata/141203.nextera_5prime_libraries.csv"
    exprname <- "data/expressiondata/fetal_JUNE_AUG_expression.txt"
    #groupN <- read.csv(filename)
    sam1 <- unlist(read.table("data/expressiondata/5pr.list")) #groupN[groupN[,4]=="5pr",1]
    sam2 <- unlist(read.table("data/expressiondata/nextera.list")) #groupN[groupN[,4]=="Nextera",1]
    
    exprM0 <- read.delim(exprname,check.names=F,row.names=1)
    exprM <- exprM0
    samplesN <- colnames(exprM)
    # delete duplicated samples
    subs <- grepl("CHD",samplesN)
    exprM <- exprM[,!subs]
    samplesN <- samplesN[!subs]
    
    a <- 1:length(samplesN)
    for(i in 1:length(samplesN)){
        a[i] <- substr(samplesN[i],1,7)
    }
    a[(197:216)-16] <- c("NFH12","NFH13","NHF82","NFH20","NFH9","NFH12","NFH13","NFH8","NFH20","NFH9","NFH12","NFH13","NFH8","NFH20","NFH9","NFH12","NFH13","NFH8","NFH20","NFH9")
    subs <- which(a%in% union(sam1,sam2))
    a <- a[subs]
    samplesN <- samplesN[subs]
    
    exprM <- exprM[,samplesN]
    exprM <- exprM[rowSums(abs(exprM))>0,]
    exprM <- log(exprM+1)
    datExpr1 <- exprM[,samplesN[a %in% sam1]]
    datExpr2 <- exprM[,samplesN[a %in% sam2]]
  
    write.table(datExpr1,file="data/expressiondata/PCGCbatch1.txt",quote=FALSE,sep="\t")
    write.table(datExpr2,file="data/expressiondata/PCGCbatch2.txt",quote=FALSE,sep="\t")   
   
    # =======================================================
    # PCGC new born expression data: FPKM values; one batches
    # log(FPKM+1) transform
    filename <- "data/expressiondata/PCGC_newBorn_FPKM.txt"
    exprM0 <- read.delim(filename,check.names=F,row.names=1)
    exprM <- exprM0
    exprM <- exprM[rowSums(abs(exprM))>0,]
    exprM <- log(exprM+1)
    write.table(exprM,file="data/expressiondata/PCGCnewborn.txt",quote=FALSE,sep="\t")
    
}

mouse_prune_exp <- function(){
    options(stringsAsFactors=FALSE)
    idmap <- read.delim("data/expressiondata/HMD_HumanPhenotype.rpt",header=FALSE)[,c(1,4)] ## got from MGI:http://www.informatics.jax.org/homology.shtml
    mouseexp <- read.csv("data/expressiondata/MOUSE_HEART_EXPRESSION_DATS.csv")
    mousegenes <- mouseexp[mouseexp[,"Mean"] > 0.6459,"NimbleGen"] 
    mousegenes <- unique(mousegenes)
    #humangenes <- idmap[match(toupper(mousegenes),toupper(idmap[,1])),1]
    humangenes <- mousegenes
    mapT <- as.matrix(read.delim("Fmap0121.txt",header=FALSE,sep="\t"))
    
    source("Network_analysis.R")
    humangenes <- mapping_to(humangenes)
    humangenes <- unique(humangenes)
    
    humangenes1 <- mapT[match(humangenes,mapT[,2]),1]
    humangenes1 <- humangenes1[!is.na(humangenes1)]
    
    dirstr <- "data/expressiondata/"
    netfiles <- c("PCGCbatch1.txt","PCGCbatch2.txt")
    netfiles <- paste(dirstr,netfiles,sep="")
    
    netfiles1 <- c("PCGCbatch1_truned.txt","PCGCbatch2_truned.txt")
    netfiles1 <- paste(dirstr,netfiles1,sep="")
    
    for(i in 1:length(netfiles)){
        datExp <- read.delim(netfiles[i],sep="\t",row.names=1)
        if(grepl("GTEx",netfiles[i])){
            subs <- rownames(datExp) %in% humangenes
        }else{
            subs <- rownames(datExp) %in% humangenes1
        }
        datExp <- datExp[subs,]
        print(dim(datExp))
        write.table(datExp,file=netfiles1[i],sep="\t",quote=FALSE)
    }
    
}

build_module <- function(){
    library(WGCNA)
    dirstr <- "data/expressiondata/"
    filenames <- c("data/expressiondata/GTExbatch1.txt","data/expressiondata/GTExbatch2.txt","data/expressiondata/PCGCbatch1.txt","data/expressiondata/PCGCbatch2.txt","data/expressiondata/PCGCnewborn.txt","data/expressiondata/PCGCbatch1_truned.txt","data/expressiondata/PCGCbatch2_truned.txt")
    strname <- c("data/expressiondata/GTExbatch1module","data/expressiondata/GTExbatch2module","data/expressiondata/PCGCbatch1module","data/expressiondata/PCGCbatch2module","data/expressiondata/PCGCnewbornmodule","data/expressiondata/PCGCbatch1module_t","data/expressiondata/PCGCbatch2module_t")
    pows <- c(6,14,20) # 6:1; 7: 7
    for(i in 1:3){
        filename <- filenames[i]
        datExpr <- read.delim(filename,check.names=F,row.names=1)
        distgene = adjacency(t(datExpr), type = "unsigned", power=1);
    
        powers = c(c(1:10), seq(from = 12, to=20, by=2))
        sft = pickSoftThreshold.fromSimilarity(distgene, powerVector = powers, verbose = 5)
        softPower <- sft$fitIndices[sft$fitIndices[,2]==max(sft$fitIndices[,2]),1] #sft$fitIndices[,1],-sign(sft$fitIndices[,3]*sft$fitIndices[,2]), #sft$fitIndices[,5]

        minModuleSize = 100; # set the minimum module size relatively high:
        #diag(distgene) <- 1
        adjacency = adjacency.fromSimilarity(distgene, type = "unsigned", power=softPower)
        TOM = TOMsimilarity(adjacency);
        dissTOM = 1-TOM
        geneTree = flashClust(as.dist(dissTOM), method = "average"); # Call the hierarchical clustering function
        dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize); # Module identification using dynamic tree cut
        dynamicColors = labels2colors(dynamicMods) # Convert numeric lables into colors
        moduleColors = dynamicColors
        colorOrder = c("grey", standardColors());
        Labels = match(moduleColors, colorOrder)-1;
        module <- cbind(rownames(distgene),Labels)
    
        module[module[,2]==0,2] <- max(as.numeric(module[,2])) + 1
        save(module,file=strname[i])
    }
    
    
    source("Network_analysis.R")
    filenames <- c("data/expressiondata/GTExbatch1module","data/expressiondata/GTExbatch2module")
    load(filenames[1])
    module[,1] <- mapping_to(module[,1])
    save(module,file="data/expressiondata/GTExbatch1modulem")
    rm(module)
    load(filenames[2])
    module[,1] <- mapping_to(module[,1])
    save(module,file="data/expressiondata/GTExbatch2modulem")
    
}

build_net <- function(){
    options(stringsAsFactors=FALSE)
    source("ASD_data_set.R")
    source("Network_analysis.R")
    
    dirstr <- "data/expressiondata/"
    filenames <- c("GTExbatch1.txt","GTExbatch2.txt","PCGCbatch1.txt","PCGCbatch2.txt","PCGCnewborn.txt","PCGCbatch1_truned.txt","PCGCbatch2_truned.txt")
    strname <- c("GTExbatch1module","GTExbatch2module","PCGCbatch1module","PCGCbatch2module","PCGCnewbornmodule","PCGCbatch1module_t","PCGCbatch2module_t")
    netfiles <- c("GTExbatch1net.txt","GTExbatch2net.txt","PCGCbatch1net.txt","PCGCbatch2net.txt","PCGCnewbornnet.txt","PCGCbatch1net_t.txt","PCGCbatch2net_t.txt")
    filenames <- paste(dirstr,filenames,sep="")
    strname <- paste(dirstr,strname,sep="")
    netfiles <- paste(dirstr,netfiles,sep="")
    
    for(k in 1:length(filenames)){
        filename <- filenames[k]
        
        if(k >1) rm(module)
        
        data <- read.delim(filename,check.names=F,row.names=1)
        load(strname[k])
       
        data <- data[rownames(data) %in% module[,1],]
        Labels <- unique(module[,2])
        allnet<- c()    
        for(i in Labels){
            modgenes <- module[module[,2]==i,1]
            modexp <- data[match(modgenes,rownames(data)),]
            corm=abs(cor(t(modexp),use='pair',nThreads=3))
            diag(corm) <- 0
        
            corm1 <- select_para(corm,d=5)$corm1
            edges <- which(corm1==1,arr.ind=TRUE)
            net <- cbind(modgenes[edges[,1]],modgenes[edges[,2]])
            net <- cbind(net,corm[edges])
            allnet <- rbind(allnet,net)
            print(i)
        }
    
        if(grepl("GTEx",filename)){
            genes <- union(allnet[,1],allnet[,2])
            genes1 <- mapping_to(genes)
            allnet[,1] <- genes1[match(allnet[,1],genes)]
            allnet[,2] <- genes1[match(allnet[,2],genes)]
        }
       
        write.table(allnet,file=netfiles[k],quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
    }
    
}

mouse_prune_net <- function(){
    options(stringsAsFactors=FALSE)
    idmap <- read.delim("data/expressiondata/HMD_HumanPhenotype.rpt",header=FALSE)[,c(1,4)] ## got from MGI:http://www.informatics.jax.org/homology.shtml
    mouseexp <- read.csv("data/expressiondata/MOUSE_HEART_EXPRESSION_DATS.csv")
    mousegenes <- mouseexp[mouseexp[,"Mean"] > 0.6459,"NimbleGen"] 
    mousegenes <- unique(mousegenes)
    #humangenes <- idmap[match(toupper(mousegenes),toupper(idmap[,1])),1]
    humangenes <- mousegenes
    mapT <- as.matrix(read.delim("Fmap0121.txt",header=FALSE,sep="\t"))
    
    source("Network_analysis.R")
    humangenes <- mapping_to(humangenes)
    humangenes <- unique(humangenes)
    
    humangenes1 <- mapT[match(humangenes,mapT[,2]),1]
    humangenes1 <- humangenes1[!is.na(humangenes1)]
    
    dirstr <- "data/expressiondata/"
    netfiles <- c("GTExbatch1net.txt","GTExbatch2net.txt","PCGCbatch1net.txt","PCGCbatch2net.txt","PCGCnewbornnet.txt")
    netfiles <- paste(dirstr,netfiles,sep="")
    
    netfiles1 <- c("GTExbatch1net_truned.txt","GTExbatch2net_truned.txt","PCGCbatch1net_truned.txt","PCGCbatch2net_truned.txt","PCGCnewbornnet_truned.txt")
    netfiles1 <- paste(dirstr,netfiles1,sep="")
    
    for(i in 1:length(netfiles)){
         allnet <- read.table(netfiles[i])
         print(length(union(allnet[,1],allnet[,2])))
         if(grepl("GTEx",netfiles[i])){
            subs <- allnet[,1] %in% humangenes & allnet[,2] %in% humangenes
         }else{
            subs <- allnet[,1] %in% humangenes1 & allnet[,2] %in% humangenes1
         }
         allnet <- allnet[subs,]
         print(length(union(allnet[,1],allnet[,2])))
         write.table(allnet,file=netfiles1[i],row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)
    }
    
}

dat_normalized1 <- function(){
    # log transform
    filename <- "141203.nextera_5prime_libraries.csv"
    exprname <- "fetal_JUNE_AUG_expression.txt"
    groupN <- read.csv(filename)
    tmp <- paste(groupN[,1],groupN[,3],sep="")
    tmp <- gsub("-","",tmp)
    n <- dim(groupN)[1]
    for(i in 1:n){
        if(substr(tmp[i],1,3)=="CHD"){
            tmp[i] <- substr(tmp[i],1,7)
        }else{
            tmp[i] <- substr(tmp[i],1,6)
        }
    }
    groupN  <- cbind(groupN,tmp)
    set1 <- groupN[groupN[,4]=="5pr",5]
    
    #========
    exprM0 <- read.delim(exprname,check.names=F,row.names=1)
    exprM <- exprM0
    samplesN <- colnames(exprM)
    # delete duplicated samples
    subs <- grepl("CHD",samplesN)
    exprM <- exprM[,!subs]
    # delete no expression genes
    exprM <- exprM[rowSums(abs(exprM))>0,]
    
    #exprM <- exprM + 1 # for log transform
    #exprM <- log(exprM)
    
    #========
    samplesN <- colnames(exprM)
    samplesN <- gsub("-","",samplesN)
    n <- length(samplesN)
    for(i in 1:n){
        if(substr(samplesN[i],1,3)=="CHD"){
            samplesN[i] <- substr(samplesN[i],1,7)
        }else{
            samplesN[i] <- substr(samplesN[i],1,6)
        }
    }
    
    Label <- 1:n
    for(i in 1:n){
        if(samplesN[i] %in% set1){
            Label[i] <- "A"
        }else if(substr(samplesN[i],1,5)=="fetal"){
            Label[i] <- "B"
        }else{
            Label[i] <- "B"
        }    
    }
    
    #======= MA plot for two batches
    R <- rowMeans(exprM[,Label=="A"])
    G <- rowMeans(exprM[,Label=="B"])
    edgeR::maPlot(R,G)
    exprMn <- log(exprM+1)
    R <- rowMeans(exprMn[,Label=="A"])
    G <- rowMeans(exprMn[,Label=="B"])
    edgeR::maPlot(R,G)
    
#     M <- log2(R/G)
#     A <- 0.5*(log2(R)+log2(G))
#     plot(A,M)
    
    batch1 <- exprMn[,Label=="A"]
    batch1m <- matrix(R,ncol=dim(batch1)[2],nrow=dim(batch1)[1])
    yR <- rowSums((batch1-0.6*batch1m)<0)/dim(batch1)[2]
    plot(R,yR)

    batch2 <- exprMn[,Label=="B"]
    batch2m <- matrix(G,ncol=dim(batch2)[2],nrow=dim(batch2)[1])
    yG <- rowSums((batch2-0.6*batch2m)<0)/dim(batch2)[2]
    plot(G,yG)

    #=== DESeq normalization and MA plot agains
    brainexp <- read.csv("genes_matrix_csv/expression_matrix.csv",header=FALSE,row.names=1)
    # delete no expression genes
    brainexp <- brainexp[rowSums(abs(brainexp))>0,]

    R <- rowMeans(brainexp)
    expm <- matrix(R,ncol=dim(brainexp)[2],nrow=dim(brainexp)[1])
    yR <- rowSums((brainexp-0.6*expm)<0)/dim(brainexp)[2]
    plot(R,yR)
    plot(R[R<10],yR[R<10])
    plot(R[R<0.01],yR[R<0.01])
    #Labels <- read.csv("genes_matrix_csv/columns_metadata.csv")[,"structure_acronym"]

    ## count table
    dataFile <- c("GTEx.heart.countstable.txt")
    batchFile <- "renameFileWithHeartTissue1125.list"
    datExpr <- read.delim(dataFile,sep="\t",row.names=1)
    datExpr <- datExpr[rowSums(abs(datExpr))>0,]
    
    # get the batch information
    samplenames <- as.matrix(read.delim(batchFile,header=FALSE,sep="\t"))
    samplenames <- cbind(samplenames,paste(samplenames[,1],1,sep=""))
    samplenames <- cbind(samplenames,paste(samplenames[,1],1,sep="_"))
    # gene name mapping
    genes <- rownames(datExpr)
    datExpr1 <- datExpr[,samplenames[samplenames[,4]=="Heart - Left Ventricle",5]]
    datExpr2 <- datExpr[,samplenames[samplenames[,4]=="Heart - Atrial Appendage",5]]
    R <- rowMeans(datExpr1)
    G <- rowMeans(datExpr2)
    edgeR::maPlot(R,G)
    library(DESeq)
    condition <- rep(1,157)
    condition[115:157] <- 2
    cds = newCountDataSet( cbind(datExpr1,datExpr2), condition )
    cds = estimateSizeFactors( cds )
    exprM = counts( cds, normalized=TRUE )
    R <- rowMeans(exprM[,1:114])
    G <- rowMeans(exprM[,115:157])
    edgeR::maPlot(R,G)

    R <- rowMeans(exprM)
    expm <- matrix(R,ncol=dim(exprM)[2],nrow=dim(exprM)[1])
    yR <- rowSums((exprM-0.6*expm)<0)/dim(exprM)[2]
    plot(R[R<10],yR[R<10])
    
    cds = estimateDispersions(cds, method="blind")
    exprM = getVarianceStabilizedData(cds)
    #exprM = varianceStabilizingTransformation( cds )
    R <- rowMeans(exprM)
    expm <- matrix(R,ncol=dim(exprM)[2],nrow=dim(exprM)[1])
    yR <- rowSums((exprM-0.6*expm)<0)/dim(exprM)[2]
    plot(R,yR)
   

}

dat_normalized <- function(){
# log transform
filename <- "141203.nextera_5prime_libraries.csv"
exprname <- "fetal_JUNE_AUG_expression.txt"
groupN <- read.csv(filename)
tmp <- paste(groupN[,1],groupN[,3],sep="")
tmp <- gsub("-","",tmp)
n <- dim(groupN)[1]
for(i in 1:n){
    if(substr(tmp[i],1,3)=="CHD"){
        tmp[i] <- substr(tmp[i],1,7)
    }else{
        tmp[i] <- substr(tmp[i],1,6)
    }
}
groupN  <- cbind(groupN,tmp)
set1 <- groupN[groupN[,4]=="5pr",5]

#========
exprM0 <- read.delim(exprname,check.names=F,row.names=1)
exprM <- exprM0
samplesN <- colnames(exprM)
# delete duplicated samples
subs <- grepl("CHD",samplesN)
exprM <- exprM[,!subs]
# delete outlier genes
n.genes <- dim(exprM)[1]
n.sample <- dim(exprM)[2]
for(i in 1:n.sample){
    boxout <- boxplot(exprM[exprM[,i] > 1,i],plot=FALSE)
    tmp <- exprM[,i] > boxout$stats[5]
    exprM[tmp,i] <- 0
    #exprM[tmp,] <- 0
}
# delete no expression genes
exprM <- exprM[rowSums(abs(exprM))>0,]
exprM <- exprM + 1 # for log transform
exprM <- log(exprM)

#========
samplesN <- colnames(exprM)
samplesN <- gsub("-","",samplesN)
n <- length(samplesN)
for(i in 1:n){
    if(substr(samplesN[i],1,3)=="CHD"){
        samplesN[i] <- substr(samplesN[i],1,7)
    }else{
        samplesN[i] <- substr(samplesN[i],1,6)
    }
}

Label <- 1:n
for(i in 1:n){
    if(samplesN[i] %in% set1){
        Label[i] <- "A"
    }else if(substr(samplesN[i],1,5)=="fetal"){
        Label[i] <- "B"
    }else{
        Label[i] <- "B"
    }	
}

#=======
library(edgeR)
datExpr <- removeBatchEffect(exprM,batch=Label)
loc <- prcomp(t(datExpr)) # PCA plot
plot(loc$x[,1],loc$x[,2],main="PCA plot",xlab="PC1",ylab="PC2")

# remove the batch effect again
Label <- rep("C",dim(datExpr)[2])
Label[loc$x[,1] < -60] <- "A"
Label[loc$x[,2] < -40] <- "B"
datExpr <- removeBatchEffect(datExpr,batch=Label)
loc <- prcomp(t(datExpr))
plot(loc$x[,1],loc$x[,2],main="PCA plot",xlab="PC1",ylab="PC2")

datExpr <- t(datExpr)
# clustering and outlier delete; WCGNA delete the outlier data
library(WGCNA)
gsg <- goodSamplesGenes(datExpr,verbose=3)    # print(gsg$allOK)
if (!gsg$allOK){ datExpr = datExpr[gsg$goodSamples, gsg$goodGenes];}
sampled <- dist(datExpr)
sampleTree = flashClust(sampled,method="average")
plot(sampleTree,main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)

write.table(datExpr,file="PCGCall.txt",quote=FALSE,sep="\t",row.names=F)

}

module_label <- function(){
    #==============
    distflag <- 2
    source("Network_analysis.R")
    #check_power(distgene,datExpr,2)
    softPower <- 5
    flag <- 1
    Labels <- WGCNA_module(distgene,datExpr,distflag,flag,para1=softPower,para3=0.01,minModuleSize = 10)
    genes <- colnames(datExpr)
    module <- cbind(genes,Labels)
    write.table(module,file="co-express_modules_0.01_10.txt",quote=F,row.names=F,col.names=F,sep="\t")
    

}

inference_net <- function(){
    
source("Network_analysis.R")
datExpr <- as.matrix(read.delim("PCGCall.txt",sep="\t",row.names=1))
datExpr <- t(datExpr)
#distgene = adjacency(datExpr, type = "unsigned", power=1) # co-expression network
distgene <- 0
module <- as.matrix(read.delim("co-express_modules_0.01_10.txt",sep="\t",header=F))
nodeflag <- 1
nodesim <- integrate_node(nodeflag)
#==========
#source("Bayes_build.R")
#Bayes_build(nodesim,module,datExpr,distgene)

# GRNFile <- "GRN_net_map.txt"
# GRNTable <- as.matrix(read.delim(GRNFile,sep="\t",header=FALSE))
# GRNTable <- GRNTable[GRNTable[,5]=="all",]
source("CRF_build.R")
strn="PCGC"
CRF_build(nodesim,module,datExpr,strn,nodeflag)

}

integrate_node <- function(flag = 1){

    #===============
    # TADA information
    mapT <- as.matrix(read.delim("Fmap.txt",header=FALSE,sep="\t"))
    genes <- mapT[,1]
    
    TADAFile <- "TADA_hardfilterednew0113.csv"
    TADAFile <- "TADA_lofmis1202.csv"
    geneinfo <- as.matrix(read.csv(TADAFile))
    genelist <- geneinfo[,1]
    genelist <- mapping_to(genelist)
    genelist <- mapT[match(genelist,mapT[,2]),1]
    nodesim <- as.matrix((as.numeric(geneinfo[,"BF"])*modratio)/(as.numeric(geneinfo[,"BF"]*modratio) + 1))
    names(nodesim) <- genelist
    subs <- !is.na(genelist)
    nodesim <- nodesim[subs]
    #write.table(cbind(names(nodesim),nodesim),file="TADAinfo.txt",quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")
    if(flag==1){
    nodeinfo <- matrix(,length(genes),3,dimnames=list(genes,1:3))
    nodeinfo[,1] <- 0.1 ## related genes percent estimate
    nodeinfo[,2:3] <- 0
    nodeinfo[names(nodesim),1] <- nodesim
    
    # CNV information; filtered: some CNVs affect a very large number of genes. double counting closely related genes like "UBXN10" and "UBXN10-AS1" 
    CNVfile <- "case.de.novo.cnv.genes.bed"
    cnvinfo <- read.delim(CNVfile,header=FALSE,sep="\t")
    ts <- table(as.vector(cnvinfo[,5]))
    boxout <- boxplot(ts,plot=FALSE)
    tmp <- names(ts)[ts <= boxout$stats[5]]
    cnvinfo <- cnvinfo[cnvinfo[,5] %in% tmp,]
    cnvinfo[,4] <- mapping_to(cnvinfo[,4],"")
    cnvgenes <- unique(cnvinfo[,4])
    cnvgenes <- mapT[match(cnvgenes,mapT[,2]),1]
    cnvgenes <- cnvgenes[!is.na(cnvgenes)]
    
    nodeinfo[cnvgenes,2] <-1
    
    # Mendelian CHD genes # SNP-associated genes # Environmental target genes
    Otherfile <- "Othergenes.txt"
    otherinfo <- read.delim(Otherfile,header=FALSE,sep="\t")
    othergenes <- otherinfo[,2]
    othergenes <- mapping_to(othergenes)
    othergenes <- mapT[match(othergenes,mapT[,2]),1]
    othergenes <- othergenes[!is.na(othergenes)]
    othergenes <- unique(othergenes)
    
    nodeinfo[othergenes,3] <- 1 
    }else{
        nodeinfo <- matrix(,length(genes),1,dimnames=list(genes,1))
        nodeinfo[,1] <- 0.5 ## related genes percent estimate
        nodeinfo[names(nodesim),1] <- nodesim
    }
    nodeinfo
}

mapping_PPI <- function(){
    options(stringsAsFactors=FALSE)
    genes <- as.matrix(read.table("GENES"))
    
    strPPI <- as.matrix(read.table("STRINGnet",header=FALSE))
    mapPPI <- as.matrix(read.delim("ENSG_SP.txt",sep="\t"))
    #mapPPI <- mapPPI[mapPPI[,1] %in% genes,]
    nodes <- union(strPPI[,1],strPPI[,2])
    nodes1 <- mapPPI[match(nodes,mapPPI[,2]),1]
    mapnod <- cbind(nodes,nodes1)
    strPPInew <- strPPI
    strPPInew[,1] <- mapnod[match(strPPInew[,1],mapnod[,1]),2]
    strPPInew[,2] <- mapnod[match(strPPInew[,2],mapnod[,1]),2]
    strPPInew <- strPPInew[!is.na(strPPInew[,1]) & !is.na(strPPInew[,2]),]
    write.table(strPPInew,file="STRINGnetmap.txt",quote=F,col.names=F,row.names=F,sep="\t")
    
    funPPI <- as.matrix(read.table("FIs_043009.txt",header=FALSE))
    mapPPI <- as.matrix(read.delim("ENSG_Uniprot.txt",sep="\t"))
    mapPPI <- mapPPI[mapPPI[,1] %in% genes,]
    nodes <- union(funPPI[,1],funPPI[,2])
    nodes1 <- mapPPI[match(nodes,mapPPI[,2]),1]
    mapnod <- cbind(nodes,nodes1)
    funPPInew <- funPPI[,1:2]
    funPPInew[,1] <- mapnod[match(funPPInew[,1],mapnod[,1]),2]
    funPPInew[,2] <- mapnod[match(funPPInew[,2],mapnod[,1]),2]
    funPPInew <- funPPInew[!is.na(funPPInew[,1]) & !is.na(funPPInew[,2]),]
    write.table(funPPInew,file="FUNCTIONnetmap.txt",quote=F,col.names=F,row.names=F,sep="\t")

}

inference_net_brain <- function(){
    
    source("Network_analysis.R")
    datExpr <- as.matrix(read.delim("genes_matrix_csv/Brain_expression_Ref.txt",sep="\t",row.names=1))
    datExpr <- t(datExpr)
    #distgene = adjacency(datExpr, type = "unsigned", power=1) # co-expression network
    distgene <- 0
    
    module <- read.table("brain_modules_0.01_10.txt",sep="\t",header=F,strip.white=TRUE)
    #===============
    mapT <- as.matrix(read.delim("Fmap.txt",header=FALSE,sep="\t"))
    TADAFile <- "NDD_TADA1205.csv"
    geneinfo <- as.matrix(read.csv(TADAFile))
    genelist <- geneinfo[,1]
    genelist <- mapping_to(genelist)
    genelist <- mapT[match(genelist,mapT[,2]),1]
    #mutatedrate <- as.numeric(geneinfo[,2])
    #mutatednum <- as.numeric(geneinfo[,"dn.LoF"])+as.numeric(geneinfo[,"dn.mis3"])
    nodesim <- as.numeric(geneinfo[,"BF"])/(as.numeric(geneinfo[,"BF"]) + 1)
    nodesim <- cbind(nodesim,1-nodesim)
    rownames(nodesim) <- genelist
    subs <- !is.na(genelist)
    nodesim <- nodesim[subs,]
    #==========
    #source("Bayes_build.R")
    #Bayes_build(nodesim,module,datExpr,distgene)
    source("CRF_build.R")
    strn="brainNDD"
    CRF_build(nodesim,module,datExpr,strn)
    
    TADAFile <- "nonNDD_TADA1205.csv"
    geneinfo <- as.matrix(read.csv(TADAFile))
    genelist <- geneinfo[,1]
    genelist <- mapping_to(genelist)
    genelist <- mapT[match(genelist,mapT[,2]),1]
    nodesim <- as.numeric(geneinfo[,"BF"])/(as.numeric(geneinfo[,"BF"]) + 1)
    nodesim <- cbind(nodesim,1-nodesim)
    rownames(nodesim) <- genelist
    subs <- !is.na(genelist)
    nodesim <- nodesim[subs,]
    source("CRF_build.R")
    strn="brainnonNDD"
    CRF_build(nodesim,module,datExpr,strn)    
}

mcl_R <- function(net){
    adj <- net$matrix
    diag(adj) <- 1    
    adj <- sweep(adj, 2, colSums(adj), FUN = "/")
    while(TRUE){
        adj0 <- adj
        
        # expansion
        while(TRUE){
            adjtmp <- adj
            adj <- adj %*% adj
            tmp <- adj - adjtmp
            if(norm(tmp,"O") < 10^-3 ) break;
        }
         
        # inflation
        power <- 2
        for(i in 1:(power-1)){
            adj <- adj * adj
        }
        adj <- sweep(adj, 2, colSums(adj), FUN="/")
        
        
        tmp <- adj - adj0
        if(norm(tmp,"O") < 10^-3 ) break;
    }
    
    # merge clusters    

}