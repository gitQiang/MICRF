allcheck <- function(){

expression_data=datExpr;module_data=module; pvalue_data=pvalue_data;cor_thres1=thres1;cor_thres2=thres2;nASDthres=0.5;fdrthres=0.1;flag_q=distflag;filename=filename


genestmp <- module[module[,2]==5,1]
datExpr <- datExpr[genestmp,]
module <- module[module[,2]==5,]
module[,2] <- 1
pvalue_data <- pvalue_data[pvalue_data[,1]%in% genestmp,]
distgene <- distgene[genestmp,genestmp]


#modulef8
genename=genenamecor
thres1 <- 1-cor_thres1
thres2 <- cor_thres2

#main_f9
#graph1,zscore,Istart,100,mustart,seedindex,sigma1^2,sigma2^2
graph <- graph1
iter <- 100
mui <- mustart
sigmas1 <- sigma1^2
sigmas2 <- sigma2^2

#emf7_node
#zscore3,newclass,theta,0,1.5,1,1
zscore <- zscore3
kclass <- newclass
mu1<- 0
mu2<- 1.5
sigma1<- 1
sigma2<- 1

# kpvalue=save1$keepgenep[keepi] #finalop
# 
# finalgene, save1$keepgene
# finalkclass save1$geneclass
# ,finalop,
# finalposter,
# finaltheta,
# flocalfdr

}

check_TADA <- function(){

TADAFile <- "TADA_lofmis1121.csv"
pvalue_data <- read.csv(TADAFile)[,c(1,10)]
pvalue_data[,1] <- mapping_to(pvalue_data[,1])
if(fileflag!=2)  pvalue_data[,1] <- mapT[match(pvalue_data[,1],mapT[,2]),1]; 
pvalue_data <- pvalue_data[!(is.na(pvalue_data[,1]) | is.na(pvalue_data[,2])),]
pvalue_data <- pvalue_data[pvalue_data[,2]< 1,]
pvalue_data <- pvalue_data[!duplicated(pvalue_data[,1]),]
pvalue_data0 <- pvalue_data
print(dim(pvalue_data0))

TADAFile <- "TADA_lofmis1202.csv"
pvalue_data <- read.csv(TADAFile)[,c(1,10)]
pvalue_data[,1] <- mapping_to(pvalue_data[,1])
if(fileflag!=2)  pvalue_data[,1] <- mapT[match(pvalue_data[,1],mapT[,2]),1]; 
pvalue_data <- pvalue_data[!(is.na(pvalue_data[,1]) | is.na(pvalue_data[,2])),]
pvalue_data <- pvalue_data[pvalue_data[,2]< 1,]
pvalue_data <- pvalue_data[!duplicated(pvalue_data[,1]),]
print(dim(pvalue_data))
for(i in 1:18){
gene1121 <- module[module[,2]==i,1]
p1121 <- pvalue_data0[match(gene1121,pvalue_data0[,1]),2]
p1202 <- pvalue_data[match(gene1121,pvalue_data[,1]),2]
subs <- !is.na(p1121)
subs1 <- !is.na(p1202)
print(all(subs==subs1))

p1121 <- p1121[subs]
p1202 <- p1202[subs]
plot(p1121,p1202)
}

}

check_RBFOX2 <- function(){
source("Network_analysis.R")
fileflag <- 4
distflag <- 1
cflag <- 0
flag <- 1
softpower <- 6
para3 <- 0.5
thres1 <- 0.65
thres2 <- 0.5
dataFile <- c("GTEx_heart_expression.txt","GTEx.heart.countstable.txt","PCGC_newBorn_FPKM.txt","fetal_JUNE_AUG_expression.txt")
batchFile <- "renameFileWithHeartTissue1125.list"
datExpr <- step_1_0_normalize(dataFile[fileflag],batchFile,fileflag)
load("dist_4_1")
distgene[is.na(distgene)] <- 0
distgene <- abs(distgene) # this for DAWN
genes <- colnames(datExpr)
dimnames(distgene) <- list(genes,genes)
distgene4_1 <- distgene
rm(distgene)

fileflag <- 1
datExpr <- step_1_0_normalize(dataFile[fileflag],batchFile,fileflag)
distfile <- paste("dist",fileflag,distflag,sep="_")
load(distfile)
distgene[is.na(distgene)] <- 0
distgene <- abs(distgene) # this for DAWN
genes <- colnames(datExpr)
dimnames(distgene) <- list(genes,genes)
distgene1_1 <- distgene
rm(distgene)

genes <- intersect(colnames(distgene4_1),colnames(distgene1_1))
n.genes <- length(genes)
b <- matrix(0,n.genes,4)
for(i in 1:n.genes){
    a1 <- distgene1_1[genes[i],genes]
    a2 <- distgene4_1[genes[i],genes]
    subs <- pmax(a1,a2) > 0.5
    b[i,1] <- sum(a1[subs] > 0.5)
    b[i,2] <- sum(a2[subs] > 0.5)
    b[i,3] <- sum(a2[subs] - 20 * a1[subs] >0 ) # appear
    b[i,4] <- sum(a1[subs] - 20 * a2[subs] >0 ) # disappear
}
plot(b[,1],b[,2])
plot(b[,3],b[,4])
plot(b[,3]/b[,4])
plot(b[,4]/b[,3])
result <- read.csv("enrichment/ComNet14_1_1_0.65_0.5_0.15.csv")
riskgenes <- result[result[,4]==1,2]

rsubs <- match(riskgenes,genes)


which(colnames(distgene4_1)=="ENSG00000100320")
which(colnames(distgene1_1)=="ENSG00000100320")
a1 <- distgene1_1[2190,genes]
a2 <- distgene4_1[2191,genes]
subs <- pmax(a1,a2) > 0.5
genemax <- genes[subs]

subnet1_1 <- distgene1_1[genemax,genemax]
subnet1_1[lower.tri(subnet1_1,diag=TRUE)] <- 0
subnet4_1 <- distgene4_1[genemax,genemax]
subnet4_1[lower.tri(subnet4_1,diag=TRUE)] <- 0

dis_app <- union(genemax[a2[subs]/a1[subs] > 20],genemax[a1[subs]/a2[subs] > 20])
dis_app <- union(dis_app,"ENSG00000100320")
subnet1_1 <- subnet1_1[dis_app,dis_app]
subnet4_1 <- subnet4_1[dis_app,dis_app]

subnet <- pmax(subnet1_1,subnet4_1)
edges <- which(subnet > 0.5,arr.ind=TRUE)
n <- dim(edges)[1]
net <- matrix(0,n,4)
weights <- subnet[edges]
net[,1] <- dis_app[edges[,1]]
net[,2] <- dis_app[edges[,2]]
net[,3] <- weights
net[subnet1_1[edges] > 0.5 & subnet4_1[edges] <= 0.5,4] <- 1
net[subnet1_1[edges] <= 0.5 & subnet4_1[edges] > 0.5,4] <- 2
net[subnet1_1[edges] > 0.5 & subnet4_1[edges] > 0.5,4] <- 3
net[,1] <- mapT[match(net[,1],mapT[,1]),2]
net[,2] <- mapT[match(net[,2],mapT[,1]),2]
net <- net[!is.na(net[,1]) & !is.na(net[,2]),]
colnames(net) <- c("Node1","Node2","weight","Condition")
write.table(net,file="RBFOX2Net.txt",quote=F,row.names=F,sep="\t")

genemax <- mapT[match(genemax,mapT[,1]),2]
topapp <- genemax[a2[subs]/a1[subs] > 20]
topdis <- genemax[a1[subs]/a2[subs] > 20]
write.table(topdis,file="disappeargenes",quote=F,row.names=F,col.names=F)
write.table(topapp,file="appeargenes",quote=F,row.names=F,col.names=F)


plot(a1[subs],a2[subs])
plot(a1[subs]/a2[subs]) # edge disappear
abline(h=200)
plot(a2[subs]/a1[subs]) # edge appeared
abline(h=20)

# edges1_1 <- which(subnet1_1>0,arr.ind=TRUE)
# weights <- subnet1_1[edges1_1]
# net1_1 <- cbind(genemax[edges1_1[,1]],genemax[edges1_1[,2]])
# net1_1 <- cbind(net1_1,weights)
# net1_1[,1] <- mapT[match(net1_1[,1],mapT[,1]),2]
# net1_1[,2] <- mapT[match(net1_1[,2],mapT[,1]),2]
# net1_1 <- net1_1[!is.na(net1_1[,1]) & !is.na(net1_1[,2]),]
# write.table(net1_1,file="RBFOX2normal.txt",quote=F,row.names=F,col.names=F,sep="\t")
# 
# edges4_1 <- which(subnet4_1 > 0,arr.ind=TRUE)
# weights <- subnet4_1[edges4_1]
# net4_1 <- cbind(genemax[edges4_1[,1]],genemax[edges4_1[,2]])
# net4_1 <- cbind(net4_1,weights)
# net4_1[,1] <- mapT[match(net4_1[,1],mapT[,1]),2]
# net4_1[,2] <- mapT[match(net4_1[,2],mapT[,1]),2]
# net4_1 <- net4_1[!is.na(net4_1[,1]) & !is.na(net4_1[,2]),]
# write.table(net4_1,file="RBFOX2PCGC.txt",quote=F,row.names=F,col.names=F,sep="\t")
}

check_DAWN_TADA <- function(){
    source("Network_analysis.R")
    library(WGCNA)
    mapT <- as.matrix(read.delim("Fmap.txt",header=FALSE,sep="\t"))
    # compare the DAWN rank and TADA rank
    DAWNFile <- "D1121result_4_1_0.65_0.45_0.5_1.csv"
    TADAFile <- "TADA_lofmis1202.csv"
    # only compare DAWN FDR<1 genes
    DAWNresult <- as.matrix(read.csv(DAWNFile))
    subs <- !is.na(DAWNresult[,"FDR"])
    DAWNfdr <- DAWNresult[subs,"FDR"]
    DAWNgenes <- DAWNresult[subs,"Gene"]
    TADAp <- DAWNresult[subs,"original_pvalue"]
    n <- length(DAWNfdr)
    DAWNsub <- sort(DAWNfdr,index.return=TRUE)$ix
    DAWNrank <- 1:n
    DAWNrank[DAWNsub] <- 1:n
#     TADAsub <- sort(TADAp,index.return=TRUE)$ix
#     TADArank <- 1:n
#     TADArank[TADAsub] <- 1:n
#     plot(DAWNrank,TADArank,main="DAWN rank VS TADA rank",xlab="DAWN rank",ylab="TADA rank")
#     subs <- which(DAWNresult[,3]=="EP300")
#     lines(DAWNrank[subs],TADArank[subs],type="p",col="red",lwd=10)

    pvalue_data <- read.csv(TADAFile)[,c(1,10)]
    pvalue_data[,1] <- mapping_to(pvalue_data[,1]) 
    pvalue_data <- pvalue_data[!(is.na(pvalue_data[,1]) | is.na(pvalue_data[,2])),]
    pvalue_data <- pvalue_data[pvalue_data[,2]< 1,]
    pvalue_data <- pvalue_data[!duplicated(pvalue_data[,1]),]
    TADAsub <- sort(as.numeric(pvalue_data[,2]),index.return=TRUE)$ix
    TADArank <- 1:length(TADAsub)
    TADArank[TADAsub] <- 1:length(TADAsub)
    TADAr <- TADArank[match(DAWNgenes,pvalue_data[,1])]
    plot(DAWNrank,TADAr,main="DAWN rank VS TADA rank",xlab="DAWN rank",ylab="TADA rank")

    load("dist_4_1")
    distgene[is.na(distgene)] <- 0
    distgene <- abs(distgene) 
    diag(distgene) <- 0 
    dataFile <- c("GTEx_heart_expression.txt","GTEx.heart.countstable.txt","PCGC_newBorn_FPKM.txt","fetal_JUNE_AUG_expression.txt","Brain_expression_Ref.txt")
    batchFile <- "renameFileWithHeartTissue1125.list"
    fileflag <- 4
    datExpr <- step_1_0_normalize(dataFile[fileflag],batchFile,fileflag)
    genes <- colnames(datExpr)
    dimnames(distgene) <- list(genes,genes)
    module <- read.delim("FPKM4_1_module.txt",sep="\t",header=F)

#     ONEset <- DAWNgenes[DAWNrank <= 30 & TADAr > 100]
#     TWOset <- DAWNgenes[DAWNrank > 80 & TADAr <= 80]
    subsone <- which(TADAr/DAWNrank > 1.5 & TADAr <= 300 & TADAr > 30)
    substwo <- which(DAWNrank/TADAr > 1.5 & DAWNrank <= 300 & DAWNrank > 30)
    ONEset <- DAWNgenes[subsone]
    TWOset <- DAWNgenes[substwo]

    ONEs1 <- mapT[match(ONEset,mapT[,2]),1]
    TWOs1 <- mapT[match(TWOset,mapT[,2]),1]

    plot(DAWNrank[1:300],TADAr[1:300],main="DAWN rank VS TADA rank",xlab="DAWN rank",ylab="TADA rank",xlim=c(1,300),ylim=c(1,300))
    abline(0,1,col=1,lwd=1)
    abline(0,1.5,col=2,lwd=1)
    abline(0,2/3,col=2,lwd=1)
    
    subs <- which(DAWNresult[,3] %in% as.character(ONEset))
    lines(DAWNrank[subs],TADAr[subs],type="p",col="red",lwd=5)
    subs <- which(DAWNresult[,3] %in% as.character(TWOset))
    lines(DAWNrank[subs],TADAr[subs],type="p",col="red",lwd=5)
    
    n.one <- length(ONEset)
    n.two <- length(TWOset)
    totalM <- matrix(0,n.one+n.two,7)
    for(i in 1:(n.one+n.two)){
        if(i <= n.one){
            P <- check_gene(module,distgene,ONEs1[i],as.character(ONEset[i]),DAWNrank,DAWNgenes,TADAr,DAWNresult)
        }else{
            P <- check_gene(module,distgene,TWOs1[i-n.one],as.character(TWOset[i-n.one]),DAWNrank,DAWNgenes,TADAr,DAWNresult)
        }
        totalM[i,] <- P
    }
    colnames(totalM) <- c("Gene","DAWN","TADA","risk_N","mean_cor","percent_risk > 0.1","cor_significant")
    write.table(totalM,file="DAWN_TADA.txt",quote=F,row.names=F,sep="\t")
    
# i <- 6
# P <- check_onegene(module,distgene,ONEs1[i],as.character(ONEset[i]))
# i<- 7
# P <- check_onegene(module,distgene,ONEs1[i],as.character(ONEset[i]))
# print(P)
# 
# i <- 1
# P <- check_onegene(module,distgene,TWOs1[i],as.character(TWOset[i]))
# print(P)
# i <- 2
# P <- check_onegene(module,distgene,TWOs1[i],as.character(TWOset[i]))
# print(P)


}

check_onegene <- function(module,distgene,gene,gene1){
    modi <- module[module[,1] %in% gene,2]
    modgenes <- module[module[,2]==modi,1]
    moddist <- distgene[modgenes,modgenes]
    subs <- which(colnames(moddist)==gene)
    rcut <- 0.45
    neiborgenes <- modgenes[moddist[subs,] > rcut]

    # disensity plot
    tmp0 <- density(moddist[upper.tri(moddist)])
    tmp <- density(moddist[subs,])
    plot(tmp0,main=paste("Correlation distribution and ",gene1,sep=""),xlab="Correlation",col=1,ylim=c(0,max(c(tmp$y,tmp0$y))+0.2))
    lines(density(moddist[subs,]),col=2)
    abline(v=tmp$x[tmp$y==max(tmp$y)])
    legend("topright",col=1:2,legend=c("ALL genes",gene1),lty=1)
    C <- tmp$x[tmp$y==max(tmp$y)]
    text(C,max(tmp$y)/2, paste(gene1,":",substr(as.character(C),1,5),sep=""), col = 2)
    #print(ks.test(moddist[upper.tri(moddist)],moddist[subs,])$p.value)
    abline(v=0.45,col=3)
    # degree plot 
    n.tmp <- dim(moddist)[1]
    x <- rep(0,n.tmp)
    for (i in 1:n.tmp){
        x[i] <- sum(moddist[i,] > rcut)
    }
    plot(density(x),main="Degree density",xlab="Degree")
    abline(v=length(neiborgenes))
    tmp<- density(x)
    C <- length(neiborgenes)
    text(C,max(tmp$y)/2, paste(gene1,":",C,sep=""), col = 2)

    #
    DAWNresult <- read.csv("D1121result_4_1_0.65_0.45_0.5.csv")
    subfdr <- !is.na(DAWNresult[,"FDR"])
    dawngenes <- DAWNresult[subfdr,2]

    modgenes <- setdiff(modgenes,gene)
    neiborgenes <- setdiff(neiborgenes,gene)

    m <- length(intersect(modgenes,dawngenes))
    k <- length(neiborgenes) 
    n <- length(modgenes) - m
    x <- length(intersect(dawngenes,neiborgenes))
    P <- phyper(x-1, m, n, k, lower.tail = FALSE, log.p = FALSE) 
    P <- P * (n+m)
    P
}

check_gene <- function(module,distgene,gene,gene1,DAWNrank,DAWNgenes,TADAr,DAWNresult){
    modi <- module[module[,1] %in% gene,2]
    modgenes <- module[module[,2]==modi,1]
    moddist <- distgene[modgenes,modgenes]
    subs <- which(colnames(moddist)==gene)
    tmp0 <- density(moddist[upper.tri(moddist)])
    x0 <- tmp0$x[tmp0$y==max(tmp0$y)]
    tmp <- density(moddist[subs,])
    x1 <- tmp$x[tmp$y==max(tmp$y)]
  
    P <- 1:7
    P[7] <- t.test(moddist[upper.tri(moddist)],moddist[subs,])$p.value
    
    P[1] <- gene1
    subs <- which(DAWNgenes %in% gene1)
    n.tmp <- dim(DAWNresult)[2]
    P[2] <- DAWNrank[subs]
    P[3] <- TADAr[subs]
    prisk <- as.numeric(DAWNresult[subs,n.tmp-1])/as.numeric(DAWNresult[subs,n.tmp])
    P[4] <- paste(DAWNresult[subs,n.tmp-1],"(",substr(as.character(prisk),1,5),")",sep="")
    P[5] <- paste(substr(as.character(x1),1,5),"(",substr(as.character(x0),1,5),")",sep="")
    P[6] <- ifelse(prisk>0.1,"Yes","No")
    
    
    P
    
}

check_coex_RBFOX2 <- function(flag=1){
    
    library(WGCNA)
    source("Network_analysis.R")
    mapT <- as.matrix(read.delim("Fmap.txt",header=FALSE,sep="\t"))
    fileflag <- 4
    distflag <- 1
    softPower <- 6

    dataFile <- c("GTEx_heart_expression.txt","GTEx.heart.countstable.txt","PCGC_newBorn_FPKM.txt","fetal_JUNE_AUG_expression.txt","Brain_expression_Ref.txt")
    batchFile <- "renameFileWithHeartTissue1125.list"
    datExpr <- step_1_0_normalize(dataFile[fileflag],batchFile,fileflag)
    distfile <- paste("dist",fileflag,distflag,sep="_")
    load(distfile)
    distgene[is.na(distgene)] <- 0
    genes <- colnames(datExpr)
    dimnames(distgene) <- list(genes,genes)
   
    if(flag == 1){
        MEDissThres = 0.2
        minModuleSize = 200
        adjacency = adjacency(datExpr, type="signed", power = softPower);
    }else if(flag==2){
        MEDissThres = 0.5
        minModuleSize = 30
        adjacency = adjacency(datExpr, type="unsigned", power = softPower);    
    }
    TOM = TOMsimilarity(adjacency);
    dissTOM = 1-TOM
    geneTree = flashClust(as.dist(dissTOM), method = "average");
    dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize);
    dynamicColors = labels2colors(dynamicMods)
    merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
    mergedColors = merge$colors;
    moduleColors = mergedColors
    colorOrder = c("grey", standardColors());
    moduleLabels = match(moduleColors, colorOrder)-1;
    
    MEList = moduleEigengenes(datExpr, colors = moduleColors)
    MEs = MEList$eigengenes
    
    TADAFile <- "TADA_lofmis1202.csv"
    pvalue_data <- read.csv(TADAFile)[,c(1,10)]
    pvalue_data[,1] <- mapping_to(pvalue_data[,1])
    if(fileflag!=2)  pvalue_data[,1] <- mapT[match(pvalue_data[,1],mapT[,2]),1]; 
    pvalue_data <- pvalue_data[!(is.na(pvalue_data[,1]) | is.na(pvalue_data[,2])),]
    pvalue_data <- pvalue_data[pvalue_data[,2]< 1,]
    pvalue_data <- pvalue_data[!duplicated(pvalue_data[,1]),]
    moduleLabels <- delete_module(moduleLabels,n=5000,pvalue_data[,1],genes)
    print(table(moduleLabels))
    check_indir_RBFOX(moduleLabels,distgene,datExpr,MEs,mapT)

}

check_indir_RBFOX <- function(moduleLabels,distgene,datExpr,MEs,mapT){
    subs <- which(colnames(distgene)=="ENSG00000100320")
    #hist(distgene[subs,])
    rbt <- distgene[subs,]
    a1 <-  quantile(rbt,0.95)
    a2 <-  quantile(rbt,0.05)
    rbe <- datExpr[,subs]
    tmp <- rep(0,dim(MEs)[2])
    for(i in 1:dim(MEs)[2]){
        tmp[i] <- cor(MEs[,i],rbe)
    }
    smod0 <- moduleLabels[subs]
    genes <- colnames(datExpr)
    genemod0 <- genes[moduleLabels==smod0]
    #hist(distgene[subs,genemod0])
    smod <- which(tmp < a2 | tmp >a1)
    
    if(length(smod)>0){
        for(i in 1:length(smod)){
            genemod <- genes[moduleLabels==smod[i]]
            geneset1 <- intersect(genemod,genes[distgene[subs,] > a1])
            geneset2 <- intersect(genemod,genes[distgene[subs,] < a2])
            geneset1 <- setdiff(geneset1,genes[subs])
            geneset2 <- setdiff(geneset2,genes[subs])
            geneset1 <- mapT[match(geneset1,mapT[,1]),2]
            geneset2 <- mapT[match(geneset2,mapT[,1]),2]
            geneset1 <- geneset1[!is.na(geneset1)]
            geneset2 <- geneset2[!is.na(geneset2)]
            write.table(geneset1,file=paste("RBupre",i,".txt",sep=""),quote=F,row.names=F,col.names=F)
           # write.table(geneset2,file="RBdownre.txt",quote=F,row.names=F,col.names=F)
    
            grn <- read.delim("enrichment/GRNallnetmap.txt",sep="\t")
            grn <- as.matrix(grn)
            grns <- unique(grn[,2])
            grns2 <- intersect(grns,geneset1)
            geneindir <- unique(grn[grn[,2] %in% grns2,1])
            write.table(geneindir,file=paste("geneindirRBFOX2",i,".txt",sep=""),quote=F,col.names=F,row.names=F)
        }
    }
    
    # downregulatory
#     grns3 <- intersect(grns,geneset2)
#     write.table(grns3,file="DownRBFOX2.txt",quote=F,col.names=F,row.names=F)
#     
#     
#     ##===========
#     davidb <- mapT[match(unique(grn[,1]),mapT[,2]),1]
#     davidb <- davidb[!is.na(davidb)]
#     write.table(davidb,file="backgdavidb.txt",quote=F,col.names=F,row.names=F)
#     
#     module <- read.delim("FPKM4_1_module.txt",sep="\t",header=F)
#     genes <- c("OSBPL8","HSPD1","IRF2","BANP","IGF1R","PTK2")
#     intersect(genes,grn[,2])
    
}