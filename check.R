check <- function(gene,netflag=1,cut=0){
    source("Multi_net.R")
    fileexp <- "PCGCall.txt"
    allnet <- build_net(netflag,fileexp)
    load(paste("module",netflag,sep=""))
   
    mapT <- as.matrix(read.delim("Fmap.txt",header=FALSE,sep="\t"))
    modLab <- setdiff(unique(module[,2]),0)
    n.module <- length(modLab)
    n.genes <- length(module[module[,2]!=0,1])
    tada <- as.matrix(read.table("CRF_PCGC.txt"))
    nodeinfo <- build_node(module,1,"CRF_PCGC.txt")
    cutn=15000
    i <- which(as.numeric(modLab) == as.numeric(module[module[,1]==gene,2]))
    modgenes <- as.vector(module[module[,2]==modLab[i],1])
    net <- discrete_net(allnet,modgenes,netflag,cutn,nodeinfo)
         
    ind1 <- which(net$node==gene)  
    ne1 <- net$node[net$matrix[ind1,]>0]
    se1 <- as.numeric(tada[match(ne1,tada[,1]),2])
    se1[is.na(se1)] <- 0
    s1 <- as.numeric(tada[tada[,1]==gene,2])
    degs1 <- rowSums(net$matrix[match(ne1,net$node),])

    
    #(sum( pmax(1/degs1,1/length(se1)) * (se1[se1>cut] + s1)/2) *s1) / (sum( pmax(1/degs1,1/length(se1)) * (1 - se1[se1>cut] + 1-s1)/2) *(1-s1))
    #sum( pmin(1/degs1,1/length(se1)) * (se1[se1>cut] + s1)/2) *s1
    
}

tmptest <- function(){
    a1 <- matrix(c(0.01,0.005,1-0.01,1-0.005),2,2)
    a2 <- matrix(c(1/20,1/40,19/20,39/40),2,2); a2 <- t(a2)
    a3 <- matrix(c(1/15,1/30,14/15,29/30),2,2)
    1/(abs(a1-a2) + 1)
    1/(abs(a2-a3) + 1)
    
    # ENSG00000183072    NKX2-5 NKX2-5 netflag=1, deg=4; netflag=2, deg=8; netflag=3 deg=18;netflag=4 deg=46;
    #ENSG00000136754    ABI1
    a1 <- which(module[,1]=="ENSG00000038358") # EDC4
    a2 <- which(module[,1]=="ENSG00000167548") # KMT2D
    
    ind1 <- which(net$node=="ENSG00000038358") # EDC4 110
    ind2 <- which(net$node=="ENSG00000167548")# KMT2D 15; netflag=4 deg=3
    ind3 <- which(net$node=="ENSG00000175634")# RPS6KB2 2
    ind4 <- which(net$node=="ENSG00000153162")# BMP6 0
    # ENSG00000100320 # RBFOX2 0
    # ENSG00000185015 # CA13 92
    # ENSG00000005379 # BZRAP1 69
    # ENSG00000175387 # SMAD2 52
    # ENSG00000124151   # NCOA3 5
    # ENSG00000100393  #  EP300 41
    # ENSG00000108854    #SMURF2 75
    # ENSG00000100393    #EP300
    # ENSG00000065883    #CDK13
    # ENSG00000175387    #SMAD2
    #ENSG00000150991    #UBC
    #ENSG00000070831    #CDC42
    #ENSG00000036257 #CUL3
    tada <- as.matrix(read.table("TADAinfo0121.txt"))
    ne1 <- net$node[net$matrix[ind1,]>0]
    ne2 <- net$node[net$matrix[ind2,]>0]
    ne3 <- net$node[net$matrix[ind3,]>0]
    se1 <- as.numeric(tada[match(ne1,tada[,1]),2])
    se2 <- as.numeric(tada[match(ne2,tada[,1]),2])
    se3 <- as.numeric(tada[match(ne3,tada[,1]),2])
    
    se1[is.na(se1)] <- 0
    se2[is.na(se2)] <- 0
    se3[is.na(se3)] <- 0
    
    s1 <- as.numeric(tada[tada[,1]=="ENSG00000038358",2])
    s2 <- as.numeric(tada[tada[,1]=="ENSG00000167548",2])
    s3 <- as.numeric(tada[tada[,1]=="ENSG00000175634",2])
    
    cut <- 0.1
    degs1 <- rowSums(net$matrix[match(ne1[se1>cut],net$node),])
    degs2 <- rowSums(net$matrix[match(ne2[se2>cut],net$node),])
    
    (sum( pmax(1/degs1,1/length(se1)) * (se1[se1>cut] + s1)/2) *s1) / (sum( pmax(1/degs1,1/length(se1)) * (1 - se1[se1>cut] + 1-s1)/2) *(1-s1))
    (sum( pmax(1/degs2,1/length(se2)) * (se2[se2>cut] + s2)/2) *s2) / (sum( pmax(1/degs2,1/length(se2)) * (1 - se2[se2>cut] + 1-s2)/2) *(1-s2))
    
    sum( pmin(1/degs1,1/length(se1)) * (se1[se1>cut] + s1)/2) *s1
    sum( pmin(1/degs2,1/length(se2)) * (se2[se2>cut] + s2)/2) *s2

}

check_DAWN_TADA <- function(){
    source("Network_analysis.R")
    
    mapT <- as.matrix(read.delim("Fmap0121.txt",header=FALSE,sep="\t"))
    # compare the DAWN rank and TADA rank
    DAWNFile <- "PCGCresult.csv"
    # only compare DAWN FDR<1 genes
    DAWNresult <- as.matrix(read.csv(DAWNFile))
    
    ntop <- 500
    DAWNfdr <- DAWNresult[1:ntop,"DAWN.FDR"]
    DAWNgenes <- DAWNresult[1:ntop,"Gene"]
    
    
    DAWNsub <- sort(DAWNfdr,index.return=TRUE)$ix
    DAWNrank <- 1:ntop
    DAWNrank[DAWNsub] <- 1:ntop
    
    TADAFile <- "TADAinfo0130.txt"
    genes <- mapT[match(DAWNgenes,mapT[,2]),1]
    TADAresult <- as.matrix(read.table(TADAFile))
    TADAr <- match(genes,TADAresult[,1])
    
    subsone <- which(TADAr/DAWNrank > 1.5 & TADAr <= ntop & TADAr > 30)
    substwo <- which(DAWNrank/TADAr > 1.5 & DAWNrank <= ntop & DAWNrank > 30)
    ONEset <- DAWNgenes[subsone]
    TWOset <- DAWNgenes[substwo]
        
    plot(DAWNrank,TADAr,main="DAWN rank VS TADA rank",xlab="DAWN rank",ylab="TADA rank",xlim=c(1,ntop),ylim=c(1,ntop))
    abline(0,1,col=1,lwd=1)
    abline(0,1.5,col=2,lwd=1)
    abline(0,2/3,col=2,lwd=1)
    
    subs <- which(DAWNresult[,"Gene"] %in% as.character(ONEset))
    lines(DAWNrank[subs],TADAr[subs],type="p",col="red",lwd=3)
    subs <- which(DAWNresult[,"Gene"] %in% as.character(TWOset))
    lines(DAWNrank[subs],TADAr[subs],type="p",col="red",lwd=3)
    
    #subs <- which(DAWNgenes=="ENSG00000100393")
    subs <- which(DAWNgenes=="EP300")
    lines(DAWNrank[subs],TADAr[subs],type="p",col="blue",lwd=6)
   
    print(length(ONEset))
    print(length(TWOset))
}

check_HotNet2_TADA <- function(){

source("Network_analysis.R")
mapT <- as.matrix(read.delim("Fmap0121.txt",header=FALSE,sep="\t"))
#DAWNFile <- "result/RWR_8_resultall.txt"
DAWNFile <- "result/nonNDD0210_50c_3LBP_31.txt"
DAWNresult <- as.matrix(read.table(DAWNFile))

ntop <- 500
DAWNgenes <- DAWNresult[1:ntop,1]
DAWNrank <- 1:ntop

TADAFile <- "TADAinfo0210nonNDD.txt"
genes <- mapT[match(DAWNgenes,mapT[,2]),1]
TADAresult <- as.matrix(read.table(TADAFile))
TADAr <- match(genes,TADAresult[,1])

subsone <- which(TADAr/DAWNrank > 1.5 & TADAr <= ntop & TADAr > 30)
substwo <- which(DAWNrank/TADAr > 1.5 & DAWNrank <= ntop & DAWNrank > 30)
ONEset <- DAWNgenes[subsone]
TWOset <- DAWNgenes[substwo]

plot(DAWNrank,TADAr,main="HotNet2 rank VS TADA rank",xlab="HotNet2 rank",ylab="TADA rank",xlim=c(1,ntop),ylim=c(1,ntop))
abline(0,1,col=1,lwd=1)
abline(0,1.5,col=2,lwd=1)
abline(0,2/3,col=2,lwd=1)

subs <- which(DAWNresult[,1] %in% as.character(ONEset))
lines(DAWNrank[subs],TADAr[subs],type="p",col="red",lwd=2)
subs <- which(DAWNresult[,1] %in% as.character(TWOset))
lines(DAWNrank[subs],TADAr[subs],type="p",col="red",lwd=2)

#subs <- which(DAWNgenes=="ENSG00000100393")
subs <- which(DAWNgenes=="EP300")
lines(DAWNrank[subs],TADAr[subs],type="p",col="blue",lwd=6,pch=19)

print(length(ONEset))
print(length(TWOset))

}