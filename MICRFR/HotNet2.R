HotNet2 <- function(){
    # generate heat
    TADAFile  <- "TADA_resultsPCGC_2_4.csv"
    geneinfo <- read.csv(TADAFile)
    pi0 <- 0.94 #!!!!!!
    pi <- 1-pi0
    posP <- (geneinfo[,"BF"]*pi)/(geneinfo[,"BF"]*pi + 1-pi)
    write.table(cbind(geneinfo[,1],posP),file="HotNetheat.txt",quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")
    
    genes <- read.table("hotnet/iref_index_genes")
    edges <- read.table("hotnet/iref_edge_list")
    net <- edges
    net[,1] <- genes[match(edges[,1] , genes[,1]),2]
    net[,2] <- genes[match(edges[,2] , genes[,1]),2]
    write.table(net,file="hotnet/iRefIndex.txt",quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")
    
    net0 <- mapping_net(net)
    write.table(net0,file="hotnet/iRefIndexmap.txt",quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")
      
    library(R.matlab)
    genes <- read.table("ASD/iref_index_genes")[,2]
    Ma <- readMat("ASD/iref_ppr_0.55.mat")
    PPR <- Ma$PPR
    save(PPR,file="iref_PPR")
    
    source("HotNet2.R")
    load("iref_PPR")
    genes <- read.table("ASD/iref_index_genes")[,2]
    strname <- c("ASDall","ASD13772","ASD13908","ASD932")
    filenames <- c("ASD/TADAresult/hotnet_inputall.txt","ASD/TADAresult/hotnet_input13772.txt","ASD/TADAresult/hotnet_input13908.txt","ASD/TADAresult/hotnet_input932.txt")
    for(i in 1:length(filenames)){
        mutp <- read.table(filenames[i])
        p1 <- matrix(0,length(genes),1,dimnames=list(genes,1))
        g1 <- intersect(mutp[,1],genes)
        p1[match(g1,genes)] <- mutp[match(g1,mutp[,1]),2]
        p1new <- t(p1) %*% t(PPR)
        p1new <- matrix(c(genes,p1new),length(genes),2)
        p1new <- p1new[order(-as.numeric(p1new[,2])),]
        filename <- paste("result/hotnet_result_",strname[i],".txt",sep="")
        write.table(p1new,file=filename,quote=FALSE,col.names=FALSE,sep="\t",row.names=FALSE)
        addFDR1(filename,strname[i],"result/")
    }   

    source("HotNet2.R")
    load("iref_PPR")
    genes <- read.table("ASD/iref_index_genes")[,2]
    filenames <- c("hotnet_PCGC.txt","random_samples/hotnet_rand1.txt","random_samples/hotnet_rand2.txt","random_samples/hotnet_rand3.txt")
    strname <- c("PCGC","rand1","rand2","rand3")
    HotNet2_1(PPR,genes,filenames,strname)
    
    
    source("Network_analysis.R")
    source("HotNet2.R")
    load("iref_PPR")
    genes <- read.table("ASD/iref_index_genes")[,2]
    genes <- mapping_to(genes)
    
    dirstr <- "result/PCGC/"
    dirstr1 <- "result/PCGC/"
    strname <- "PCGC"
    filename <- paste(dirstr,"hotnet_input",strname,".txt",sep="")
    HotNet2_1(PPR,genes,filename,strname,dirstr1,flag=2)
    
    HotNet2R(1,strname,paste(dirstr,"CRF_input",strname,".txt",sep=""),"PCGCall.txt",dirstr,cutn=0.50,idm=TRUE)
    HotNet2R(3,strname,paste(dirstr,"CRF_input",strname,".txt",sep=""),"PCGCall.txt",dirstr,cutn=0.50,idm=TRUE)
    HotNet2R(8,strname,paste(dirstr,"hotnet_input",strname,".txt",sep=""),"PCGCall.txt",dirstr,cutn=0.50,idm=FALSE)
    source("Multi_net.R")
    source("CRF_build.R")
    for(netflag in c(9,10,11,12,13)){
        if(netflag==3 | netflag==1 | netflag>10){ instr <- "CRF_input";idm=TRUE; }else{ instr <- "hotnet_input"; idm=FALSE;}
        HotNet2R(netflag,strname,paste(dirstr,instr,strname,".txt",sep=""),"PCGCall.txt",dirstr,cutn=0.50,idm=idm)
    }
 
    
    
    source("Network_analysis.R")
    source("HotNet2.R")
    load("iref_PPR")
    genes <- read.table("ASD/iref_index_genes")[,2]
    genes <- mapping_to(genes)
    
    dirstr <- "result/"
    dirstr1 <- "result/"
    strname <- "control"
    filename <- paste(dirstr,"hotnet_input",strname,".txt",sep="")
    HotNet2_1(PPR,genes,filename,strname,dirstr1,flag=2)
    

    HotNet2R(3,strname,paste(dirstr,"CRF_input",strname,".txt",sep=""),"PCGCall.txt",dirstr,cutn=0.50,idm=TRUE)
    HotNet2R(7,strname,paste(dirstr,"hotnet_input",strname,".txt",sep=""),"PCGCall.txt",dirstr,cutn=0.50,idm=FALSE)
    HotNet2R(8,strname,paste(dirstr,"hotnet_input",strname,".txt",sep=""),"PCGCall.txt",dirstr,cutn=0.50,idm=FALSE)
    
}

HotNet2_randset <- function(){
    source("Network_analysis.R")
    source("HotNet2.R")
    load("iref_PPR")
    genes <- read.table("ASD/iref_index_genes")[,2]
    genes <- mapping_to(genes)
    
    dirstr <- "result/randset/"
    dirstr1 <- "result/randresult/iRef/"
    
    strname <- "ASDall"
    filename <- paste(dirstr,"hotnet_input",strname,".txt",sep="")
    HotNet2_1(PPR,genes,filename,strname,dirstr1,flag=2)
    
    ## different sample size and exclude sample sets
    for(j in 2:9){
        for(i in 1:5){
            strname <- paste("part",j,"_",i,sep="")
            filename <- paste(dirstr,"hotnet_input",strname,".txt",sep="")
            HotNet2_1(PPR,genes,filename,strname,dirstr1,flag=2)
            
            strname <- paste("rest",j,"_",i,sep="")
            filename <- paste(dirstr,"hotnet_input",strname,".txt",sep="")
            HotNet2_1(PPR,genes,filename,strname,dirstr1,flag=2)           
        }
    }

dirstr <- "result/randset2/"
dirstr1 <- "result/randresult2/iRef/"

## different sample size and more samples
for(j in 1:3){
    for(i in 1:5){
        strname <- paste("part",j,"_",i,sep="")
        filename <- paste(dirstr,"hotnet_input",strname,".txt",sep="")
        HotNet2_1(PPR,genes,filename,strname,dirstr1,flag=2)         
    }
}

}

HotNet2_randset2 <- function(){
    source("Network_analysis.R")
    source("HotNet2.R")
    load("iref_PPR")
    genes <- read.table("ASD/iref_index_genes")[,2]
    genes <- mapping_to(genes)
    
    dirstr <- "result/leaveone3/"
    dirstr1 <- "result/leaveone3result/iRef/"

    ## leave one genes with FDR<=0.2
    for(i in 1:251){
            strname <- paste("rand2_de1_",i,sep="")
            filename <- paste(dirstr,"hotnet_input",strname,".txt",sep="")
            HotNet2_1(PPR,genes,filename,strname,dirstr1,flag=2)
    }
    
}

HotNet2_1 <- function(PPR,genes,filenames,strname,dirstr="result/",flag=1){
    
    for(i in 1:length(filenames)){
        mutp <- read.table(filenames[i])
        p1 <- matrix(0,length(genes),1,dimnames=list(genes,1))
        g1 <- intersect(mutp[,1],genes)
        p1[match(g1,genes)] <- mutp[match(g1,mutp[,1]),2]
        p1new <- t(p1) %*% t(PPR)
        p1new <- matrix(c(genes,p1new),length(genes),2)
        p1new <- p1new[order(-as.numeric(p1new[,2])),]
        filename <- paste(dirstr,"hotnet_result_",strname[i],".txt",sep="")
        write.table(p1new,file=filename,quote=FALSE,col.names=FALSE,sep="\t",row.names=FALSE)
        addFDR3(filename,strname[i],dirstr,flag)
    }
}

batch_hotnet2 <- function(){
    
    source("CRF_build.R")
    source("HotNet2.R")
    source("Network_analysis.R")
    source("Multi_net.R")
    options(stringsAsFactors=FALSE)
    fileexp=""
    netstr <- c("STRING/","iRef/","coexp/","Infer/","GNAT/")
    k = 1
    for(netflag in c(3,6,7,8,21)){
        
        if(netflag==3){ instr <- "CRF_input";idm=TRUE; }else{ instr <- "hotnet_input"; idm=FALSE;}
                         
        dirin <- "result/randset_1/"
        dirstr <- paste("result/randresult_1/",netstr[k],sep="")
    
        strname <- "ASD4_16"
        filename <- paste(dirin,instr,strname,".txt",sep="")
        HotNet2R(netflag,strname,filename,fileexp,dirstr,cutn=0.50,idm)
        
        
        ## random set 1: different sample size and exclude sample sets
        for(j in 2:9){
            for(i in 1:20){
                strname <- paste("part",j,"_",i,sep="")
                filename <- paste(dirin,instr,strname,".txt",sep="")
                HotNet2R(netflag,strname,filename,fileexp,dirstr,cutn=0.50,idm)
            }
        }
        
        ### random set 2: leaveone mutation
        dirin <- "result/leaveone4/"
        dirstr <- paste("result/leaveone4result/",netstr[k],sep="")
    
        for(i in 1:187){
            strname <- paste("rand2_1_",i,sep="")
            filename <- paste(dirin,instr,strname,".txt",sep="")
            HotNet2R(netflag,strname,filename,fileexp,dirstr,cutn=0.50,idm)
        }
  
     
        
        ## random set 3: rand set 4
        dirin <- "result/randset4/"
        dirstr <- paste("result/randresult4/",netstr[k],sep="")
        
        ## random set 3: different sample size and exclude sample sets
        for(j in 3){
            for(i in 1:20){
                strname <- paste("part",j,"_",i,sep="")
                filename <- paste(dirin,instr,strname,".txt",sep="")
                HotNet2R(netflag,strname,filename,fileexp,dirstr,cutn=0.50,idm)
            }
        }
        
        
        k <- k + 1
        
    }
    
}

HotNet2R <- function(netflag,strn,filenames,fileexp,dirstr,cutn,idm=FALSE){
    
    allnet <- build_net(netflag,fileexp)
    allnet$matrix[allnet$matrix > 0] <- 1
    
    nodeinfo <- as.matrix(read.table(filenames[1]))
    POSpro <- build_rwr(nodeinfo,allnet,netflag,strn,dirstr,cutn)
    
    if(idm){
        mapT <- as.matrix(read.delim("Fmap0121.txt",header=FALSE,sep="\t"))
        pros1 <- POSpro
        pros1[,1] <- mapT[match(pros1[,1],mapT[,1]),2]
        pros1 <- pros1[!is.na(pros1[,1]),]
        write.table(pros1,file=paste(dirstr,"hotnetresult1",strn,netflag,"1.txt",sep=""),quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")
    }
    
}

build_rwr <- function(nodeinfo,allnet,netflag,strn,dirstr,cutn=0.5){

    alpha <- cutn
    net <- allnet#build_net1(allnet,netflag)
    # column normalized
    W <- net$matrix
    W <- W/matrix(colSums(W),net$size,net$size,byrow=TRUE)
    
    # run the RWR algorithm
    p0 <- matrix(0,net$size,1,dimnames=list(net$node,1))
    tmpg <- intersect(nodeinfo[,1],net$node)
    p0[match(tmpg,net$node)] <- as.numeric(nodeinfo[match(tmpg,nodeinfo[,1]),2]) 
    eng <- sum(p0)
    p0 <- p0/eng
    d <- p0
    p1 <- RWR(W,p0,alpha,d)
    p1 <- p1*eng  
    p2 <- p1
    p1 <- fit_mix_nor(p1,k=1)
    
    POSpro <- matrix(0,length(p1),5)
    POSpro[,1] <- net$node
    POSpro[,2] <- p2
    FDR <- estimate.FDR(p1)
    POSpro[,3] <- FDR
    POSpro[,5] <- p.adjust(p1,method="fdr")
    POSpro[,4] <- 1 - as.numeric(POSpro[,5])
    
    POSpro <- POSpro[order(as.numeric(POSpro[,2]),decreasing = TRUE),]
    write.table(POSpro,file=paste(dirstr,"hotnetresult1",strn,netflag,".txt",sep=""),quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")
    
    POSpro
}

build_net1 <- function(allnet,netflag){
    
    net <- allnet
    if(netflag==1){
        cutf <- quantile(net$matrix[upper.tri(net$matrix)],0.95)
        net$matrix[net$matrix >= cutf] <- 1
        net$matrix[net$matrix < cutf] <- 0
    }
    
    if(netflag==2){
        neta <- minet::aracne(net$matrix)
        net$matrix <- neta
        net$matrix[net$matrix>0] <- 1
    } 
    
    if(netflag==3 | netflag==4 | netflag==0 | netflag==5){
        net$matrix[net$matrix > 0] <- 1
    }
    
    # eliminate isolated nodes
    tmp <- allnet$matrix
    mindeg <- 2
    subs <- which(colSums(net$matrix) == 0)
    if(length(subs)>0){
        for (i in 1:length(subs)){
            net$matrix[subs[i],sort(tmp[subs[i],],decreasing=TRUE,index.return=TRUE)$ix[1:mindeg]] <- 1
            net$matrix[sort(tmp[subs[i],],decreasing=TRUE,index.return=TRUE)$ix[1:mindeg],subs[i]] <- 1
        }
    }
    
    net
}

RWR <- function(W,p0,alpha,d){
    # RWR: random walk with restart 
    
    # parameters:
    # W: the edge weights matrix
    # p0 : the start point for iteration
    # alpha: the restart probability
    # d: starting vector
    
    p1 <- (1-alpha)*W%*%p0 + alpha*d;
    iter <- 1
    while (sum(abs(p0-p1))>10^(-10)){
        p0=p1;
        p1 <- (1-alpha)*W%*%p0 + alpha*d;
        iter <- iter + 1
        #print(iter)
    }
    
    #print(iter)
    p1
    
}

mapping_net <- function(net){
    source("Network_analysis.R")
    mapT <- as.matrix(read.delim("Fmap0121.txt",header=FALSE,sep="\t"))
    genes <- union(net[,1],net[,2])
    genesm <- mapping_to(genes)
    genemap <- cbind(genes,genesm)
    
    net0 <- net
    net0[,1] <- genemap[match(net[,1],genes),2]
    net0[,2] <- genemap[match(net[,2],genes),2]
    
    net0[,1] <- mapT[match(net0[,1],mapT[,2]),1]
    net0[,2] <- mapT[match(net0[,2],mapT[,2]),1]
    
    net0 <- net0[!is.na(net0[,1]) & !is.na(net0[,2]),]
    net0
}

build_module <- function(){
    
    ## gene id mapping here 
    source("Network_analysis.R")
    genes <- read.table("hotnet/iref_index_genes")
    genes[,2] <- mapping_to(genes[,2])
    
    edges <- read.table("hotnet/iref_edge_list")
    net <- edges
    net[,1] <- genes[match(edges[,1] , genes[,1]),2]
    net[,2] <- genes[match(edges[,2] , genes[,1]),2]
    write.table(net,file="hotnet/iRefIndexm.txt",quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")
    
    source("CRF_build.R")
    net.text <- as.matrix(read.table("hotnet/iRefIndexm.txt",sep="\t",header=FALSE))
    net.text <- rbind(net.text,net.text[,c(2,1,3)])
    net <- read_net(net.text)
    distgene <- net$matrix
    library(WGCNA)
    softPower=1
    minModuleSize = 200
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
    module[module[,2]==0,2] <- 10
    
    save(module,file="iRefmodule")
    
}

addFDR <- function(){
    mapT <- as.matrix(read.delim("Fmap0121.txt",header=FALSE,sep="\t"))
    
    filenames <- c("ASD/hotnet2_13772.txt","ASD/hotnet2_13908.txt","ASD/hotnet2_932.txt","ASD/hotnet2_932p.txt")
    #genes <- read.table("ASD/iref_index_genes")[,2]
    netflag=6
    strn0="result/ASD_hotnet"
    strname <- c("13772","13908","932","932p")
    for(i in 1:length(filenames)){
        strn <- paste(strn0,strname[i],sep="")
        tescore <- read.table(filenames[i])
        tescore[,1] <- mapping_to(tescore[,1])
        p1 <- as.numeric(tescore[,2])
        n <- length(p1)
        p1 <- fit_mix_nor(p1,k=1)
        
        POSpro <- matrix(0,n,4)
        POSpro[,1] <- tescore[,1]
        POSpro[,3] <- pmin(p1*n,1)
        POSpro[,2] <- 1 - as.numeric(POSpro[,3])
        
        
        POSpro <- POSpro[order(as.numeric(POSpro[,2]),decreasing = TRUE),]
        POSpro[,4] <- as.numeric(POSpro[,3]) * dim(POSpro)[1]
        POSpro[as.numeric(POSpro[,4])>1,4] <- 1
        b <- as.numeric(POSpro[,3])
        BF <- as.numeric(POSpro[,2])/b
        BF[BF==Inf] <- 2*max(BF[BF<Inf])
        FDR <- Bayesian.FDR(BF)$FDR ####!!!!!!!
        POSpro <- cbind(POSpro,FDR)
        POSpro <- POSpro[order(as.numeric(POSpro[,5]),decreasing=FALSE),]
        
        pros1 <- POSpro
        pros1[,1] <- mapT[match(pros1[,1],mapT[,2]),1]
        pros1 <- pros1[!is.na(pros1[,1]),]
        write.table(POSpro,file=paste(strn,"LBP_",netflag,"1.txt",sep=""),quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")
        write.table(pros1,file=paste(strn,"LBP_",netflag,".txt",sep=""),quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")
    }
    
}

addFDR1 <- function(filename,str,dirstr,flag){
    source("Network_analysis.R")
    source("Multi_net.R")
    library(fdrtool)
    
    mapT <- as.matrix(read.delim("Fmap0121.txt",header=FALSE,sep="\t"))
    netflag=6
    strn <- str
    tescore <- read.table(filename)
    if(flag==1){
        tescore[,1] <- mapping_to(tescore[,1])
    }
    p1 <- as.numeric(tescore[,2])
    n <- length(p1)
    p0 <- p1
    p1 <- fit_mix_nor(p1,k=1)
    
    POSpro <- matrix(0,n,4)
    POSpro[,1] <- tescore[,1]
    POSpro[,2] <- p0
    POSpro[,4] <- pmin(p1*n,1)
    POSpro[,3] <- 1 - as.numeric(POSpro[,4])
    
    POSpro <- POSpro[order(as.numeric(POSpro[,2]),decreasing = TRUE),]
    b <- as.numeric(POSpro[,4])
    BF <- as.numeric(POSpro[,3])/b
    BF[BF==Inf] <- 2*max(BF[BF<Inf])
    FDR <- Bayesian.FDR(BF)$FDR ####!!!!!!!
    POSpro <- cbind(POSpro,FDR)
    #POSpro <- POSpro[order(as.numeric(POSpro[,2]),decreasing=TRUE),]
    
    pros1 <- POSpro
    pros1[,1] <- mapT[match(pros1[,1],mapT[,2]),1]
    pros1 <- pros1[!is.na(pros1[,1]),]
    write.table(POSpro,file=paste(dirstr,"hotnetresult1",strn,netflag,".txt",sep=""),quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")
    #write.table(pros1,file=paste(dirstr,"hotnetresult",strn,netflag,".txt",sep=""),quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")
    
}

addFDR2 <- function(filename,str,dirstr){
    source("Network_analysis.R")
    source("Multi_net.R")
    library(fdrtool)
    
    mapT <- as.matrix(read.delim("Fmap0121.txt",header=FALSE,sep="\t"))
    netflag=6
    strn <- paste(dirstr,str,sep="")
    tescore <- read.table(filename)
    tescore[,1] <- mapping_to(tescore[,1])
    p1 <- as.numeric(tescore[,2])
    n <- length(p1)
    fdrest <- fdrtool((p1-mean(p1))/sd(p1),statistic="normal",plot=FALSE)
    
    POSpro <- matrix(0,n,5)
    POSpro[,1] <- tescore[,1]
    POSpro[,2] <- p1
    POSpro[,3] <- fdrest$pval
    POSpro[,4] <- as.numeric(POSpro[,3]) * dim(POSpro)[1]
    POSpro[as.numeric(POSpro[,4])>1,4] <- 1
    POSpro[,5] <- fdrest$qval
    POSpro <- POSpro[order(as.numeric(POSpro[,5]),decreasing=FALSE),]
    
    pros1 <- POSpro
    pros1[,1] <- mapT[match(pros1[,1],mapT[,2]),1]
    pros1 <- pros1[!is.na(pros1[,1]),]
    write.table(POSpro,file=paste("hotnetresult1",strn,netflag,".txt",sep=""),quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")
    #write.table(pros1,file=paste("hotnetresult",strn,netflag,".txt",sep=""),quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")
    
}

addFDR3 <- function(filename,str,dirstr,flag,idm=FALSE){
    source("Network_analysis.R")
    source("Multi_net.R")
    
    mapT <- as.matrix(read.delim("Fmap0121.txt",header=FALSE,sep="\t"))
    netflag=6
    strn <- str
    tescore <- read.table(filename)
    if(flag==1){
        tescore[,1] <- mapping_to(tescore[,1])
    }
    p1 <- as.numeric(tescore[,2])
    n <- length(p1)
    p0 <- p1
    p1 <- fit_mix_nor(p1,k=1)
    
    POSpro <- matrix(0,n,5)
    POSpro[,1] <- tescore[,1]
    POSpro[,2] <- p0
    FDR <- estimate.FDR(p1)
    POSpro[,3] <- FDR
    POSpro[,5] <- p.adjust(p1,method="fdr")
    POSpro[,4] <- 1 - as.numeric(POSpro[,5])
    
    POSpro <- POSpro[order(as.numeric(POSpro[,2]),decreasing = TRUE),]
    write.table(POSpro,file=paste(dirstr,"hotnetresult1",strn,netflag,".txt",sep=""),quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")
    
    if(idm){
        pros1 <- POSpro
        pros1[,1] <- mapT[match(pros1[,1],mapT[,2]),1]
        pros1 <- pros1[!is.na(pros1[,1]),]
    }
    
    #POSpro
}

HotNet2Rpa <- function(){
    source("CRF_build.R")
    source("Multi_net1.R")
    source("Network_analysis.R")
    fileexp <- "PCGCall.txt"
    for(netflag in 1){
        allnet <- build_net(netflag,fileexp) 
        net <- build_net1(allnet,netflag)
        W <- net$matrix
        W <- W/matrix(colSums(W),net$size,net$size,byrow=TRUE)
        
        ps <- seq(0.1,0.99,0.1)
        pa1 <- matrix(0,net$size,length(ps))
        pa2 <- matrix(0,net$size,length(ps))
        k <- 1
        for(alpha in ps){
            Fm <- matrix(0,net$size,net$size)
            diag(Fm) <- 1
            P <- (1-alpha)*W
            while(max(P) >= 10e-10){
                Fm <- Fm + P
                P <- P%*%P
            }
            Fm <- alpha*Fm
            
            vals <- heatdis(Fm,W)
            pa1[,k] <- vals[,2]/vals[,3]
            pa2[,k] <- vals[,2]/(vals[,1] + vals[,3])
            print(k)
            k <- k + 1
            write.csv(vals,file=paste("heat",netflag,alpha,".csv",sep=""),row.names=FALSE)
        }
        print(sort(apply(pa1,2,mean))[length(ps)])
        print(ps[sort(apply(pa1,2,mean),index.return=TRUE)$ix[length(ps)]])
        
        print(sort(apply(pa2,2,mean))[length(ps)])
        print(ps[sort(apply(pa2,2,mean),index.return=TRUE)$ix[length(ps)]])
    }
}

heatdis <- function(Fm,W){
    
    vals <- matrix(0,dim(W)[1],3)
    vals[,1] <- diag(Fm)
    diag(Fm) <- 0
    a1 <- matrix(0,dim(W)[1],dim(W)[2])
    a2 <- a1
    a1[W>0] <- Fm[W>0]
    a2[W==0] <- Fm[W==0]
    vals[,2] <- rowSums(a1) # i --> its neighbors  #colSums(a1) j --> i
    vals[,3] <- rowSums(a2) # i --> not neighbors  #colSums(a2) j --> i
    
    vals
}
