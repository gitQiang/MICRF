Multi_net <- function(netflag,strn="LBP",filenames,fileexp){
    #netflag <- 1 # 1: expression; 2: regulatory ARACNE; 3: PPI; 4: Phenotype based
    #filenames <- c("TADAinfo.txt","CNVinfo.txt","otherinfo.txt")
    
    allnet <- build_net(netflag,fileexp)
    if(netflag ==1 | netflag==2){
        if(fileexp=="Brain_expression_Ref.txt"){
            load(paste("module_brain",1,sep=""))
        }else{
            load(paste("module",netflag,sep=""))
        }
    }else{
        load(paste("module",netflag,sep=""))
    }
    #module <- build_module(allnet$matrix)
    nodeinfo <- build_node(module,1,filenames)
    build_CRF(nodeinfo,allnet,module,netflag,strn,cutn=15000)
}

batch_Multi_net <- function(){
    source("CRF_build.R")
    source("MixNet.R")
    source("Network_analysis.R")
    source("Multi_net.R")
    options(stringsAsFactors=FALSE)
    
    strn0 <- "result/Newtest" #strn: the label for different data set names
    fnode <- "NewtestPheno_score_map.txt"
    #fnode <- c("TADAinfo0113.txt","CNVinfo.txt","otherinfo.txt")
    #fnode <- c("random_samples/CHDnaturePheno_score_map.txt","CNVinfo.txt","otherinfo.txt")
    fileexp <- "PCGCall.txt" 
    #fileexp <- "Brain_expression_Ref.txt"
    filenames <- c()
    ngene=300
    for(i in 1:4){
        strn <- paste(strn0,i,sep="_")
        Multi_net(i,strn,fnode,fileexp);
        rank_combine1(paste(strn,"LBP_",i,"1.txt",sep=""),strn,i,ngene)
        filenames <- union(filenames,paste(strn,"LBP_",i,".txt",sep=""))
    }
    flag=6
    rank_combine(filenames,strn0,flag,ngene)
    
    a <- read.table(paste(strn0,flag,"resultall.txt",sep="_"))
    testgenes <- c("RBFOX2","KMT2D","EP300","NAA16","KDM5B","NCOA3","SMAD2","PTPN11","KDM5B","NAA15","POGZ","AHNAK","MYH6","NOTCH1","SUV420H1","MASTL","CHD7","BRAF","SMARCC1","CDKL2","JAG1","FLT4","PCNXL2","MKRN2","MROH7","TLK2","ARID1B","CTNNB1","ETS1","CUL3")
    sum(match(testgenes,a[,1])<422)
    print(match(testgenes,a[,1]))
    a1 <- match(testgenes,a[,1])
}

build_node <- function(module,flag=1,filenames=c("TADAinfo.txt","CNVinfo.txt","otherinfo.txt")){

    mapT <- as.matrix(read.delim("Fmap.txt",header=FALSE,sep="\t")) ###!!! 可以进一步精细
    genes <- module[,1]
    mapT <- mapT[mapT[,1] %in% genes,]
    
    if(flag==1){
        nodeinfo <- matrix(0,length(genes),1,dimnames=list(genes,1))
        #nodeinfo[,1] <- 0.1 ## related genes percent estimate
        tmp <- read.table(filenames[1])
        tmp <- tmp[tmp[,1] %in% genes,]
        nodeinfo[match(tmp[,1],rownames(nodeinfo)),1] <- tmp[,2]
    }else{
        n.info <- length(filenames)
        nodeinfo <- matrix(0,length(genes),n.info,dimnames=list(genes,1:n.info))
        #nodeinfo[,1] <- 0.1 ## related genes percent estimate
        #nodeinfo[,2:n.info] <- 0
        #nodeinfo[match(names(nodesim),rownames(nodeinfo)),1] <- nodesim
        for(i in 1:n.info){
            tmp <- read.table(filenames[i])
            tmp <- tmp[tmp[,1] %in% genes,]
            nodeinfo[match(tmp[,1],rownames(nodeinfo)),i] <- tmp[,2]
        }
    }
    nodeinfo
}

build_net <- function(netflag,fileexp="PCGCall.txt"){
    # net: a list: size node matrix (systemtic and zero diag)
    net <- list()
    if(netflag==1 | netflag==2){
        datExpr <- as.matrix(read.delim(fileexp,sep="\t",row.names=1))
        datExpr <- t(datExpr)
        if(netflag==1){library(WGCNA); distM = adjacency(datExpr, type = "unsigned", power=1);}
        if(netflag==2){library(minet); distM <- build.mim(datExpr,estimator="spearman"); distM <- distM/max(distM);}
        genes <- colnames(datExpr)
        net$matrix <- distM
        diag(net$matrix) <- 0
        net$size <- length(genes)
        net$node <- genes
    }else if(netflag==3 | netflag==4){
        if(netflag==3){ filename <- "STRINGnetmap.txt";}
        if(netflag==4){ filename <- "phenotype_net.txt";}
        net.text <- as.matrix(read.table(filename,sep="\t",header=FALSE))
        net.text <- rbind(net.text,net.text[,c(2,1,3)])
        net.text[,3] <- as.numeric(net.text[,3])/max(as.numeric(net.text[,3]))
        net <- read_net(net.text)
        diag(net$matrix) <- 0
    }else if(netflag==0){
        load("net0")
        net <- net0
    }
    
    net
}

build_module <- function(distgene){
    library(WGCNA)
    powers = c(c(1:10), seq(from = 12, to=20, by=2))
    sft = pickSoftThreshold.fromSimilarity(distgene, powerVector = powers, verbose = 5)
    softPower <- sft$fitIndices[sft$fitIndices[,2]==max(sft$fitIndices[,2]),1] #sft$fitIndices[,1],-sign(sft$fitIndices[,3]*sft$fitIndices[,2]), #sft$fitIndices[,5]

    # WGCNA to module relationships
    MEDissThres = quantile(distgene[upper.tri(distgene)],0.2)
    minModuleSize = 200; # set the minimum module size relatively high:
    diag(distgene) <- 1
    adjacency = adjacency.fromSimilarity(distgene, type = "unsigned", power=softPower)
    TOM = TOMsimilarity(adjacency);
    dissTOM = 1-TOM
    geneTree = flashClust(as.dist(dissTOM), method = "average"); # Call the hierarchical clustering function
    dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize); # Module identification using dynamic tree cut
    dynamicColors = labels2colors(dynamicMods) # Convert numeric lables into colors
    #merge = mergeCloseModules(distgene, dynamicColors, cutHeight = MEDissThres, verbose = 3) # Call an automatic merging function
    #moduleColors = merge$colors; # The merged module colors
    moduleColors = dynamicColors
    colorOrder = c("grey", standardColors());
    Labels = match(moduleColors, colorOrder)-1;
    module <- cbind(rownames(distgene),Labels)
    module
}

build_CRF <- function(nodeinfo,allnet,module,netflag,strn,cutn=20000){
    library(CRF)
    library(Corbi)
    mapT <- as.matrix(read.delim("Fmap.txt",header=FALSE,sep="\t"))
    modLab <- setdiff(unique(module[,2]),0)
    n.module <- length(modLab)
    n.genes <- length(module[module[,2]!=0,1])
    
    POSpro <- matrix(0,n.genes,3)
    POSpro[,1] <- as.vector(module[module[,2]!=0,1]) # gene, probability; FDR, node score, edge score
    for (i in 1:n.module){
        modgenes <- as.vector(module[module[,2]==modLab[i],1])
        if(length(netflag)==1){
            net <- discrete_net(allnet,modgenes,netflag,cutn)
        }else{
            net <- combine_net(allnet,modgenes,netflag)
        }
        print(sum(net$matrix)/2)
        
        modnodeinfo <- nodeinfo[match(modgenes,rownames(nodeinfo)),]
        model <- build_model(modnodeinfo,net,bflag=1)
        crfresult <- infer_crf(model, query.type=4)
        POSpro[match(modgenes,POSpro[,1]),2:3] <- crfresult
    }
    
    POSpro <- POSpro[order(as.numeric(POSpro[,2]),decreasing = TRUE),]
    POSpro <- as.matrix(cbind(POSpro,as.numeric(POSpro[,3]) * dim(POSpro)[1]))
    POSpro[as.numeric(POSpro[,4])>1,4] <- 1
    
    pros1 <- POSpro
    pros1[,1] <- mapT[match(pros1[,1],mapT[,1]),2]
    pros1 <- pros1[!is.na(pros1[,1]),]
    write.table(POSpro,file=paste(strn,"LBP_",netflag,".txt",sep=""),quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")
    write.table(pros1,file=paste(strn,"LBP_",netflag,"1.txt",sep=""),quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")
    
}

rank_combine <- function(filenames,strn,flag=1,n){
    source("CRF_build.R")
    source("Network_analysis.R")
    
    mapT <- as.matrix(read.delim("Fmap.txt",header=FALSE,sep="\t"))
    n.net <- length(filenames)
    genes <- c()
    for(i in 1:n.net){
        tmp <- read.table(filenames[i])
        genes <- union(genes,tmp[,1])
    }
    n.genes <- length(genes)
    
    scoreM <- matrix(0,n.genes,n.net+1,dimnames=list(genes,1:(n.net+1)))
    pM <- matrix(1,n.genes,n.net+2,dimnames=list(genes,1:(n.net+2)))
    for(i in 1:n.net){
        tmp <- read.table(filenames[i])
        scoreM[match(tmp[,1],genes),i] <- as.numeric(tmp[,2])
        scoreM[is.na(scoreM[,i]),i] <- 0
        
        pM[match(tmp[,1],genes),i] <- as.numeric(tmp[,3])
        pM[pM[,i]==0,i] <- 1e-50
        pM[is.na(pM[,i]),i] <- 1
    }
    
    if(flag==1){
        y <- sapply(1:n.genes,function(i) { 
            if(any(pM[i,1:n.net]< 1)){
                exp(-(prod(-log(pM[i,pM[i,]< 1]))^(1/n.net)))
            }else{
                1;
            }
        })
    }else if(flag==2){
        y <- sapply(1:n.genes,function(i) { 
            if(any(pM[i,1:n.net]< 0.5)){
                exp(-(prod(-log(pM[i,pM[i,]< 0.5]))^(1/n.net)))
            }else{
                exp(-(prod(-log(pM[i,1:n.net]))^(1/n.net)))
            }
        })
    }else if(flag==3){
        y <- sapply(1:n.genes,function(i) { 
            if(any(pM[i,1:n.net]< 0.5)){
                exp(-(prod(-log(pM[i,pM[i,]< 0.5]))^(1/sum(pM[i,1:n.net]< 0.5))))
            }else{
                exp(-(prod(-log(pM[i,1:n.net]))^(1/n.net)))
            }
        })    
    }else if(flag==4){
        cut <- mean(sapply(1:n.net,function(i) quantile(scoreM[scoreM[,i]>0,i],0.99)))
        x <- sapply(1:n.genes,function(i) ifelse(sum(scoreM[i, scoreM[i,] > cut])==0,max(scoreM[i,]),sum(scoreM[i, scoreM[i,] > cut])))
        y1 <- plogis(x, location = 0.5, scale = 1, lower.tail = TRUE, log.p = FALSE)
        y <- 1-y1
    }else if(flag==5){
        cut <- mean(sapply(1:n.net,function(i) quantile(scoreM[scoreM[,i]>0,i],0.99)))
        x <- sapply(1:n.genes,function(i) ifelse(sum(scoreM[i, scoreM[i,] > cut])==0,max(scoreM[i,]),sum(scoreM[i, scoreM[i,] > cut])))
        y1 <- plogis(x, location = 0.5, scale = 1, lower.tail = TRUE, log.p = FALSE)
        
        y <- sapply(1:n.genes,function(i) { 
            if(sum(scoreM[i, scoreM[i,] > cut])==0){
                exp(-(prod(-log(pM[i,1:n.net]))^(1/n.net)))
            }else{
                pM[i,which(scoreM[i,]==max(scoreM[i,]))][1]
            }
        })
    }else if(flag==6){
        x <- rowSums(scoreM)
        y1 <- plogis(x, location = 0.5, scale = 1, lower.tail = TRUE, log.p = FALSE)
        y <- 1-y1
    }else if(flag==7){
        y1 <- sapply(1:n.genes,function(i) mean(scoreM[i,scoreM[i,]>0]))
        y <- 1-y1    
    }
    
    pM[,n.net+1] <- y
    pM[,n.net+2] <- y ### adjust????
    pM[pM[,n.net+2]>1,n.net+2] <- 1
    
    if(flag<=3){
        scoreM[,n.net+1] <- 1-pM[,n.net+1]
    }else if(flag>=4){
        scoreM[,n.net+1] <- y1
    }
    result <- as.matrix(cbind(scoreM,pM))
    result <- result[order(as.numeric(pM[,n.net+1]),decreasing=FALSE),]
    result0 <- result
    
    rownames(result) <- mapT[match(rownames(result),mapT[,1]),2]
    result <- result[!is.na(rownames(scoreM)),]

    filename <- paste(strn,flag,"resultall.txt",sep="_")
    filename0 <- paste(strn,flag,"resultall0.txt",sep="_")
    testname <- paste(strn,flag,"enrichall","1.txt",sep="_")    
    filetopgene <- paste(strn,flag,"topgeneall.txt",sep="_")
    write.table(result0,file=filename0,quote=FALSE,col.names=FALSE,sep="\t")
    write.table(result,file=filename,quote=FALSE,col.names=FALSE,sep="\t")
    
    Output_file1(filename,mapT,testname,netflag,n)
    testgenes <- rownames(result)[result[,n.net+2] < 0.05]
    write.table(testgenes,file=filetopgene,quote=FALSE,col.names=FALSE,row.names=FALSE)
}

rank_combinetmp <- function(filename,strn,netflag){
    source("CRF_build.R")
    source("Network_analysis.R")
    
    mapT <- as.matrix(read.delim("Fmap.txt",header=FALSE,sep="\t"))
    n.net <- length(filename)
    scoreM <- matrix(1,n.genes,n.net+2,dimnames=list(genes,1:(n.net+2)))
    
    tmp <- read.table(filename)
    genes <- unique(tmp[,1])
    n.genes <- length(genes)
    scoreM[match(tmp[,1],genes),1] <- as.numeric(tmp[,2])
    scoreM[match(tmp[,1],genes),2] <- as.numeric(tmp[,3])
    scoreM[is.na(scoreM[,1]),1] <- 0
    scoreM[is.na(scoreM[,2]),2] <- 1
    scoreM[scoreM[,2]==0,2] <- min(as.numeric(tmp[,3])[as.numeric(tmp[,3]) >0 ])/100
    scoreM[,n.net+2] <- scoreM[,n.net+1]*dim(scoreM)[1]    
    scoreM[scoreM[,n.net+2]>1,n.net+2] <- 1
    scoreM <- scoreM[order(as.numeric(scoreM[,n.net+1]),decreasing=TRUE),]

    result <- scoreM
    rownames(result) <- mapT[match(rownames(result),mapT[,1]),2]
    result <- result[!is.na(rownames(scoreM)),]
    
    testname <- paste(strn,"enrich",netflag,"1.txt",sep="_")
    filename <- paste(strn,"LBP_",netflag,"1.txt",sep="")
    filetopgene <- paste(strn,"topgene",netflag,"1.txt",sep="_")
    #write.table(scoreM,file=paste("result",strn,netflag,"10.txt",sep="_"),quote=FALSE,col.names=FALSE,sep="\t")
    #write.table(result[,(n.net+1):(n.net+2)],file=filename,quote=FALSE,col.names=FALSE,sep="\t")
    
    Output_file1(filename,mapT,testname,netflag)
    testgenes <- rownames(result)[result[,n.net+2] < 0.05]
    write.table(testgenes,file=filetopgene,quote=FALSE,col.names=FALSE,row.names=FALSE)
}

rank_combine1 <- function(filename,strn,netflag,n){
  
    mapT <- as.matrix(read.delim("Fmap.txt",header=FALSE,sep="\t"))
    testname <- paste(strn,"enrich",netflag,"1.txt",sep="_")    
    Output_file1(filename,mapT,testname,netflag,n)

}

#====================
# auxiliary function
#====================
# tmp used from CRF_build.R

combine_net <- function(allnet,modgenes,netflag){


}

discrete_net <- function(allnet,modgenes,netflag,cutn=20000){
    net <- list()
    net$node <- modgenes
    net$size <- length(modgenes)
    net$matrix <- allnet$matrix[modgenes,modgenes]
    
    if(netflag==1){
        cutf <- quantile(net$matrix[upper.tri(net$matrix)],0.95)
        net$matrix[net$matrix >= cutf] <- 1
        net$matrix[net$matrix < cutf] <- 0
    }
    
    if(netflag==2){
        neta <- minet::aracne(net$matrix)
        net$matrix <- neta[modgenes,modgenes]
        net$matrix[net$matrix>0] <- 1
    } 
    
    if(netflag==3 | netflag==4 | netflag==0){
        net$matrix[net$matrix > 0] <- 1
    }
    
    if(sum(net$matrix==1) > 2*cutn){ # memory limit
        net$matrix[lower.tri(net$matrix,diag=TRUE)] <- 0
        edges <- which(net$matrix==1,arr.ind=TRUE)
        tmp <- allnet$matrix[modgenes,modgenes]
        edges <- as.matrix(cbind(edges,tmp[edges]))
        edges <- edges[order(as.numeric(edges[,3]),decreasing=TRUE),]
        edges <- edges[1:cutn,1:2]
        net$matrix[upper.tri(net$matrix,diag=TRUE)] <- 0
        net$matrix[edges] <- 1
        net$matrix[edges[,c(2,1)]] <- 1        
    }
    
    diag(net$matrix) <- 0
    
    net
}

combine_allnet <- function(){
    source("Multi_net.R")
    library(WGCNA)
    library(minet)
    fileexp <- c("PCGCall.txt","Brain_expression_Ref.txt")
    k=1
    # net: a list: size node matrix (systemtic and zero diag)
    netlist <- list()
    nodes <- c()
   
    netlist[[1]] <- build_net(1,fileexp[k])
    cutf <- quantile(netlist[[1]]$matrix[upper.tri(netlist[[1]]$matrix)],0.999)
    netlist[[1]]$matrix[netlist[[1]]$matrix >= cutf] <- 1
    netlist[[1]]$matrix[netlist[[1]]$matrix < cutf] <- 0
    print("Network1")
    
    netlist[[2]] <- build_net(2,fileexp[k])
    neta <- minet::aracne(netlist[[2]]$matrix)
    netlist[[2]]$matrix <- neta
    netlist[[2]]$matrix[netlist[[2]]$matrix > 0] <- 1
    print("Network2")
    
    netlist[[3]] <- build_net(3,fileexp[k])
    netlist[[3]]$matrix[netlist[[3]]$matrix > 0] <- 1
    netlist[[4]] <- build_net(4,fileexp[k])
    netlist[[4]]$matrix[netlist[[4]]$matrix > 0] <- 1
    print("Network4")
    
    nodes <- union(nodes,netlist[[1]]$node)
    nodes <- union(nodes,netlist[[2]]$node)
    nodes <- union(nodes,netlist[[3]]$node)
    nodes <- union(nodes,netlist[[4]]$node)
      

    net0 <- list()
    net0$node <- nodes
    net0$size <- length(nodes)
    net0$matrix <- matrix(0,net0$size,net0$size,dimnames=list(nodes,nodes))
    for(i in 1:4){
        net0$matrix[netlist[[i]]$node,netlist[[i]]$node] <- net0$matrix[netlist[[i]]$node,netlist[[i]]$node] + netlist[[i]]$matrix
    }
    rm(netlist)
    print("RM") 
    net0$matrix[net0$matrix < 2] <- 0
    subs <- rowSums(net0$matrix) > 0
    net0$node <- net0$node[subs]
    net0$size <- sum(subs)
    net0$matrix <- net0$matrix[subs,subs]
    net0$matrix <- net0$matrix/max(net0$matrix)
    print("Saving")
    save(net0,file="net0");
    
    
    Labels <- build_module(net0$matrix)
    module <- cbind(net0$node,Labels)
    print("Module done")
    save(module,file="module0");
    
}

allrun <- function(){
    source("CRF_build.R")
    source("MixNet.R")
    source("Network_analysis.R")
    source("Multi_net.R")
    options(stringsAsFactors=FALSE)
    
    strn0 <- c("result/PCGClog","result/TADAnew","result/NDD","result/nonNDD","result/other","result/CHDnature","result/CHDother","result/P2","result/T2","result/P3","result/T3") 
    fileexp <- c("PCGCall.txt","PCGCall.txt","Brain_expression_Ref.txt","Brain_expression_Ref.txt","Brain_expression_Ref.txt","PCGCall.txt","PCGCall.txt","PCGCall.txt","PCGCall.txt","PCGCall.txt","PCGCall.txt")
    fnode <- c("NewlogPheno_score_map.txt","TADAinfo0113.txt","random_samples/NDDPheno_score_map.txt","random_samples/nonNDDPheno_score_map.txt","random_samples/otherPheno_score_map.txt","random_samples/CHDnaturePheno_score_map.txt","random_samples/CHDotherPheno_score_map.txt","random_samples/P2Pheno_score_map.txt","random_samples/T2Pheno_score_map.txt","random_samples/P3Pheno_score_map.txt","random_samples/T3Pheno_score_map.txt")
    for(k in 8:length(fnode)){
        filenames <- c()
        ngene=300
        for(i in 1:4){
            strn <- paste(strn0[k],i,sep="_")
            Multi_net(i,strn,fnode[k],fileexp[k]);
            rank_combine1(paste(strn,"LBP_",i,"1.txt",sep=""),strn,i,ngene)
            filenames <- union(filenames,paste(strn,"LBP_",i,".txt",sep=""))
        }
        flag=6
        rank_combine(filenames,strn0[k],flag,ngene)
    }
    
    
    strn0 <- "PAH/PAH"
    fileexp <- "PCGCall.txt"
    fnode <- "PAH/PAHPheno_score_map.txt"

        filenames <- c()
        ngene=300
        for(i in 3:4){
            strn <- paste(strn0,i,sep="_")
            Multi_net(i,strn,fnode,fileexp);
            rank_combine1(paste(strn,"LBP_",i,"1.txt",sep=""),strn,i,ngene)
            filenames <- union(filenames,paste(strn,"LBP_",i,".txt",sep=""))
        }
        flag=6
        rank_combine(filenames,strn0,flag,ngene)
    
}
