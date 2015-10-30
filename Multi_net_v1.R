Multi_net <- function(netflag,strn="LBP"){
    #netflag <- 1 # 1: expression; 2: regulatory ARACNE; 3: PPI; 4: Phenotype based
    #filenames <- c("TADA_lofmis1202.csv","case.de.novo.cnv.genes.bed","Othergenes.txt")
    #fileexp <- "PCGCall.txt"
    
    filenames <- c("NDD_TADA1205.csv","case.de.novo.cnv.genes.bed","Othergenes.txt")
    fileexp <- "Brain_expression_Ref.txt"
    allnet <- build_net(netflag,fileexp)
    module <- build_module(allnet$matrix)
    nodeinfo <- build_node(module,2,filenames)
    build_CRF(nodeinfo,allnet,module,netflag,strn)
}

batch_Multi_net <- function(){
    source("CRF_build.R")
    source("MixNet.R")
    source("Network_analysis.R")
    source("Multi_net.R")
    options(stringsAsFactors=FALSE)
    
    strn0 <- "NDD" #strn: the label for different data set names
    filenames <- c()
    for(i in 1:4){
        strn <- paste(strn0,i,sep="_")
        Multi_net(i,strn);
        filename <- paste("LBP_",strn,"_",i,".txt",sep="")
        rank_combine(filename,strn,i,0)
        filenames <- union(filenames,filename)
    }
    
    rank_combine(filenames,strn0,1,1)
    
}

build_node <- function(module,flag=1,filenames=c("TADA_lofmis1202.csv","case.de.novo.cnv.genes.bed","Othergenes.txt")){
    #===============
    # TADA information
    mapT <- as.matrix(read.delim("Fmap.txt",header=FALSE,sep="\t")) ###!!! 可以进一步精细
    genes <- module[,1]
    mapT <- mapT[mapT[,1] %in% genes,]
    
    TADAFile <- filenames[1]
    #TADAFile <- "TADA_lofmis1202.csv"
    geneinfo <- as.matrix(read.csv(TADAFile))
    genelist <- geneinfo[,1]
    genelist <- mapping_to(genelist)
    genelist <- mapT[match(genelist,mapT[,2]),1]
    nodesim <- as.matrix(as.numeric(geneinfo[,"BF"])/(as.numeric(geneinfo[,"BF"]) + 1))
    names(nodesim) <- genelist
    subs <- !is.na(genelist)
    nodesim <- nodesim[subs]
    
    if(flag==1){
        nodeinfo <- matrix(,length(genes),3,dimnames=list(genes,1:3))
        nodeinfo[,1] <- 0.1 ## related genes percent estimate
        nodeinfo[,2:3] <- 0
        nodeinfo[names(nodesim),1] <- nodesim
        
        # CNV information; filtered: some CNVs affect a very large number of genes. double counting closely related genes like "UBXN10" and "UBXN10-AS1" 
        CNVfile <- filenames[2]
        #CNVfile <- "case.de.novo.cnv.genes.bed"
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
        Otherfile <- filenames[3]
        #Otherfile <- "Othergenes.txt"
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

build_CRF <- function(nodeinfo,allnet,module,netflag,strn){
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
            net <- discrete_net(allnet,modgenes,netflag)
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
    pros1 <- mapT[match(POSpro[,1],mapT[,1]),2]
    pros1 <- cbind(pros1,POSpro[,2])
    write.table(POSpro,file=paste("LBP_",strn,"_",netflag,".txt",sep=""),quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")
    write.table(pros1,file=paste("LBP",strn,netflag,"1.txt",sep="_"),quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")
    
    ### tmp used for debug
#     testname <- paste("enrich",strn,netflag,"1.txt",sep="_")
#     filename <- paste("LBP",strn,netflag,"1.txt",sep="_")
#     Output_file1(filename,mapT,testname,netflag)
#     print("The known genes ranks:")
#     print(match(c("RBFOX2","KMT2D","EP300","NAA16","KDM5B","NCOA3","SMAD2","PTPN11","KDM5B","NAA15","POGZ","AHNAK","MYH6","NOTCH1","SUV420H1","MASTL","CHD7","IRG1","BRAF","SMARCC1","LRIG1","CDKL2","JAG1","FLT4","PCNXL2","MKRN2","MROH7","TLK2","ARID1B","CTNNB1","ETS1","CUL3"),pros1[,1]))
}

rank_combine <- function(filenames,strn,netflag,flag=0){
    source("CRF_build.R")
    source("Network_analysis.R")
    
    mapT <- as.matrix(read.delim("Fmap.txt",header=FALSE,sep="\t"))
    n.net <- length(filenames)
    genes <- c()
    for(i in 1:n.net){
        #filename <- paste("LBP_",strn,"_",i,"_",i,".txt",sep="")
        tmp <- read.table(filenames[i])
        genes <- union(genes,tmp[,1])
    }
    n.genes <- length(genes)
    scoreM <- matrix(0,n.genes,n.net+2,dimnames=list(genes,1:(n.net+2)))
    for(i in 1:n.net){
        #filename <- paste("LBP_",strn,"_",i,"_",i,".txt",sep="")
        tmp <- read.table(filenames[i])
        scoreM[match(tmp[,1],genes),i] <- as.numeric(tmp[,2])
    }
    scoreM[is.na(scoreM)] <- 0
    x <- rowSums(scoreM)
    if(flag==0){
        y <- x
    }else{
        y <- plogis(x, location = 0.5, scale = 1/n.net, lower.tail = TRUE, log.p = FALSE)
    }
    scoreM[,n.net+1] <- y
    scoreM <- scoreM[order(as.numeric(scoreM[,n.net+1]),decreasing=TRUE),]
    #scoreM <- cbind(scoreM,(1-scoreM[,n.net+1])*dim(scoreM)[1])
    scoreM[,n.net+2] <- (1-scoreM[,n.net+1])*dim(scoreM)[1]
    scoreM[scoreM[,n.net+2]>1,n.net+2] <- 1
    
    result <- scoreM
    rownames(result) <- mapT[match(rownames(result),mapT[,1]),2]
    result <- result[!is.na(rownames(scoreM)),]
    
    testgenes <- rownames(result)[result[,n.net+2] < 0.05]
    
    if(flag==0){
        testname <- paste("enrich",strn,netflag,"1.txt",sep="_")
        filename <- paste("result",strn,netflag,"1.txt",sep="_")
        filetopgene <- paste("topgene",strn,netflag,"1.txt",sep="_")
        #write.table(scoreM,file=paste("result",strn,netflag,"10.txt",sep="_"),quote=FALSE,col.names=FALSE,sep="\t")
        write.table(result[,(n.net+1):(n.net+2)],file=filename,quote=FALSE,col.names=FALSE,sep="\t")
    }else if(flag==1){
        filename <- paste(strn,"LBPresult.txt",sep="")
        filename0 <- paste(strn,"LBPresult0.txt",sep="")
        testname <- paste(strn,"enrichall","1.txt",sep="_")    
        filetopgene <- paste(strn,"topgenes.txt",sep="")
        write.table(scoreM,file=filename0,quote=FALSE,col.names=FALSE,sep="\t")
        write.table(result[,(n.net+1):(n.net+2)],file=filename,quote=FALSE,col.names=FALSE,sep="\t")
    }
    
    Output_file1(filename,mapT,testname,netflag)
    write.table(testgenes,file=filetopgene,quote=FALSE,col.names=FALSE,row.names=FALSE)
}
#====================
# auxiliary function
#====================
# tmp used from CRF_build.R

combine_net <- function(allnet,modgenes,netflag){


}

discrete_net <- function(allnet,modgenes,netflag){
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
    
    if(netflag==3 | netflag==4){
        net$matrix[net$matrix > 0] <- 1
    }
    
    if(sum(net$matrix==1) > 40000){ # memory limit
        net$matrix[lower.tri(net$matrix,diag=TRUE)] <- 0
        edges <- which(net$matrix==1,arr.ind=TRUE)
        tmp <- allnet$matrix[modgenes,modgenes]
        edges <- as.matrix(cbind(edges,tmp[edges]))
        edges <- edges[order(as.numeric(edges[,3]),decreasing=TRUE),]
        edges <- edges[1:20000,1:2]
        net$matrix[upper.tri(net$matrix,diag=TRUE)] <- 0
        net$matrix[edges] <- 1
        net$matrix[edges[,c(2,1)]] <- 1        
    }
    
    diag(net$matrix) <- 0
    
    net
}
