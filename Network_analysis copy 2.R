# Network analysis 
Network_analysis_batch <- function(){
  source("Network_analysis.R")
   for(i in 1){
     for(j in 1:2){
       Network_analysis(i, j);
     }
   }
  
}

Network_analysis <- function(fileflag = 1, distflag = 1,cflag=0){

# initial packages and files
# -----
  library(RUVSeq)
  library(stringr)
  library(WGCNA)
  #library(parmigene)
  #library(entropy)
  library(infotheo)  
  mapT <- as.matrix(read.delim("Fmap.txt",header=FALSE,sep="\t"))
# ----

# step 1: network construction
# ---- 
    if(cflag ==0){
  # 0) PCA analysis and normalized methods
  dataFile <- c("GTEx_heart_expression.txt","GTEx.heart.countstable.txt","PCGC_newBorn_FPKM.txt","fetal_JUNE_AUG_expression.txt")
  batchFile <- "renameFileWithHeartTissue1125.list"
  # if use, housekeeping genes File "HK_genes.txt"
  datExpr <- step_1_0_normalize(dataFile[fileflag],batchFile,fileflag)
    if(fileflag<=2){
        check_batch(datExpr,batchFile)
    }
  # 1) co-expression based on pearson coefficient 
  # 2) based on mutual information
  if(distflag==1){
    distfile <- paste("distgene",fileflag,distflag,sep="_")
    load(distfile)
    distgene[is.na(distgene)] <- 0
    distgene <- abs(distgene) # this for DAWN
    diag(distgene) <- 0 # to delete isolated nodes
  }else if(distflag==2){
    distfile <- paste("dist",fileflag,distflag,sep="_")
    load(distfile)  
    distgene <- distM
    rm(distM)
    distgene[is.na(distgene)] <- 0
    diag(distgene) <- 0 
    distgene <- distgene/max(distgene)
  }

  genes <- colnames(datExpr)
  dimnames(distgene) <- list(genes,genes)
    }else if(cflag==1){
        tmp <- combine_batch_samples(distflag)
        datExpr <- tmp$datExpr
        distgene <- tmp$distgene
        rm(tmp)
        flag <- 3 #!!!!
    }else if(cflag==2){
        tmp <- combine_batch_networks(distflag)
        datExpr <- tmp$datExpr
        distgene <- tmp$distgene
        rm(tmp)   
        flag <- 4 #!!!!
    }


    # delete isolated node
    x <- colSums(distgene>0)
    distgene <- distgene[x>0,x>0]
    genes <- colnames(distgene)
    datExpr <- datExpr[,genes]
# ---- 

# step 2: network check
# ----
  knownFile <- "known_genes.csv"
  # plot network degree distribution
  # power_plot(distgene)
  # known gene connected
  if(fileflag==2) {knowngenes <- as.matrix(read.csv(knownFile,header=FALSE))[,1];} else {knowngenes <- as.matrix(read.csv(knownFile,header=FALSE))[,2];}
  # to check
  print("These are known genes to check:")
  #as.matrix(read.csv(knownFile,header=FALSE))
  knowngenes <- intersect(knowngenes,genes)
  n.known <- length(knowngenes)
#----
  #colSums(distgene[,knowngenes]>0)
#   tmp <- matrix(0,n.known,n.known)
#   meantmp <- matrix(0,n.known,n.known)
#   for(i in 1:n.known){
#     for(j in 1:n.known){
#       tmp[i,j] <- t.test(as.vector(distgene[,knowngenes[i]]),as.vector(distgene[,knowngenes[j]]))$p.value
#       meantmp[i,j] <- mean(as.vector(distgene[,knowngenes[i]]))/mean(as.vector(distgene[,knowngenes[j]]))
#     }
#     }
#   tmp[c(1,3,2,8),c(1,3,2,8)]
#   meantmp[c(1,3,2,8),c(1,3,2,8)]
#----
    density_plot(knowngenes[c(1,3,2,8)],distgene)


  # known TF degree
    TFFile <- "TF-co-regulate.csv"
  knownTF <- read.csv(TFFile,header=FALSE)
  knownTF <- union(knownTF[,1],knownTF[,2])
  knownTF <- mapping_to(knownTF)
  if(fileflag!=2)  knownTF <- mapT[match(knownTF,mapT[,2]),1]; 
  knownTF <- unique(knownTF[!is.na(knownTF)])
  #deg <- colSums(distgene>0)
  #x <- deg[genes %in% knownTF]
  #y <- deg[!(genes %in% knownTF)]
  #t.test(x,y)
  x <- genes %in% knownTF
  y <- !(genes %in% knownTF)
  hist_plot(as.vector(distgene[,x]),as.vector(distgene[,y]))
    

# ----

# step 3: DAWN run
# ----
  # !!! random slect half of mutation to run and whole to test significant
  TADAFile <- "TADA_lofmis1121.csv"
  pvalue_data <- read.csv(TADAFile)[,c(1,10)]
  pvalue_data[,1] <- mapping_to(pvalue_data[,1])
  if(fileflag!=2)  pvalue_data[,1] <- mapT[match(pvalue_data[,1],mapT[,2]),1]; 
  pvalue_data <- pvalue_data[!(is.na(pvalue_data[,1]) | is.na(pvalue_data[,2])),]
  pvalue_data <- pvalue_data[pvalue_data[,2]< 1,]
  pvalue_data <- pvalue_data[!duplicated(pvalue_data[,1]),]
  
  moduleLabels <- WGCNA_module(distgene,datExpr,distflag,flag)
  # delete the module larger than 5000
  moduleLabels <- delete_module(moduleLabels,n=5000,pvalue_data[,1],genes)
  module <- data.frame(GENE=genes,MODULE=moduleLabels)
  diag(distgene) <- 1
  source("DAWN_sourcefile.R")
  datExpr <- t(datExpr)
  filename <- paste("DAWNresult_",fileflag,"_",distflag,"_",cflag,".csv",sep="")
    if(flag == 1){
    DAWN_package(expression_data=datExpr,module_data=module, pvalue_data=pvalue_data,cor_thres1=0.6,cor_thres2=0.5,nASDthres=0.5,fdrthres=0.1,flag_q=distflag,filename=filename,distgene=distgene)
    }else if(flag==2){
     DAWN_package(expression_data=datExpr,module_data=module, pvalue_data=pvalue_data,cor_thres1=0.19,cor_thres2=0.13,nASDthres=0.5,fdrthres=0.1,flag_q=distflag,filename=filename,distgene=distgene)
    }else if(flag ==3){
    
    }else if(flag==4){
        DAWN_package(expression_data=datExpr,module_data=module, pvalue_data=pvalue_data,cor_thres1=0.6,cor_thres2=0.5,nASDthres=0.5,fdrthres=0.1,flag_q=distflag,filename=filename,distgene=distgene)
        filename <- filename <- paste("DAWNresult_",fileflag,"_",distflag,"_",cflag,"_0.75",".csv",sep="")
        DAWN_package(expression_data=datExpr,module_data=module, pvalue_data=pvalue_data,cor_thres1=0.75,cor_thres2=0.7,nASDthres=0.5,fdrthres=0.1,flag_q=distflag,filename=filename,distgene=distgene)
    }
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

# step 4: result analysis
# ----
  # histone modification genes  
  # Networks and pathways (transcriptional regulation, TGF-beta, cell adhesion)
  # GATA4, TBX5, Nkx2-5 etc
  # alternative splicing gene: RBFOX2
  HMGFile <- "HMGs.csv"
  TGFFile <- "TGF-beta.txt"
  #filename <-  "DAWN_analysis_result.csv" # for test  
  HMGgene <- read.csv(HMGFile,header=FALSE)[,2] # there is an unkonwn in 152 HMGgene
  TGFgene <- read.delim(TGFFile,header=FALSE,sep="\t")[,2]
  HMGgene <- mapping_to(HMGgene)
  TGFgene <- mapping_to(TGFgene)

  result <- as.matrix(read.csv(filename))
  if(fileflag!=2){  result[,2] <- mapT[match(result[,2],mapT[,1]),2];
  write.csv(result,file=paste(unlist(strsplit(filename,"\\."))[1],"_1.csv",sep=""));}
  rASDgene <- result[result[,4]==1,2]
  
 
  # check details
  # 1) whether includes RBFOX2
  print("RBFOX2 is: ")
  "RBFOX2" %in% toupper(rASDgene)

  # 2) enrich or not: HMG genes and TGF genes
  # n.genes <- dim(result)[1]
  print("HMG genes test p value is:")
  hyper_test(HMGgene,rASDgene,result[,2])
  print("TGF-beta gene test p value is:")
  hyper_test(TGFgene,rASDgene,result[,2])

  # 3) known genes in or not
  knownFile <- "known_genes.csv"
  knowngenes <- as.matrix(read.csv(knownFile,header=FALSE))[,1]
  knowngenes <- mapping_to(knowngenes)
  print("Total known genes and in rASD genes: ")  
  length(knowngenes)
  which(toupper(knowngenes) %in% toupper(rASDgene))
#----

}

# ----
step_1_0_normalize <- function(dataFile,batchFile,flag=1){
    
    # delete no expression genes
    datExpr <- read.delim(dataFile,sep="\t")
    genes <- datExpr[,1]
    samples <- colnames(datExpr)[-1]
    datExpr <- as.matrix(datExpr[,-1])
    datExpr <- t(datExpr)
    colnames(datExpr) <- genes  
    datExpr <- datExpr[,colSums(abs(datExpr))>0]
    
    if(flag == 2){
      # get the batch information
        samplenames <- as.matrix(read.delim(batchFile,header=FALSE,sep="\t"))
        samplenames <- cbind(samplenames,paste(samplenames[,1],1,sep=""))
        samplenames <- cbind(samplenames,paste(samplenames[,1],1,sep="_"))
      # gene name mapping
      genes <- colnames(datExpr)
      genes <- mapping_to(genes,"fi",mapfile="GENE_SYNONYMS")
      colnames(datExpr) <- genes
      # count table RUVSeq normalize
      datExpr <- RUV_normalized(datExpr,samplenames)
    }
    # PCA plot
    sampled <- dist(datExpr)
    loc <- cmdscale(sampled,k=2)
    plot(loc[,1],loc[,2],main="MDS plot",xlab="PC1",ylab="PC2")
    loc <- prcomp(datExpr)
    plot(loc$x[,1],loc$x[,2],main="PCA plot",xlab="PC1",ylab="PC2")
    
    # clustering and outlier delete
    # WCGNA delete the outlier data
    gsg <- goodSamplesGenes(datExpr,verbose=3)    # print(gsg$allOK)
    if (!gsg$allOK){ datExpr = datExpr[gsg$goodSamples, gsg$goodGenes];}
    
    sampleTree = flashClust(sampled,method="average")
    
    if(flag==1) {height <- 35000;} else if(flag==2) {height <- 1500000;} else if(flag==3) {height <- 85000;}else if(flag==4){height <- 120000;}
    plot(sampleTree,main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
    abline(h=height,col="red")
    
    clust = cutreeStatic(sampleTree, cutHeight = height, minSize = 10) ### choose the cutHeight
    
    # split into two group expression values
    keepSamples <-  clust == 1
    datExpr <- datExpr[keepSamples,]
    datExpr <- datExpr[,colSums(abs(datExpr))>0]
    
    datExpr
}

RUV_normalized <- function(datExpr,samplenames){
  
  datExpr1 <- datExpr[samplenames[samplenames[,4]=="Heart - Left Ventricle",5],]
  subs <- colSums(abs(datExpr1))>0
  datExpr1 <- datExpr1[,subs]
  # dim(datExpr1)
  datExpr2 <- datExpr[samplenames[samplenames[,4]=="Heart - Atrial Appendage",5],]
  subs <- colSums(abs(datExpr2))>0
  datExpr2 <- datExpr2[,subs]  
  # dim(datExpr2)
  
  genes <- intersect(colnames(datExpr1),colnames(datExpr2))
  datExpr1 <- datExpr1[,genes]
  datExpr2 <- datExpr2[,genes]
  # normalized two group by RUV
  # house keeping genes
  HK_g <- as.matrix(read.table("HK_genes.txt",sep="\t",header=FALSE))
  HK_g[,1] <- str_trim(HK_g[,1])
  HK_g <- mapping_to(HK_g,"fi",mapfile="GENE_SYNONYMS")
  HK_g <- HK_g[HK_g %in% genes]
  
  seqRUVg1 <- RUVg(t(datExpr1), HK_g, k=3)
  datExprN1 <- t(seqRUVg1$normalizedCounts)
  
  seqRUVg2 <- RUVg(t(datExpr2), HK_g, k=3)
  datExprN2 <- t(seqRUVg2$normalizedCounts)
  
  datExprnew <- rbind(datExprN1,datExprN2)
  
  datExprnew
}

check_batch <- function(datExpr,batchFile){
    samplenames <- as.matrix(read.delim(batchFile,header=FALSE,sep="\t"))
    samplenames <- cbind(samplenames,paste(samplenames[,1],1,sep=""))
    samplenames <- cbind(samplenames,paste(samplenames[,1],1,sep="_"))
    Exprnames <- rownames(datExpr)
    n0 <- length(Exprnames)
    print(n0)
    if(sum(samplenames[,6]%in% Exprnames)==0){
        subbatch <- samplenames[,5]%in%Exprnames
    }else {subbatch <- samplenames[,6]%in% Exprnames;}
    tmp <- samplenames[subbatch,]
    table(tmp[,4])
    #tmp[tmp[,4]=="Heart - Atrial Appendage",2]
   
}

power_plot <- function(distgene){
  
  x <- colSums(distgene>0)
  hist(x,main="Histgram of degree",xlab="Degree",ylab="Frequency")
   
  library(igraph)
  library(poweRlaw)
  #Fit the power-laws
  m1 = displ$new(x)
  m1$setXmin(estimate_xmin(m1))
  
  #Plot and add lines as normal:
  plot(m1, xlab="Degree", ylab="Percentage", main="Degree distribution")
  lines(m1, col=2)

}

density_plot <- function(subgenes,distgene){
    # plot density for GATA4, NKX2-5, TBX1, TBX5
    n.genes <- length(subgenes)
    plot(density(distgene[,subgenes[1]]),main="Density for GATA4, NKX2-5, TBX1, TBX5",xlab="correlation",col=1)
    for(i in 2:n.genes){
        lines(density(distgene[,subgenes[i]]),col=i)
    }
    lines(density(distgene[upper.tri(distgene)]),col=n.genes+1)
    legend("topright",col=1:(n.genes+1),legend=c("GATA4","NKX2-5","TBX5","TBX1","ALL genes"),lty=1)
    
    
    # plot GATA4 and NKX2-5 correlation significant
    n <- dim(distgene)[1]
    C <- distgene[subgenes[1],subgenes[2]]
    P <- sum(distgene[upper.tri(distgene)] > C )*2/(n*(n-1))
    allden <- density(distgene[upper.tri(distgene)])
    plot(allden,col=1,main="Gene pairs correlation distribution",xlab="correlation")
    abline(v=C,col=2)   
    text(C,max(allden$y)/2, paste("GATA4-NKX2-5:",C,", p-vlaue is: ",P,sep=""), col = 2)

}

hist_plot <- function(X,Y){
    plot(density(X),main="Correlations distributions for Known TFs  and Other genes",xlab="correlation",col=1)
    lines(density(Y),col=2)
    legend("topright",col=1:2,legend=c("TFs","Other genes"),lty=1)
}

WGCNA_module <- function(distgene,datExpr,distflag,flag){
  
  # WGCNA to module relationships
  # when use the pearson correlation, it should be changed to power 6 for WGCNA
    
  if(flag ==1 ){
    softPower = 6;
    adjacency = adjacency(datExpr, power = softPower);  
    MEDissThres = 0.5 ## A cut off value, We choose a height cut of 0.25, corresponding to correlation of 0.75, to merge
  }else if(flag==2){
    softPower = 3;
    diag(distgene) <- 1
    adjacency = adjacency.fromSimilarity(distgene, type = "unsigned", power=softPower)
    MEDissThres = 0.01
  }else if(flag==3){
    softPower = 7;
    adjacency = adjacency(datExpr, power = softPower);
    MEDissThres = 0.1
  }else if(flag==4){
      softPower = 9;
      diag(distgene) <- 1
      adjacency = adjacency.fromSimilarity(distgene, type = "unsigned", power=softPower)
      MEDissThres = 0.05  
  }
  
  minModuleSize = 30; # set the minimum module size relatively high:
  
  
  TOM = TOMsimilarity(adjacency);
  dissTOM = 1-TOM
  # Call the hierarchical clustering function
  geneTree = flashClust(as.dist(dissTOM), method = "average");
  # Module identification using dynamic tree cut:
  dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize);

  # Convert numeric lables into colors
  dynamicColors = labels2colors(dynamicMods)
  
  # Should we change this for ????
  # Calculate eigengenes
  MEList = moduleEigengenes(datExpr, colors = dynamicColors)
  MEs = MEList$eigengenes
  # Calculate dissimilarity of module eigengenes
  if(distflag==1) { MEDiss = 1-cor(MEs);
  }else if(distflag==2){
    tmp <-discretize(MEs)
    muthq <- mutinformation(tmp)
    muthq <- muthq/max(muthq)
    diag(muthq) <- 1
    MEDiss <- 1- muthq
  }
  
  # Cluster module eigengenes
  METree = flashClust(as.dist(MEDiss), method = "average");
  
  # Call an automatic merging function
  merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
  # The merged module colors
  mergedColors = merge$colors;
  
  moduleColors = mergedColors
  colorOrder = c("grey", standardColors());
  moduleLabels = match(moduleColors, colorOrder)-1;
  
  moduleLabels
  
}

delete_module <- function(moduleLabels,n,pvalue_genes,genes){
  
  tmp <- table(moduleLabels)
  tmpna <- names(table(moduleLabels))
  for(i in 1:length(tmp)){
    if(tmp[i] > n & as.numeric(tmpna[i])!=0){
      moduleLabels[moduleLabels==as.numeric(tmpna[i])] <- 0
    }
    subgenes <- genes[moduleLabels==as.numeric(tmpna[i])]
    if(as.numeric(tmpna[i])!=0 & length(intersect(subgenes,pvalue_genes))==0){
      moduleLabels[moduleLabels==as.numeric(tmpna[i])] <- 0
    }
  }
  
  tmp <- table(moduleLabels)
  tmpna <- names(table(moduleLabels))
  subs <- as.numeric(tmpna) > 0
  tmpna <- tmpna[subs]
  tmp <- tmp[subs]
  for(i in 1:length(tmp)){
      moduleLabels[moduleLabels==as.numeric(tmpna[i])] <- i
  }
  
  moduleLabels
}

hyper_test <- function(geneset,rgeneset,genes){
  m <- length(intersect(genes,geneset))
  k <- length(rgeneset)
  n <- length(genes) - m
  x <- length(intersect(rgeneset,geneset))
  P <- phyper(x-1, m, n, k, lower.tail = FALSE, log.p = FALSE) 
  P
}

mapping_to <- function(Tab,filename,mapfile="GENE_SYNONYMS"){
  
  Tab <- as.matrix(Tab)
  map <- as.matrix(read.delim(mapfile,sep="\t",header=FALSE))
  map0 <- map
  map <- toupper(map)
  tmptab <- paste("|",map[,2],"|",sep="")
  tmpmap <- paste("|",map[,1],"|",sep="")
  #tmptab <- toupper(tmptab)
  
  map <- map0
  
  library(stringr)
  tmptabname <- paste("|",toupper(as.character(Tab[,1])),"|",sep="")
  for(i in 1:dim(Tab)[1]){
    subs <- which(str_detect(tmpmap,fixed(tmptabname[i])))
    #if(length(subs)==1){
    #  Tab[i,1] <- map[subs,1]	
    #}else 
    if(length(subs)==0){
      subs1 <- which(str_detect(tmptab,fixed(tmptabname[i])))
      if(length(subs1)==1){
        Tab[i,1] <- map[subs1,1]	
      }
    }
  }
  
  #write.table(Tab,file=filename,quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")
  
  Tab
}

make_batchfile <- function(){
    batchFile <- "renameFileWithHeartTissue.list"
    GtexFile <- "GTEx_Data_2014-01-17_Annotations_SampleAttributesDS.txt"
    batch <- as.matrix(read.delim(batchFile,header=FALSE,sep="."))
    gtex <- as.matrix(read.delim(GtexFile,sep="\t"))[,c("SAMPID","SMTS","SMTSD")]
    subs <- match(batch[,2],gtex[,"SAMPID"])
    batchnew <- batch[,1:2]
    batchnew <- cbind(batchnew,gtex[subs,2:3])
    batchnew[89,3:4] <- batchnew[1,3:4]
    write.table(batchnew,file="renameFileWithHeartTissue1125.list",quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")
}

check_count_FPKM <- function(){
    fileflag <- 1
    distflag <- 1
    
    dataFile <- c("GTEx_heart_expression.txt","GTEx.heart.countstable.txt","PCGC_newBorn_FPKM.txt","fetal_JUNE_AUG_expression.txt")
    batchFile <- "renameFileWithHeartTissue1125.list"
    # if use, housekeeping genes File "HK_genes.txt"
    datExpr <- step_1_0_normalize(dataFile[fileflag],batchFile,fileflag)
    datExpr1 <- datExpr
    
    fileflag <- 2
    distflag <- 1
    
    datExpr <- step_1_0_normalize(dataFile[fileflag],batchFile,fileflag)
    colnames(datExpr) <- mapT[match(colnames(datExpr),mapT[,2]),1]
    
    genes <- intersect(colnames(datExpr),colnames(datExpr1))
    sample <- rownames(datExpr)
    sample1 <- rownames(datExpr1)
    sample1 <- gsub("_","",sample1)
    rownames(datExpr1) <- sample1
    sam <- intersect(sample,sample1)
    
    datExpr1 <- datExpr1[sam,genes]
    datExpr <- datExpr[sam,genes]
    
    a <- sapply(1:17423,function(i) cor(datExpr[,i],datExpr1[,i]))
    
    hist(a)
}

check_WGCNA <- function(distgene,datExpr,distflag,MEDissThres = 0.15){
    
    if(distflag ==1 ){
        softPower = 6;
        adjacency = adjacency(datExpr, power = softPower);  
    }else if(distflag==2){
        softPower = 3;
        diag(distgene) <- 1
        adjacency = adjacency.fromSimilarity(distgene, type = "unsigned", power=softPower)
    }else if(distflag==3){
        softPower = 7;
        adjacency = adjacency(datExpr, power = softPower);
    }else if(distflag==4){
        softPower = 9;
        adjacency = adjacency(datExpr, power = softPower);
    }
    
    minModuleSize = 30; # set the minimum module size relatively high:
    
    TOM = TOMsimilarity(adjacency);
    dissTOM = 1-TOM
    # Call the hierarchical clustering function
    geneTree = flashClust(as.dist(dissTOM), method = "average");
    # Module identification using dynamic tree cut:
    dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize);
    # Convert numeric lables into colors
    dynamicColors = labels2colors(dynamicMods)
    # Calculate eigengenes
    MEList = moduleEigengenes(datExpr, colors = dynamicColors)
    MEs = MEList$eigengenes
    # Calculate dissimilarity of module eigengenes
    if(distflag==1) { MEDiss = 1-cor(MEs);
    }else if(distflag==2){
        tmp <-discretize(MEs)
        muthq <- mutinformation(tmp)
        muthq <- muthq/max(muthq)
        diag(muthq) <- 1
        MEDiss <- 1- muthq
    }
    # Cluster module eigengenes
    METree = flashClust(as.dist(MEDiss), method = "average");
    # Call an automatic merging function
    merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
    # The merged module colors
    mergedColors = merge$colors;
    moduleColors = mergedColors
    colorOrder = c("grey", standardColors());
    moduleLabels = match(moduleColors, colorOrder)-1;
       
    genes <- colnames(datExpr)
    HMGFile <- "HMGs.csv"
    TGFFile <- "TGF-beta.txt"
    HMGgene <- read.csv(HMGFile,header=FALSE)[,2] 
    TGFgene <- read.delim(TGFFile,header=FALSE,sep="\t")[,2]
    HMGgene <- mapping_to(HMGgene)
    TGFgene <- mapping_to(TGFgene)
    knownFile <- "known_genes.csv"
    knowngenes <- as.matrix(read.csv(knownFile,header=FALSE))[,2]
    #knowngenes <- as.matrix(read.csv(knownFile,header=FALSE))[,1]
    HMGgene <- mapT[match(HMGgene,mapT[,2]),1]
    TGFgene <- mapT[match(TGFgene,mapT[,2]),1]
    
    n.mo <- max(moduleLabels)
    a1 <- 1
    a2 <- 1
    a3 <- 1
    
    for(i in 1:n.mo){
        modulegene <- genes[moduleLabels==i]
        a1 <- min(a1,hyper_test(HMGgene,modulegene,genes))
        a2 <- min(a2,hyper_test(TGFgene,modulegene,genes))
        a3 <- min(a3,hyper_test(knowngenes,modulegene,genes))   
    }
    
    
    a <- c(a1,a2,a3)
    
    a
    
}

check_power <- function(distgene,datExpr,distflag){
    powers = c(c(1:10), seq(from = 12, to=20, by=2))
    if(distflag ==1) {
        sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
    }else if(distflag == 2){
        sft = pickSoftThreshold.fromSimilarity(distgene, powerVector = powers, verbose = 5)
    }
    
    par(mfrow=c(1,2))
    cex1=0.9
    plot(sft$fitIndices[,1],-sign(sft$fitIndices[,3]*sft$fitIndices[,2]), ylim=c(0,1), xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n", main = paste("Scale independence"));
    text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers,cex=cex1,col="red");
    abline(h=0.8,col="red")
    
    plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))
    text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
    
    par(mfrow=c(1,1))

}

combine_batch_samples <- function(distflag){

    datExpr <- combine_samples()
    if(distflag==1){
        distgene <- adjacency(datExpr, power = 1, type="unsigned")
        distgene[is.na(distgene)] <- 0
        diag(distgene) <- 0
    }else if(distflag==2){
        distfile <- paste("distc",distflag,sep="_")
        load(distfile)  
        distgene <- distM
        rm(distM)
        distgene[is.na(distgene)] <- 0
        diag(distgene) <- 0 
        distgene <- distgene/max(distgene)
    }
    
    genes <- colnames(datExpr)
    dimnames(distgene) <- list(genes,genes)
    
    list(datExpr=datExpr,distgene=distgene)
}

combine_batch_networks <- function(distflag){

    dataFile <- c("GTEx_heart_expression.txt","GTEx.heart.countstable.txt","PCGC_newBorn_FPKM.txt","fetal_JUNE_AUG_expression.txt")
    batchFile <- "renameFileWithHeartTissue1125.list"
    datExpr1 <- step_1_0_normalize(dataFile[1],batchFile,1)
    datExpr2 <- step_1_0_normalize(dataFile[3],batchFile,3)
    
    if(distflag==1){
    load("distgene_1_1")
    distgene[is.na(distgene)] <- 0
    distgene <- abs(distgene) 
    diag(distgene) <- 0 
    distgene1 <- distgene
    rm(distgene)
    genes <- colnames(datExpr1)
    dimnames(distgene1) <- list(genes,genes)
    
    load("distgene_3_1")
    distgene[is.na(distgene)] <- 0
    distgene <- abs(distgene) 
    diag(distgene) <- 0 
    distgene2 <- distgene
    rm(distgene)
    genes <- colnames(datExpr2)
    dimnames(distgene2) <- list(genes,genes)
    }
    
    genes <- intersect(colnames(datExpr1),colnames(datExpr2))
    distgene1 <- distgene1[genes,genes]
    distgene2 <- distgene2[genes,genes]
    
    distgene <- pmax(distgene1,distgene2)
    rm(distgene1)
    rm(distgene2)
    datExpr <- rbind(datExpr1[,genes],datExpr2[,genes])
    
    list(datExpr=datExpr,distgene=distgene)
}

combine_samples <- function(){

    dataFile <- c("GTEx_heart_expression.txt","GTEx.heart.countstable.txt","PCGC_newBorn_FPKM.txt","fetal_JUNE_AUG_expression.txt")
    
    # delete no expression genes
    datExpr <- read.delim(dataFile[1],sep="\t")
    genes <- datExpr[,1]
    samples <- colnames(datExpr)[-1]
    datExpr <- as.matrix(datExpr[,-1])
    datExpr <- t(datExpr)
    colnames(datExpr) <- genes  
    datExpr <- datExpr[,colSums(abs(datExpr))>0]
    datExpr1 <- datExpr
    
    datExpr <- read.delim(dataFile[3],sep="\t")
    genes <- datExpr[,1]
    samples <- colnames(datExpr)[-1]
    datExpr <- as.matrix(datExpr[,-1])
    datExpr <- t(datExpr)
    colnames(datExpr) <- genes  
    datExpr <- datExpr[,colSums(abs(datExpr))>0]   
    datExpr2 <- datExpr
    
    genes <- intersect(colnames(datExpr1),colnames(datExpr2))
    datExpr <- rbind(datExpr1[,genes],datExpr2[,genes])
    
    # PCA plot
    sampled <- dist(datExpr)
    loc <- prcomp(datExpr)
    plot(loc$x[,1],loc$x[,2],main="PCA plot",xlab="PC1",ylab="PC2")
    
    gsg <- goodSamplesGenes(datExpr,verbose=3)    # print(gsg$allOK)
    if (!gsg$allOK){ datExpr = datExpr[gsg$goodSamples, gsg$goodGenes];}
    
    sampleTree = flashClust(sampled,method="average")
    
    height <- 52000
    plot(sampleTree,main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
    abline(h=height,col="red")
    
    clust = cutreeStatic(sampleTree, cutHeight = height, minSize = 10) ### choose the cutHeight
    
    # split into two group expression values
    keepSamples <-  clust == 1
    datExpr <- datExpr[keepSamples,]
    datExpr <- datExpr[,colSums(abs(datExpr))>0]
    
    datExpr
}