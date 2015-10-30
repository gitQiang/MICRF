# Network analysis 
Network_analysis_batch <- function(){
  source("Network_analysis.R")
#   for(i in 1:4){
#     j <- 1
#       Network_analysis(i, j);
#   }
  for(i in 3:4){
    j=2
    Network_analysis(i, j);
  }
  
}

Network_analysis <- function(fileflag = 1, distflag = 1){

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
  # 0) PCA analysis and normalized methods
  dataFile <- c("GTEx_heart_expression.txt","GTEx.heart.countstable.txt","PCGC_newBorn_FPKM.txt","fetal_JUNE_AUG_expression.txt")
  batchFile <- "renameFileWithHeartTissue.list"
  # if use, housekeeping genes File "HK_genes.txt"
  datExpr <- step_1_0_normalize(dataFile[fileflag],batchFile,fileflag)

  # 1) co-expression based on pearson coefficient 
  # 2) based on mutual information
  if(distflag==1){
    distfile <- paste("distgene",fileflag,distflag,sep="_")
    load(distfile)
    distgene[is.na(distgene)] <- 0
    distgene <- abs(distgene) # this for DAWN
  }else if(distflag==2){
    distfile <- paste("dist",fileflag,distflag,sep="_")
    load(distfile)  
    distgene <- distM
    rm(distM)
    distgene[is.na(distgene)] <- 0
  }

  genes <- colnames(datExpr)
  dimnames(distgene) <- list(genes,genes)
# ---- 

# step 2: network check
# ----
  TFFile <- "TF-co-regulate.csv"
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
  #colSums(distgene[,knowngenes]>0)
  tmp <- matrix(0,n.known,n.known)
  meantmp <- matrix(0,n.known,n.known)
  for(i in 1:n.known){
    for(j in 1:n.known){
      tmp[i,j] <- t.test(as.vector(distgene[,knowngenes[i]]),as.vector(distgene[,knowngenes[j]]))$p.value
      meantmp[i,j] <- mean(as.vector(distgene[,knowngenes[i]]))/mean(as.vector(distgene[,knowngenes[j]]))
    }
    }
  tmp[c(1,3,2,8),c(1,3,2,8)]
  meantmp[c(1,3,2,8),c(1,3,2,8)]
  distgene[knowngenes,knowngenes]

  # known TF degree
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
  t.test(as.vector(distgene[,x]),as.vector(distgene[,y]))


if(distflag == 1){
  cut <- 0.6 # ????
  distgene[distgene < cut] <- 0
  diag(distgene) <- 0 
}else if(distflag == 2){
  cut <- 0.126
  distgene[distgene < cut] <- 0
  diag(distgene) <- 0 
  distgene <- distgene/max(distgene)
}

# delete isolated node
x <- colSums(distgene>0)
distgene <- distgene[x>0,x>0]
genes <- colnames(distgene)
datExpr <- datExpr[,genes]

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
  
  moduleLabels <- WGCNA_module(distgene,datExpr,distflag)
  # delete the module larger than 5000
  moduleLabels <- delete_module(moduleLabels,n=5000,pvalue_data[,1],genes)
  module <- data.frame(GENE=genes,MODULE=moduleLabels)
  diag(distgene) <- 1
  source("DAWN_sourcefile.R")
  datExpr <- t(datExpr)
  filename <- paste("DAWNresult_",fileflag,"_",distflag,".csv",sep="")
  DAWN_package(expression_data=datExpr,module_data=module, pvalue_data=pvalue_data,cor_thres1=0.75,cor_thres2=0.7,nASDthres=0.5,fdrthres=0.1,flag_q=distflag,filename=filename,distgene=distgene)
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
      batch <- as.matrix(read.delim(batchFile,header=FALSE,sep="."))
      samplenames <- batch[,c(1,2,5)]
      samplenames[samplenames[,3]=="Ventricle",3] <- 1
      samplenames[samplenames[,3]=="Appendage",3] <- 2
      samplenames[,1] <- paste(samplenames[,1],1,sep="")
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
  
  datExpr1 <- datExpr[samplenames[samplenames[,3]==1,1],]
  subs <- colSums(abs(datExpr1))>0
  datExpr1 <- datExpr1[,subs]
  # dim(datExpr1)
  datExpr2 <- datExpr[samplenames[samplenames[,3]==2,1],]
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

WGCNA_module <- function(distgene,datExpr,distflag){
  
  # WGCNA to module relationships
  # when use the pearson correlation, it should be changed to power 6 for WGCNA
  if(distflag ==1 ){
    softPower = 6;
    adjacency = adjacency(datExpr, power = softPower);  
  }else if(distflag==2){
    diag(distgene) <- 1
    adjacency = adjacency.fromSimilarity(distgene, type = "distance")
  }
  
  # set the minimum module size relatively high:
  minModuleSize = 30;
  MEDissThres = 0.15 ## A cut off value, We choose a height cut of 0.25, corresponding to correlation of 0.75, to merge
  
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
    if(tmp[i] > 5000 & as.numeric(tmpna[i])!=0){
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
