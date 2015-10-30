addresult <- function(){
    ### network degree distributions
    source("Multi_net.R")
    source("CRF_build.R")

    plotf=2    
    pdf(file="plot/degree_dis.pdf",width=7,height=7)
    
    k <- 1
    for(netflag in c(3,6,7,8,24)){
        
    net3 <- build_net(netflag,fileexp="")
    dev3 <- colSums(net3$matrix>0)
    
    kmax <- max(dev3)
    kmin <- min(dev3)
    de <- kmin:kmax
    Fv3 <- rep(0,kmax-kmin+1)
    for(i in kmin:kmax){
        Fv3[i-kmin+1] <- sum(dev3==i)
    }
    if(plotf==1){
        if(k==1){
            plot(Fv3[Fv3>1],col=k,type="l",xlab="degree",ylab="Frequency")
            k <- k +1
        }else{
            lines(Fv3[Fv3>1],col=k,type="l")
            k <- k +1
        }
    }else if(plotf==2){
        x <- de[Fv3>1]
        if(k==1){
            plot(log10(x),log10(Fv3[Fv3>1]),col=k,type="l",ylim=c(0,4),xlab="log10(degree)",ylab="log10(Frequency)")
            k <- k +1
        }else{
            lines(log10(x),log10(Fv3[Fv3>1]),col=k,type="l")
            k <- k +1
        }
    }
    
    }
    if(plotf==1){
        legend("topright",legend=c("STRING","iRefindex","Co-exp","Infer","GNAT"),lty=rep(1,5),lwd=rep(1,5),col=1:5)
    }else if(plotf==2){
        legend("bottomleft",legend=c("STRING","iRefindex","Co-exp","Infer","GNAT"),lty=rep(1,5),lwd=rep(1,5),col=1:5)
    }
    
    dev.off()
    
    
    ### GATA4, TBX5, and NKX2-5
    
    ##a. chromatin modification genes and their targets.  Similarly, transcription factors and their targets. 
    ##b. more detailed mRNA processing network, including spliceosome interaction and regulator-targets relationship. 
    k <- 1
    dlist <- list()
    glist <- list()
    nlist <- list()
    for(netflag in c(3,6,7,8)){
        HMGgene <- read.csv("data/genedata/HMGs.csv",header=FALSE)[,2] 
        chs <- c("GATA4","GATA6","TBX5","NKX2-5")
        if(netflag==3){
            mapT <- as.matrix(read.delim("Fmap0121.txt",header=FALSE,sep="\t"))
            HMGgene <- mapT[match(HMGgene,mapT[,2]),1]
            HMGgene <- HMGgene[!is.na(HMGgene)]
            chs <- mapT[match(chs,mapT[,2]),1]
        }
        net3 <- build_net(netflag,fileexp="")
        print(net3$size)
        chs <- intersect(chs,net3$node)
        #print(net3$matrix[chs,chs])
        dev3 <- colSums(net3$matrix>0)
    
        subs <- match(HMGgene,net3$node)
        subs <- subs[!is.na(subs)]
        
        dlist[[k]] <- dev3[subs]
        glist[[k]] <- net3$node[subs]
        k <- k+1
        
        targs <- net3$node[colSums(net3$matrix[subs,]>0)>=1]
        glist[[k]] <- targs
        dlist[[k]] <- dev3[match(targs,net3$node)]
        nlist[[k/2]] <- net3$node
        k <- k +1
    }
    
    Mtar <- rep(0,4)
    for(i in 1:4){
        Mtar[i] <- length(dlist[[2*i]])/length(dlist[[2*i-1]])
    }
    pdf(file="plot/hmgs_target.pdf",width=7,height=7)
    boxplot(list(dlist[[1]],dlist[[3]],dlist[[5]],dlist[[7]]),ylim=c(0,100),names=c("STRING","iRefindex","Co-exp","Infer"),ylab="The number of target genes",main="The number of target genes for HMG")
    dev.off()
    
    
    glist[[1]] <- mapT[match(glist[[1]],mapT[,1]),2]
    glist[[2]] <- mapT[match(glist[[2]],mapT[,1]),2]
    nlist[[1]] <- mapT[match(nlist[[1]],mapT[,1]),2]
    for(i in 1:4){
        write.table(glist[[2*i]],file=paste("target",i,sep=""),quote=FALSE,col.names=F,row.names=F)
    }
    brainlist <- c()
    source("Network_analysis.R")
    brainpagen <- read.delim("data/brain_specific/PaGen.txt")
    brainpagen <- unique(brainpagen[grepl("GDS707",brainpagen[,2]),1])
    brainpagen <- mapping_to(brainpagen)
    brainlist[[1]] <- brainpagen
    
    
    specond <- read.csv("data/brain_specific/SpeCond_Supple.csv")
    specond <- specond[specond[,"Adrenal_cortex"]>0,1] ##Whole_brain
    specond <- mapT[match(specond,mapT[,1]),2]
    specond <- specond[!is.na(specond)]
    specond <- mapping_to(specond)
    brainlist[[2]] <- specond

    TiGER1 <- read.delim("data/brain_specific/TiGER_genes.txt")
    TiGER1 <- unique(TiGER1[,"Gene_Symbol"])
    TiGER1 <- mapping_to(TiGER1)
    brainlist[[3]] <- TiGER1
    
    TiGER2 <- read.delim("data/brain_specific/TiGER_tf.txt")
    TiGER2 <- union(TiGER2[,1],TiGER2[,2])
    TiGER2 <- mapping_to(TiGER2)
    brainlist[[4]] <- TiGER2
    
    P <- matrix(0,4,4)
    for(i in 1:4){
        for(j in 1:4){
            P[i,j] <- (length(intersect(brainlist[[i]],glist[[2*j]]))/length(glist[[2*j]])) / (length(intersect(brainlist[[i]],nlist[[j]]))/length(nlist[[j]]))
        }
    }

  
}

gold_set <- function(){
    source("enrichana_6_12.R")
    options(stringsAsFactors=FALSE)
    TPcut=0.1
    options(stringsAsFactors=FALSE)
    netstr <- c("STRING/","iRef/","coexp/","Infer/","GNAT/")
    netnum <- c(31,6,7,8,21)
    TADAFile <- "ASD/TADAresult/TADAdenovo_ASD4_16.csv"
    asdall <- read.csv(TADAFile)
    Tset <- asdall[asdall[,"qvalue.dn"]<=TPcut,1] 
    
    g1 <- read.csv("ASD/nature10989/topriskg.csv",header=FALSE)
    g1 <- unique(g1[,3])
    g1 <- gsub(" ","",g1)
    
    g2 <- read.csv("ASD/nature13772/TADA107.csv")
    g2 <- g2[g2[,"qvalue"] < 0.3,1]
    
    g3 <- read.csv("ASD/cell/pASD.csv")
    g3 <- g3[g3[,"category"]=="pASD","gene"]
    
    g4 <- read.csv("ASD/cell/hcASD.csv")
    g4 <- unique(g4[,"gene"])
    g4 <- setdiff(g4,"")
    
    vg <- union(g1,union(g2,union(g3,g4)))
    
    length(intersect(Tset,vg))
    length(intersect(Tset,g1))
    length(intersect(Tset,g2))
    length(intersect(Tset,union(g3,g4)))
    
    ### recurrent in DDD 
    metam1 <- read.csv("DDD_mutations/datasheet/mutation_4_15_table1.csv",skip=1)
    metam1 <- metam1[,c(1,6,7,10,11)]
    metam1[,2] <- metam1[,2] + metam1[,4]
    metam1[,3] <- metam1[,3] + metam1[,5]
    metam1 <- metam1[,1:3]
    metam2 <- read.csv("DDD_mutations/datasheet/DDDmutation_4_15_table1.csv",skip=1)
    genes <- intersect(metam2[,1],metam1[,1])
    metam1[match(genes,metam1[,1]),2] <- metam1[match(genes,metam1[,1]),2] + metam2[match(genes,metam2[,1]),2]
    metam1[match(genes,metam1[,1]),3] <- metam1[match(genes,metam1[,1]),3] + metam2[match(genes,metam2[,1]),3]
    metam <- metam1
    
    pg <- setdiff(Tset,vg)
    tmp <- rowSums(metam[match(pg,metam[,1]),2:3])>0
    sum(tmp)
    
    vg1 <- pg[tmp]
    vg2 <- pg[!tmp]
    
    write.table(vg2,file="vg2.txt",quote=FALSE,col.names=FALSE,row.names=FALSE)
    
    source("enrichana.R")
    tmpl <- test_genes()
    testgenes <- tmpl$testgenes
    
    a<- rep(0,length(testgenes))
    for(i in 1:length(testgenes)){
        a[i] <- length(intersect(testgenes[[i]],vg2))
        print(intersect(testgenes[[i]],vg2))
    }
    
    vg3 <- setdiff(vg2,c("CACNA1S","RYR1","KMT2C","DSCAM","KMT2E","TNRC6B","ZC3H4","DSCAM"))
    
    #g2 <- read.csv("ASD/nature13772/TADA107.csv")
    nu1 <- rowSums(asdall[match(vg3,asdall[,1]),c("dn.LoF","dn.mis3")])
    #nu2 <- rowSums(g2[match(vg3,g2[,1]),c("dn.LoF","dn.mis3")])
    #nu1 <- nu1 - nu2 
    nu1
    
}

# figure 4: co-exp better than PPIs 
co_exp_PPIs <- function(){
    source("~/.Rprofile")
    source("enrichana_6_12.R")
    library(ROCR)
    library(pROC)
    
    load("BestP_6_10")
    betaV <- BestP[,3]
    
    filenames <- filenamesf(0,0,betaV=betaV,"","")
    allP <- allPf(TPcut=0.1)
    Tset <- allP
    
    DAWNflag <- 1
    DAWNnames <-  readLines(con <- file("result/randresult_1/DAWN/TADAdenovo_randset_1.txt","r"))   
    close(con)
    Maginames <-  readLines(con <- file("result/randresult_1/MAGI/MAGIfilenames.txt","r"))   
    close(con)
    
    TADAg <- list()
    NetPg <- list()
    NetiRef <- list()
    TADAgs <- c()
    NetPgs <- c()
    NetiRefg <- c()
    
    for(j in 2){
        for(i in 1:20){   
            filenames <- filenamesf(j,i,betaV=betaV,DAWNnames,Maginames)
            tmp <- gene_select(filenames,Tset)
            genes <- tmp$genes
            allP <- tmp$allP
            allN <- tmp$allN
            
            if(DAWNflag==1){
                tmp <- read.csv(filenames[2])
                genes <- intersect(genes,tmp[!is.na(tmp[,"FDR"]),"Gene"])
            }
            
            for(k in c(1,3,10,11)){
                filename <- filenames[k]
                cutoff <- total_genes_cutoff(filename,Tset,muts="",alpha=0.5)
                
                if(grepl(".csv",filename)){ result <- read.csv(filename); 
                }else if(grepl(".txt",filename)){result <- read.table(filename);
                }else{result <- read.table(filename);}
                
                if(k==1){
                    ones <- result[result[,"qvalue.dn"] <= 1-cutoff,1]
                    TADAg[[i]] <- ones
                    TADAgs <- union(TADAgs,ones)
                }else if(k==11){
                    ones <- result[result[,2]  >= cutoff,1]
                    NetPg[[i]] <- ones
                    NetPgs <- union(NetPgs,ones)
                }else if(k==10){
                    ones <- result[result[,2]  >= cutoff,1]
                    NetiRef[[i]] <- ones
                    NetiRefg <- union(NetiRefg,ones)            
                }
            }
            
        }
    }
    
    a <- list()
    a[[1]] <- rep(0,20)
    a[[2]] <- rep(0,20)
    a[[3]] <- rep(0,20)
    for(i in 1:20){
        a[[1]][i] <- 100*length(intersect(TADAg[[i]],Tset))/length(TADAg[[i]])
        a[[2]][i] <- 100*length(intersect(NetiRef[[i]],Tset))/length(NetiRef[[i]])
        a[[3]][i] <- 100*length(intersect(NetPg[[i]],Tset))/length(NetPg[[i]])
    }
    
    b <- matrix(0,20,3)
    for(i in 1:20){
        b[i,1] <- length(intersect(TADAg[[i]],Tset))/length(TADAg[[i]])
        b[i,2] <- length(intersect(NetiRef[[i]],Tset))/length(NetiRef[[i]])
        b[i,3] <- length(intersect(NetPg[[i]],Tset))/length(NetPg[[i]])
    }
    
    pdf(file="plot/TsetFrac.pdf",width=10,height=7)
    main="Fraction of gold risk genes"
    boxplot(a,names=c("TADA","NetP-iRef","Net-coexp"),main="",ylab="Percentage",outline=FALSE)
    dev.off()
    
    #b <- matrix(c(length(intersect(Tset,NetPgs)),length(NetPgs),length(intersect(Tset,TADAgs)),length(TADAgs)),2,2)
    #fisher.test(b,alternative="greater")
    
    tmp <- matrix(0,20,4)
    for(i in 1:20){
        tmp[i,] = c(length(intersect(NetPg[[i]],Tset)),length(intersect(NetiRef[[i]],Tset)),a[[3]][i],  a[[2]][i]);
    }
    subs <- which(tmp[,1]> tmp[,2] & tmp[,3] >=  tmp[,4])
    
    i=20
    x1 <- intersect(NetPg[[i]],Tset)
    x2 <- intersect(NetiRef[[i]],Tset)
    ge = setdiff(x1,x2);
    qwt(setdiff(x1,x2),file="plot/addedgenesL.txt")
    qwt(cbind(setdiff(x1,x2),"=",3),file="plot/addedgenes.txt")
    
    input <- read.csv("ASD/TADAresult/randset_1/TADAdenovo_ASDrand1_2_20.csv")
    tmp <- cbind(input[,1],"=",-log10(input[,"pval.TADA.dn"]))
    tmp[tmp[,1] %in% ge,3] <- 10
    
    net  <- read.table("data/network_inference/ASDexp_net_top5.txt")
    net1  <- read.table("data/hotnet/iRefIndexm.txt")
    net <- rbind(net,net1)
    netg <- union(net[,1],net[,2])
    tg <- setdiff(netg,tmp[,1])
    tg <- cbind(tg," = ",0)
    tmp <- rbind(tmp,tg)
    qwt(tmp,file="plot/input2_20.txt")
    
    ## single edge for co-expression
    source("CRF_build.R")
    a <- read.table("data/network_inference/ASDexp_net_top5.txt")
    net <- read_net(a)
    
    adj <- net$matrix
    adj[adj > 0] <- 1
    adj[lower.tri(adj,diag=TRUE)] <- 0
    edges <- which(adj>0,arr.ind=TRUE)
    net1 <- cbind(net$node[edges[,1]],net$node[edges[,2]],net$matrix[edges])
    qwt(net1,file="data/network_inference/ASDexp_net_top5_single.txt")
    
    
    #==================================================================================
    #####
    ### average network distances
    library(igraph)
    source("Multi_net.R")
    source("CRF_build.R")
            
    net <- build_net(netflag=7,fileexp="")
    adj <- net$matrix 
    adj[adj>0] <- 1
    g1 <- graph.adjacency(adj, mode="undirected", weighted=NULL, diag=FALSE) 
    D1 <- shortest.paths(g1)
    
    net6 <- build_net(netflag=6,fileexp="")
    adj6 <- net6$matrix
    adj6[adj6>0 ] <- 1
    g6 <- graph.adjacency(adj6, mode="undirected", weighted=NULL, diag=FALSE)
    D6 <- shortest.paths(g6)
    
    qwt(cbind(intersect(NetPgs,Tset),"=",1),file="plot/coexp_risk.txt")
    qwt(cbind(intersect(NetiRefg,Tset),"=",1),file="plot/iRef_risk.txt")    
    
    load("iRefmodule")
    qwt(cbind(module[,1],"=",module[,2]),file="plot/iRef_module.txt") 
    net.text <- read.table("data/hotnet/iRefIndexm.txt")
    modL <- matrix(0,dim(net.text)[1],2)
    modL[,1] <- module[match(net.text[,1], module[,1]),2]
    modL[,2] <- module[match(net.text[,2], module[,1]),2]
    
    subs <- modL[,1]==modL[,2]
    net.text <- net.text[subs,]
    qwt(net.text,file="plot/modulenet_iRef.txt")
    
    load("DAWNmodule")
    ####
    
    dlist <- list()
    tg1 <- intersect(Tset,net$node)
    subs1 <- match(tg1,net$node)
    dg1 <- D[subs1,subs1]
    dlist[[1]] <- as.vector(dg1[dg1>0 & dg1<Inf])
    
    dlist[[2]] <- rep(0,20)
    dlist[[3]] <- rep(0,20)
    for(i in 1:20){
    tg2 <- intersect(TADAg[[i]],net$node)
    subs2 <- match(tg2,net$node)
    dg2 <- D[subs2,subs2]
    dlist[[2]][i] <- mean(as.vector(dg2[dg2>0 & dg2<Inf]))
    
    tg3 <- intersect(NetPg[[i]],net$node)
    subs3 <- match(tg3,net$node)
    dg3 <- D[subs3,subs3]
    dlist[[3]][i] <- mean(as.vector(dg3[dg3>0 & dg3<Inf]))
    }
    
    a <- list()
    a[[1]] <- dlist[[2]]
    a[[2]] <- dlist[[3]]
    boxplot(a)
   
    #####====
    ###Distances for Tsets
    library(igraph)
    source("Multi_net.R")
    source("CRF_build.R")
    mapT <- as.matrix(read.delim("Fmap0121.txt",header=FALSE,sep="\t"))
    Dlist <- list()
    netflags <- c(3,6,7,8)
    for(i in 1:4){
        net <- build_net(netflag=netflags[i],fileexp="")
        adj <- net$matrix 
        adj[adj>0] <- 1
        g1 <- graph.adjacency(adj, mode="undirected", weighted=NULL, diag=FALSE) 
        D <- shortest.paths(g1)
        
        if(i==1){
            net$node <- mapT[match(net$node,mapT[,1]),2]
        }
        
        tg1 <- intersect(Tset,net$node)
        subs1 <- match(tg1,net$node)
        dg1 <- D[subs1,subs1]
        dg1[upper.tri(dg1,diag=TRUE)] <- 0
        Dlist[[i]] <- as.vector(dg1[dg1>0 & dg1<Inf])
    }
    pdf(file="plot/Tset_distances.pdf",width=7,height=7)
    boxplot(Dlist,names=c("STRING","iRefindex","Co-exp","Infer"),ylab="Distances")
    dev.off()
    #####
    #==========================
}

degree_be <- function(){
    source("Multi_net.R")
    source("CRF_build.R")
    source("enrichana_6_12.R")
    library(ROCR)
    library(pROC)
    
    load("BestP_6_10")
    betaV <- BestP[,3]
    
    filenames <- filenamesf(0,0,betaV=betaV,"","")
    allP <- allPf(TPcut=0.1)
    Tset <- allP
    
    k <- 1
    dlist <- list()
    netflag=7

    net3 <- build_net(netflag,fileexp="")
    dev3 <- colSums(net3$matrix>0)
    
    subs <- match(Tset,net3$node)
    subs <- subs[!is.na(subs)]
    
    dlist[[k]] <- dev3[subs]   
    k <- k +1
    
    targs <- setdiff(net3$node,Tset)
    dlist[[k]] <- dev3[match(targs,net3$node)]
    pdf(file="plot/degree_Tset.pdf",width=7,height=7)
    boxplot(dlist,ylim=c(0,40),names=c("Gold risk","Others"),ylab="Degree")    
    dev.off()
    
    d1 <- list()
    be3 <- read.table("data/Network_betweenness/Betweenness_node_coexp.txt")
    d1[[1]] <- be3[match(intersect(Tset,be3[,1]),be3[,1]),2]
    d1[[2]] <- be3[match(setdiff(be3[,1],Tset),be3[,1]),2]
    pdf(file="plot/between_Tset.pdf",width=7,height=7)
    boxplot(d1,ylim=c(0,5000),names=c("Gold risk","Others"),ylab="Betweenness centrality")
    dev.off()
}

# figure 1: sparse the mutation
global_map <- function(){

    ### 3953 ASD trios data and ASD co-expression network
    mut <- read.csv("ASD/TADAresult/TADAdenovo_ASD4_16.csv")
    num <- mut[,"dn.LoF"] + mut[,"dn.mis3"] + mut[,"dn.mis"]
    mutnum1 <- cbind(mut[,1],"=",num)
    
    net <- read.table("data/network_inference/ASDexp_net_top5.txt")
    genes1 <- union(net[,1],net[,2])
    genes1 <- setdiff(genes1,mut[,1])
    mutnum2 <- cbind(genes1,"=",0)

    qwt(rbind(mutnum1,mutnum2),file="ASD/muation_num.txt")
    
    pdf(file="plot/hist_mut.pdf",width=12,height=7)
    hist(num[num>0],xlab="De novo count",ylab="The number of genes",col=2)
    #axis(1, at=0:15, labels=0:15)
    dev.off()
    
    ### nature 13908 and iREF PPI network
    deset2 <- read.csv("ASD/nature13908/Supplementary Table 2.csv")
    inchild <- c("pF","pFsF","pFsM","pM","pMsF","pMsM")
    deset2 <- deset2[deset2[,"inChild"] %in% inchild,]
    lofset2 <- c("frame-shift","nonsense","splice-site") # LGD mutations: likely gene disrupting: nonsense, frameshift and splice site
    misset2 <- c("missense")
    deset2 <- deset2[deset2[,"effectType"] %in% union(lofset2,misset2),]
    deset2[,"effectGene"] <- sapply(deset2[,"effectGene"],function(s) {unlist(strsplit(s,","))[1]})
    #source("Network_analysis.R")
    #deset2[,"effectGene"] <- mapping_to(deset2[,"effectGene"])
    
    tmp <- table(deset2[,"effectGene"])
    mutnum1 <- cbind(names(tmp),"=",tmp)
    
    a <- rep(0,8)
    for(i in 1:8){a[i] <- sum(tmp==i);}
    
    pdf(file="plot/sparsity.pdf",width=12,height=10)
    #hist(tmp,xlab="De novo count",ylab="The number of genes",col=2)
    par(cex.lab=1.7,cex.main=2,mai=c(2,2,1,1))
    mp <- barplot(a,space=0.4, col=2,cex.axis=1.6,xlab="De novo count",ylab="The number of genes",main=substitute(paste("Sparsity of ", italic('de novo'), " mutations" )))
    axis(1, at=mp, labels=1:8,cex.axis=1.6)
    dev.off()   
   
    net <- read.table("plot/modulenet_iRef.txt")
    #mapT <- as.matrix(read.delim("Fmap0121.txt",header=FALSE,sep="\t"))
    #net[,1] <- mapT[match(net[,1],mapT[,1]),2]
    #net[,2] <- mapT[match(net[,2],mapT[,1]),2]
    #net <- net[!is.na(net[,1]) & !is.na(net[,2]),]
    
    #load("module3")
    #module[,1] <- mapT[match(module[,1],mapT[,1]),2]
    #module <- module[!is.na(module[,1]),]
    #qwt(paste(module[,1]," = ", module[,2], sep=""),file="ASD/SPRING_module.txt")
    
    #tgene <- mutnum1[,1]
    #mgene <- module[module[,2]==1,1]
    #togene <- union(tgene,mgene)
    #togene <- tgene
    
    #net <- net[net[,1] %in% togene | net[,2] %in% togene,]
    #qwt(net,file="ASD/netiref_nature13908.txt")
    
    genes1 <- union(net[,1],net[,2])
    genes1 <- setdiff(genes1,mutnum1[,1])
    mutnum2 <- cbind(genes1,"=",0)
    source("~/.Rprofile")
    qwt(rbind(mutnum1,mutnum2),file="ASD/muation_num_nature13908.txt")
}

# figure 2: power of the statistical test methods
statis_number <- function(){
    source("Poisson_test_hq.R")
    filename <- "ASD/TADAresult/TADAdenovo_ASD4_16.csv"
    ntrio <- 3953
    tmp <- Poisson_test_hq(filename,ntrio)
    tmp <- tmp[order(tmp[,"min_p"]),]
    tmp <- cbind(tmp,dim(tmp)[1]*3*tmp[,"min_p"])
    
    mut <- read.csv("ASD/TADAresult/TADAdenovo_ASD4_16.csv")
    
    cutV <- seq(0.0001,0.05,0.0005)
    y1 <- rep(0,length(cutV))
    y2 <- rep(0,length(cutV))
    y3=y1
    y4=y1
    
    for(i in 1:length(cutV)){
        y1[i] <- sum(mut[,"qvalue.dn"] <= cutV[i])
        y2[i] <- sum(tmp[,19] <= cutV[i])
        y3[i] <- sum(mut[,"pval.TADA.dn"] <= cutV[i])
        y4[i] <- sum(tmp[,"min_p"] <= cutV[i])
    }
    
    #plot
    pdf(file="plot/statis_number.pdf",width=13,height=6)
    par(mfrow=c(1,2))
    plot(-log10(cutV),y1,xlab="-log10 (p-value)",ylab="Number of significant genes",ylim=c(0,max(c(y1,y2))),main="Significance genes (Corrected p-value)",type="l")
    lines(-log10(cutV),y2,col=2,type="l")
    #abline(v=2)
    legend("topright",legend=c("TADA","Poisson Test"),lwd=c(1,1),lty=c(1,1),col=1:2)
    
    plot(-log10(cutV),y3,xlab="-log10 (p-value)",ylab="Number of significant genes",ylim=c(0,max(c(y3,y4))),main="Significance genes (p-value)",type="l")
    lines(-log10(cutV),y4,col=2,type="l")
    #abline(v=2)
    legend("topright",legend=c("TADA","Poisson Test"),lwd=c(1,1),lty=c(1,1),col=1:2)
    dev.off()
}

# figure 3: first simulation performance
first_performance <- function(){
    source("enrichana_6_12.R")
    #### selected plot
    nmeth=7
    subs <- c(1,2,3,5,6,10,11)
    cols <- c("Black","Blue","Cyan","Brown","Coral","DarkViolet","Green")
    
    
    pdf(file="plot/samplesize.pdf",width=10,height=8)
    
    load("AUC_m_6_11")
     #c("Black","Blue","Cyan","Brown","Coral","DarkViolet","Green","HotPink","Orange","Plum","deeppink")
    labels <- c("TADA","DAWN","MAGI","HotNet2-STRING","HotNet2-iRef","HotNet2-coexp","HotNet2-Infer","HotNet2-GNAT","MICRF-STRING","MICRF-iRef","MICRF-coexp","MICRF-Infer","MICRF-GNAT")
    #legend <- labels
    
    main="AUC performance" # for paper
    main="" # for slides
    xlab="The number of trios";ylab="AUC";
    ymint <- 1;ymaxt <- 0;
    ymint <- min(AUC);ymaxt <- max(AUC);
    ymin=max(min(ymint)+0.1,0);ymax = min(ymaxt+0.1,1);
    x <- seq(0.2,0.9,0.1)
    Y1 <- rep(0,nmeth)
    
    nn <- subs ###1:nmeth
    
    k=1
    for(i in nn){
        y <- rowMeans(AUC[,,i])
        sd <- apply( AUC[,,i], 1, sd)/2 
        plot_err(x[1:5],y[1:5],sd[1:5],k,main,xlab,ylab,ylim=c(ymin,ymax),xlim=c(0.0,0.65))
        #plot_err(x,y,sd,k,main,xlab,ylab,ylim=c(ymin,ymax),xlim=c(-0.05,0.95))
        Y1[k] <- y[1]
        k <- k+1
    }
    ## nmeth=5
    Y1[3] <- 0.851
    Y1[4] <- 0.845
    Y1[5] <- Y1[5] + 0.003
    ## nmeth=7
    Y1[6] <- Y1[6] - 0.0065
    
    text(rep(0.2,nmeth),Y1,labels=labels[nn],pos=2,col=cols[1:nmeth],cex=1.5)
    #legend("bottomright",col=cols[1:nmeth],legend=legend[nn],lty=rep(1,nmeth),lwd=rep(3,nmeth),cex=1.2,y.intersp=0.8)
    dev.off()

}

# Recurrent in two ASD data sets
ASD_re <- function(){
    
    source("CHD_normalized.R")
    source("Network_analysis.R")
    source("CRF_build.R")
    source("Multi_net.R")
    options(stringsAsFactors=FALSE)
    fileexp <- "PCGCall.txt" 
    
    netstr <- c("STRING/","iRef/","Infer/","coexp/","GNAT/")
    k = 4
    modulefiles <- c("module3","iRefmodule","DAWNmodule","DAWNmodule","data/GNATnet/module_brain")
    
    netflag=7
    if(netflag==3){ instr <- "CRF_input";idm=TRUE; }else{ instr <- "hotnet_input"; idm=FALSE;}
    modulefile <- modulefiles[k]
    
    dirstr <- "result/"
    strname <- "ASDnat"
    filename <- paste(dirstr,instr,strname,".txt",sep="")
    strn <- paste("result/randresult_1/",netstr[k],"CRFresult_",strname,sep="")
    Multi_net(netflag,strn,filename,fileexp,modulefile,cutn=200000,idm,beta=0.2,net_e=4);     
       
    source("enrichana_6_12.R")
    rfile <- "result/randresult_1/coexp/CRFresult_ASDnatLBP_7.txt"
    Tset <- unlist(read.table("autism.FDR-0.3.txt"))
    cutoff <- total_genes_cutoff(rfile,Tset,muts="",alpha=0.5)
    
    rtable <- read.table(rfile)
    pg <- rtable[rtable[,2]>=cutoff,1]
    
    
    ## nature 13908 unique
    options(stringsAsFactors=FALSE)
    ## read the families used in 13908, supplementary table 1
    fam13908 <- read.csv("ASD/nature13908/Supplementary Table 1.csv")[,1]
    fam13908 <- unique(fam13908)
    ## read the families used in other papers, supplementary table 10
    #fam13908subset <- read.csv("ASD/nature13908/Supplementary Table 10.csv") ## 208 only 1 less
    ## read the families from lossifov, 343 families
    famlos <- read.csv("ASD/Neuron/mmc2.csv")[,1]
    famlos <- unique(famlos)
    ## read the families from sanders nature 10945, 238 families
    famsan <- read.csv("ASD/nature10945/nature10945-s2.csv")[,1]
    famsan <- unique(famsan)
    ## read the familes from O' Roak nature 10989, 189 + 20  families
    deroa <- read.csv("ASD/nature10989/nature10989-s1.csv")
    deroa <- deroa[deroa[,"type"] %in% "SSC189",]
    famroa <- deroa[,1]
    famroa <- sapply(1:length(famroa),function(i) {unlist(strsplit(famroa[i],"\\."))[1]})
    famroa <- unique(famroa)
    
    fam13908F <- setdiff(fam13908,union(famlos,union(famsan,famroa)))
    length(fam13908F) ## 1751
    newfile <- "DDD_mutations/datasheet/anno_mutationlist_4_15_all.csv"
    allr <- read.csv(newfile)
    allr <- allr[allr[,"From"]=="nature 13908" & allr[,"mutation_type_metaSVM"] %in% c("LOF","dMIS") & allr[,"sampleID"] %in% fam13908F,]
    sum(pg%in% allr[,"Gene"])
    length(unique(allr[,"Gene"]))
    
}

# meta ana 
meta_ana <- function(){

    poifile <- "DDD_mutations/datasheet/Meta-analysis_result_7_3.csv"
    allresult <- read.csv(poifile)
    
    NetPg <- allresult[allresult[,2]>0.7,1]
    TADAg <- allresult[allresult[,"TADA_qvalue"] <= 0.05,1]
    Poig <- allresult[as.numeric(allresult[,"min_p"]) * 3 <= 0.05,1]
    
    knowngenes <- unlist(read.table("DDD_mutations/nature14135-s2/Genes.txt"))[1:34]
    knowngenes <- unique(knowngenes)
    length(knowngenes)
    
    length(intersect(NetPg,TADAg))
    length(intersect(NetPg,Poig))
    length(intersect(Poig,TADAg))
    length(intersect(Poig,intersect(NetPg,TADAg)))
    
    length(intersect(NetPg,knowngenes))
    length(intersect(TADAg,knowngenes))
    length(intersect(Poig,knowngenes))
    
    
    length(NetPg)
    length(TADAg)
    length(Poig)
    
    source("DDD.R")
    test_genesDDD(TADAg)
    
    a <- read.table("result/control_6_12/v4/coexp/CRFresult_0.2Meta_7_3LBP_7.txt")
    genelist <- a[a[,2]>0.7,1]
    qwt(genelist,file="plot/Meta_risk_genes.txt")
    
}