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
    
    pdf(file="plot/sparsity_10_8.pdf",width=12,height=10)
    #hist(tmp,xlab="De novo count",ylab="The number of genes",col=2)
    par(cex.lab=1.7,cex.main=2,mai=c(2,2,1,1))
    #mp <- barplot(a,space=0.4, col=2,cex.axis=1.6,xlab="De novo count",ylab="The number of genes",main=substitute(paste("Sparsity of ", italic('de novo'), " mutations" )))
    mp <- barplot(a,space=0.4, col=2,cex.axis=1.6,xlab="De novo count",ylab="Number of genes",main="")
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
    #filename <- "ASD/TADAdenovo_nat8_20.csv"
    filename <- "ASD/TADAresult/TADAdenovo_ASD4_16.csv"
    ntrio <- 3953
    tmp <- Poisson_test_hq(filename,ntrio)
    tmp <- tmp[order(tmp[,"min_p"]),]
    tmp <- cbind(tmp,p.adjust(tmp[,"min_p"],method="fdr"))
    
    mut <- read.csv("ASD/TADAdenovo_nat8_20.csv")
    mut[is.na(mut[,"qvalue.dn"]),"qvalue.dn"] <- 1
    cutV <- seq(0.01,0.3,0.005)
    y1 <- rep(0,length(cutV))
    y2 <- rep(0,length(cutV))
    y3=y1
    y4=y1
    
    for(i in 1:length(cutV)){
        y1[i] <- sum(mut[,"qvalue.dn"] <= cutV[i])
        y2[i] <- sum(tmp[,19] <= cutV[i])
        #y3[i] <- sum(mut[,"pval.TADA.dn"] <= cutV[i])
        #y4[i] <- sum(tmp[,"min_p"] <= cutV[i])
    }
    
    #plot
    pdf(file="plot/statis_number.pdf",width=8,height=8)
    par(mfrow=c(1,1),mai=c(2,1,1,1))
    plot(-log10(cutV),y1,xlab="p-value",ylab="Number of genes",ylim=c(0,max(c(y1,y2))),main="Significant genes",type="l",xaxt="n",cex.axis=1.6,cex.lab=1.7,cex.main=2)
    lines(-log10(cutV),y2,col=2,type="l")
    #abline(v=-log10(0.05))
    legend("topright",legend=c("TADA","Poisson Test"),lwd=c(1,1),lty=c(1,1),col=1:2,cex=1.3)
    axis(1, at=-log10(cutV), labels=cutV,cex.axis=1.6)
         
#     plot(-log10(cutV),y3,xlab="p-value",ylab="Number of significant genes",ylim=c(0,max(c(y3,y4))),main="Significance genes (p-value)",type="l",xaxt="n")
#     lines(-log10(cutV),y4,col=2,type="l")
#     #abline(v=-log10(0.05))
#     legend("topright",legend=c("TADA","Poisson Test"),lwd=c(1,1),lty=c(1,1),col=1:2)
#     axis(1, at=-log10(cutV), labels=cutV)
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

# Nature 13772 
Nat13772 <- function(){
    
#     source("Network_analysis.R")
#     a <- read.csv("ASD/TADAresult/TADA_results.csv")
#     colnames(a)[15] <- "qvalue.dn"
#     a[,1] <- mapping_to(a[,1])
#     write.csv(a,file="ASD/TADAresult/TADA_results_1.csv",row.names=FALSE)
    
#         source("Network_analysis.R")
#         a <- read.csv("ASD/TADAresult/TADAdenovo_Nat13772.csv")
#         a[,1] <- mapping_to(a[,1])
#         write.csv(a,file="ASD/TADAresult/TADAdenovo_Nat13772_1.csv",row.names=FALSE)
    
    source("enrichana_6_12.R")
    ### 107 gold risk genes
    #T107 <- unlist(read.table("ASD/nature13772/autism.FDR-0.3.txt"))
    
    tmp <- read.csv("ASD/TADAresult/TADAdenovo_Nat13772_1.csv")
    Tset <- tmp[tmp[,"qvalue.dn"] < 0.1,1]
    #Tset <- T107
    source("Network_analysis.R")
    Tset <- mapping_to(Tset)
    
    
    filenames <- c("ASD/TADAresult/TADAdenovo_Nat13772_1.csv","result/Nat1377_meta/DAWN/DAWN_nat13772.csv","result/Nat1377_meta/MAGI/RandomGeneList.4","result/Nat1377_meta/RWR/coexp/hotnetresult1ASDnat137727.txt","result/Nat1377_meta/MICRF/v4/coexp/CRFresult_0.2ASDnat13772LBP_7.txt")
    #filenames <- c("ASD/TADAresult/TADA_results_1.csv","result/Nat1377_meta/DAWN/DAWN_nat13772_0.csv","result/Nat1377_meta/MAGI/RandomGeneList.4","result/Nat1377_meta/RWR/coexp/hotnetresult1nat13772_07.txt","result/Nat1377_meta/MICRF/v4/coexp/CRFresult_0.2nat13772_0LBP_7.txt")
    
    
    geneT <- list()
    genes <- c()
    for(i in 1:length(filenames)){
        filename <- filenames[i]
        cutoff <- total_genes_cutoff(filename,Tset,muts="",alpha=0.5)
        
        if(grepl(".csv",filename)){ result <- read.csv(filename); 
        }else{result <- read.table(filename);}
        
        if(grepl("DAWN",basename(filename))){
            result <- result[!is.na(result[,"FDR"]),]
        }
        
        if(i==1){ geneT[[i]] <- result[result[,"qvalue.dn"] <= 1-cutoff,1];}else if(i==2){ geneT[[i]] <- result[result[,"Stage2_posterior"] <= 1-cutoff,1];
        }else{
            geneT[[i]] <- result[result[,2]>= cutoff,1]
        } 
        
        genes <- c(genes,geneT[[i]])
    }
    
    ##========== recurrent in ASD
    n13908 <- read.csv("result/Nat1377_meta/nature13908.csv")
    a <- 1:5
    b <- 1:5
    for(i in 1:5){
        a[i] <- sum(geneT[[i]] %in% n13908[,"Gene"])/length(geneT[[i]])
        b[i] <- sum(intersect(geneT[[i]],Tset) %in% n13908[,"Gene"])/sum(Tset %in% n13908[,"Gene"])
    }
    
    a
    tmpg <- setdiff(geneT[[5]],geneT[[2]])
    tmpg1 <- setdiff(geneT[[2]],geneT[[5]])
    sum(tmpg %in% n13908[,"Gene"])/length(tmpg)
    sum(tmpg1 %in% n13908[,"Gene"])/length(tmpg1)
    
    
    ## poisson test for Nature 13772
    poiT <- read.csv("ASD/TADAresult/PtestNat13772.csv")
    gpoi <- poiT[poiT[,"p_value"] < 0.05,1]
    
    length(intersect(tmpg,gpoi))
    length(intersect(tmpg1,gpoi))
    
    length(tmpg)
    length(tmpg1)
    length(intersect(geneT[[2]],geneT[[5]]))
    
    ##=============================================
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
    
    tmp <- rowSums(metam[match(tmpg1,metam[,1]),2:3])>0
    sum(tmp)
    
    
    # showing in Cytoscape ======================== 
    source("~/.Rprofile")
    mes <- c("TADA","DAWN","MAGI","RWR","MICRF")
    tn <- table(genes)
    tg <- names(table(genes))
    
    tT <- c()
    for(i in 1:length(tg)){
        if(tn[i]==1){
            subs <- rep(FALSE,5)
            for(j in 1:5){
                subs[j] <- tg[i] %in% geneT[[j]]
            }
            tme <- mes[subs]
        }else{
            tme <- "Overlap"
        }
        tmp <- c(tg[i]," = ",tme)
        tT <- rbind(tT,tmp)
    }
    #qwt(tT,file="result/Nat1377_meta/gene_method.txt")
    
    input <- read.csv("ASD/TADAresult/TADA_results.csv")
    tmp <- cbind(input[,1],"=",-log10(input[,"pval.TADA"]))
    #input <- read.csv("ASD/TADAresult/TADAdenovo_Nat13772.csv")
    #tmp <- cbind(input[,1],"=",-log10(input[,"pval.TADA.dn"]))
    
    net  <- read.table("data/network_inference/ASDexp_net_top5_single.txt")
    netg <- union(net[,1],net[,2])
    tg <- setdiff(netg,tmp[,1])
    tg <- cbind(tg," = ",0)
    tmp <- rbind(tmp,tg)
    #qwt(tmp,file="plot/inputnat13772.txt")
    
    ### check the connected in co-expression network
    source("CRF_build.R")
    net.text <- read.table("data/network_inference/ASDexp_net_top5.txt")
    net <- read_net(net.text)
    deV <- colSums(net$matrix>0)
    deV[match(tmpg,net$node)]
    
    deV[match(tmpg1,net$node)]
        
}

batch_poisson_hq <- function(){
    
    source("Poisson_test_hq.R")
    filename <- "ASD/TADAresult/TADAdenovo_Nat13772.csv"
    N <- 2270
    tmp <- read.csv(filename)
    n <- dim(tmp)[1]
    pV <- sapply(1:n, function(i){
        poisson.test(x=tmp[i,"dn.LoF"]+ tmp[i,"dn.mis3"], T = N *2* tmp[i,"mut.rate"], alternative = "greater", conf.level = 0.95)$p.value
    })
    tmp[,"p_value"] <- pV 
    tmp[,"adjust_p"] <- 1
    tmp[,"adjust_p"] <- p.adjust(tmp[,"p_value"],method="fdr") ###!!!!!! n test number
    
    write.csv(tmp,file="ASD/TADAresult/PtestNat13772.csv",row.names=FALSE)
    
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
    
    write.csv(allr,file="result/Nat1377_meta/nature13908.csv",row.names=FALSE)
}

## output performance: high ranked gene lists
high_ranked <- function(){
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
    
    resultg <- list()
    m=1
    for(j in 2){
        for(i in 1:20){   
            filenames <- filenamesf(j,i,betaV=betaV,DAWNnames,Maginames)
            
            for(k in c(1,2,3,6,11)){
                filename <- filenames[k]
                cutoff <- total_genes_cutoff(filename,Tset,muts="",alpha=0.5)
                
                if(grepl(".csv",filename)){ result <- read.csv(filename); 
                }else{result <- read.table(filename);}
                
                if(k==1){
                    ones <- result[result[,"qvalue.dn"] <= 1-cutoff,1]
                }else if(k==2){
                    subs <- is.na(result[,"FDR"])
                    result <- result[!subs,]
                    ones <- result[result[,"Stage2_posterior"]  <= 1-cutoff,1]
                }else if(k>2){
                    ones <- result[result[,2]  >= cutoff,1]          
                }
             
             resultg[[m]] <- ones   
             m  <- m+1
            }
        }
    }

    TL <- matrix(0,20,5)
    for(k in 1:100){
        i <- floor((k-1)/20) + 1
        j <- k %% 20
        if(j==0) j=20
        TL[j,i] <- length(intersect(resultg[[k]],Tset))/length(resultg[[k]])
    }
    
    boxplot(TL)
    
    ## example: i=5
    exL <- list()
    exL1 <- list()
    for(i in 1:5){
        exL[[i]] <- setdiff(intersect(resultg[[(i-1)*20+5]],Tset), intersect(resultg[[85]],Tset))
        exL1[[i]] <- setdiff(intersect(resultg[[85]],Tset),intersect(resultg[[(i-1)*20+5]],Tset))
    }
    
    

}

## figure 3
leaveone_recurrent <- function(){
    source("enrichana_6_12.R")
    options(stringsAsFactors=FALSE)
    TPcut=0.1
    options(stringsAsFactors=FALSE)
    netstr <- c("STRING/","iRef/","HPRD/","corr/","coexp5/")
    netnum <- c(31,6,20,26,27)
    TADAFile <- "../TADA_DAWN/result/TADAdenovo_meta_dmis.csv"
    asdall <- read.csv(TADAFile)
    asdall[is.na(asdall[,"qvalue.dn"]),"qvalue.dn"] <- 1
    Tset <- asdall[asdall[,"qvalue.dn"]<=TPcut,1]
    
    load("BestP_8_19_4")
    BestP <- BestP[c(1,2,4,6,7),]
    ### leave one mutation ==========================================
    nLeo <- matrix(0,13,2)
    
    oneg <- c()
    ### TADA method
    cfile <- "../TADA_DAWN/result/TADAdenovo_control.csv"
    ctmp <- read.csv(cfile)
    for(i in 1:length(Tset)){
        filename <- paste("../TADA_DAWN/result/leaveone4/TADAdenovo_rand2_",i,".csv",sep="")
        cutoff <- total_genes_cutoff(filename,Tset,muts="",alpha=0.5)
        tmp <- read.csv(filename)
        tmp[is.na(tmp[,"qvalue.dn"]),"qvalue.dn"] <- 1
        oner <- tmp[match(Tset[i],tmp[,1]),"qvalue.dn"]  <= 1-cutoff
        nLeo[1,1] <- nLeo[1,1] + oner
        #if(oner & (ctmp[match(Tset[i],ctmp[,1]),2] <= 1-cutoff)) nLeo[1,2] <- nLeo[1,2] + 1
        if( ctmp[match(Tset[i],ctmp[,1]),"qvalue.dn"] <= 1-cutoff ) nLeo[1,2] <- nLeo[1,2] + 1
    }
    
    ## DAWN method 
    ### control data result  
    cfile <- "../TADA_DAWN/DAWN_package/DAWNcontrol.csv" ###  checked
    ctmp <- read.csv(cfile)
    ctmp[is.na(ctmp[,"Stage2_posterior"]),"Stage2_posterior"] <- 1
    
    for(i in 1:length(Tset)){
        filename <- paste("../TADA_DAWN/DAWN_package/result/leaveone4/DAWN_leaveone4_",i,".csv",sep="")
        cutoff <- total_genes_cutoff(filename,Tset,muts="",alpha=0.5)
        tmp <- read.csv(filename)
        tmp[is.na(tmp[,"Stage2_posterior"]),"Stage2_posterior"] <- 1
        oner <- tmp[match(Tset[i],tmp[,1]),"Stage2_posterior"]  <= 1-cutoff
        
        if(!is.na(oner)){
            nLeo[2,1] <- nLeo[2,1] + oner
            if(oner) oneg <- union(oneg,Tset[i])
        }

        if( Tset[i] %in% ctmp[,1] ){
            if(ctmp[match(Tset[i],ctmp[,1]),"Stage2_posterior"] <= 1-cutoff) nLeo[2,2] <- nLeo[2,2] + 1
        }
        
    }
    DAWNg <- oneg
    
    oneg <- c()
    ###  MAGI method
    cfile <- "../MAGI_V0.1_Alpha/mydata/control_0/RandomGeneList.1"
    ctmp <- read.table(cfile)
    
    for(i in 1:length(Tset)){
        filename <- paste("../MAGI_V0.1_Alpha/mydata/leaveone4/RandomGeneList.",i,sep="")
        cutoff <- total_genes_cutoff(filename,Tset,muts="",alpha=0.5)
        tmp <- read.table(filename)
        oner <- tmp[match(Tset[i],tmp[,1]),2]  >= cutoff

        if(!is.na(oner)){
            nLeo[3,1] <- nLeo[3,1] + oner
            if(oner) oneg <- union(oneg,Tset[i])
        }
        if(Tset[i] %in% ctmp[,1]){
        if( ctmp[match(Tset[i],ctmp[,1]),2] >= cutoff ) nLeo[3,2] <- nLeo[3,2] + 1;}
    }
    MAGIg <- oneg
    
    oneg <- c()
    ### hotnet2 method 
    for(j in 1:5){
        cfile <- paste("result/control/hotnet/",netstr[j],"/hotnetresult1control",netnum[j],".txt",sep="")
        ctmp <- read.table(cfile)
        
        for(i in 1:length(Tset)){
            filename <- paste("result/leaveone4result/",netstr[j],"hotnetresult1rand2_",i,netnum[j],".txt",sep="")
            cutoff <- total_genes_cutoff(filename,Tset,muts="",alpha=0.5)
            tmp <- read.table(filename)
            oner <- tmp[match(Tset[i],tmp[,1]),2]  >= cutoff

            if(!is.na(oner)){
                if(oner) oneg <- union(oneg,Tset[i])
                nLeo[j+3,1] <- nLeo[j+3,1] + oner
            }
            if(Tset[i] %in% ctmp[,1]){
            if( ctmp[match(Tset[i],ctmp[,1]),2] >= cutoff ) nLeo[j+3,2] <- nLeo[j+3,2] + 1;}
        }       
    }
    RWRg <- oneg
    
    oneg <- c()
    usoneg <- c()
    ### our method
    for(j in 1:5){
        cfile <- paste("result/control/v",BestP[j,2],"/",netstr[j],"CRFresult_",BestP[j,3],"controlLBP_",netnum[j],".txt",sep="")
        ctmp <- read.table(cfile)
        for(i in 1:length(Tset)){
            filename <- paste("result/leaveone4result_",BestP[j,2],"/",netstr[j],"CRFresult_",BestP[j,3],"rand2_",i,"LBP_",netnum[j],".txt",sep="")
            cutoff <- total_genes_cutoff(filename,Tset,muts="",alpha=0.5)
            tmp <- read.table(filename)
            oner <- tmp[match(Tset[i],tmp[,1]),2]  >= cutoff
            if(!is.na(oner)){
                if(oner) {oneg <- union(oneg,Tset[i]); usoneg <- union(usoneg,Tset[i]);}
                nLeo[j+8,1] <- nLeo[j+8,1] + oner
            }
            if(Tset[i] %in% ctmp[,1]){
            if( ctmp[match(Tset[i],ctmp[,1]),2] >= cutoff ) nLeo[j+8,2] <- nLeo[j+8,2] + 1;}
        }       
    }
    MICRFg <- oneg
    save(nLeo,file="nLeo_8_16")
    
    ### plot figures
    
    pdf(file="plot/leaveone_8_19_4.pdf",width=10,height=10)
    par(mai=c(4,2,2,2))
    legend <- c("TADA","DAWN","MAGI","RWR-STRING","RWR-iRef","RWR-HPRD","RWR-corr","RWR-coexp","MICRF-STRING","MICRF-iRef","MICRF-HPRD","MICRF-corr","MICRF-coexp")
    ##ords <- c(1,2,3,4,9,5,10,6,11,7,12,8,13)
    ords <- c(1,2,3,4,9,5,10,6,11,8,13)
    #cols <- c("red","green","blue","orange","orange","yellow","yellow","pink","pink","cyan","cyan","purple","purple")
    
    barplot(t(nLeo)[,ords], names.arg=legend[ords], beside=T, ylab="#gold risk genes",cex.names=1.5, las=2, ylim=c(0,max(nLeo)), col=c("red","darkblue"),cex.lab=1.5,cex.axis=1.25)
    legend("topleft",legend=c("ASD case","SSC control"),cex=1.3,fill=c("red","darkblue")) 
    dev.off()    
    
}

meta_poisson_hq <- function(){
    
    source("Poisson_test_hq.R")
    filename <- "result/Nat1377_meta/TADAdenovo_Meta_dmis.csv"
    N <- 5542
    tmp <- read.csv(filename)
    n <- dim(tmp)[1]
    pV <- sapply(1:n, function(i){
        a1 <- poisson.test(x=tmp[i,"dn.LoF"], T = N *2* tmp[i,"LOF"], alternative = "greater", conf.level = 0.95)$p.value
        a2 <- poisson.test(x=tmp[i,"dn.mis3"], T = N *2*tmp[i,"dmis"], alternative = "greater", conf.level = 0.95)$p.value
        a3 <- poisson.test(x=tmp[i,"dn.LoF"]+ tmp[i,"dn.mis3"]+tmp[i,"dn.mis"], T = N *2* (tmp[i,"LOF"]+tmp[i,"mis"]) , alternative = "greater", conf.level = 0.95)$p.value
        c(a1,a2,a3)
    })
    pV <- t(pV)
    tmp[,"min_p"] <- apply(pV,1,min) 
    tmp[,"adjust_p"] <- p.adjust(tmp[,"min_p"],n=3*length(tmp[,"min_p"]),method="fdr") ##* 3
    #tmp[,"adjust_p"] <- p.adjust(tmp[,"min_p"],method="fdr") ##******* 3
    
    write.csv(tmp,file="result/Nat1377_meta/Ptest_meta_dmis.csv",row.names=FALSE)
    
}

meta_nbinom_hq_test <- function(){
    
    source("Poisson_test_hq.R")
    filename <- "result/Nat1377_meta/TADAdenovo_Meta_dmis.csv"
    N <- 5542
    tmp <- read.csv(filename)
    
    a <- rowSums(tmp[,c("dn.mis3","dn.LoF","dn.mis")])
    a <- a[a>0]
    b <- rowSums(tmp[,c("dn.mis3","dn.LoF")])
    b <- b[b>0]
    
    library(MASS)
    pas1 <- fitdistr(tmp[tmp[,"dn.LoF"]>0,"dn.LoF"],densfun="negative binomial")
    pas2 <- fitdistr(tmp[tmp[,"dn.mis3"]>0,"dn.mis3"],densfun="negative binomial")
    pas3 <- fitdistr(a,densfun="negative binomial",lower=c(0.01,0.01))
    pas4 <- fitdistr(b,densfun="negative binomial",lower=c(0.01,0.01))

    a1 <- dnbinom(x=tmp[,"dn.LoF"], size=pas1$estimate[1],mu=pas1$estimate[2])/dpois(x=tmp[,"dn.LoF"], 2*N* tmp[,"LOF"])
    a2 <- dnbinom(x=tmp[,"dn.mis3"], size=pas2$estimate[1],mu=pas2$estimate[2])/dpois(x=tmp[,"dn.mis3"], 2*N* tmp[,"dmis"])
    a3 <- dnbinom(x=rowSums(tmp[,c("dn.mis3","dn.LoF","dn.mis")]), size=pas3$estimate[1],mu=pas3$estimate[2])/dpois(x=rowSums(tmp[,c("dn.mis3","dn.LoF","dn.mis")]), 2*N*rowSums(tmp[,c("LOF","mis")]))
    a4 <- dnbinom(x=rowSums(tmp[,c("dn.mis3","dn.LoF")]), size=pas4$estimate[1],mu=pas4$estimate[2])/dpois(x=rowSums(tmp[,c("dn.mis3","dn.LoF")]), 2*N*tmp[,"LOFDmis"])
    
#     b1 <- pnbinom(apply(cbind(tmp[,"dn.LoF"]-1,0),1,max), size=pas1$estimate[1],mu=pas1$estimate[2],lower.tail=FALSE)
#     b2 <- pnbinom(apply(cbind(tmp[,"dn.mis3"]-1,0),1,max), size=pas2$estimate[1],mu=pas2$estimate[2],lower.tail=FALSE)
#     b3 <- pnbinom(apply(cbind(rowSums(tmp[,c("dn.mis3","dn.LoF","dn.mis")])-1,0),1,max), size=pas3$estimate[1],mu=pas3$estimate[2],lower.tail=FALSE)
#     
#     pV1 <- apply(cbind(b1,b2,b3),1,min)
#     pd <- p.adjust(pV1,n=3*length(pV1),method="BH")
#     sum(pd<0.1)
    
    pV <- cbind(a1,a2,a3,a4)
    tmp <- cbind(tmp,pV)
    
    source("Multi_net.R")
    # FDR estimation
    fdr1  <- Bayesian.FDR(pV[,1])$FDR
    fdr2  <- Bayesian.FDR(pV[,2])$FDR
    fdr3  <- Bayesian.FDR(pV[,3])$FDR
    fdr4  <- Bayesian.FDR(pV[,4])$FDR
    
    tmp <- cbind(tmp,fdr1,fdr2,fdr3,fdr4)
    write.csv(tmp,file="result/Nat1377_meta/Ptest_meta_BFs.csv",row.names=FALSE)
    
}

meta_nbinom_hq <- function(){
    
    source("Poisson_test_hq.R")
    filename <- "result/Nat1377_meta/TADAdenovo_Meta_dmis.csv"
    N <- 5542
    tmp <- read.csv(filename)
    tmp <- tmp[,c(1:10,15:17)]
    
    b <- rowSums(tmp[,c("dn.mis3","dn.LoF")])
    b <- b[b>0]
    
    library(MASS)
    pas4 <- fitdistr(b,densfun="negative binomial",lower=c(0.01,0.01))
    BFp <- dnbinom(x=rowSums(tmp[,c("dn.mis3","dn.LoF")]), size=pas4$estimate[1],mu=pas4$estimate[2])/dpois(x=rowSums(tmp[,c("dn.mis3","dn.LoF")]), 2*N*tmp[,"LOFDmis"])
    tmp <- cbind(tmp,BFp)
    
    source("Multi_net.R")
    # FDR estimation
    FDR  <- Bayesian.FDR(BFp)$FDR
    tmp <- cbind(tmp,FDR)
    write.csv(tmp,file="result/Nat1377_meta/Ptest_meta_BFp.csv",row.names=FALSE)
    
}

# Meta-data
Meta <- function(){
    
    source("enrichana_6_12.R")
    source("addresult.R")
    source("~/.Rprofile")
    
    library(VennDiagram)
    filenames <- c("result/Nat1377_meta/TADAdenovo_Meta_dmis.csv","result/Nat1377_meta/DAWN/DAWN_Meta_dmis.csv","result/Nat1377_meta/MAGI/RandomGeneList.3","result/Nat1377_meta/RWR/coexp/hotnetresult1meta_dmis7.txt","result/Nat1377_meta/MICRF/v4/coexp/CRFresult_0.2meta_dmisLBP_7.txt")
    
    geneT <- list()
    for(i in 1:length(filenames)){
        filename <- filenames[i]
        if(grepl(".csv",filename)){ result <- read.csv(filename); 
        }else{result <- read.table(filename);}
        
        if(i==1){
            geneT[[i]] <- result[rank(result[,"qvalue.dn"])<=100,1]
        }else if(i==2){
            geneT[[i]] <- result[rank(result[,"FDR"]) <= 100,1]
        }else{
            geneT[[i]] <- result[rank(-result[,2]) <= 100,1]
        }
    }
 
   # venn.diagram(x=list(TADA=geneT[[1]],DAWN=geneT[[2]],MAGI=geneT[[3]],RWR=geneT[[4]],MICRF=geneT[[5]]),file="test.tiff",margin=c(0.2,0.2,0.2,0.2))
    
    methods <- c("TADA","DAWN","MAGI","RWR","MICRF")
    allg <- c(geneT[[1]],geneT[[2]],geneT[[3]],geneT[[4]],geneT[[5]])
    uniT <- list()
    for(i in 1:5){
        uniT[[i]] <- setdiff(geneT[[i]],allg[-(((i-1)*100+1):(i*100))])
       # qwt(uniT[[i]],file=paste("plot/uniqueg",methods[i],".txt",sep=""))
    }
    
   dis <- read.delim("../human_disease_textmining_full.tsv",header=FALSE)
   for(i in 1:5){
        tmp <- dis[dis[,2]%in% uniT[[i]],]
        write.csv(tmp,file=paste("diseaseA",methods[i],".csv",sep=""),row.names=FALSE)
   }
   
    poifile <- "result/Nat1377_meta/Ptest_meta_dmis.csv"
    gp <- read.csv(poifile)  
    gpoi <- gp[gp[,"adjust_p"] < 0.05,1]
    
    Ltest <- test_genes()
    testgs <- Ltest$testgenes
    testgs[[12]] <- gpoi
    
    NM <- matrix(0,length(testgs),length(filenames))
    for(i in 1:length(testgs)){
        for(j in 1:length(filenames)){
            NM[i,j] <- length(intersect(testgs[[i]],geneT[[j]]))
        }
    }

    ##-------------------------------------------------------------------
    ## show in Cytoscape
    
    qwt(geneT[[5]],file="result/Nat1377_meta/MICRF_meta_dmis.txt")
    
    
    input <- read.csv("result/Nat1377_meta/TADAdenovo_Meta_dmis.csv")
    tmp <- cbind(input[,1],"=",-log10(input[,"pval.TADA.dn"]))
    
    net  <- read.table("data/network_inference/ASDexp_net_top5_single.txt")
    netg <- union(net[,1],net[,2])
    tg <- setdiff(netg,tmp[,1])
    tg <- cbind(tg," = ",0)
    tmp <- rbind(tmp,tg)
    qwt(tmp,file="plot/inputmeta.txt")
    
    genecol <- cbind(tmp[,1]," = ", 0)
    genecol[genecol[,1] %in% geneT[[5]] ,3] <- 1
    
    bnet <- net[net[,1] %in% geneT[[5]] | net[,2] %in% geneT[[5]],]
    coreg <- union(bnet[,1],bnet[,2])
    cnet <- net[,1] %in% coreg & net[,2] %in% coreg
    
    cnete <- net[cnet,]
    qwt(cnete,file="plot/subnetworks.txt")
    
    edgea <- cbind(net[,1], "(pp)", net[,2], "=", 0)
    edgea[cnet,5] <- 1
    qwt(edgea,file="plot/edges_meta1.txt",sep=" ")
    
    dg <- setdiff(coreg,geneT[[5]])
    genecol[genecol[,1] %in% dg,3] <- 2
    qwt(genecol ,file="plot/genecol_Cyto.txt")
    
    
    geneid <- cbind(tmp[,1]," = ", "")
    geneid[match(coreg,geneid[,1]),3] <- coreg
    qwt(geneid ,file="plot/geneid_Cyto.txt")

    
    showCyto(geneT)
}

# repeat DDD analysis
repDDDf <- function(){
    
    filenames <- c("../TADA_DAWN/result/TADAdenovo_repDDD_8_24.csv","../TADA_DAWN/DAWN_package/result/DAWN_repDDD.csv","../MAGI_V0.1_Alpha/mydata/control/RandomGeneList.2","result/control/hotnet/coexp5/hotnetresult1repDDD27.txt","result/control/v4/coexp5/CRFresult_0.6repDDDLBP_27.txt")
    
    geneT <- list()
    n <- 246
    for(i in 1:length(filenames)){
        filename <- filenames[i]
        if(grepl(".csv",filename)){ result <- read.csv(filename); 
        }else{result <- read.table(filename);}
        
        if(i==1){
            geneT[[i]] <- result[rank(result[,"qvalue.dn"])<= n,1]
        }else if(i==2){
            geneT[[i]] <- result[rank(result[,"FDR"]) <= n,1]
        }else{
            geneT[[i]] <- result[rank(-result[,2]) <= n,1]
        }
    }
    
    ##--------------------------------------------------------------------
    ## meta recurrece genes
    metad <- read.csv("../TADA_DAWN/result/TADAdenovo_meta_dmis.csv")
    metad[is.na(metad[,"qvalue.dn"]),"qvalue.dn"] <- 1
    
    genes <- metad[metad[,"qvalue.dn"] <= 0.1,1]
    genes1 <- metad[rowSums(metad[,c("dn.LoF","dn.mis3")]) > 0,1]
    
    nm1 <- matrix(0,2,5)
    for(i in 1:length(geneT)){
        nm1[1,i] <- length(intersect(geneT[[i]],genes))
        nm1[2,i] <- length(intersect(geneT[[i]],genes1))
    }

    nm1
}

## show in Cytoscape with iRefIndex PPI
showCyto <- function(geneT){

    netfile <- "data/hotnet/iRefIndexm.txt"
    str <- "iRef"
    
    ##-------------------------------------------------------------------
    ## show in Cytoscape
    input <- read.csv("result/Nat1377_meta/TADAdenovo_Meta_dmis.csv")
    tmp <- cbind(input[,1],"=",-log10(input[,"pval.TADA.dn"]))
    
    net  <- read.table(netfile)
    netg <- union(net[,1],net[,2])
    tg <- setdiff(netg,tmp[,1])
    tg <- cbind(tg," = ",0)
    tmp <- rbind(tmp,tg)
    qwt(tmp,file=paste("plot/inputmeta",str,".txt",sep=""))
    
    genecol <- cbind(tmp[,1]," = ", 0)
    genecol[genecol[,1] %in% geneT[[5]] ,3] <- 1
    
    bnet <- net[net[,1] %in% geneT[[5]] | net[,2] %in% geneT[[5]],]
    coreg <- union(bnet[,1],bnet[,2])
    cnet <- net[,1] %in% coreg & net[,2] %in% coreg
    
    cnete <- net[cnet,]
    qwt(cnete,file=paste("plot/subnetworks",str,".txt",sep=""))
    
    edgea <- cbind(net[,1], "(pp)", net[,2], "=", 0)
    edgea[cnet,5] <- 1
    qwt(edgea,file=paste("plot/edges_meta1",str,".txt",sep=""),sep=" ")
    
    dg <- setdiff(coreg,geneT[[5]])
    genecol[genecol[,1] %in% dg,3] <- 2
    qwt(genecol ,file=paste("plot/genecol_Cyto",str,".txt",sep=""))
    
    
    geneid <- cbind(tmp[,1]," = ", "")
    geneid[match(coreg,geneid[,1]),3] <- coreg
    qwt(geneid ,file=paste("plot/geneid_Cyto",str,".txt",sep=""))

}

## module function for co-expression networks
modulef <- function(){
    load("DAWNmodule")
    
    for(i in 1:max(module[,2])){
        genes <- module[module[,2]==i,1]
        qwt(genes,file=paste("moduleCo/module",i,".txt",sep=""))
    }
}

test_genes <- function(){
    
    source("Network_analysis.R")
    Rbfoxt <- read.table("ASD/Targets/RBFOXt/RBFOX_target_genes.txt")
    Rbfoxt <- mapping_to(Rbfoxt) # 597
    
    FMRPt1 <- read.csv("ASD/Targets/FMRP2/mmc2.csv",skip=1)[,"Gene.Symbol"]
    FMRPt1 <- mapping_to(FMRPt1) # 842
    FMRPt1 <- unique(FMRPt1)
    
    FMRPt2 <- read.csv("ASD/Targets/FMRP2/mmc5.csv",skip=1)[,"Gene.Symbol"]
    FMRPt2 <- mapping_to(FMRPt2)  # 67
    FMRPt2 <- unique(FMRPt2)
    
    Sfarigene <- read.csv("ASD/Targets/SFARIgene/gene-S-H-Strong.csv")[,1]
    Sfarigene <- mapping_to(Sfarigene) # 71
    Sfarigene <- unique(Sfarigene)
    
    cancer1 <- read.csv("ASD/Targets/cancergenes/cancer_gene_census.csv")[,"Gene.Symbol"]
    cancer1 <- mapping_to(cancer1)
    cacer1 <- unique(cancer1)
    
    cancer2 <- read.csv("ASD/Targets/cancergenes/1235122TablesS2A.csv")[,"Gene.Symbol"]
    cancer2 <- mapping_to(cancer2)
    cancer2 <- unique(cancer2)
    
    CHDgenes <- read.table("ASD/Targets/CHDnonslientgenes.txt")
    CHDgenes <- mapping_to(CHDgenes)
    CHDgenes <- unique(CHDgenes)
    
    IDgenes <- read.table("ASD/Targets/IDgenes.txt")
    IDgenes <- mapping_to(IDgenes)
    IDgenes <- unique(IDgenes)
    
    DDDgenes <- read.table("ASD/Targets/DDDgenes.txt")
    DDDgenes <- mapping_to(DDDgenes)
    DDDgenes <- unique(DDDgenes)
    
    HMGgene <- read.csv("data/genedata/HMGs.csv",header=FALSE)[,2] # there is an unkonwn in 152 HMGgene
    TGFgene <- read.delim("data/genedata/TGF-beta.txt",header=FALSE,sep="\t")[,2]
    HMGgene <- mapping_to(HMGgene)
    HMGgene <- HMGgene[HMGgene!=""]
    HMGgene <- unique(HMGgene)
    TGFgene <- mapping_to(TGFgene)
    TGFgene <- unique(TGFgene)
    
    testgenes <- list()
    testgenes[[1]] <- Rbfoxt
    testgenes[[2]] <- FMRPt1
    testgenes[[3]] <- Sfarigene
    testgenes[[4]] <- cancer1
    testgenes[[5]] <- cancer2
    testgenes[[6]] <- CHDgenes
    testgenes[[7]] <- IDgenes
    testgenes[[8]] <- DDDgenes
    testgenes[[9]] <- FMRPt2
    testgenes[[10]] <- HMGgene
    testgenes[[11]] <- TGFgene
    
    Lsets <- c("RBFOX target","FMRP 842 target","SFARI autism genes","COMIC cancer census","Cancer driver genes","CHD genes","ID genes","DDD genes","FMRP autism target","HMG genes","TGF beta genes")
    
    list(testgenes=testgenes,Lsets=Lsets)
}

## two gene lists for Yufeng
twolists <- function(){
    
    filenames <- c("ASD/TADAresult/TADAdenovo_Meta.csv","result/Nat1377_meta/DAWN/DAWN_Meta.csv","result/Nat1377_meta/MAGI/RandomGeneList.3","result/Nat1377_meta/RWR/coexp/hotnetresult1Meta_7_37.txt","result/Nat1377_meta/MICRF/v4/coexp/CRFresult_0.2Meta_7_3LBP_7.txt")
 
        #filenames <- c("result/Nat1377_meta/TADAdenovo_nat13772.csv","result/Nat1377_meta/DAWN/DAWN_nat13772.csv","result/Nat1377_meta/MAGI/RandomGeneList.4","result/Nat1377_meta/RWR/coexp/hotnetresult1nat137727.txt","result/Nat1377_meta/MICRF/CRFresult_0.2nat13772LBP_7.txt")
        MICRFr <- read.table(filenames[5])
        MICRFg <- MICRFr[MICRFr[,2] > 0.5,1]
        n <- length(MICRFg)
        
    TADAr <- read.csv(filenames[1])
    TADAg <- TADAr[TADAr[,"qvalue.dn"] < 0.5,1]
    
    DAWNr <- read.csv(filenames[2])
    DAWNr <- DAWNr[!is.na(DAWNr[,"FDR"]),]
    DAWNg <- DAWNr[DAWNr[,"FDR"] < 0.5,1]
    
    MAGIr <- read.table(filenames[3])
    MAGIr <- MAGIr[order(- MAGIr[,2]),]
    MAGIg <- MAGIr[1:n,1]
    
    RWRr <- read.table(filenames[4])
    RWRr <- RWRr[order(-RWRr[,2]),]
    RWRg <- RWRr[1:n,1]
    
 
    uni <- setdiff(MICRFg, union(TADAg,union(DAWNg,union(MAGIg,RWRg))))
    qwt(uni,file="uniMICRF_Meta.txt")
  
    uni1 <- setdiff(TADAg, union(MICRFg,union(DAWNg,union(MAGIg,RWRg))))
    qwt(uni1,file="uniTADA_Meta.txt")
}

## complement networks
comple <- function(){

    source("enrichana_6_12.R")
    
    filenames <- c("result/Nat1377_meta/TADAdenovo_Meta_dmis.csv","result/Nat1377_meta/DAWN/DAWN_Meta_dmis.csv","result/Nat1377_meta/MAGI/RandomGeneList.3","result/Nat1377_meta/RWR/coexp1/hotnetresult1meta_dmis21.txt","result/Nat1377_meta/RWR/HPRD/hotnetresult1meta_dmis20.txt","result/Nat1377_meta/RWR/STRING/hotnetresult1meta_dmis31.txt","result/Nat1377_meta/RWR/iRef/hotnetresult1meta_dmis6.txt","result/Nat1377_meta/RWR/coexp/hotnetresult1meta_dmis7.txt","result/Nat1377_meta/MICRF/v4/coexp/CRFresult_0.2meta_dmisLBP_7.txt") 
    
    ## meta data
    
    # coexp: 204  STRING:231 HPRD:192  IREF: 227
    b <- read.csv(filenames[1])
    b <- b[!duplicated(b[,1]),]
    noncov <- b[b[,"qvalue.dn"]<0.3,1] # 0.1:  178; 0.3 400
    
    for(i in 2:9){
        filename <- filenames[i]
        if(grepl(".csv",filename)){ result <- read.csv(filename); 
        }else{result <- read.table(filename);}
        
        noncov <- setdiff(noncov,result[,1])
        print(length(noncov))
    }
    

    source("~/.Rprofile")
    qwt(noncov,file="TADA_meta_3.txt")
    
    net <- read.csv("/Users/qh2159/Desktop/Papers/Network analysis method/data/CellNet/Human_Big_GRN_032014.csv")
    netg <- union(net[,1],net[,2])
    
    sum(noncov %in% netg) ## 86
    
    Rbfoxt <- read.table("ASD/Targets/RBFOXt/RBFOX_target_genes.txt")
    Rbfoxt <- mapping_to(Rbfoxt) # 597
    
    FMRPt1 <- read.csv("ASD/Targets/FMRP2/mmc2.csv",skip=1)[,"Gene.Symbol"]
    FMRPt1 <- mapping_to(FMRPt1) # 842
    FMRPt1 <- unique(FMRPt1)
    
    sum(noncov %in% Rbfoxt) # 3
    sum(noncov %in% FMRPt1) # 5
    
    fbrain <- read.table("../cell6312mmc4/fBrain/genes-regulate-genes.txt")
    netg1 <- union(fbrain[,1],fbrain[,2])
    sum(noncov %in% netg1) # 0
    
   
    RBPs <- read.csv("../splicingnet/023010-2.csv")[,1]
    sum(noncov %in% RBPs)
     
}

## supplement AUC for all networks and methods 
AUC_supple <- function(){
    nmeth=13
    pdf(file="plot/samplesize_allS.pdf",width=12,height=10)    
    
    load("AUC_all_7_29")
    cols <- c("black","darkred","red","orange","green","blue","blue","orange","yellow4","violet","darkorange","coral","violet")
    legend <- c("TADA","DAWN","MAGI","RWR-STRING","RWR-iRef","RWR-coexp","RWR-HPRD","RWR-corr","MICRF-STRING","MICRF-iRef","MICRF-coexp","MICRF-HPRD","MICRF-corr")
    
    main="AUC performance for different trios"
    xlab="The number of trios";ylab="AUC";
    ymint <- 1;ymaxt <- 0;
    ymint <- min(AUC);ymaxt <- max(AUC);
    ymin=max(min(ymint)+0.1,0);ymax = min(ymaxt+0.1,1);
    x <- seq(0.2,0.9,0.1)
    Y1 <- rep(0,nmeth)
    subs <- 1:nmeth
    nn <- subs
    k=1
    for(i in subs){
        y <- rowMeans(AUC[,,i])
        sd <- apply( AUC[,,i], 1, sd)/2 
        plot_err(x[1:5],y[1:5],sd,i,main,xlab,ylab,ylim=c(0,ymax),xlim=c(0.15,0.62),cols=cols)
        Y1[k] <- y[1]
        k <- k+1
    } 
    Y2 <- round(Y1*100)/100
    labels <- c("1,","2,","3","4","5,","6,","7,","8,","9","10","11","12","13")
    tx <- c(0.18,0.192,0.2,0.2,0.162,0.172,0.182,0.192,0.2,0.2,0.2,0.2,0.2)
    ty <- -sort(-Y2)
    ty[5] <- ty[4]
    ty[1:3] <- mean(ty[1:3]) 
    ty[10:12] <- ty[10:12] - 0.005
    ty[5:9] <- mean(ty[5:9])
    ty[4] <- ty[4] + 0.005
    ords <- sort(-Y1,index.return=T)$ix
    text(tx,ty,labels=labels,pos=2,col=cols[ords],cex=1.3)
    
    #ords <- c(1,2,3,4,7,5,8,6,9,12,10,13,11)
    legend("bottomright",col=cols[ords],legend=paste(subs,legend[ords],sep=":"),lty=rep(1,length(subs)),lwd=rep(2,length(subs)),cex=1.5,y.intersp=1)
    
    dev.off()

}

### not cover genes by network methods 
outlier_gf <- function(){
    Tadaf <- "../TADA_DAWN/result/TADAdenovo_meta_dmis.csv"
    tmp <- read.csv(Tadaf)
    tmp[is.na(tmp[,"qvalue.dn"]),"qvalue.dn"] <- 1
    tadag <-  tmp[tmp[,"qvalue.dn"]<=0.3,1] ## FDR <= 0.3
    
    n <- length(tadag)
    dawng <- read.csv("../TADA_DAWN/DAWN_package/result/DAWN_meta_dmis.csv")[1:n,1]
    magir <- read.table("../MAGI_V0.1_Alpha/mydata/control_0/RandomGeneList.3")
    magir <- magir[order( - as.numeric(magir[,2])), ]
    magig <- magir[1:n,1]
    rwrg <- read.table("result/control/hotnet/coexp5/hotnetresult1meta27.txt")[1:n,1]
    rwrg2 <- read.table("result/control/hotnet/STRING/hotnetresult1meta31.txt")[1:n,1]
    rwrg3 <- read.table("result/control/hotnet/iRef/hotnetresult1meta6.txt")[1:n,1]
    rwrg4 <- read.table("result/control/hotnet/HPRD/hotnetresult1meta20.txt")[1:n,1]
    
    micrfg <- read.table("result/control/v4/coexp5/CRFresult_0.6metaLBP_27.txt")[1:n,1]
    micrfg2 <- read.table("result/control/v4/STRING/CRFresult_0.4metaLBP_31.txt")[1:n,1]
    micrfg3 <- read.table("result/control/v4/iRef/CRFresult_0.6metaLBP_6.txt")[1:n,1]
    micrfg4 <- read.table("result/control/v4/HPRD/CRFresult_0.4metaLBP_20.txt")[1:n,1]
    
    
    a <- tadag
    a <- setdiff(a,dawng)
    a <- setdiff(a,magig)
    a <- setdiff(a,rwrg)
    a <- setdiff(a,rwrg2)
    a <- setdiff(a,rwrg3)
    a <- setdiff(a,rwrg4)
    a <- setdiff(a,micrfg)
    a <- setdiff(a,micrfg2)
    a <- setdiff(a,micrfg3)
    a <- setdiff(a,micrfg4)
    
    a
    write.table(a,file="TADAonly_0.3.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)
}

### pathways genes
pathways_genes <- function(){
    source("Network_analysis.R")
    HMGgene <- read.csv("data/genedata/HMGs.csv",header=FALSE)[,2] # there is an unkonwn in 152 HMGgene
    TGFgene <- read.delim("data/genedata/TGF-beta.txt",header=FALSE,sep="\t")[,2]
    #HMGgene <- mapping_to(HMGgene)
    HMGgene <- HMGgene[HMGgene!=""]
    HMGgene <- unique(HMGgene)
    #TGFgene <- mapping_to(TGFgene)
    TGFgene <- unique(TGFgene)
    length(TGFgene)
    length(HMGgene)
    length(intersect(TGFgene,micrfg))
    length(intersect(HMGgene,micrfg))
    
    filenames <- c("data/pathways/Hedgehog_pathway.txt","data/pathways/Wnt_pathway.txt","data/pathways/Notch_pathway.txt","data/pathways/MAPK_signaling_pathway.txt","data/pathways/Mtor_pathway.txt","data/pathways/Ras_pathway.txt")
    for(i in 1:length(filenames)){
        a <- read.delim(filenames[i],header=FALSE,sep="\t",comment.char="#")
        genes <- sapply(1:dim(a)[1], function(i) {tmp <- unlist(strsplit(a[i,1],";"))[1]; unlist(strsplit(tmp,","))[1];})
        print(length(genes))
        genes <- mapping_to(genes)
        print(length(intersect(genes,micrfg)))
    }

}
