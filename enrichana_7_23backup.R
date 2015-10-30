randset_1_AUC_6_3 <- function(){
    
    source("enrichana_6_12.R")
    
    load("BestP_6_10")
    betaV <- BestP[,3]
    filenames <- filenamesf(0,0,betaV=betaV,"","")
    allP <- allPf(TPcut=0.1)
    Tset <- allP
    
#     tmp <- gene_select(filenames,allP)
#     genes <- tmp$genes
#     allP <- tmp$allP
#     allN <- tmp$allN
#     tmp <- new_auc(filenames,allP,genes,flag=2,iplot=FALSE);
#     
#     tmp <- read.csv("result/randresult_1/DAWN/DAWN_randset1.csv")
#     geneDAWN <- intersect(genes,tmp[!is.na(tmp[,"FDR"]),"Gene"])
#     tmp <- new_auc(filenames,allP,geneDAWN,flag=2,iplot=FALSE);
    
    DAWNflag <- 1
    DAWNnames <-  readLines(con <- file("result/randresult_1/DAWN/TADAdenovo_randset_1.txt","r"))   
    close(con)
    Maginames <-  readLines(con <- file("result/randresult_1/MAGI/MAGIfilenames.txt","r"))   
    close(con)

    ### ROC example plot ================================================
    #pdf(file="plot/ROC_1.pdf",width=12,height=12)
    par(mai=c(1.2,1.2,1,1))
    j=2
    i=20
    subs =  c(1,10,11)
    filenames <- filenamesf(j,i,betaV=betaV,DAWNnames,Maginames)
    #tmp <- gene_select(filenames[subs],Tset)
    tmp <- gene_select(filenames,Tset)
    genes <- tmp$genes
    allP <- tmp$allP
    allN <- tmp$allN
   
    if(DAWNflag==1){
        tmp <- read.csv(filenames[2])
        genes <- intersect(genes,tmp[!is.na(tmp[,"FDR"]),"Gene"])
    }

    tmp <- ROC_ex(filenames[subs],allP,genes,flag=2)
    #tmp <- ROC_ex(filenames[subs],allP,genes,flag=2,xlim=c(0,0.1),ylim=c(0.3,0.6))
    #tmp <- ROC_ex(filenames[subs],allP,genes,flag=2,xlim=c(0.05,0.3),ylim=c(0.8,1))
    #dev.off()
    ### ==================================================================
    
    nnet <- 5
    AUC <- array(0,dim=c(8,20,3+nnet*2)) 
    FPR <- array(0,dim=c(8,20,3+nnet*2))
    TPR <- array(0,dim=c(8,20,3+nnet*2))

    for(j in 2:9){
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
            
            tmp <- new_auc(filenames,allP,genes,flag=2,iplot=FALSE,DAWNflag=DAWNflag)
            AUC[j-1,i,] <- tmp$auc
            
            for(kk in 1:length(filenames)){
                TPR[j-1,i,kk] <- tmp$TPR[[kk]][which.min(abs(tmp$FPR[[kk]]-0.2))]
                FPR[j-1,i,kk] <- tmp$TPR[[kk]][which.min(abs(tmp$FPR[[kk]]-0.1))]
            }
        }
    }
    
    save(AUC,file="AUC_m_6_11")
    save(TPR,file="TPR_m_6_11")
    save(FPR,file="FRP_m_6_11")
    
    ## plot auc values
    ## Figure 1: AUC values for eight random sets
    
    #nmeth <- length(filenames)

    nmeth=13
    pdf(file="plot/AUCs6_112.pdf",width=10,height=12)
    
    load("AUC_m_6_11")
    cols <- 1:13 #c("Black","Blue","Cyan","Brown","Coral","DarkViolet","Green","HotPink","Orange","Plum","deeppink")
    labels <- c("TADA","DAWN","MAGI","RWR-STRING","RWR-iRef","RWR-coexp","RWR-Infer","RWR-GNAT","MICRF-STRING","MICRF-iRef","MICRF-coexp","MICRF-Infer","MICRF-GNAT")
    legend <- c("TADA","DAWN","MAGI","RWR-STRING","RWR-iRef","RWR-coexp","RWR-Infer","RWR-GNAT","MICRF-STRING","MICRF-iRef","MICRF-coexp","MICRF-Infer","MICRF-GNAT")
    
    main="AUC performance for different trios"
    xlab="The number of trios";ylab="AUC";
    ymint <- 1;ymaxt <- 0;
    ymint <- min(AUC);ymaxt <- max(AUC);
    ymin=max(min(ymint)+0.1,0);ymax = min(ymaxt+0.1,1);
    x <- seq(0.2,0.9,0.1)
    
    Y1 <- rep(0,nmeth)
    
    nn <- nmeth  ##!!!!!
    
    for(i in 1:nn){
        y <- rowMeans(AUC[,,i])
        sd <- apply( AUC[,,i], 1, sd)/2 
        plot_err(x,y,sd,i,main,xlab,ylab,ylim=c(ymin,ymax),xlim=c(-0.05,0.95))
        Y1[i] <- y[1]
    }
    text(rep(0.2,nn),Y1,labels=labels[1:nn],pos=2,col=cols,cex=1.5)
    legend("bottomright",col=cols[1:nmeth],legend=legend[1:nn],lty=rep(1,nmeth),lwd=rep(2,nmeth),cex=1.2,y.intersp=0.8)
    
    dev.off()
    
    #### selected plot
    nmeth=7
    subs <- c(1,2,3,6,8,11,13)
    pdf(file="plot/AUCs6_112.pdf",width=10,height=12)
    
    load("AUC_m_6_11")
    cols <- 1:13 #c("Black","Blue","Cyan","Brown","Coral","DarkViolet","Green","HotPink","Orange","Plum","deeppink")
    labels <- c("TADA","DAWN","MAGI","RWR-STRING","RWR-iRef","RWR-coexp","RWR-Infer","RWR-GNAT","MICRF-STRING","MICRF-iRef","MICRF-coexp","MICRF-Infer","MICRF-GNAT")
    legend <- c("TADA","DAWN","MAGI","RWR-STRING","RWR-iRef","RWR-coexp","RWR-Infer","RWR-GNAT","MICRF-STRING","MICRF-iRef","MICRF-coexp","MICRF-Infer","MICRF-GNAT")
    
    main="AUC performance for different trios"
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
        plot_err(x,y,sd,k,main,xlab,ylab,ylim=c(ymin,ymax),xlim=c(-0.05,0.95))
        Y1[k] <- y[1]
        k <- k+1
    }
    text(rep(0.2,nmeth),Y1,labels=labels[nn],pos=2,col=cols[1:nmeth],cex=1.5)
    legend("bottomright",col=cols[1:nmeth],legend=legend[nn],lty=rep(1,nmeth),lwd=rep(2,nmeth),cex=1.2,y.intersp=0.8)
    
    dev.off()
    
    
}

batch_randset4_6_2 <- function(){
    source("enrichana_6_12.R")
    betaV <- c(0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0)
    nmeth <- 4
    
    TPcut=0.1
    TADAFile <- "ASD/TADAresult/TADAdenovo_ASD4_16.csv"
    asdall <- read.csv(TADAFile)
    Tset <- asdall[asdall[,"qvalue.dn"]<=TPcut,1]
    
    netstr <- c("STRING/","iRef/","Infer/","coexp/","GNAT/")
    netflags <- c(31,6,8,7,21)
    nodeN <- c(12969,12124,9157,9448,9998)
    
    PerM <- array(0,c(5,4,12))
    ReVL <- array(0,c(5,4,12))
    for(netflag in netflags){
        for(meth in 1:nmeth){
            for(beta in betaV){
                
                netsub <- which(netflags==netflag)
                ## recurrent
                path <- paste("result/randresult4_",meth,"/",netstr[netsub],sep="")
                ReM <- randset4_one(netflag,meth,beta,path)
                ReV <- colMeans(ReM)
                
                ## leave-one rank
                path <- paste("result/leaveone4result_",meth,"/",netstr[netsub],sep="")
                HiV <- matrix(0,1,nodeN[netsub])
                HiV <- leaveone4_one(netflag,meth,beta,Tset,path,HiV)
                
                if(FALSE){
                    pdf(file="plot/ParaSel.pdf",width=9,height=9)
                    plot(ReV/max(ReV),type="l",col=1,xlab="Gene list rank",ylab="Fraction of recurrent (leave-one)", main="Parameter Selection")
                    lines(HiV/max(HiV),type="l",col=2)
                    lines(ReV/max(ReV) * HiV/max(HiV) ,type="l",col=3)
                    a <- ReV/max(ReV) * HiV/max(HiV)
                    #abline(h=a[which.max(a)],lty=2) 
                    #abline(v=which.max(a),lty=2)
                    lines(which.max(a),a[which.max(a)],type="p",pch=19,col=3)
                    legend("topright",legend=c("Recurrent","Leave-one","Combined"),col=1:3,lwd=rep(1,3),lty=rep(1,3))
                    dev.off()
                }
                
                tmp <- ReV*HiV
                
                PerM[which(netflag==netflags),meth,which(betaV==beta)] <- max(tmp[10:length(tmp)])
                #ReVL[which(netflag==netflags),meth,which(betaV==beta)] <- max(ReV*HiV)
            }
        }
    } 
    save(PerM,file="PerM_6_10")
    
    BestP <- matrix(0,5,3)
    for(i in 1:5){
        BestP[i,1] <- netflags[i]
        tmp <- which(PerM[i,,]==max(PerM[i,,]),arr.ind=TRUE)
        BestP[i,2] <- tmp[1]
        BestP[i,3] <- betaV[tmp[2]]
    }
    BestP <- BestP[c(1,2,4,3,5),]
    save(BestP,file="BestP_6_10")
}

batch_randset4_7_22 <- function(){
    source("enrichana_6_12.R")
    betaV <- c(0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0)
    nmeth <- 4
    
    TPcut=0.1
    TADAFile <- "ASD/TADAresult/TADAdenovo_ASD4_16.csv"
    asdall <- read.csv(TADAFile)
    Tset <- asdall[asdall[,"qvalue.dn"]<=TPcut,1]
    
    netstr <- c("STRING/","iRef/","coexp/","HPRD/","coexp1/")
    netflags <- c(31,6,7,20,21)
    nodeN <- c(12969,12124,9448,10526,8782)
    
    PerM <- array(0,c(5,4,12))
    for(netflag in netflags){
        for(meth in 4){
            for(beta in betaV){
                
                netsub <- which(netflags==netflag)
                ## recurrent
                path <- paste("result/randresult4_",meth,"/",netstr[netsub],sep="")
                ReM <- randset4_one(netflag,meth,beta,path)
                ReV <- colMeans(ReM)
                
                ## leave-one rank
                path <- paste("result/leaveone4result_",meth,"/",netstr[netsub],sep="")
                HiV <- matrix(0,1,nodeN[netsub])
                HiV <- leaveone4_one(netflag,meth,beta,Tset,path,HiV)
                     
                tmp <- ReV*HiV
                
                PerM[which(netflag==netflags),meth,which(betaV==beta)] <- max(tmp[10:length(tmp)])
            }
        }
    } 
    
    save(PerM,file="PerM_7_23")
    
    BestP <- matrix(0,5,3)
    for(i in 1:5){
        BestP[i,1] <- netflags[i]
        tmp <- which(PerM[i,,]==max(PerM[i,,]),arr.ind=TRUE)
        BestP[i,2] <- tmp[1]
        BestP[i,3] <- betaV[tmp[2]]
    }
    save(BestP,file="BestP_7_23")
}

randset4_one <- function(netflag,meth,beta,path){
    
    ReM <- matrix(0,20,nodeN[netsub])
    for(i in 1:20){
        filename <- paste(path,"CRFresult_",beta,"part3_",i,"LBP_",netflag,".txt",sep="")
        ReM[i,] <- one_recurrent(filename,i)
    }
    
    ReM
    
}

one_recurrent <- function(filename,i){
    
    options(stringsAsFactors=FALSE)
    j=3
    onefile <- paste("ASD/TADAresult/randset4/TADAdenovo_ASD2sub",j,"derest",i,".csv",sep="")
    oneresult <- read.csv(onefile)
    
    TADAfile <- paste("ASD/TADAresult/randset4/TADAdenovo_ASD1sub",j,"depart",i,".csv",sep="")
    subresult <- read.csv(TADAfile)
    
    result <- read.table(filename)
    
    inrtmp <- rep(0, dim(result)[1])
    subs <- result[,1] %in% oneresult[,1]
    genes <- result[subs,1]
    inrtmp[subs] <- rowSums(oneresult[match(genes,oneresult[,1]),c("dn.LoF","dn.mis3")])>0
    
    tmp <- rep(0,length(inrtmp))
    for(m in 1:dim(result)[1]) tmp[m] <- sum(inrtmp[1:m]) / m 
    
    tmp
}

leaveone4_one <- function(netflag,meth,beta,Tset,path,HiV){
    
    n <- length(Tset)
    for(i in 1:n){
        filename <- paste(path,"CRFresult_",beta,"rand2_",i,"LBP_",netflag,".txt",sep="")
        result <- read.table(filename)
        subs <- match(Tset[i],result[,1])
        if(!is.na(subs)){HiV[subs] <- HiV[subs] + 1;}
    }
    
    tmp <- rep(0,length(HiV))
    for(m in 1:length(HiV)) tmp[m] <- sum(HiV[1:m])^2 / (m*n) 
    
    tmp
}

best_performance_6_11 <- function(){
    
    source("enrichana_6_12.R")
    options(stringsAsFactors=FALSE)
    TPcut=0.1
    options(stringsAsFactors=FALSE)
    netstr <- c("STRING/","iRef/","coexp/","Infer/","GNAT/")
    netnum <- c(31,6,7,8,21)
    TADAFile <- "ASD/TADAresult/TADAdenovo_ASD4_16.csv"
    asdall <- read.csv(TADAFile)
    Tset <- asdall[asdall[,"qvalue.dn"]<=TPcut,1]
    
    ### leave one mutation ==========================================
    
    nLeo <- matrix(0,13,2)
    
    oneg <- c()
    
    ### TADA method
    cfile <- "ASD/TADAresult/TADAdenovo_control.csv"
    ctmp <- read.csv(cfile)
    for(i in 1:length(Tset)){
        filename <- paste("ASD/TADAresult/leaveone4/TADAdenovo_rand2_",i,".csv",sep="")
        cutoff <- total_genes_cutoff(filename,Tset,muts="",alpha=0.5)
        tmp <- read.csv(filename)
        oner <- tmp[match(Tset[i],tmp[,1]),"qvalue.dn"]  <= 1-cutoff
        nLeo[1,1] <- nLeo[1,1] + oner
        if(oner & (ctmp[match(Tset[i],ctmp[,1]),2] <= 1-cutoff)) nLeo[1,2] <- nLeo[1,2] + 1
    }
    
    ## DAWN method 
    ### control data result  
    cfile <- "result/control_6_12/DAWN/DAWNcontrol.csv"
    ctmp <- read.csv(cfile)
    ctmp[is.na(ctmp[,"Stage2_posterior"]),"Stage2_posterior"] <- 1
    
    for(i in 1:length(Tset)){
        filename <- paste("result/leaveone4result/DAWN/DAWN_leaveone4_",i,".csv",sep="")
        cutoff <- total_genes_cutoff(filename,Tset,muts="",alpha=0.5)
        tmp <- read.csv(filename)
        oner <- tmp[match(Tset[i],tmp[,1]),"Stage2_posterior"]  <= 1-cutoff
        if(!is.na(oner)){
            nLeo[2,1] <- nLeo[2,1] + oner
            if(oner) oneg <- union(oneg,Tset[i])
            if(oner & (Tset[i] %in% ctmp[,1])){
                if(ctmp[match(Tset[i],ctmp[,1]),"Stage2_posterior"] <= 1-cutoff) nLeo[2,2] <- nLeo[2,2] + 1
            }
        }
    }
    
    ###  MAGI method
    cfile <- "result/control_6_12/MAGI/RandomGeneList.1"
    ctmp <- read.table(cfile)
    
    for(i in 1:length(Tset)){
        filename <- paste("result/leaveone4result/MAGI/RandomGeneList.",i,sep="")
        cutoff <- total_genes_cutoff(filename,Tset,muts="",alpha=0.5)
        tmp <- read.table(filename)
        oner <- tmp[match(Tset[i],tmp[,1]),2]  >= cutoff
        if(!is.na(oner)){
            nLeo[3,1] <- nLeo[3,1] + oner
            if(oner) oneg <- union(oneg,Tset[i])
            if(oner & (ctmp[match(Tset[i],ctmp[,1]),2] >= cutoff)) nLeo[3,2] <- nLeo[3,2] + 1
        }
    }
    
    ### hotnet2 method 
    for(j in 1:5){
        cfile <- paste("result/control_6_12/hotnet/hotnetresult1control",netnum[j],".txt",sep="")
        ctmp <- read.table(cfile)
        
        for(i in 1:length(Tset)){
            filename <- paste("result/leaveone4result/",netstr[j],"hotnetresult1rand2_",i,netnum[j],".txt",sep="")
            cutoff <- total_genes_cutoff(filename,Tset,muts="",alpha=0.5)
            tmp <- read.table(filename)
            oner <- tmp[match(Tset[i],tmp[,1]),2]  >= cutoff
            if(!is.na(oner)){
                if(oner) oneg <- union(oneg,Tset[i])
                nLeo[j+3,1] <- nLeo[j+3,1] + oner
                if(oner & (ctmp[match(Tset[i],ctmp[,1]),2] >= cutoff)) nLeo[j+3,2] <- nLeo[j+3,2] + 1
            }
        }       
    }
    
    usoneg <- c()
    ### our method
    load("BestP_6_10")    
    for(j in 1:5){
        cfile <- paste("result/control_6_12/v",BestP[j,2],"/",netstr[j],"CRFresult_",BestP[j,3],"controlLBP_",netnum[j],".txt",sep="")
        ctmp <- read.table(cfile)
        for(i in 1:length(Tset)){
            filename <- paste("result/leaveone4result_",BestP[j,2],"/",netstr[j],"CRFresult_2rand2_",i,"LBP_",netnum[j],".txt",sep="")
            cutoff <- total_genes_cutoff(filename,Tset,muts="",alpha=0.5)
            tmp <- read.table(filename)
            oner <- tmp[match(Tset[i],tmp[,1]),2]  >= cutoff
            if(!is.na(oner)){
                if(oner) {oneg <- union(oneg,Tset[i]); usoneg <- union(usoneg,Tset[i]);}
                nLeo[j+8,1] <- nLeo[j+8,1] + oner
                if(oner & (ctmp[match(Tset[i],ctmp[,1]),2] >= cutoff)) nLeo[j+8,2] <- nLeo[j+8,2] + 1
            }
        }       
    }
    
    length(oneg)
    write.table(oneg,file="leaveone_genes_all.txt",quote=FALSE,col.names=FALSE,row.names=FALSE)
    write.table(usoneg,file="leaveone_genes_us.txt",quote=FALSE,col.names=FALSE,row.names=FALSE)
    ### genes in co-expression network degrees and betweenness 
    source("Multi_net.R")
    source("CRF_build.R")
    netex <- build_net(7,fileexp="PCGCall.txt",net_e=FALSE)
    deV <- colSums(netex$matrix)
    beV <- read.table("data/Network_betweenness/Betweenness_node_coexp.txt")
    
    pdf(file="plot/MeaNet.pdf",width=10,height=7)
    a1 <- max(deV[netex$node %in% Tset])
    a2 <- max(beV[beV[,1] %in% Tset,2])
    nulist <- list()
    nulist[[1]] <- deV[netex$node %in% oneg]
    nulist[[2]] <- deV[netex$node %in% setdiff(Tset,oneg)]
    
    nu1list <- list()
    nu1list[[1]] <- log10(beV[beV[,1] %in% oneg,2])
    nu1list[[2]] <- log10(beV[beV[,1] %in% setdiff(Tset,oneg),2])
    
    par(mfrow=c(1,2))
    boxplot(nulist,names=c("CTP","NTP"),main = "Gene interactions", xlab = "", ylab = "Degree")
    text(1.5,35,"t test pvalue: 0.0024")
    boxplot(nu1list,names=c("CTP","NTP"),main= "Gene Betweenness", xlab="", ylab="Betweenness (log10)")
    text(1.5,1,"t test pvalue: 0.00088")
    dev.off()   
    
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
    sum(rowSums(metam[metam[,1] %in% oneg,2:3])>0)
    sum(rowSums(metam[metam[,1] %in% setdiff(Tset,oneg),2:3])>0)
    sum(rowSums(metam[metam[,1] %in% Tset,2:3])>0)
    a <- matrix(c(14,23,22,56),2,2)
    fisher.test(t(a))
    
    
    
    
    
    ### plot figures
    pdf(file="plot/FDRs6_12.pdf",width=10,height=10)
    subs <- c(1,2,3,5,6,10,11)
    par(mai=c(4,2,2,2))
    legend=c("TADA","DAWN","MAGI","HotNet2-STRING","HotNet2-iRef","HotNet2-coexp","HotNet2-Infer","HotNet2-GNAT","MICRF-STRING","MICRF-iRef","MICRF-coexp","MICRF-Infer","MICRF-GNAT") 
    #legend=c("TADA","DAWN","MAGI","RWR-STRING","RWR-iRef","RWR-coexp","RWR-Infer","RWR-GNAT","MICRF-STRING","MICRF-iRef","MICRF-coexp","MICRF-Infer","MICRF-GNAT")  
    barplot(t(nLeo)[,subs], names.arg=legend[subs], beside=T, ylab="#gold risk genes", 
            cex.names=1.5, las=2, ylim=c(0,max(nLeo)+3), col=c("red","darkblue"),cex.lab=1.5,cex.axis=1.25)
    #box(bty="l")
    legend("topleft",col=c("darkblue","red"),legend=c("ASD case","SSC control"),cex=1.3,fill=c("red","darkblue")) 
    dev.off()    
    

    #######  independent samples  ===================================
    source("enrichana_6_12.R")
    ReM <- matrix(0,20,13)
    load("BestP_6_10") 
    j=3
    for(i in 1:20){
        onefile <- paste("ASD/randset4/ASD2sub",j,"derest",i,".csv",sep="")
        oneresult <- read.csv(onefile)
        
        subs <- rowSums(oneresult[match(Tset,oneresult[,1]),c("dn.LoF","dn.mis3")])>0
        reTset <- Tset[subs]
        
        ### TADA method
        TADAfile <- paste("ASD/TADAresult/randset4/TADAdenovo_ASD1sub",j,"depart",i,".csv",sep="")
        subresult <- read.csv(TADAfile)
        cutoff <- total_genes_cutoff(TADAfile,Tset,muts="",alpha=0.5)
        tmpg <- subresult[subresult[,"qvalue.dn"] <= 1-cutoff,1]
        ReM[i,1] <- length(intersect(tmpg,reTset))/length(reTset)##sum(rowSums(oneresult[match(tmpg,oneresult[,1]),c("dn.LoF","dn.mis3")])>0)/length(tmpg)
        
        ### DAWN method
        DAWNfile <- paste("result/randresult4/DAWN/DAWN_randset4_",i,".csv",sep="")
        cutoff <- total_genes_cutoff(DAWNfile,Tset,muts="",alpha=0.5)
        tmp <- read.csv(DAWNfile)
        tmp[is.na(tmp[,"Stage2_posterior"]),"Stage2_posterior"] <- 1
        tmpg <- tmp[tmp[,"Stage2_posterior"] <= 1-cutoff,"Gene"]
        ReM[i,2] <- length(intersect(tmpg,reTset))/length(reTset)##sum(rowSums(oneresult[match(tmpg,oneresult[,1]),c("dn.LoF","dn.mis3")])>0)/length(tmpg)
        
        ### MAGI method 
        MAGIfile <- paste("result/randresult4/MAGI/RandomGeneList.",i,sep="")
        cutoff <- total_genes_cutoff(MAGIfile,Tset,muts="",alpha=0.5)
        tmp <- read.table(MAGIfile)
        tmpg <- tmp[tmp[,2] >= cutoff,1]
        ReM[i,3] <- length(intersect(tmpg,reTset))/length(reTset)##sum(rowSums(oneresult[match(tmpg,oneresult[,1]),c("dn.LoF","dn.mis3")])>0)/length(tmpg)
            
        ### hotnet2 method and our method
        for(netk in 1:5){
            ### hotnet2
            hotfile <- paste("result/randresult4/",netstr[netk],"hotnetresult1part3_",i,netnum[netk],".txt",sep="")
            cutoff <- total_genes_cutoff(hotfile,Tset,muts="",alpha=0.5)
            tmp <- read.table(hotfile)
            tmpg <- tmp[tmp[,2] >= cutoff,1]
            tmpg <- intersect(tmpg,oneresult[,1])
            ReM[i,netk+3] <- length(intersect(tmpg,reTset))/length(reTset)##sum(rowSums(oneresult[match(tmpg,oneresult[,1]),c("dn.LoF","dn.mis3")])>0)/length(tmpg)
          
            ### our method
            myfile <- paste("result/randresult4_",BestP[netk,2],"/",netstr[netk],"CRFresult_",BestP[netk,3],"part3_",i,"LBP_",netnum[netk],".txt",sep="")
            cutoff <- total_genes_cutoff(myfile,Tset,muts="",alpha=0.5)
            tmp <- read.table(myfile)
            tmpg <- tmp[tmp[,2] >= cutoff,1]
            ReM[i,netk+8] <- length(intersect(tmpg,reTset))/length(reTset)##sum(rowSums(oneresult[match(tmpg,oneresult[,1]),c("dn.LoF","dn.mis3")])>0)/length(tmpg)
        }
    }
    
    
    source("enrichana_6_12.R")
    load("BestP_6_10") 
    j=3
    a <- rep(0,20)
    for(i in 1:20){
        onefile <- paste("ASD/randset4/ASD2sub",j,"derest",i,".csv",sep="")
        oneresult <- read.csv(onefile)
        TADAfile <- paste("ASD/TADAresult/randset4/TADAdenovo_ASD1sub",j,"depart",i,".csv",sep="")
        subresult <- read.csv(TADAfile)
        
        subs <- rowSums(subresult[,c("dn.LoF","dn.mis3")])>0
        genes <- subresult[subs,1]
        
        subs <- rowSums(oneresult[,c("dn.LoF","dn.mis3")])>0
        genes1 <- oneresult[subs,1]
        
        a[i] <- length(intersect(genes1,genes))/length(genes)
    }
    
    
    ### plot figures ===================================
    pdf(file="plot/sub6_12_f.pdf",width=10,height=7)
    subs <- c(1,2,3,5,6,10,11)
    par(mai=c(2,1,1,1))
    legend=c("TADA","DAWN","MAGI","HotNet2-STRING","HotNet2-iRef","HotNet2-coexp","HotNet2-Infer","HotNet2-GNAT","MICRF-STRING","MICRF-iRef","MICRF-coexp","MICRF-Infer","MICRF-GNAT")
    #legend=c("TADA","DAWN","MAGI","RWR-STRING","RWR-iRef","RWR-coexp","RWR-Infer","RWR-GNAT","MICRF-STRING","MICRF-iRef","MICRF-coexp","MICRF-Infer","MICRF-GNAT")
    boxplot(ReM[,subs],col=2:(length(subs)+1),names=legend[subs],las=2,ylab="Fraction of recurrent genes",cex.axis=1.5,cex.lab=2,cex.names=2)
    dev.off()
    
}

allPf <- function(TPcut,TADAFile="ASD/TADAresult/TADAdenovo_ASD4_16.csv"){
    
    source("enrichana_6_12.R")
    TADAall <- read.csv(TADAFile)
    Tset <- TADAall[TADAall[,"qvalue.dn"]<= TPcut,1]
    
    Tset
}

filenamesf <- function(j,i,betaV,DAWNnames,Maginames){
    
    netstr <- c("STRING/","iRef/","coexp/","Infer/","GNAT/")
    netnum <- c(31,6,7,8,21)
    
    tadastr <- "ASD/TADAresult/randset_1/TADAdenovo_ASD"
    DAWNstr <- "DAWN/DAWN_randset"
    MAGIstr <- "MAGI/RandomGeneList."
    resultstr <- "result/randresult_1/"
    
    nuk <- c(0,0)
    if(i==0 & j==0){
        tmpstr <- "ASD4_16"
        tmpstr1 <- "4_16"
        nuk <- c(1,1)
    }else{
        tmpstr <- paste("rand1",j,"_",i,sep="")
        tmpstr1 <- paste("rand1_",j,"_",i,sep="")
        nuk[1] <- which(DAWNnames==paste("TADAdenovo_ASD",tmpstr1,".csv",sep=""))
        nuk[2] <- which(Maginames==paste("ASD",tmpstr1,".csv",sep=""))
    }
    
    TADAfile <- paste(tadastr,tmpstr1,".csv",sep="")
    DAWNfile <- paste(resultstr,DAWNstr,nuk[1],".csv",sep="")
    MAGIfile <- paste(resultstr,MAGIstr,nuk[2],sep="") 
    
    nnet <- length(netstr)
    hotnetfile <- 1:nnet
    misfile <- 1:nnet
    for(k in 1:nnet){
        hotnetfile[k] <- paste(resultstr,netstr[k],"hotnetresult1",tmpstr,netnum[k],".txt",sep="")
        misfile[k] <- paste(resultstr,netstr[k],"CRFresult_",betaV[k],tmpstr,"LBP_",netnum[k],".txt",sep="")
    }
    
    filenames <- c(TADAfile,DAWNfile,MAGIfile,hotnetfile,misfile)
    
    filenames
    
}

gene_select <- function(filenames, allP){
    
    Tset <- allP
    ## all positive samples could be predicted by all methods
    for(k in 1:length(filenames)){
        filename <- filenames[k]
        if(grepl(".csv",filename)){ result <- read.csv(filename,stringsAsFactors=FALSE); 
        }else if(grepl(".txt",filename)){result <- read.table(filename);}
        allP <- intersect(allP,result[,1])
        #print(length(allP))
    } 
    ## all negative samples could be predicted by all methods
    allN <- read.csv(filenames[1])[,1]
    for(k in 1:length(filenames)){
        filename <- filenames[k]
        if(grepl(".csv",filename)){ result <- read.csv(filename); 
        }else if(grepl(".txt",filename)){result <- read.table(filename);}
        allN <- intersect(allN,result[,1])
        #print(length(allN))
    } 
    allN <- setdiff(allN,Tset)
    genes <- union(allP,allN)
    
    list(genes=genes,allP=allP,allN=allN)
}

ROC_ex <- function(filenames,allP,genes,flag=2,xlim=c(0,1),ylim=c(0,1)){
    
    #  flag==1: only rank information is used
    #  flag==2: FDR information is used
    
    library(ROCR)
    library(pROC)
    N <- length(genes)
    predictions <- list()
    labels <- list()
    ped0 <- 1 - (1:N)/N
    
    for(i in 1:length(filenames)){
        filename <- filenames[i]
        if(grepl(".csv",filename)){ result <- read.csv(filename); 
        }else if(grepl(".txt",filename)){result <- read.table(filename);
        }else{result <- read.table(filename);}
        
        result <- result[!duplicated(result[,1]),]
        result <- result[result[,1] %in% genes,]
        
        if(flag==1){
            predictions[[i]] <- ped0
        }else if(flag==2){
            if(grepl(".csv",filename)){
                predictions[[i]] <- 1 - result[,"qvalue.dn"]
            }else{
                #predictions[[i]] <- 1 - result[,5]
                predictions[[i]] <- result[,2]/max(result[,2])
            }
        }
        
        labels[[i]] <- rep(-1,length(predictions[[i]]))
        labels[[i]][match(allP,result[,1])] <- 1
    }

    pred <- prediction(predictions, labels)
    auc <- unlist(performance(pred, "auc")@y.values)
    perf <- performance(pred, measure = "tpr", x.measure = "fpr") 

        cols <- c("Black","Blue","Cyan","Brown","Coral","DarkViolet","Green")
        plot(perf,main="",xlab="False positive rate",ylab="True positive rate",cex.lab=2,cex.main=2,xlim=xlim,ylim=ylim,xaxt="n",yaxt="n")
        for(i in 1:length(filenames)){
            if(i==1){
                lines(perf@x.values[[i]], perf@y.values[[i]], col = cols[i],lwd=2)
            }else{
                lines(perf@x.values[[i]], perf@y.values[[i]], col = cols[i],lwd=1)
            }
        }
        lines(x=c(0,1),y=c(0,1),col=i+1,lty=2)
        #axis(1,at=seq(0,1,0.1))
        #axis(2,at=seq(0,1,0.1))
        legend=c("TADA","MICRF-iRef","MICRF-coexp") ###
        legend <- paste(legend,format(auc,digits=3),sep="--AUC:")
    if(max(xlim)>=0.5){
        legend("bottomright",col=cols[1:length(filenames)],legend=legend,lty=rep(1,length(filenames)),lwd=rep(2,length(filenames)),cex=2,y.intersp=0.7) 
        a=0; b=0.1; c=0.3; d=0.6;
        red_box(a,b,c,d);
        a=0.05; b=0.3; c=0.8; d=1;
        red_box(a,b,c,d)
    }
        
} 

red_box <- function(a,b,c,d){
    lines(c(a,a),c(c,d),col="red",type="l")
    lines(c(a,b),c(d,d),col="red",type="l")
    lines(c(b,b),c(c,d),col="red",type="l")
    lines(c(a,b),c(c,c),col="red",type="l")
}

new_auc <- function(filenames,allP,genes,flag=2,iplot=FALSE,DAWNflag=0){
    
    #  flag==1: only rank information is used
    #  flag==2: FDR information is used
    
    library(ROCR)
    library(pROC)
    N <- length(genes)
    predictions <- list()
    labels <- list()
    ped0 <- 1 - (1:N)/N
    
    DAWNr <- read.csv(filenames[2])
    subs <- is.na(DAWNr[,"FDR"])
    DAWNr1 <- DAWNr[!subs,]
    DAWNr2 <- DAWNr[subs,]
    if(DAWNflag==0){
        DAWNr2 <- DAWNr2[order(DAWNr2[,"Stage1_posterior"]),]
        DAWNr <- rbind(DAWNr1,DAWNr2)
        DAWNr <- DAWNr[!duplicated(DAWNr[,1]),]
        DAWNr <- DAWNr[DAWNr[,"Gene"] %in% genes,]
        DAWNr[is.na(DAWNr[,"FDR"]),"FDR"] <- 1
        DAWNr[is.na(DAWNr[,"Stage2_posterior"]),"Stage2_posterior"] <- 1
    }else if(DAWNflag==1){
        DAWNr <- DAWNr1
        DAWNr <- DAWNr[!duplicated(DAWNr[,1]),]
        DAWNr <- DAWNr[DAWNr[,"Gene"] %in% genes,]
    }
    
    for(i in c(1,3:length(filenames))){
        filename <- filenames[i]
        if(grepl(".csv",filename)){ result <- read.csv(filename); 
        }else if(grepl(".txt",filename)){result <- read.table(filename);
        }else{result <- read.table(filename);}
        
        result <- result[!duplicated(result[,1]),]
        result <- result[result[,1] %in% genes,]
        
        if(flag==1){
            predictions[[i]] <- ped0
        }else if(flag==2){
            if(grepl(".csv",filename)){
                predictions[[i]] <- 1 - result[,"qvalue.dn"]
            }else{
                #predictions[[i]] <- 1 - result[,5]
                predictions[[i]] <- result[,2]/max(result[,2])
            }
        }
        
        labels[[i]] <- rep(-1,length(predictions[[i]]))
        labels[[i]][match(allP,result[,1])] <- 1
    }
    
    i=2
    if(flag==1){
        predictions[[i]] <- ped0
    }else if(flag==2){
        predictions[[i]] <- 1 - DAWNr[,"Stage2_posterior"] #1 - DAWNr[,"FDR"]
    }
    labels[[i]] <- rep(-1,length(predictions[[i]]))
    labels[[i]][match(allP,DAWNr[,1])] <- 1
    
    pred <- prediction(predictions, labels)
    auc <- unlist(performance(pred, "auc")@y.values)
    perf <- performance(pred, measure = "tpr", x.measure = "fpr") 
    FPR <- list()
    TPR <- list()
    for(i in 1:length(filenames)){
        FPR[[i]] <- perf@x.values[[i]]
        TPR[[i]] <- perf@y.values[[i]]   
    }
    
    if(iplot){
        cols <- c("Black","Blue","Cyan","Brown","Coral","DarkViolet","Green")
        plot(perf,main="Comparison of Performance",xlab="False positive rate",ylab="True positive rate",cex.lab=2,cex.main=2,xlim=c(0,1),ylim=c(0,1))
        for(i in 1:length(filenames)){
            if(i==1){
                lines(perf@x.values[[i]], perf@y.values[[i]], col = cols[i],lwd=2)
            }else{
                lines(perf@x.values[[i]], perf@y.values[[i]], col = cols[i],lwd=1)
            }
        }
        lines(x=c(0,1),y=c(0,1),col=i+1,lty=2)
        axis(1,at=seq(0,1,0.1))
        axis(2,at=seq(0,1,0.1))
        legend=c("TADA","DAWN","MAGI","HotNet2-iRef","HotNet2-coexp","MICRF-iRef","MICRF-coexp") ###
        
        ##legend=c("TADA","DAWN","MAGI","RWR-STRING","RWR-iRef","RWR-coexp","RWR-Infer","RWR-GNAT","MICRF-STRING","MICRF-iRef","MICRF-coexp","MICRF-Infer","MICRF-GNAT") ###
        legend <- paste(legend,format(auc,digits=3),sep="--AUC:")
        legend("bottomright",col=cols[1:length(filenames)],legend=legend,lty=rep(1,length(filenames)),lwd=rep(2,length(filenames)),cex=2,y.intersp=0.7) 
        
        #         for(i in 1:length(filenames)){
        #             if(i==1) {
        #                 plot.roc(labels[[i]],predictions[[i]])
        #             }else{
        #                 plot.roc(labels[[i]],predictions[[i]], add=TRUE)
        #             }
        #         }  
        #plot.new()
        #myplot(predictions,labels)
    }
    
    list(auc=auc,FPR=FPR,TPR=TPR)
}

plot_err <- function(x,y,sd,flag=1,main="ROC",xlab="Samples",ylab="AUC",xlim=c(0,1),ylim=c(0,1)){
    cols <- c("Black","Blue","Cyan","Brown","Coral","DarkViolet","Green")
    if(flag==1){
        par(mar=c(6,6,3,2))
        plot(x, y,xaxt="n", type="b",col=cols[flag], main=main, ylab=ylab,xlab=xlab,xlim=xlim,ylim=ylim,cex.lab=2,cex.main=2,lwd=3,cex.axis=1.8,pch=10)
        axis(1,at=seq(0.2,0.9,0.1)[1:length(x)],labels=floor(3953*seq(0.2,0.9,0.1))[1:length(x)],font=1,cex.axis=1.5)
    }else{
        lines(x, y, type="b",col=cols[flag],xaxt="n",lwd=3)
    }
    
    #axis(2,at=seq(ylim[1],ylim[2],0.1))
    #segments(x,y-sd,x,y+sd)
    #epsilon <- 0.004
    #segments(x-epsilon,y-sd,x+epsilon,y-sd)
    #segments(x-epsilon,y+sd,x+epsilon,y+sd)
    
}

total_genes_cutoff <- function(filenames,Tset,muts,alpha=0.5){
    
    library(ROCR)
    predictions <- list()
    labels <- list()
    result <- list()
    n <- length(filenames)
    cutoff <- 1:n
    
    for(i in 1:length(filenames)){
        filename <- filenames[i]
        if(grepl(".csv",filename)){ tmp <- read.csv(filename); 
        }else if(grepl(".txt",filename)){tmp <- read.table(filename);}else{
            tmp <- read.table(filename)
        }
        
        tmpn <- length(unlist(strsplit(filename,"/")))
        subfilename <- unlist(strsplit(filename,"/"))[tmpn]
        
        if(grepl("DAWN",subfilename)){
            subs <- is.na(tmp[,"FDR"])
            result[[i]] <- tmp[!subs,]
        }else{
            #result[[i]] <- tmp[!duplicated(tmp[,1]),]
            result[[i]] <- tmp       
        }
        
        if(grepl("TADA",subfilename)){
            predictions[[i]] <- 1 - result[[i]][,"qvalue.dn"]
        }else if(grepl("DAWN",subfilename)){
            predictions[[i]] <- 1 - result[[i]][,"Stage2_posterior"]     
        }else if(grepl("Ptest",subfilename)){
            predictions[[i]] <- result[[i]][,"score"]
        }else{
            predictions[[i]] <- result[[i]][,2] #1 - result[[i]][,5]
        }
        
        labels[[i]] <- rep(-1,length(predictions[[i]]))
        labels[[i]][match(Tset,result[[i]][,1])] <- 1
        
        pred <- prediction(predictions[[i]], labels[[i]])
        Fmea <- performance(pred, "f",alpha=alpha)@y.values
        Cutoffs <- performance(pred, "f",alpha=alpha)@x.values
        Fmea[[1]][is.na(Fmea[[1]])] <- 0
        subs <- which.max(Fmea[[1]])
        cutoff[i] <- Cutoffs[[1]][subs]
    }
    
    cutoff
}
