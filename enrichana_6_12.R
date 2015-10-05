randset_1_AUC_6_3 <- function(){
    
    source("enrichana_6_12.R")
    
    load("BestP_8_19_4")
    betaV <- BestP[,3]
    TADAFile="../TADA_DAWN/result/TADAdenovo_meta_dmis.csv"
    Tset <- allPf(TPcut=0.1,TADAFile)

    DAWNnames <-  readLines(con <- file("../TADA_DAWN/DAWN_package/TADAdenovo_randset_1.txt","r"))   
    close(con)
    Maginames <-  readLines(con <- file("../MAGI_V0.1_Alpha/mydata/randset_1/MAGIfilenames.txt","r"))   
    close(con)

    
    nnet <- 8
    AUC <- array(0,dim=c(8,20,3+nnet*2)) 
    FPR <- array(0,dim=c(8,20,3+nnet*2))
    TPR <- array(0,dim=c(8,20,3+nnet*2))

    DAWNflag=1
    for(j in 2:9){
        for(i in 1:20){   
            filenames <- filenamesf(j,i,betaV=betaV,DAWNnames,Maginames)
            tmp <- gene_select_8_13(filenames,Tset)
            genes <- tmp$genes
            allP <- tmp$allP
            allN <- tmp$allN
            
            if(DAWNflag==1){
                tmp <- read.csv(filenames[2])
                genes <- intersect(genes,tmp[!is.na(tmp[,"FDR"]),"Gene"])
                allP <- intersect(allP,genes)
            }
            
            tmp <- new_auc_8_13(filenames,allP,genes)
            AUC[j-1,i,] <- tmp$auc
            
            for(kk in 1:length(filenames)){
                TPR[j-1,i,kk] <- tmp$TPR[[kk]][which.min(abs(tmp$FPR[[kk]]-0.2))]
                FPR[j-1,i,kk] <- tmp$TPR[[kk]][which.min(abs(tmp$FPR[[kk]]-0.1))]
            }
        }
    }
    
    save(AUC,file="AUC_m_8_19_4")
    save(TPR,file="TPR_m_8_19_4")
    save(FPR,file="FRP_m_8_19_4")

}

plotAUC <- function(aucfile,aucf){

    ## plot auc values
    ## Figure 1: AUC values for eight random sets
        source("enrichana_6_12.R")
       
        pdf(file=aucf,width=12,height=10)    
        load(aucfile)
        #cols <- c("black","red","green","aquamarine3","bisque4","blueviolet","brown","chartreuse","chocolate4","cyan","coral4","aquamarine3","bisque4","blueviolet","brown","chartreuse","chocolate4","cyan","coral4")
        cols <- c("black","red","green","mediumorchid1","orchid","blueviolet","brown","chartreuse","chocolate4","cyan","coral4","orangered","blue","blue","darkslategray1","","chocolate4","coral4","")
        #legend <- c("TADA","DAWN","MAGI","RWR-STRING","RWR-iRef","RWR-coexp","RWR-HPRD","RWR-corr","RWR-corr1","RWR-coexp5","RWR-coexp7","MICRF-STRING","MICRF-iRef","MICRF-coexp","MICRF-HPRD","MICRF-corr","MICRF-corr1","MICRF-coexp5","MICRF-coexp7")
        
        #subs <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19)
        cols <- cols[c(1,2,3,4,5,7,10,12,13,15,18)]
        subs <- 1:11 #c(1,2,3,4,5,7,10,12,13,15,18)
        nn <- subs
        nmeth=length(subs)
        legend <- c("TADA","DAWN","MAGI","RWR-STRING","RWR-iRef","RWR-HPRD","RWR-coexp","MICRF-STRING","MICRF-iRef","MICRF-HPRD","MICRF-coexp")
        
        main="AUC performance for different trios"
        xlab="The number of trios";ylab="AUC";
        ymint <- 1;ymaxt <- 0;
        ymint <- min(AUC);ymaxt <- max(AUC);
        ymin=max(min(ymint)+0.1,0);ymax = min(ymaxt+0.1,1);
        x <- seq(0.2,0.9,0.1)
        Y1 <- rep(0,nmeth)
        
        k=1
        for(i in subs){
            y <- rowMeans(AUC[,,i])
            sd <- apply( AUC[,,i], 1, sd)/2 
            plot_err(x[1:5],y[1:5],sd,i,main,xlab,ylab,ylim=c(0,ymax),xlim=c(0.15,0.62),cols=cols)
            Y1[k] <- y[1]
            k <- k+1
        } 
        Y2 <- round(Y1*100)/100
        labels <- 1:length(nn)
    
        tx <- c(0.193,0.2,0.193,0.2,0.193,0.2,0.193,0.2,0.188,0.2,0.2)
        ty <- -sort(-Y2)
        ty[5] = ty[6] + 0.001
        ty[6] = ty[6] + 0.001
#         ty[3] = ty[3] - 0.01
#         ty[1] = ty[1] + 0.01
#         ty[4] = ty[4] + 0.01
#         ty[6] = ty[6] - 0.01
#         ty[11] = ty[11] - 0.01
#         ty[12] = ty[12] - 0.01
        
        cols <- cols[subs]
        ords <- sort(-Y1,index.return=T)$ix
        text(tx,ty,labels=labels,pos=2,col=cols[ords],cex=1.3)
        legend("bottomright",col=cols[ords],legend=paste(labels,legend[ords],sep=":"),lty=rep(1,length(subs)),lwd=rep(2,length(subs)),cex=1.5,y.intersp=1)
    
        dev.off()

}

randset_1_AUC_7_29 <- function(){
    source("enrichana_6_12.R")
    
    load("BestP_8_13")
    betaV <- BestP[,3]
    TADAFile="../TADA_DAWN/result/TADAdenovo_meta_dmis.csv"
    Tset <- allPf(TPcut=0.1,TADAFile)
  
    DAWNnames <-  readLines(con <- file("../TADA_DAWN/DAWN_package/TADAdenovo_randset_1.txt","r"))   
    close(con)
    Maginames <-  readLines(con <- file("../MAGI_V0.1_Alpha/mydata/randset_1/MAGIfilenames.txt","r"))   
    close(con)
    
    nnet <- 8
    AUC <- array(0,dim=c(8,20,3+nnet*2)) 
    for(j in 2:9){
        for(i in 1:20){   
            filenames <- filenamesf(j,i,betaV=betaV,DAWNnames,Maginames)
            tmp <- new_auc_7_29(filenames,Tset)
            AUC[j-1,i,] <- tmp
        }
    }
    
    save(AUC,file="AUC_all_8_14")
    
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
    save(PerM,file="PerM_7_23")
    
    BestP <- matrix(0,5,3)
    for(i in 1:5){
        BestP[i,1] <- netflags[i]
        tmp <- which(PerM[i,,]==max(PerM[i,,]),arr.ind=TRUE)
        BestP[i,2] <- tmp[1]
        BestP[i,3] <- betaV[tmp[2]]
    }
    BestP <- BestP[c(1,2,4,3,5),]
    save(BestP,file="BestP_7_23")
}

batch_randset4_7_22 <- function(){
    source("enrichana_6_12.R")
    betaV <- c(0,0.2,0.4,0.5,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0)
    nmeth <- 4
    
    TPcut=0.1
    TADAFile <- "../TADA_DAWN/result/TADAdenovo_meta_dmis.csv"
    asdall <- read.csv(TADAFile)
    asdall[is.na(asdall[,"qvalue.dn"]),"qvalue.dn"] <- 1
    Tset <- asdall[asdall[,"qvalue.dn"]< TPcut,1]
    
    netstr <- c("STRING/","iRef/","coexp/","HPRD/","coexp1/","corr/","coexp5/","coexp7/")
    netflags <- c(31,6,7,20,21,26,27,28)
    nodeN <- c(12969,12124,9448,10526,8782,9801,17587,17587)
    
    PerM <- array(0,c(8,4,12))
    
    ReVM <- list()
    HiVM <- list()
    
    k=1
    for(netflag in netflags){
        ReVM[[k]] <- array(0,c(nmeth,length(betaV),nodeN[k]))
        HiVM[[k]] <- array(0,c(nmeth,length(betaV),nodeN[k]))
        for(meth in 1:nmeth){
            for(beta in betaV){
                netsub <- which(netflags==netflag)
                ## recurrent
                path <- paste("result/randresult4_",meth,"/",netstr[netsub],sep="")
                ReM <- matrix(0,20,nodeN[netsub])
                ReM <- randset4_one(netflag,meth,beta,path,ReM)
                ReV <- colMeans(ReM)
                
                ## leave-one rank
                path <- paste("result/leaveone4result_",meth,"/",netstr[netsub],sep="")
                HiV <- matrix(0,1,nodeN[netsub])
                HiV <- leaveone4_one(netflag,meth,beta,Tset,path,HiV)
                     
                tmp <- ReV*HiV
                
                ReVM[[k]][meth,which(betaV==beta),] <- ReV
                HiVM[[k]][meth,which(betaV==beta),] <- HiV
                
                PerM[which(netflag==netflags),meth,which(betaV==beta)] <- max(tmp[10:length(tmp)])
            }
        }
        k=k+1
    } 
    
    save(ReVM,file="ReVM")
    save(HiVM,file="HiVM")
    save(PerM,file="PerM_8_17")

    BestP <- matrix(0,8,3)
    for(i in 1:8){
        BestP[i,1] <- netflags[i]
        tmp <- which(PerM[i,,]==max(PerM[i,,]),arr.ind=TRUE)
        BestP[i,2] <- tmp[1]
        BestP[i,3] <- betaV[tmp[2]]
    }
    save(BestP,file="BestP_8_17")

}

BestP_8_17 <- function(k){
    source("enrichana_6_12.R")
    betaV <- c(0,0.2,0.4,0.5,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0)
    nmeth <- 4
    
    TPcut=0.1
    TADAFile <- "../TADA_DAWN/result/TADAdenovo_meta_dmis.csv"
    asdall <- read.csv(TADAFile)
    asdall[is.na(asdall[,"qvalue.dn"]),"qvalue.dn"] <- 1
    Tset <- asdall[asdall[,"qvalue.dn"]< TPcut,1]
    
    netstr <- c("STRING/","iRef/","coexp/","HPRD/","coexp1/","corr/","coexp5/","coexp7/")
    netflags <- c(31,6,7,20,21,26,27,28)
    nodeN <- c(12969,12124,9448,10526,8782,9801,17587,17587)
    
    PerM <- array(0,c(8,4,12))
        
    netflag=netflags[k]
    ReVM <- array(0,c(nmeth,length(betaV),nodeN[k]))
    HiVM <- array(0,c(nmeth,length(betaV),nodeN[k]))
        for(meth in 1:nmeth){
            for(beta in betaV){
                netsub <- which(netflags==netflag)
                ## recurrent
                path <- paste("result/randresult4_",meth,"/",netstr[netsub],sep="")
                ReM <- matrix(0,20,nodeN[netsub])
                ReM <- randset4_one(netflag,meth,beta,path,ReM)
                ReV <- colMeans(ReM)
                
                ## leave-one rank
                path <- paste("result/leaveone4result_",meth,"/",netstr[netsub],sep="")
                HiV <- matrix(0,1,nodeN[netsub])
                HiV <- leaveone4_one(netflag,meth,beta,Tset,path,HiV)
                
                tmp <- ReV*HiV
                
                ReVM[meth,which(betaV==beta),] <- ReV
                HiVM[meth,which(betaV==beta),] <- HiV
                
                PerM[which(netflag==netflags),meth,which(betaV==beta)] <- max(tmp[10:length(tmp)])
            }
        }
    
    save(ReVM,file=paste("ReVM",k,sep="_"))
    save(HiVM,file=paste("HiVM",k,sep="_"))
    save(PerM,file=paste("PerM_8_17",k,sep="_"))

}

dealBestP <- function(){
    source("enrichana_6_12.R")
    betaV <- c(0,0.2,0.4,0.5,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0)
    nmeth <- 4
    
    TPcut=0.1
    TADAFile <- "../TADA_DAWN/result/TADAdenovo_meta_dmis.csv"
    asdall <- read.csv(TADAFile)
    asdall[is.na(asdall[,"qvalue.dn"]),"qvalue.dn"] <- 1
    Tset <- asdall[asdall[,"qvalue.dn"]< TPcut,1]
    
    netstr <- c("STRING/","iRef/","coexp/","HPRD/","coexp1/","corr/","coexp5/","coexp7/")
    netflags <- c(31,6,7,20,21,26,27,28)
    nodeN <- c(12969,12124,9448,10526,8782,9801,17587,17587)
    
    PerM <- array(0,c(8,4,12))
    
    for(k in 1:length(netflags)){
        netflag=netflags[k]
        load(paste("ReVM",k,sep="_"))
        load(paste("HiVM",k,sep="_"))
        
        for(meth in 1:nmeth){
            for(beta in betaV){
            ReV <- ReVM[meth,which(betaV==beta),] ## recurrent
            HiV <- HiVM[meth,which(betaV==beta),]## leave-one rank
            tmp <- ReV*HiV
            
            PerM[which(netflag==netflags),meth,which(betaV==beta)] <- max(tmp[10:300])
            }
        }
    }

    BestP <- matrix(0,8,3)
    for(i in 1:8){
        BestP[i,1] <- netflags[i]
        tmp <- which(PerM[i,,]==max(PerM[i,,]),arr.ind=TRUE)
        BestP[i,2] <- tmp[1,1]
        BestP[i,3] <- betaV[tmp[1,2]]
    }
    save(BestP,file="BestP_8_13")

}

randset4_one <- function(netflag,meth,beta,path,ReM){
    
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
    
#     TADAfile <- paste("ASD/TADAresult/randset4/TADAdenovo_ASD1sub",j,"depart",i,".csv",sep="")
#     subresult <- read.csv(TADAfile)
    
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
    netstr <- c("STRING/","iRef/","HPRD/","corr/","coexp5/")
    netnum <- c(31,6,20,26,27)
    TADAFile <- "../TADA_DAWN/result/TADAdenovo_meta_dmis.csv"
    asdall <- read.csv(TADAFile)
    asdall[is.na(asdall[,"qvalue.dn"]),"qvalue.dn"] <- 1
    Tset <- asdall[asdall[,"qvalue.dn"]<=TPcut,1]
    
    ### leave one mutation ==========================================
    ## see addresult leaveone_recurrent

    #######  independent samples  ===================================
    source("enrichana_6_12.R")
    ReM <- matrix(0,20,13)
    load("BestP_8_19_4") 
    BestP <- BestP[c(1,2,4,6,7),]
    j=3
    for(i in 1:20){
        onefile <- paste("DDD_mutations/randset4/ASD2sub",j,"derest",i,".csv",sep="")
        oneresult <- read.csv(onefile)
        
        subs <- rowSums(oneresult[match(Tset,oneresult[,1]),c("dn.LoF","dn.mis3")])>0
        reTset <- Tset[subs]
        
        ### TADA method
        TADAfile <- paste("../TADA_DAWN/result/randset4/TADAdenovo_ASD1sub",j,"depart",i,".csv",sep="")
        subresult <- read.csv(TADAfile)
        cutoff <- total_genes_cutoff(TADAfile,Tset,muts="",alpha=0.5)
        tmpg <- subresult[subresult[,"qvalue.dn"] <= 1-cutoff,1]
        ReM[i,1] <- length(intersect(tmpg,reTset))/length(reTset)##sum(rowSums(oneresult[match(tmpg,oneresult[,1]),c("dn.LoF","dn.mis3")])>0)/length(tmpg)
        
        ### DAWN method
        DAWNfile <- paste("../TADA_DAWN/DAWN_package/result/randset4/DAWN_randset4_",i,".csv",sep="")
        cutoff <- total_genes_cutoff(DAWNfile,Tset,muts="",alpha=0.5)
        tmp <- read.csv(DAWNfile)
        tmp[is.na(tmp[,"Stage2_posterior"]),"Stage2_posterior"] <- 1
        tmpg <- tmp[tmp[,"Stage2_posterior"] <= 1-cutoff,"Gene"]
        ReM[i,2] <- length(intersect(tmpg,reTset))/length(reTset)##sum(rowSums(oneresult[match(tmpg,oneresult[,1]),c("dn.LoF","dn.mis3")])>0)/length(tmpg)
        
        ### MAGI method 
        MAGIfile <- paste("../MAGI_V0.1_Alpha/mydata/randset4/RandomGeneList.",i,sep="")
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

    ### plot figures ===================================
    pdf(file="plot/recurr_8_19_4.pdf",width=10,height=7)
    #subs <- c(1,2,3,4,9,7,12,5,10,8,13,6,11)
    ##subs <- c(1,2,3,4,9,5,10,6,11,7,12,8,13)
    subs <- c(1,2,3,4,9,5,10,6,11,8,13)
    par(mai=c(2.3,1,1,1))
    cols <- c("red","green","dodgerblue","orange","orange","yellow","yellow","pink","pink","cyan","cyan","purple","purple")
    legend=c("TADA","DAWN","MAGI","RWR-STRING","RWR-iRef","RWR-HPRD","RWR-corr","RWR-coexp","MICRF-STRING","MICRF-iRef","MICRF-HPRD","MICRF-corr","MICRF-coexp")
    boxplot(ReM[,subs],col=cols,names=legend[subs],las=2,ylab="Fraction of recurrent genes",cex.axis=1.5,cex.lab=1.7,cex.names=2,outline=FALSE)
    dev.off()
    
}

allPf <- function(TPcut,TADAFile="ASD/TADAresult/TADAdenovo_ASD4_16.csv"){
    
    TADAall <- read.csv(TADAFile)
    TADAall[is.na(TADAall[,"qvalue.dn"]),"qvalue.dn"] <- 1
    Tset <- TADAall[TADAall[,"qvalue.dn"]< TPcut,1]
    
    Tset
}

filenamesf <- function(j,i,betaV,DAWNnames,Maginames){
    
    netstr <- c("STRING/","iRef/","coexp/","HPRD/","coexp1/","corr/","coexp5/","coexp7/")
    netnum <- c(31,6,7,20,21,26,27,28)
    
    tadastr <- "../TADA_DAWN/result/randset_1/TADAdenovo_ASD"
    DAWNstr <- "../TADA_DAWN/DAWN_package/result/randset_1/DAWN_randset"
    MAGIstr <- "../MAGI_V0.1_Alpha/mydata/randset_1/RandomGeneList."
    resultstr <- "result/randresult_1/"
    
    nuk <- c(0,0)
    if(i==0 & j==0){
        tmpstr <- ""
        tmpstr1 <- ""
        nuk <- c(1,1)
    }else{
        tmpstr <- paste("rand1",j,"_",i,sep="")
        tmpstr1 <- paste("rand1_",j,"_",i,sep="")
        nuk[1] <- which(DAWNnames==paste("TADAdenovo_ASD",tmpstr1,".csv",sep=""))
        nuk[2] <- which(Maginames==paste("ASD",tmpstr1,".csv",sep=""))
    }
    
    TADAfile <- paste(tadastr,tmpstr1,".csv",sep="")
    DAWNfile <- paste(DAWNstr,nuk[1],".csv",sep="")
    MAGIfile <- paste(MAGIstr,nuk[2],sep="") 
    
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

gene_select_8_13 <- function(filenames, Tset){
    
    genes <- read.csv(filenames[1])[,1]
    for(k in 2:length(filenames)){
        filename <- filenames[k]
        if(grepl(".csv",filename)){ result <- read.csv(filename); 
        }else if(grepl(".txt",filename)){result <- read.table(filename);}
        genes <- intersect(genes,result[,1])
    } 
    allN <- setdiff(genes,Tset)
    allP <- intersect(genes,Tset)
    
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
        plot(x=c(0,1),y=c(0,1),col="black",lty=2,main="",xlab="False positive rate",ylab="True positive rate",cex.axis=2,cex.lab=2,cex.main=2,xlim=xlim,ylim=ylim,type="l")
        for(i in 1:length(filenames)){
            if(i==1){
                lines(perf@x.values[[i]], perf@y.values[[i]], col = cols[i],lwd=2,xaxt="n")
            }else{
                lines(perf@x.values[[i]], perf@y.values[[i]], col = cols[i],lwd=1,xaxt="n")
            }
        }

        legend=c("TADA","MICRF-STRING","MICRF-iRef","MICRF-coexp","MICRF-HPRD","MICRF-corr") ###
        ords <- c(1,2,5,3,6,4)
        legend <- paste(legend,format(auc,digits=3),sep="--AUC:")
    if(max(xlim)>=0.5){
        legend("bottomright",col=cols[ords],legend=legend[ords],lty=rep(1,length(filenames)),lwd=rep(2,length(filenames)),cex=1.8,y.intersp=0.7) 
        a=0; b=0.1; c=0.3; d=0.6;
        red_box(a,b,c,d);
        #a=0.05; b=0.3; c=0.8; d=1;
        #red_box(a,b,c,d)
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

new_auc_8_13 <- function(filenames,allP,genes){
    
    library(ROCR)
    library(pROC)
    N <- length(genes)
    predictions <- list()
    labels <- list()
   
    for(i in 1:length(filenames)){
        filename <- filenames[i]
        if(grepl(".csv",filename)){ result <- read.csv(filename); 
        }else{result <- read.table(filename);}
        
        result <- result[!duplicated(result[,1]),]
        result <- result[result[,1] %in% genes,]
        
        if(i==1){
            predictions[[i]] <- 1 - result[,"qvalue.dn"]
        }else if(i==2){
            result[is.na(result[,"FDR"]),"FDR"] <- 1
            result[is.na(result[,"Stage2_posterior"]),"Stage2_posterior"] <- 1
            predictions[[i]] <- 1 - result[,"Stage2_posterior"]
        }else{
            predictions[[i]] <- result[,2]/max(result[,2])
        }
        
        labels[[i]] <- rep(-1,length(predictions[[i]]))
        labels[[i]][match(allP,result[,1])] <- 1
    }
    #-------------------------------------------------------------------
    pred <- prediction(predictions, labels)
    auc <- unlist(performance(pred, "auc")@y.values)
    perf <- performance(pred, measure = "tpr", x.measure = "fpr") 
    FPR <- list()
    TPR <- list()
    for(i in 1:length(filenames)){
        FPR[[i]] <- perf@x.values[[i]]
        TPR[[i]] <- perf@y.values[[i]]   
    }
    
    list(auc=auc,FPR=FPR,TPR=TPR)
}

new_auc_7_29 <- function(filenames,Tset){
    
    #  flag==1: only rank information is used
    #  flag==2: FDR information is used
    
    library(ROCR)
    library(pROC)
    predictions <- list()
    labels <- list()
    
    for(i in 1:length(filenames)){
        filename <- filenames[i]
        if(grepl(".csv",filename)){ result <- read.csv(filename); 
        }else{result <- read.table(filename);}
        
        result <- result[!duplicated(result[,1]),]
        
        if(i==1){
            predictions[[i]] <- 1 - result[,"qvalue.dn"]
        }else if(i==2){
            result[is.na(result[,"FDR"]),"FDR"] <- 1
            result[is.na(result[,"Stage2_posterior"]),"Stage2_posterior"] <- 1
            predictions[[i]] <- 1 - result[,"Stage2_posterior"]
        }else{
            predictions[[i]] <- result[,2]/max(result[,2])
        }
        
        labels[[i]] <- rep(-1,length(predictions[[i]]))
        labels[[i]][result[,1] %in% Tset] <- 1
    }
    
    pred <- prediction(predictions, labels)
    auc <- unlist(performance(pred, "auc")@y.values)

    auc
}

plot_err <- function(x,y,sd,flag=1,main="ROC",xlab="Samples",ylab="AUC",xlim=c(0,1),ylim=c(0,1),cols){
    
    if(flag==1){
        par(mar=c(6,6,3,2))
        plot(x, y,xaxt="n", type="b",col=cols[flag], main=main, ylab=ylab,xlab=xlab,xlim=xlim,ylim=ylim,cex.lab=2,cex.main=2,lwd=2,cex.axis=1.8)
        axis(1,at=seq(0.2,0.9,0.1)[1:length(x)],labels=floor(5542*seq(0.2,0.9,0.1))[1:length(x)],font=1,cex.axis=1.5)
    }else{
        lines(x, y, type="b",col=cols[flag],xaxt="n",lwd=2)
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
        
        subfilename <- basename(filename)
        
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
