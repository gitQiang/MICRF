BestP_8_17 <- function(k){
    ##source("BestP_8_17.R")
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
    
    PerM <- array(0,c(12,2))
    
    k1 <- k %% length(netflags)
    if(k1==0) k1=8
    netflag=netflags[k1]
    meth <- floor((k-1)/8) + 1

        for(beta in betaV){
            netsub <- which(netflags==netflag)
            ## recurrent
            path <- paste("result/randresult4_",meth,"/",netstr[netsub],sep="")
            ReM <- matrix(0,20,nodeN[netsub])
            ReM <- randset4_one(netflag,meth,beta,path,ReM)
            
            ## leave-one rank
            path <- paste("result/leaveone4result_",meth,"/",netstr[netsub],sep="")
            HiV <- matrix(0,1,nodeN[netsub])
            HiV <- leaveone4_one(netflag,meth,beta,Tset,path,HiV)
            
            a1 <- sum(HiV/(1:nodeN[netsub]))
            tmp <- sapply(1:20,function(ki) sum(ReM[ki,]/(1:nodeN[netsub])))
            a2 <- median(tmp)
            
            PerM[which(betaV==beta),] <- c(a1,a2)
        }

    save(PerM,file=paste("PerM_8_19",k,sep="_"))
    
}

dealBestP <- function(){
    source("BestP_8_17.R")
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
#     
#     BestP <- matrix(0,8,3)
#     for(i in 1:8){
#         tmp1 <- c()
#         tmp2 <- c()
#         for(meth in 1:4){
#             k= (meth-1)*8+i
#             load(paste("PerM_8_19",k,sep="_"))
#             tmp1 <- cbind(tmp1,PerM[,1])
#             tmp2 <- cbind(tmp2,PerM[,2])
#         }
#         
#         ##oneM <- matrix(rank(tmp1),ncol=4) + matrix(rank(tmp2),ncol=4)
#         tmp1 <- (tmp1 - min(tmp1))/(max(tmp1)-min(tmp1));
#         tmp2 <- (tmp2 - min(tmp2))/(max(tmp2)-min(tmp2));
#         ##oneM <- tmp1 + tmp2
#         BestP[i,1] <- netflags[i]
#         tmp <- which(oneM==max(oneM),arr.ind=TRUE)
#         BestP[i,2] <- tmp[1,2]
#         BestP[i,3] <- betaV[tmp[1,1]]
#     }
#     save(BestP,file="BestP_8_19")
    
    netflags <- c(31,6,7,20,21,26,27,28)
    #betaV <- c(0,0.2,0.4,0.5,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0)
    betaV <- c(0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0)
    BestP <- matrix(0,8,3)
    meth = 4 
    for(i in 1:8){       
        k= (meth-1)*8+i
        load(paste("PerM_8_19",k,sep="_"))
        ##oneM <- rank(PerM[,1]) + rank(PerM[,2])
        ##oneM <- (PerM[,1]-min(PerM[,1]))/(max(PerM[,1])-min(PerM[,1])) + (PerM[,2]-min(PerM[,2]))/(max(PerM[,2])-min(PerM[,2]))
        tmp <- cbind( (PerM[,1]-min(PerM[,1]))/(max(PerM[,1])-min(PerM[,1])) , (PerM[,2]-min(PerM[,2]))/(max(PerM[,2])-min(PerM[,2])) )
        ##tmp <- cbind( rank(PerM[,1]) , rank(PerM[,2]))
        oneM <- apply(tmp,1,min)
        BestP[i,1] <- netflags[i]
        BestP[i,2] <- 4
        BestP[i,3] <- betaV[which.max(oneM[-4])[1]]
    }

    save(BestP,file="BestP_8_19_4")
      
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
        
    result <- read.table(filename)
    
    inrtmp <- rep(0, dim(result)[1])
    subs <- result[,1] %in% oneresult[,1]
    genes <- result[subs,1]
    inrtmp[subs] <- rowSums(oneresult[match(genes,oneresult[,1]),c("dn.LoF","dn.mis3")])>0
    
    tmp <- inrtmp
    #tmp <- rep(0,length(inrtmp))
    #for(m in 1:dim(result)[1]) tmp[m] <- sum(inrtmp[1:m]) /m 
    
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
    
    tmp <- HiV
    #tmp <- rep(0,length(HiV))
    #for(m in 1:length(HiV)) tmp[m] <- sum(HiV[1:m]) /m
    
    tmp
}

AUC_faster <- function(j){
    
    source("enrichana_6_12.R")
    
    load("BestP_8_19_4")
    betaV <- BestP[,3]
    TADAFile="../TADA_DAWN/result/TADAdenovo_meta_dmis.csv"
    Tset <- allPf(TPcut=0.1,TADAFile)
    
    DAWNnames <-  readLines(con <- file("../TADA_DAWN/DAWN_package/TADAdenovo_randset_1.txt","r"))   
    close(con)
    Maginames <-  readLines(con <- file("../MAGI_V0.1_Alpha/mydata/randset_1/MAGIfilenames.txt","r"))   
    close(con)
    
    nnet <- 4
    AUC0 <- array(0,dim=c(20,3+nnet*2))     
    
    DAWNflag=1
    subs <- c(1,2,3,4,5,7,10,12,13,15,18)
    
        for(i in 1:20){   
            filenames <- filenamesf(j,i,betaV=betaV,DAWNnames,Maginames)
            filenames <- filenames[subs]
            
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
            AUC0[i,] <- tmp$auc
            
        }
    
    ##save(AUC0,file=paste("AUC0",j,sep="_"))
    save(AUC0,file=paste("AUC_1",j,sep="_"))
    
}

dealAUC <- function(){
    
    ### do not consider the DAWN predicted results
    
    AUC <- array(0,dim=c(8,20,11)) 
    for(j in 2:9){
        load(paste("AUC0",j,sep="_"))
        AUC[j-1,,] <- AUC0
    }
    
    save(AUC,file="AUC_m_8_19_4")
    
    
    ### overlapped with DAWN predicted results
    AUC <- array(0,dim=c(8,20,11)) 
    for(j in 2:9){
        load(paste("AUC_1",j,sep="_"))
        AUC[j-1,,] <- AUC0
    }
    
    save(AUC,file="AUC_m_8_23_4")   
    
}