batch_randset4_6_2 <- function(){
    source("BestP_6_13.R")
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
                ReM <- randset4_one(netflag,meth,beta,path,Tset)
                ReV <- colMeans(ReM)
                
                ## leave-one rank
                path <- paste("result/leaveone4result_",meth,"/",netstr[netsub],sep="")
                HiV <- matrix(0,1,nodeN[netsub])
                HiV <- leaveone4_one(netflag,meth,beta,Tset,path,HiV)
                
                tmp <- ReV*HiV
                
                PerM[which(netflag==netflags),meth,which(betaV==beta)] <- max(tmp[10:length(tmp)])
                #ReVL[which(netflag==netflags),meth,which(betaV==beta)] <- max(ReV*HiV)
            }
        }
    } 
    save(PerM,file="PerM_6_13")
    
    BestP <- matrix(0,5,3)
    for(i in 1:5){
        BestP[i,1] <- netflags[i]
        tmp <- which(PerM[i,,]==max(PerM[i,,]),arr.ind=TRUE)
        BestP[i,2] <- tmp[1]
        BestP[i,3] <- betaV[tmp[2]]
    }
    BestP <- BestP[c(1,2,4,3,5),]
    save(BestP,file="BestP_6_13")
}

randset4_one <- function(netflag,meth,beta,path,Tset){
    
    ReM <- matrix(0,20,nodeN[netsub])
    for(i in 1:20){
        filename <- paste(path,"CRFresult_",beta,"part3_",i,"LBP_",netflag,".txt",sep="")
        ReM[i,] <- one_recurrent(filename,i,Tset)
    }
    
    ReM
    
}

one_recurrent <- function(filename,i,Tset){
    
    options(stringsAsFactors=FALSE)
    j=3
    onefile <- paste("ASD/TADAresult/randset4/TADAdenovo_ASD2sub",j,"derest",i,".csv",sep="")
    oneresult <- read.csv(onefile)
    
    result <- read.table(filename)
    
    inrtmp <- rep(0, dim(result)[1])
    subs <- result[,1] %in% oneresult[,1]
    genes <- intersect(result[subs,1],Tset)
    inrtmp[match(genes,result[,1])] <- rowSums(oneresult[match(genes,oneresult[,1]),c("dn.LoF","dn.mis3")])>0
    
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
