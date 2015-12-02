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
        Tset <- intersect(Tset,result[[i]][,1])
        
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

new_auc_11_9 <- function(filenames,Tset,gf){
    
    library(ROCR)
    predictions <- list()
    labels <- list()
    
    if(gf){ tmp <- read.csv(filenames[1])
        genes <- tmp[rowSums(tmp[,c("dn.LoF","dn.mis3")]) > 0,1]
    }  
    
    for(i in 1:length(filenames)){
        filename <- filenames[i]
        if(grepl(".csv",filename)){ result <- read.csv(filename); 
        }else{result <- read.table(filename);}
        
        if(gf) result <- result[result[,1] %in% genes,]
        
        if(i==1){
            predictions[[i]] <-  -1* result[,"qvalue.dn"]
        }else if(i==2){
            #result <- result[!is.na(result[,"Stage2_posterior"]), ]
            subs <- is.na(result[,"Stage2_posterior"])
            result[subs,"Stage2_posterior"] <- result[subs, "Stage1_posterior"]
            predictions[[i]] <-  -1* result[,"Stage2_posterior"]
        }else{
            predictions[[i]] <- result[,2]
        }
        
        allP <- intersect(Tset,result[,1])
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

allPf <- function(TPcut,TADAFile="ASD/TADAresult/TADAdenovo_ASD4_16.csv"){
    
    TADAall <- read.csv(TADAFile)
    TADAall[is.na(TADAall[,"qvalue.dn"]),"qvalue.dn"] <- 1
    if(TPcut<=1){
        Tset <- TADAall[TADAall[,"qvalue.dn"]< TPcut,1]
    }else{
        Tset <- TADAall[1:TPcut,1]
    }
    
    Tset
}

onefnew <- function(j,i,DAWNnames,Maginames){
    
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
    
#     nnet <- length(netstr)
#     hotnetfile <- 1:nnet
#     misfile <- 1:nnet
#     for(k in 1:nnet){
#         hotnetfile[k] <- paste(resultstr,netstr[k],"hotnetresult1",tmpstr,netnum[k],".txt",sep="")
#         if(k==1){
#             misfile[k] <- paste("result/randresult_5/",netstr[k],"CRFresult_",3,tmpstr,"LBP_",netnum[k],".txt",sep="")
#         }else{
#             misfile[k] <- paste("result/randresult_5/",netstr[k],"CRFresult_",netnum[k],tmpstr,"LBP_",netnum[k],".txt",sep="")
#         }
#     }
#     
#     filenames <- c(TADAfile,DAWNfile,MAGIfile,hotnetfile,misfile)

    misfile <- paste("result/randresult_5/MICRFresult_",2:6,"_",j,"_",i,".txt",sep="")
    filenames <- c(TADAfile,DAWNfile,MAGIfile,misfile)

    filenames
    
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
        
        #result <- result[!duplicated(result[,1]),]
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
