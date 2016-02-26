Poisson_test_hq <- function(filename,ntrio){
    tmp <- read.csv(filename)
    pV <- poisson_test(tmp,ntrio)
    tmp <- cbind(tmp,pV)
    tmp[,"min_p"] <- apply(pV,1,min)
    
    subs <- rowSums(tmp[,c("dn.LoF","dn.mis3","dn.mis")]) > 0
    nc <- sum(subs)
    tmp[,"adjust_p"] <- tmp[,"min_p"]
    tmp[subs,"adjust_p"] <- p.adjust(tmp[subs,"min_p"],n=dim(pV)[2]*nc,method="fdr")
    
    #apply(cbind(tmp[,"min_p"] *3 * dim(tmp)[1],1),1,min) ### there should use the number of genes with one mutation at least ????
    ##tmp[,"score"] <- 1- tmp[,"adjust_p"]
    
    tmp
}

meta_poisson_hq <- function(){

    source("ASD_data_set.R")
    mutrate0 <- read.csv("ASD/GeneMutaRatem.csv")
    mutrate0[which(mutrate0[,1]=="KBTBD4"),"dmis"] <- 0.00000127
    mutrate1 <- read.csv("ASD/MutationRatem.csv")
    mutrate0[which(mutrate0[,1]=="HIST1H2AE"),"dmis"] <- mutrate1[which(mutrate1[,1]=="HIST1H2AE"),"mut.rate"]
    mutrate0[which(mutrate0[,1]=="HIST1H2AL"),"dmis"] <- mutrate1[which(mutrate1[,1]=="HIST1H2AL"),"mut.rate"]

    mutrate <- cbind(mutrate0,mutrate1[match(mutrate0[,1],mutrate1[,1]),"mut.rate"])
    colnames(mutrate)[6] <- "mut.rate"
    mutrate <- cbind(mutrate,0,0,0,0)
    colnames(mutrate) <- c("Gene","dmis","mis","LOF","LOFDmis","mut.rate","dn.LoF","dn.mis3","dn.mis","dn.syn")
    
    data <- read.csv("DDD_mutations/datasheet/mutation_4_15_table1.csv",skip=1)
    data1 <- read.csv("DDD_mutations/datasheet/DDDmutation_4_15_table1.csv",skip=1)
    ntrio <- 4409 + 1133
    
    tmpg <- intersect(mutrate[,1],data1[,1])
    subs <- match(tmpg,mutrate[,1])
    mutrate[subs,7:10] <- mutrate[subs,7:10] + data1[match(tmpg,data1[,1]),2:5]
    
    tmpg <- intersect(mutrate[,1],data[,1])
    subs <- match(tmpg,mutrate[,1])
    mutrate[subs,7:10] <- mutrate[subs,7:10] + data[match(tmpg,data[,1]),14:17]
    write.csv(mutrate,file="DDD_mutations/datasheet/Meta_analysis_4_28.csv",row.names=FALSE)
    
    source("Poisson_test_hq.R")
    filename <- "DDD_mutations/datasheet/TADAdenovo_DDD_mis.csv"
    ntrio <- 4409 + 1133
    tmp <- Poisson_test_hq(filename,ntrio)
    tmp <- tmp[order(tmp[,"min_p"]),]
    write.csv(tmp,file="DDD_mutations/datasheet/Ptest_7_30.csv",row.names=FALSE)

    
    source("Poisson_test_hq.R")
    filename <- "ASD/ASDLOF_dMIS_4_16.csv"
    ntrio <- 3953
    tmp <- Poisson_test_hq(filename,ntrio)
    tmp <- tmp[,c(1:13,16)]
    tmp[,"min_p"] <- apply(tmp[,11:13],1,min)
    tmp <- tmp[order(tmp[,"min_p"]),]
    write.csv(tmp,file="ASD/Ptest_4_29.csv",row.names=FALSE)
}

batch_poisson_hq <- function(){
    source("Poisson_test_hq.R")
    filename <- "ASD/ASDLOF_dMIS_4_16.csv"
    ntrio <- 3953
    tmp <- Poisson_test_hq(filename,ntrio)
    write.csv(tmp,file="ASD/TADAresult/randset4/PtestASDall.csv",row.names=FALSE)
    TADAinput2(tmp,"ASD_4_16",dirstr="result/randset4/")
    
    source("ASD_data_set.R")    
    n <- 3953
    for(j in 3){
        for(i in 1:20){
            filename <- paste("ASD/randset4/ASD1sub",j,"depart",i,".csv",sep="")
            ntrio <- round(j*n/10)
            tmp <- Poisson_test_hq(filename,ntrio)
            write.csv(tmp,file=paste("ASD/TADAresult/randset4/Ptest_",j,"depart",i,".csv",sep=""),row.names=FALSE)
            TADAinput2(tmp,paste("part",j,"_",i,sep=""),dirstr="result/randset4/")
            
            filename <- paste("ASD/randset4/ASD2sub",j,"derest",i,".csv",sep="")
            ntrio <- n - ntrio
            tmp <- Poisson_test_hq(filename,ntrio)
            write.csv(tmp,file=paste("ASD/TADAresult/randset4/Ptest_",j,"derest",i,".csv",sep=""),row.names=FALSE)
            TADAinput2(tmp,paste("rest",j,"_",i,sep=""),dirstr="result/randset4/")            
        }
    }

    filename <- "ASD/TADAresult/randset4/PtestASDall.csv"
    asddata <- read.csv(filename)
    genes <- asddata[asddata[,"adjust_p"]<=0.5,"Gene"]
    
    ntrio <- 3953
    for(i in 1:length(genes)){
        tmp <- asddata
        subs <- which(tmp[,"Gene"]==genes[i])
        data <- tmp[subs[1],1:10,drop=FALSE]
        pV <- poisson_test(data,ntrio)
        tmp[subs,colnames(pV)] <- pV
        tmp[subs,"min_p"] <- min(pV)
        
        subt <- rowSums(tmp[,c("dn.LoF","dn.mis3","dn.mis")]) > 0
        nc <- sum(subt)
        tmp[subs,"adjust_p"] <- p.adjust(tmp[subs,"min_p"],n=dim(pV)[2]*nc,method="fdr")
        if(tmp[subs,"adjust_p"] <1){
            print(genes[i])
        }
        write.csv(tmp,file=paste("ASD/TADAresult/leaveone4/Rand2_",i,".csv",sep=""),row.names=FALSE)
        TADAinput2(tmp,paste("rand2",i,sep=""),dirstr="result/leaveone4/")        
    }
    
    
}

poisson_test <- function(data,N){
    
    n <- dim(data)[1]
    data[is.na(data[,"mut.rate"]),"mut.rate"] <- 1
    
    logpM <- sapply(1:n, function(i){
        a1 <- poisson.test(x=data[i,"dn.LoF"], T = N *2* data[i,"LOF"], alternative = "greater", conf.level = 0.95)$p.value
        a2 <- poisson.test(x=data[i,"dn.mis3"] + data[i,"dn.mis"], T = N *2* data[i,"mis"], alternative = "greater", conf.level = 0.95)$p.value
        a3 <- poisson.test(x=data[i,"dn.LoF"]+ data[i,"dn.mis3"] +  data[i,"dn.mis"], T = N *2* (data[i,"LOF"]+ data[i,"mis"]), alternative = "greater", conf.level = 0.95)$p.value
        #a4 <- poisson.test(x=data[i,"dn.mis"]+ data[i,"dn.mis3"], T = N *2* data[i,"mis"], alternative = "greater", conf.level = 0.95)$p.value #!!!!
        #a5 <- poisson.test(x=data[i,"dn.LoF"]+ data[i,"dn.mis3"]+ data[i,"dn.mis"], T = N *2* data[i,"mut.rate"], alternative = "greater", conf.level = 0.95)$p.value
        c(a1,a2,a3)
        #c(a1,a2,a3,a4,a5)
    })
    logpM <- t(logpM)
    #colnames(logpM) <- c("p_LOF","p_dmis","p_both","p_mis","p_all")
    colnames(logpM) <- c("p_LOF","p_mis","p_both")
    
    logpM
}

TADAinput2 <- function(geneinfo,strname,dirstr="random_samples/"){
    
    mapT <- as.matrix(read.delim("Fmap0121.txt",header=FALSE,sep="\t"))
    infofile <- paste(dirstr,"hotnet_input",strname,".txt",sep="")
    write.table(geneinfo[,c("Gene","score")],file=infofile,quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")
        
    genelist <- geneinfo[,1]
    genelist <- mapT[match(genelist,mapT[,2]),1]
    subs <- !is.na(genelist)
    infofile <- paste(dirstr,"CRF_input",strname,".txt",sep="")
    write.table(cbind(genelist[subs],geneinfo[subs,"score"]),file=infofile,quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")
    
}

addtoCRF <- function(){
    CRFr <- read.table("result/control_6_12/v4/coexp/CRFresult_0.2Meta_7_3LBP_7.txt")
    colnames(CRFr) <- c("Gene","Probability","Alternative probability","Pvalue_zscore","Adjust_p")
    posi <- read.csv("DDD_mutations/datasheet/Ptest_7_01.csv")
    tadar <- read.csv("ASD/TADAresult/TADAdenovo_Meta.csv")
    
    #tg <- intersect(posi[,1],CRF[,1])
    result <- cbind(CRFr,posi[match(CRFr[,1],posi[,1]),2:14])
    result <- cbind(result,1,1)
    subs <- match(intersect(CRFr[,1],tadar[,1]),result[,1]) 
    result[subs,19] <- tadar[ match(intersect(CRFr[,1],tadar[,1]),tadar[,1]) ,"pval.TADA.dn"]
    result[subs,20] <- tadar[ match(intersect(CRFr[,1],tadar[,1]),tadar[,1]) ,"qvalue.dn"]
    colnames(result)[19:20] <- c("TADA_pvalue","TADA_qvalue")
    
    result[is.na(result[,"min_p"]),"min_p"] <- 1
    result[is.na(result[,"TADA_pvalue"]),"TADA_pvalue"] <- 1
    result[is.na(result[,"TADA_qvalue"]),"TADA_qvalue"] <- 1
    
    write.csv(result,file="DDD_mutations/datasheet/Meta-analysis_result_7_3.csv",row.names=FALSE)
}