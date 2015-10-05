prenode <- function(){
options(stringsAsFactors=FALSE)
source("Network_analysis.R") # mapping_to
source("prenode.R")
# mutrate <- read.csv("TADA_lofmis1202.csv")[c("Gene","mut.rate")]
# mutrate[,1] <- mapping_to(mutrate[,1])
mutrate <- 0

inherifile <- "harvard_PCGC_inherited1017.csv"
denovofile <- "harvard_PCGC_denovo1030.csv"
denovo <- as.matrix(read.csv(denovofile))
denovo[is.na(denovo[,"VarClass"]),"VarClass"] <- "frameshiftdeletion"
inheri <- as.matrix(read.csv(inherifile))
# AFa <- sapply(cases[,"AF"], function(x){max(as.numeric(unlist(strsplit(x,","))))})
# subs <- AFa <= 0.05 # the individual allele frequency in cases  5%
# cases[,"AF"] <- AFa
# cases <- cases[subs,]
                                     
header <- colnames(denovo)
denovo <- cbind(denovo,"denovo")
inheri <- cbind(inheri,"inherited")
samples <- rbind(denovo,inheri)
colnames(samples) <- c(header,"DENOVO")

#subs <- !(is.na(samples[,"ESPfreq"]) & is.na(samples[,"X1KGfreq"]) & is.na(samples[,"esp.MAF"]))
#samples <- samples[subs,]
# samples[is.na(samples[,"ESPfreq"]),"ESPfreq"] <- 0
# samples[is.na(samples[,"X1KGfreq"]),"X1KGfreq"] <- 0
# samples[is.na(samples[,"esp.MAF"]),"esp.MAF"] <- 0
# tmp <- sapply(samples[,"esp.MAF"], function(x){max(as.numeric(unlist(strsplit(x,","))))})
# tmp <- tmp/100
# PF <- apply(as.matrix(cbind(tmp,cbind(as.numeric(samples[,"ESPfreq"]),as.numeric(samples[,"X1KGfreq"])))),1,max)
# samples <- samples[PF <= 0.001,] #0.1%

trun <- c("frameshiftdeletion","frameshiftinsertion","stopgainSNV","stoplossSNV","frameshiftsubstitution")
mis3 <- c("nonframeshiftdeletion","nonframeshiftinsertion","nonframeshiftsubstitution","nonsynonymousSNV")
samples <- samples[samples[,"VarClass"] %in% union(trun,mis3),]

strscore <- c("SIFTscr","PP2scr","MutTscr","GERP")
for(i in 1:length(strscore)){
    samples[,strscore[i]] <- change(samples[,strscore[i]],-2)
}
tmp <- sapply(1:dim(samples)[1],function(i) mean(as.numeric( samples[i, strscore[as.numeric(samples[i,strscore])>0] ] )))
samples <- cbind(samples,tmp)
colnames(samples)[dim(samples)[2]] <- "SCORE"
score <- tmp
samples <- as.matrix(samples)

tmp <- samples[,"GeneName"]
tmp <- sapply(tmp, function(x){unlist(strsplit(x,","))[1]})
tmp <- sapply(tmp, function(x){unlist(strsplit(x,"\\("))[1]})
genes <- unique(tmp)
gmap <- cbind(genes,genes)
gmap[,2] <- mapping_to(gmap[,2])
samples[,"GeneName"] <- gmap[match(tmp,gmap[,1]),2]
genes <- unique(samples[,"GeneName"])

# deal with NA
label <- 1:dim(samples)[1]
label[(samples[,"DENOVO"]=="denovo") & (samples[,"VarClass"] %in% trun)] <- 1
label[(samples[,"DENOVO"]=="denovo") & (samples[,"VarClass"] %in% mis3)] <- 2
label[(samples[,"DENOVO"]=="inherited") & (samples[,"VarClass"] %in% trun)] <- 3
label[(samples[,"DENOVO"]=="inherited") & (samples[,"VarClass"] %in% mis3)] <- 4
subs <- !is.na(score)
a <- 1:4
tmp <- (samples[,"DENOVO"]=="denovo") & (samples[,"VarClass"] %in% trun) & subs
a[1] <- mean(as.numeric(score[tmp]))
tmp <- (samples[,"DENOVO"]=="denovo") & (samples[,"VarClass"] %in% mis3) & subs
a[2] <- mean(as.numeric(score[tmp]))
tmp <- (samples[,"DENOVO"]=="inherited") & (samples[,"VarClass"] %in% trun) & subs
a[3] <- mean(as.numeric(score[tmp]))
tmp <- (samples[,"DENOVO"]=="inherited") & (samples[,"VarClass"] %in% mis3) & subs
a[4] <- mean(as.numeric(score[tmp]))
tmp <- dealNA(samples,label,trun,mis3,a,score)
samples[,"SCORE"] <- tmp
write.csv(samples,file="PCGC_node_SCORE.csv",row.names=FALSE)

prenode0(samples,trun,mis3,genes,mutrate)

}

prenode1 <- function(){
    options(stringsAsFactors=FALSE)
    source("Network_analysis.R") # mapping_to
    source("prenode.R")
    mutrate <- 0
    
    inherifile <- "PCGCinherited.csv"
    denovofile <- "PCGCdenovo.csv"
    inheri <- as.matrix(read.csv(inherifile))
    denovo <- as.matrix(read.csv(denovofile))  # exclude "Polyphen2_HVAR_score"
    
    varfunc <- c("exonic","exonic;exonic","exonic;splicing","intronic","ncRNA_exonic","splicing")
    varmap <- c("exonic","exonic","exonic,splicing","intronic","ncRNA_exonic","splicing")
    denovo[,"Func.refgene"] <- varmap[match(denovo[,"Func.refgene"],varfunc)]
    denovo <- denovo[denovo[,"Func.refgene"] %in% c("exonic","exonic","exonic,splicing","splicing"),]
    
    # make the same denovo and inheri col names
    mutype1 <- c("","frameshift deletion","frameshift insertion","nonframeshift deletion","nonframeshift insertion","nonsynonymous SNV","nonsynonymous SNV;nonsynonymous SNV", "stopgain","stoploss","synonymous SNV","unknown") 
    mutmap <- c("frameshiftdeletion","frameshiftdeletion","frameshiftinsertion","nonframeshiftdeletion","nonframeshiftinsertion","nonsynonymousSNV","nonsynonymousSNV", "stopgainSNV","stoplossSNV","synonymousSNV","unknown") ##!!!!!!! 20 splicing mutations are seen as LOF in denovo.
    denovo[,"ExonicFunc.refgene"] <- mutmap[match(denovo[,"ExonicFunc.refgene"],mutype1)]
    denovo[denovo[,"SIFT_score"]==".","SIFT_score"] = ""
    denovo[denovo[,"Polyphen2_HDIV_score"]==".","Polyphen2_HDIV_score"] = ""
    denovo[denovo[,"MutationTaster_score"]==".","MutationTaster_score"] = ""
    denovo[denovo[,"GERP.._RS"]==".","GERP.._RS"] = ""
    denovo[denovo[,"SIFT_score"]=="","SIFT_score"] = NA
    denovo[denovo[,"Polyphen2_HDIV_score"]=="","Polyphen2_HDIV_score"] = NA
    denovo[denovo[,"MutationTaster_score"]=="","MutationTaster_score"] = NA
    denovo[denovo[,"GERP.._RS"]=="","GERP.._RS"] = NA
    
    title1 <- c("CHR","POS","REF","ALT","Subject","parent_min","esp6500si_all","Gene","Func.refgene","ExonicFunc.refgene","AAChange.refgene","X1000g2014oct_all","SIFT_score","Polyphen2_HDIV_score","Polyphen2_HDIV_pred","MutationTaster_score","GERP.._RS")            
    title2 <- c("CHROM","POS","REF","ALT","proband.ID.GT.AD.DP.GQ.PL.","Inherited.from","esp.MAF","GeneName","VarFunc","VarClass","AAChange","X1KGfreq","SIFTscr","PP2scr","PP2prd","MutTscr","GERP") 
    denovo <- denovo[,title1]
    inheri <- inheri[,title2]        
    header <- title2
    denovo <- cbind(denovo,"denovo")
    inheri <- cbind(inheri,"inherited")
    samples <- rbind(denovo,inheri)
    colnames(samples) <- c(header,"DENOVO")
    
    trun <- c("frameshiftdeletion","frameshiftinsertion","stopgainSNV","stoplossSNV","frameshiftsubstitution")
    mis3 <- c("nonframeshiftdeletion","nonframeshiftinsertion","nonframeshiftsubstitution","nonsynonymousSNV")
    samples <- samples[samples[,"VarClass"] %in% union(trun,mis3),]
    
    subs1 <- samples[,"VarClass"] %in% trun
    subs2 <- samples[,"VarClass"] %in% mis3 &  samples[,"PP2prd"] %in% c("D","P")
    samples <- samples[subs1 | subs2,]
    
    strscore <- c("SIFTscr","PP2scr","MutTscr","GERP")
    for(i in 1:length(strscore)){
        samples[,strscore[i]] <- change(samples[,strscore[i]],-2)
    }
    tmp <- sapply(1:dim(samples)[1],function(i) mean(as.numeric( samples[i, strscore[as.numeric(samples[i,strscore])>0] ] )))
    samples <- cbind(samples,tmp)
    colnames(samples)[dim(samples)[2]] <- "SCORE"
    score <- tmp
    samples <- as.matrix(samples)
    
    tmp <- samples[,"GeneName"]
    tmp <- sapply(tmp, function(x){unlist(strsplit(x,","))[1]})
    tmp <- sapply(tmp, function(x){unlist(strsplit(x,"\\("))[1]})
    genes <- unique(tmp)
    gmap <- cbind(genes,genes)
    gmap[,2] <- mapping_to(gmap[,2])
    samples[,"GeneName"] <- gmap[match(tmp,gmap[,1]),2]
    genes <- as.vector(unique(samples[,"GeneName"]))
    
    # deal with NA
    label <- 1:dim(samples)[1]
    label[(samples[,"DENOVO"]=="denovo") & (samples[,"VarClass"] %in% trun)] <- 1
    label[(samples[,"DENOVO"]=="denovo") & (samples[,"VarClass"] %in% mis3)] <- 2
    label[(samples[,"DENOVO"]=="inherited") & (samples[,"VarClass"] %in% trun)] <- 3
    label[(samples[,"DENOVO"]=="inherited") & (samples[,"VarClass"] %in% mis3)] <- 4
    subs <- !is.na(score)
    a <- 1:4
    tmp <- (samples[,"DENOVO"]=="denovo") & (samples[,"VarClass"] %in% trun) & subs
    a[1] <- mean(as.numeric(score[tmp]))
    tmp <- (samples[,"DENOVO"]=="denovo") & (samples[,"VarClass"] %in% mis3) & subs
    a[2] <- mean(as.numeric(score[tmp]))
    tmp <- (samples[,"DENOVO"]=="inherited") & (samples[,"VarClass"] %in% trun) & subs
    a[3] <- mean(as.numeric(score[tmp]))
    tmp <- (samples[,"DENOVO"]=="inherited") & (samples[,"VarClass"] %in% mis3) & subs
    a[4] <- mean(as.numeric(score[tmp]))
    tmp <- dealNA(samples,label,trun,mis3,a,score)
    samples[,"SCORE"] <- tmp
    
    write.csv(samples,file="PCGC1_node_SCORE.csv",row.names=FALSE)
    
    samples <- read.csv("PCGC1_node_SCORE.csv")
    prenode0(samples,trun,mis3,genes,mutrate,str="New",flag=2)
    
}

prenode2 <- function(filenames,flag=1,str="PCGC"){
    options(stringsAsFactors=FALSE)
    source("Network_analysis.R") # mapping_to
    source("prenode.R")
    mutrate <- 0
    trun <- c("frameshiftdeletion","frameshiftinsertion","stopgainSNV","stoplossSNV","frameshiftsubstitution")
    mis3 <- c("nonframeshiftdeletion","nonframeshiftinsertion","nonframeshiftsubstitution","nonsynonymousSNV")
    #filenames <- c("cases0120.csv","PCGCinherited.csv")
    datalist <- predeal(trun,mis3,filenames)
    denovo <- datalist[[1]]
    inheri <- datalist[[2]]        
    denovo <- cbind(denovo,"denovo")
    inheri <- cbind(inheri,"inherited")
    samples <- rbind(denovo,inheri)
    colnames(samples)[dim(samples)[2]] <- "DENOVO"
    #str <- "PCGC"
    strscore <- c("SIFTscr","PP2scr","MutTscr","GERP")
    prestep2(strscore,samples,trun,mis3,mutrate=0,str,flag)

}

filter2 <- function(){
    source("prenode.R")
    
    # all samples 
    filenames <- c("data/cases0120.csv","data/PCGCinherited.csv")
    flag <- 2
    prenode2(filenames,flag)
    
    
    strss <- c("cade","cain")
    flag <- 1
    # NDD, non-NDD and others
    strV <- c("random_samples/NDD_","random_samples/nonNDD_","random_samples/other_")
    for(i in 1:3){
        filenames <- paste(strV[i],strss,".csv",sep="")
        prenode2(filenames,flag,strV[i])    
    }
    # CHD naturre and Others
    strV <- c("random_samples/CHD_","random_samples/CHDo_")
    for(i in 1:2){
        filenames <- paste(strV[i],strss,".csv",sep="")
        prenode2(filenames,flag,strV[i])    
    }
    # random 2 and 3 
    strV <- c("random_samples/T2","random_samples/T3","random_samples/P2","random_samples/P3")
    for(i in 1:4){
        filenames <- paste(strV[i],strss,".csv",sep="")
        prenode2(filenames,flag,strV[i])    
    }    
    
}

prenode0 <- function(samples,trun,mis3,genes,mutrate,str="",flag=1){

    tablelist <- list()
    tablelist[[1]] <- samples[(samples[,"DENOVO"]=="denovo") & (samples[,"VarClass"] %in% trun) ,]
    tablelist[[2]] <- samples[(samples[,"DENOVO"]=="denovo") & (samples[,"VarClass"] %in% mis3) ,]
    tablelist[[3]] <- samples[(samples[,"DENOVO"]=="inherited") & (samples[,"VarClass"] %in% trun) ,]
    tablelist[[4]] <- samples[(samples[,"DENOVO"]=="inherited") & (samples[,"VarClass"] %in% mis3) ,]
    
    if(flag==2){
        n.de <- dim(samples[ samples[,"DENOVO"]=="denovo" ,])[1]
        n.in <- dim(samples[ samples[,"DENOVO"]=="inherited" ,])[1]
        n.tru <- sum(samples[,"VarClass"] %in% trun)
        n.mis <- sum(samples[,"VarClass"] %in% mis3)
        W1 <- 1:2;W2=W1;W3=W1;
        W1[1] <- (length(unique(samples[ samples[,"DENOVO"]=="denovo" ,"GeneName"]))/n.de) / (length(unique(samples[ samples[,"DENOVO"]=="inherited" , "GeneName"]))/n.in); 
        W1[2] <- 1;
        W2[1] <- (length(unique(samples[samples[,"VarClass"] %in% trun, "GeneName"]))/n.tru) / (length(unique(samples[samples[,"VarClass"] %in% mis3, "GeneName"]))/n.mis); 
        W2[2] <- 1; 
        W3 <- W2;
    }else if(flag==1){
        #W1 <- c(6.68,1)
        #W2 <- c(3.97,1)
        #W3 <- c(3.97,1)
        W1 <- c(4.544202, 1)
        W2 <- c(2.80892, 1)
        W3 <- c(2.80892, 1)
    }
    
    print(W1)
    print(W2)
    print(W3)
    
    sM <- batchscoref(genes,tablelist,mutrate,W1,W2,W3)
    result <- cbind(genes,sM)
    result <- result[order(as.numeric(result[,2]),decreasing=TRUE),]
    
    mapT <- as.matrix(read.delim("Fmap0121.txt",header=FALSE,sep="\t"))
    result1 <- result
    result1[,1] <- mapT[match(result[,1],mapT[,2]),1]
    result1 <- result1[!is.na(result1[,1]),]
    
    write.table(result,file=paste(str,"Pheno_score.txt",sep=""),quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")
    write.table(result1,file=paste(str,"Pheno_score_map.txt",sep=""),quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")

}

batchscoref <- function(genes,tablelist,mutrate,W1,W2,W3){
    tmp <- sapply(1:length(genes), scoref,genes,tablelist)
    tmp <- t(tmp)
    tmp <- tmp/matrix(apply(tmp,2,max),dim(tmp)[1],dim(tmp)[2],byrow=TRUE)
    tmp1 <- tmp[,1:2] %*% W2
    tmp2 <- tmp[,3:4] %*% W3
    #tmp <- cbind(tmp1/max(tmp1),tmp2/max(tmp2)) %*% W1
    tmp <- cbind(tmp1/sum(W2),tmp2/sum(W3)) %*% W1
    
    #tmp <- tmp/max(tmp)
    #tmp <- change(tmp,-2)
    
#     subs <- tmp > 0.5
#     tmp1 <- tmp[subs]
#     tmp1 <- plogis(tmp1, location = 0.5, scale = 1, lower.tail = TRUE, log.p = FALSE)
#     tmp[subs] <- tmp1

    #library(MASS)
    #pa <- fitdistr(tmp,"gamma")
    #tmp <- pgamma(tmp, shape=pa$estimate[1], rate = pa$estimate[2], lower.tail = TRUE, log.p = FALSE)
    #tmp <- pgamma(tmp, shape=0.5, rate = 1, lower.tail = TRUE, log.p = FALSE)

    k=2
    library(nor1mix)
    re <- norMixEM(tmp,k)
    paM <- matrix(re[1:(3*k)],k,3)
    vaM <- matrix(c(dnorm(tmp,mean=paM[1,1],sd=paM[1,2]) * paM[1,3],dnorm(tmp,mean=paM[2,1],sd=paM[2,2]) * paM[2,3]),length(tmp),k)
    tmp <- pnorm(tmp,mean=re[k],sd=re[2*k])*exp(log(vaM[,k])-log(rowSums(vaM)))
    a <- sort(tmp,decreasing=TRUE)
    print(a[1]-a[50])
    print(sum(a>0.5))
    
    tmp
}

change <- function(Aa,flag=1){
    # flag: 2, >=; 1, >; 0, ==; -1, <; -2, <=;
    subs <- is.na(Aa)
    subs1 <- !subs
    A <- as.numeric(Aa[subs1])
    #n <- sum(A>cutt)
    n <- length(A)
    if(flag==2) 	tmp0 <- sapply(A, function(x){sum(A>=x)/n})
    if(flag==1) 	tmp0 <- sapply(A, function(x){sum(A>x)/n})
    if(flag==0) 	tmp0 <- sapply(A, function(x){sum(A==x)/n})
    if(flag==-1) 	tmp0 <- sapply(A, function(x){sum(A<x)/n})
    if(flag==-2) 	tmp0 <- sapply(A, function(x){sum(A<=x)/n})
    
    tmp <- Aa
    tmp[subs] <- -1
    tmp[subs1] <- tmp0
    
    tmp
}

scoref <- function(k,genes,tablelist){
    gene <- genes[k]
    nn <- 1:4
    sn <- 1:4    
    for(i in 1:4){
        tmp <- tablelist[[i]][ tablelist[[i]][,"GeneName"]==gene, "SCORE"]
        nn[i] <- length(tmp)
        sn[i] <- sum(as.numeric(tmp))
    }
    #W1 <- matrix(c(7,1),2,1) # weights for denovo and inherited
    #W2 <- matrix(c(2,1),2,1) # weights for LOF and missence for denovo
    #W3 # weights for LOF and missence for inherited
    #s <- sum(W1 * c(sum(W2*sn[1:2]), sum(W3*sn[3:4])))
    #s <- 1/(nn[1]+nn[2]+1) * (1/(nn[1]+1) * sn[1]  + 1/(nn[2]+1) * sn[2] ) + 1/(nn[3]+nn[4]+1) *(1/(nn[3]+1) * sn[3] + 1/(nn[4]+1) * sn[4])
    #s <- max(s/(sum(nn) * as.numeric(mutrate[mutrate[,1]==gene,2])))
    #s <- s/(sum(nn)
        
    sn
}

dealNA <- function(samples,label,trun,mis3,a,score){
    #score <- samples[,"SCORE"]
    subs <- !is.na(score)
    subs1 <- which(is.na(score))
    
    for(i in 1:length(subs1)){
        if(label[subs1[i]]==1){
        tmp <- (samples[,"GeneName"]==samples[subs1[i],"GeneName"]) & (samples[,"DENOVO"]=="denovo") & (samples[,"VarClass"] %in% trun) & subs
        }else if(label[subs1[i]]==2){
        tmp <- (samples[,"GeneName"]==samples[subs1[i],"GeneName"]) & (samples[,"DENOVO"]=="denovo") & (samples[,"VarClass"] %in% mis3) & subs
        }else if(label[subs1[i]]==3){
        tmp <- (samples[,"GeneName"]==samples[subs1[i],"GeneName"]) & (samples[,"DENOVO"]=="inherited") & (samples[,"VarClass"] %in% trun) & subs
        }else if(label[subs1[i]]==4){
        tmp <- (samples[,"GeneName"]==samples[subs1[i],"GeneName"]) & (samples[,"DENOVO"]=="inherited") & (samples[,"VarClass"] %in% mis3) & subs
        }
        b <- mean(as.numeric(score[tmp]))
        score[subs1[i]] <- ifelse(is.na(b),a[label[subs1[i]]],b)
    }
    score
}

controlM <- function(){
    options(stringsAsFactors=FALSE)
    source("Network_analysis.R") # mapping_to
    mutrate <- 0
    
    inherifile <- "harvard_ssc_inherited1017.csv"
    denovofile <- "harvard_ssc_denovo1030.csv"
    denovo <- as.matrix(read.csv(denovofile))
    inheri <- as.matrix(read.csv(inherifile))

    header <- colnames(denovo)
    denovo <- cbind(denovo,"denovo")
    inheri <- cbind(inheri,"inherited")
    samples <- rbind(denovo,inheri)
    colnames(samples) <- c(header,"DENOVO")
      
    trun <- c("frameshiftdeletion","frameshiftinsertion","stopgainSNV","stoplossSNV","frameshiftsubstitution")
    mis3 <- c("nonframeshiftdeletion","nonframeshiftinsertion","nonframeshiftsubstitution","nonsynonymousSNV")
    samples <- samples[samples[,"VarClass"] %in% union(trun,mis3),]
    
    strscore <- c("SIFTscr","PP2scr","MutTscr","GERP")
    for(i in 1:length(strscore)){
        samples[,strscore[i]] <- change(samples[,strscore[i]],-2)
    }
    tmp <- sapply(1:dim(samples)[1],function(i) mean(as.numeric( samples[i, strscore[as.numeric(samples[i,strscore])>0] ] )))
    tmp[is.na(tmp)] <- median(tmp[!is.na(tmp)])
    samples <- cbind(samples,tmp)
    colnames(samples)[dim(samples)[2]] <- "SCORE"
    samples <- as.matrix(samples)
    
    tmp <- samples[,"GeneName"]
    tmp <- sapply(tmp, function(x){unlist(strsplit(x,","))[1]})
    tmp <- sapply(tmp, function(x){unlist(strsplit(x,"\\("))[1]})
    genes <- unique(tmp)
    gmap <- cbind(genes,genes)
    gmap[,2] <- mapping_to(gmap[,2])
    samples[,"GeneName"] <- gmap[match(tmp,gmap[,1]),2]
    genes <- unique(samples[,"GeneName"])
    
    write.csv(samples,file="PCGC_node_SCORE.csv",row.names=FALSE)
    tablelist <- list()
    tablelist[[1]] <- samples[(samples[,"DENOVO"]=="denovo") & (samples[,"VarClass"] %in% trun) ,]
    tablelist[[2]] <- samples[(samples[,"DENOVO"]=="denovo") & (samples[,"VarClass"] %in% mis3) ,]
    tablelist[[3]] <- samples[(samples[,"DENOVO"]=="inherited") & (samples[,"VarClass"] %in% trun) ,]
    tablelist[[4]] <- samples[(samples[,"DENOVO"]=="inherited") & (samples[,"VarClass"] %in% mis3) ,]
    
    n.de <- dim(samples[ samples[,"DENOVO"]=="denovo" ,])[1]
    n.in <- dim(samples[ samples[,"DENOVO"]=="inherited" ,])[1]
    W1 <- 1:2;W2=W1;W3=W1;
    W1[1] <- (n.in/length(unique(samples[ samples[,"DENOVO"]=="inherited" , "GeneName"]))) / (n.de/length(unique(samples[ samples[,"DENOVO"]=="denovo" ,"GeneName"]))); W1[2] <- 1;
    W2[1] <- (dim(tablelist[[2]])[1]/length(unique(tablelist[[2]][, "GeneName"]))) / (dim(tablelist[[1]])[1]/length(unique(tablelist[[1]][, "GeneName"]))); W2[2] <- 1; 
    W3[1] <- (dim(tablelist[[4]])[1]/length(unique(tablelist[[4]][, "GeneName"]))) / (dim(tablelist[[3]])[1]/length(unique(tablelist[[3]][, "GeneName"]))); W3[2] <- 1; 
       


}

NDD_nonNDD <- function(){

    filename <- "NDD_1121.csv"
    ID <- as.matrix(read.csv(filename))
    ID[is.na(ID[,2]),2] <- 2
    nonNDDid <- ID[as.numeric(ID[,2])==0,1]
    NDDid <- ID[as.numeric(ID[,2])==1,1]
    
    sfile <- "PCGC_node_SCORE.csv"
    samples <- read.csv(sfile)
    n <- dim(samples)[1]
    subs <- sapply(1:n,function(i) {
        b <- unlist(strsplit(gsub("[(.);]","\t",samples[i,"proband.ID.GT.AD.DP.GQ.PL."]),"\t"))
        c(any(b[seq(1,length(b),3)] %in% NDDid),any(b[seq(1,length(b),3)] %in% nonNDDid))
        })
    subs <- t(subs)
    NDDsamples <- samples[subs[,1],]
    nonNDDsamples <- samples[subs[,2],]
    str <- "NDD"
    str1 <- "nonNDD"
    trun <- c("frameshiftdeletion","frameshiftinsertion","stopgainSNV","stoplossSNV","frameshiftsubstitution")
    mis3 <- c("nonframeshiftdeletion","nonframeshiftinsertion","nonframeshiftsubstitution","nonsynonymousSNV")
    mutrate <- 0
    
    flag <- 2
    if(flag==1){
        chdfile <- "pheno6topgenes.txt"
        chdgene <- as.matrix(read.table(chdfile))
        NDDsamples <- NDDsamples[!(NDDsamples[,"GeneName"] %in% chdgene),]
        nonNDDsamples <- nonNDDsamples[!(nonNDDsamples[,"GeneName"] %in% chdgene),]
        str <- "NDD1"
        str1 <- "nonNDD1"
    }else if(flag==2){
        chdfile1 <- "NDDtopgenes.txt"
        chdgene1 <- as.matrix(read.table(chdfile1))
        NDDsamples <- NDDsamples[!(NDDsamples[,"GeneName"] %in% chdgene1),]
        chdfile2 <- "nonNDDtopgenes.txt"
        chdgene2 <- as.matrix(read.table(chdfile2))        
        nonNDDsamples <- nonNDDsamples[!(nonNDDsamples[,"GeneName"] %in% chdgene2),]
        str <- "NDD2"
        str1 <- "nonNDD2"   
    }
    write.csv(NDDsamples,file=paste(str,"_node_SCORE.csv",sep=""),row.names=FALSE)
    write.csv(nonNDDsamples,file=paste(str1,"_node_SCORE.csv",sep=""),row.names=FALSE)
    
    genes <- as.vector(unique(NDDsamples[,"GeneName"]))
    prenode0(NDDsamples,trun,mis3,genes,mutrate,str)
  
    genes <- as.vector(unique(nonNDDsamples[,"GeneName"]))
    prenode0(nonNDDsamples,trun,mis3,genes,mutrate,str1)
}

random_samples <- function(){
    mutrate <- 0
    sfile <- "PCGC1_node_SCORE.csv"
    samples <- read.csv(sfile)
    trun <- c("frameshiftdeletion","frameshiftinsertion","stopgainSNV","stoplossSNV","frameshiftsubstitution")
    mis3 <- c("nonframeshiftdeletion","nonframeshiftinsertion","nonframeshiftsubstitution","nonsynonymousSNV")
    
    tablelist <- list()
    tablelist[[1]] <- samples[(samples[,"DENOVO"]=="denovo") & (samples[,"VarClass"] %in% trun) ,]
    tablelist[[2]] <- samples[(samples[,"DENOVO"]=="denovo") & (samples[,"VarClass"] %in% mis3) ,]
    tablelist[[3]] <- samples[(samples[,"DENOVO"]=="inherited") & (samples[,"VarClass"] %in% trun) ,]
    tablelist[[4]] <- samples[(samples[,"DENOVO"]=="inherited") & (samples[,"VarClass"] %in% mis3) ,]
    
    ## random select mutations
    Kvec <- 2:3
    n.col <- dim(samples)[2]
    for(K in Kvec){
        samrandt <- matrix(,0,n.col)
        samrandp <- matrix(,0,n.col)
        strt <- "random_samples/T"
        strp <- "random_samples/P"
        for(i in 1:4){
            n.mutation <- dim(tablelist[[i]])[1]
            samsplit <- sample.cross(n.mutation,K)
            samrandt <- rbind(samrandt,tablelist[[i]][samsplit$train[[1]],])
            samrandp <- rbind(samrandp,tablelist[[i]][samsplit$pred[[1]],])
        }
        write.csv(samrandt,file=paste(strt,K,"PCGC1_node_SCORE.csv",sep=""),row.names=FALSE)
        write.csv(samrandp,file=paste(strp,K,"PCGC1_node_SCORE.csv",sep=""),row.names=FALSE)
        
        genes <- as.vector(unique(samrandt[,"GeneName"]))
        prenode0(samrandt,trun,mis3,genes,mutrate,paste(strt,K,sep=""))
        genes <- as.vector(unique(samrandp[,"GeneName"]))
        prenode0(samrandp,trun,mis3,genes,mutrate,paste(strp,K,sep=""))
    }
    
    ## random samples: NDD;  non-NDD;  others
    filename <- "NDD_1121.csv"
    ID <- as.matrix(read.csv(filename))
    ID[is.na(ID[,2]),2] <- 2
    nonNDDid <- ID[as.numeric(ID[,2])==0,1]
    NDDid <- ID[as.numeric(ID[,2])==1,1]
    otherid <- ID[as.numeric(ID[,2])==2,1]
    
    n <- dim(samples)[1]
    subs <- sapply(1:n,function(i) {
        b <- unlist(strsplit(gsub("[(.);]","\t",samples[i,"proband.ID.GT.AD.DP.GQ.PL."]),"\t"))
        c(any(b[seq(1,length(b),3)] %in% NDDid),any(b[seq(1,length(b),3)] %in% nonNDDid),any(b[seq(1,length(b),3)] %in% otherid))
    })
    subs <- t(subs)
    K <- 3
    strV <- c("random_samples/NDD","random_samples/nonNDD","random_samples/other")
    for(i in 1:K){
        splitsams <- samples[subs[,i],]
        write.csv(splitsams,file=paste(strV[i],"_node_SCORE.csv",sep=""),row.names=FALSE)
        genes <- as.vector(unique(splitsams[,"GeneName"]))
        prenode0(splitsams,trun,mis3,genes,mutrate,strV[i])
    }
    ## CHD Nature used; others
    CHDid <- as.matrix(read.table("random_samples/CHDnature.txt"))
    CHDid <- unique(CHDid)
    n <- dim(samples)[1]
    subs <- sapply(1:n,function(i) {
        b <- unlist(strsplit(gsub("[(.);]","\t",samples[i,"proband.ID.GT.AD.DP.GQ.PL."]),"\t"))
        any(b[seq(1,length(b),3)] %in% CHDid)
    })
    strV <- c("random_samples/CHDnature","random_samples/CHDother")
    for(i in 1:2){
        if(i==1){splitsams <- samples[subs,];}
        if(i==2){splitsams <- samples[!subs,];}
        write.csv(splitsams,file=paste(strV[i],"_node_SCORE.csv",sep=""),row.names=FALSE)
        genes <- as.vector(unique(splitsams[,"GeneName"]))
        prenode0(splitsams,trun,mis3,genes,mutrate,strV[i])
    }
    
}

random_samples1 <- function(){
    
    filenames <-  c("data/cases0120.csv","data/PCGCinherited.csv","data/PCGCnontransmitted0120.csv") #,"data/controls0120.csv","data/sscinherited0120.csv")
    trun <- c("frameshiftdeletion","frameshiftinsertion","stopgainSNV","stoplossSNV","frameshiftsubstitution")
    mis3 <- c("nonframeshiftdeletion","nonframeshiftinsertion","nonframeshiftsubstitution","nonsynonymousSNV")
    strss <- c("cade","cain","cano")
    
    ## random samples: NDD;  non-NDD;  others
    filename <- "random_samples/NDD_1121.csv"
    ID <- as.matrix(read.csv(filename))
    ID[is.na(ID[,2]),2] <- 2
    nonNDDid <- ID[as.numeric(ID[,2])==0,1]
    NDDid <- ID[as.numeric(ID[,2])==1,1]
    otherid <- ID[as.numeric(ID[,2])==2,1]
    
    for(i in 1:length(filenames)){
        samples <- read.csv(filenames[i])
        n <- dim(samples)[1]
        if("proband.ID.GT.AD.DP.GQ.PL." %in% colnames(samples)){
            subs <- sapply(1:n,function(i) {
                b <- unlist(strsplit(gsub("[();]","\t",samples[i,"proband.ID.GT.AD.DP.GQ.PL."]),"\t"))
                c(any(b[seq(1,length(b),3)] %in% NDDid),any(b[seq(1,length(b),3)] %in% nonNDDid),any(b[seq(1,length(b),3)] %in% otherid))
        })
        }else{
            subs <- matrix(c(samples[,"Subject"] %in% NDDid,samples[,"Subject"] %in% nonNDDid,samples[,"Subject"] %in% otherid),3,n,byrow=TRUE)
        }
        subs <- t(subs)
        K <- 3
        strV <- c("random_samples/NDD_","random_samples/nonNDD_","random_samples/other_")
        for(k in 1:K){
            splitsams <- samples[subs[,k],]
            write.csv(splitsams,file=paste(strV[k],strss[i],".csv",sep=""),row.names=FALSE)    
        }
    }
    
    ## CHD Nature used; others
    CHDid <- as.matrix(read.table("random_samples/CHDnature.txt"))
    CHDid <- unique(CHDid)
    
    for(i in 1:length(filenames)){
        samples <- read.csv(filenames[i])
        n <- dim(samples)[1]
        if("proband.ID.GT.AD.DP.GQ.PL." %in% colnames(samples)){
            subs <- sapply(1:n,function(i) {
                b <- unlist(strsplit(gsub("[(.);]","\t",samples[i,"proband.ID.GT.AD.DP.GQ.PL."]),"\t"))
                any(b[seq(1,length(b),3)] %in% CHDid)
            })
        }else{
            subs <- samples[,"Subject"] %in% CHDid
        }
        
        strV <- c("random_samples/CHD_","random_samples/CHDo_")
        for(k in 1:2){
            if(k==1){splitsams <- samples[subs,];}
            if(k==2){splitsams <- samples[!subs,];}
            write.csv(splitsams,file=paste(strV[k],strss[i],".csv",sep=""),row.names=FALSE)
        }
    }

    ## random select mutations
    strt <- "random_samples/T"
    strp <- "random_samples/P"
    for(i in 1:length(filenames)){
        samples <- read.csv(filenames[i])
        n <- dim(samples)[1]
        k=2
        samsplit <- sample.cross(n,k)
        samrandt <- samples[samsplit$train[[1]],]
        samrandp <- samples[samsplit$pred[[1]],]
        write.csv(samrandt,file=paste(strt,k,strss[i],".csv",sep=""),row.names=FALSE)
        write.csv(samrandp,file=paste(strp,k,strss[i],".csv",sep=""),row.names=FALSE)
        
        k=3
        samsplit <- sample.cross(n,k)
        samrandt <- samples[samsplit$train[[1]],]
        samrandp <- samples[samsplit$pred[[1]],]
        write.csv(samrandt,file=paste(strt,k,strss[i],".csv",sep=""),row.names=FALSE)
        write.csv(samrandp,file=paste(strp,k,strss[i],".csv",sep=""),row.names=FALSE)
    }
    
}

random_samples2 <- function(){
    filenames <-  c("data/cases0120.csv","data/PCGCinherited.csv","data/PCGCnontransmitted0120.csv") 
    trun <- c("frameshiftdeletion","frameshiftinsertion","stopgainSNV","stoplossSNV","frameshiftsubstitution")
    mis3 <- c("nonframeshiftdeletion","nonframeshiftinsertion","nonframeshiftsubstitution","nonsynonymousSNV")
    strss <- c("cade","cain","cano")
    ## random select samples
    
    allcase <- c()
    for(i in 1:length(filenames)){
        samples <- read.csv(filenames[i])
        n <- dim(samples)[1]
        if("proband.ID.GT.AD.DP.GQ.PL." %in% colnames(samples)){
            samplenames <- c()
            for(k in 1:n){
                b <- unlist(strsplit(gsub("[();]","\t",samples[k,"proband.ID.GT.AD.DP.GQ.PL."]),"\t"))
                samplenames <- union(b[seq(1,length(b),3)],samplenames)    
            }
        }else{
            samplenames <- unique(samples[,"Subject"])
        }
        allcase <- union(allcase,samplenames)
    }
    samplenames <- allcase
    n.sample <- length(samplenames)
    n.set <- 3
    for(j in 1:n.set){
        subs <- sample(n.sample,round(n.sample/2))
        subsamples <- samplenames[subs]
        for(i in 1:length(filenames)){
            samples <- read.csv(filenames[i])
            if("proband.ID.GT.AD.DP.GQ.PL." %in% colnames(samples)){
                samsubs <- sapply(1:n,function(k){ 
                    b <- unlist(strsplit(gsub("[();]","\t",samples[k,"proband.ID.GT.AD.DP.GQ.PL."]),"\t"))
                    if(any(b[seq(1,length(b),3)] %in% subsamples)){
                        paste(intersect(b[seq(1,length(b),3)],subsamples),collapse='();')
                    }else{FALSE;}
                })
                samrand <- samples[samsubs!=FALSE,]
                samrand[,"proband.ID.GT.AD.DP.GQ.PL."] <- samsubs[samsubs!=FALSE]
                write.csv(samrand,file=paste("random_samples/rand",j,strss[i],".csv",sep=""),row.names=FALSE)
            }else{
                samrand <- samples[samples[,"Subject"] %in% subsamples,]
                write.csv(samrand,file=paste("random_samples/rand",j,strss[i],".csv",sep=""),row.names=FALSE)
            }
        }
    }
    
}

sample.cross <- function(nsample,K){
    
    train_sample <- list()
    pred_sample <- list()
    
    nk <- floor(nsample/K)
    sam <- sample(nsample)
    
    for(i in 1:(K-1)){
        pred_sample[[i]] <- sam[((i-1)*nk+1):(i*nk)]
        train_sample[[i]] <- setdiff(sam,pred_sample[[i]])
    }
    pred_sample[[K]] <- sam[((K-1)*nk+1):nsample]
    train_sample[[K]] <- setdiff(sam,pred_sample[[K]])
    
    list(train=train_sample,pred=pred_sample)
}

PAH_samples <- function(){
    source("Network_analysis.R")
    fpath <- "/Users/qh2159/Dropbox (Personal)/genet_Pulmonary Hypertension/FilteringJan2015/trios"
    sufVec <- c(".denovo.annot.tsv$","maternal.annot.tsv$","paternal.annot.tsv$")
    
    PAHlist <- list()
    for(k in 1:3){
        pathfiles <- list.files(path = fpath, pattern=sufVec[k], full.names = TRUE, include.dirs=TRUE, recursive = TRUE)
        n <- length(pathfiles)
        samnames <- 1:n
        pahmut <- c()
        for(i in 1:n){
            tmp <- unlist(strsplit(pathfiles[i],"/"))
            samnames[i] <- unlist(strsplit(tmp[length(tmp)],"\\."))[1]
             
            tmp <- read.delim(pathfiles[i],sep="\t")
            colnames(tmp)[22:24] <- c("proband","Mother","Father") 
            tmp <- as.matrix(cbind(tmp,samnames[i]))
            pahmut <- as.matrix(rbind(pahmut,tmp))
        }
        PAHlist[[k]] <- pahmut
    }
    denovo <- PAHlist[[1]]
    inheri <- rbind(PAHlist[[2]],PAHlist[[3]])
    write.csv(denovo,file="PAH/PAHdenovo.csv",row.names=FALSE)
    write.csv(inheri,file="PAH/PAHinherited.csv",row.names=FALSE)

    denovo <- read.csv("PAH/PAHdenovo.csv")
    inheri <- read.csv("PAH/PAHinherited.csv")
    title1 <- c("Chromosome","Position","ID","REF","ALT","Gene","VariantFunction","VariantClass","AAchange","AlleleFrequency.ExAC","AlleleFrequencyKG","AlleleFrequency.ESP","SIFT","PolyPhen2","Mutation.Assessor","Mutation.Taster","GERP..","CADD.score","proband","Mother","Father","CLINVAR","SYNONYMS","RVIS_SCORE","TOLERANCE_ALL_DALY","TOLERANCE.specific")
    denovo <- denovo[,title1]
    inheri <- inheri[,title1]
    #====
    header <- colnames(denovo)
    denovo <- cbind(denovo,"denovo")
    colnames(denovo) <- c(header,"DENOVO")
    inheri <- cbind(inheri,"inherited")
    colnames(inheri) <- c(header,"DENOVO")
    samples <- rbind(denovo,inheri)
    colnames(samples)[colnames(samples)=="Gene"] <- "GeneName"
    colnames(samples)[colnames(samples)=="VariantFunction"] <- "VarFunc"
    colnames(samples)[colnames(samples)=="VariantClass"] <- "VarClass"
    
    samples <- as.matrix(samples)
    samples[samples[,"VarClass"] %in% "stopgain","VarClass"] <- "stopgainSNV"
    samples[samples[,"VarClass"] %in% "stoploss","VarClass"] <- "stoplossSNV"
    samples[samples=="."] <- NA
    samples[samples==""] <- NA

    trun <- c("frameshiftdeletion","frameshiftinsertion","stopgainSNV","stoplossSNV","frameshiftsubstitution")
    mis3 <- c("nonframeshiftdeletion","nonframeshiftinsertion","nonframeshiftsubstitution","nonsynonymousSNV")
    samples <- samples[samples[,"VarClass"] %in% union(trun,mis3),]
    samples <- samples[samples[,"VarFunc"] %in% c("exonic","exonic,splicing","splicing"),]
    
    strscore <- c("GERP..","CADD.score")
    str <- "PAH/PAH"
    flag=1
    prenode2(strscore,samples,trun,mis3,mutrate=0,str,flag)
}

prestep2 <- function(strscore,samples,trun,mis3,mutrate=0,str,flag=1){

    for(i in 1:length(strscore)){
        samples[,strscore[i]] <- change(samples[,strscore[i]],-2)
    }
    tmp <- sapply(1:dim(samples)[1],function(i) mean(as.numeric( samples[i, strscore[as.numeric(samples[i,strscore])>0] ] )))
    samples <- cbind(samples,tmp)
    colnames(samples)[dim(samples)[2]] <- "SCORE"
    score <- tmp
    
    tmp <- samples[,"GeneName"]
    tmp <- sapply(tmp, function(x){unlist(strsplit(x,","))[1]})
    tmp <- sapply(tmp, function(x){unlist(strsplit(x,"\\("))[1]})
    genes <- unique(tmp)
    gmap <- cbind(genes,genes)
    gmap[,2] <- mapping_to(gmap[,2])
    samples[,"GeneName"] <- gmap[match(tmp,gmap[,1]),2]
    genes <- unique(samples[,"GeneName"])
    
    # deal with NA
    label <- 1:dim(samples)[1]
    label[(samples[,"DENOVO"]=="denovo") & (samples[,"VarClass"] %in% trun)] <- 1
    label[(samples[,"DENOVO"]=="denovo") & (samples[,"VarClass"] %in% mis3)] <- 2
    label[(samples[,"DENOVO"]=="inherited") & (samples[,"VarClass"] %in% trun)] <- 3
    label[(samples[,"DENOVO"]=="inherited") & (samples[,"VarClass"] %in% mis3)] <- 4
    subs <- !is.na(score)
    a <- 1:4
    tmp <- (samples[,"DENOVO"]=="denovo") & (samples[,"VarClass"] %in% trun) & subs
    a[1] <- mean(as.numeric(score[tmp]))
    tmp <- (samples[,"DENOVO"]=="denovo") & (samples[,"VarClass"] %in% mis3) & subs
    a[2] <- mean(as.numeric(score[tmp]))
    tmp <- (samples[,"DENOVO"]=="inherited") & (samples[,"VarClass"] %in% trun) & subs
    a[3] <- mean(as.numeric(score[tmp]))
    tmp <- (samples[,"DENOVO"]=="inherited") & (samples[,"VarClass"] %in% mis3) & subs
    a[4] <- mean(as.numeric(score[tmp]))
    tmp <- dealNA(samples,label,trun,mis3,a,score)
    samples[,"SCORE"] <- tmp
    
    write.csv(samples,file=paste(str,"_node_SCORE.csv",sep=""),row.names=FALSE)
    
    prenode0(samples,trun,mis3,genes,mutrate,str,flag)

}

write.csv3 <- function(samples, ...){
    x <- samples[,"GeneName"]
    x <- gsub("^((JAN|FEB|MAR|APR|MAY|JUN|JUL|AUG|SEP|OCT|NOV|DEC)\\d+)$", "=\"\\1\"",x, ignore.case=T)
    samples[,"GeneName"] <- x
    write.csv(samples, ...)
}

read.csv3 <- function(...){
    samples <- read.csv(...)
    x <- samples[,"GeneName"]
    x <- gsub("^=\"((JAN|FEB|MAR|APR|MAY|JUN|JUL|AUG|SEP|OCT|NOV|DEC)\\d+)\"$", "\\1", toupper(x))
    samples[,"GeneName"] <- x
    samples
}

predeal <- function(LOF,Mis,filenames){
    
    gene_alias = matrix(c('AGO1','EIF2C1','BRINP3','FAM5C','CEP162','KIAA1009','CIPC','KIAA1737','CLUH','KIAA0664','KMT2D','MLL2','KMT2C','MLL3', 'TENM4','ODZ4','TRMT13','CCDC761','UVSSA','KIAA1530','VWA9','C15orf44','WDR45B','WDR45L','ZNF106','ZFP106','COLGALT1','GLT25D1','DPH6','ATPBD4','EMC1','KIAA0090','ENKD1','C16orf48', 'EQTN','C9orf11','GSE1','KIAA0182','HID1','C17orf28','KIAA1549L','C11orf41','MAP3K19','YSK4','MROH7','HEATR8','MYRF','C11orf9','NRROS','LRRC33','SBSPON','C8orf84','SPDYE2B','SPDYE2L','SZRD1','C1orf144'), 28,2,byrow=TRUE);
    
    n <- length(filenames)	
    title1 <- c("CHR","POS","REF","ALT","Subject","parent_min","esp6500si_all","Gene","Func.refgene","ExonicFunc.refgene","AAChange.refgene","X1000g2014oct_all","SIFT_score","Polyphen2_HDIV_score","Polyphen2_HDIV_pred","MutationTaster_score","GERP.._RS")
    title2 <- c("CHROM","POS","REF","ALT","proband.ID.GT.AD.DP.GQ.PL.","Inherited.from","esp.MAF","GeneName","VarFunc","VarClass","AAChange","X1KGfreq","SIFTscr","PP2scr","PP2prd","MutTscr","GERP")
    genes <- c()
    datalist <- list()
    for(nf in 1:n){
        dfile <- filenames[nf]
        samples <- as.matrix(read.csv(dfile))
        if(all(title2%in% colnames(samples) ) ){
            samples[is.na(samples[,"ESPfreq"]),"ESPfreq"] <- 0
            samples[is.na(samples[,"X1KGfreq"]),"X1KGfreq"] <- 0
            subs1 <- as.numeric(samples[,"X1KGfreq"]) < 0.0001 & as.numeric(samples[,"ESPfreq"]) < 0.0001
            samples[is.na(samples[,"esp.MAF"]),"esp.MAF"] <- "0,0,0"
            subs2 <- sapply(samples[,"esp.MAF"], function(val) as.numeric(unlist(strsplit(val,",")))[3] ) < 0.1
            subs <- subs1 & subs2
            samples <- samples[subs,]
            
            samples[samples[,"VarFunc"]=="splicing","VarClass"] <- "frameshiftdeletion"  ### splicing as LOF
            samples <- samples[samples[,"VarFunc"] %in% c("exonic","exonic,splicing","splicing"),]
            
            ###### Difference ===================
            ###samples <- samples[samples[,"VarClass"] %in% union(LOF,Mis),]	
            subs1 <- samples[,"VarClass"] %in% LOF
            subs2 <- samples[,"VarClass"] %in% Mis & samples[,"PP2prd"] %in% c("D","P")
            samples <- samples[subs1 | subs2,]
            
            samples <- samples[,title2]
        }else{		
            samples[samples[,"SIFT_score"]==".","SIFT_score"] = ""
            samples[samples[,"Polyphen2_HDIV_score"]==".","Polyphen2_HDIV_score"] = ""
            samples[samples[,"MutationTaster_score"]==".","MutationTaster_score"] = ""
            samples[samples[,"GERP.._RS"]==".","GERP.._RS"] = ""
            samples[samples[,"SIFT_score"]=="","SIFT_score"] = NA
            samples[samples[,"Polyphen2_HDIV_score"]=="","Polyphen2_HDIV_score"] = NA
            samples[samples[,"MutationTaster_score"]=="","MutationTaster_score"] = NA
            samples[samples[,"GERP.._RS"]=="","GERP.._RS"] = NA
            
            samples[is.na(samples[,"esp6500si_all"]),"esp6500si_all"] <- 0
            samples[is.na(samples[,"X1000g2014oct_all"]),"X1000g2014oct_all"] <- 0
            subs <- as.numeric(samples[,"esp6500si_all"]) < 0.0001 & as.numeric(samples[,"X1000g2014oct_all"]) < 0.0001
            samples <- samples[subs,]
            
            subs <- as.numeric(samples[,"parent_min"]) >= 20
            samples <- samples[subs,]
            
            varfunc <- c("exonic","exonic;exonic","exonic;splicing","intronic","ncRNA_exonic","splicing")
            varmap <- c("exonic","exonic","exonic,splicing","intronic","ncRNA_exonic","splicing")
            samples[,"Func.refgene"] <- varmap[match(samples[,"Func.refgene"],varfunc)]
            samples <- samples[samples[,"Func.refgene"] %in% c("exonic","exonic,splicing","splicing"),]
            mutype1 <- c("","frameshift deletion","frameshift insertion","nonframeshift deletion","nonframeshift insertion","nonsynonymous SNV","nonsynonymous SNV;nonsynonymous SNV", "stopgain","stoploss","synonymous SNV","unknown")
            mutmap <- c("frameshiftdeletion","frameshiftdeletion","frameshiftinsertion","nonframeshiftdeletion","nonframeshiftinsertion","nonsynonymousSNV","nonsynonymousSNV", "stopgainSNV","stoplossSNV","synonymousSNV","unknown") ##!!!!!!! splicing mutations are seen as LOF in samples.
            samples[,"ExonicFunc.refgene"] <- mutmap[match(samples[,"ExonicFunc.refgene"],mutype1)]
            samples <- samples[,title1]
            samples <- samples[samples[,"ExonicFunc.refgene"] %in% union(LOF,Mis),]
            colnames(samples) <- title2
            
            subs1 <- samples[,"VarClass"] %in% LOF
            subs2 <- samples[,"VarClass"] %in% Mis &  samples[,"PP2prd"] %in% c("D","P")
            samples <- samples[subs1 | subs2,]
        }
        tmp <- samples[,"GeneName"]
        tmp <- sapply(tmp, function(x){unlist(strsplit(x,","))[1]})
        tmp <- sapply(tmp, function(x){unlist(strsplit(x,"\\("))[1]})
        
        tmpg <- intersect(gene_alias[,1],tmp)
        tmp[match(tmpg,tmp)] <- gene_alias[match(tmpg,gene_alias[,1]),2]
        
        samples[,"GeneName"] <- tmp
        
        datalist[[nf]] <- samples
        #genes <- union(genes,as.vector(toupper(samples[,"GeneName"])))	
        genes <- union(genes,as.vector(samples[,"GeneName"]))	
    }
    datalist[[n+1]] <- genes
    
    datalist
}