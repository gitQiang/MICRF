enrichana <- function(){
    source("Network_analysis.R")
    source("enrichana.R")
#     genes <- read.table("ASD/iref_index_genes")[,2] ##!!!!
#     genes <- mapping_to(genes)
#     cutV <- seq(0.1,0.5,0.1)
#     TADAfiles <- c("ASD/Mutations13772.csv","ASD/Mutations13772.csv","ASD/Mutations13772.csv")
#     NDDlist <- list()    
#     NDDlist[[1]] <- c("ASD/Mutations13772.csv","result/ASD_hotnet13772LBP_61.txt","result/ASD13772_6LBP_61.txt")
#     NDDlist[[2]] <- c("ASD/Muta13908.csv","result/ASD_hotnet13908LBP_61.txt","result/ASD13908_6LBP_61.txt")
#     NDDlist[[3]] <- c("ASD/TADA_resultsMutap932.csv","result/ASD_hotnet932pLBP_61.txt","result/ASD932p_6LBP_61.txt")
#     dataname <- c("data set 1","data set 2","data set 3")
#     a1 <- enric3(NDDlist,TADAfiles,dataname,genes,flag=1)

    TADAfiles <- c("TADA_resultsPCGC_2_4.csv","TADA_resultsPCGC_2_4.csv","TADA_resultsPCGC_2_4.csv","TADA_resultsPCGC_2_4.csv","TADA_resultsPCGC_2_4.csv")
    NDDlist <- list()    
    NDDlist[[1]] <- c("random_samples/TADA_resultsrand1.csv","result/hotnetresult1rand16.txt","result/CRFresult_rand1_6LBP_61.txt","result/DAWNresultrand1.txt")
    NDDlist[[2]] <- c("random_samples/TADA_resultsrand2.csv","result/hotnetresult1rand26.txt","result/CRFresult_rand2_6LBP_61.txt","result/DAWNresultrand2.txt")
    NDDlist[[3]] <- c("random_samples/TADA_resultsrand3.csv","result/hotnetresult1rand36.txt","result/CRFresult_rand3_6LBP_61.txt","result/DAWNresultrand3.txt")
    NDDlist[[4]] <- c("TADA_resultsPCGC_2_4.csv","result/hotnetresult1PCGC6.txt","result/CRFresult_PCGC_6LBP_61.txt","result/DAWNresultallPCGC.txt")
    dataname <- c("rand set 1","rand set 2","rand set 3","ALL samples")
    #a3 <- enric3(NDDlist,TADAfiles,dataname,genes,flag=1)
    enric5(NDDlist,TADAfiles,dataname,genes,flag=1,N=0.3)

    source("Network_analysis.R")
    source("enrichana.R")    
    TADAfiles <- c("ASD/TADAresult/TADA_resultsall.csv","ASD/TADAresult/TADA_resultsall.csv","ASD/TADAresult/TADA_resultsall.csv","ASD/TADAresult/TADA_resultsall.csv")
    NDDlist <- list()    
    NDDlist[[1]] <- c("ASD/TADAresult/TADA_results13772.csv","result/hotnetresult1ASD137726.txt","result/CRFresult_ASD13772_6LBP_61.txt","result/DAWNresult13772.txt")
    NDDlist[[2]] <- c("ASD/TADAresult/TADA_results13908.csv","result/hotnetresult1ASD139086.txt","result/CRFresult_ASD13908_6LBP_61.txt","result/DAWNresult13908.txt")
    NDDlist[[3]] <- c("ASD/TADAresult/TADA_results932.csv","result/hotnetresult1ASD9326.txt","result/CRFresult_ASD932_6LBP_61.txt","result/DAWNresult932.txt")
    NDDlist[[4]] <- c("ASD/TADAresult/TADA_results932.csv","result/hotnetresult1denovo9326.txt","result/CRFresult_ASDdenovo932_6LBP_61.txt","result/DAWNresult932.txt")
    dataname <- c("data set 1 with 2270 trios","data set 2 with 1735 trios","data set 3 with 932 trios","data set 3 with 932 trios")
    genes=""
    enric5(NDDlist,TADAfiles,dataname,genes,flag=1,N=0.1)
    
    source("enrichana.R")
    NDDlist <- list()
    NDDlist[[1]] <- c("ASD/TADAresult/hotnet_inputall.txt","result/hotnet_result_ASDall.txt","result/CRFresult_ASDall_6LBP_61.txt")
    NDDlist[[2]] <- c("ASD/TADAresult/hotnet_input13772.txt","result/hotnet_result_ASD13772.txt","result/CRFresult_ASD13772_6LBP_61.txt")
    NDDlist[[3]] <- c("ASD/TADAresult/hotnet_input13908.txt","result/hotnet_result_ASD13908.txt","result/CRFresult_ASD13908_6LBP_61.txt")
    NDDlist[[4]] <- c("ASD/TADAresult/hotnet_input932.txt","result/hotnet_result_ASD932.txt","result/CRFresult_ASD932_6LBP_61.txt")
    for(i in 1:4){
        plotecdf(NDDlist[[i]]);
    }
}

enric3 <- function(NDDlist,TADAfiles,dataname,genes,flag=1,N=c(12129,12129,12129)){
    cutV <- seq(0.1,0.5,0.1)
    sMlist <- list()
    cMlist <- list()
    Llist <- list()
    tMlist <- list()
    nmeth <- length(NDDlist[[1]])
    ncut <- length(cutV)
    nst <- length(TADAfiles)
    
    TADAall <- read.csv(TADAfiles[1])
    TADAall[,1] <- mapping_to(TADAall[,1])
    TADAall <- TADAall[TADAall[,1] %in% genes,]
    for(cut in cutV){
        j <- which(cutV %in% cut)
        sM <- matrix(0,nmeth,nst)
        cM <- matrix(0,nmeth,nst)
        Lv <- matrix(0,nmeth,nst) 
        tM <- matrix(0,nmeth,nst)
        
        Tset <- TADAall[TADAall[,"qvalue"]< 0.3,1]
        #Tset <- mapping_to(Tset)
        TADAset <- TADAall[TADAall[,"qvalue"]< cut,1]
        #TADAset <- mapping_to(TADAset)
        n.all <- length(TADAset)
        for(d in 1:nst){     
            for(i in 1:nmeth){
                NDDfiles <- NDDlist[[d]]    
                if(i==1){
                    result <- read.csv(NDDfiles[i]) 
                    if(flag==1){
                        NDDset <- result[result[,"qvalue"]< cut,1]
                    }else{
                        NDDset <- result[1:n.all,1]
                    }
                    if(length(NDDset)>0 ){
                        NDDset <- mapping_to(NDDset)
                    }
                }else{
                    result <- read.table(NDDfiles[i])
                    if(flag==1){
                        NDDset <- result[result[,5]< cut,1]
                    }else{
                        NDDset <- result[1:n.all,1]
                    }
                }
                
                m <- n.all
                n <- N[i] - m
                k <- length(NDDset)
                q <- length(intersect(NDDset,TADAset))
                sM[i,d] <- phyper(q-1, m, n, k, lower.tail = FALSE, log.p = FALSE)
                cM[i,d] <- q
                Lv[i,d] <- k
                tM[i,d] <- length(intersect(NDDset,Tset))
            }
        }
        sMlist[[j]] <- sM
        cMlist[[j]] <- cM
        Llist[[j]] <- Lv  
        tMlist[[j]] <- tM
        barplot(cM, main=paste("FDR < ",cut,sep=""), xlab="Groups of data sets", ylab="The number of overlap genes", col=1:nmeth, legend = c("TADA","RWR","Our"), beside=TRUE,args.legend=list(x = "bottomright"),xaxt="n")
        axis(side = 1, at = seq(2.5,2.5+4*(nst-1),4), labels = dataname)
        
        barplot(tM, main=paste(" TADA (FDR < 0.3) overlap with FDR < ",cut,sep=""), xlab="Groups of data sets", ylab="The number of overlap genes", col=1:nmeth, legend = c("TADA","RWR","Our"), beside=TRUE,args.legend=list(x = "bottomright"),xaxt="n")
        axis(side = 1, at = seq(2.5,2.5+4*(nst-1),4), labels = dataname)
        
    }
    
    a <- matrix(0,nst*nmeth,ncut*3)
    for(i in 1:ncut){
        a[,seq(1,ncut*3,3)[i]] <- as.vector(cMlist[[i]])
        a[,seq(2,ncut*3,3)[i]] <- as.vector(Llist[[i]])
        a[,seq(3,ncut*3,3)[i]] <- as.vector(tMlist[[i]])
    }
    
    a
    
}

enric4 <- function(NDDlist,TADAfiles,dataname,genes,flag=1,N=c(12129,12129,12129)){
    
    nmeth <- length(NDDlist[[1]])
    nst <- length(TADAfiles)
    TADAall <- read.csv(TADAfiles[1])
    #TADAall[,1] <- mapping_to(TADAall[,1])
    #TADAall <- TADAall[TADAall[,1] %in% genes,];
    
    resultT <- list()
    cutV <- seq(0,0.5,0.01) 
    ncut <- length(cutV)
    Tset <- TADAall[TADAall[,"qvalue"]<= 0.3,1]
    n.all <- length(Tset)
    Fset <- setdiff(genes,Tset)
    
        for(d in 1:nst){    
            NDDfiles <- NDDlist[[d]] 
            sM <- matrix(0,4*nmeth,ncut)
            for(i in 1:nmeth){
                if(i==1){ result <- read.csv(NDDfiles[i]);
                }else{ result <- read.table(NDDfiles[i]);}
                for(cut in cutV){
                    j <- which(cutV %in% cut)
                    if(i==1){
                        NDDset <- result[result[,"qvalue"]<= cut,1]
                        #if(length(NDDset)>0 ){ NDDset <- mapping_to(NDDset);}
                    }else{ NDDset <- result[result[,5]<= cut,1];}
                    NDDFset <- setdiff(genes,NDDset)
                    sM[((i-1)*4+1),j] <- length(intersect(NDDset,Tset))
                    sM[((i-1)*4+2),j] <- length(intersect(NDDset,Fset))
                    sM[((i-1)*4+3),j] <- length(intersect(NDDFset,Tset))
                    sM[((i-1)*4+4),j] <- length(intersect(NDDFset,Fset))
                }
                x <- sM[((i-1)*4+2),]/(sM[((i-1)*4+2),]+sM[((i-1)*4+4),])
                y <- sM[((i-1)*4+1),]/(sM[((i-1)*4+1),]+sM[((i-1)*4+3),])
                if(i==1){
                    plot(x,y,main=paste("ROC plot on ",dataname[d],sep=""),xlab="1-specificity",ylab="sensitivity",col=i,type="l",xlim=c(0,1),ylim=c(0,1));
                }else{lines(x,y,col=i,type="l");}
            }
            lines(x=c(0,1),y=c(0,1),col=i+1,lty=2)
            legend("bottomright",col=1:nmeth,legend=c("TADA","RWR","Our"),lty=c(1,1,1),lwd=c(1,1,1))
            resultT[[d]] <- sM
        }
    
    resultT
    
}

enric5 <- function(NDDlist,TADAfiles,dataname,genes,flag=1,N=0.3,qval="qvalue"){
    library(ROCR)
    nmeth <- length(NDDlist[[1]])
    nst <- length(NDDlist)
    TADAall <- read.csv(TADAfiles[1])
    Tset <- TADAall[TADAall[,qval]<= N,1]

    for(d in 1:nst){    
        NDDfiles <- NDDlist[[d]]
        predictions <- list()
        labels <- list()
        for(i in 1:nmeth){
            if(grepl(".csv",NDDfiles[i]) & grepl("DAWN",NDDfiles[i])){ result <- read.csv(NDDfiles[i]); predictions[[i]] <- result[,"FDR"]; predictions[[i]][is.na(predictions[[i]])] <- 1;
            }else if(grepl(".csv",NDDfiles[i]) & grepl("TADA",NDDfiles[i])){ result <- read.csv(NDDfiles[i]); predictions[[i]] <- result[,"qvalue.dn"];
            }else if(grepl(".txt",NDDfiles[i])){result <- read.table(NDDfiles[i]); predictions[[i]] <- result[,5]}
            
            if(flag!=1){result <- result[1:flag,]; predictions[[i]] <- predictions[[i]][1:flag];}
            
            predictions[[i]] <- 1 - predictions[[i]] 
            labels[[i]] <- rep(-1,length(predictions[[i]]))
            labels[[i]][match(Tset,result[,1])] <- 1
        }
        pred <- prediction(predictions, labels)
        auc <- unlist(performance(pred, "auc")@y.values)
        perf <- performance(pred, measure = "tpr", x.measure = "fpr") 
        plot(perf,main=paste("ROC plot on ",dataname[d],sep=""),xlab="1-specificity",ylab="sensitivity")
        for(i in 1:nmeth){
                lines(perf@x.values[[i]], perf@y.values[[i]], col = i)
        }
        lines(x=c(0,1),y=c(0,1),col=i+1,lty=2)
        if(nmeth==3){legend=c("TADA","RWR","MIS");}else if(nmeth==4){legend=c("TADA","DAWN","RWR","MIS");}
        legend <- paste(legend,format(auc,digits=3),sep="--AUC:")
        legend("bottomright",col=(1:nmeth),legend=legend,lty=rep(1,nmeth),lwd=rep(1,nmeth))
    }
        
}

plotecdf <- function(NDDfiles){

    for(i in 1:length(NDDfiles)){
        p1 <- read.table(NDDfiles[i])[,2];
        if(i==1){
            plot(ecdf(p1-min(p1)),col=i,main="Empirical distribution for gene scores.")
        }else{lines(ecdf(p1-min(p1)),col=i);}
    }
    #abline(x=c(0,1),y=c(0,1),col=i+1,lty=2)
    legend("bottomright",col=1:length(NDDfiles),lty=c(1,1,1),lwd=c(1,1,1),legend=c("TADA","RWR","Our"))

}

net_EP300 <- function(){
    if(netflag==3){ filename <- "STRINGnetmap.txt";}
    net.text <- as.matrix(read.table(filename,sep="\t",header=FALSE))
    subs <- net.text[,1]=="ENSG00000100393" | net.text[,2]=="ENSG00000100393"
    genes <- union(net.text[subs,1],net.text[subs,2])
    #write.table(genes,file="nodeneighs.txt",quote=FALSE,col.names=FALSE,row.names=FALSE)
    #subs <- net.text[,1] %in% genes & net.text[,2] %in% genes
    subnet <- net.text[subs,]
    write.table(subnet,file="ep300net.txt",quote=FALSE,col.names=FALSE,row.names=FALSE)
    
    nodeinfo <- as.matrix(read.table("TADAinfo0206.txt"))
    a <- paste(nodeinfo[,1]," = ", nodeinfo[,2],sep="")
    write.table(a,file="nodeTADA.txt",quote=FALSE,col.names=FALSE,row.names=FALSE)

}

net_UBC <- function(){
    if(netflag==3){ filename <- "STRINGnetmap.txt";}
    net.text <- as.matrix(read.table(filename,sep="\t",header=FALSE))
    subs <- net.text[,1]=="ENSG00000150991" | net.text[,2]=="ENSG00000150991"
    genes <- union(net.text[subs,1],net.text[subs,2])
    #write.table(genes,file="nodeneighs.txt",quote=FALSE,col.names=FALSE,row.names=FALSE)
    #subs <- net.text[,1] %in% genes & net.text[,2] %in% genes
    subnet <- net.text[subs,]
    write.table(subnet,file="UBCnet.txt",quote=FALSE,col.names=FALSE,row.names=FALSE)
    
    nodeinfo <- as.matrix(read.table("TADAinfo0206.txt"))
    a <- paste(nodeinfo[,1]," = ", nodeinfo[,2],sep="")
    write.table(a,file="nodeTADA.txt",quote=FALSE,col.names=FALSE,row.names=FALSE)
    
}

wirtetmp <- function(){

###  section 3.1
tmpfunc <- function(){
mapT <- as.matrix(read.delim("Fmap0121.txt",header=FALSE,sep="\t"))
heart1 <- as.matrix(read.table("GRN_net_heart.txt"))
heart2 <- as.matrix(read.table("function_net_heart.txt"))
genes1 <- union(heart1[,1],heart1[,2])
genes2 <- union(heart2[,1],heart2[,2])
genes <- union(genes1,genes2)
genesm <- mapping_to(genes)
genesmap <- cbind(genes,genesm)

net <- rbind(heart1[,1:2],heart2[,1:2])
net[,1] <- genesmap[match(net[,1],genes),2]
net[,2] <- genesmap[match(net[,2],genes),2]

genesm1 <- mapT[match(genesm,mapT[,2]),1]
genemap1 <- cbind(genesm1,genesm)
net[,1] <- genemap1[match(net[,1],genesm),1]
net[,2] <- genemap1[match(net[,2],genesm),1]
net <- net[!is.na(net[,1]) & !is.na(net[,2]),]
write.table(net,file="heart_net.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)


mapT <- as.matrix(read.delim("Fmap0121.txt",header=FALSE,sep="\t"))
net <- read.csv("Human_Big_GRN_032014.csv")
genes <- union(net[,1],net[,2])
genesm <- mapping_to(genes)
genesmap <- cbind(genes,genesm)
net[,1] <- genesmap[match(net[,1],genes),2]
net[,2] <- genesmap[match(net[,2],genes),2]
genes <- union(net[,1],net[,2])
genesm1 <- mapT[match(genes,mapT[,2]),1]
genemap1 <- cbind(genesm1,genes)
net[,1] <- genemap1[match(net[,1],genes),1]
net[,2] <- genemap1[match(net[,2],genes),1]
net <- net[!is.na(net[,1]) & !is.na(net[,2]),]
write.table(net[,c(1,2,4)],file="GRN_net.txt",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")

}

# add FDR to our score 
Our <- as.matrix(read.table("NewPheno_score.txt"))
BF <- as.numeric(Our[,2])/(1-as.numeric(Our[,2]))
FDR <- Bayesian.FDR(BF)$FDR
Our <- cbind(Our,FDR)
Our <- Our[order(FDR,decreasing=FALSE),]
write.table(Our,file="NewPheno_score_FDR.txt",col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)

mapT <- as.matrix(read.delim("Fmap0121.txt",header=FALSE,sep="\t"))
Output_file1("result/RWR_1LBP_11.txt",mapT,"result/RWR1_20.txt",1,flag=20)
Output_file1("result/RWR_1LBP_11.txt",mapT,"result/RWR1_500.txt",1,flag=500)

Output_file1("result/RWR_3LBP_31.txt",mapT,"result/RWR3_20.txt",1,flag=20)
Output_file1("result/RWR_3LBP_31.txt",mapT,"result/RWR3_500.txt",1,flag=500)

Output_file1("result/RWR_4LBP_41.txt",mapT,"result/RWR4_20.txt",1,flag=20)
Output_file1("result/RWR_4LBP_41.txt",mapT,"result/RWR4_500.txt",1,flag=500)

Output_file1("result/RWR_8_resultall.txt",mapT,"result/RWR8_20.txt",1,flag=20)
Output_file1("result/RWR_8_resultall.txt",mapT,"result/RWR8_500.txt",1,flag=500)



tadafin <- read.csv("TADA_resultsPCGC.csv")
rASDgene <- as.vector(tadafin[tadafin[,"qvalue"] < 0.5,1])
Our <- as.matrix(read.table("NewPheno_score_FDR.txt"))
Topgenes <- Our[as.numeric(Our[,3])<0.5,1]
length(intersect(Topgenes,rASDgene))
genes <- intersect(Topgenes,rASDgene)
x <- match(genes,Topgenes)
y <- match(genes,rASDgene)
wilcox.test(x,y)
t.test(x,y)

## FDR < 0.5, compare with DAWN
dawnfin <- read.csv("PCGCresult.csv")
dawnfin[is.na(dawnfin[,"DAWN.FDR"]),"DAWN.FDR"] <- 1
rASDgene <- as.vector(dawnfin[as.numeric(dawnfin[,"DAWN.FDR"])<0.5,1])
#dawnfin[is.na(dawnfin[,"DAWN.Stage2_posterior"]),"DAWN.Stage2_posterior"] <- 1
#rASDgene <- as.vector(dawnfin[as.numeric(dawnfin[,"DAWN.Stage2_posterior"])<0.5,1])
fileflag <- 1
testname <- "DAWNFinal0127.txt"
mapT <- as.matrix(read.delim("Fmap0121.txt",header=FALSE,sep="\t"))
result <- dawnfin
flag <- 0.5

# RWR 


tadafin <- read.csv("TADA_resultsPCGC.csv")
tadafin[is.na(tadafin[,"qvalue"]),"qvalue"] <- 1
rASDgene <- as.vector(tadafin[as.numeric(tadafin[,"qvalue"])<0.5,1])
fileflag <- 1
testname <- "tadaFinal0127.txt"
mapT <- as.matrix(read.delim("Fmap0121.txt",header=FALSE,sep="\t"))
result <- tadafin

resultS <- as.matrix(read.table("NewPheno_score_FDR.txt"))
resultS[is.na(resultS[,3]),3] <- 1
rASDgene <- as.vector(resultS[as.numeric(resultS[,3])<0.5,1])
fileflag <- 1
testname <- "ScoreFinal0127.txt"
mapT <- as.matrix(read.delim("Fmap0121.txt",header=FALSE,sep="\t"))
result <- resultS

resultO <- as.matrix(read.table("result/Newscore_6_resultall.txt"))
rASDgene <- as.vector(resultO[as.numeric(resultO[,12]) <0.5,1])
fileflag <- 1
testname <- "result/OurscoreFinal0127.txt"
mapT <- as.matrix(read.delim("Fmap0121.txt",header=FALSE,sep="\t"))
result <- resultO

###=============== gene specific score based ========================
resultO <- as.matrix(read.table("result/Newscore_6_resultall0.txt"))
One <- as.vector(resultO[as.numeric(resultO[,13]) < 0.1,1])
fnode <- c("random_samples/NDDPheno_score_map.txt","random_samples/nonNDDPheno_score_map.txt","random_samples/otherPheno_score_map.txt","random_samples/CHDnaturePheno_score_map.txt","random_samples/CHDotherPheno_score_map.txt","random_samples/P2Pheno_score_map.txt","random_samples/T2Pheno_score_map.txt","random_samples/P3Pheno_score_map.txt","random_samples/T3Pheno_score_map.txt")

# for(i in 1:9){
# Our <- as.matrix(read.table(fnode[i]))
# BF <- as.numeric(Our[,2])/(1-as.numeric(Our[,2]))
# FDR <- Bayesian.FDR(BF)$FDR
# Our <- cbind(Our,FDR)
# Our <- Our[order(FDR,decreasing=FALSE),]
# write.table(Our,file=gsub("_map","_FDR",fnode[i]),col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)
# }
library(exactRankTests) 
flag=6
strn0 <- c("result/NDD","result/nonNDD","result/other","result/CHDnature","result/CHDother","result/P2","result/T2","result/P3","result/T3")
pvs <- matrix(-1,9,3)
for(i in 1:9){
    tmp1 <- as.matrix(read.table(gsub("_map","_FDR",fnode[i])))
    #genes <- tmp1[as.numeric(tmp1[,dim(tmp1)[2]]) < 0.1,1]
    filename <- paste(strn0[i],flag,"resultall0.txt",sep="_")
    tmp2 <- as.matrix(read.table(filename))
    genes <- tmp2[as.numeric(tmp2[,dim(tmp2)[2]]) < 0.1,1]
    x <- match(genes, tmp1[,1])
    y <- match(genes, tmp2[,1])
    z <- match(genes,One)
    pvs[i,1] <- wilcox.exact(y,x,alternative="less",paired=TRUE,exact=TRUE)$p.value
    pvs[i,2] <- wilcox.exact(z,x,alternative="less",paired=TRUE,exact=TRUE)$p.value
    pvs[i,3] <- wilcox.exact(z,y,alternative="less",paired=TRUE,exact=TRUE)$p.value
}

fnode <- c("random_samples/NDDPheno_score_map.txt","random_samples/nonNDDPheno_score_map.txt","random_samples/otherPheno_score_map.txt","random_samples/CHDnaturePheno_score_map.txt","random_samples/CHDotherPheno_score_map.txt","random_samples/P2Pheno_score_map.txt","random_samples/T2Pheno_score_map.txt","random_samples/P3Pheno_score_map.txt","random_samples/T3Pheno_score_map.txt")

for(i in 1:9){
Our <- as.matrix(read.table(gsub("_map","",fnode[i])))
BF <- as.numeric(Our[,2])/(1-as.numeric(Our[,2]))
FDR <- Bayesian.FDR(BF)$FDR
Our <- cbind(Our,FDR)
Our <- Our[order(FDR,decreasing=FALSE),]
write.table(Our,file=gsub("_map","_FDR0",fnode[i]),col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)
}
mapT <- as.matrix(read.delim("Fmap0121.txt",header=FALSE,sep="\t"))
fileflag = 1
flag = 0.5
str <- c("NDD","nonNDD","other","CHD","chdother","P2","T2","P3","T3")
for(i in 1:9){
    testname <- paste("GeneScore",str[i],".txt",sep="")
    Output_file1(gsub("_map","_FDR0",fnode[i]),mapT,testname,fileflag,flag)
}


###=============== TADA score based ========================
resultO <- as.matrix(read.table("result/Newtada_6_resultall0.txt"))
One <- as.vector(resultO[as.numeric(resultO[,13]) < 0.1,1])
str <- c("CHD","chdother","NDD","nonNDD","other","P2","T2","P3","T3")
fnode <- paste("TADAinfo",str,".txt",sep="")
library(exactRankTests) 
flag=6
strn0 <- c("result/NDD1","result/nonNDD1","result/other1","result/CHDnature1","result/CHDother1","result/P21","result/T21","result/P31","result/T31")
pvs <- matrix(-1,9,3)
for(i in 1:9){
    tmp1 <- as.matrix(read.table(fnode[i]))
    genes <- intersect(tmp1[,1],One)
    filename <- paste(strn0[i],flag,"resultall0.txt",sep="_")
    tmp2 <- as.matrix(read.table(filename))
    genes <- intersect(genes,tmp2[,1])
    x <- match(genes, tmp1[,1])
    y <- match(genes, tmp2[,1])
    z <- match(genes,One)
    pvs[i,1] <- wilcox.exact(y,x,alternative="less",paired=TRUE,exact=TRUE)$p.value
    pvs[i,2] <- wilcox.exact(z,x,alternative="less",paired=TRUE,exact=TRUE)$p.value
    pvs[i,3] <- wilcox.exact(z,y,alternative="less",paired=TRUE,exact=TRUE)$p.value
}

pvs[,1]

L <- c(437,299,357,310,181,64,89,41,130)
strn0 <- c("result/NDD1","result/nonNDD1","result/other1","result/CHDnature1","result/CHDother1","result/P21","result/T21","result/P31","result/T31") 
for(k in 1:length(strn0)){
    filenames <- c()
    for(i in 1:4){
        strn <- paste(strn0[k],i,sep="_")
        filenames <- union(filenames,paste(strn,"LBP_",i,".txt",sep=""))
    }
    flag=6
    rank_combine(filenames,strn0[k],flag,L[k])
}  


fnode <- c("TADA_resultsCHD.csv","TADA_resultschdother.csv","TADA_resultsNDD.csv","TADA_resultsnonNDD.csv","TADA_resultsother.csv","TADA_resultsP2.csv","TADA_resultsT2.csv","TADA_resultsP3.csv","TADA_resultsT3.csv")
fnode <- paste("result/",fnode,sep="")
mapT <- as.matrix(read.delim("Fmap0121.txt",header=FALSE,sep="\t"))
fileflag = 1
flag = 0.5
str <- c("NDD","nonNDD","other","CHD","chdother","P2","T2","P3","T3")
for(i in 1:9){
    testname <- paste("TADAenr",str[i],".txt",sep="")
    Output_file1(fnode[i],mapT,testname,fileflag,flag)
}



### TADA BF to probability
mapT <- as.matrix(read.delim("Fmap0121.txt",header=FALSE,sep="\t"))
TADAFile <- "TADA_resultsPCGC_2_4.csv"
geneinfo <- as.matrix(read.csv(TADAFile))
genelist <- geneinfo[,1]
genelist <- mapping_to(genelist)
genelist <- mapT[match(genelist,mapT[,2]),1]
pi0 <- 0.85 #!!!!!!
pi <- 1-pi0
nodesim <- as.matrix((as.numeric(geneinfo[,"BF"])*pi)/(as.numeric(geneinfo[,"BF"])*pi + 1-pi))
names(nodesim) <- genelist
subs <- !is.na(genelist)
nodesim <- nodesim[subs]
write.table(cbind(names(nodesim),nodesim),file="TADAinfo0206.txt",quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")

nodes  <- as.matrix(read.table("TADAinfo0206.txt"))
nodesV <- as.numeric(nodes[,2])
k=2
library(nor1mix)
re <- norMixEM(nodesV,k)
paM <- matrix(re[1:(3*k)],k,3)
vaM <- matrix(c(dnorm(nodesV,mean=paM[1,1],sd=paM[1,2]) * paM[1,3],dnorm(nodesV,mean=paM[2,1],sd=paM[2,2]) * paM[2,3]),length(nodesV),k)
P <- pnorm(nodesV,mean=re[k],sd=re[2*k])*exp(log(vaM[,k])-log(rowSums(vaM)))
a <- sort(P,decreasing=TRUE)
a[1]-a[50]
sum(a>0.5)
nodes2 <- nodes
nodes2[,2] <- P
write.table(nodes2,file="TADAinfo0206_c.txt",quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")

# library(MASS)
# pa <- fitdistr(nodesV[nodesV>0.05],"lognormal")
# nodesV <- plnorm(nodesV, meanlog=pa$estimate[1], sdlog = pa$estimate[2], lower.tail = TRUE, log.p = FALSE)
# a <- sort(nodesV,decreasing=TRUE)
# a[1]-a[50]
# pa <- fitdistr(nodesV[nodesV>0.05],"weibull")
# nodesV <- pweibull(nodesV, shape=pa$estimate[1], scale = pa$estimate[2], lower.tail = TRUE, log.p = FALSE)
# a <- sort(nodesV,decreasing=TRUE)
# a[1]-a[50]
# sum(a>0.5)
# # zpV <- (nodesV-mean(nodesV))/sd(nodesV)
# # pV <- pnorm(zpV)
# # score1 <- 1-pV
# # nodes2 <- nodes
# nodes2[,2] <- nodesV
# write.table(nodes2,file="TADAinfo0204_c.txt",quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")



filenames <- c("TADA_resultsCHD.csv","TADA_resultschdother.csv","TADA_resultsNDD.csv","TADA_resultsnonNDD.csv","TADA_resultsother.csv","TADA_resultsP2.csv","TADA_resultsT2.csv","TADA_resultsP3.csv","TADA_resultsT3.csv")
filenames <- paste("result/",filenames,sep="")
mapT <- as.matrix(read.delim("Fmap0121.txt",header=FALSE,sep="\t"))
str <- c("CHD","chdother","NDD","nonNDD","other","P2","T2","P3","T3")
for(i in 1:9){
    TADAFile <- filenames[i]
    geneinfo <- as.matrix(read.csv(TADAFile))
    genelist <- geneinfo[,1]
    genelist <- mapping_to(genelist)
    genelist <- mapT[match(genelist,mapT[,2]),1]
    pi0 <- 0.94
    pi <- 1-pi0
    nodesim <- as.matrix((as.numeric(geneinfo[,"BF"])*pi)/(as.numeric(geneinfo[,"BF"])*pi + 1-pi))
    names(nodesim) <- genelist
    subs <- !is.na(genelist)
    nodesim <- nodesim[subs]
    write.table(cbind(names(nodesim),nodesim),file=paste("TADAinfo",str[i],".txt",sep=""),quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")
}

for(i in 1:9){
    TADAFile <- filenames[i]
    geneinfo <- as.matrix(read.csv(TADAFile))
    genelist <- geneinfo[,1]
    genelist <- mapping_to(genelist)
    genelist <- mapT[match(genelist,mapT[,2]),1]
    pi0 <- 0.94
    pi <- 1-pi0
    nodesim <- as.matrix((as.numeric(geneinfo[,"BF"])*pi)/(as.numeric(geneinfo[,"BF"])*pi + 1-pi))
    names(nodesim) <- genelist
    FDR <- as.numeric(geneinfo[,"qvalue"])
    subs <- !is.na(genelist)
    nodesim <- nodesim[subs]
    FDR <- FDR[subs]
    tadar <- data.frame(names(nodesim),nodesim,FDR)
    write.table(tadar,file=paste("TADAinfoFDR",str[i],".txt",sep=""),quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")
}




#======== 01 26 before ========
## Figure1
a1 <- c(2.7E-7,1.6E-7,7.0E-7,0.001,2.09E-9)
a2 <- c(0.01,0.011,0.0093,0.0066,0.00088)

x<- 1:5
plot(x,-log10(a2),col=1,xaxt="n",xlab="Network",ylab="-log p values",main="HMG and TGF beta pathways enrichment",type="b",ylim=c(1,10))
lines(x,-log10(a1),col=2,type="b")
axis(1,at=x,labels=c("Co-expression","GRN","PPI","Function","Combined"))
legend("topright",col=c(1,2),lty=1,legend=c("HMG","TGF-beta"))

## TGF genes in NDD and non NDD cases:: Table 6
source("prenode.R")
fk <- function(tablelist,k,genes){
    gene <- genes[k]
    nn <- 1:4
    sn <- 1:4    
    for(i in 1:4){
        tmp <- tablelist[[i]][ tablelist[[i]][,"GeneName"]==gene, "SCORE"]
        nn[i] <- length(tmp)
        sn[i] <- sum(as.numeric(tmp))
    }
    print(nn)
    print(sn)
}
trun <- c("frameshiftdeletion","frameshiftinsertion","stopgainSNV","stoplossSNV","frameshiftsubstitution")
mis3 <- c("nonframeshiftdeletion","nonframeshiftinsertion","nonframeshiftsubstitution","nonsynonymousSNV")
mutrate <- 0


samples <- read.csv("PCGC1_node_SCORE.csv")
samples <- read.csv("random_samples/NDD_node_SCORE.csv")
samples <- read.csv("random_samples/nonNDD_node_SCORE.csv")

genes <- as.vector(unique(samples[,"GeneName"]))
tablelist <- list()
tablelist[[1]] <- samples[(samples[,"DENOVO"]=="denovo") & (samples[,"VarClass"] %in% trun) ,]
tablelist[[2]] <- samples[(samples[,"DENOVO"]=="denovo") & (samples[,"VarClass"] %in% mis3) ,]
tablelist[[3]] <- samples[(samples[,"DENOVO"]=="inherited") & (samples[,"VarClass"] %in% trun) ,]
tablelist[[4]] <- samples[(samples[,"DENOVO"]=="inherited") & (samples[,"VarClass"] %in% mis3) ,]

k = which(genes=="SMAD4")
fk(tablelist,k,genes)
k = which(genes=="SMAD2")
fk(tablelist,k,genes)
k = which(genes=="SMURF2")
fk(tablelist,k,genes)
k = which(genes=="PITX2")
fk(tablelist,k,genes)


### other samples 
NDD <- as.matrix(read.table("result/NDD_6_resultall.txt"))
nonNDD <- as.matrix(read.table("result/nonNDD_6_resultall.txt"))
other <- as.matrix(read.table("result/other_6_resultall.txt"))
our <- as.matrix(read.table("result/PCGClog_6_resultall.txt"))
ngene <- 300
overlapg <- function(S1,S2,S3,our,ngene){
    o <- 1:4
    o[1] <- length(intersect(S1[1:ngene,1],our[1:ngene,1]))
    o[2] <- length(intersect(S2[1:ngene,1],our[1:ngene,1]))
    o[3] <- length(intersect(S3[1:ngene,1],our[1:ngene,1]))
    o[4] <- length(intersect(S1[1:ngene,1],intersect(S2[1:ngene,1],S3[1:ngene,1])))
    o
}
o=overlapg(NDD,nonNDD,other,our,ngene)

a<- matrix(c(o[1:3],387,393,469),2,3,byrow=TRUE)
chisq.test(a)



### PAH enrichment analysis
source("Network_analysis.R")
knowngenes <- as.matrix(read.table("PAH/knowngenes.txt"))
knowngenes <- mapping_to(knowngenes)
result <- read.table("PAH/PAH_6_resultall.txt")
rans <- match(knowngenes,result[,1])

### Nature paper and others
chdgenes <- read.table("CHDnonslientgenes.txt")
chdgenes <- as.matrix(chdgenes)
chdgenes <- unique(chdgenes)
source("Network_analysis.R")
chdgenes <- mapping_to(chdgenes)
length(chdgenes)

chdresult <- read.table("result/CHDnature_6_resultall.txt")
chdother <- read.table("result/CHDother_6_resultall.txt")
our <- read.table("result/PCGClog_6_resultall.txt")

ngene=231
a1 <- sum(chdgenes %in% chdresult[1:ngene,1])
a2 <- sum(chdgenes %in% chdother[1:ngene,1])
#a3 <- sum(chdgenes %in% our[1:ngene,1])

a <- matrix(c(a1,217,a2,1032),2,2)
fisher.test(a)

#phyper(203,231,21232-231,243,lower.tail = FALSE, log.p = FALSE)*243
#phyper(13,231,21232-231,243,lower.tail = FALSE, log.p = FALSE)*243

o1 <- length(intersect(chdresult[1:ngene,1],our[1:ngene,1]))
o2 <- length(intersect(chdother[1:ngene,1],our[1:ngene,1]))

a <- matrix(c(o1,217,o2,1032),2,2)
fisher.test(a)

o3 <- length(intersect(chdresult[1:ngene,1],chdother[1:ngene,1]))
phyper(o3-1,231,21232-231,231,lower.tail = FALSE, log.p = FALSE)



samid <- as.matrix(read.table("random_samples/CHDnature.txt"))
samid <- unique(samid)
length(samid)



### random mutations 
overlapg <- function(S1,S2,our,ngene){
    o <- 1:3
    o[1] <- length(intersect(S1[1:ngene,1],our[1:ngene,1]))
    o[2] <- length(intersect(S2[1:ngene,1],our[1:ngene,1]))
    o[3] <- length(intersect(S1[1:ngene,1],S2[1:ngene,1]))
    o
}

p2 <- as.matrix(read.table("result/P2_6_resultall.txt"))
t2 <- as.matrix(read.table("result/T2_6_resultall.txt"))
our <- read.table("result/PCGClog_6_resultall.txt")
ngene <- 300

o=overlapg(p2,t2,our,ngene)
a<- matrix(c(68125,68127,o[1:2]),2,2)
fisher.test(a)
o[3]

phyper(o[3]-1,300,21232-300,300,lower.tail = FALSE, log.p = FALSE)

p3 <- as.matrix(read.table("result/P3_6_resultall.txt"))
t3 <- as.matrix(read.table("result/T3_6_resultall.txt"))
our <- read.table("result/PCGClog_6_resultall.txt")
ngene <- 300
o=overlapg(p3,t3,our,ngene)
o[3]
a<- matrix(c(45417,90835,o[1:2]),2,2)
fisher.test(a)

phyper(o[3]-1,300,21232-300,300,lower.tail = FALSE, log.p = FALSE)



## total samples used
n <- dim(samples)[1]
ids <- c()
for(i in 1:n){
    b <- unlist(strsplit(gsub("[(.);]","\t",samples[i,"proband.ID.GT.AD.DP.GQ.PL."]),"\t"))
    ids <- union(ids,b[seq(1,length(b),3)])
}


}

randset_enri <- function(){
    source("enrichana.R")
    library(ROCR)
    TPcut=0.1
    TADAFile <- "ASD/TADAresult/TADAdenovo_ASDall.csv"
    TADAall <- read.csv(TADAFile)
    Tset <- TADAall[TADAall[,"qvalue.dn"]<= TPcut,1]
    dirstr <- "result/randresult/"
    netflag <- 6      
    DAWNnames <-  readLines(con <- file("result/randresult/TADArandset.txt","r"))   
    close(con)
    
    AUC <- list()
    AUC$TADA <- matrix(0,8,5);AUC$DAWN <- matrix(0,8,5);AUC$RWR2 <- matrix(0,8,5); AUC$MIS <- matrix(0,8,5);
    AUCfdr <- list()
    AUCfdr$TADA <- matrix(0,8,5);AUCfdr$DAWN <- matrix(0,8,5);AUCfdr$RWR2 <- matrix(0,8,5); AUCfdr$MIS <- matrix(0,8,5);
    AUCnum <- list()
    AUCnum$TADA <- matrix(0,8,5);AUCnum$DAWN <- matrix(0,8,5);AUCnum$RWR2 <- matrix(0,8,5); AUCnum$MIS <- matrix(0,8,5);
    fixSEN <- list()
    fixSEN$TADA <- matrix(0,8,5);fixSEN$DAWN <- matrix(0,8,5);fixSEN$RWR2 <- matrix(0,8,5); fixSEN$MIS <- matrix(0,8,5);
    comgene <- list()
    comgene$TADA <- matrix(0,8,5);comgene$DAWN <- matrix(0,8,5);comgene$RWR2 <- matrix(0,8,5); comgene$MIS <- matrix(0,8,5); 
    
    ## ALL results
    Fixspe=0.8
    FixFDR=0.3
    FixNUM=500
    for(j in 2:9){
        for(i in 1:5){
            ### TADA method
            TADAfile <- paste("ASD/TADAresult/TADAdenovo_part",j,"_",i,".csv",sep="")
            tmp <- one_auc(TADAfile,Tset,Fixspe=Fixspe,FixFDR=FixFDR, FixNUM=FixNUM,pl=FALSE)
            AUC$TADA[j-1,i] <- tmp$auc;AUCfdr$TADA[j-1,i] <- tmp$aucfdr;fixSEN$TADA[j-1,i] <- tmp$sensitivity;AUCnum$TADA[j-1,i] <- tmp$aucnum;
            tmpfile <- paste("ASD/TADAresult/TADAdenovo_rest",j,"_",i,".csv",sep="")
            tmp1 <- one_auc(tmpfile,Tset,Fixspe=Fixspe,FixFDR=FixFDR, FixNUM=FixNUM,pl=FALSE)
            comgene$TADA[j-1,i] <- length(intersect(tmp$geneset,tmp1$geneset))
            
            ### DAWN method
            k <- which(DAWNnames==paste("TADAdenovo_part",j,"_",i,".csv",sep=""))
            DAWNfile=paste(dirstr,"DAWN_randset",k,".csv",sep="")
            tmp <- one_auc(DAWNfile,Tset,Fixspe=Fixspe,FixFDR=FixFDR, FixNUM=FixNUM,pl=FALSE)
            AUC$DAWN[j-1,i] <- tmp$auc;AUCfdr$DAWN[j-1,i] <- tmp$aucfdr;fixSEN$DAWN[j-1,i] <- tmp$sensitivity;AUCnum$DAWN[j-1,i] <- tmp$aucnum;
            k <- which(DAWNnames==paste("TADAdenovo_rest",j,"_",i,".csv",sep=""))
            DAWNfile=paste(dirstr,"DAWN_randset",k,".csv",sep="")
            tmp1 <- one_auc(DAWNfile,Tset,Fixspe=Fixspe,FixFDR=FixFDR, FixNUM=FixNUM,pl=FALSE)
            comgene$DAWN[j-1,i] <- length(intersect(tmp$geneset,tmp1$geneset))
            
            ### RWR method
            hotnetfile <- paste(dirstr,"hotnetresult1part",j,"_",i,"6.txt",sep="")
            tmp <- one_auc(hotnetfile,Tset,Fixspe=Fixspe,FixFDR=FixFDR, FixNUM=FixNUM,pl=FALSE)
            AUC$RWR2[j-1,i] <- tmp$auc;AUCfdr$RWR2[j-1,i] <- tmp$aucfdr;fixSEN$RWR2[j-1,i] <- tmp$sensitivity;AUCnum$RWR2[j-1,i] <- tmp$aucnum;
            hotnetfile <- paste(dirstr,"hotnetresult1rest",j,"_",i,"6.txt",sep="")
            tmp1 <- one_auc(hotnetfile,Tset,Fixspe=Fixspe,FixFDR=FixFDR, FixNUM=FixNUM,pl=FALSE)
            comgene$RWR2[j-1,i] <- length(intersect(tmp$geneset,tmp1$geneset))
            
            ### modified IS model, shorted as, MIS
            MISfile <- paste(dirstr,"CRFresult_part",j,"_",i,"LBP_61.txt",sep="")
            tmp <- one_auc(MISfile,Tset,Fixspe=Fixspe,FixFDR=FixFDR, FixNUM=FixNUM,pl=FALSE)
            AUC$MIS[j-1,i] <- tmp$auc;AUCfdr$MIS[j-1,i] <- tmp$aucfdr;fixSEN$MIS[j-1,i] <- tmp$sensitivity;AUCnum$MIS[j-1,i] <- tmp$aucnum;
            MISfile <- paste(dirstr,"CRFresult_rest",j,"_",i,"LBP_61.txt",sep="")
            tmp1 <- one_auc(MISfile,Tset,Fixspe=Fixspe,FixFDR=FixFDR, FixNUM=FixNUM,pl=FALSE)
            comgene$MIS[j-1,i] <- length(intersect(tmp$geneset,tmp1$geneset))
        }
    }
    
    x <- seq(0.2,0.9,0.1)
    ## Figure 1: AUC values for eight random sets
    main="Average AUC based on different fractions of 3962 trios"
    xlab="The fractions of trios";ylab="AUC";
    plot_one(x,main,xlab,ylab,AUC,fig=1,nmeth=4)
    
    ## Figure 2: AUC values for eight random sets with a fix FDR
    main="Average AUC values based on FDR<=0.3"
    xlab="The fraction of samples";ylab="AUC"  
    plot_one(x,main,xlab,ylab,AUCfdr,fig=2,nmeth=4)
    

    main=paste("Average AUC values based on top ",FixNUM," genes",sep="")
    xlab="The fraction of samples";ylab="AUC" 
    plot_one(x,main,xlab,ylab,AUCnum,fig=2,nmeth=4)
    
    ## Figure 3: sensitivity for given specificity = 0.8
    main="Average sensitivities given specificity=0.8"
    xlab="The fraction of samples";ylab="Sensitivity"    
    plot_one(x,main,xlab,ylab,fixSEN,fig=2,nmeth=4)
    
    ## Figure 4: the fraction of recurrent genes in two complementary sets
    main="The fraction of recurrent genes"
    xlab="The fraction of samples";ylab="The number of recurrent genes"  
    plot_one(x,main,xlab,ylab,comgene,fig=2,nmeth=4)
    
}

randset2_enri <- function(){

    source("enrichana.R")
    library(ROCR)
    TPcut=0.2
    TADAFile <- "ASD/TADAresult/TADAdenovo_ASDall.csv"
    TADAall <- read.csv(TADAFile)
    Tset0 <- TADAall[TADAall[,"qvalue.dn"]<= TPcut,1]
    
    Tset <- Tset0
    i=1
    DAWNr <- read.csv(paste("result/leaveoneresult/DAWN_randset2_",i,".csv",sep=""))
    hotnetr <- read.table(paste("result/leaveoneresult/hotnetresult1rand2_",i,"6.txt",sep=""))
    MISr <- read.table(paste("result/leaveoneresult/CRFresult_rand2_",i,"LBP_61.txt",sep=""))
    Tset <- intersect(Tset,DAWNr[,1])
    Tset <- intersect(Tset,hotnetr[,1])
    Tset <- intersect(Tset,MISr[,1])
    
    ALLset <- TADAall[,1]
    ALLset <- intersect(ALLset,DAWNr[,1])
    ALLset <- intersect(ALLset,hotnetr[,1])
    ALLset <- intersect(ALLset,MISr[,1])
    
    subs <- match(Tset,Tset0)
    dirstr <- "result/leaveoneresult/"
    netflag <- 6      
    nmeth <- 4
    RankM <- matrix(0,length(Tset),nmeth) 
    FDRM <- matrix(0,length(Tset),nmeth)
    ## ALL results
    k=1
    for(i in subs){
        ### TADA method
        TADAfile <- paste("ASD/TADAresult/leaveone/TADAdenovo_rand2_",i,".csv",sep="")
        TADAr <- read.csv(TADAfile)
        TADAr <- TADAr[TADAr[,1]%in% ALLset,]
        RankM[k,1] <- which(TADAr[,1]==Tset0[i])
        FDRM[k,1] <- -log(TADAr[RankM[k,1],"qvalue.dn"])
        
        ### DAWN method
        DAWNr <- read.csv(paste("result/leaveoneresult/DAWN_randset2_",i,".csv",sep=""))
        DAWNr <- DAWNr[DAWNr[,1]%in% ALLset,]
        DAWNr[is.na(DAWNr[,"FDR"]),"FDR"] <- 1
        RankM[k,2] <- which(DAWNr[,"Gene"]==Tset0[i])
        FDRM[k,2] <- -log(DAWNr[RankM[k,2],"FDR"])            
        
        ### RWR method
        hotnetr <- read.table(paste("result/leaveoneresult/hotnetresult1rand2_",i,"6.txt",sep=""))
        hotnetr <- hotnetr[hotnetr[,1]%in% ALLset,]
        RankM[k,3] <- which(hotnetr[,1]==Tset0[i])
        FDRM[k,3] <- -log(hotnetr[RankM[k,3],5])           
        
        ### modified IS model, shorted as, MIS
        MISr <- read.table(paste("result/leaveoneresult/CRFresult_rand2_",i,"LBP_61.txt",sep=""))
        MISr <- MISr[MISr[,1]%in% ALLset,]
        RankM[k,4] <- ifelse(Tset0[i] %in% MISr[,1],which(MISr[,1]==Tset0[i]),0)
        FDRM[k,4] <-  ifelse(Tset0[i] %in% MISr[,1],-log(MISr[RankM[k,4],5]),1) 
        
        k <- k + 1
    }
    
    main="Gene rank with leave-one process"
    xlab="The original ranks";ylab="Prediction rank";
    library(vioplot)
    vioplot(RankM[,1],RankM[,2],RankM[,3],RankM[,4],names=c("TADA","DAWN","RWR","MIS"))
    #bar_err(RankM,main,xlab,ylab)
    
    fM <- RankM
    Tcut <- seq(1,length(ALLset),length.out = dim(RankM)[1])
    for(i in 1:nmeth){
        for(j in 1:dim(RankM)[1]){
            fM[j,i] <- sum(RankM[1:j,i] <= Tcut[j])
        }
    }
    main="The number of cumulative genes"
    xlab="Ranks";ylab="The number of genes";
    plot_two(Tcut,main,xlab,ylab,fM,pos="bottomright")
   
    FDRM[FDRM>1] <- 1
    main="Gene predicted FDR with leave-one process"
    xlab="The original ranks";ylab="Prediction FDR";
    plot_two(subs,main,xlab,ylab,FDRM,pos="bottomright")   
    
}

randset3_enri <- function(){
    
    source("enrichana.R")
    dirstr <- "result/randresult3/"
    DAWNnames <-  readLines(con <- file("result/randresult3/TADArandset3.txt","r"))   
    close(con)
    Maginames <-  readLines(con <- file("result/randresult3/MAGI/MAGIfilenames.txt","r"))   
    close(con)
    netstr <- c("STRING/","iRef/","coexp/","Infer/")
    netnum <- c(31,6,7,8)
    nnet=4;TPcut=0.1
    
    filenames <- c("ASD/TADAresult/TADAdenovo_ASDall.csv","result/randresult3/coexp/DAWN_randset1.csv","result/randresult3/MAGI/RandomGeneList.1","result/randresult3/STRING/hotnetresult1ASDall31.txt","result/randresult3/iRef/hotnetresult1ASDall6.txt","result/randresult3/coexp/hotnetresult1ASDall7.txt","result/randresult3/Infer/hotnetresult1ASDall8.txt","result/randresult3/STRING/CRFresult_ASDallLBP_31.txt","result/randresult3/iRef/CRFresult_ASDallLBP_6.txt","result/randresult3/coexp/CRFresult_ASDallLBP_7.txt","result/randresult3/Infer/CRFresult_ASDallLBP_8.txt")
    
    tmp <- gene_select(j=0,i=0,TPcut=0.1,nnet,datdir="randresult3/",tadadir="randset3/",Dawnfile="TADArandset3.txt",Magifile="MAGIfilenames.txt",filens=filenames)
    genes <- tmp$genes
    allP <- tmp$allP
    allN <- tmp$allN
    
    tmp <- new_auc(filenames,allP,genes,flag=2,iplot=TRUE);
    tmp <- read.csv("result/randresult3/coexp/DAWN_randset1.csv")
    geneDAWN <- intersect(genes,tmp[!is.na(tmp[,"FDR"]),"Gene"])
    tmp <- new_auc(filenames,allP,geneDAWN,flag=2,iplot=TRUE);
    
    AUC <- array(0,dim=c(8,10,3+nnet*2)) 
    FPR <- array(0,dim=c(8,10,3+nnet*2))
    TPR <- array(0,dim=c(8,10,3+nnet*2))
    
    DAWNflag <- 1
    datdir="randresult3/";
    tadadir="randset3/"
    for(j in 2:9){
        for(i in 1:10){    
            tmp <- gene_select(j=j,i=i,TPcut=0.1,nnet,datdir="randresult3/",tadadir="randset3/",Dawnfile="TADArandset3.txt")
            genes <- tmp$genes
            allP <- tmp$allP
            allN <- tmp$allN
            
            TADAfile <- paste("ASD/TADAresult/randset3/TADAdenovo_part",j,"_",i,".csv",sep="")
            k <- which(DAWNnames==paste("TADAdenovo_part",j,"_",i,".csv",sep=""))
            DAWNfile=paste(dirstr,"coexp/DAWN_randset",k,".csv",sep="")
            
            
            netstr <- c("STRING/","iRef/","coexp/","Infer/")
            netnum <- c(31,6,7,8)
            
            #dirstr <- paste("result/",datdir,sep="")
            ## all positive samples could be predicted by all methods
            #j=2;i=1;
            TADAfile <- paste("ASD/TADAresult/",tadadir,"TADAdenovo_part",j,"_",i,".csv",sep="")
            k <- which(DAWNnames==paste("TADAdenovo_part",j,"_",i,".csv",sep=""))
            DAWNfile=paste("result/",datdir,"coexp/DAWN_randset",k,".csv",sep="")
            
            k <- which(Maginames==paste("ASDrand",j,"depart",i,".csv",sep=""))
            MAGIfile=paste("result/",datdir,"MAGI/RandomGeneList.",k,sep="")
            
            
            hfiles <- 1:nnet
            mfiles <- 1:nnet
            for(k in 1:nnet){
                hfiles[k] <- paste(dirstr,netstr[k],"hotnetresult1part",j,"_",i,netnum[k],".txt",sep="")
                mfiles[k] <- paste(dirstr,netstr[k],"CRFresult_part",j,"_",i,"LBP_",netnum[k],".txt",sep="")
            }
            filenames <- c(TADAfile,DAWNfile,MAGIfile,hfiles,mfiles)
            if(DAWNflag==1){
                tmp <- read.csv(DAWNfile)
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
    
    save(AUC,file="AUC_m")
    save(TPR,file="TPR_m")
    save(FPR,file="FRP_m")

    ## plot auc values
    ## Figure 1: AUC values for eight random sets
    
    pdf(file="plot/AUCs.pdf",width=10,height=8)
    load("AUC")
    cols <- c("Black","Blue","Cyan","Brown","Coral","DarkViolet","Green","HotPink","Orange","Plum","deeppink")
    main="AUC performance for different trios"
    xlab="The number of trios";ylab="AUC";
    ymint <- 1;ymaxt <- 0;
    ymint <- min(AUC);ymaxt <- max(AUC);
    ymin=max(min(ymint)+0.2,0);ymax = min(ymaxt+0.1,1);
    x <- seq(0.2,0.9,0.1)
    nmeth <- length(filenames)
    Y1 <- rep(0,nmeth)
    for(i in 1:nmeth){
        y <- rowMeans(AUC[,,i])
        sd <- apply( AUC[,,i], 1, sd)/2 
        plot_err(x,y,sd,i,main,xlab,ylab,ylim=c(ymin,ymax),xlim=c(-0.05,0.95))
        Y1[i] <- y[1]
    }
    Y1[6] <- 0.86
    Y1[10] <- 0.87
    Y1[1] <- 0.745
    Y1[3] <- 0.82
    Y1[4] <- 0.775
    Y1[5] <- 0.785
    
    text(rep(0.2,length(filenames)),Y1,labels=c("TADA","DAWN","MAGI","RWR-STRING","RWR-iRef","RWR-coexp","RWR-Infer","MIS-STRING","MIS-iRef","MIS-coexp","MIS-Infer"),pos=2,col=cols,cex=1.5)
    dev.off()
    
    legend("bottomright",col=cols[1:nmeth],legend=c("TADA","DAWN","MAGI","RWR-STRING","RWR-iRef","RWR-coexp","RWR-Infer","MIS-STRING","MIS-iRef","MIS-coexp","MIS-Infer"),lty=rep(1,nmeth),lwd=rep(2,nmeth),cex=1.2,y.intersp=0.8)

}

plot_TPR <- function(TPR,nmeth=10){
    
    cols <- c("Black","Blue","Cyan","Brown","Coral","DarkViolet","Green","HotPink","Orange","Plum")
    main="True positive rate with fixed FPR=0.1"
    xlab="The number of trios";ylab="TPR";
    
    ymint <- 1; ymaxt <- 0;
    ymint <- min(TPR);ymaxt <- max(TPR);
    ymin=max(min(ymint)+0.2,0);ymax = min(ymaxt+0.1,1);
    xlim=c(-0.05,0.95)
    ylim=c(ymin,ymax)
    
    x <- seq(0.2,0.9,0.1)
    #nmeth <- length(filenames)
    Y1 <- rep(0,nmeth)
    for(i in 1:nmeth){
        y <- rowMeans(TPR[,,i])
        sd <- apply( TPR[,,i], 1, sd)/2     
        if(i==1){
            par(mai=c(2,2,1,1))
            plot(x, y,xaxt="n", type="b",col=cols[i], main=main, ylab=ylab,xlab=xlab,xlim=xlim,ylim=ylim,cex.lab=2,cex.main=2,lwd=2,cex.axis=1.8)
            axis(1,at=seq(0.2,0.9,0.1),labels=floor(3953*seq(0.2,0.9,0.1)),font=1,cex.axis=1.5)
        }else{
            lines(x, y, type="b",col=cols[i],xaxt="n")
        }
        Y1[i] <- y[1]
    }
    Y1[6] <- 0.8
    Y1[5] <- 0.79
    Y1[7] <- 0.74
    Y1[8] <- 0.73
    text(rep(0.2,10),Y1,labels=c("TADA","DAWN","RWR-STRING","RWR-iRef","RWR-coexp","RWR-Infer","MIS-STRING","MIS-iRef","MIS-coexp","MIS-Infer"),pos=2,col=cols,cex=1)
    
    legend("bottomright",col=cols[1:nmeth],legend=c("TADA","DAWN","RWR-STRING","RWR-iRef","RWR-coexp","RWR-Infer","MIS-STRING","MIS-iRef","MIS-coexp","MIS-Infer"),lty=rep(1,nmeth),lwd=rep(2,nmeth),cex=1.2,y.intersp=0.8)

}

randset4_enri <- function(){
    source("enrichana.R")
    library(ROCR)
    TPcut=0.2
    TADAFile <- "ASD/TADAresult/TADAdenovo_ASDall.csv"
    TADAall <- read.csv(TADAFile)
    Tset <- TADAall[TADAall[,"qvalue.dn"]<= TPcut,1]
    Tset0 <- Tset
    dirstr <- "result/leaveone1result/"
    netstr <- c("STRING/","iRef/","coexp/","Infer/")
    netnum <- c(31,6,7,8)
    
    ## all positive samples could be predicted by all methods
    i=1;
    TADAfile <- paste("ASD/TADAresult/leaveone1/TADAdenovo_rand2_1_",i,".csv",sep="")
    DAWNfile <- paste("result/leaveone1result/coexp/DAWN_randset2_1_",i,".csv",sep="")
    hfiles <- 1:2
    mfiles <- 1:2
    for(k in 1:2){
        hfiles[k] <- paste(dirstr,netstr[k],"hotnetresult1rand2_1_",i,netnum[k],".txt",sep="")
        mfiles[k] <- paste(dirstr,netstr[k],"CRFresult_rand2_1_",i,"LBP_",netnum[k],".txt",sep="")
    }
    filenames <- c(TADAfile,DAWNfile,hfiles,mfiles)
    allP <- TADAall[,1]
    for(k in 1:length(filenames)){
        filename <- filenames[k]
        if(grepl(".csv",filename)){ result <- read.csv(filename); 
        }else if(grepl(".txt",filename)){result <- read.table(filename);}
        allP <- intersect(allP,result[,1])
        print(length(allP))
    } 
    Tset <- intersect(Tset,allP)
    genes <- allP
    
    RANK <- array(0,dim=c(length(Tset),2,length(filenames)))
    kf = 1
    subs <- match(Tset,Tset0)
    for(i in subs){
        TADAfile <- paste("ASD/TADAresult/leaveone1/TADAdenovo_rand2_1_",i,".csv",sep="")
        DAWNfile <- paste("result/leaveone1result/coexp/DAWN_randset2_1_",i,".csv",sep="")
        hfiles <- 1:2
        mfiles <- 1:2
        for(k in 1:2){
            hfiles[k] <- paste(dirstr,netstr[k],"hotnetresult1rand2_1_",i,netnum[k],".txt",sep="")
            mfiles[k] <- paste(dirstr,netstr[k],"CRFresult_rand2_1_",i,"LBP_",netnum[k],".txt",sep="")
        }
        filenames <- c(TADAfile,DAWNfile,hfiles,mfiles)
        tmp <- new_rank(filenames,Tset0[i],genes)
        RANK[kf,,] <- t(tmp)
        kf <- kf + 1
    }

    ## plot rank and FDR distributions
    library(caroline)
    legend=c("TADA","DAWN","RWR-STRING","RWR-iRef","MIS-STRING","MIS-iRef")  
    tmp1 <- data.frame(RANK[,1,])
    par(mai=c(2,2,2,2),cex.lab=2,cex.main=2)
    violins(tmp1,names=legend,col=c('purple','lightblue','lightgreen','red','orange','yellow'), main="Violin plot on rank",ylab="Rank",deciles=FALSE,connect="")
    
    par(mai=c(2,2,2,2),cex.lab=2,cex.main=2)
    tmp2 <- data.frame(RANK[,2,])
    violins(tmp2,names=legend,col=c('purple','lightblue','lightgreen','red','orange','yellow'), main="Violin plot on FDR",ylab="FDR",deciles=FALSE,connect="")
    
}

randset5_enri <- function(){
    source("enrichana.R")
    library(ROCR)
    TPcut=0.2
    tmp <- gene_select(TPcut=0.2,nnet=2)
    allgene <- tmp$genes
    TADAFile <- "ASD/TADAresult/TADAdenovo_ASDall.csv"
    TADAall <- read.csv(TADAFile)
    Tset <- TADAall[TADAall[,"qvalue.dn"]<= TPcut,1]
    
    dirstr <- "result/randresult/"
    DAWNnames <-  readLines(con <- file("result/randresult/TADArandset.txt","r"))   
    close(con)
    netstr <- c("STRING/","iRef/","coexp/","Infer/")
    netnum <- c(31,6,7,8)
    
    flag=1
    AUC <- array(0,dim=c(8,5,6))    
    for(i in 1:5){
        partL <- list()
        restL <- list()
        for(j in 2:9){       
            mutafile1 <- paste("ASD/randset/ASDrand",j,"depart",i,".csv",sep="")
            mutafile2 <- paste("ASD/randset/ASDrand",j,"denrest",i,".csv",sep="")
            muta1 <- read.csv(mutafile1)
            muta2 <- read.csv(mutafile2)
            genes <- intersect(Tset,intersect(muta1[rowSums(muta1[,3:4])>0,"Gene"], muta2[rowSums(muta2[,3:4])>0,"Gene"]))
            
            filenames <- matrix(,6,2)
            filenames[1,1] <- paste("ASD/TADAresult/randset/TADAdenovo_part",j,"_",i,".csv",sep="")
            filenames[1,2] <- paste("ASD/TADAresult/randset/TADAdenovo_rest",j,"_",i,".csv",sep="")
            
            k <- which(DAWNnames==paste("TADAdenovo_part",j,"_",i,".csv",sep=""))
            filenames[2,1] <- paste(dirstr,"coexp/DAWN_randset",k,".csv",sep="")
            k <- which(DAWNnames==paste("TADAdenovo_rest",j,"_",i,".csv",sep=""))
            filenames[2,2] <- paste(dirstr,"coexp/DAWN_randset",k,".csv",sep="")

            for(k in 1:2){
                filenames[3+k-1,1] <- paste(dirstr,netstr[k],"hotnetresult1part",j,"_",i,netnum[k],".txt",sep="")
                filenames[3+k-1,2] <- paste(dirstr,netstr[k],"hotnetresult1rest",j,"_",i,netnum[k],".txt",sep="")
                filenames[5+k-1,1] <- paste(dirstr,netstr[k],"CRFresult_part",j,"_",i,"LBP_",netnum[k],".txt",sep="")
                filenames[5+k-1,2] <- paste(dirstr,netstr[k],"CRFresult_rest",j,"_",i,"LBP_",netnum[k],".txt",sep="")
            }
            genes <- intersect(genes,allgene)
            tmp <- new_recurrent(filenames,genes,allgene,flag=flag)
            partL[[j-1]] <- tmp[,1,]
            restL[[j-1]] <- tmp[,2,]
        }
        main="Violin plot"; xlab=""; if(flag==1){ylab="Rank";}else{ylab="FDR";}
        new_vioplot(partL,restL,main,xlab,ylab)
    }

}

randset6_enri <- function(){
    options(stringsAsFactors=FALSE)
    source("enrichana.R")
    library(ROCR)
    TPcut=0.1
    TADAFile <- "ASD/TADAresult/TADAdenovo_ASDall.csv"
    TADAall <- read.csv(TADAFile)
    Tset <- TADAall[TADAall[,"qvalue.dn"]<= TPcut,1]
    
    dirstr <- "result/randresult2/"
    DAWNnames <-  readLines(con <- file("result/randresult2/TADArandset2.txt","r"))   
    close(con)
    netstr <- c("STRING/","iRef/","coexp/","Infer/")
    netnum <- c(31,6,7,8)
    
    TPcut = 0.1;
    nnet=4;
    partL <- list()
    m <- 1
    
    sM <- matrix(0,10,10)
    cutoff <- c(0.09937977,0.05662114,0.3968514,0.3616888,0.4189705,0.4202216,0.6476619,0.6398611,0.6978669,0.6823814)
    sML <- array(0,dim=c(10,10,5))
    for(i in 1:5){
        for(j in 3){
            onefile <- paste("ASD/TADAresult/randset2/TADAdenovo_ASD2sub",j,"denrest",i,".csv",sep="")
            oneresult <- read.csv(onefile)
            
            filenames <- 1:2
            filenames[1] <- paste("ASD/TADAresult/randset2/TADAdenovo_ASD1sub",j,"depart",i,".csv",sep="")
            subresult <- read.csv(filenames[1])
            
            k <- which(DAWNnames==paste("TADAdenovo_ASD1sub",j,"depart",i,".csv",sep=""))
            filenames[2] <- paste(dirstr,"coexp/DAWN_randset2_3_",k,".csv",sep="")
            hfiles <- 1:nnet
            mfiles <- 1:nnet
            for(k in 1:nnet){
                hfiles[k] <- paste(dirstr,netstr[k],"hotnetresult1part",j,"_",i,netnum[k],".txt",sep="")
                mfiles[k] <- paste(dirstr,netstr[k],"CRFresult_part",j,"_",i,"LBP_",netnum[k],".txt",sep="")
            }
            filenames <- c(filenames,hfiles,mfiles)
            
            result <- list()
            for(k in 1:length(filenames)){
                
                if(k==2){
                    DAWNr <- read.csv(filenames[2])
                    subs <- is.na(DAWNr[,"FDR"])
                    result[[2]] <- DAWNr[!subs,]
                    genes <- result[[2]][result[[2]][,"FDR"] <= TPcut,1]
                }else{
                filename <- filenames[k]
                if(grepl(".csv",filename)){ tmp <- read.csv(filename); 
                }else if(grepl(".txt",filename)){tmp <- read.table(filename);}
                result[[k]] <- tmp[!duplicated(tmp[,1]),]
                if(grepl(".csv",filename)){
                    genes <- result[[k]][result[[k]][,"qvalue.dn"] <= cutoff[k],1] ## TPcut,1]
                }else{ genes <- result[[k]][result[[k]][,2] >= cutoff[k],1] }##  !!!<= TPcut,1] }
                }
                
                partL[[k]] <- oneresult[match(genes,oneresult[,1]),"qvalue.dn"]
                partL[[k]][is.na(partL[[k]])] <- 1
                oneg1 <- oneresult[oneresult[,"qvalue.dn"] <=0.1,1]
                oneg2 <- oneresult[oneresult[,"qvalue.dn"] <=0.2,1]
                tmpgenes <- intersect(genes,oneg1)
                sM[k,1] <- length(intersect(genes,oneg1))
                sM[k,2] <- length(intersect(genes,oneg2))
                sM[k,3] <- phyper(sM[k,1]-1, length(intersect(oneg1,result[[k]][,1])), dim(result[[k]])[1] - length(intersect(oneg1,result[[k]][,1])), length(genes), lower.tail = FALSE, log.p = FALSE)
                sM[k,4] <- phyper(sM[k,2]-1, length(intersect(oneg2,result[[k]][,1])), dim(result[[k]])[1] - length(intersect(oneg2,result[[k]][,1])), length(genes), lower.tail = FALSE, log.p = FALSE)
                sM[k,5] <- length(genes)
                sM[k,6] <- sum(oneresult[match(tmpgenes,oneresult[,1]),3:4])
                sM[k,7] <- sum(oneresult[match(tmpgenes,oneresult[,1]),3:4])- sum(subresult[match(tmpgenes,subresult[,1]),3:4])
                sM[k,8] <- (sum(oneresult[match(tmpgenes,oneresult[,1]),3:4])- sum(subresult[match(tmpgenes,subresult[,1]),3:4]))/length(tmpgenes)
                sM[k,9] <- (sum(oneresult[match(genes,oneresult[,1]),3:4])- sum(subresult[match(genes,subresult[,1]),3:4]))/length(genes)
                sM[k,10] <- sum(rowSums(oneresult[match(genes,oneresult[,1]),3:4]- subresult[match(genes,subresult[,1]),3:4])>0)/length(genes)
            }
        }
        
        sML[,,i] <- sM
    }
    
    k=10
    a <- cbind(sML[,k,1],sML[,k,2],sML[,k,3],sML[,k,4],sML[,k,5]) 
    pdf(file="plot/subsams.pdf",width=10,height=7)
    legend=c("TADA","DAWN","RWR-STRING","RWR-iRef","RWR-coexp","RWR-Infer","MIS-STRING","MIS-iRef","MIS-coexp","MIS-Infer")  
    cols <- c("Black","Blue","Cyan","Brown","Coral","DarkViolet","Green","HotPink","Orange","Plum")
    par(mai=c(2,1,1,1))
    boxplot(t(a),col=cols,names=legend,las=2,ylab="The fraction of genes",cex.axis=1.5,cex.lab=2,cex.names=2)
    
    dev.off()
    
    
    
    pdf(file="plot/subsams.pdf",width=10,height=10)
    legend=c("TADA","DAWN","RWR-STRING","RWR-iRef","RWR-coexp","RWR-Infer","MIS-STRING","MIS-iRef","MIS-coexp","MIS-Infer")  
    cols <- c("Black","Blue","Cyan","Brown","Coral","DarkViolet","Green","HotPink","Orange","Plum")
    par(mai=c(3,1,1,1))
    Y <- sM[,6]#(sM[,1]/sM[,5])/(sM[1,1]/sM[1,5])
    barplot(Y,col=cols,names.arg=legend,beside=T,las=2,ylim=c(0,max(1.2*Y)),ylab="De novo mutations",cex.axis=1.2,cex.lab=1.5,cex.names=1.5)
    #barplot(t(a), names.arg=legend, beside=T, ylab="The number of genes (FDR<=0.1)", 
            #cex.names=1.5, las=2, ylim=c(0,max(a)), col=c("darkblue","red"),cex.lab=1.5,cex.axis=1.25)
    #abline(h=1,lty=2)
    #text(seq(0.8,12,1.2)-0.1,Y,sM[,1],pos=3,cex=2)
    #text(seq(0.8,12,1.2)-0.1,sM[,1],c("**","","***","**","***","***","***","***","***","***"),pos=3,cex=2)
    #text(rep(0,3),c(31,33,35),c("*    p < 0.01","**   p < 1e-30","***  p < 1e-40"),pos=4,cex=1.2)
    dev.off()
    
    
    library(caroline)
    plot.new()
    main="Violin plot for rank distribution"; xlab=""; ylab="Rank";
    legend=c("TADA","DAWN","RWR-STRING","RWR-iRef","MIS-STRING","MIS-iRef")  
    legend <- rep(legend,3)
    par(bg=NA,mai=c(4,2,2,2))
    violins(partL,names=legend,col=rep(c('purple','lightblue','lightgreen','red','orange','yellow'),3), main=main,xlab=xlab,ylab=ylab,deciles=FALSE,connect="",at=c(1:6,8:13,15:20))
    abline(v=7,col=1,lwd=3)
    abline(v=14,col=1,lwd=3)
    mtext(text=c("395 -> 1186 trios","791 -> 2372 trios","1186 -> 3558 trios"), side=3, font=2, line=0.5,at=c(3.5,10.5,17.5))
    
}

randset4_boxplot_4_15 <- function(){
    options(stringsAsFactors=FALSE)
    source("enrichana.R")
    dirstr <- "result/randresult4/"
    netstr <- c("STRING/","iRef/","coexp/","Infer/")
    netnum <- c(31,6,7,8)
    
    filenames <- c("ASD/TADAresult/randset4/PtestASDall.csv","result/randresult4/STRING/CRFresult_ASD_4_16LBP_31.txt","result/randresult4/iRef/CRFresult_ASD_4_16LBP_6.txt","result/randresult4/coexp/CRFresult_ASD_4_16LBP_7.txt","result/randresult4/Infer/CRFresult_ASD_4_16LBP_8.txt")
    asdall <- read.csv(filenames[1])
    Tset <- asdall[asdall[,"min_p"]<=0.001,1]
    
    cutoff <- total_genes_cutoff(filenames,Tset,muts="score")
    #cutoff[1] <- 1 - cutoff[1]
    #cutoff <- c(0.09937977,0.05662114,0.3968514,0.3616888,0.4189705,0.4202216,0.6476619,0.6398611,0.6978669,0.6823814)
    sM <- matrix(-1,20,5)
    nnet <- 4
    for(i in 1:20){
        for(j in 3){
            onefile <- paste("ASD/TADAresult/randset4/Ptest_",j,"derest",i,".csv",sep="")
            oneresult <- read.csv(onefile)
            TADAfile <- paste("ASD/TADAresult/randset4/Ptest_",j,"depart",i,".csv",sep="")
            subresult <- read.csv(TADAfile)
            
            mfiles <- 1:nnet
            for(k in 1:nnet){
                mfiles[k] <- paste(dirstr,netstr[k],"CRFresult_part",j,"_",i,"LBP_",netnum[k],".txt",sep="")
            }
            filenames <- c(TADAfile,mfiles)
            
            result <- list()
            for(k in 1:length(filenames)){
                    filename <- filenames[k]
                    if(grepl(".csv",filename)){ tmp <- read.csv(filename); 
                    }else if(grepl(".txt",filename)){tmp <- read.table(filename);}
                    result[[k]] <- tmp[!duplicated(tmp[,1]),]
                    
                    if(grepl(".csv",filename)){
                        genes <- result[[k]][result[[k]][,"score"] >= cutoff[k],1] ## TPcut,1]
                    }else{ genes <- result[[k]][result[[k]][,2] >= cutoff[k],1] }##  !!!<= TPcut,1] }
                
                    #oneg1 <- oneresult[oneresult[,"qvalue.dn"] <=0.1,1]
                    sM[i,k] <- sum(rowSums(oneresult[match(genes,oneresult[,1]),c("dn.LoF","dn.mis3","dn.mis")]- subresult[match(genes,subresult[,1]),c("dn.LoF","dn.mis3","dn.mis")])>0)/length(genes)
            }
        }
    }

    subs <- !is.na(sM[,1])
    sM <- sM[subs,]
    #pdf(file="plot/subsams.pdf",width=10,height=7)
    legend=c("TADA","MIS-STRING","MIS-iRef","MIS-coexp","MIS-Infer")  
    cols <- c("White","Blue","Cyan","Brown","Coral")
    par(mai=c(2,1,1,1))
    boxplot(sM,col=cols,names=legend,las=2,ylab="The fraction of genes",cex.axis=1.5,cex.lab=2,cex.names=2)
    #dev.off()   
}

randset7_enri <- function(){
    source("enrichana.R")
    library(ROCR)
    TPcut=0.2
    TADAFile <- "ASD/TADAresult/TADAdenovo_ASDall.csv"
    TADAall <- read.csv(TADAFile)
    Tset <- TADAall[TADAall[,"qvalue.dn"]<= TPcut,1]
    Tset0 <- Tset
    dirstr <- "result/leaveone2result/"
    netstr <- c("STRING/","iRef/","coexp/","Infer/")
    netnum <- c(31,6,7,8)
    
    ## all positive samples could be predicted by all methods
    i=1;
    TADAfile <- paste("ASD/TADAresult/leaveone2/TADAdenovo_rand2_de1_",i,".csv",sep="")
    DAWNfile <- paste("result/leaveone2result/coexp/DAWN_randset2_de1_",i,".csv",sep="")
    hfiles <- 1:2
    mfiles <- 1:2
    for(k in 1:2){
        hfiles[k] <- paste(dirstr,netstr[k],"hotnetresult1rand2_de1_",i,netnum[k],".txt",sep="")
        mfiles[k] <- paste(dirstr,netstr[k],"CRFresult_rand2_de1_",i,"LBP_",netnum[k],".txt",sep="")
    }
    filenames <- c(TADAfile,DAWNfile,hfiles,mfiles)
    allP <- TADAall[,1]
    for(k in 1:length(filenames)){
        filename <- filenames[k]
        if(grepl(".csv",filename)){ result <- read.csv(filename); 
        }else if(grepl(".txt",filename)){result <- read.table(filename);}
        allP <- intersect(allP,result[,1])
        print(length(allP))
    } 
    Tset <- intersect(Tset,allP)
    genes <- allP
    
    RANK <- array(0,dim=c(length(Tset),2,length(filenames)))
    kf = 1
    subs <- match(Tset,Tset0)
    for(i in subs){
        TADAfile <- paste("ASD/TADAresult/leaveone2/TADAdenovo_rand2_de1_",i,".csv",sep="")
        DAWNfile <- paste("result/leaveone2result/coexp/DAWN_randset2_de1_",i,".csv",sep="")
        hfiles <- 1:2
        mfiles <- 1:2
        for(k in 1:2){
            hfiles[k] <- paste(dirstr,netstr[k],"hotnetresult1rand2_de1_",i,netnum[k],".txt",sep="")
            mfiles[k] <- paste(dirstr,netstr[k],"CRFresult_rand2_de1_",i,"LBP_",netnum[k],".txt",sep="")
        }
        filenames <- c(TADAfile,DAWNfile,hfiles,mfiles)
        tmp <- new_rank(filenames,Tset0[i],genes)
        RANK[kf,,] <- t(tmp)
        kf <- kf + 1
    }
    
    ## plot rank and FDR distributions
    library(caroline)
    plot.new()
    par(bg=NA,mai=c(5,2,2,2))
    legend=c("TADA","DAWN","RWR-STRING","RWR-iRef","MIS-STRING","MIS-iRef")  
    tmp1 <- data.frame(RANK[,1,])
    #tmp1[tmp1>2000] <- 2000
    par(mai=c(2,2,2,2),cex.lab=2,cex.main=2)
    violins(tmp1,names=legend,col=c('purple','lightblue','lightgreen','red','orange','yellow'), main="Violin plot on rank",ylab="Rank",deciles=FALSE,connect="")
    plot.new()
    par(mai=c(2,2,2,2),cex.lab=2,cex.main=2)
    tmp2 <- data.frame(RANK[,2,])
    violins(tmp2,names=legend,col=c('purple','lightblue','lightgreen','red','orange','yellow'), main="Violin plot on FDR",ylab="FDR",deciles=FALSE,connect="")
    
}

randset8_enri <- function(){
    source("enrichana.R")
    library(ROCR)
    TPcut=0.3;
    nnet=4;
    TADAFile <- "ASD/TADAresult/TADAdenovo_ASDall.csv"
    TADAall <- read.csv(TADAFile)
    Tset <- TADAall[TADAall[,"qvalue.dn"]<= TPcut,1]
    Tset0 <- Tset
    dirstr <- "result/leaveone3result/"
    netstr <- c("STRING/","iRef/","coexp/","Infer/")
    netnum <- c(31,6,7,8)
    
    ## all positive samples could be predicted by all methods
    i=1;
    TADAfile <- paste("ASD/TADAresult/leaveone3/TADAdenovo_rand2_de1_",i,".csv",sep="")
    DAWNfile <- paste("result/leaveone3result/coexp/DAWN_leaveone3_",i,".csv",sep="")
    hfiles <- 1:nnet
    mfiles <- 1:nnet
    for(k in 1:nnet){
        hfiles[k] <- paste(dirstr,netstr[k],"hotnetresult1rand2_de1_",i,netnum[k],".txt",sep="")
        mfiles[k] <- paste(dirstr,netstr[k],"CRFresult_rand2_de1_",i,"LBP_",netnum[k],".txt",sep="")
    }
    filenames <- c(TADAfile,DAWNfile,hfiles,mfiles)
    allP <- TADAall[,1]
    for(k in 1:length(filenames)){
        filename <- filenames[k]
        if(grepl(".csv",filename)){ result <- read.csv(filename); 
        }else if(grepl(".txt",filename)){result <- read.table(filename);}
        allP <- intersect(allP,result[,1])
        print(length(allP))
    } 
    Tset <- intersect(Tset,allP)
    genes <- allP
    
    RANK <- array(0,dim=c(length(Tset),2,length(filenames)))
    kf = 1
    subs <- match(Tset,Tset0)
    for(i in subs){
        TADAfile <- paste("ASD/TADAresult/leaveone3/TADAdenovo_rand2_de1_",i,".csv",sep="")
        DAWNfile <- paste("result/leaveone3result/coexp/DAWN_leaveone3_",i,".csv",sep="")
        hfiles <- 1:nnet
        mfiles <- 1:nnet
        for(k in 1:nnet){
            hfiles[k] <- paste(dirstr,netstr[k],"hotnetresult1rand2_de1_",i,netnum[k],".txt",sep="")
            mfiles[k] <- paste(dirstr,netstr[k],"CRFresult_rand2_de1_",i,"LBP_",netnum[k],".txt",sep="")
        }
        filenames <- c(TADAfile,DAWNfile,hfiles,mfiles)
        tmp <- new_rank(filenames,Tset0[i],genes)
        RANK[kf,,] <- t(tmp)
        kf <- kf + 1
    }
    
    ## plot rank and FDR distributions
#     plot.new()
#     pchs <- c(4,2,6,24,25,5,23,0,22)
#     for(i in 2:length(filenames)){
#         if(i==2){
#             plot(RANK[,2,1],RANK[,2,i],type="p",pch=pchs[i-1])
#         }else{
#             lines(RANK[,2,1],RANK[,2,i],type="p",pch=pchs[i-1])
#         }
#     }
#     lines(c(0,1),c(0,1),col=1,lty=2)
    
    par(bg=NA)
    pdf(file="plot/FDRs.pdf",width=12,height=12)
    load("RANK")
    cols <- c("Black","Blue","Cyan","Brown","Coral","DarkViolet","Green","HotPink","Orange","Plum")
    nmeth=length(filenames)-1
    subs <- c(1:3,5:length(filenames))
    Y <- rep(0,10)
    X <- rep(0,10)
    legend=c("TADA","DAWN","RWR-STRING","RWR-iRef","RWR-coexp","RWR-Infer","MIS-STRING","MIS-iRef","MIS-coexp","MIS-Infer")  
    for(i in subs){
        if(i==1){
            par(mar=c(4,6,4,2))
            plot(density(RANK[,2,i]),type="l",col=cols[i],lwd=2,xlim=c(0,0.5),xaxs="i",yaxs="i",xlab="FDR",ylab="Density",main="FDR distribution for risk genes",cex.lab=2,cex.main=2)
            #axis(2, at=seq(0, 0.5, by=0.1), labels = FALSE)
            #text(y = seq(0, 0.5, by=0.1), par("usr")[1], labels = seq(0, 0.5, by=0.1), srt = 0, pos = 2, xpd = TRUE)
        }else{
            lines(density(RANK[,2,i]),type="l",col=cols[i],lwd=1)
        }
        Y[i] <- max(density(RANK[,2,i])$y[density(RANK[,2,i])$x <= 0.2])
        X[i] <- density(RANK[,2,i])$x[which.max(density(RANK[,2,i])$y[density(RANK[,2,i])$x <= 0.2])]
    }
    X[3] <- 0.08
    X[5] <- 0.09
    X[6] <- 0.07
    X[7] <- 0.05
    X[8] <- 0.04
    X[9] <- 0.04
    X[10] <- 0.04
    Y[7] <- 1.2
    Y[3] <- 0.65
    Y[5] <- 0.8
    Y[6] <- 0.94
    Y[9] <- 1.5
    Y[10] <- 1.35
    abline(v=0.1,lty=2)
    abline(v=0.2,lty=2)
    abline(v=0.3,lty=2)
    text(X[subs],Y[subs],labels=legend[subs],col=cols[subs],cex=2)
    dev.off()    

#     legend("bottomright",col=cols[subs],legend=c("TADA","DAWN","RWR-STRING","RWR-coexp","RWR-Infer","MIS-STRING","MIS-iRef","MIS-coexp","MIS-Infer"),lty=rep(1,nmeth),lwd=rep(2,nmeth),cex=1.2,y.intersp=0.8)
#     
#     library(caroline)
#     plot.new()
#     par(bg=NA,mai=c(5,2,2,2))
#     legend=c("TADA","DAWN","RWR-STRING","RWR-iRef","RWR-coexp","RWR-Infer","MIS-STRING","MIS-iRef","MIS-coexp","MIS-Infer")  
#     cols <- c("Black","Blue","Cyan","Brown","Coral","DarkViolet","Green","HotPink","Orange","Plum")
#     tmp1 <- data.frame(RANK[,2,])
#     par(mai=c(2,2,1,1),cex.lab=1.2,cex.main=1.2)
#     violins(tmp1,names=legend,col=cols, main="Violin plot",ylab="FDR",deciles=FALSE,connect="")
   
    
}

leave_one_control <- function(){

    source("enrichana.R")
    library(ROCR)
    TPcut=0.3;
    nnet=4;
    TADAFile <- "ASD/TADAresult/TADAdenovo_ASDall.csv"
    TADAall <- read.csv(TADAFile)
    Tset <- TADAall[TADAall[,"qvalue.dn"]<= TPcut,1]
    Tset0 <- Tset
    dirstr <- "result/leaveone3result/"
    netstr <- c("STRING/","iRef/","coexp/","Infer/")
    netnum <- c(31,6,7,8)
    
    ## all positive samples could be predicted by all methods
    i=1;
    TADAfile <- paste("ASD/TADAresult/leaveone3/TADAdenovo_rand2_de1_",i,".csv",sep="")
    DAWNfile <- paste("result/leaveone3result/coexp/DAWN_leaveone3_",i,".csv",sep="")
    hfiles <- 1:nnet
    mfiles <- 1:nnet
    for(k in 1:nnet){
        hfiles[k] <- paste(dirstr,netstr[k],"hotnetresult1rand2_de1_",i,netnum[k],".txt",sep="")
        mfiles[k] <- paste(dirstr,netstr[k],"CRFresult_rand2_de1_",i,"LBP_",netnum[k],".txt",sep="")
    }
    filenames <- c(TADAfile,DAWNfile,hfiles,mfiles)
    allP <- TADAall[,1]
    for(k in 1:length(filenames)){
        filename <- filenames[k]
        if(grepl(".csv",filename)){ result <- read.csv(filename); 
        }else if(grepl(".txt",filename)){result <- read.table(filename);}
        allP <- intersect(allP,result[,1])
        print(length(allP))
    } 
    Tset <- intersect(Tset,allP)
    genes <- allP
    
    RANK <- array(0,dim=c(length(Tset),3,length(filenames)))
    kf = 1
    subs <- match(Tset,Tset0)
    for(i in subs){
        TADAfile <- paste("ASD/TADAresult/leaveone3/TADAdenovo_rand2_de1_",i,".csv",sep="")
        DAWNfile <- paste("result/leaveone3result/coexp/DAWN_leaveone3_",i,".csv",sep="")
        hfiles <- 1:nnet
        mfiles <- 1:nnet
        for(k in 1:nnet){
            hfiles[k] <- paste(dirstr,netstr[k],"hotnetresult1rand2_de1_",i,netnum[k],".txt",sep="")
            mfiles[k] <- paste(dirstr,netstr[k],"CRFresult_rand2_de1_",i,"LBP_",netnum[k],".txt",sep="")
        }
        filenames <- c(TADAfile,DAWNfile,hfiles,mfiles)
        tmp <- new_rank(filenames,Tset0[i],genes)
        RANK[kf,,] <- t(tmp)
        kf <- kf + 1
    }
    
    filenames <- c("result/control/TADAdenovo_control.csv","result/control/DAWNcontrol.csv","result/control/hotnetresult1control31.txt","result/control/hotnetresult1control6.txt","result/control/hotnetresult1control7.txt","result/control/hotnetresult1control8.txt","result/control/STRINGCRFresult_controlLBP_31.txt","result/control/iRefCRFresult_controlLBP_6.txt","result/control/coexpCRFresult_controlLBP_7.txt","result/control/InferCRFresult_controlLBP_8.txt")
    conRank <- matrix(0,length(Tset),length(filenames))
    for (i in 1:length(filenames)){
        if(i==1){
            tmp <- read.csv(filenames[i])
            conRank[,i] <- tmp[match(Tset,tmp[,1]),"qvalue.dn"]
        }else if(i==2){
            tmp <- read.csv(filenames[i])
            conRank[,i] <- tmp[match(Tset,tmp[,1]),"FDR"]     
            conRank[is.na(conRank[,i]),i] <- 1
        }else{
            tmp <- read.table(filenames[i])
            conRank[,i] <- tmp[match(Tset,tmp[,1]),5]
        }
    }
    
    #cutoff <- c(0.1,0.1,0.3968514,0.3616888,0.4189705,0.4202216,0.6476619,0.6398611,0.6978669,0.6823814)
    cutoff <- rep(0.1,10)
    load("RANK1")
    cases10 <- 1:dim(RANK)[3]
    contr10 <- 1:dim(RANK)[3]
    for(i in 1:dim(RANK)[3]){
        if(i <=2){
            cases10[i] <- sum(RANK[,2,i] <= cutoff[i])
            contr10[i] <- sum(conRank[,i] <= cutoff[i])
        }else{
            cases10[i] <- sum(RANK[,2,i] <= cutoff[i])
            contr10[i] <- sum(conRank[,i] <= cutoff[i])        
        }
    }
    
    a <- cbind(cases10,contr10)
    pdf(file="plot/FDRs1.pdf",width=11,height=10)
    par(mai=c(4,2,2,2))
    legend=c("TADA","DAWN","RWR-STRING","RWR-iRef","RWR-coexp","RWR-Infer","MIS-STRING","MIS-iRef","MIS-coexp","MIS-Infer")  
    barplot(t(a), names.arg=legend, beside=T, ylab="The number of genes (FDR<=0.1)", 
            cex.names=1.5, las=2, ylim=c(0,max(a)), col=c("darkblue","red"),cex.lab=1.5,cex.axis=1.25)
    box(bty="l")
    legend("topleft",col=c("darkblue","red"),legend=c("ASD case","SSC control"),cex=1.3,fill=c("darkblue","red")) 
    dev.off()    
    
  
    

}

randset_newnet_enri <- function(){
    
    source("enrichana.R")
    dirstr <- "result/randresult3/"
    DAWNnames <-  readLines(con <- file("result/randresult3/TADArandset3.txt","r"))   
    close(con)
    netstr <- c("STRING/","iRef/","coexp/","Infer/","Overlapexp/","Jointexp/","Overlapexp_PPI/","Jointexp_PPI/","Overlap/")
    netnum <- c(31,6,7,8,15:19)
    nnet=9;TPcut=0.1
    
    AUC <- array(0,dim=c(8,10,nnet+2)) 
    DAWNflag <- 0
    for(j in 2:9){
        for(i in 1:10){    
            tmp <- gene_select(j=j,i=i,TPcut=0.1,4,datdir="randresult3/",tadadir="randset3/",Dawnfile="TADArandset3.txt")
            genes <- tmp$genes
            allP <- tmp$allP
            allN <- tmp$allN
            
            TADAfile <- paste("ASD/TADAresult/randset3/TADAdenovo_part",j,"_",i,".csv",sep="")
            k <- which(DAWNnames==paste("TADAdenovo_part",j,"_",i,".csv",sep=""))
            DAWNfile=paste(dirstr,"coexp/DAWN_randset",k,".csv",sep="")

            mfiles <- 1:nnet
            for(k in 1:nnet){
                mfiles[k] <- paste(dirstr,netstr[k],"CRFresult_part",j,"_",i,"LBP_",netnum[k],".txt",sep="")
            }
            filenames <- c(TADAfile,DAWNfile,mfiles)
            if(DAWNflag==1){
                tmp <- read.csv(DAWNfile)
                genes <- intersect(genes,tmp[!is.na(tmp[,"FDR"]),"Gene"])
            }
            print(length(genes))
            tmp <- new_auc(filenames,allP,genes,flag=2,iplot=FALSE,DAWNflag=DAWNflag)
            AUC[j-1,i,] <- tmp$auc
        }
    }
    
    ## plot auc values
    ## Figure 1: AUC values for eight random sets

    cols <- c("Black","Blue","Cyan","Brown","Coral","DarkViolet","Green","HotPink","Orange","Plum",11)
    main="AUC for different number of trios"
    xlab="The number of trios";ylab="AUC";
    ymint <- 1;ymaxt <- 0;
    ymint <- min(AUC);ymaxt <- max(AUC);
    ymin=max(min(ymint),0);ymax = min(ymaxt,1);
    x <- seq(0.2,0.9,0.1)
    nmeth <- length(filenames)
    Y1 <- rep(0,nmeth)
    for(i in 1:nmeth){
        y <- rowMeans(AUC[,,i])
        sd <- apply( AUC[,,i], 1, sd)/2 
        plot_err(x,y,sd,i,main,xlab,ylab,ylim=c(ymin,ymax),xlim=c(0.05,0.95))
        Y1[i] <- y[1]
    }
    text(rep(0.2,11),Y1,labels=c("TADA","DAWN","MIS-STRING","MIS-iRef","MIS-coexp","MIS-Infer","overlapexp","jointexp","overlapexp_PPI","jointexp_PPI","Overlap_ppi"),pos=2,col=cols,cex=1.3)
 
    
    
    source("enrichana.R")
    TADAFile <- "ASD/TADAresult/TADAdenovo_ASDall.csv"
    TADAall <- read.csv(TADAFile)
    Tset <- TADAall[TADAall[,"qvalue.dn"]<= TPcut,1]
    allP <- Tset
    genes <- TADAall[,1]
    
    nnet=4
    #netstr <- c("STRING/","iRef/","coexp/","Infer/","Overlapexp/","Jointexp/","Overlapexp_PPI/","Jointexp_PPI/","Overlap/")
    #netnum <- c(31,6,7,8,15:19)
    netstr <- c("coexp/","Jointexp/","Overlapexp_PPI/","Jointexp_PPI/")
    netnum <- c(7,16,17,18)
    AUC1 <- array(0,dim=c(8,10,nnet)) 
    for(j in 2:9){
        for(i in 1:10){    
            filenames <- 1:nnet
            for(kk in 1:nnet){
                filenames[kk] <- paste(dirstr,netstr[kk],"CRFresult_part",j,"_",i,"LBP_",netnum[kk],".txt",sep="")
                tmp1 <- read.table(filenames[kk]);
                genes <- intersect(tmp1[,1],genes);
                #print(length(genes))
            }
            
            allP <- intersect(genes,allP)
            predictions <- list()
            labels <- list()
            for(kk in 1:nnet){
                filename <- filenames[kk]
                result <- read.table(filename);
                
                result <- result[!duplicated(result[,1]),]
                #result <- result[result[,1] %in% genes,]
               
                predictions[[kk]] <- result[,2]/max(result[,2])
                labels[[kk]] <- rep(-1,length(predictions[[kk]]))
                labels[[kk]][match(allP,result[,1])] <- 1
            }
            pred <- prediction(predictions, labels)
            auc <- unlist(performance(pred, "auc")@y.values)
            
            AUC1[j-1,i,] <- auc
        }
        print(j)
    }
    
    for(i in 1:nnet){
        y <- rowMeans(AUC1[,,i])
       if(i==1){ plot(y,col=1,type="b");
       }else{ lines(y,col=i,type="b");} 
    }
   
}

gene_select <- function(j=2,i=1,TPcut=0.1,nnet=2,datdir="randresult3/",tadadir="randset3/",Dawnfile="TADArandset3.txt",Magifile="MAGIfilenames.txt",filens=""){

    source("enrichana.R")
    TADAFile <- "ASD/TADAresult/TADAdenovo_ASDall.csv"
    TADAall <- read.csv(TADAFile)
    Tset <- TADAall[TADAall[,"qvalue.dn"]<= TPcut,1]

    DAWNnames <-  readLines(con <- file(paste("result/",datdir,Dawnfile,sep=""),"r"))   
    close(con)
    Maginames <-  readLines(con <- file(paste("result/",datdir,"MAGI/",Magifile,sep=""),"r"))   
    close(con)
    
    netstr <- c("STRING/","iRef/","coexp/","Infer/")
    netnum <- c(31,6,7,8)
    
    dirstr <- paste("result/",datdir,sep="")
    ## all positive samples could be predicted by all methods
    #j=2;i=1;
    TADAfile <- paste("ASD/TADAresult/",tadadir,"TADAdenovo_part",j,"_",i,".csv",sep="")
    k <- which(DAWNnames==paste("TADAdenovo_part",j,"_",i,".csv",sep=""))
    DAWNfile=paste("result/",datdir,"coexp/DAWN_randset",k,".csv",sep="")
    
    k <- which(Maginames==paste("ASDrand",j,"depart",i,".csv",sep=""))
    MAGIfile=paste("result/",datdir,"MAGI/RandomGeneList.",k,sep="")
    
    hfiles <- 1:nnet
    mfiles <- 1:nnet
    for(k in 1:nnet){
        hfiles[k] <- paste(dirstr,netstr[k],"hotnetresult1part",j,"_",i,netnum[k],".txt",sep="")
        mfiles[k] <- paste(dirstr,netstr[k],"CRFresult_part",j,"_",i,"LBP_",netnum[k],".txt",sep="")
    }
    filenames <- c(TADAfile,DAWNfile,MAGIfile,hfiles,mfiles)
    
    if(i==0 | j==0){
        filenames <- filens
    }
    
    allP <- Tset
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

new_recurrent <- function(filenames,genes,allgene,flag=2){
    
    #  flag==1: only rank information is used
    #  flag==2: FDR information is used  
    n <- dim(filenames)[1]
    Mlist <- array(0,dim=c(length(genes),2,n))
    for(i in 1:n){
        for(j in 1:2){
            filename <- filenames[i,j]
            if(grepl(".csv",filename)){ result <- read.csv(filename); 
            }else if(grepl(".txt",filename)){result <- read.table(filename);}
            if(i==2){
                DAWNr <- read.csv(filename)
                subs <- is.na(DAWNr[,"FDR"])
                DAWNr1 <- DAWNr[!subs,]
                DAWNr2 <- DAWNr[subs,]
                DAWNr2 <- DAWNr2[order(DAWNr2[,"Stage1_posterior"]),]
                DAWNr <- rbind(DAWNr1,DAWNr2)
                DAWNr <- DAWNr[!duplicated(DAWNr[,1]),]
                DAWNr[is.na(DAWNr[,"FDR"]),"FDR"] <- 1
                result <- DAWNr
            }
            
            result <- result[result[,1] %in% allgene,]
            
            subs <-  which(result[,1] %in% genes)
            if(flag==1){
                Mlist[,j,i] <- subs
            }else if(flag==2){
                if(i==1){Mlist[,j,i] <- result[subs,"qvalue.dn"];}else if(i==2){Mlist[,j,i] <- result[subs,"FDR"];
                }else{Mlist[,j,i] <- result[subs,5];}
            }
            
        }
    }
    
    Mlist
    
}

new_vioplot <- function(partL,restL,main="",xlab="",ylab=""){
    ## plot rank and FDR distributions
    library(caroline)
    legend=c("TADA","DAWN","RWR-STRING","RWR-iRef","MIS-STRING","MIS-iRef")  
    legend <- rep(legend,3)
    L <- c(1,4,7)
    tmp1 <- list()
    for(i in 1:18){
        k <- i%%6
        if(k==0){k=6;}
        subs <- (i-1)%/% 6 + 1
        tmp1[[i]] <- partL[[L[subs]]][,k]
        #print(k)
        #print(L[subs])
    }
    violins(tmp1,names=legend,col=rep(c('purple','lightblue','lightgreen','red','orange','yellow'),3), main=main,xlab=xlab,ylab=ylab,deciles=FALSE,connect="",at=c(1:6,8:13,15:20))
    abline(v=7,col=1,lwd=3)
    abline(v=14,col=1,lwd=3)
    mtext(text=c("20% trios","50% trios","80% trios"), side=3, font=2, line=1,at=c(3.5,10.5,17.5))
    
    tmp2 <- list()
    for(i in 1:18){
        k <- i%%6
        if(k==0){k=6;}
        subs <- (i-1)%/% 6 + 1
        tmp2[[i]] <- restL[[L[subs]]][,k]
        #print(k)
        #print(L[subs])
    }
    violins(tmp2,names=legend,col=rep(c('purple','lightblue','lightgreen','red','orange','yellow'),3), main=main,ylab=ylab,deciles=FALSE,connect="",at=c(1:6,8:13,15:20))
    abline(v=7,col=1,lwd=3)
    abline(v=14,col=1,lwd=3)    
    mtext(text=c("80% trios","50% trios","20% trios"), side=3, font=2, line=1, at=c(3.5,10.5,17.5))
}

new_rank <- function(filenames,gene,genes){
    #  flag==1: only rank information is used
    #  flag==2: FDR information is used
    tmp <- matrix(0,length(filenames),3)
    
    DAWNr <- read.csv(filenames[2],stringsAsFactors=FALSE)
    subs <- is.na(DAWNr[,"FDR"])
    DAWNr1 <- DAWNr[!subs,]
    DAWNr2 <- DAWNr[subs,]
    DAWNr2 <- DAWNr2[order(DAWNr2[,"Stage1_posterior"]),]
    DAWNr <- rbind(DAWNr1,DAWNr2)
    DAWNr <- DAWNr[!duplicated(DAWNr[,1]),]
    DAWNr <- DAWNr[DAWNr[,"Gene"] %in% genes,]
    DAWNr[is.na(DAWNr[,"FDR"]),"FDR"] <- 1
    tmp[2,1] <- which(DAWNr[,1]==gene) ## rank
    tmp[2,2] <- DAWNr[tmp[2,1],"FDR"] ## FDR
    tmp[2,3] <- DAWNr[tmp[2,1],"Stage2_posterior"] ## Stage2_posterior
    
    for(k in c(1,3:length(filenames))){
        filename <- filenames[k]
        if(grepl(".csv",filename)){ result <- read.csv(filename,stringsAsFactors=FALSE); 
        }else if(grepl(".txt",filename)){result <- read.table(filename,stringsAsFactors=FALSE);}
        
        result <- result[!duplicated(result[,1]),]
        result <- result[result[,1] %in% genes,]
        
        tmp[k,1] <- which(result[,1]==gene) ## rank
        if(grepl(".csv",filename)){
            tmp[k,2] <- result[tmp[k,1],"qvalue.dn"]
            tmp[k,3] <- result[tmp[k,1],"BF.dn"]
        }else{
            tmp[k,2] <- result[tmp[k,1],5]
            tmp[k,3] <- result[tmp[k,1],2]
        }
        
    } 
    
    tmp
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
        cols <- c("Black","Blue","Cyan","Brown","Coral","DarkViolet","Green","HotPink","Orange","Plum","deeppink")
        plot(perf,main="Comparison of Performance",xlab="False positive rate",ylab="True positive rate",cex.lab=1.2,cex.main=1.2,xlim=c(0,1),ylim=c(0,1))
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
        legend=c("TADA","DAWN","MAGI","RWR-STRING","RWR-iRef","RWR-coexp","RWR-Infer","MIS-STRING","MIS-iRef","MIS-coexp","MIS-Infer")
        legend <- paste(legend,format(auc,digits=3),sep="--AUC:")
        legend("bottomright",col=cols[1:length(filenames)],legend=legend,lty=rep(1,length(filenames)),lwd=rep(2,length(filenames)),cex=1.2,y.intersp=0.5) 

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

myplot <- function(predictions,labels){

    n <- length(predictions);
    cutV <- seq(0,1,0.001)
    kc <- length(cutV)
    for(i in 1:n){
        x <- sapply(1:kc, function(k) sum(labels[[i]][1 - predictions[[i]] <= cutV[k]] == -1))
        x <- x/sum(labels[[i]]==-1)
        y <- sapply(1:kc, function(k) sum(labels[[i]][1 - predictions[[i]] <= cutV[k]] ==1) )
        y <- y/sum(labels[[i]]==1)
                    
        if(i==1){
            plot(x,y,col=i,type="l")
        }else{
            lines(x,y,col=i,type="l")
        }
    }

}

plot_two <- function(x,main,xlab,ylab,RankM,pos="topright",nmeth=4){
    
    ymint <- 1;ymaxt <- 0;
    for(i in 1:nmeth){ymint <- min(min(RankM[,i]),ymint);ymaxt <- max(max(RankM[,i]),ymaxt);}
    ymin=max(min(ymint)-1,0);ymax=max(ymaxt) + 1;
    for(i in 1:nmeth){
        if(i==1){
            plot(x,RankM[,1],main=main,xlab=xlab,ylab=ylab,col=1,type="b",ylim=c(ymin,ymax))
        }else{
            lines(x,RankM[,i],col=i,type="b")
        }
    }
    legend(pos,col=(1:nmeth),legend=c("TADA","DAWN","RWR","MIS"),lty=rep(1,nmeth),lwd=rep(1,nmeth))
    
}

plot_one <- function(x,main,xlab,ylab,result,fig=1,nmeth=4){
    
    ymint <- 1;ymaxt <- 0;
    for(i in 1:nmeth){ymint <- min(min(result[[i]]),ymint);ymaxt <- max(max(result[[i]]),ymaxt);}
    ymin=max(min(ymint)-0.1,0);ymax=max(ymaxt)+0.1
    for(i in 1:nmeth){
        y <- rowMeans(result[[i]])
        sd <- apply( result[[i]], 1, sd) 
        plot_err(x,y,sd,i,main,xlab,ylab,ylim=c(ymin,ymax))
    }
    #if(fig==1){ lines(c(0,1),c(0,1),col=1,lty=2,type="l");}
    legend("bottomright",col=(1:nmeth),legend=c("TADA","DAWN","RWR","MIS"),lty=rep(1,nmeth),lwd=rep(1,nmeth))

}

one_auc <- function(filename,Tset,Fixspe=0.8,FixFDR=0.3,FixNUM=500,pl=FALSE){
    
    if(grepl(".csv",filename) & grepl("DAWN",filename)){ result <- read.csv(filename); predictions <- result[,"FDR"]; predictions[is.na(predictions)] <- 1;
    }else if(grepl(".csv",filename) & grepl("TADA",filename)){ result <- read.csv(filename); predictions <- result[,"qvalue.dn"];
    }else if(grepl(".txt",filename)){result <- read.table(filename); predictions <- result[,5]}
    
    geneset <- result[predictions<=FixFDR,1]
    
    predictions <- 1 - predictions 
    labels <- rep(-1,length(predictions))
    labels[match(Tset,result[,1])] <- 1
    
    pred <- prediction(predictions, labels)
    auc <- unlist(performance(pred, "auc")@y.values)
    perf <- performance(pred, measure = "tpr", x.measure = "fpr") 
    if(pl){
        plot(perf,main="ROC plot",xlab="1-specificity",ylab="sensitivity")
    }
    
    fixX <- 1-Fixspe
    senspe <- perf@y.values[[1]][which.min(abs(perf@x.values[[1]]-fixX))]
 
    fixF <- 1 - FixFDR  
    predictions1 <- predictions[predictions>=fixF]
    labels1 <- labels[predictions>=fixF]
    if(length(unique(labels1))==2){
        pred1 <- prediction(predictions1, labels1)
        aucfdr <- unlist(performance(pred1, "auc")@y.values)   
    }else{ aucfdr <- 1;}#!!!!!
    
    predictions2 <- predictions[1:FixNUM]
    labels2 <- labels[1:FixNUM]
    pred2 <- prediction(predictions2, labels2)
    aucnum <- unlist(performance(pred2, "auc")@y.values) 
    
    predictions3 <- 1 - ((1:FixNUM)/FixNUM)
    labels3 <- labels[1:FixNUM]
    pred3 <- prediction(predictions3, labels3)
    aucrank <- unlist(performance(pred3, "auc")@y.values)    
    
    list(auc=auc,aucfdr=aucfdr,aucnum=aucnum,aucrank=aucrank,sensitivity=senspe,geneset=geneset)
}

plot_err <- function(x,y,sd,flag=1,main="ROC",xlab="Samples",ylab="AUC",xlim=c(0,1),ylim=c(0,1)){
    cols <- c("Black","Blue","Cyan","Brown","Coral","DarkViolet","Green","HotPink","Orange","Plum")
    if(flag==1){
        par(mar=c(6,6,3,2))
        plot(x, y,xaxt="n", type="b",col=cols[flag], main=main, ylab=ylab,xlab=xlab,xlim=xlim,ylim=ylim,cex.lab=2,cex.main=2,lwd=2,cex.axis=1.8)
        axis(1,at=seq(0.2,0.9,0.1),labels=floor(3953*seq(0.2,0.9,0.1)),font=1,cex.axis=1.5)
    }else{
        lines(x, y, type="b",col=cols[flag],xaxt="n")
    }
    
    #axis(2,at=seq(ylim[1],ylim[2],0.1))
    #segments(x,y-sd,x,y+sd)
    #epsilon <- 0.004
    #segments(x-epsilon,y-sd,x+epsilon,y-sd)
    #segments(x-epsilon,y+sd,x+epsilon,y+sd)

}

bar_err <- function(G,main,xlab,ylab,groups=c("TADA","DAWN","RWR","MIS")){

    nmeth <- dim(G)[2]
    means <- rep(0,nmeth)
    stDevs <- rep(0,nmeth)
    for(i in 1:nmeth){
        means[i] <- mean(G[,i])
        stDevs[i] <- sd(G[,i])
    }
    ymin <- min(G); ymax <- max(G)
    # Get the mean values of both groups
    # Create the plot (without axes, we will add those later)
    # Range of y-axis is set to interval (0, 60) to give space 
    # to the error bars
    # Midpoints of the bars will be saved in variable mp
    mp <- barplot(means, axes=FALSE, axisnames=FALSE, ylim=c(ymin, ymax),
                  col=1:nmeth, main=main, xlab=xlab, ylab=ylab)
    
    
    # The x-axis with labels for each group
    # The labels are drawn at the bar midpoints
    axis(1, labels=groups, at = mp)
    # The y-axis with the age going from 0 to 60 
    axis(2, at=seq(ymin , ymax,length.out=10))
    # Put the plot in a box
    box()
    
    # 4. Add the error bars

    # Get standard deviation of each group
    # The standard deviations are saved in a matrix of same size 
    # as the matrix with midpoints, this is useful for plotting 
    # the error bars
    
    # Plot the vertical lines of the error bars
    # The vertical bars are plotted at the midpoints
    segments(mp, means - stDevs, mp, means + stDevs, lwd=2)
    # Now plot the horizontal bounds for the error bars
    # 1. The lower bar
    segments(mp - 0.1, means - stDevs, mp + 0.1, means - stDevs, lwd=2)
    # 2. The upper bar
    segments(mp - 0.1, means + stDevs, mp + 0.1, means + stDevs, lwd=2)

}

landscape<- function(){
    
    source("enrichana.R")
    filenames <- c("ASD/TADAresult/TADAdenovo_ASDall.csv","result/randresult3/coexp/DAWN_randset1.csv","result/randresult3/STRING/hotnetresult1ASDall31.txt","result/randresult3/iRef/hotnetresult1ASDall6.txt","result/randresult3/coexp/hotnetresult1ASDall7.txt","result/randresult3/Infer/hotnetresult1ASDall8.txt","result/randresult3/STRING/CRFresult_ASDallLBP_31.txt","result/randresult3/iRef/CRFresult_ASDallLBP_6.txt","result/randresult3/coexp/CRFresult_ASDallLBP_7.txt","result/randresult3/Infer/CRFresult_ASDallLBP_8.txt")
    
    TPcut <- 0.1
    TADAFile <- "ASD/TADAresult/TADAdenovo_ASDall.csv"
    TADAall <- read.csv(TADAFile)
    Tset <- TADAall[TADAall[,"qvalue.dn"]<= TPcut,1]
    
    # enrichment of FMRP and RBFOX targets; SFARI autism genes with syndromic S, High confidence 1, Strong Candidate 2,
    tmp <- test_genes()
    testgenes <- tmp$testgenes
    Lsets <- tmp$Lsets

    tmp <- total_genes(filenames,Tset,flag=2,alpha=0.1)
    Tgenes <- tmp$topgenes
    result <- tmp$result
 
    testresult <- enrich_genes(Tgenes,result,testgenes,Lsets)
    
    cols <- c("Black","Blue","Cyan","Brown","Coral","DarkViolet","Green","HotPink","Orange","Plum")
    legend=c("TADA","DAWN","RWR-STRING","RWR-iRef","RWR-coexp","RWR-Infer","MIS-STRING","MIS-iRef","MIS-coexp","MIS-Infer") 
    #att <- cytoscape_attri(Tgenes,testresult,Lsets,legend,netflag=7,strn="Topgene")
    
    plot.new()
    par(bg=NA,mai=c(5,2,2,2))
    Foldc <- testresult$Foldc
    n <- dim(Foldc)[2]
    for(i in 1:dim(Foldc)[1]){
        if(i==1){
            plot(1:n,Foldc[i,],xaxt="n",col=cols[i],type="b",main="Fold changes",xlab="",ylab="Fold changes")
        }else{
            lines(1:n,Foldc[i,],col=cols[i],type="b")
        }
    }
    nmeth=length(filenames)
    legend("topright",col=cols,legend=legend,lty=rep(1,nmeth),lwd=rep(1,nmeth))
    axis(1, at=1:n, labels=Lsets,las=2)
}

cytoscape_attri <- function(Tgenes,testresult,Lsets,legend,netflag=7,strn="Topgene"){
    
    n <- length(testresult$Lgenes)
    Lgenes <- testresult$Lgenes
    Ugenes <- c()
    for(i in 1:n){
        Ugenes <- union(Tgenes[[i]],Ugenes)
    }
    
    att <- matrix(0.0,length(Ugenes),length(legend)+length(Lsets)+2,dimnames=list(Ugenes,c(Lsets,legend,"ALL","GeneSets")))
    for(i in 1:n){
        tmpg <- rownames(Lgenes[[i]])
        for(j in 1:length(Lgenes[[i]])){
            tmpn <- setdiff(unlist(strsplit(Lgenes[[i]][j],"\\|")),"")
            tmpn <- union(tmpn,legend[i])
            #att[tmpg[j],tmpn] <- tmpn
            att[tmpg[j],tmpn] <- 1.1
        }
    }
    att[,"ALL"] <- rowSums(att[,legend])
    att[,"GeneSets"] <- rowSums(att[,Lsets])
    write.csv(att,file=paste(strn,"attri.csv",sep=""),quote=FALSE)
   
    ### wirte subnet to show 
    if(netflag==7) {filename <- "data/hotnet/iRefIndexm.txt";}
    net.text <- as.matrix(read.table(filename,sep="\t",header=FALSE))
    subs <- net.text[,1] %in% Ugenes & net.text[,2] %in% Ugenes
    subnet <- net.text[subs,]
    isonodes <- setdiff(Ugenes,union(subnet[,1],subnet[,2]))
    isonodes <- cbind(isonodes,"gap",1)
    subnet <- rbind(subnet,isonodes)
    write.table(subnet,file=paste(strn,"net.txt",sep=""),row.names=FALSE,quote=FALSE,col.names=FALSE,sep="\t")

    att
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

total_genes <- function(filenames,Tset,flag=1,alpha=0.1){
    
    library(ROCR)
    predictions <- list()
    labels <- list()
    result <- list()
    topgenes <- list()
    n <- length(filenames)
    
    if(flag==1){
    DAWNr <- read.csv(filenames[2])
    subs <- is.na(DAWNr[,"FDR"])
    result[[2]] <- DAWNr[!subs,]
    i=2
    predictions[[i]] <- 1 - result[[i]][,"Stage2_posterior"] #result[[i]][,"FDR"] #
    labels[[i]] <- rep(-1,length(predictions[[i]]))
    labels[[i]][match(Tset,result[[i]][,1])] <- 1
    
    for(i in c(1,3:length(filenames))){
        filename <- filenames[i]
        if(grepl(".csv",filename)){ tmp <- read.csv(filename); 
        }else if(grepl(".txt",filename)){tmp <- read.table(filename);}
        result[[i]] <- tmp[!duplicated(tmp[,1]),]
        if(grepl(".csv",filename)){
            predictions[[i]] <- 1 - result[[i]][,"qvalue.dn"]
        }else{
            predictions[[i]] <- result[[i]][,2] #1 - result[[i]][,5]
        }
        
        labels[[i]] <- rep(-1,length(predictions[[i]]))
        labels[[i]][match(Tset,result[[i]][,1])] <- 1
    }

    pred <- prediction(predictions, labels)
    Fmea <- performance(pred, "f",alpha=alpha)@y.values
    Cutoffs <- performance(pred, "f")@x.values
    cutoff <- 1:n
    for(i in 1:n){
        tmp <- Fmea[[i]]
        tmp[is.na(tmp)] <- 0
        subs <- which.max(tmp)
        cutoff[i] <- Cutoffs[[i]][subs]
        topgenes[[i]] <- result[[i]][predictions[[i]] >= cutoff[i],1]
    }
    cutoff <- 1 - cutoff
    }else if(flag==2){
        TPcut <- alpha
        for(k in 1:length(filenames)){
            if(k==2){
                DAWNr <- read.csv(filenames[2])
                subs <- is.na(DAWNr[,"FDR"])
                result[[2]] <- DAWNr[!subs,]
                genes <- result[[2]][result[[2]][,"FDR"] <= TPcut,1]
            }else{
                filename <- filenames[k]
                if(grepl(".csv",filename)){ tmp <- read.csv(filename); 
                }else if(grepl(".txt",filename)){tmp <- read.table(filename);}
                result[[k]] <- tmp[!duplicated(tmp[,1]),]
                if(grepl(".csv",filename)){
                    genes <- result[[k]][result[[k]][,"qvalue.dn"] <= TPcut,1];
                }else{ genes <- result[[k]][result[[k]][,5]<= TPcut,1]; }
            }
            topgenes[[k]] <- genes
        }
        
        cutoff=alpha;
    }
    
    
    list(result=result,topgenes=topgenes,cutoff=cutoff)
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
            }else if(grepl(".txt",filename)){tmp <- read.table(filename);}

            if(grepl("DAWN",filename)){
                subs <- is.na(tmp[,"FDR"])
                result[[i]] <- tmp[!subs,]
            }else{
                #result[[i]] <- tmp[!duplicated(tmp[,1]),]
                result[[i]] <- tmp       
            }
            
            tmpn <- length(unlist(strsplit(filename,"/")))
            subfilename <- unlist(strsplit(filename,"/"))[tmpn]
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
            Cutoffs <- performance(pred, "f")@x.values
            Fmea[[1]][is.na(Fmea[[1]])] <- 0
            subs <- which.max(Fmea[[1]])
            cutoff[i] <- Cutoffs[[1]][subs]
        }
   
    cutoff
}

enrich_genes <- function(Tgenes,result,testgenes,Lsets){
    
    n <- length(result)
    n.t  <- length(testgenes)
    
    Overlap <- matrix(0,n,n.t)
    Foldc <- matrix(0,n,n.t)
    fisherp <- matrix(0,n,n.t)
    Lgenes <- list()
    
    N <- sapply(1:n,function(i) dim(result[[i]])[1])
    N.t <- sapply(1:n.t,function(i) length(testgenes[[i]]))
    NT <- sapply(1:n,function(i) length(Tgenes[[i]]))
    
    for(i in 1:n){
        Lgenes[[i]] <- matrix(rep("",NT[i]),nrow=NT[i],ncol=1)
        rownames(Lgenes[[i]]) <- Tgenes[[i]]
        for(j in 1:n.t){
            lgenes <- intersect(Tgenes[[i]],testgenes[[j]])
            Lgenes[[i]][lgenes,1] <- paste(Lgenes[[i]][lgenes,1],Lsets[j],sep="|") 
            Overlap[i,j] <- length(lgenes)
            Foldc[i,j] <- (Overlap[i,j]/NT[i]) / (N.t[j]/N[i])
            fisherp[i,j] <- fisher.test( matrix(c(Overlap[i,j],NT[i],N.t[j],N[i]),2,2),alternative="greater")$p.value
        }
    }
    
    list(Lgenes=Lgenes,Overlap=Overlap,Foldc=Foldc,fisherp=fisherp)
}

overlap_genes <- function(Tgenes,result,netflag=7,legend=c("TADA","DAWN","RWR-STRING","RWR-iRef","MIS-STRING","MIS-iRef")){
    
    n <- length(Tgenes)
    Ugenes <- c()
    for(i in 1:n){
        Ugenes <- union(Tgenes[[i]],Ugenes)
    }
    
    Lgenes <- matrix("",length(Ugenes),1)
    rownames(Lgenes) <- Ugenes
    for(i in 1:n){
        tmp <- intersect(Ugenes,Tgenes[[i]])
        Lgenes[tmp,1] <- paste(Lgenes[tmp,1],legend[i],sep="|") 
    }
    
    write.csv(cbind(Ugenes,Lgenes),file="Topgenes.csv",row.names=FALSE)
    write.table(paste(Ugenes," = ", Lgenes, sep=""),file="Topgenes.txt",row.names=FALSE,quote=FALSE,col.names=FALSE)
    
    if(netflag==7) {filename <- "hotnet/iRefIndexm.txt";}
    
    net.text <- as.matrix(read.table(filename,sep="\t",header=FALSE))
    subs <- net.text[,1] %in% Ugenes & net.text[,2] %in% Ugenes
    subnet <- net.text[subs,]
    write.table(subnet,file="Topgenenet.txt",row.names=FALSE,quote=FALSE,col.names=FALSE,sep="\t")

}

case_gene <- function(){
    source("CRF_build.R")
    source("Multi_net.R")
    allnet <- read.table("hotnet/iRefIndex.txt")
    result <- list()
    result[[1]] <- read.csv("ASD/TADAresult/TADAdenovo_ASDall.csv")
    result[[2]] <- read.csv("result/randresult/coexp/DAWN_randset1.csv")
    result[[4]] <- read.table("result/randresult/iRef/hotnetresult1ASDall6.txt")
    result[[6]] <- read.table("result/randresult/iRef/CRFresult_ASDallLBP_6.txt")
    result[[3]] <- read.table("result/randresult/STRING/hotnetresult1ASDall31.txt")
    result[[5]] <- read.table("result/randresult/STRING/CRFresult_ASDallLBP_31.txt")
    
    genes=c("EP300","KDM5B","CUL3")
    P <- list()
    for(k in 1:length(genes)){
        gene <- genes[k]
        subs <- allnet[,1] %in% gene | allnet[,2] %in% gene
        neighbors <- union(allnet[subs,1],allnet[subs,2])
        #neighbors <- union(neighbors,gene)
        neighbors <- setdiff(neighbors,gene)
    
        ranks <- 1:6
        for(i in 1:6){
         ranks[i] <- which(result[[i]][,1]==gene)
        }
        pi=0.06
        BF <- result[[1]][result[[1]][,1] %in% neighbors,"BF.dn"]
        p <- pi*BF/(1-pi+pi*BF)
        P[[k]] <- p
        print(ranks)
        print(length(neighbors))
        print(sum(p>=0.5))
    }
    
    vioplot(P[[1]],P[[2]],P[[3]],names=genes)
    
    
    ### ALL top genes
    TPcut=0.1
    TADAFile <- "ASD/TADAresult/TADAdenovo_ASDall.csv"
    TADAall <- read.csv(TADAFile)
    Tset0 <- TADAall[TADAall[,"qvalue.dn"]<= TPcut,1]
    
    Tset <- Tset0
    i=1
    DAWNr <- read.csv(paste("result/leaveoneresult/DAWN_randset2_",i,".csv",sep=""))
    hotnetr <- read.table(paste("result/leaveoneresult/hotnetresult1rand2_",i,"6.txt",sep=""))
    MISr <- read.table(paste("result/leaveoneresult/CRFresult_rand2_",i,"LBP_61.txt",sep=""))
    Tset <- intersect(Tset,DAWNr[,1])
    Tset <- intersect(Tset,hotnetr[,1])
    Tset <- intersect(Tset,MISr[,1])
    
    Tset <- intersect(Tset,allnet[,1])
    Tset <- intersect(Tset,allnet[,2])
    load("geneexp.Rdata")
    module <- read.csv("module.csv")
    library(WGCNA)
    Fp <- 1:length(Tset)
    Fp1 <- Fp
    for(k in 1:length(Tset)){
        gene <- Tset[k]
        subs <- allnet[,1] %in% gene | allnet[,2] %in% gene
        neighbors <- union(allnet[subs,1],allnet[subs,2])
        neighbors <- setdiff(neighbors,gene)
        
        pi=0.06
        BF <- result[[1]][result[[1]][,1] %in% neighbors,"BF.dn"]
        p <- pi*BF/(1-pi+pi*BF)
        b <- length(neighbors)
        a <- sum(p>=0.5)
        
        label <- module[module[,1]==gene,2]
        modgenes <- module[module[,2]==label,1]
        adj = adjacency(t(data[rownames(data) %in% modgenes,]), type = "unsigned", power=1)
        subs <- rownames(adj)==gene
        subs <- adj[subs,] >= 0.7
        nei1 <- rownames(adj)[subs]
        pi=0.06
        BF <- result[[1]][result[[1]][,1] %in% nei1,"BF.dn"]
        p <- pi*BF/(1-pi+pi*BF)
        d <- length(nei1)
        c <- sum(p>=0.5)        
        
        Fp[k] <- (a/b) #fisher.test(matrix(c(a,b,c,d),2,2),alternative="less")$p.value
        Fp1[k] <- (c/d)
    }
    
    plot(Fp1,type="b",col=1,ylim=c(0,max(Fp1)),ylab="Fractions of risk genes")
    lines(Fp,type="b",col=2)
    legend("topleft",col=(1:2),legend=c("co-expression","PPI"),lty=rep(1,2),lwd=rep(1,2))
    
    
    TADAfiles <- c("ASD/TADAresult/TADAdenovo_ASDall.csv","ASD/TADAresult/TADAdenovo_ASDall.csv","ASD/TADAresult/TADAdenovo_ASDall.csv","ASD/TADAresult/TADAdenovo_ASDall.csv")
    NDDlist <- list()    
    NDDlist[[1]] <- c("ASD/TADAresult/TADAdenovo_ASDall.csv","result/randresult/DAWN_randset1.csv","result/randresult/hotnetresult1ASDall6.txt","result/randresult/CRFresult_ASDallLBP_61.txt")
    NDDlist[[2]] <- c("ASD/TADAresult/TADAdenovo_part9_5.csv","result/randresult/DAWN_randset41.csv","result/randresult/hotnetresult1part9_56.txt","result/randresult/CRFresult_part9_5LBP_61.txt")
    NDDlist[[3]] <- c("ASD/TADAresult/TADAdenovo_part5_5.csv","result/randresult/DAWN_randset26.csv","result/randresult/hotnetresult1part5_56.txt","result/randresult/CRFresult_part5_5LBP_61.txt")
    dataname <- c("ALL samples","90 percent samples","50 percent samples")
    genes=""
    enric5(NDDlist,TADAfiles,dataname,genes,flag=1,N=0.3,qval="qvalue.dn")
}

auc_bug <- function(){

predictions1 <- 1- (1:100/100)
predictions2 <- 1- (seq(30,100,length.out=100)/100)
max(predictions2)

max(predictions1)

pred1 <- prediction(predictions1, labels)
pred2 <- prediction(predictions2, labels)
unlist(performance(pred1, "auc")@y.values)

unlist(performance(pred2, "auc")@y.values)


}

CHD_results <- function(){
#     source("Network_analysis.R")
#     source("enrichana.R")
#     source("CRF_build.R")
#     options(stringsAsFactors=FALSE)
#     filenames<- c("ASD/TADAresult/TADAdenovo_resultsPCGC.csv","result/PCGC/coexpCRFresult_PCGCLBP_121.txt")
#     
#     a <- read.table(filenames[2])
#     b <- read.csv(filenames[1])
#     topgene <- a[a[,2]>=0.6,1]
#     mtable <- cbind(topgene,b[match(topgene,b[,1]),c("dn.LoF","dn.mis3")])
#     write.csv(mtable,file="functions/topgenes.csv",row.names=FALSE)
#     
#     strname <- c("TADA","MIS")
#     genelist <- function_genes()
#     k <- length(genelist)
#     TPcut=0.1
#     
#     for(i in 1:length(filenames)){
#         filename <- filenames[i]
#         if(grepl(".csv",filename)){
#             tmp <- read.csv(filename);genes <- tmp[tmp[,"qvalue.dn"]<=TPcut,1]
#         }else{
#             tmp <- read.table(filename);genes <- tmp[tmp[,5]<=TPcut,1]
#         }
#         genelist[[k+i]] <- genes
#     }
#     Lsets <- c("HMG","TGF","Notch","Wnt","Chromation","TADA","MIS")
#     Ugenes <- c()
#     for(i in 1:length(genelist)){
#         genelist[[i]] <- mapping_to(unique(genelist[[i]]))
#         Ugenes <- union(Ugenes,genelist[[i]])
#         write.table(genelist[[i]],file=paste("functions/",Lsets[i],"gene.txt",sep=""),quote=FALSE,col.names=FALSE,row.names=FALSE)
#         write.table(paste(genelist[[i]]," = ", i,sep=""),file=paste("functions/",Lsets[i],"gene.txt",sep=""),quote=FALSE,col.names=Lsets[i],row.names=FALSE)
#     }
#     n <- length(Ugenes)
#    
#     att <- matrix(0.0,n,length(genelist),dimnames=list(Ugenes,Lsets))
#     for(i in 1:length(Lsets)){
#         att[Ugenes %in% genelist[[i]] ,i] <- 1.1
#     }
#     write.csv(att,file=paste("functions/","PCGCatt.csv",sep=""),quote=FALSE)
#     
#     mapT <- as.matrix(read.delim("Fmap0121.txt",header=FALSE,sep="\t"))
#     filename <- "data/expressiondata/PCGCbatch2net_t.txt";
#     tmpgenes <- union(genelist[[6]],genelist[[7]])
    
    
    source("Network_analysis.R")
    source("enrichana.R")
    source("CRF_build.R")
    options(stringsAsFactors=FALSE)
    filenames<- c("ASD/TADAresult/TADAdenovo_resultsPCGC.csv","result/PCGC/coexpCRFresult_PCGCLBP_121.txt","result/PCGC/CRF_inputPCGC.txt")
    mapT <- as.matrix(read.delim("Fmap0121.txt",header=FALSE,sep="\t"))
    
    a <- read.table(filenames[2])
    b <- read.csv(filenames[1])
    topgene <- a[a[,2]>=0.6,1]
    mtable <- cbind(topgene,b[match(topgene,b[,1]),c("dn.LoF","dn.mis3")])
    genes1 <- b[b[,"qvalue.dn"]<=0.1,1]
    genes2 <- a[a[,2]>=0.6,1]
    inputs <- read.table(filenames[3])
    tmp <- inputs[inputs[,2]>=0.21,1] # quantile 0.99
    tmp <- mapT[match(tmp,mapT[,1]),2]
    genes3 <- union(a[a[,5]<=0.3,1],tmp)
        
    filename <- "data/expressiondata/PCGCbatch2net_t.txt";
    tmpgenes <- union(genes1,genes2)
    genes <- mapT[match(tmpgenes,mapT[,2]),1]
    net <- read.table(filename)
    subs <- net[,1] %in% genes | net[,2] %in% genes
    subnet <- net[subs,]
    subnet[,1] <- mapT[match(subnet[,1],mapT[,1]),2]
    subnet[,2] <- mapT[match(subnet[,2],mapT[,1]),2]
    write.table(subnet,file="functions/subnetPCGC.txt",sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)
 
    net.text <- as.matrix(read.table("functions/subnetPCGC.txt",sep="\t",header=FALSE))
    net.text <- rbind(net.text,net.text[,c(2,1,3)])
    a <- read_net(net.text)
    a$matrix[lower.tri(a$matrix,diag=TRUE)] <- 0
    edges <- which(a$matrix>0,arr.ind=TRUE)
    subnet <- cbind(a$node[edges[,1]],a$node[edges[,2]],1)
    write.table(subnet,file="functions/subnetPCGC1.txt",sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)
    
    a <- read_net(net.text)
    nodes <- cbind(a$node," = ",rowSums(a$matrix>0))
    nodes[!(nodes[,1] %in% tmpgenes),3] <- 1
    write.table(nodes,file="functions/genesize.txt",quote=FALSE,col.names=FALSE,row.names=FALSE) 
    nodes1 <- nodes
    nodes1[,3] <- 0
    nodes1[nodes1[,1] %in% genes1,3] <- 1
    nodes1[nodes1[,1] %in% setdiff(genes2,genes1),3] <- 2
    nodes1[nodes1[,1] %in% setdiff(genes3,genes2),3] <- 3
    write.table(nodes1,file="functions/genecolor.txt",quote=FALSE,col.names=FALSE,row.names=FALSE)
    nodes2 <- nodes
    nodes2[,3] <- ""
    nodes2[ match(tmpgenes,nodes2[,1]),3] <- tmpgenes
    grgene <- intersect(setdiff(genes3,genes2),nodes2[,1])
    nodes2[ match(grgene,nodes2[,1]),3] <- grgene
    write.table(nodes2,file="functions/genenames.txt",quote=FALSE,col.names=FALSE,row.names=FALSE)
    
#     filename <- "data/hotnet/iRefIndexm.txt"
#     net <- read.table(filename)
#     subs <- net[,1] %in% Ugenes & net[,2] %in% Ugenes
#     subnet <- net[subs,]
#     write.table(subnet,file="functions/subnetPCGCiRef.txt",sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)
#     
#     filename <- "data/hotnet/iRefIndexm.txt"
#     net <- read.table(filename)
#     subs <- net[,1] %in% union(genelist[[6]],genelist[[7]]) & net[,2] %in% union(genelist[[6]],genelist[[7]])
#     subnet <- net[subs,]
#     write.table(subnet,file="functions/subnetPCGCMIS.txt",sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)
    
}

function_genes <- function(){

    filenames <- c("functions/pathways/HMG_pathway.txt","functions/pathways/TGF_beta_pathway.txt","functions/pathways/Notch_pathway.txt","functions/pathways/Wnt_pathway.txt","functions/pathways/REACTOME/Chromatin_organization.txt")
    genelist <- list()
    genelist[[1]] <- read.delim(filenames[1],sep="\t",skip=2,header=FALSE)[,2] ## HMG genes
    
    for(k in 2:4){
    a <- read.delim(filenames[k],sep="\t",skip=1,header=FALSE)
    tmp <- sapply(seq(2,dim(a)[1],2),function(i) unlist(strsplit(a[i,1],";"))[1])
    genelist[[k]] <- tmp
    }
    
    a <- read.delim(filenames[5],sep=" ",skip=1,header=FALSE)
    genelist[[5]] <- a[,2]

    genelist
}

PCGC_results <- function(){
    options(stringsAsFactors=FALSE)
    source("Multi_net.R")
    source("CRF_build.R")
    source("Network_analysis.R")
    filename <- "result/PCGC/coexpCRFresult_PCGCLBP_121.txt"
    TADAfile <- "ASD/TADAresult/TADAdenovo_resultsPCGC.csv"
    
    countT <- read.csv("ASD/PCGCmutation.csv")
    coexpnet <- build_net(netflag=12,"")
    mapT <- as.matrix(read.delim("Fmap0121.txt",header=FALSE,sep="\t"))
    coexpnet$node <- mapT[match(coexpnet$node,mapT[,1]),2]
    dimnames(coexpnet$matrix) <- list(coexpnet$node,coexpnet$node)

    CRFr <- read.table(filename)
    TADAr <- read.csv(TADAfile)
    TADAr[,1] <- mapping_to(TADAr[,1])
    
    genes <- intersect(CRFr[,1],TADAr[,1])
    CRFr <- CRFr[CRFr[,1] %in% genes,]
    TADAr <- TADAr[TADAr[,1] %in% genes,]
    
    genes1 <- CRFr[CRFr[,5]<=0.1,1]
    genes2 <- TADAr[TADAr[,"qvalue.dn"]<=0.3,1]
    
    subs <- match(CRFr[,1],TADAr[,1])
    resultT <- cbind(CRFr[,1],TADAr[subs,c("dn.LoF","dn.mis3")])
    #resultT[is.na(resultT)] <- 0
    
    source("enrichana.R")
    poipvalue <- poisson_testhq(resultT,N=1228)
    
    resultT <- cbind(resultT,poipvalue,TADAr[subs,"qvalue.dn"])
    resultT <- cbind(resultT,CRFr[,c(2,5)])
    
    adj <- coexpnet$matrix[genes,genes]
    adj[adj>0] <- 1
    resultT <- cbind(resultT,rowSums(adj))
    
    resultT <- cbind(resultT,rowSums(adj[,genes1]))
    resultT <- cbind(resultT,rowSums(adj[,genes2]))
    
    colnames(resultT) <- c("Gene","dn.LoF","dn.mis3","poission_test_LOF_pvalue","poission_test_MIS_pvalue","poission_test_LOFMIS_pvalue","TADA_FDR","MIS_score","MIS_FDR","Degree_in_network","edges_in_FDR<=0.1_MISresult","edges_in_FDR<=0.3_TADAresult")
    
    write.csv(resultT,file="result/PCGC/resultTable.csv",row.names=FALSE)

}

poisson_testhq <- function(data,N,mutfile="ASD/GeneMutaRatem.csv"){
    
    mutrate <- read.csv(mutfile,stringsAsFactors=FALSE)
    n <- dim(data)[1]
    logpM <- sapply(1:n, function(i){
        sub <- match(data[i,1],mutrate[,1])
        if(!is.na(sub)){
        a1 <- poisson.test(x=data[i,2], T = N *2* as.numeric(mutrate[sub,"LOF"]), alternative = "greater", conf.level = 0.95)$p.value
        a2 <- poisson.test(x=data[i,3], T = N *2* as.numeric(mutrate[sub,"dmis"]), alternative = "greater", conf.level = 0.95)$p.value
        a3 <- poisson.test(x=data[i,2]+ data[i,3], T = N *2* as.numeric(mutrate[sub,"LOFDmis"]), alternative = "greater", conf.level = 0.95)$p.value
        }else{
        a1=1;a2=1;a3=1;
        }
        c(a1,a2,a3)
        })
    logpM <- t(logpM)
    
    logpM
}