randset_1_AUC_6_3 <- function(){
    
    source("enrichana.R")
    
    load("BestP_6_10")
    betaV <- BestP[,3]
    filenames <- filenamesf(0,0,betaV=betaV,"","")
    allP <- allPf(TPcut=0.1)
    Tset <- allP
    
    tmp <- gene_select(filenames,allP)
    genes <- tmp$genes
    allP <- tmp$allP
    allN <- tmp$allN
    tmp <- new_auc(filenames,allP,genes,flag=2,iplot=TRUE);
    
    tmp <- read.csv("result/randresult_1/DAWN/DAWN_randset1.csv")
    geneDAWN <- intersect(genes,tmp[!is.na(tmp[,"FDR"]),"Gene"])
    tmp <- new_auc(filenames,allP,geneDAWN,flag=2,iplot=TRUE);
    
    nnet <- 5
    AUC <- array(0,dim=c(8,20,3+nnet*2)) 
    FPR <- array(0,dim=c(8,20,3+nnet*2))
    TPR <- array(0,dim=c(8,20,3+nnet*2))
    
    DAWNflag <- 1
    DAWNnames <-  readLines(con <- file("result/randresult_1/DAWN/TADAdenovo_randset_1.txt","r"))   
    close(con)
    Maginames <-  readLines(con <- file("result/randresult_1/MAGI/MAGIfilenames.txt","r"))   
    close(con)
    
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
    
    nmeth <- length(filenames)
    nmeth=13
    pdf(file="plot/AUCs11.pdf",width=10,height=12)
    
    load("AUC_m_6_11")
    cols <- 1:13 #c("Black","Blue","Cyan","Brown","Coral","DarkViolet","Green","HotPink","Orange","Plum","deeppink")
    labels <- c("TADA","DAWN","MAGI","RWR-STRING","RWR-iRef","RWR-coexp","RWR-Infer","RWR-GNAT","NetP-STRING","NetP-iRef","NetP-coexp","NetP-Infer","NetP-GNAT")
    legend <- c("TADA","DAWN","MAGI","RWR-STRING","RWR-iRef","RWR-coexp","RWR-Infer","RWR-GNAT","NetP-STRING","NetP-iRef","NetP-coexp","NetP-Infer","NetP-GNAT")
    
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
    text(rep(0.2,10),Y1,labels=c("TADA","DAWN","RWR-STRING","RWR-iRef","RWR-coexp","RWR-Infer","NetP-STRING","NetP-iRef","NetP-coexp","NetP-Infer"),pos=2,col=cols,cex=1)
    
    legend("bottomright",col=cols[1:nmeth],legend=c("TADA","DAWN","RWR-STRING","RWR-iRef","RWR-coexp","RWR-Infer","NetP-STRING","NetP-iRef","NetP-coexp","NetP-Infer"),lty=rep(1,nmeth),lwd=rep(2,nmeth),cex=1.2,y.intersp=0.8)

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
    legend=c("TADA","DAWN","RWR-STRING","RWR-iRef","RWR-coexp","RWR-Infer","NetP-STRING","NetP-iRef","NetP-coexp","NetP-Infer")  
    cols <- c("Black","Blue","Cyan","Brown","Coral","DarkViolet","Green","HotPink","Orange","Plum")
    par(mai=c(2,1,1,1))
    boxplot(t(a),col=cols,names=legend,las=2,ylab="The fraction of genes",cex.axis=1.5,cex.lab=2,cex.names=2)
    
    dev.off()
    
    
    
    pdf(file="plot/subsams.pdf",width=10,height=10)
    legend=c("TADA","DAWN","RWR-STRING","RWR-iRef","RWR-coexp","RWR-Infer","NetP-STRING","NetP-iRef","NetP-coexp","NetP-Infer")  
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
    legend=c("TADA","DAWN","RWR-STRING","RWR-iRef","NetP-STRING","NetP-iRef")  
    legend <- rep(legend,3)
    par(bg=NA,mai=c(4,2,2,2))
    violins(partL,names=legend,col=rep(c('purple','lightblue','lightgreen','red','orange','yellow'),3), main=main,xlab=xlab,ylab=ylab,deciles=FALSE,connect="",at=c(1:6,8:13,15:20))
    abline(v=7,col=1,lwd=3)
    abline(v=14,col=1,lwd=3)
    mtext(text=c("395 -> 1186 trios","791 -> 2372 trios","1186 -> 3558 trios"), side=3, font=2, line=0.5,at=c(3.5,10.5,17.5))
    
}

batch_randset4_4_24 <- function(){
    source("enrichana.R")
    betaV <- c(0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0)
    result <- list()
    nmeth <- 4
    legend=c("STRING","iRef","coexp","Infer")
    for(netk in 1:4){
        result[[netk]] <- randset4_all(netk,betaV,TPcut=0.1,flag=1)
    }
    
    result1 <- list()
    nmeth <- 4
    legend=c("STRING","iRef","coexp","Infer")
    for(netk in 1:4){
        result1[[netk]] <- leaveone4_all(netk,betaV,TPcut=0.1)
    }

    load("Recurrent")
    load("Leave-one")  
    source("enrichana.R")
    betaV <- c(0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0)
    nmeth <- 4
    legend=c("STRING","iRef","coexp","Infer")
    N <- c(153,134,115,109) ## FDR 0.3
    ## N <- c(72,73,67,66) ## FDR 0.1
    for(netk in 1:4){ 
        b <- result[[netk]]$sM
        a <- matrix(0,length(betaV),nmeth)
        for(i in 1:nmeth){
            a[,i] <- sqrt(colMeans(b[,,i]) * (result1[[netk]][,i]/N[netk]))
        }
        
        for(i in 1:nmeth){
            ymax <- max(a)
            ymin <- min(a)
            if(i==1){
                plot(a[,i],col=i,type="b",main="",xlab="beta",ylab="fraction of recurrent",ylim=c(ymin,ymax),xaxt="n")
            }else{
                lines(a[,i],col=i,type="b")
            }
        }
        axis(1, at=1:length(betaV), labels = betaV)
        legend("bottomleft",col=1:nmeth,legend=paste(legend[netk],1:4,sep="_"),lty=rep(1,nmeth),lwd=rep(2,nmeth),cex=1.2,y.intersp=0.8)
    }
    
    
    for(netk in 1:4){ 
        filename <- paste("plot/recurrent",legend[netk],netk,".pdf",sep="")
        pdf(file=filename,width=12,height=5)
        b <- result[[netk]]$sM
        a <- matrix(0,length(betaV),nmeth)
        for(i in 1:nmeth){
            a[,i] <- colMeans(b[,,i]) 
        }
        
        for(i in 1:nmeth){
            ymax <- max(a)
            ymin <- min(a)
            if(i==1){
                plot(a[,i],col=i,type="b",main="",xlab="beta",ylab="fraction of recurrent",ylim=c(ymin,ymax),xaxt="n")
            }else{
                lines(a[,i],col=i,type="b")
            }
        }
        axis(1, at=1:length(betaV), labels = betaV)
        legend("bottomleft",col=1:nmeth,legend=paste(legend[netk],1:4,sep="_"),lty=rep(1,nmeth),lwd=rep(2,nmeth),cex=1.2,y.intersp=0.8)
        dev.off()
    }
    
    for(netk in 1:4){
        filename <- paste("plot/leave-one",legend[netk],netk,".pdf",sep="")
        pdf(file=filename,width=12,height=5)
        a <- result1[[netk]]/N[netk]
        for(i in 1:nmeth){
            ymax <- max(a)
            ymin <- min(a)
            if(i==1){
                plot(a[,i],col=i,type="b",main="",xlab="beta",ylab="number of genes",ylim=c(ymin,ymax),xaxt="n")
            }else{
                lines(a[,i],col=i,type="b")
            }
        }
        axis(1, at=1:length(betaV), labels = betaV)
        legend("bottomleft",col=1:nmeth,legend=paste(legend[netk],1:4,sep="_"),lty=rep(1,nmeth),lwd=rep(2,nmeth),cex=1.2,y.intersp=0.8)
        dev.off()
    }

    
    ### network genes and not 
    source("enrichana.R")
    betaV <- c(0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0)
    retada <- list()
    for(netk in 1:4){
        retada[[netk]] <- randset4_tada(netk,betaV)
        print(netk)
    }
    save(retada,file="Recurrent_tada")
    
    load("Recurrent_tada") ### FDR cut 0.1 there is some NAN
    netstr <- c("STRING","iRef","coexp","Infer")
    for(netk in 1:4){
        filename <- paste("plot/recurrent",netstr[netk],".pdf",sep="")
        pdf(file=filename,width=10,height=7)
        a <- matrix(0,11,8)
        k=1
        for(i in 1:4){ # methods
            a[,k] <- colMeans(retada[[netk]]$TADAs1[,,i],na.rm = TRUE)
            k <- k+1
            a[,k] <- colMeans(retada[[netk]]$TADAs2[,,i],na.rm = TRUE)
            k <- k+1
            #a[,k] <- colMeans(retada[[netk]]$TADAs3[,,i])
            #k <- k+1
            #a[,k] <- colMeans(retada[[netk]]$TADAs4[,,i])
            #k <- k+1
        }
        boxplot(a,names=c("Version1","TADA only","Version2","TADA only","Version3","TADA only","Version4","TADA only"),main=netstr[netk],ylab="Fraction of recurrent genes")
        abline(v=c(2.5,4.5,6.5))
        
        dev.off()
    }
    
    source("enrichana.R")
    betaV <- c(0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0)
    retada <- list()
    for(netk in 1:4){
        retada[[netk]] <- randset4_tada(netk,betaV,TPcut=0.1)
        print(netk)
    }
    save(retada,file="Recurrent_tada_0.3")
    
    if(FALSE){
        b <- sM
        a <- matrix(0,length(betaV),nnet)
        for(i in 1:nnet){
            a[,i] <- colMeans(b[,,i])
            ymax <- max(b)
            ymin <- min(b)
            if(i==1){
                plot(a[,i],col=i,type="b",main="",xlab="beta",ylab="recurrent",ylim=c(ymin,ymax),xaxt="n")
            }else{
                lines(a[,i],col=i,type="b")
            }
        }
        axis(1, at=1:length(betaV), labels = betaV)
        legend("bottomleft",col=1:nnet,legend=c("NetP-STRING","NetP-iRef","NetP-coexp","NetP-Infer"),lty=rep(1,nnet),lwd=rep(2,nnet),cex=1.2,y.intersp=0.8)
        
        
        
        b <- TADAs2
        a <- matrix(0,10,nnet)
        for(i in 1:nnet){
            a[,i] <- colMeans(b[,,i])
            ymax <- 0.5
            ymin <- min(b)
            lines(a[,i],col=i+nnet,type="p",lwd=3)
        }
        legend("bottomleft",col=1:(2*nnet),legend=c("NetP-STRING","NetP-iRef","NetP-coexp","NetP-Infer","NetP-STRING_e","NetP-iRef_e","NetP-coexp_e","NetP-Infer_e"),lty=rep(1,2*nnet),lwd=rep(2,2*nnet),cex=1.2,y.intersp=0.8)
        
        
        b <- TADAs1
        b1 <- TADAs2
        legend=c("NetP-STRING","NetP-iRef","NetP-coexp","NetP-Infer")
        for(i in 1:nnet){
            a1 <- colMeans(b[,,i])
            a2 <- colMeans(b1[,,i])
            sd1 <- apply(b[,,i],2,sd)
            sd2 <- apply(b1[,,i],2,sd)
            ymax <- 0.4
            ymin <- min(min(b),min(b1))
            plot(1:10,a1,col=i,type="b",main=legend[i],xlab="beta",ylab="recurrent",ylim=c(ymin,ymax),xaxt="n")
            x <- 1:10
            y <- a1
            sd <- sd1
            segments(x,y-sd,x,y+sd)
            epsilon <- 0.004
            segments(x-epsilon,y-sd,x+epsilon,y-sd)
            segments(x-epsilon,y+sd,x+epsilon,y+sd)
            
            x <- 1:10
            y <- a2
            sd <- sd2
            lines(1:10,a2,col=i,type="b") 
            segments(x,y-sd,x,y+sd)
            epsilon <- 0.004
            segments(x-epsilon,y-sd,x+epsilon,y-sd)
            segments(x-epsilon,y+sd,x+epsilon,y+sd)
        }
        axis(1, at=1:10, labels = betaV)
    }
}

best_parameter_6_7 <- function(){

    ### network genes and not 
    source("enrichana.R")
    betaV <- c(0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0)
    retada <- list()
    for(netk in 1:5){
        retada[[netk]] <- best_para(netk,betaV,TPcut=0.1)
        ##retada[[netk]] <- randset4_tada(netk,betaV,TPcut=0.1)
        print(netk)
    }
    save(retada,file="Recurrent_tada_6_7")
       
    netstr <- c("STRING","iRef","coexp","Infer","GNAT")
    methods <- c("Degree-Translation","Degree-Scaling","Edge-betweenness-Translation","Edge-betweenness-Scaling")
    for(netk in 1:5){
        filename <- paste("plot/resub",netstr[netk],".pdf",sep="")
        pdf(file=filename,width=12,height=12)
        a <- matrix(0,11,8)
        b <- matrix(0,11,8)
        par(mfrow=c(2,2),oma=c(1,1,2,1))
        k=1
        for(i in 1:4){ # methods
            a[,k] <- colMeans(retada[[netk]]$TADAs1[,,i],na.rm = TRUE)
            b[,k] <- apply(retada[[netk]]$TADAs1[,,i],2,sd,na.rm=TRUE)
            k <- k+1
            k <- K+1
            plot(a[,k-2],col=1,type="b",ylim=c(min(a[,(k-2):(k-1)]),max(a[,(k-2):(k-1)])),xaxt="n",xlab="Beta",ylab="The fraction of recurrent",main=methods[i])
            axis(1, at=1:length(betaV), labels = betaV)
        }
        title(main=netstr[netk],outer=T)
        dev.off()
    }
    
    ## output the best parameters with the largest recurrent values in Tset
    load("Recurrent_tada_6_7")
    betaV <- c(0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0)
  
    nnet <- 5
    nmeth <- 4
    BestP <- array(0,c(nnet,nmeth,3))
    BePair <- matrix(0,nnet,2)
    
    for(netk in 1:5){    
        a <- matrix(0,11,4)
        b <- matrix(0,11,4)
        for(i in 1:4){ # methods
            a[,i] <- colMeans(retada[[netk]]$TADAs1[,,i],na.rm = TRUE)
            b[,i] <- apply(retada[[netk]]$TADAs1[,,i],2,sd,na.rm=TRUE)
         
            c <- a[,i]/b[,i]
            BestP[netk,i,1] <- betaV[which.max(a[,i])]
            BestP[netk,i,2] <- max(c)
            BestP[netk,i,3] <- max(a[,i])
        }
        BePair[netk,1] <- which.max(BestP[netk,,3])## method
        BePair[netk,2] <- BestP[netk,BePair[netk,1],1] ## beta
    }
    
    

}

batch_randset4_6_2 <- function(){
    source("enrichana.R")
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

leaveone4_one_6_9 <- function(netflag,meth,beta,Tset,path,HiV){
    
    n <- length(Tset)
    tmp <- 0
    for(i in 1:n){
        filename <- paste(path,"CRFresult_",beta,"rand2_",i,"LBP_",netflag,".txt",sep="")
        result <- read.table(filename)
        cutoff <- total_genes_cutoff(filename,Tset)
        subs <- Tset[i] %in% result[result[,2]>=cutoff,1]
        tmp <- tmp + subs
    }
  
    tmp
}

add_sdbar <- function(x,y,sd,col){
    
    segments(x,y-sd,x,y+sd,col=col)
    epsilon <- 0.004
    segments(x-epsilon,y-sd,x+epsilon,y-sd,col=col)
    segments(x-epsilon,y+sd,x+epsilon,y+sd,col=col)
    
}

best_performance_6_11 <- function(){
    
    source("enrichana.R")
    
    TPcut=0.1
    options(stringsAsFactors=FALSE)
    netstr <- c("STRING/","iRef/","coexp/","Infer/","GNAT/")
    netnum <- c(31,6,7,8,21)
    TADAFile <- "ASD/TADAresult/TADAdenovo_ASD4_16.csv"
    asdall <- read.csv(TADAFile)
    Tset <- asdall[asdall[,"qvalue.dn"]<=TPcut,1]
    
    ### leave one mutation ==========================================
    
    nLeo <- matrix(0,13,2)
    
    ### TADA method
    ntada <- 0
    for(i in 1:length(Tset)){
        filename <- paste("ASD/TADAresult/leaveone4/TADAdenovo_rand2_",i,".csv",sep="")
        tmp <- read.csv(filename)
        oner <- tmp[match(Tset[i],tmp[,1]),"qvalue.dn"]  <= TPcut
        ntada <- ntada + oner
    }
    ### control data result  
    ntadac <- 0
    filename <- "ASD/TADAresult/TADAdenovo_control.csv"
    tmp <- read.csv(filename)
    subs <- tmp[match(Tset,tmp[,1]),"qvalue.dn"] <= TPcut
    subs[is.na(subs)] <- FALSE
    ntadac <- sum(subs)
    
    nLeo[1,] <- c(ntada,ntadac) 
    
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
                nLeo[j+3,1] <- nLeo[j+3,1] + oner
                if(oner & (ctmp[match(Tset[i],ctmp[,1]),2] >= cutoff)) nLeo[j+3,2] <- nLeo[j+3,2] + 1
            }
        }       
    }
    
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
                nLeo[j+8,1] <- nLeo[j+8,1] + oner
                if(oner & (ctmp[match(Tset[i],ctmp[,1]),2] >= cutoff)) nLeo[j+8,2] <- nLeo[j+8,2] + 1
            }
        }       
    }
    
    ### plot figures
    pdf(file="plot/FDRs6_12.pdf",width=8,height=10)
    par(mai=c(4,2,2,2))
    legend=c("TADA","DAWN","MAGI","RWR-STRING","RWR-iRef","RWR-coexp","RWR-Infer","RWR-GNAT","NetP-STRING","NetP-iRef","NetP-coexp","NetP-Infer","NetP-GNAT")  
    barplot(t(nLeo), names.arg=legend, beside=T, ylab="#genes (FDR<=0.1)", 
            cex.names=1.5, las=2, ylim=c(0,max(nLeo)+3), col=c("red","darkblue"),cex.lab=1.5,cex.axis=1.25)
    #box(bty="l")
    legend("topleft",col=c("darkblue","red"),legend=c("ASD case","SSC control"),cex=1.3,fill=c("red","darkblue")) 
    dev.off()    
    
    
    
    #######  independent samples  ===================================
    source("enrichana.R")
    TPcut=0.1
    sM1 <- matrix(0,20,4)
    options(stringsAsFactors=FALSE)
    netstr <- c("STRING/","iRef/","coexp/","Infer/")
    netnum <- c(31,6,7,8)
    TADAFile <- "ASD/TADAresult/TADAdenovo_ASD4_16.csv"
    asdall <- read.csv(TADAFile)
    Tset <- asdall[asdall[,"qvalue.dn"]<=TPcut,1]
    betaV <- rep(0.5,4)
    methV <- rep(1,4)
    j=3
    
    for(i in 1:20){
        for(netk in 1:4){
            beta <- betaV[netk]
            k <- methV[netk]
            
            dirstr <- paste("result/randresult4_",k,"/",sep="")
            filename <- paste(dirstr,netstr[netk],"CRFresult_",beta,"part",j,"_",i,"LBP_",netnum[netk],".txt",sep="")
            cutoff <- total_genes_cutoff(filename,Tset)
            allP <- asdall[,1]
            tmp <- read.table(filename)
            allP <- intersect(allP,tmp[,1])
            
            onefile <- paste("ASD/TADAresult/randset4/TADAdenovo_ASD2sub",j,"derest",i,".csv",sep="")
            oneresult <- read.csv(onefile)
            TADAfile <- paste("ASD/TADAresult/randset4/TADAdenovo_ASD1sub",j,"depart",i,".csv",sep="")
            subresult <- read.csv(TADAfile)
            
            genes <- tmp[tmp[,2] >= cutoff ,1]
            genes <- intersect(genes,allP)
            tmpg <- intersect(genes,oneresult[,1])
            sM1[i,netk] <- sum(rowSums(oneresult[match(tmpg,oneresult[,1]),c("dn.LoF","dn.mis3")])>0)/length(genes)
        }
    }
    
    ntada <- rep(0,20)
    j=3
    for(i in 1:20){
        onefile <- paste("ASD/TADAresult/randset4/TADAdenovo_ASD2sub",j,"derest",i,".csv",sep="")
        oneresult <- read.csv(onefile)
        TADAfile <- paste("ASD/TADAresult/randset4/TADAdenovo_ASD1sub",j,"depart",i,".csv",sep="")
        subresult <- read.csv(TADAfile)
        
        tmpg <- subresult[subresult[,"qvalue.dn"] <= TPcut,1]
        ntada[i] <- sum(rowSums(oneresult[match(tmpg,oneresult[,1]),c("dn.LoF","dn.mis3")])>0)/length(tmpg)
    }
    
    sM1 <- cbind(ntada,sM1)
    pdf(file="plot/subsams4_28.pdf",width=7,height=7)
    legend=c("TADA","NetP-STRING","NetP-iRef","NetP-coexp","NetP-Infer")  
    cols <- c("Black","Blue","Cyan","Brown","Coral")
    par(mai=c(2,1,1,1))
    boxplot(sM1,col=cols,names=legend,las=2,ylab="The fraction of genes",cex.axis=1.5,cex.lab=2,cex.names=2)
    dev.off()

}

randset4_all <- function(netk,betaV,TPcut=0.3,flag=1){
    
    options(stringsAsFactors=FALSE)
    netstr <- c("STRING/","iRef/","coexp/","Infer/")
    netnum <- c(31,6,7,8)
#     TADAFile <- "ASD/TADAresult/TADAdenovo_ASD4_16.csv"
#     asdall <- read.csv(TADAFile)
#     Tset <- asdall[asdall[,"qvalue.dn"]<=TPcut,1]
    
    sM <- array(-1,dim=c(20,length(betaV),4))
    sM1 <- sM
    TADAs1 <- sM
    TADAs2 <- sM
    nnet <- 4
    
    for(kk1 in 1:length(betaV)){
        beta <- betaV[kk1]
#         filenames <- 1:nnet
#         for(kk in 1:nnet){
#             dirstr <- paste("result/randresult4_",kk,"/",sep="")
#             filenames[kk] <- paste(dirstr,netstr[netk],"CRFresult_",beta,"ASD4_16LBP_",netnum[netk],".txt",sep="")
#         }
#         
# #         if(flag==1){
# #             cutoff <- total_genes_cutoff(filenames,Tset,muts="score")
# #         }else if(flag==2){
# #             cutoff <- rep(0.5,4)
# #         }
#         
#         allP <- asdall[,1]
#         for(kk in 1:nnet){
#             tmp <- read.table(filenames[kk])
#             allP <- intersect(allP,tmp[,1])
#         }
        j=3
        for(i in 1:20){
                onefile <- paste("ASD/TADAresult/randset4/TADAdenovo_ASD2sub",j,"derest",i,".csv",sep="")
                oneresult <- read.csv(onefile)
                
                TADAfile <- paste("ASD/TADAresult/randset4/TADAdenovo_ASD1sub",j,"depart",i,".csv",sep="")
                subresult <- read.csv(TADAfile)
            
                Tset <- subresult[subresult[,"qvalue.dn"] <= TPcut,1]
                
                result <- list()
                for(k in 1:nnet){
                    dirstr <- paste("result/randresult4_",k,"/",sep="")
                    filename <- paste(dirstr,netstr[netk],"CRFresult_",beta,"part",j,"_",i,"LBP_",netnum[netk],".txt",sep="")
                    cutoff <- total_genes_cutoff(filename,Tset)
                    
                    tmp <- read.table(filename)
                    result[[k]] <- tmp[!duplicated(tmp[,1]),]
                    genes <- result[[k]][result[[k]][,2] >= cutoff,1]
                    
                    tmpg <- intersect(genes,oneresult[,1])
                    sM[i,kk1,k] <- sum(rowSums(oneresult[match(tmpg,oneresult[,1]),c("dn.LoF","dn.mis3")])>0)/length(genes)
                    sM1[i,kk1,k] <- sum(rowSums(oneresult[match(tmpg,oneresult[,1]),c("dn.LoF","dn.mis3")])>0)
                }
        }
    }
    
    list(sM=sM,sM1=sM1)
    
}

best_para <- function(netk,betaV,TPcut=0.1,flag=1){
    
    options(stringsAsFactors=FALSE)
    netstr <- c("STRING/","iRef/","coexp/","Infer/","GNAT/")
    netnum <- c(31,6,7,8,21)
    nmeth <- 4
    
    TADAs1 <- array(-1,dim=c(20,length(betaV),nmeth))
    
    for(kk1 in 1:length(betaV)){
        beta <- betaV[kk1]  
        j=3
        for(i in 1:20){
            onefile <- paste("ASD/TADAresult/randset4/TADAdenovo_ASD2sub",j,"derest",i,".csv",sep="")
            oneresult <- read.csv(onefile)
            
            TADAfile <- paste("ASD/TADAresult/randset4/TADAdenovo_ASD1sub",j,"depart",i,".csv",sep="")
            subresult <- read.csv(TADAfile)
            Tset <- subresult[subresult[,"qvalue.dn"] <= TPcut,1]
            
            result <- list()
            for(k in 1:nmeth){
                dirstr <- paste("result/randresult4_",k,"/",sep="")
                filename <- paste(dirstr,netstr[netk],"CRFresult_",beta,"part",j,"_",i,"LBP_",netnum[netk],".txt",sep="")
                cutoff <- total_genes_cutoff(filename,Tset)
                
                tmp <- read.table(filename)
                result[[k]] <- tmp[!duplicated(tmp[,1]),]
                
                genes <- result[[k]][result[[k]][,2] >= cutoff,1]
                g1 <- intersect(genes,oneresult[,1])
                
                TADAs1[i,kk1,k] <- sum(rowSums(oneresult[match(g1,oneresult[,1]),c("dn.LoF","dn.mis3")])>0)/length(g1)
            }
        }
    }
    
    list(TADAs1=TADAs1)
    
}

randset4_tada <- function(netk,betaV,TPcut=0.1,flag=1){
    
    options(stringsAsFactors=FALSE)
    netstr <- c("STRING/","iRef/","coexp/","Infer/","GNAT/")
    netnum <- c(31,6,7,8,21)
    nmeth <- 4
    
    TADAs1 <- array(-1,dim=c(20,length(betaV),nmeth))
    TADAs2 <- TADAs1
    TADAs3 <- TADAs1
    TADAs4 <- TADAs1
    
    for(kk1 in 1:length(betaV)){
        beta <- betaV[kk1]  
        j=3
        for(i in 1:20){
            onefile <- paste("ASD/TADAresult/randset4/TADAdenovo_ASD2sub",j,"derest",i,".csv",sep="")
            oneresult <- read.csv(onefile)
            
            TADAfile <- paste("ASD/TADAresult/randset4/TADAdenovo_ASD1sub",j,"depart",i,".csv",sep="")
            subresult <- read.csv(TADAfile)
            Tset <- subresult[subresult[,"qvalue.dn"] <= TPcut,1]
            
            result <- list()
            for(k in 1:nmeth){
                dirstr <- paste("result/randresult4_",k,"/",sep="")
                filename <- paste(dirstr,netstr[netk],"CRFresult_",beta,"part",j,"_",i,"LBP_",netnum[netk],".txt",sep="")
                cutoff <- total_genes_cutoff(filename,Tset)
                
                tmp <- read.table(filename)
                result[[k]] <- tmp[!duplicated(tmp[,1]),]
                
                genes <- result[[k]][result[[k]][,2] >= cutoff,1]
                
                #genes <- intersect(genes,subresult[,1])
                g1 <- intersect(genes,Tset)
                #supn <- max(match(g1,subresult[,1]))
                #g2 <- setdiff(subresult[1:supn,1],g1)
                g2 <- setdiff(Tset,g1)
                g3 <- intersect(g2,result[[k]][,1])
                g4 <- intersect(g2,g3)
                TADAs1[i,kk1,k] <- sum(rowSums(oneresult[match(g1,oneresult[,1]),c("dn.LoF","dn.mis3")])>0)/length(g1)
                TADAs2[i,kk1,k] <- sum(rowSums(oneresult[match(g2,oneresult[,1]),c("dn.LoF","dn.mis3")])>0)/length(g2)
                TADAs3[i,kk1,k] <- sum(rowSums(oneresult[match(g3,oneresult[,1]),c("dn.LoF","dn.mis3")])>0)/length(g3)
                TADAs4[i,kk1,k] <- sum(rowSums(oneresult[match(g4,oneresult[,1]),c("dn.LoF","dn.mis3")])>0)/length(g4)
            }
        }
    }
    
    list(TADAs1=TADAs1,TADAs2=TADAs2,TADAs3=TADAs3,TADAs4=TADAs4)
    
}

leaveone4_all <- function(netk,betaV,TPcut=0.3){
    options(stringsAsFactors=FALSE)
    source("enrichana.R")
    netstr <- c("STRING/","iRef/","coexp/","Infer/")
    netnum <- c(31,6,7,8)
    TADAFile <- "ASD/TADAresult/TADAdenovo_ASD4_16.csv"
    asdall <- read.csv(TADAFile)
    Tset0 <- asdall[asdall[,"qvalue.dn"]<=TPcut,1]
    Tset <- Tset0
    
    sM <- array(0,dim=c(length(betaV),4))
    nmeth <- 4
    
    filenames <- 1:nmeth
    for(kk in 1:nmeth){
        filenames[kk] <- paste("result/leaveone4result_",kk,"/",netstr[netk],"CRFresult_",0,"rand2","_",1,"LBP_",netnum[netk],".txt",sep="")
        tmp <- read.table(filenames[kk])
        Tset <- intersect(Tset,tmp[,1])
    }
    subs <- match(Tset,Tset0)
    print(length(subs))
    
    rank1 <- array(0,dim=c(length(betaV),nmeth,length(subs)))
    for(kk1 in 1:length(betaV)){
        beta <- betaV[kk1]
        ksub <- 1
        for(i in subs){                  
            result <- list()
            for(k in 1:nmeth){
                dirstr <- paste("result/leaveone4result_",k,"/",sep="")
                filename <- paste(dirstr,netstr[netk],"CRFresult_",beta,"rand2","_",i,"LBP_",netnum[netk],".txt",sep="")
                
                cutoff <- total_genes_cutoff(filename,Tset)
                
                if(grepl(".csv",filename)){ tmp <- read.csv(filename); 
                }else if(grepl(".txt",filename)){tmp <- read.table(filename);}
                result[[k]] <- tmp[!duplicated(tmp[,1]),]
                
                #tmp <- result[[k]][match(Tset0[i],result[[k]][,1]),2]  >= cutoff
                tmp <- match(Tset0[i],result[[k]][,1])  <= max(subs)
                sM[kk1,k] <-  sM[kk1,k] + tmp
                
                rank1[kk1,k,ksub] <- match(Tset0[i],result[[k]][,1])
            }
            ksub <- ksub + 1
        }
        print(kk1)
    }
    
  list(sM=sM,rank1=rank1,subs=subs)

}

randset4_boxplot_4_15 <- function(dirstr,betaV){
    
    options(stringsAsFactors=FALSE)
    source("enrichana.R")
    netstr <- c("STRING/","iRef/","coexp/","Infer/")
    netnum <- c(31,6,7,8)
    TADAFile <- "ASD/TADAresult/TADAdenovo_ASD4_16.csv"
    asdall <- read.csv(TADAFile)
    Tset <- asdall[asdall[,"qvalue.dn"]<=0.3,1]
    
    #betaV <- c(0,0.1,0.2,0.3,0.4,0.5,1,2,5,10)
    #betaV <- c(0.0001,0.005,0.01,0.02,0.05,0.1,0.2,0.5,1,2,10)
    #betaV <- c(0.0001,0.005,0.01,0.02,0.05,0.1,0.2,0.5,1,10,100,1000)
    #betaV <- c(0.01,0.02,0.05,0.1,0.5,1,2,3,4,5,6,7,10)
    sM <- array(-1,dim=c(20,length(betaV),4))
    sM1 <- sM
    TADAs1 <- sM
    TADAs2 <- sM
    nnet <- 4
    
    
    for(kk1 in 1:length(betaV)){
        beta <- betaV[kk1]
        filenames <- 1:nnet
        for(kk in 1:nnet){
            filenames[kk] <- paste(dirstr,netstr[kk],"CRFresult_",beta,"ASD4_16LBP_",netnum[kk],".txt",sep="")
        }
        cutoff <- total_genes_cutoff(filenames,Tset,muts="score")
        #cutoff <- rep(0.51,4)
        allP <- asdall[,1]
        for(kk in 1:nnet){
            tmp <- read.table(filenames[kk])
            allP <- intersect(allP,tmp[,1])
        }
        
    for(i in 1:20){
        for(j in 3){
            
            onefile <- paste("ASD/TADAresult/randset4/TADAdenovo_ASD2sub",j,"derest",i,".csv",sep="")
            oneresult <- read.csv(onefile)
            
            TADAfile <- paste("ASD/TADAresult/randset4/TADAdenovo_ASD1sub",j,"depart",i,".csv",sep="")
            subresult <- read.csv(TADAfile)
            
            mfiles <- 1:nnet
            for(k in 1:nnet){
                mfiles[k] <- paste(dirstr,netstr[k],"CRFresult_",beta,"part",j,"_",i,"LBP_",netnum[k],".txt",sep="")
            }
            #filenames <- c(TADAfile,mfiles)
            filenames <- mfiles
            
            result <- list()
            for(k in 1:length(filenames)){
                    filename <- filenames[k]
                    if(grepl(".csv",filename)){ tmp <- read.csv(filename); 
                    }else if(grepl(".txt",filename)){tmp <- read.table(filename);}
                    result[[k]] <- tmp[!duplicated(tmp[,1]),]
                    
                    if(grepl(".csv",filename)){
                        genes <- result[[k]][result[[k]][,"score"] >= cutoff[k],1] ## TPcut,1]
                    }else{ genes <- result[[k]][result[[k]][,2] >= cutoff[k],1] } ##  !!!<= TPcut,1] }
                    
                    tmpgenes <- intersect(genes,subresult[,1])
                    g1 <- intersect(tmpgenes,subresult[,1])
                    g2 <- setdiff(subresult[1:max(match(tmpgenes,subresult[,1])),1],tmpgenes)
                    TADAs1[i,kk1,k] <- sum(rowSums(oneresult[match(g1,oneresult[,1]),c("dn.LoF","dn.mis3")])>0)/length(g1)
                    TADAs2[i,kk1,k] <- sum(rowSums(oneresult[match(g2,oneresult[,1]),c("dn.LoF","dn.mis3")])>0)/length(g2)
                    
                    genes <- intersect(genes,allP)
                    
                    #oneg1 <- oneresult[oneresult[,"qvalue.dn"] <=0.1,1]
                    tmpg <- intersect(genes,oneresult[,1])
                    sM[i,kk1,k] <- sum(rowSums(oneresult[match(tmpg,oneresult[,1]),c("dn.LoF","dn.mis3")])>0)/length(genes)
                    sM1[i,kk1,k] <- sum(rowSums(oneresult[match(tmpg,oneresult[,1]),c("dn.LoF","dn.mis3")])>0)
            }
        }
    }
    }
    

    

    list(sM=sM,sM1=sM1,TADAs1=TADAs1,TADAs2=TADAs2)
    
}

leaveone4_boxplot_4_20 <- function(){
    
    options(stringsAsFactors=FALSE)
    source("enrichana.R")
    
    #dirstr <- "result/leaveone4result/"
    dirstr <- "result/leaveone4result_1/"
    netstr <- c("STRING/","iRef/","coexp/","Infer/")
    netnum <- c(31,6,7,8)
    TADAFile <- "ASD/TADAresult/TADAdenovo_ASD4_16.csv"
    asdall <- read.csv(TADAFile)
    Tset0 <- asdall[asdall[,"qvalue.dn"]<=0.3,1]
    Tset <- Tset0
    
    sM <- array(0,dim=c(10,4))
    #betaV <- c(0,0.1,0.2,0.3,0.4,0.5,1,2,5,10)
    #betaV <- c(0.0001,0.005,0.01,0.02,0.05,0.1,0.2,0.5,1,2)
    betaV <- c(0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0)
    nnet <- 4
    
    filenames <- 1:nnet
    for(kk in 1:nnet){
        filenames[kk] <- paste("result/randresult4/",netstr[kk],"CRFresult_",0,"ASD4_16LBP_",netnum[kk],".txt",sep="")
        tmp <- read.table(filenames[kk])
        Tset <- intersect(Tset,tmp[,1])
    }
    subs <- match(Tset,Tset0)
    
    
    for(kk1 in 1:10){
        beta <- betaV[kk1]
#         filenames <- 1:nnet
#         for(kk in 1:nnet){
#             filenames[kk] <- paste("result/randresult4/",netstr[kk],"CRFresult_",beta,"ASD4_16LBP_",netnum[kk],".txt",sep="")
#         }
#         cutoff <- total_genes_cutoff(filenames,Tset,muts="")
        
        for(i in subs){                  
                result <- list()
                for(k in 1:nnet){
                    filename <- paste(dirstr,netstr[k],"CRFresult_",beta,"rand2","_",i,"LBP_",netnum[k],".txt",sep="")
                    if(grepl(".csv",filename)){ tmp <- read.csv(filename); 
                    }else if(grepl(".txt",filename)){tmp <- read.table(filename);}
                    result[[k]] <- tmp[!duplicated(tmp[,1]),]
                    
                    tmp <- result[[k]][match(Tset0[i],result[[k]][,1]),5]  <= 0.1 #>= cutoff[k]
                    #print(tmp)
                    sM[kk1,k] <-  sM[kk1,k] + tmp
                }
        }
        print(kk1)
    }
    
    a <- sM
    for(i in 1:nnet){
        ymax <- max(a)
        ymin <- min(a)
        if(i==1){
            plot(a[,i],col=i,type="b",main="",xlab="beta",ylab="recurrent",ylim=c(ymin,ymax),xaxt="n")
        }else{
            lines(a[,i],col=i,type="b")
        }
    }
    axis(1, at=1:10, labels = betaV)
    legend("bottomleft",col=1:nnet,legend=c("NetP-STRING","NetP-iRef","NetP-coexp","NetP-Infer"),lty=rep(1,nnet),lwd=rep(2,nnet),cex=1.2,y.intersp=0.8)
    
   
}

betweenness_look <- function(){
    source("Multi_net.R")
    source("CRF_build.R")
    TADAFile <- "ASD/TADAresult/TADAdenovo_ASD4_16.csv"
    tadar <- read.csv(TADAFile)
    
    #tadag <- tadar[1:79,1]
    #tadao <- setdiff(tadar[,1],tadag)
    
    tadag <- tadar[rowSums(tadar[,c("dn.LoF","dn.mis3")]) > 0,1]
    #tadao <- tadar[rowSums(tadar[,c("dn.LoF","dn.mis3")]) == 0 & tadar[,"dn.mis"] > 0 ,1]
    tadao <- setdiff(tadar[,1],tadag)
    
    ### node betweenness
    filenames <- c("data/Network_betweenness/Betweenness_node_STRINGm.txt","data/Network_betweenness/Betweenness_node_iRef.txt","data/Network_betweenness/Betweenness_node_Infer.txt","data/Network_betweenness/Betweenness_node_coexp.txt","data/Network_betweenness/Betweenness_node_MAGI.txt")
    
    netflags <- c(3,6,8,7,20)
    bes1 <- list()
    bes2 <- list()
    
    edn1 <- list()
    edn2 <- list()
    
    bes <- list()
    edn <- list()
    netmr <- list()
    netnr <- list()
    for(i in 1:length(filenames)){
        tmp <- read.table(filenames[i])
        betw <- tmp[,2]#/max(tmp[,2])
        netmr[[i]] <- tmp
        
        g1 <- intersect(tadag,tmp[,1])
        den1 <- betw[match(g1,tmp[,1])]
        g2 <- intersect(tadao,tmp[,1])
        den2 <- betw[match(g2,tmp[,1])]
        
        bes1[[i]] <- den1
        bes2[[i]] <- den2
        
        netflag <- netflags[i]
        net <- build_net(netflag,"")
        net$matrix[net$matrix > 0 ] <- 1
        ednw <- rowSums(net$matrix)
        
        if(netflag==3){
            mapT <- as.matrix(read.delim("Fmap0121.txt",header=FALSE,sep="\t"))
            net$node <- mapT[match(net$node,mapT[,1]),2]
        }
        edn1[[i]] <- ednw[match(g1,net$node)]
        edn2[[i]] <- ednw[match(g2,net$node)]
        
        bes[[i]] <- betw
        edn[[i]] <- ednw
        netnr[[i]] <- net$node
    }
    
    ### TADA score and betweeness number
    cormT <- matrix(0,5,3)
    for(i in 1:5){
        tmpg <- intersect(netnr[[i]],tadar[,1])
        #print("The number of genes ")
        cormT[i,3] <- length(tmpg)
        x <- 1 - tadar[match(tmpg,tadar[,1]),"qvalue.dn"]
        y <- bes[[i]][match(tmpg,netmr[[i]][,1])]
        #plot(x,y,ylim=c(0,10000))
        cormT[i,1] <- cor(x,y)
        
        z <- edn[[i]][match(tmpg,netnr[[i]])]
        #print("degree")
        cormT[i,2] <- cor(x,z)
        #plot(x,z,ylim=c(0,100))
    }
    colnames(cormT) <- c("tada_be","tada_de","n.gene")
    rownames(cormT) <- c("STRING","iRef","Infer","coexp","MAGI")
    
    ### cor difference
    corn <- matrix(0,5,4)
    for(i in 1:5){
        tmpg <- intersect(netnr[[i]],tadar[,1])
        tmpg1 <- intersect(tmpg,tadag)
        tmpg2 <- intersect(tmpg,tadao)
        
        x1 <- 1 - tadar[match(tmpg1,tadar[,1]),"qvalue.dn"]
        y1 <- bes[[i]][match(tmpg1,netmr[[i]][,1])]
        corn[i,1] <- cor(x1,y1)
        
        x2 <- 1 - tadar[match(tmpg2,tadar[,1]),"qvalue.dn"]
        y2 <- bes[[i]][match(tmpg2,netmr[[i]][,1])]
        corn[i,2] <- cor(x2,y2)
        
        x3 <- 1 - tadar[match(tmpg1,tadar[,1]),"qvalue.dn"]
        y3 <- edn[[i]][match(tmpg1,netnr[[i]])]
        corn[i,3] <- cor(x3,y3)
        
        x4 <- 1 - tadar[match(tmpg2,tadar[,1]),"qvalue.dn"]
        y4 <- edn[[i]][match(tmpg2,netnr[[i]])]
        corn[i,4] <- cor(x4,y4)        
    }
    colnames(corn) <- c("ta_be_lof_dmis","ta_be_oth","ta_de_lof_dmis","ta_de_oth")
    rownames(corn) <- c("STRING","iRef","Infer","coexp","MAGI")

    ### mean difference test
    meanT <- matrix(0,5,6)
    for(i in 1:5){
        a <- list()
        a[[1]] <- bes1[[i]]
        a[[2]] <- bes2[[i]]
        tmp <- t.test(a[[1]],a[[2]])
        meanT[i,1] <- tmp$p.value
        meanT[i,2] <- tmp$estimate[1]
        meanT[i,3] <- tmp$estimate[2]
        
        b <- list()
        b[[1]] <- edn1[[i]]
        b[[2]] <- edn2[[i]]
        tmp <- t.test(b[[1]],b[[2]])
        meanT[i,4] <- tmp$p.value
        meanT[i,5] <- tmp$estimate[1]
        meanT[i,6] <- tmp$estimate[2]        
    }    
    
    colnames(meanT) <- c("pvalue1","be_lof_dmis","be1","pvalue2","de_lof_dmis","de2")
    rownames(meanT) <- c("STRING","iRef","Infer","coexp","MAGI")
    
    
    a <- read.delim("data/Network_betweenness/STR_centralities",sep=" ")
    a <- a[,c(2,5)]
    source("Network_analysis.R")
    a[,1] <- mapping_to(a[,1])
    
    tmpg <- intersect(netmr[[1]][,1],a[,1])
    x <- a[match(tmpg,a[,1]),2] ## large
    y <- bes[[1]][match(tmpg,netmr[[1]][,1])] # small
    cor(x,y)
    
    
    b <- read.delim("data/Network_betweenness/STRING_node_degree.txt",sep=" ")
    b <- b[,2:3]
    b <- b[!is.na(b[,1]),]
    tmpg <- intersect(b[,1],netnr[[1]])
    x <- b[match(tmpg,b[,1]),2]
    y <- edn[[1]][match(tmpg,netnr[[1]])]
    
        
    
    
    
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
    legend=c("TADA","DAWN","RWR-STRING","RWR-iRef","RWR-coexp","RWR-Infer","NetP-STRING","NetP-iRef","NetP-coexp","NetP-Infer")  
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

#     legend("bottomright",col=cols[subs],legend=c("TADA","DAWN","RWR-STRING","RWR-coexp","RWR-Infer","NetP-STRING","NetP-iRef","NetP-coexp","NetP-Infer"),lty=rep(1,nmeth),lwd=rep(2,nmeth),cex=1.2,y.intersp=0.8)
#     
#     library(caroline)
#     plot.new()
#     par(bg=NA,mai=c(5,2,2,2))
#     legend=c("TADA","DAWN","RWR-STRING","RWR-iRef","RWR-coexp","RWR-Infer","NetP-STRING","NetP-iRef","NetP-coexp","NetP-Infer")  
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
    legend=c("TADA","DAWN","RWR-STRING","RWR-iRef","RWR-coexp","RWR-Infer","NetP-STRING","NetP-iRef","NetP-coexp","NetP-Infer")  
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
    text(rep(0.2,11),Y1,labels=c("TADA","DAWN","NetP-STRING","NetP-iRef","NetP-coexp","NetP-Infer","overlapexp","jointexp","overlapexp_PPI","jointexp_PPI","Overlap_ppi"),pos=2,col=cols,cex=1.3)
 
    
    
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

allPf <- function(TPcut,TADAFile="ASD/TADAresult/TADAdenovo_ASD4_16.csv"){
    
    source("enrichana.R")
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
    legend=c("TADA","DAWN","RWR-STRING","RWR-iRef","NetP-STRING","NetP-iRef")  
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
        legend=c("TADA","DAWN","MAGI","RWR-STRING","RWR-iRef","RWR-coexp","RWR-Infer","RWR-GNAT","NetP-STRING","NetP-iRef","NetP-coexp","NetP-Infer","NetP-GNAT") ###
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
    legend(pos,col=(1:nmeth),legend=c("TADA","DAWN","RWR","NetP"),lty=rep(1,nmeth),lwd=rep(1,nmeth))
    
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
    legend("bottomright",col=(1:nmeth),legend=c("TADA","DAWN","RWR","NetP"),lty=rep(1,nmeth),lwd=rep(1,nmeth))

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
    cols <- 1:20 ##c("Black","Blue","Cyan","Brown","Coral","DarkViolet","Green","HotPink","Orange","Plum")
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

bar_err <- function(G,main,xlab,ylab,groups=c("TADA","DAWN","RWR","NetP")){

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
    legend=c("TADA","DAWN","RWR-STRING","RWR-iRef","RWR-coexp","RWR-Infer","NetP-STRING","NetP-iRef","NetP-coexp","NetP-Infer") 
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

overlap_genes <- function(Tgenes,result,netflag=7,legend=c("TADA","DAWN","RWR-STRING","RWR-iRef","NetP-STRING","NetP-iRef")){
    
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