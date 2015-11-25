leaveone_Rank <- function(Rcut){
    
    options(stringsAsFactors=FALSE)
    TPcut=0.1
    TADAFile <- "../TADA_DAWN/result/TADAdenovo_meta_dmis.csv"
    asdall <- read.csv(TADAFile)
    asdall[is.na(asdall[,"qvalue.dn"]),"qvalue.dn"] <- 1
    Tset <- asdall[asdall[,"qvalue.dn"]<=TPcut,1]
    ### leave one mutation ==========================================
    nLeo <- matrix(0,7,2)
    
    ### TADA method
    cfile <- "../TADA_DAWN/result/TADAdenovo_control1911.csv"
    ctmp <- read.csv(cfile)
    for(i in 1:length(Tset)){
        filename <- paste("../TADA_DAWN/result/leaveone4/TADAdenovo_rand2_",i,".csv",sep="")
        tmp <- read.csv(filename)
        nLeo[1,1] <- nLeo[1,1] + sum(match(Tset[i],tmp[,1]) <= Rcut)
        nLeo[1,2] <- nLeo[1,2] + sum(match(Tset[i],ctmp[,1]) <= Rcut)
    }
    
    ## DAWN method 
    cfile <- "../TADA_DAWN/DAWN_package/DAWNcontrol1911.csv" ###  checked
    ctmp <- read.csv(cfile)
    ctmp[is.na(ctmp[,"Stage2_posterior"]),"Stage2_posterior"] <- 1
    for(i in 1:length(Tset)){
        filename <- paste("../TADA_DAWN/DAWN_package/result/leaveone4/DAWN_leaveone4_",i,".csv",sep="")
        tmp <- read.csv(filename)
        if(Tset[i] %in% tmp[,1]) nLeo[2,1] <- nLeo[2,1] + sum(match(Tset[i],tmp[,1]) <= Rcut)
        if(Tset[i] %in% ctmp[,1]) nLeo[2,2] <- nLeo[2,2] + sum(match(Tset[i],ctmp[,1]) <= Rcut)
    }
    
    ###  MAGI method
    cfile <- "../MAGI_V0.1_Alpha/mydata/control/RandomGeneList.3"
    ctmp <- read.table(cfile)
    ctmp <- ctmp[order(-as.numeric(ctmp[,2])),] 
    for(i in 1:length(Tset)){
        filename <- paste("../MAGI_V0.1_Alpha/mydata/leaveone4/RandomGeneList.",i,sep="")
        tmp <- read.table(filename)
        tmp <- tmp[order(-as.numeric(tmp[,2])), ]
        if(Tset[i] %in% tmp[,1]) nLeo[3,1] <- nLeo[3,1] + sum(match(Tset[i],tmp[,1]) <= Rcut)
        if(Tset[i] %in% ctmp[,1]) nLeo[3,2] <- nLeo[3,2] + sum(match(Tset[i],ctmp[,1]) <= Rcut)
    }
    
    ### our method
    outputstr='/ifs/scratch/c2b2/ys_lab/qh2159/Mutations/CHD/MIS/result/control/v5/MICRFresult_';
    for(j in 1:4){        
        cfile=paste(outputstr,j+1,'.txt',sep=""); #!!!
        ctmp <- read.table(cfile)
        ctmp <- ctmp[order(-as.numeric(ctmp[,2])),] 
        for(i in 1:length(Tset)){
            filename <- paste("result/leaveone4result_5/MICRFresult_",j+1,"_",i,".txt",sep="") ##!!!!!!!!
            tmp <- read.table(filename)
            tmp <- tmp[order(-as.numeric(tmp[,2])),]
            if(Tset[i] %in% tmp[,1]) nLeo[j+3,1] <- nLeo[j+3,1] + sum(match(Tset[i],tmp[,1]) <= Rcut)
            if(Tset[i] %in% ctmp[,1]) nLeo[j+3,2] <- nLeo[j+3,2] + sum(match(Tset[i],ctmp[,1]) <= Rcut)
        }       
    }
    
    nLeo
    
}

leaveone_Rankes <- function(){
    source("plot_MICRF.R")
    Vrank <- c(100,200,300,500,1000,5000)
    TPM <- matrix(0,7,length(Vrank))
    FPM <- TPM
    for(i in 1:length(Vrank)){
        tmp <- leaveone_Rank(Vrank[i]);
        TPM[,i] <- tmp[,1]
        FPM[,i] <- tmp[,2]
    }
    leaveone_plots(TPM,FPM)
}

leaveone_plots <- function(TPM,FPM){
   
    Vrank <- c(100,200,300,500,1000,5000) 
    legend <- c("TADA","DAWN","MAGI","MICRF-iRef","MICRF-HPRD","MICRF-corr","MICRF-coexp")
    cols <- c("red","green","blue","orange","darkorange","yellow","pink")
    
    pdf(file="plot/leaveone_rankes_sim.pdf",width=24,height=10)
    par(mfrow=c(1,2),mar=c(6,6,3,2))
    plot(log10(Vrank),TPM[1,],xaxt="n", type="b",col=cols[1],main="DDD cases",ylim=c(0,max(TPM)+5),xlab="rank",ylab="number of risk genes",cex.lab=2,cex.main=2,lwd=2,cex.axis=1.8)
    axis(1,at=log10(Vrank),labels=Vrank,font=1,cex.axis=1.5)
    for(i in 2:7){
        lines(log10(Vrank),TPM[i,], type="b",col=cols[i],xaxt="n",lwd=2)
    }
    legend("topleft",col=cols,legend=legend,lty=rep(1,15),lwd=rep(2,15),cex=1,y.intersp=1)
    
    plot(log10(Vrank),FPM[1,],xaxt="n", type="b",col=cols[1],main="SSC controls",ylim=c(0,max(FPM)+5),xlab="rank",ylab="number of risk genes",cex.lab=2,cex.main=2,lwd=2,cex.axis=1.8)
    axis(1,at=log10(Vrank),labels=Vrank,font=1,cex.axis=1.5)
    for(i in 2:7){
        lines(log10(Vrank),FPM[i,], type="b",col=cols[i],xaxt="n",lwd=2)
    }
    legend("topleft",col=cols,legend=legend,lty=rep(1,15),lwd=rep(2,15),cex=1,y.intersp=1)
    
    dev.off()

    ### plot figures
    for(i in 1:length(Vrank)){
        nLeo <- cbind(TPM[,i],FPM[,i])
        pdf(file=paste("plot/leaveone_",Vrank[i],".pdf",sep=""),width=10,height=10)
        par(mai=c(4,2,2,2))
        legend <- c("TADA","DAWN","MAGI","MICRF-iRef","MICRF-HPRD","MICRF-corr","MICRF-coexp")
        ords <- 1:7
        mp <- barplot(t(nLeo)[,ords], beside=T, ylab="#gold risk genes",cex.names=1.5, las=2, ylim=c(0,max(nLeo)), col=c("red","darkblue"),cex.lab=1.5,cex.axis=1.25,axes = FALSE, axisnames = FALSE)
        text(colMeans(mp), par("usr")[3], labels = legend[ords], srt = 45, adj = c(1.1,1.1), xpd = TRUE,cex=1.5)
        axis(2)
        legend("topleft",legend=c("DDD case","SSC control"),cex=1.3,fill=c("red","darkblue")) 
        dev.off()
    }
}

leaveone4_one <- function(){
    
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


## try missing data
recurrent_Rank <- function(Rcut,flag=1,TPcut=0.3){
    ## Rcut: rank cutoff
    ## flag: only risk gene or not
    options(stringsAsFactors=FALSE)
    TADAFile <- "../TADA_DAWN/result/TADAdenovo_meta_dmis.csv"
    asdall <- read.csv(TADAFile)
    asdall[is.na(asdall[,"qvalue.dn"]),"qvalue.dn"] <- 1
    Tset <- asdall[asdall[,"qvalue.dn"]<=TPcut,1]
    
    #######  independent samples  ===================================
    ReM <- matrix(0,20,8)
    j=3
    for(i in 1:20){
        #onefile <- paste("DDD_mutations/randset4/ASD2sub",j,"derest",i,".csv",sep="")
        onefile <- paste("DDD_mutations/randset4/ASD1sub",j,"depart",i,".csv",sep="")
        oneresult <- read.csv(onefile)
        
        ### TADA method
        #TADAfile <- paste("../TADA_DAWN/result/randset4/TADAdenovo_ASD1sub",j,"depart",i,".csv",sep="")
        TADAfile <- paste("../TADA_DAWN/result/randset4/TADAdenovo_ASD2sub",j,"derest",i,".csv",sep="")
        subresult <- read.csv(TADAfile)
        tmpg <- subresult[1:Rcut,1]
        if(flag == 2) tmpg <- intersect(tmpg,Tset)
        ReM[i,1] <- sum(rowSums(oneresult[match(tmpg,oneresult[,1]),c("dn.LoF","dn.mis3")])>0)
        
        ### DAWN method
        #DAWNfile <- paste("../TADA_DAWN/DAWN_package/result/randset4/DAWN_randset4_",i,".csv",sep="")
        DAWNfile <- paste("../TADA_DAWN/DAWN_package/result/randset4/DAWN_randset4rest_",i,".csv",sep="")
        tmp <- read.csv(DAWNfile)
        tmpg <- tmp[1:Rcut,"Gene"]
        if(flag == 2) tmpg <- intersect(tmpg,Tset)
        ReM[i,2] <- sum(rowSums(oneresult[match(tmpg[tmpg %in% oneresult[,1]],oneresult[,1]),c("dn.LoF","dn.mis3")])>0)
        
        ### MAGI method 
        #MAGIfile <- paste("../MAGI_V0.1_Alpha/mydata/randset4/RandomGeneList.",i,sep="")
        MAGIfile <- paste("../MAGI_V0.1_Alpha/mydata/randset4rest/RandomGeneList.",i,sep="")
        tmp <- read.table(MAGIfile)
        tmp <- tmp[order(-as.numeric(tmp[,2])), ]
        tmpg <- tmp[1:Rcut,1]
        if(flag == 2) tmpg <- intersect(tmpg,Tset)
        ReM[i,3] <- sum(rowSums(oneresult[match(tmpg[tmpg %in% oneresult[,1]],oneresult[,1]),c("dn.LoF","dn.mis3")])>0)
        
        ### our method
        path = "/ifs/scratch/c2b2/ys_lab/qh2159/Mutations/CHD/MIS/result/randresult4_5/"
        for(netk in 1:5){
            myfile <- paste(path,"MICRFresult_",netk,"_",i,".txt",sep="")
            if(file.exists(myfile)){
                tmp <- read.table(myfile)
                tmp <- tmp[order(-as.numeric(tmp[,2])), ]
                tmpg <- tmp[1:Rcut,1]
                if(flag == 2) tmpg <- intersect(tmpg,Tset)
                ReM[i,netk+3] <- sum(rowSums(oneresult[match(tmpg[tmpg %in% oneresult[,1]],oneresult[,1]),c("dn.LoF","dn.mis3")])>0)
            }
        }
    }
     
    ReM
}

recurrent_Rankes <- function(){
    setwd("/ifs/scratch/c2b2/ys_lab/qh2159/Mutations/CHD/MIS/")
    source("plot_MICRF.R")
    Vrank <- c(50,100,200,500,1000,5000)
    Rlist <- list()
    for(i in 1:length(Vrank)){
        tmp <- recurrent_Rank(Vrank[i]);
        Rlist[[i]] <- tmp
    }
    recurrent_plots(Rlist)
    
    Vrank <- c(50,100,200,500,1000,5000)
    Rlist2 <- list()
    for(i in 1:length(Vrank)){
        tmp <- recurrent_Rank(Vrank[i],2,0.1);
        Rlist2[[i]] <- tmp
    }
    recurrent_plots(Rlist2,2,0.1)  
    
}

recurrent_plots <- function(Rlist,flag=1,TPcut=0.3){
#     load("ReMlist_rank_11_22")
#     Vrank <- c(50,100,200,500,1000,5000)
#     ReM <- matrix(0,15,length(Vrank))
#     for(i in 1:length(Rlist)){
#         y <- colMeans(Rlist[[i]],na.rm=TRUE)
#         ReM[,i] <- y
#     }
#     legend <- c("TADA","DAWN","MAGI","RWR-STRING","RWR-iRef","RWR-HPRD","RWR-corr","RWR-coexp","MICRF-STRING","MICRF-iRef","MICRF-HPRD","MICRF-corr","MICRF-coexp","MICRF-PrePPI","MICRF-Co-PrePPI")
#     cols <- c("red","green","blue","orange","darkorange","yellow","yellow4","pink","deeppink","cyan","darkcyan","purple","purple4","deepskyblue","darkgreen")
#     
#     pdf(file="plot/recurrent_rankes_tmpg.pdf",width=12,height=10)
#     par(mar=c(6,6,3,2))
#     plot(log10(Vrank),ReM[1,],xaxt="n", type="b",col=cols[1],ylim=c(0,max(ReM)+0.05),xlab="rank",ylab="fraction of recurrent risk genes",cex.lab=2,cex.main=2,lwd=2,cex.axis=1.8)
#     axis(1,at=log10(Vrank),labels=Vrank,font=1,cex.axis=1.5)
#     for(i in 2:15){
#         lines(log10(Vrank),ReM[i,], type="b",col=cols[i],xaxt="n",lwd=2)
#     }
#     legend("bottomright",col=cols,legend=legend,lty=rep(1,15),lwd=rep(2,15),cex=1,y.intersp=1)
#     dev.off()
    
    
    Vrank <- c(50,100,200,500,1000,5000)
    ### plot figures
    for(i in 1:length(Vrank)){
        ReM <- Rlist[[i]]
        pdfname <- paste("plot/recurr_",Vrank[i],".pdf",sep="")
        if(flag > 1) pdfname <- paste("plot/recurr_",Vrank[i],"_",flag,"_",TPcut,".pdf",sep="")
        pdf(file=pdfname,width=12,height=7)
        subs <- 
        par(mai=c(2.3,1,1,1))
        cols <- c("red","green","blue","orange","yellow","pink","cyan")
        legend=c("TADA","DAWN","MAGI","MICRF-iRef","MICRF-HPRD","MICRF-corr","MICRF-coexp")
        boxplot(ReM[,subs],col=cols,ylab="Fraction of recurrent genes",cex.axis=1.5,cex.lab=1.7,cex.names=0,outline=FALSE, axisnames = FALSE,xaxt="n")
        axis(1,at=1:length(cols),labels=FALSE)
        text(1:length(cols), par("usr")[3], labels = legend[subs], srt = 45, adj = c(1.1,1.1), xpd = TRUE,cex=1.7)
        abline(h=median(ReM[,subs[1]]),lty=2)
        dev.off()
    }
    
}


randset4_one <- function(tmpg,oneresult){
      
    inrtmp <- rep(0, length(tmpg))
    subs <- tmpg %in% oneresult[,1]
    genes <- tmpg[subs]
    inrtmp[subs] <- rowSums(oneresult[match(genes,oneresult[,1]),c("dn.LoF","dn.mis3")])>0
    tmp <- sapply(1:length(tmpg), function(i) sum(inrtmp[1:i])) 
    
    tmp
}


## ROC plots
plot_AUC <- function(){
    source("misc_output.R")
    AUC <- dealAUC()
    
    aucf="plot/samplesize11_11_22.pdf"
    subs <- 1:7
    nmeth=length(subs)
    cols <- c("red","green","blue","orange","yellow","pink","cyan")
    legend <- c("TADA","DAWN","MAGI","MICRF-iRef","MICRF-HPRD","MICRF-corr","MICRF-coexp")
    
    ## above to change
    main="AUC performance for different trios"
    xlab="The number of trios";ylab="AUC";
    x <- seq(0.2,0.9,0.1)
    Y1 <- rep(0,nmeth)
    pdf(file=aucf,width=12,height=10)
    k=1
    for(i in subs){
        y <- rowMeans(AUC[,,i])
        sd <- apply( AUC[,,i], 1, sd)/2
        
        x <- x[1:5]
        y <- y[1:5]
        if(i==1){
            par(mar=c(6,6,3,2))
            plot(x, y,xaxt="n", type="b",col=cols[i], main=main, ylab=ylab,xlab=xlab,xlim=c(0.16,0.12+length(x)/10),ylim=c(0,1),cex.lab=2,cex.main=2,lwd=2,cex.axis=1.8)
            axis(1,at=seq(0.2,0.9,0.1)[1:length(x)],labels=floor(5542*seq(0.2,0.9,0.1))[1:length(x)],font=1,cex.axis=1.5)
        }else{
            lines(x, y, type="b",col=cols[i],xaxt="n",lwd=2)
        }
        Y1[k] <- y[1]
        k <- k+1
    } 
    Y2 <- round(Y1*100)/100
    labels <- 1:length(subs)
    tx <- rep(0.19,length(subs))
    ty <- -sort(-Y2)
    cols <- cols[subs]
    ords <- sort(-Y1,index.return=T)$ix
    text(tx,ty,labels=labels,pos=2,col=cols[ords],cex=1.3)
    legend("bottomright",col=cols[ords],legend=paste(labels,legend[ords],sep=":"),lty=rep(1,length(subs)),lwd=rep(2,length(subs)),cex=1.5,y.intersp=1)
    dev.off()

}

getAUC <- function(j,TPcut){
    
    source("misc_output.R")
    TADAFile="../TADA_DAWN/result/TADAdenovo_meta_dmis.csv"
    Tset <- allPf(TPcut,TADAFile)
    
    DAWNnames <-  readLines(con <- file("../TADA_DAWN/DAWN_package/TADAdenovo_randset_1.txt","r"))   
    close(con)
    Maginames <-  readLines(con <- file("../MAGI_V0.1_Alpha/mydata/randset_1/MAGIfilenames.txt","r"))   
    close(con)
    
    nnet <- 4
    AUC0 <- array(0,dim=c(20,3+nnet))
    for(i in 1:20){   
        filenames <- onefnew(j,i,DAWNnames,Maginames) 
        tmp <- new_auc_11_9(filenames,Tset)
        AUC0[i,] <- tmp$auc
    }
    
    AUC0
}

dealAUC <- function(){
    
    AUC <- array(0,dim=c(8,20,7)) 
    for(j in 2:9){
        AUC0 <- getAUC(j,TPcut=0.1)
        AUC[j-1,,] <- AUC0
    }
    
    AUC
}

oneROC <- function(){
    
    TPcut <- 0.1; j <- 2; i <- 10; fprc <- 1
    subs <- 1:7
    cols <- c("red","green","blue","orange","yellow","pink","cyan")
    legend <- c("TADA","DAWN","MAGI","MICRF-iRef","MICRF-HPRD","MICRF-corr","MICRF-coexp")
    
    #### above to change 
    source("misc_output.R")
    TADAFile="../TADA_DAWN/result/TADAdenovo_meta_dmis.csv"
    Tset <- allPf(TPcut,TADAFile)
    DAWNnames <-  readLines(con <- file("../TADA_DAWN/DAWN_package/TADAdenovo_randset_1.txt","r"))   
    close(con)
    Maginames <-  readLines(con <- file("../MAGI_V0.1_Alpha/mydata/randset_1/MAGIfilenames.txt","r"))   
    close(con)
    
    filenames <- onefnew(j,i,DAWNnames,Maginames) 
    tmp <- new_auc_11_9(filenames,Tset)
    TPR <- tmp$TPR
    FPR <- tmp$FPR
    
    pdf(file=paste("plot/one_ROC",fprc,".pdf",sep=""),width=12,height=10)
    par(mar=c(6,6,3,2))
    x <- FPR[[subs[1]]]
    y <- TPR[[subs[1]]]
    plot(x[x<=fprc],y[x<=fprc],type="l",col=cols[1],ylim=c(0,1),xlab="1-specificity",ylab="sensitivity",cex.lab=2,cex.main=2,lwd=2,cex.axis=1.8)
    for(i in 2:length(subs)){
        x <- FPR[[subs[i]]]
        y <- TPR[[subs[i]]]
        lines(x[x<=fprc],y[x<=fprc],type="l",col=cols[i],xaxt="n",lwd=2)
    }
    legend("bottomright",col=cols,legend=legend,lty=rep(1,15),lwd=rep(2,15),cex=1,y.intersp=1)
    lines(c(0,1),c(0,1),type="l",lwd=2,lty=2)
    dev.off()
    
}
