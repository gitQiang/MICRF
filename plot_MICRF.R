leaveone_Rank <- function(Rcut,flag=1){
    
    options(stringsAsFactors=FALSE)
    TPcut=0.1
    TADAFile <- "../TADA_DAWN/result/TADAdenovo_meta_dmis.csv"
    asdall <- read.csv(TADAFile)
    asdall[is.na(asdall[,"qvalue.dn"]),"qvalue.dn"] <- 1
    Tset <- asdall[asdall[,"qvalue.dn"]<=TPcut,1]
    ### leave one mutation ==========================================
    if(flag==1) nLeo <- array(0,dim=c(8,2,1))
    if(flag==2) nLeo <- array(0,dim=c(8,2,Rcut)) 

    ### contorls data
    cfiles <- c("../TADA_DAWN/result/TADAdenovo_control1911.csv","../TADA_DAWN/DAWN_package/DAWNcontrol1911.csv","../MAGI_V0.1_Alpha/mydata/control/RandomGeneList.3",paste('/ifs/scratch/c2b2/ys_lab/qh2159/Mutations/CHD/MIS/result/control/v5/MICRFresult_',2:6,".txt",sep=""))
    for(i in 1:length(cfiles)){
        nLeo[i,2,] <- onefileC(Tset,cfiles[i],i,Rcut,flag)
    }
    
    n.t <- length(Tset)
    tmp <- array(0,dim=c(8,n.t))
    outputstr='/ifs/scratch/c2b2/ys_lab/qh2159/Mutations/CHD/MIS/result/control/v5/MICRFresult_';
    for(i in 1:n.t){
        ### TADA method 
        filename <- paste("../TADA_DAWN/result/leaveone4/TADAdenovo_rand2_",i,".csv",sep="")  
        tmp[1,i] <- onefileL(Tset[i],filename,1,Rcut,flag)
        ## DAWN method 
        filename <- paste("../TADA_DAWN/DAWN_package/result/leaveone4/DAWN_leaveone4_",i,".csv",sep="")
        tmp[2,i] <- onefileL(Tset[i],filename,2,Rcut,flag)
        ###  MAGI method
        filename <- paste("../MAGI_V0.1_Alpha/mydata/leaveone4/RandomGeneList.",i,sep="")
        tmp[3,i] <- onefileL(Tset[i],filename,3,Rcut,flag)
        ### MICRF method
        for(j in 1:5){
            filename <- paste("result/leaveone4result_5/MICRFresult_",j+1,"_",i,".txt",sep="");
            tmp[3+j,i] <- onefileL(Tset[i],filename,3+j,Rcut,flag)
        }
    }
    
    for(i in 1:8){
        if(flag==1)  nLeo[i,1,] <- sum(tmp[i,]) 
        if(flag==2){
            inrtmp <- rep(0,Rcut);
            #inrtmp[tmp[i,tmp[i,] > 0]] <- 1; ## occur in the same position
            ttmp <- table(tmp[i,])
            subs <- as.numeric(names(ttmp)) > 0
            inrtmp[as.numeric(names(ttmp))[subs]] <- ttmp[subs];
            nLeo[i,1,] <- sapply(1:Rcut, function(kk) sum(inrtmp[1:kk]))
        }
    }
    
    nLeo          
                  
}

leaveone_Rankes <- function(){
    source("plot_MICRF.R")
    Vrank <- c(100,200,300,500,1000,5000)
    TPM <- matrix(0,8,length(Vrank))
    FPM <- TPM
    for(i in 1:length(Vrank)){
        tmp <- leaveone_Rank(Vrank[i]);
        TPM[,i] <- tmp[,1,1]
        FPM[,i] <- tmp[,2,1]
    }
    leaveone_plots(TPM,FPM)
    
    source("plot_MICRF.R")
    Rcut=500
    tmp <- leaveone_Rank(Rcut,2);
    str <- c("case","control")
    for (j in 1:2){
    atmp <- tmp[,j,]
    pdf(file=paste("plot/leaveone",str[j],".pdf",sep=""),width=10,height=8)
    par(mai=c(2,1,1,1))
    main=""
    xlab="Rank"
    ylab="Number of risk genes"
    cols <- c("black","red","blue","brown","blueviolet","deeppink","darkgreen","darkcyan")
    legend=c("TADA","DAWN","MAGI","MICRF-iRef","MICRF-HPRD","MICRF-corr","MICRF-coexp","MICRF-CoPrePPI")
    plot_matrix(atmp,cols,main,xlab,ylab)
    legend("topleft",legend=legend,col=cols,lwd=2,lty=1,cex=1.25,y.intersp=1)
    dev.off()
    }
}

leaveone_plots <- function(TPM,FPM){
    Vrank <- c(100,200,300,500,1000,5000)
    ### plot figures
    for(i in 1:length(Vrank)){
        nLeo <- cbind(TPM[,i],FPM[,i])
        pdf(file=paste("plot/leaveone_",Vrank[i],".pdf",sep=""),width=10,height=10)
        par(mai=c(4,2,2,2))
        legend=c("TADA","DAWN","MAGI","MICRF-STRING","MICRF-iRef","MICRF-HPRD","MICRF-corr","MICRF-coexp")
        ords <- 1:dim(nLeo)[1]
        mp <- barplot(t(nLeo)[,ords], beside=T, ylab="#gold risk genes",cex.names=1.5, las=2, ylim=c(0,max(nLeo)), col=c("red","darkblue"),cex.lab=1.5,cex.axis=1.25,axes = FALSE, axisnames = FALSE)
        text(colMeans(mp), par("usr")[3], labels = legend[ords], srt = 45, adj = c(1.1,1.1), xpd = TRUE,cex=1.5)
        axis(2)
        legend("topleft",legend=c("DDD case","SSC control"),cex=1.3,fill=c("red","darkblue")) 
        dev.off()
    }
}

onefileL <- function(gene,filename,i,Rcut,flag){
    if(grepl(".csv",filename)){ tmp <- read.csv(filename); }else{ tmp <- read.table(filename);}
    if(i>2) tmp <- tmp[order(-as.numeric(tmp[,2])), ]
    oner <- gene %in% tmp[1:Rcut,1];
    if(flag==2){if(oner){oner=which(tmp[1:Rcut,1] %in% gene);}}
    
    oner
}

onefileC <- function(Tset,cfile,i,Rcut,flag){
    if(grepl(".csv",cfile)){ ctmp <- read.csv(cfile); }else{ ctmp <- read.table(cfile);}
    if(i>2) ctmp <- ctmp[order(-as.numeric(ctmp[,2])),]
    
    tmpg <- ctmp[1:Rcut,1]
    if(flag==1){
        oner <- sum(Tset %in% tmpg)
    }else if(flag==2){
        inrtmp <- rep(0,Rcut)
        inrtmp[tmpg %in% Tset] <- 1
        oner <- sapply(1:Rcut, function(kk) sum(inrtmp[1:kk])) 
    }
    
    oner
}

## try missing data
recurrent_Rank <- function(Rcut,flag=1,TPcut=0.1){
    ## Rcut: rank cutoff
    ## flag: only risk gene or not
    options(stringsAsFactors=FALSE)
    TADAFile <- "../TADA_DAWN/result/TADAdenovo_meta_dmis.csv"
    asdall <- read.csv(TADAFile)
    asdall[is.na(asdall[,"qvalue.dn"]),"qvalue.dn"] <- 1
    Tset <- asdall[asdall[,"qvalue.dn"]<=TPcut,1]
    
    #######  independent samples  ===================================
    if(flag <=2) ReM <- array(0,dim=c(20,8,1))
    if(flag ==3) ReM <- array(0,dim=c(20,8,Rcut))
    j=3
    for(i in 1:20){
        #onefile <- paste("DDD_mutations/randset4/ASD2sub",j,"derest",i,".csv",sep="")
        onefile <- paste("DDD_mutations/randset4/ASD1sub",j,"depart",i,".csv",sep="")
        oneresult <- read.csv(onefile)
        
        ### TADA method
        #TADAfile <- paste("../TADA_DAWN/result/randset4/TADAdenovo_ASD1sub",j,"depart",i,".csv",sep="")
        TADAfile <- paste("../TADA_DAWN/result/randset4/TADAdenovo_ASD2sub",j,"derest",i,".csv",sep="")
        ReM[i,1,] <- onefile(TADAfile,Rcut,Tset,oneresult,flag,0);
        
        ### DAWN method
        #DAWNfile <- paste("../TADA_DAWN/DAWN_package/result/randset4/DAWN_randset4_",i,".csv",sep="")
        DAWNfile <- paste("../TADA_DAWN/DAWN_package/result/randset4/DAWN_randset4rest_",i,".csv",sep="")
        ReM[i,2,] <- onefile(DAWNfile,Rcut,Tset,oneresult,flag,0);
        
        ### MAGI method 
        #MAGIfile <- paste("../MAGI_V0.1_Alpha/mydata/randset4/RandomGeneList.",i,sep="")
        MAGIfile <- paste("../MAGI_V0.1_Alpha/mydata/randset4rest/RandomGeneList.",i,sep="")
        ReM[i,3,] <- onefile(MAGIfile,Rcut,Tset,oneresult,flag,1);
        
        ### our method
        path = "/ifs/scratch/c2b2/ys_lab/qh2159/Mutations/CHD/MIS/result/randresult4_5/"
        for(netk in 1:5){
            myfile <- paste(path,"MICRFresult_",netk+1,"_",i,".txt",sep="")
            ReM[i,3+netk,] <- onefile(myfile,Rcut,Tset,oneresult,flag,1);
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
    setwd("/ifs/scratch/c2b2/ys_lab/qh2159/Mutations/CHD/MIS/")
    source("plot_MICRF.R")
    Rcut=500
    tmp <- recurrent_Rank(Rcut,3,0.1);
    atmp <- matrix(0,8,Rcut);
    for(i in 1:8) atmp[i,] <- colMeans(tmp[,i,])
    pdf(file="plot/reranks.pdf",width=10,height=8)
    par(mai=c(2,1,1,1))
    main=""
    xlab="Rank"
    ylab="Number of recurrent genes"
    cols <- c("black","red","blue","brown","blueviolet","deeppink","darkgreen","darkcyan")
    legend=c("TADA","DAWN","MAGI","MICRF-iRef","MICRF-HPRD","MICRF-corr","MICRF-coexp","MICRF-CoPrePPI")
    plot_matrix(atmp,cols,main,xlab,ylab)
    legend("bottomright",legend=legend,col=cols,lwd=2,lty=1,cex=1.25,y.intersp=1)
    dev.off()
}

recurrent_plots <- function(Rlist,flag=1,TPcut=0.3){    
    Vrank <- c(50,100,200,500,1000,5000)
    cols <- c("red","green","blue","brown","orange","yellow","pink","cyan")
    legend=c("TADA","DAWN","MAGI","MICRF-STRING","MICRF-iRef","MICRF-HPRD","MICRF-corr","MICRF-coexp")
    ### plot figures
    for(i in 1:length(Vrank)){
        ReM <- Rlist[[i]]
        pdfname <- paste("plot/recurr_",Vrank[i],".pdf",sep="")
        if(flag > 1) pdfname <- paste("plot/recurr_",Vrank[i],"_",flag,"_",TPcut,".pdf",sep="")
        pdf(file=pdfname,width=12,height=8)
        subs <- 1:dim(ReM)[2]
        par(mai=c(2.3,1,1,1))
        boxplot(ReM[,subs,1],col=cols,ylab="Fraction of recurrent genes",cex.axis=1.5,cex.lab=1.8,cex.names=0,outline=FALSE, axisnames = FALSE,xaxt="n")
        axis(1,at=1:length(cols),labels=FALSE)
        text(1:length(cols), par("usr")[3], labels = legend[subs], srt = 45, adj = c(1.1,1.1), xpd = TRUE,cex=1.8)
        abline(h=median(ReM[,subs[1],1]),lty=2)
        dev.off()
    }
    
}

onefile <- function(filename,Rcut,Tset,oneresult,flag,k){
    if(grepl(".csv",filename)){ tmp <- read.csv(filename); }else{ tmp <- read.table(filename);}
    if(k==1){ tmp <- tmp[order(-as.numeric(tmp[,2])), ];}
    tmpg <- tmp[1:Rcut,1]
    if(flag == 2) tmpg <- intersect(tmpg,Tset)
    if(flag <= 2){
        oner <- sum(rowSums(oneresult[match(tmpg[tmpg %in% oneresult[,1]],oneresult[,1]),c("dn.LoF","dn.mis3")])>0)
    }else if(flag==3){
        oner <- randset4_one(tmpg,oneresult)
    }
    oner
}

randset4_one <- function(tmpg,oneresult){
      
    inrtmp <- rep(0, length(tmpg))
    subs <- tmpg %in% oneresult[,1]
    genes <- tmpg[subs]
    inrtmp[subs] <- rowSums(oneresult[match(genes,oneresult[,1]),c("dn.LoF","dn.mis3")])>0
    tmp <- sapply(1:length(tmpg), function(i) sum(inrtmp[1:i])) 
    
    tmp
}

plot_matrix <- function(tmp,cols,main,xlab,ylab,byrow=TRUE){
    # plot each row
    if(!byrow) tmp <- t(tmp)
    n <- dim(tmp)[1]
    ymax <- max(tmp)
    plot(1:dim(tmp)[2],tmp[1,],type='l',col=cols[1],main=main,xlab=xlab,ylab=ylab,ylim=c(0,ymax),xaxs="i",yaxs="i",cex.lab=1.5,cex.axis=1.25,lwd=2)
    for(i in 2:n){
        lines(1:dim(tmp)[2],tmp[i,],type='l',col=cols[i],lwd=2)
    }
}

## ROC plots
plot_AUC <- function(){
    source("plot_MICRF.R")
    source("misc_output.R")
    AUC <- array(0,dim=c(8,20,8)) 
    TPcut <- 0.1
    for(j in 2:9){
        AUC0 <- getAUC(j,TPcut)
        AUC[j-1,,] <- AUC0
    }
    
    aucf=paste("plot/AUC_11_26_",TPcut,".pdf",sep="")
    subs <- 1:8
    nmeth=length(subs)
    cols <- c("black","red","blue","brown","blueviolet","deeppink","darkgreen","darkcyan")
    legend=c("TADA","DAWN","MAGI","MICRF-iRef","MICRF-HPRD","MICRF-corr","MICRF-coexp","MICRF-CoPrePPI")
    
    ## above to change
    pdf(file=aucf,width=12,height=10)
    par(mar=c(6,6,3,2))
    main="AUC performance"
    xlab="Number of trios";ylab="AUC";
    x <- seq(0.2,0.9,0.1)
    Y1 <- rep(0,nmeth)
    for(i in subs){
        y <- rowMeans(AUC[,,i])
        if(i==1){
            plot(x, y,xaxt="n", type="b",col=cols[i], main=main, ylab=ylab,xlab=xlab,xlim=c(0.16,0.12+length(x)/10),ylim=c(0,1),cex.lab=2,cex.main=2,lwd=2,cex.axis=1.8)
            axis(1,at=seq(0.2,0.9,0.1)[1:length(x)],labels=floor(5542*seq(0.2,0.9,0.1))[1:length(x)],font=1,cex.axis=1.5)
        }else{
            lines(x, y, type="b",col=cols[i],xaxt="n",lwd=2)
        } 
    } 
    legend("bottomright",col=cols,legend=legend,lty=rep(1,nmeth),lwd=rep(2,nmeth),cex=1.5,y.intersp=1)
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
    
    nnet <- 5
    AUC0 <- array(0,dim=c(20,3+nnet))
    for(i in 1:20){   
        filenames <- onefnew(j,i,DAWNnames,Maginames) 
        tmp <- new_auc_11_9(filenames,Tset)
        AUC0[i,] <- tmp$auc
    }
    
    AUC0
}

oneROC <- function(){
    
    TPcut <- 0.1; j <- 3; i <- 10; fprc <- 0.2
    subs <- 1:8
    cols <- c("black","red","blue","brown","blueviolet","deeppink","darkgreen","darkcyan")
    legend=c("TADA","DAWN","MAGI","MICRF-iRef","MICRF-HPRD","MICRF-corr","MICRF-coexp","MICRF-CoPrePPI")
    nmeth <- length(subs)
    
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
    for(i in 2:nmeth){
        x <- FPR[[subs[i]]]
        y <- TPR[[subs[i]]]
        lines(x[x<=fprc],y[x<=fprc],type="l",col=cols[i],xaxt="n",lwd=2)
    }
    legend("bottomright",col=cols,legend=legend,lty=rep(1,nmeth),lwd=rep(2,nmeth),cex=1,y.intersp=1)
    lines(c(0,1),c(0,1),type="l",lwd=2,lty=2)
    dev.off()
    
}
