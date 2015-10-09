## step 1: fig1: AUC performance; compute AUC values
plot_AUC <- function(){
    ## plot auc values
    ## Figure 1: AUC values for eight random sets
    
    aucfile="AUC_10_8_2"
    aucf="plot/samplesize10_8.pdf"
    source("enrichana_6_12.R")
    load(aucfile)
    
    pdf(file=aucf,width=12,height=10)
    #cols <- c("black","red","green","aquamarine3","bisque4","blueviolet","brown","chartreuse","chocolate4","cyan","coral4","aquamarine3","bisque4","blueviolet","brown","chartreuse","chocolate4","cyan","coral4")
    cols <- c("black","red","green","mediumorchid1","orchid","blueviolet","brown","chartreuse","chocolate4","cyan","coral4","orangered","blue","blue","darkslategray1","","chocolate4","coral4","")
    #legend <- c("TADA","DAWN","MAGI","RWR-STRING","RWR-iRef","RWR-coexp","RWR-HPRD","RWR-corr","RWR-corr1","RWR-coexp5","RWR-coexp7","MICRF-STRING","MICRF-iRef","MICRF-coexp","MICRF-HPRD","MICRF-corr","MICRF-corr1","MICRF-coexp5","MICRF-coexp7")
    
    #subs <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19)
    cols <- cols[c(1,2,3,4,5,7,9,10,12,13,15,17,18)]
    subs <- c(1,2,3,4,5,7,9,10,12,13,15,17,18)
    nn <- subs
    nmeth=length(subs)
    legend <- c("TADA","DAWN","MAGI","RWR-STRING","RWR-iRef","RWR-HPRD","RWR-corr","RWR-coexp","MICRF-STRING","MICRF-iRef","MICRF-HPRD","MICRF-corr","MICRF-coexp")
    
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
    #ty[5] = ty[6] + 0.001
    #ty[6] = ty[6] + 0.001
    #ty[3] = ty[3] - 0.01
#     ty[1] = ty[1] + 0.01
#     ty[4] = ty[4] + 0.01
#     ty[6] = ty[6] - 0.01
#     ty[11] = ty[11] - 0.01
#     ty[12] = ty[12] - 0.01
    
    cols <- cols[subs]
    ords <- sort(-Y1,index.return=T)$ix
    text(tx,ty,labels=labels,pos=2,col=cols[ords],cex=1.3)
    legend("bottomright",col=cols[ords],legend=paste(labels,legend[ords],sep=":"),lty=rep(1,length(subs)),lwd=rep(2,length(subs)),cex=1.5,y.intersp=1)
    
    dev.off()
    

}

getAUC <- function(j){
    DAWNflag=2
    
    source("enrichana_6_12.R")
    TADAFile="../TADA_DAWN/result/TADAdenovo_meta_dmis.csv"
    Tset <- allPf(TPcut=0.1,TADAFile)
    
    DAWNnames <-  readLines(con <- file("../TADA_DAWN/DAWN_package/TADAdenovo_randset_1.txt","r"))   
    close(con)
    Maginames <-  readLines(con <- file("../MAGI_V0.1_Alpha/mydata/randset_1/MAGIfilenames.txt","r"))   
    close(con)
    
    nnet <- 8
    AUC0 <- array(0,dim=c(20,3+nnet*2))
    for(i in 1:20){   
        filenames <- onefnew(j,i,DAWNnames,Maginames)
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
    
    save(AUC0,file=paste("AUC_10_8",j,DAWNflag,sep="_"))
        
}

dealAUC <- function(){

    ### overlapped with DAWN predicted results
    DAWNflag=2
    AUC <- array(0,dim=c(8,20,19)) 
    for(j in 2:9){
        #load(paste("AUC_10_7_",j,sep="_"))
        load(paste("AUC_10_8",j,DAWNflag,sep="_"))
        AUC[j-1,,] <- AUC0
    }
    
    save(AUC,file=paste("AUC_10_8",DAWNflag,sep="_"))

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
    
    nnet <- length(netstr)
    hotnetfile <- 1:nnet
    misfile <- 1:nnet
    for(k in 1:nnet){
        hotnetfile[k] <- paste(resultstr,netstr[k],"hotnetresult1",tmpstr,netnum[k],".txt",sep="")
        if(k==1){
            misfile[k] <- paste("result/randresult_5/",netstr[k],"CRFresult_",3,tmpstr,"LBP_",netnum[k],".txt",sep="")
        }else{
            misfile[k] <- paste("result/randresult_5/",netstr[k],"CRFresult_",netnum[k],tmpstr,"LBP_",netnum[k],".txt",sep="")
        }
    }
    
    filenames <- c(TADAfile,DAWNfile,MAGIfile,hotnetfile,misfile)
    
    filenames
    
}

## step 2: fig2: recurrent performance 
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
    ReM <- matrix(0,20,13)
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
            if(netk==1){
                myfile <- paste("result/randresult4_5/",netstr[netk],"CRFresult_",3,"part3_",i,"LBP_",netnum[netk],".txt",sep="")
            }else{
                myfile <- paste("result/randresult4_5/",netstr[netk],"CRFresult_",netnum[netk],"part3_",i,"LBP_",netnum[netk],".txt",sep="")
            }
            cutoff <- total_genes_cutoff(myfile,Tset,muts="",alpha=0.5)
            tmp <- read.table(myfile)
            tmpg <- tmp[tmp[,2] >= cutoff,1]
            ReM[i,netk+8] <- length(intersect(tmpg,reTset))/length(reTset)##sum(rowSums(oneresult[match(tmpg,oneresult[,1]),c("dn.LoF","dn.mis3")])>0)/length(tmpg)
        }
    }
    
    save(ReM,file="ReM_10_8")
    
    ### plot figures ===================================
    load("ReM_10_8")
    pdf(file="plot/recurr_10_7.pdf",width=10,height=7)
    #subs <- c(1,2,3,4,9,7,12,5,10,8,13,6,11)
    subs <- c(1,2,3,4,9,5,10,6,11,7,12,8,13)
    #subs <- c(1,2,3,4,9,5,10,6,11,8,13)
    par(mai=c(2.3,1,1,1))
    cols <- c("red","green","dodgerblue","orange","orange","yellow","yellow","pink","pink","cyan","cyan","purple","purple")
    legend=c("TADA","DAWN","MAGI","RWR-STRING","RWR-iRef","RWR-HPRD","RWR-corr","RWR-coexp","MICRF-STRING","MICRF-iRef","MICRF-HPRD","MICRF-corr","MICRF-coexp")
    boxplot(ReM[,subs],col=cols,names=legend[subs],las=2,ylab="Fraction of recurrent genes",cex.axis=1.5,cex.lab=1.7,cex.names=2,outline=FALSE)
    abline(h=median(ReM[,subs[1]]),lty=2)
    dev.off()
    
}

## step 3: fig3: leaveon performance
leaveone_recurrent <- function(){
    source("enrichana_6_12.R")
    options(stringsAsFactors=FALSE)
    TPcut=0.1
    netstr <- c("STRING/","iRef/","HPRD/","corr/","coexp5/")
    netnum <- c(31,6,20,26,27)
    TADAFile <- "../TADA_DAWN/result/TADAdenovo_meta_dmis.csv"
    asdall <- read.csv(TADAFile)
    asdall[is.na(asdall[,"qvalue.dn"]),"qvalue.dn"] <- 1
    Tset <- asdall[asdall[,"qvalue.dn"]<=TPcut,1]
    
    ### leave one mutation ==========================================
    nLeo <- matrix(0,13,2)
    
    oneg <- c()
    ### TADA method
    cfile <- "../TADA_DAWN/result/TADAdenovo_control1911.csv"
    ctmp <- read.csv(cfile)
    for(i in 1:length(Tset)){
        filename <- paste("../TADA_DAWN/result/leaveone4/TADAdenovo_rand2_",i,".csv",sep="")
        cutoff <- total_genes_cutoff(filename,Tset,muts="",alpha=0.5)
        tmp <- read.csv(filename)
        tmp[is.na(tmp[,"qvalue.dn"]),"qvalue.dn"] <- 1
        oner <- tmp[match(Tset[i],tmp[,1]),"qvalue.dn"]  <= 1-cutoff
        nLeo[1,1] <- nLeo[1,1] + oner
        #if(oner & (ctmp[match(Tset[i],ctmp[,1]),2] <= 1-cutoff)) nLeo[1,2] <- nLeo[1,2] + 1
        if( ctmp[match(Tset[i],ctmp[,1]),"qvalue.dn"] <= 1-cutoff ) nLeo[1,2] <- nLeo[1,2] + 1
    }
    
    ## DAWN method 
    ### control data result  
    cfile <- "../TADA_DAWN/DAWN_package/DAWNcontrol1911.csv" ###  checked
    ctmp <- read.csv(cfile)
    ctmp[is.na(ctmp[,"Stage2_posterior"]),"Stage2_posterior"] <- 1
    
    for(i in 1:length(Tset)){
        filename <- paste("../TADA_DAWN/DAWN_package/result/leaveone4/DAWN_leaveone4_",i,".csv",sep="")
        cutoff <- total_genes_cutoff(filename,Tset,muts="",alpha=0.5)
        tmp <- read.csv(filename)
        tmp[is.na(tmp[,"Stage2_posterior"]),"Stage2_posterior"] <- 1
        oner <- tmp[match(Tset[i],tmp[,1]),"Stage2_posterior"]  <= 1-cutoff
        
        if(!is.na(oner)){
            nLeo[2,1] <- nLeo[2,1] + oner
            if(oner) oneg <- union(oneg,Tset[i])
        }
        
        if( Tset[i] %in% ctmp[,1] ){
            if(ctmp[match(Tset[i],ctmp[,1]),"Stage2_posterior"] <= 1-cutoff) nLeo[2,2] <- nLeo[2,2] + 1
        }
        
    }
    DAWNg <- oneg
    
    oneg <- c()
    ###  MAGI method
    cfile <- "../MAGI_V0.1_Alpha/mydata/control/RandomGeneList.3"
    ctmp <- read.table(cfile)
    
    for(i in 1:length(Tset)){
        filename <- paste("../MAGI_V0.1_Alpha/mydata/leaveone4/RandomGeneList.",i,sep="")
        cutoff <- total_genes_cutoff(filename,Tset,muts="",alpha=0.5)
        tmp <- read.table(filename)
        oner <- tmp[match(Tset[i],tmp[,1]),2]  >= cutoff
        
        if(!is.na(oner)){
            nLeo[3,1] <- nLeo[3,1] + oner
            if(oner) oneg <- union(oneg,Tset[i])
        }
        if(Tset[i] %in% ctmp[,1]){
            if( ctmp[match(Tset[i],ctmp[,1]),2] >= cutoff ) nLeo[3,2] <- nLeo[3,2] + 1;}
    }
    MAGIg <- oneg
    
    oneg <- c()
    ### hotnet2 method 
    for(j in 1:5){
        cfile <- paste("result/control/hotnet/",netstr[j],"/hotnetresult1control1911",netnum[j],".txt",sep="")
        ctmp <- read.table(cfile)
        
        for(i in 1:length(Tset)){
            filename <- paste("result/leaveone4result/",netstr[j],"hotnetresult1rand2_",i,netnum[j],".txt",sep="")
            cutoff <- total_genes_cutoff(filename,Tset,muts="",alpha=0.5)
            tmp <- read.table(filename)
            oner <- tmp[match(Tset[i],tmp[,1]),2]  >= cutoff
            
            if(!is.na(oner)){
                if(oner) oneg <- union(oneg,Tset[i])
                nLeo[j+3,1] <- nLeo[j+3,1] + oner
            }
            if(Tset[i] %in% ctmp[,1]){
                if( ctmp[match(Tset[i],ctmp[,1]),2] >= cutoff ) nLeo[j+3,2] <- nLeo[j+3,2] + 1;}
        }       
    }
    RWRg <- oneg
    
    oneg <- c()
    usoneg <- c()
    ### our method
    for(j in 1:5){
        if(j==1){
            cfile <- paste("result/control/v5/",netstr[j],"CRFresult_",3,"control1911_2LBP_",netnum[j],".txt",sep="")
            }else{
            cfile <- paste("result/control/v5/",netstr[j],"CRFresult_",netnum[j],"control1911_2LBP_",netnum[j],".txt",sep="")
        }
        ctmp <- read.table(cfile)
        for(i in 1:length(Tset)){
            if(j==1){
                filename <- paste("result/leaveone4result_5/",netstr[j],"CRFresult_",3,"rand2_",i,"LBP_",netnum[j],".txt",sep="")
            }else{
                filename <- paste("result/leaveone4result_5/",netstr[j],"CRFresult_",netnum[j],"rand2_",i,"LBP_",netnum[j],".txt",sep="")
            }
            
            cutoff <- total_genes_cutoff(filename,Tset,muts="",alpha=0.5)
            tmp <- read.table(filename)
            oner <- tmp[match(Tset[i],tmp[,1]),2]  >= cutoff
            if(!is.na(oner)){
                if(oner) {oneg <- union(oneg,Tset[i]); usoneg <- union(usoneg,Tset[i]);}
                nLeo[j+8,1] <- nLeo[j+8,1] + oner
            }
            if(Tset[i] %in% ctmp[,1]){
                if( ctmp[match(Tset[i],ctmp[,1]),2] >= cutoff ) nLeo[j+8,2] <- nLeo[j+8,2] + 1;}
        }       
    }
    MICRFg <- oneg
    save(nLeo,file="nLeo_10_7")
    
    ### plot figures
    load("nLeo_10_7")
    pdf(file="plot/leaveone_10_7.pdf",width=10,height=10)
    par(mai=c(4,2,2,2))
    legend <- c("TADA","DAWN","MAGI","RWR-STRING","RWR-iRef","RWR-HPRD","RWR-corr","RWR-coexp","MICRF-STRING","MICRF-iRef","MICRF-HPRD","MICRF-corr","MICRF-coexp")
    ords <- c(1,2,3,4,9,5,10,6,11,7,12,8,13)
    #ords <- c(1,2,3,4,9,5,10,6,11,8,13)
    #cols <- c("red","green","blue","orange","orange","yellow","yellow","pink","pink","cyan","cyan","purple","purple")
    
    barplot(t(nLeo)[,ords], names.arg=legend[ords], beside=T, ylab="#gold risk genes",cex.names=1.5, las=2, ylim=c(0,max(nLeo)), col=c("red","darkblue"),cex.lab=1.5,cex.axis=1.25)
    legend("topleft",legend=c("DDD case","SSC control"),cex=1.3,fill=c("red","darkblue")) 
    dev.off()    
    
}
