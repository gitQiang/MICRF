# check randset1 not done
checkrandset1 <- function(){
        nSim=160;
        inputpath='/ifs/scratch/c2b2/ys_lab/qh2159/Mutations/CHD/MIS/result/randset_1/';
        inputstr='rand1';
        outputstr='/ifs/scratch/c2b2/ys_lab/qh2159/Mutations/CHD/MIS/result/randresult_5/MICRFresult_';
        filenames = 1:640
        for(kk in 161:800){
                netj=kk %% nSim;
                if(netj == 0)netj=nSim;
                
                netj_i = floor((netj-1)/20)+2;
                netj_j = netj %% 20;
                if (netj_j == 0) netj_j=20;
                
                netflag=floor((kk-1)/nSim)+1;
                if (netflag == 1){
                        instr = 'CRF_input';
                        idm=1;
                }else{
                        instr = 'hotnet_input'; 
                        idm=0;
                }
                
                outputfile=paste(outputstr,netflag,'_',netj_i,'_',netj_j,'.txt',sep="");
                filenames[kk-160] <- basename(outputfile)
        }
        
        runfile <- list.files("/ifs/scratch/c2b2/ys_lab/qh2159/Mutations/CHD/MIS/result/randresult_5/","*.txt$")
        
        subs <- which(!(filenames %in% runfile)) + 160
        qwt(subs,file="randset1Rest.txt")
}

## check the additional genes finding by coexp and PrePPI networks
addig <- function(kk=4){
        source("misc_output.R")
        nmeth=9
        j=2;TPcut=0.1;
        subs <- 1:nmeth
        gf=TRUE
        ncut=100
        
        TADAFile="../TADA_DAWN/result/TADAdenovo_meta_dmis.csv"
        Tset <- allPf(TPcut,TADAFile)
        DAWNnames <-  readLines(con <- file("../TADA_DAWN/DAWN_package/TADAdenovo_randset_1.txt","r"))   
        close(con)
        Maginames <-  readLines(con <- file("../MAGI_V0.1_Alpha/mydata/randset_1/MAGIfilenames.txt","r"))   
        close(con)
        
        filenames <- onefnew(j,kk,DAWNnames,Maginames) 
        if(gf){ tmp <- read.csv(filenames[1])
        genes <- tmp[rowSums(tmp[,c("dn.LoF","dn.mis3")]) > 0,1]
        }  
        
        glist <- list()
        nvec <- rep(0,length(filenames))
        aa <- c()
        for(i in 1:length(filenames)){
                filename <- filenames[i]
                if(grepl(".csv",filename)){ result <- read.csv(filename); 
                }else{result <- read.table(filename);}
                result <- cbind(result,1:dim(result)[1])
                
                if(gf) result <- result[result[,1] %in% genes,]
                
                if(i==3) result <- result[order(-result[,2]),]
                
                glist[[i]] <- intersect(result[1:ncut,1],Tset)
                nvec[i] <- length(glist[[i]])
                if(i<9) aa <- c(aa,glist[[i]])
        }
        
        ones <- setdiff(glist[[9]],glist[[8]])
        uni1 <- setdiff(glist[[9]],aa)
        
        print(kk)
        #print(ones)
        print(uni1)
        
}

## one gene example to show the CoPrePPI effect
oneGeneEx <- function(){
        source("misc_output.R")
        j=2;TPcut=0.1;i=16;
        
        TADAFile="../TADA_DAWN/result/TADAdenovo_meta_dmis.csv"
        Tset <- allPf(TPcut,TADAFile)
        DAWNnames <-  readLines(con <- file("../TADA_DAWN/DAWN_package/TADAdenovo_randset_1.txt","r"))   
        close(con)
        Maginames <-  readLines(con <- file("../MAGI_V0.1_Alpha/mydata/randset_1/MAGIfilenames.txt","r"))   
        close(con)
        
        filenames <- onefnew(j,i,DAWNnames,Maginames) 
        ## gene edges in Coexp and PrePPI
        gene="SH3D19" # PTEN
        #coexp <- read.table("data/network_inference/brainspan_net_top5.txt")
        #subs <- coexp[,1]==gene | coexp[,2]==gene
        #subco <- coexp[subs,]
        #subco[,3] <- 1
        preppi <- read.table("data/PrePPI.txt")
        subs <- preppi[,1]==gene | preppi[,2]==gene
        subpre <- preppi[subs,]
        subpre[,3] <- 2
        #subC <- rbind(subco,subpre)
        subC <- subpre #!!!!!
        #qwt(subC,file=paste("CoPrePPI",gene,"net.txt",sep=""))
        
        genes <- union(subC[,1],subC[,2])
        a <- read.csv(filenames[1])
        s <- a[match(genes,a[,1]),"BF.dn"]
        s[is.na(s)] <- 1
        pi=0.06
        s <- (s*pi)/(s*pi + 1-pi)
        qwt(paste(genes,"=",s,sep=" "),file=paste("nodeS",gene,"_s.txt",sep=""))
        
        preppi <- read.table("data/PrePPI.txt")
        subs <- preppi[,1] %in% genes & preppi[,2] %in% genes
        subpre <- preppi[subs,]
        subpre[,3] <- 2
        subC <- subpre
        qwt(subC,file=paste("CoPrePPI",gene,"edges.txt",sep=""))
        
        lab <- rep(0,length(genes))
        lab[genes %in% Tset] <- 1
        rgene <- paste(genes," = ",lab,sep="")
        names(rgene) <- "Adjrisk"
        qwt(rgene,file=paste("AdjRisk",gene,"g.txt",sep=""),flag=2)
        
        ####================
        ## local change node color
        aa <- read.table(paste("CoPrePPI",gene,"net.txt",sep=""))
        genes <- unique(c(aa[,1],aa[,2]))
        col <- rep(0,length(genes))
        a1 <- union(aa[aa[,3]==1,1],aa[aa[,3]==1,2])
        a2 <- union(aa[aa[,3]==2,1],aa[aa[,3]==2,2])
        a3 <- intersect(a1,a2)
        col[genes %in% a1] <- 1
        col[genes %in% a2] <- 2
        col[genes %in% a3] <- 3
        col[genes %in% gene] <- 0
        genes <- paste(genes," = ", col, sep="")
        qwt(genes,file=paste(gene,"geneCols.txt",sep=""),flag=2)
        
        dd <- read.csv("DDD_mutations/datasheet/anno_DDDmutationlist_4_15_all.csv")
        ndd <- read.csv("DDD_mutations/datasheet/anno_mutationlist_4_15_all.csv")
        ndd <- rbind(dd,ndd)
        ndd <- ndd[ndd[,"mutation_type_metaSVM"] %in% c("LOF","dMIS"),]
        aa <- ndd[ndd[,"Gene"]==gene,]
        
        
}

## check the additional risk genes identify by our methods
MICRFaddg <- function(){
        source("enrichana_6_12.R")
        source("addresult.R")
        source("~/.Rprofile")
        ncut <- 246 # TADA FDR 0.3
        geneT <- geneList9(ncut)
        
        #load("GeneT")
        aa <- c()
        for(i in 1:3){
                aa <- c(aa,geneT[[i]])
        }
        aa1 <- names(table(aa))[table(aa) >= 2]
        aa <- unique(aa)
        
        bb <- geneT[[4]]
        for(i in 5:9){
                bb <- c(bb,geneT[[i]])
        }
        bb <- names(table(bb))[table(bb)>=4]
        newg <- setdiff(bb,aa)
        allgenes <- union(aa1,bb)
        
        nodecol <- cbind(allgenes,1)
        nodecol[nodecol[,1] %in% newg,2] <- 2 
        ## annota new genes
        dd <- read.csv("DDD_mutations/datasheet/anno_DDDmutationlist_4_15_all.csv")
        ndd <- read.csv("DDD_mutations/datasheet/anno_mutationlist_4_15_all.csv")
        ndd <- rbind(dd,ndd)
        ndd <- ndd[ndd[,"mutation_type_metaSVM"] %in% c("LOF","dMIS"),]
        
        source("MICRF_net.R")
        filename <- "data/network_inference/ComCo_PrePPI.txt";
        net.text <- as.matrix(read.table(filename,sep="\t",header=FALSE))
        net1 <- build_net(36) ##CoPrePPI
        allgenes <- allgenes[allgenes %in% net1$node]
        
        deg <- colSums(net1$matrix>0)
        nodedeg <- cbind(allgenes,deg[match(allgenes,net1$node)])
        
        allgenes0=allgenes
        relg <- net1$node[colSums(net1$matrix[newg,]>0) >= 2]
        relg1 <- net1$node[colSums(net1$matrix[setdiff(allgenes,newg),]>0) >= 2]
        relg <- intersect(relg,relg1)
        allgenes <- union(allgenes,relg)
        subnet <- net1$matrix[allgenes,allgenes]
        subnet[upper.tri(subnet)] <- 0
        sube <- which(subnet > 0, arr.ind=TRUE)
        edges <- cbind(allgenes[sube[,1]],allgenes[sube[,2]])
        
        disV <- c()
        for(i in 1:length(allgenes)){
                disV[i] <- paste(unique(ndd[ndd[,"Gene"]==allgenes[i],"Disease"]),collapse=",")
        }
        nodedis <- cbind(allgenes,paste(allgenes,":",disV,sep=""))
        
        nodesize <- rep(0,length(allgenes))
        for(i in 1:length(allgenes)){
                nodesize[i] <- sum(ndd[,"Gene"]==allgenes[i])
        }
        nodesize[is.na(nodesize)] <- 0
        nodesize <- cbind(allgenes,nodesize)
        
        qwt(paste(nodecol[,1]," = ",nodecol[,2],sep=""),file="NDDnodeCol.txt",flag=2)
        qwt(paste(nodedeg[,1]," = ",nodedeg[,2],sep=""),file="NDDnodedeg.txt",flag=2)
        qwt(paste(nodedis[,1]," = ",nodedis[,2],sep=""),file="NDDnodedisease.txt",flag=2)
        nodedis[!(nodedis[,1] %in% newg),2] <- nodedis[!(nodedis[,1] %in% newg),1]
        qwt(paste(nodedis[,1]," = ",nodedis[,2],sep=""),file="NDDnodedisease_newg.txt",flag=2)
        qwt(edges,file="NDDsubnet.txt")
        qwt(paste(nodesize[,1]," = ",nodesize[,2],sep=""),file="NDDnodesize.txt",flag=2)
        qwt(newg,file="newg.txt")
        qwt(setdiff(allgenes0,newg),file="knowng.txt")
        
        #     library(VennDiagram)
        #     #venn.diagram(x=list(TADA=geneT[[1]],DAWN=geneT[[2]],MAGI=geneT[[3]],MICRF_STRING=geneT[[4]],MICRF_iRef=geneT[[5]],MICRF_HPRD=geneT[[6]],MICRF_corr=geneT[[7]],MICRF_coexp=geneT[[8]],MICRF_CoPrePPI=geneT[[9]]),file="test.png",margin=c(0.2,0.2,0.2,0.2))
        #     venn.diagram(x=list(TADA=geneT[[1]],DAWN=geneT[[2]],MAGI=geneT[[3]],MICRF=geneT[[9]]),file="test.png",margin=c(0.2,0.2,0.2,0.2))
        #     venn.diagram(x=list(iRefIndex=geneT[[5]],HPRD=geneT[[6]],COEXP=geneT[[8]],CoPrePPI=geneT[[9]]),file="test.png",margin=c(0.2,0.2,0.2,0.2))
        
        
        #     tadar <- read.csv(filenames[1])
        #     ourr <- read.table(filenames[9])
        #     
        #     tag <- tadar[tadar[,"qvalue.dn"] >= 0.3,1]
        #     tmp <- intersect(ourr[1:200,1],tag)
        #     qwt(tmp,file="MICRF_top200_TADA0.3.txt")
        
        # Ltest <- test_genes()
        # testgs <- Ltest$testgenes
        # 
        # NM <- matrix(0,length(testgs),length(filenames))
        # for(i in 1:length(testgs)){
        #     for(j in 1:length(filenames)){
        #         NM[i,j] <- length(intersect(testgs[[i]],geneT[[j]]))
        #     }
        # }
        
}

## unique gene number by each method
NumUniqueg <- function(){
        geneT <- geneList9(500)
        
        tmp <- sapply(1:ncut,function(kk){
                c(length(setdiff(geneT[[1]][1:kk],c(geneT[[2]][1:kk],geneT[[3]][1:kk],geneT[[9]][1:kk]))),  length(setdiff(geneT[[2]][1:kk],c(geneT[[1]][1:kk],geneT[[3]][1:kk],geneT[[9]][1:kk]))), length(setdiff(geneT[[3]][1:kk],c(geneT[[2]][1:kk],geneT[[1]][1:kk],geneT[[9]][1:kk]))), length(setdiff(geneT[[9]][1:kk],c(geneT[[2]][1:kk],geneT[[3]][1:kk],geneT[[1]][1:kk]))) )
        })
        
        tmp1 <- sapply(1:ncut,function(kk){
                c(length(setdiff(geneT[[2]][1:kk],c(geneT[[3]][1:kk],geneT[[9]][1:kk]))), length(setdiff(geneT[[3]][1:kk],c( geneT[[2]][1:kk],geneT[[9]][1:kk]))), length(setdiff(geneT[[9]][1:kk],c(geneT[[2]][1:kk],geneT[[3]][1:kk]))) )
        })
        
        source("plot_MICRF.R")
        main=""
        xlab="Rank"
        ylab="Number of unique genes"
        cols <- c("black","red","blue","green")
        legend <- c("TADA","DAWN","MAGI","MICRF-CoPrePPI")
        pdf(file="plot/NumUniqueg4.pdf",width=10,height=8)
        par(mai=c(2,1,1,1))
        plot_matrix(tmp,cols,main,xlab,ylab)
        legend("topleft",legend=legend,col=cols,lwd=2,lty=1,cex=1.5,y.intersp=1)
        dev.off()
        
        pdf(file="plot/NumUniqueg3.pdf",width=10,height=8)
        par(mai=c(2,1,1,1))
        plot_matrix(tmp1,cols[2:4],main,xlab,ylab)
        legend("topleft",legend=legend[2:4],col=cols[2:4],lwd=2,lty=1,cex=1.5,y.intersp=1)
        dev.off()
        
}

## final figure to show MICRF-CoPrePPI improve the risk gene prediction
Overlappedg <- function(){
        source("src.R")
        source("misc.R")
        ncut=246 ## TADA FDR 0.3
        geneT <- geneList9(ncut)
        genes <- unique(c(geneT[[1]],geneT[[2]],geneT[[3]],geneT[[9]]))
        
        geneT[[4]] <- geneT[[9]]
        mes <- c("TADA","DAWN","MAGI","MICRF")
        arr <- matrix(0,length(genes),5)
        arr[,1] <- genes
        for(i in 1:4){
                arr[genes %in% geneT[[i]],i+1] <- mes[i]
        }
        
        arr <- cbind(arr,4-rowSums(arr[,2:5]=="0"))
        arr <- cbind(arr,0)
        for(i in 1:dim(arr)[1]){
                if(arr[i,6]=="1"){
                        arr[i,7] <- match(arr[i,which(arr[i,2:5]!="0")+1],mes)
                }
        }
        
        #qwt(arr,file="Overlappedg.txt")
        #qwt(paste(arr[,1],arr[,6],sep="="),file="NumOfMes.txt")
        #qwt(paste(arr[,1],arr[,7],sep="="),file="IndexOfMes.txt")
        
        subC <- extract_subnets(genes,"Four246")
        netg <- union(subC[,1],subC[,2])
        for(i in 1:4){
                tmpg <- arr[arr[,6]==i,1]
                qwt(tmpg,file=paste("Detectedby",i,"mes.txt",sep=""))
        }
        for(i in 1:4){
                tmpg <- arr[arr[,7]==i,1]
                qwt(tmpg,file=paste("Detectedby",mes[i],"only.txt",sep=""))
                qwt(setdiff(tmpg,netg),file=paste("Isolatedby",mes[i],"only.txt",sep=""))
        }
        
}

geneList9 <- function(ncut,filenames=""){
        if(filenames==""){
                filenames <- c("../TADA_DAWN/result/TADAdenovo_meta_dmis.csv","../TADA_DAWN/DAWN_package/result/DAWN_meta_dmis.csv","../MAGI_V0.1_Alpha/mydata/control/RandomGeneList.1",paste("result/DDD5542/MICRFresult_",1:6,".txt",sep=""))
                #filenames <- c("result/DDD5542/TADAdenovo_meta_dmis.csv","result/DDD5542/DAWN_meta_dmis.csv","result/DDD5542/RandomGeneList.1",paste("result/DDD5542/MICRFresult_",1:6,".txt",sep=""))
        }
        
        geneT <- list()
        for(i in 1:length(filenames)){
                filename <- filenames[i]
                if(grepl(".csv",filename)){ result <- read.csv(filename); 
                }else{result <- read.table(filename);}
                
                if(i==1){
                        geneT[[i]] <- result[rank(result[,"qvalue.dn"])<=ncut,1]
                }else if(i==2){
                        result[is.na(result[,"Stage2_posterior"]),"Stage2_posterior"] <- 1
                        geneT[[i]] <- result[rank(result[,"Stage2_posterior"]) <= ncut,1]
                }else if(i==3){
                        geneT[[i]] <- result[rank(-result[,2]) <= ncut,1]
                }else{
                        geneT[[i]] <- result[1:ncut,1]
                }
        }
        
        geneT
}

extract_subnets <- function(genes,str){
        print(length(genes))
        ## gene edges in Coexp and PrePPI
        coexp <- read.table("data/network_inference/brainspan_net_top5.txt")
        subs <- coexp[,1] %in% genes | coexp[,2] %in% genes
        subco <- coexp[subs,]
        subco[,3] <- 1
        preppi <- read.table("data/PrePPI.txt")
        subs <- preppi[,1] %in% genes | preppi[,2] %in% genes
        subpre <- preppi[subs,]
        subpre[,3] <- 2
        subC <- rbind(subco,subpre)
        
        ## interactions in minimum spanning tree
        source("misc.R")
        net <- read_net(subC)
        library(igraph)
        print(setdiff(genes,net$node))
        genes <- intersect(genes,net$node)
        net$matrix <- net$matrix + t(net$matrix)
        g1 <- graph.adjacency(net$matrix, mode="undirected", weighted=NULL, diag=FALSE)
        mst <- minimum.spanning.tree(g1)
        edges <- as_edgelist(mst)
        edges <- cbind(edges,1)
        
        ## delete leaf and not risk genes
        print(dim(edges))
        unig <- names(table(c(edges[,1],edges[,2])))[table(c(edges[,1],edges[,2]))==1]
        unig <- setdiff(unig,genes)
        subs <- edges[,1] %in% unig | edges[,2] %in% unig
        edges <- edges[!subs,]
        subs <- edges[,1] %in% genes & edges[,2] %in% genes
        edges <- edges[!subs,]
        print(dim(edges))
        
        ## included all edges inter-risk genes
        subs <- coexp[,1] %in% genes & coexp[,2] %in% genes
        subco <- coexp[subs,]
        subco[,3] <- 1
        subs <- preppi[,1] %in% genes & preppi[,2] %in% genes
        subpre <- preppi[subs,]
        subpre[,3] <- 2
        subC <- rbind(subco,subpre)
        print(dim(subC))
        subC <- rbind(edges,subC)
        
        ## delete duplication edges
        net <- read_net(subC)
        print(net$size)
        net$matrix <- net$matrix + t(net$matrix)
        net$matrix[upper.tri(net$matrix,diag=TRUE)] <- 0
        sube <- which(net$matrix > 0, arr.ind=TRUE)
        subC <- cbind(net$node[sube[,1]],net$node[sube[,2]])
        qwt(subC,file=paste("CoPrePPI",str,"net.txt",sep=""))
        
        subC
}

### ROC example details
genesROC <- function(){
        source("misc.R")
        source("misc_output.R")
        gf=TRUE
        ncut=100
        TPcut=0.1;
        j=2;i=16;fprc=1
        #### above to change 
        TADAFile="../TADA_DAWN/result/TADAdenovo_meta_dmis.csv"
        Tset <- allPf(TPcut,TADAFile)
        DAWNnames <-  readLines(con <- file("../TADA_DAWN/DAWN_package/TADAdenovo_randset_1.txt","r"))   
        close(con)
        Maginames <-  readLines(con <- file("../MAGI_V0.1_Alpha/mydata/randset_1/MAGIfilenames.txt","r"))   
        close(con)
        
        filenames <- onefnew(j,i,DAWNnames,Maginames) 
        if(gf){ tmp <- read.csv(filenames[1])
        genes <- tmp[rowSums(tmp[,c("dn.LoF","dn.mis3")]) > 0,1]
        } 
        
        tmp <- c()
        for(i in 1:length(filenames)){
                filename <- filenames[i]
                if(grepl(".csv",filename)){ result <- read.csv(filename); 
                }else{result <- read.table(filename);}
                result <- cbind(result,1:dim(result)[1])
                if(gf) result <- result[result[,1] %in% genes,]
                if(i==3) result <- result[order(-result[,2]),]
                oner <- oneVec(result[1:ncut,1],Tset)
                tmp <- cbind(tmp,oner)
        }
        
        source("plot_MICRF.R")
        pdf(file="plot/onerocdetails.pdf",width=11,height=8)
        par(mai=c(2,1,1,1))
        plot_matrix(t(tmp),cols=c("black","red","blue","yellow4","brown","blueviolet","deeppink","darkgreen","darkcyan"),main="",xlab="",ylab="")
        dev.off()
        ## get cut off 60
        geneT <- geneList9(60,filenames)
        aa <- c()
        for(i in 1:8){
                aa <- union(aa,geneT[[i]])
        }
        
        
}

oneVec <- function(genes,Tset){
        ct <- rep(0,length(genes))
        ct[genes %in% Tset] <- 1
        oner <- sapply(1:length(genes), function(kk) sum(ct[1:kk]))
        oner
}

## network edges intersect
edgeInter <- function(){
        
        netsix <- list()
        source("MICRF_net.R")
        
        netsix[[1]] <- build_net(3) ## STRING
        netsix[[2]] <- build_net(6) ## iRefindex
        netsix[[3]] <- build_net(20) ## HPRD
        
        netsix[[4]] <- build_net(26) ## corr
        netsix[[5]] <- build_net(27) ## COEXP
        netsix[[6]] <- build_net(34) ## PrePPI
        
        nodes <- netsix[[1]]$node
        mapT <- as.matrix(read.delim("Fmap0121.txt",header=FALSE,sep="\t"))
        mapT <- mapT[!is.na(mapT[,1]) & !is.na(mapT[,2]),]
        subs <- which(nodes %in% mapT[,1])
        nodes[subs] <- mapT[match(nodes[subs],mapT[,1]),2]
        netsix[[1]]$node <- nodes
        dimnames(netsix[[1]]$matrix) <- list(nodes,nodes)
        
        inter.edges <- matrix(0,6,6)
        for(i in 1:6){
                for(j in 1:6){
                        inter.edges[i,j] <- net1_net2(netsix[[i]],netsix[[j]])
                }
        }
        save(inter.edges,file="Sixnetedges")
        
        ## pdf plot
        load("Sixnetedges")
        library(gplots)
        for(i in 1:6){
                inter.edges[,i] <- inter.edges[,i]/inter.edges[i,i]
        }
        my_palette <- colorRampPalette(c("white", "red"))(n = 1000)
        labs=c("STRING","iRefIndex","HPRD","CORR","CoEXP","PrePPI")
        pdf(file="EdgeInterSix.pdf",height=10,width=10)
        par(mai=c(2,1,1,1))
        heatmap.2(inter.edges,dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none',col=my_palette,labRow=labs,labCol=labs,margin=c(10,10))
        color.bar(my_palette,0,1)
        dev.off()
        
        
        deglist <- list()
        for(i in 1:6){
                deglist[[i]] <- colSums(netsix[[i]]$matrix > 0)
        }
        save(deglist,file="degSixNet")
        
        load("degSixNet")
        cols <- c("black","red","blue","yellow4","brown","deeppink")
        legend <- c("STRING","iRefIndex","HPRD","CORR","CoEXP","PrePPI")
        pdf(file="DegSix.pdf",height=10,width=10)
        par(mai=c(2,1,1,1))
        for(i in 1:6){
                fre <- table(deglist[[i]])
                deg <- names(fre)
                
                x <- sort(as.numeric(deg))
                y <- fre[sort(as.numeric(deg),index.return=TRUE)$ix]
                if(i==5){
                        lines(log10(x[x<=5]),log10(y[x<=5]),type='b',col=cols[i],lwd=2,lty=2)
                        subs <- x>=5
                        x <- x[subs]
                        y <- y[subs]
                }
                
                if(i==1){
                        plot(log10(x),log10(y),,col=cols[i],type="l",main="Degree distributions of biological networks",xlab="Degree (log10-scale)",ylab="Frequency (log10-scale)",cex.lab=1.5,cex.axis=1.35,cex.main=2,lwd=2,ylim=c(0,4))
                }else{
                        lines(log10(x),log10(y),type='l',col=cols[i],lwd=2)
                }
        }
        legend("topright",legend=legend,col=cols,lwd=2,lty=1,cex=1.5,y.intersp=1)
        dev.off()
        
}

net1_net2 <- function(net1,net2){
        nodes <- intersect(net1$node,net2$node)
        sum((net1$matrix[nodes,nodes] > 0) & (net2$matrix[nodes,nodes] > 0))/2
}

# Nature 13772 
Nat13772 <- function(){
        
        filenames <- c("../TADA_DAWN/result/TADAdenovo_nat8_20.csv","../TADA_DAWN/DAWN_package/result/DAWN_Nat.csv","../MAGI_V0.1_Alpha/mydata/control/RandomGeneList.4",paste("result/nat13772/MICRFresult_",1:6,".txt",sep=""))
        oneresult <-  read.csv("result/nature13908.csv")
        ncut <- 500
        
        atmp <- matrix(0,length(filenames),ncut);
        for(i in 1:length(filenames)){
                filename <- filenames[i]
                if(grepl(".csv",filename)){ result <- read.csv(filename); 
                }else{result <- read.table(filename);}
                
                if(i==3){
                        result <- result[order(-result[,2]),]
                }
                
                tmpg <- result[1:ncut,1]
                inrtmp <- rep(0, length(tmpg))
                inrtmp[tmpg %in% oneresult[,"Gene"]] <- 1 #table(oneresult[,"Gene"])[tmpg[tmpg %in% oneresult[,"Gene"]]]
                atmp[i,] <- sapply(1:length(tmpg), function(kk) sum(inrtmp[1:kk])) 
        }
        
        source("plot_MICRF.R")
        cols <- c("black","red","blue","yellow4","brown","blueviolet","deeppink","darkgreen","darkcyan")
        legend <- c("TADA","DAWN","MAGI","MICRF-STRING","MICRF-iRefIndex","MICRF-HPRD","MICRF-CORR","MICRF-CoEXP","MICRF-CoPrePPI")
        pdf(file="plot/re13908.pdf",width=11,height=8)
        par(mai=c(2,1,1,1))
        main=""
        xlab="Rank"
        ylab="Number of recurrent genes"
        plot_matrix(atmp,cols,main,xlab,ylab)
        legend("topleft",legend=legend,col=cols,lwd=2,lty=1,cex=1.3,y.intersp=0.9)
        dev.off()
}

# DDD repeat
DDD_repeat <- function(){
        ncut <- 500
        filenames <- c("../TADA_DAWN/result/TADAdenovo_meta_dmis.csv","../TADA_DAWN/DAWN_package/result/DAWN_meta_dmis.csv","../MAGI_V0.1_Alpha/mydata/control/RandomGeneList.1",paste("result/DDD5542/MICRFresult_",1:6,".txt",sep=""))
        phenos <- c("Abnormality of the cardiovascular system","Abnormality of the nervous system","Abnormality of the peripheral nervous system","Autism spectrum disorder","Growth abnormality","Multiple congenital anomalies","Seizures")  
        a <- read.delim("DDD_mutations/Table3Sup.txt")
        oneresult <- a[a[,"Primary.phenotype"] %in% phenos,]
        
        atmp <- matrix(0,length(filenames),ncut);
        for(i in 1:length(filenames)){
                filename <- filenames[i]
                if(grepl(".csv",filename)){ result <- read.csv(filename); 
                }else{result <- read.table(filename);}
                
                if(i==3){
                        result <- result[order(-result[,2]),]
                }
                
                tmpg <- result[1:ncut,1]
                inrtmp <- rep(0, length(tmpg))
                inrtmp[tmpg %in% oneresult[,"Gene"]] <- 1 #table(oneresult[,"Gene"])[tmpg[tmpg %in% oneresult[,"Gene"]]]
                atmp[i,] <- sapply(1:length(tmpg), function(kk) sum(inrtmp[1:kk])) 
        }
        
        
        genes <- tmpg[inrtmp>0]
        qwt(genes,file="reDDD78.txt")
        
        tadar <- read.csv(filenames[1])
        newgenes <- setdiff(genes,tadar[1:121,1])
        qwt(newgenes,file="newDDDreg33.txt")
        
        
        source("plot_MICRF.R")
        cols <- c("black","red","blue","yellow4","brown","blueviolet","deeppink","darkgreen","darkcyan")
        legend <- c("TADA","DAWN","MAGI","MICRF-STRING","MICRF-iRefIndex","MICRF-HPRD","MICRF-CORR","MICRF-CoEXP","MICRF-CoPrePPI")
        pdf(file="plot/reDDD.pdf",width=11,height=8)
        par(mai=c(2,1,1,1))
        main=""
        xlab="Rank"
        ylab="Number of recurrent genes"
        plot_matrix(atmp,cols,main,xlab,ylab)
        legend("bottomright",legend=legend,col=cols,lwd=2,lty=1,cex=1.5,y.intersp=1)
        dev.off()
        
        
        
        
}

## 121 and new 33 DDD genes in CoPrePPI network
subCoPrePPI <- function(){
        
        filename <- "../TADA_DAWN/result/TADAdenovo_meta_dmis.csv"
        genes <- unlist(read.table("reDDD78.txt"))
        tadar <- read.csv(filename)
        newgenes <- setdiff(genes,tadar[1:121,1])
        allgenes <- union(genes,tadar[1:121,1])
        
        nodecol <- cbind(allgenes,1)
        nodecol[nodecol[,1] %in% newgenes,2]  <- 2
        source("MICRF_net.R")
        net1 <- build_net(36)
        allgenes <- allgenes[allgenes %in% net1$node]
        
        deg <- colSums(net1$matrix)
        nodedeg <- cbind(allgenes,deg[match(allgenes,net1$node)])
        
        subnet <- net1$matrix[allgenes,allgenes]
        subnet[upper.tri(subnet)] <- 0
        sube <- which(subnet > 0, arr.ind=TRUE)
        edges <- cbind(allgenes[sube[,1]],allgenes[sube[,2]])
        
        qwt(paste(nodecol[,1]," = ",nodecol[,2],sep=""),file="NewDDDnodeCol.txt",flag=2)
        qwt(paste(nodedeg[,1]," = ",nodedeg[,2],sep=""),file="NewDDDnodedeg.txt",flag=2)
        qwt(edges,file="NewDDDsubnet.txt")
        
}

color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
        scale = (length(lut)-1)/(max-min)
        
        dev.new(width=1.75, height=5)
        plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
        axis(2, ticks, las=1)
        for (i in 1:(length(lut)-1)) {
                y = (i-1)/scale + min
                rect(0,y,1,y+1/scale, col=lut[i], border=NA)
        }
}

### network betweenness and degree centrality comparison
CentralityCom <- function(){
        source("misc_output.R")
        source("MICRF_net.R")
        source("src.R")
        source("misc.R")
        
        netfs <- c("data/STRINGnetmap.txt","data/hotnet/iRefIndexm.txt","data/StringNew_HPRD.txt","data/network_inference/brainspan_net_cor.txt","data/network_inference/brainspan_net_top5.txt","data/network_inference/ComCo_PrePPI.txt")
        netflags <- c(3,6,20,26,27,36)
        nodebefs <- c("data/Network_betweenness/Betweenness_node_STRING.txt","data/Network_betweenness/Betweenness_node_iRef.txt","data/Network_betweenness/Betweenness_node_HPRD.txt","data/Network_betweenness/Betweenness_node_corr1.txt","data/Network_betweenness/Betweenness_node_coexp5.txt","data/Network_betweenness/Betweenness_node_CoPrePPI.txt")
        nodeclfs <- c("data/Network_betweenness/Closeness_node_STRING.txt","data/Network_betweenness/Closeness_node_iRefIndex.txt","data/Network_betweenness/Closeness_node_HPRD.txt","data/Network_betweenness/Closeness_node_CORR.txt","data/Network_betweenness/Closeness_node_CoEXP.txt","data/Network_betweenness/Closeness_node_CoPrePPI.txt")
        edgebefs <- c("data/Network_betweenness/Betweenness_edge_STRING.txt","data/Network_betweenness/Betweenness_edge_iRef.txt", "data/Network_betweenness/Betweenness_edge_HPRD.txt", "data/Network_betweenness/Betweenness_edge_corr1.txt","data/Network_betweenness/Betweenness_edge_coexp5.txt","data/Network_betweenness/Betweenness_edge_Co_PrePPI.txt")
        
        alist <- list()
        for(i in 1:6){
                if(i==1){idm=TRUE;}else{idm=FALSE;}
                alist[[i]] <- OneNetAlist(netflags[i],nodebefs[i],nodeclfs[i],edgebefs[i],idm)
        }
        

        #### one boxplot example
        ncontrol=1
        gset <- Geneset(flag=2,ncontrol)
        
        k=1
        nset <- length(gset)
        degBoxL <- list()
        betBoxL <- list()
        cloBoxL <- list()
        edgebeBoxL <- list()
        for(i in 1:6){
                tmp <- OnenetList(gset,alist[[i]])
                for(j in 1:nset){
                        degBoxL[[k]] <- as.numeric(tmp[[1]][[j]]); 
                        cloBoxL[[k]] <- as.numeric(tmp[[2]][[j]]);
                        betBoxL[[k]] <- as.numeric(tmp[[3]][[j]]); 
                        edgebeBoxL[[k]] <- as.numeric(tmp[[4]][[j]]);
                        k <- k + 1        
                }
        }
        CentralityOne(degBoxL,betBoxL,cloBoxL,edgebeBoxL,wstr=paste("NetworkCentralities",ncontrol,sep=""),nset)
        
        ### 500 random non-risk gene sets mean fold-change curves
        source("src.R")
        ncontrol=1
        n.sim = 100
        TPcut=0.3
        
        TADAFile="../TADA_DAWN/result/TADAdenovo_meta_dmis.csv"
        tadar <- read.csv(TADAFile)
        tadar[is.na(tadar[,"qvalue.dn"]),"qvalue.dn"] <- 1
        Tset <- tadar[tadar[,"qvalue.dn"] < TPcut,1]
        subs <- tadar[,"qvalue.dn"] > TPcut & tadar[,"qvalue.dn"] < 0.9
        exgenes <- "" #tadar[subs,1]
        
        MsList <- list()
        for(i in 1:4){
            MsList[[i]] <-  matrix(0,n.sim,6)  
        }
        #scaMs <- matrix(0,3,6)
        for(kk in 1:n.sim){
                gset <- Geneset(flag=2,ncontrol,TPcut)
                nset <- length(gset)
                for(i in 1:6){
                        tmp <- OnenetList(gset,alist[[i]],exgenes=exgenes)
                        for(j in 1:4){
                                a1 <- as.numeric(tmp[[j]][[1]])
                                a2 <- as.numeric(tmp[[j]][[2]])
                                if(all(a1==0)) a1=1
                                if(all(a2==0)) a2=1
                                if(j==3 | j==4){ MsList[[j]][kk,i] <- median(log(a1[a1>0]))/median(log(a2[a2>0]))
                                        }else{
                                        MsList[[j]][kk,i] <- median(log(a1))/median(log(a2)) #!!!!
                                        }
                                #if(is.na(MsList[[j]][kk,i])) stop("Stoped")
                        }
                }
        }
        wstr=paste("NetworkCentralitiesFoldsnew",ncontrol,"_",n.sim,"_",TPcut,sep="")
        BatchesMeanFolds(MsList,wstr)
        
        
        #------------------------pathways
#         sixgL <- sixPathways()
#         sixstr <- c("WntPathway","TGF-betaPathway","CalciumPathway","RasPathway","NotchPathway","MAPKPathway")
#         for(kk in 1:6){
#                 gset1 <- list()
#                 for(i in 1:3){ gset1[[i]] <- intersect(gset[[i]],sixgL[[kk]]); }
#                 k=1
#                 degBoxL <- list()
#                 betBoxL <- list()
#                 cloBoxL <- list()
#                 for(i in 1:6){
#                         if(i==1){idm=TRUE;}else{idm=FALSE;}
#                         tmp <- OnenetList(gset1,netfs[i],netflags[i],nodebefs[i],nodeclfs[i],idm)
#                         degBoxL[[k]] <- tmp[[1]][[1]]; degBoxL[[k+1]] <- tmp[[1]][[2]]; degBoxL[[k+2]] <- tmp[[1]][[3]];
#                         betBoxL[[k]] <- tmp[[2]][[1]]; betBoxL[[k+1]] <- tmp[[2]][[2]]; betBoxL[[k+2]] <- tmp[[2]][[3]];
#                         cloBoxL[[k]] <- tmp[[3]][[1]]; cloBoxL[[k+1]] <- tmp[[3]][[2]]; cloBoxL[[k+2]] <- tmp[[3]][[3]]; 
#                         k=k+3;
#                 }
#                 CentralityOne(degBoxL,betBoxL,cloBoxL,wstr=sixstr[kk])
#         }
        
}

BatchesMeanFolds <- function(MsList,wstr){
        nets <- c("STRING","iRefIndex","HPRD","CORR","CoEXP","CoPrePPI")
        cols <- c("red","brown","deeppink","cyan")
        pdf(file=paste("plot/",wstr,".pdf",sep=""),width=10,height = 8)
        par(mai=c(1.2,1,1,1))
        n <- dim(MsList[[1]])[2]
        x <- 1:n
        
        natt <- length(MsList)
        yMs <- matrix(0,natt,6)
        sdMs <- matrix(0,natt,6)
        for(i in 1:natt){
                for(j in 1:dim(MsList[[i]])[2]){
                        yMs[i,j] <- mean(MsList[[i]][MsList[[i]][,j]<Inf,j])   
                        sdMs[i,j] <- sd(MsList[[i]][MsList[[i]][,j]<Inf,j])
                }
        }
        ymax <- max(yMs+sdMs)+0.05
        
        mp <- barplot(yMs, beside=T, ylab="",cex.names=1.5, ylim=c(0.999,ymax), col=cols,xaxt='n',xpd = FALSE)
        #text(colMeans(mp), par("usr")[3], labels = nets, srt = 45, pos = 1)
        #axis(1,at=colMeans(mp),labels=nets)
#         for(i in 1:natt){
#                 if(i==1){
#                         plot(x,yMs[i,],col=i,xaxt='n',xlab="",ylab="Folds",main="",ylim=c(0.9,ymax),type="b",cex.lab=1.5)
#                         axis(1,at=1:6,labels=nets,cex=2)
#                         adderr(x,yMs[i,],sdMs[i,],col=i)
#                 }else{
#                         lines(x,yMs[i,],col=i,type="b")
#                         adderr(x,yMs[i,],sdMs[i,],col=i)  
#                 }
#         }
        for(i in 1:natt){
                adderr(mp[i,],yMs[i,],sdMs[i,],col=cols[i])
        }
        legend(x=20.6,y=1.4,legend=c("Degree","Closeness","Betweenness","EdgeBetweenness"),cex=1.3,fill=cols,bty="n")
        dev.off()
}

adderr <- function(x,y,sd,col){
        segments(x,y-sd,x,y+sd,col=col,lwd=2)
        epsilon <- 0.4
        segments(x-epsilon,y-sd,x+epsilon,y-sd,col=col,lwd=2)
        segments(x-epsilon,y+sd,x+epsilon,y+sd,col=col,lwd=2)
}

CentralityOne <- function(degBoxL,betBoxL,cloBoxL,edgebeBoxL,wstr,nset=2){
        
        cols=c(2,3,5)
        degBoxL1 <- list()
        betBoxL1 <- list()
        cloBoxL1 <- list()
        edgebeBoxL1 <- list()
        for(i in 1:length(degBoxL)){
                degBoxL1[[i]] <- log(degBoxL[[i]]); 
                if(max(degBoxL1[[i]])==Inf) degBoxL1[[i]][degBoxL1[[i]]==Inf] <- max(degBoxL1[[i]][degBoxL1[[i]] < Inf]) + 1; 
                if(min(degBoxL1[[i]])==-Inf) degBoxL1[[i]][degBoxL1[[i]]== -Inf] <- min(degBoxL1[[i]][degBoxL1[[i]] > -Inf]) - 1;
                betBoxL1[[i]] <- log(betBoxL[[i]]); 
                if(max(betBoxL1[[i]])==Inf) betBoxL1[[i]][betBoxL1[[i]]==Inf] <- max(betBoxL1[[i]][betBoxL1[[i]] < Inf]) + 1; 
                if(min(betBoxL1[[i]])==-Inf) betBoxL1[[i]][betBoxL1[[i]]== -Inf] <- min(betBoxL1[[i]][betBoxL1[[i]] > -Inf]) - 1;
                cloBoxL1[[i]] <- log(cloBoxL[[i]]); 
                if(max(cloBoxL1[[i]])==Inf) cloBoxL1[[i]][cloBoxL1[[i]]==Inf] <- max(cloBoxL1[[i]][cloBoxL1[[i]] < Inf]) + 1; 
                if(min(cloBoxL1[[i]])==-Inf) cloBoxL1[[i]][cloBoxL1[[i]]== -Inf] <- min(cloBoxL1[[i]][cloBoxL1[[i]] > -Inf]) - 1;
                edgebeBoxL1[[i]] <- log(edgebeBoxL[[i]]); 
                if(max(edgebeBoxL1[[i]])==Inf) edgebeBoxL1[[i]][edgebeBoxL1[[i]]==Inf] <- max(edgebeBoxL1[[i]][edgebeBoxL1[[i]] < Inf]) + 1; 
                if(min(edgebeBoxL1[[i]])==-Inf) edgebeBoxL1[[i]][edgebeBoxL1[[i]]== -Inf] <- min(edgebeBoxL1[[i]][edgebeBoxL1[[i]] > -Inf]) - 1;
        }
        
        
        pdf(file=paste("plot/",wstr,"1.pdf",sep=""),width=15,height = 8)
        par(mfrow=c(2,6),mai=c(0.1,0.3,1,0.3))
        for(i in 1:6){
                onelist <- list()
                for(j in 1:nset){onelist[[j]] <- degBoxL1[[(i-1)*nset+j]];}
                aa=boxplot(onelist,outline=FALSE,xaxt='n',xlab='',ylab='');
                lab <- ifelse(t.test(onelist[[1]],onelist[[2]])$p.value < 0.0001, "p-value < 0.0001",paste("p-value : ",round(t.test(onelist[[1]],onelist[[2]])$p.value,digits=4),sep=""))
                #text(x=1.5,y=(max(aa$stats[4,])+max(aa$stats[5,]))/2,labels=lab,cex=1.7)
        }
        
        par(mai=c(1.5,0.3,0.1,0.3))
        for(i in 1:6){  
                onelist <- list()
                for(j in 1:nset){onelist[[j]] <- cloBoxL1[[(i-1)*nset+j]];}
                aa=boxplot(onelist,outline=FALSE,xaxt='n',xlab='',ylab='');
                lab <- ifelse(t.test(onelist[[1]],onelist[[2]])$p.value < 0.0001, "p-value < 0.0001",paste("p-value : ",round(t.test(onelist[[1]],onelist[[2]])$p.value,digits=4),sep=""))
                #text(x=1.5,y=(min(aa$stats[1,])+min(aa$stats[2,]))/2,labels=lab,cex=1.7)
        } 
        dev.off()
        
        pdf(file=paste("plot/",wstr,"2.pdf",sep=""),width=15,height = 8)
        par(mfrow=c(2,6),mai=c(0.1,0.3,1,0.3))
        for(i in 1:6){  
                onelist <- list()
                for(j in 1:nset){onelist[[j]] <- betBoxL1[[(i-1)*nset+j]];}
                aa=boxplot(onelist,outline=FALSE,xaxt='n',xlab='',ylab='');
                lab <- ifelse(t.test(onelist[[1]],onelist[[2]])$p.value < 0.0001, "p-value < 0.0001",paste("p-value : ",round(t.test(onelist[[1]],onelist[[2]])$p.value,digits=4),sep=""))
                #text(x=1.5,y=(max(aa$stats[4,])+max(aa$stats[5,]))/2,labels=lab,cex=1.7)
        }
        par(mai=c(1.5,0.3,0.1,0.3))
        for(i in 1:6){  
                onelist <- list()
                for(j in 1:nset){onelist[[j]] <- edgebeBoxL1[[(i-1)*nset+j]];}
                aa=boxplot(onelist,outline=FALSE,xaxt='n',xlab='',ylab='');
                lab <- ifelse(t.test(onelist[[1]],onelist[[2]])$p.value < 0.0001, "p-value < 0.0001",paste("p-value : ",round(t.test(onelist[[1]],onelist[[2]])$p.value,digits=4),sep=""))
                #text(x=1.5,y=(min(aa$stats[1,])+min(aa$stats[2,]))/2,labels=lab,cex=1.7)
                axis(1, at=1:nset,labels=c("Risk","Non-risk"),cex.axis=2)
        }
        dev.off()
        
}

Geneset <- function(flag=2,ncontrol=1,TPcut=0.1){
        
        TADAFile="../TADA_DAWN/result/TADAdenovo_meta_dmis.csv"
        tadar <- read.csv(TADAFile)
        gset <- list()
        if(flag==1){
                gset[[1]] <- tadar[rowSums(tadar[,c("dn.LoF","dn.mis3")]) >= 2,1]
                gset[[2]] <- tadar[rowSums(tadar[,c("dn.LoF","dn.mis3")]) == 1,1]
                gset[[3]] <- tadar[rowSums(tadar[,c("dn.LoF","dn.mis3")]) == 0,1]
        }else if(flag==2){
                ## risk genes
                #TPcut=0.2
                tadar[is.na(tadar[,"qvalue.dn"]),"qvalue.dn"] <- 1
                subs <- tadar[,"qvalue.dn"] < TPcut
                gset[[1]] <- tadar[subs,1]
                r1 <- -log10(rowSums(tadar[subs,2:3]))
                
                ## non-risk genes
                TPcut2 <- 0.9
                subs <- tadar[,"qvalue.dn"] > TPcut2
                gtmp <- tadar[subs,1]
                r2 <- -log10(rowSums(tadar[subs,2:3]))
                
                dr2 <- dnorm(r2,mean=mean(r1),sd=sd(r1))
                dr2 <- dr2/sum(dr2)
                subs <- sample(1:length(r2),ncontrol*length(r1), replace = FALSE, prob = dr2)
                gset[[2]] <- gtmp[subs]
                
                pdf("plot/risk_randon_nonrisk_mutationrateDIS.pdf",width=10,height=8)
                ymax <- max(density(r1)$y,density(r2[subs])$y)*1.01
                xmin <- min(density(r1)$x,density(r2[subs])$x)*0.99
                xmax <- max(density(r1)$x,density(r2[subs])$x)*1.01
                plot(density(r1),col=1,xlim=c(xmin,xmax),ylim=c(0,ymax),xlab="",ylab="",main="",cex.lab=1.4)
                lines(density(r2[subs]),col=2)
                legend("topright",legend=c("Risk","Non-risk"),col=1:2,lwd=rep(1,2),lty=rep(1,2),cex=2)
                dev.off()
        }
        

        
        gset
}

OneNetAlist <- function(netflag,nodebef,nodeclf,edgebef,idm=FALSE){
        
        net <- build_net(netflag)
        net$matrix <- net$matrix + t(net$matrix)
        
        deg <-  colSums(net$matrix > 0)
        be <- read.delim(nodebef,header=FALSE)
        cl <- read.delim(nodeclf,header=FALSE)
        edgebe <- read.table(edgebef)
        
        if(idm){
                net$node <- idmapf(net$node)
                be[,1] <- idmapf(be[,1])
                cl[,1] <- idmapf(cl[,1])
                edgebe[,1] <- idmapf(edgebe[,1])
                edgebe[,3] <- idmapf(edgebe[,3])
        }
        edgebe <- edgebe[,c(1,3,4)]
        deg <- cbind(net$node,deg)
        
        tmp <- edgebe[,c(2,1,3)]
        tmp <- rbind(edgebe,tmp)
        edgeNet <- read_net(tmp)
        diag(edgeNet$matrix) <- 0
        #edgebe1 <- cbind(edgeNet$node, rowSums(edgeNet$matrix))
        edgebe1 <- edgeNet$matrix
        #genes <- union(edgebe[,1],edgebe[,2])
        #edbe <- sapply(1:length(genes), function(kk) max(edgebe[edgebe[,1] %in% genes[kk] | edgebe[,2] %in% genes[kk],3]))
        #edbe[is.na(edbe)] <- 1
        #edgebe1 <- cbind(genes,edbe)
        
        list(deg,cl,be,edgebe1)
}

OnenetList <- function(gset,onenetL,exgenes=""){
        deg <- onenetL[[1]]
        cl <- onenetL[[2]]
        be <- onenetL[[3]]
        edgebe <- onenetL[[4]]
        
        degL <- list()
        beL <- list()
        clL <- list()
        edgebeL <- list()
        for(i in 1:length(gset)){
                degL[[i]] <- deg[deg[,1] %in% gset[[i]],2]
                beL[[i]] <- be[be[,1] %in% gset[[i]],2]
                clL[[i]] <- cl[cl[,1] %in% gset[[i]],2]
                #edgebeL[[i]] <- edgebe[edgebe[,1] %in% gset[[i]] & edgebe[,2] %in% gset[[i]],3]  #!!!!      
                #edgebeL[[i]] <- edgebe[edgebe[,1] %in% gset[[i]],2]  #!!!! 
                tmp <- intersect(rownames(edgebe),gset[[i]])
                edgebeL[[i]] <- rowSums(edgebe[tmp,])  
                #edgebeL[[i]] <- apply(edgebe[tmp,tmp],1,max)  
                #edgebeL[[i]] <- sapply(1:length(tmp), function(kk) mean(edgebe[tmp[kk], edgebe[tmp[kk],]>0 ])) 
                #edgebeL[[i]] <- sapply(1:length(tmp), function(kk) median(edgebe[tmp[kk], edgebe[tmp[kk],]>0 ]))
                #edgebeL[[i]] <- edgebeL[[i]][edgebeL[[i]] > 0]
        }
        #tmp1 <- intersect(rownames(edgebe),gset[[1]])
        #tmp2 <- intersect(rownames(edgebe),gset[[2]])
        #edgebeL[[1]] <- apply(edgebe[tmp1,setdiff(rownames(edgebe),c(gset[[2]],exgenes))],1,max)
        #edgebeL[[2]] <- apply(edgebe[tmp2,setdiff(rownames(edgebe),c(gset[[1]],exgenes))],1,max)
        
        list(degL,clL,beL,edgebeL)
}

idmapf <- function(nodes){
        mapT <- as.matrix(read.delim("Fmap0121.txt",header=FALSE,sep="\t"))
        mapT <- mapT[!is.na(mapT[,2]), ]
        
        subnodes <- intersect(nodes,mapT[,1])
        subs <- which(nodes %in% subnodes)
        nodes[subs] <- mapT[match(nodes[subs],mapT[,1]),2]
        nodes
}

### whether that risk genes have highly ranked centrality in their pathways to non-risk genes
PathwaysAna <- function(alist,Tset,flag=1){
        ## flag=1: Only plot risk genes in pathways and networks
        ## flag=2: plot both risk and non-risk genes in pathways and networks
        
        files <- list.files("./pathways/l5L100_121","*_pathways.txt$",full.names = TRUE)
        sixKf <- c("Wnt_pathways.txt","Calcium_pathways.txt","MAPK_pathways.txt","Notch_pathways.txt","Ras_pathways.txt","TGF-beta_pathways.txt")
        files <- files[!(basename(files) %in% sixKf)]
        
        np <- length(files)
        nets <- c("STRING","iRefIndex","HPRD","CORR","CoEXP","CoPrePPI")
        meas <- c("Degree","Closeness","Betweenness")
        
        if(flag==2){
                ### get risk genes and non-risk genes in all pathways
                Apg <- c() # ALL pathway genes
                for(i in 1:length(files)){
                        Apg <- union(Apg,unlist(read.table(files[i])))      
                }
        }
        
        for(i in 1:length(nets)){
                for(j in 1:length(meas)){
                        pdf(file=paste("plot/pathways/",nets[i],"_",meas[j],"_pathwayRank.pdf",sep=""),height=10,width = 10)
                        par(mai=c(2,1,1,1))
                        netR <- alist[[i]][[j]]
                        
                        if(flag==1) genes <- Tset[Tset %in% netR[,1]]
                        if(flag==2) genes <- Apg[Apg %in% netR[,1]]
                        
                        rankM <- matrix(-1,length(genes),np+1)
                        tmp <- rank(-as.numeric(netR[,2]))/length(netR[,2])
                        rankM[,1] <- tmp[match(genes,netR[,1])]
                        for(k in 1:length(files)){
                                tmpg <- unlist(read.table(files[k]))
                                tmpR <- netR[netR[,1] %in% tmpg,]
                                tmp <- rank(-as.numeric(tmpR[,2]))/length(tmpR[,2])
                                if(flag==1) { 
                                        tmprg <- intersect(tmpR[,1],Tset)
                                        rankM[match(tmprg,genes),k+1] <- tmp[match(tmprg,tmpR[,1])]
                                }
                                if(flag==2) rankM[match(tmpR[,1],genes),k+1] <- tmp
                        }
                        if((np+1) > 2 ) y <- sapply(1:length(genes), function(kk){aa=rankM[kk,2:(np+1)];min(aa[aa>0]);})
                        if((np+1) == 2) y <- rankM[,2];
                        subs <- !is.na(y) & (y < Inf) & (y>0)
                        ymax <- 1 #ymax <- max(c(-log(rankM[subs,1]),-log(y[subs])))
                        if(flag==2){
                                subs1 <- subs & (genes %in% Tset)
                                subs2 <- !subs1
                                krisk <- sum(subs1)
                                plot(rankM[subs1,1],y[subs1],xlab=paste("Top ranks in ",nets[i]," network",sep=""),ylab="Top ranks in pathways",main=paste(meas[j]," : ",krisk, " risk genes in ",np, " pathways",sep=""),xlim=c(0,ymax),ylim=c(0,ymax),cex.lab=2,cex.axis=2,col=1,cex.main=2)
                                lines(rankM[subs2,1],y[subs2],type="p",col=2)
                        }
                        
                        if(flag==1){
                                krisk <- sum(subs)
                                plot(rankM[subs,1],y[subs],xlab=paste("Top ranks in ",nets[i]," network",sep=""),ylab="Top ranks in pathways",main=paste(meas[j]," : ",krisk, " risk genes in ",np, " pathways",sep=""),xlim=c(0,ymax),ylim=c(0,ymax),cex.lab=2,cex.axis=2,col=1,cex.main=2)
                        }
                        #library(RColorBrewer)
                        #k <- 11
                        #my.cols <- rev(brewer.pal(k, "RdYlBu"))
                        #smoothScatter(cbind(rankM[subs,1],y[subs]), nrpoints=.3*sum(subs), colramp=colorRampPalette(my.cols), pch=19, xlab=paste("Top ", meas[j]," ranks in ",nets[i]," network",sep=""),ylab=paste("Top ", meas[j]," ranks in pathways",sep=""),main="",xlim=c(0,ymax),ylim=c(0,ymax),xaxs="i",yaxs="i",cex.lab=2,cex.axis=2)
                        lines(c(0,ymax),c(0,ymax),type="l",col=1,lty=2,lwd=2)
                        dev.off()
                }
        }
        
}

sixPathways <- function(){
        psf <- c("data/pathways/new6/Wnt_pathway.txt","data/pathways/new6/TGF-beta_pathway.txt","data/pathways/new6/Calcium_pathway.txt","data/pathways/new6/Ras_pathway.txt","data/pathways/new6/Notch_pathway.txt","data/pathways/new6/MAPK_pathway.txt") 
        sixgL <- list()
        for(i in 1:length(psf)){
                tmp <- read.delim(psf[i],header=FALSE,sep=";",quote="")[,1]
                tmp <- unique(gsub("\"","",tmp))
                genes <- sapply(1:length(tmp), function(kk) { a1 <- unlist(strsplit(tmp[kk],","))[1]; unlist(strsplit(a1," "))[1];})
                sixgL[[i]] <- genes
                qwt(genes,file=paste("pathways/",gsub("_pathway.txt","_pathways.txt",basename(psf[i])),sep=""))
        }
        #sixgL
}

## 1 violin plot of ranks in network and pathways


## 2 compare genes with damaging de novo mutations in cases vs controls in these context


