ASD_data_set <- function(){
    options(stringsAsFactors=FALSE)
    ## read the families used in 13908, supplementary table 1
    fam13908 <- read.csv("ASD/nature13908/Supplementary Table 1.csv")[,1]
    fam13908 <- unique(fam13908)
    ## read the families used in other papers, supplementary table 10
    #fam13908subset <- read.csv("ASD/nature13908/Supplementary Table 10.csv") ## 208 only 1 less
    ## read the families from lossifov, 343 families
    famlos <- read.csv("ASD/Neuron/mmc2.csv")[,1]
    famlos <- unique(famlos)
    ## read the families from sanders nature 10945, 238 families
    famsan <- read.csv("ASD/nature10945/nature10945-s2.csv")[,1]
    famsan <- unique(famsan)
    ## read the familes from O' Roak nature 10989, 189 + 20  families
    deroa <- read.csv("ASD/nature10989/nature10989-s1.csv")
    deroa <- deroa[deroa[,"type"] %in% "SSC189",]
    famroa <- deroa[,1]
    famroa <- sapply(1:length(famroa),function(i) {unlist(strsplit(famroa[i],"\\."))[1]})
    famroa <- unique(famroa)
    
    fam13908F <- setdiff(fam13908,union(famlos,union(famsan,famroa)))
    length(fam13908F) ## 1751
    ## get the denovo mutations for this subset familes 
    ## here, we do not distinguish the two type of missense mutations!!!!
    ### Only denovo mutations are considered.
    
    
    ## Neale et al., 174 trios, 18 lof and 101 mis;
    #     deNea <- read.csv("ASD/nature11011/nature11011-s2.csv")
    #     lofNea <- c("Frameshift","Nonsense","Splice Site")
    #     misNea <- c("Missense")
    #     deNea <- deNea[deNea[,"Category"] %in% union(lofNea,misNea),c("Child_ID","Gene_name","Category")]
    
    ##  Lossifov et al., 343 tiros, 
    delos1 <- read.csv("ASD/Neuron/mmc2.csv") ## SNV
    delos1 <- delos1[delos1[,"denovoScore"]>=60 & delos1[,"chi2APval"]>=0.0001 & delos1[,"pg.parentsCalledCnt"] >= 300 & delos1[,"pg.altPrcnt"] <= 0.2 & delos1[,"noise.rate.prcnt"] <=0.3,]
    inchild <- c("autF","autFsibF","autFsibM","autM","autMsibF","autMsibM")
    loflos1 <- c("nonsense","splice_site")
    mislos1 <- c("missense")
    delos1 <- delos1[delos1[,"inChild"] %in% inchild,]
    delos1 <- delos1[delos1[,"effectType"] %in% union(loflos1,mislos1),]
    subs <- !grepl("^benign:",delos1[,"PolyPhen2.prob.site.region"])
    delos1 <- delos1[subs,]
    delos1 <- delos1[,c("quadId","inChild","effectType","effectGenes")]
    
    genes <- delos1[,"effectGenes"]
    genes <- gsub(" ","",genes)
    genes <- sapply(genes, function(gene) {
        b <- unlist(strsplit(gsub("[:|]","\t",gene),"\t"))
        subs <- which(b[seq(2,length(b),2)] %in% union(loflos1,mislos1))
        b[subs[1]]
    })
    delos1[,"effectGenes"]   <- genes  
    #delos1 <- data.frame(delos1,"MutaType")
    delos1[delos1[,"effectType"] %in% loflos1,"MutaType"] <- "LOF"
    delos1[delos1[,"effectType"] %in% mislos1,"MutaType"] <- "MIS"
    delos1[,"From"] <- "Lossifov et al."
    
    delos2 <- read.csv("ASD/Neuron/mmc4.csv") ## indel
    delos2 <- delos2[delos2[,"IndelFilter"]==TRUE,]   
    inchild <- c("autF","autFsibF","autFsibM","autM","autMsibF","autMsibM")
    loflos2 <- c("frame-shift","splice-site","no-frame-shift-new-stop","no-frame-shift-new-Stop") ## no-frame-shift is no included here.
    delos2 <- delos2[delos2[,"inChild"] %in% inchild,]
    delos2 <- delos2[delos2[,"effectType"] %in% loflos2,]
#     subs <- !grepl("^benign:",delos2[,"PolyPhen2.prob.site.region"])
#     delos2 <- delos2[subs,]
    delos2 <- delos2[,c("quadId","inChild","effectType","effectGenes")] 
    
    genes <- delos2[,"effectGenes"]
    genes <- gsub(" ","",genes)
    genes <- sapply(genes, function(gene) {
        b <- unlist(strsplit(gsub("[:|]","\t",gene),"\t"))
        subs <- which(b[seq(2,length(b),2)] %in% loflos2)
        b[subs[1]]
    })
    delos2[,"effectGenes"]   <- genes
    delos2[delos2[,"effectType"] %in% loflos2,"MutaType"] <- "LOF"
    delos2[,"From"] <- "Lossifov et al."
    
    ##  Sanders et al., 238 trios,  ## two type of missense mutations !!!!
    desan <- read.csv("ASD/nature10945/nature10945-s3.csv")
    lofsan <- c("Frameshift","Nonsense","Splice_Site")
    missan <- c("Missense")
    desan <- desan[desan[,"Effect"] %in% union(lofsan,missan),]
    dam <- c("probably damaging","possibly damaging","#N/A","0","Nonsense","possiblydamaging","probablydamaging","Splice_Site","unknown")
    desan <- desan[desan[,"PolyPhen2"] %in% dam,]
    desan <- desan[,c("Child_ID","Gene","Effect")]
    desan[,1] <- sapply(desan[,1],function(s) {unlist(strsplit(s,"\\."))[1]})        
    desan[desan[,"Effect"] %in% lofsan, "MutaType"] <- "LOF"
    desan[desan[,"Effect"] %in% missan, "MutaType"] <- "MIS"
    desan[,"inChild"] <- 0
    desan <- desan[,c(1,5,3,2,4)]
    colnames(desan) <- c("quadId","inChild","effectType","effectGenes","MutaType")
    desan[,"From"] <- "Sanders et al."

    ## O' Roak et al., 189 trios, 
    deroa <- read.csv("ASD/nature10989/nature10989-s2.csv")
    lofroa <- c("Del_aa","frameshift","missense/nonsense","nonsense","splice")
    misroa <- c("missense")
    deroa <- deroa[deroa[,"Type"] %in% union(lofroa,misroa),]
    deroa <- deroa[deroa[,"Result"] %in% c("de novo","de novo, complex","de novo, weak"),]
    deroa <- deroa[,c("Person","Refseq_Symbol","Type")]
    deroa[,1] <- sapply(deroa[,1],function(s) {unlist(strsplit(s,"\\."))[1]})        
    deroa[deroa[,"Type"] %in% lofroa, "MutaType"] <- "LOF"
    deroa[deroa[,"Type"] %in% misroa, "MutaType"] <- "MIS"
    deroa[,"inChild"] <- 0
    deroa <- deroa[,c(1,5,3,2,4)]
    colnames(deroa) <- c("quadId","inChild","effectType","effectGenes","MutaType")
    deroa[,"From"] <- "O' Roak et al."

    ## nature 13772, 1445 trios: 2270-825 trios
    deset1 <- read.csv("ASD/nature13772/nature13772-s4.csv")
    lofset1 <- c("CODON_CHANGE_PLUS_CODON_DELETION","CODON_DELETION","frameshift","nonsense","splice","START_GAINED")
    misset1 <- c("missense")
    deset1 <- deset1[deset1[,"Category"] %in% union(lofset1,misset1),]
    deset1 <- deset1[,c("Child_ID","Gene_name","Category")]
    deset1[deset1[,"Category"] %in% lofset1, "MutaType"] <- "LOF"
    deset1[deset1[,"Category"] %in% misset1, "MutaType"] <- "MIS"
    deset1[,"inChild"] <- 0
    deset1 <- deset1[,c(1,5,3,2,4)]
    colnames(deset1) <- c("quadId","inChild","effectType","effectGenes","MutaType")
    deset1[,"From"] <- "nature 13772"

    ## nature 13908, 2508
    deset2 <- read.csv("ASD/nature13908/Supplementary Table 2.csv")
    #deset2 <- deset2[deset2[,1] %in% fam13908F,] ## 2517 trios
    inchild <- c("pF","pFsF","pFsM","pM","pMsF","pMsM")
    deset2 <- deset2[deset2[,"inChild"] %in% inchild,]
    lofset2 <- c("frame-shift","nonsense","splice-site") # LGD mutations: likely gene disrupting: nonsense, frameshift and splice site
    misset2 <- c("missense")
    deset2 <- deset2[deset2[,"effectType"] %in% union(lofset2,misset2),]
    deset2 <- deset2[,c("familyId","inChild","effectType","effectGene")]
    deset2[deset2[,"effectType"] %in% lofset2, "MutaType"] <- "LOF"
    deset2[deset2[,"effectType"] %in% misset2, "MutaType"] <- "MIS"
    colnames(deset2) <- c("quadId","inChild","effectType","effectGenes","MutaType")
    deset2[,"From"] <- "nature 13908"
    
    ## 2508 + 1445 = 3953 trios !!!!1751 + 1445 + 343 + 238 + 189 = 3966 trios
    data <- rbind(deset1,deset2)
    
    source("Network_analysis.R")
    data[,"effectGenes"] <- mapping_to(data[,"effectGenes"])
    write.csv(data,file="ASD/ASDmutationlists.csv",row.names=FALSE)

    mutrate <- read.csv("ASD/MutationRate.csv")
    mutrate[,1] <- mapping_to(mutrate[,1])
    write.csv(mutrate,file="ASD/MutationRatem.csv",row.names=FALSE)

    ## TADA input format and random samples
    samples <- unique(data[,"quadId"])
    n <- length(samples)
    ta <- Random_ASD(data,mutrate,samples)
    write.csv(ta,file="ASD/randset/ASDalldenovo.csv",row.names=FALSE)

    ## different sample size and exclude sample sets
    for(j in 2:9){
        for(i in 1:10){
                subsamples <- sample(samples,round(j*n/10))
                ta <- Random_ASD(data,mutrate,subsamples)
                write.csv(ta,file=paste("ASD/randset/ASDrand",j,"depart",i,".csv",sep=""),row.names=FALSE)
                
                subsamples <- setdiff(samples,subsamples)
                ta <- Random_ASD(data,mutrate,subsamples)
                write.csv(ta,file=paste("ASD/randset/ASDrand",j,"denrest",i,".csv",sep=""),row.names=FALSE)   
        }
    }


    ## different sample size and more additional samples
    for(j in 1:3){
        for(i in 1:5){
            subsamples <- sample(samples,round(j*n/10))
            ta <- Random_ASD(data,mutrate,subsamples)
            write.csv(ta,file=paste("ASD/randset2/ASD1sub",j,"depart",i,".csv",sep=""),row.names=FALSE)
        
            restsamples <- setdiff(samples,subsamples)
            subset2 <- sample(restsamples,round(2*j*n/10))
            subsets <- union(subset2,subsamples)
            ta <- Random_ASD(data,mutrate,subsets)
            write.csv(ta,file=paste("ASD/randset2/ASD2sub",j,"denrest",i,".csv",sep=""),row.names=FALSE)  
        }
    }

}

ASD_rand_set <- function(){
    
    source("ASD_data_set.R")
    mutrate0 <- read.csv("ASD/GeneMutaRatem.csv")
    mutrate1 <- read.csv("ASD/MutationRatem.csv")
    mutrate <- cbind(mutrate0,mutrate1[match(mutrate0[,1],mutrate1[,1]),"mut.rate"])
    colnames(mutrate)[6] <- "mut.rate"
    
    data <- read.csv("DDD_mutations/datasheet/anno_mutationlist_4_15.csv")
    data <- data[data[,"Disease"]=="ASD",]
    mutype <- "mutation_type_metaSVM"
    
    samples <- unique(data[,"sampleID"])
    ta <- Random_ASD_N(data,mutrate,samples,mutype)
    write.csv(ta,file="ASD/ASDLOF_dMIS_4_16.csv",row.names=FALSE)
    
    n <- length(samples)
    #### the thrid case: independent data sets
    for(j in 3){
        for(i in 1:20){
            subsamples <- sample(samples,round(j*n/10))
            ta <- Random_ASD_N(data,mutrate,subsamples,mutype)
            write.csv(ta,file=paste("ASD/randset4/ASD1sub",j,"depart",i,".csv",sep=""),row.names=FALSE)
            
            restsamples <- setdiff(samples,subsamples)
            ta <- Random_ASD_N(data,mutrate,restsamples,mutype)
            write.csv(ta,file=paste("ASD/randset4/ASD2sub",j,"derest",i,".csv",sep=""),row.names=FALSE)  
        }
    }
    
    filename="ASD/TADAresult/TADAdenovo_ASD4_16.csv"
    strname <- "ASD4_16"
    dirstr <- "result/"
    TADAinput1(filename,strname,dirstr,genelist0=1)
    
    filename="ASD/TADAresult/TADAdenovo_control.csv"
    strname <- "control_4_27"
    dirstr="result/control/"
    TADAinput1(filename,strname,dirstr,genelist0=1)
    
    source("Network_analysis.R")
    source("Multi_net.R")
    dirstr <- "result/randset4/"
    for(j in 3){
        for(i in 1:20){
            filename=paste("ASD/TADAresult/randset4/TADAdenovo_ASD1sub",j,"depart",i,".csv",sep="")
            strname <- paste("part",j,"_",i,sep="")
            TADAinput1(filename,strname,dirstr,genelist0=1)
        }
    }
    
    
    ### the second case: leave one mutation  ## changed by 8_5
    #filename <- "ASD/TADAresult/randset4/PtestASDall.csv"
    #filename <- "ASD/TADAresult/randset4/TADAdenovo_ASD4_16.csv"
    filename <- "ASD/TADAresult/TADAdenovo_ASD8_5.csv"
    asddata <- read.csv(filename)
    genes <- asddata[asddata[,"qvalue.dn"]<=0.1,"Gene"]
    
    for(i in 1:52){
        tmp <- asddata
        subs <- which(tmp[,"Gene"]==genes[i])
        if(length(subs) >1 ) print(i)
        if(tmp[subs,"dn.LoF"] > 0){
        #    tmp[subs,c("dn.LoF","dn.mis3","dn.mis")] <- c(1,0,0)  
            tmp[subs,c("dn.LoF","dn.mis3")] <- c(1,0)
        }else if(tmp[subs,"dn.LoF"] ==0 & tmp[subs,"dn.mis3"] > 0){
            #tmp[subs,c("dn.LoF","dn.mis3","dn.mis")] <- c(0,1,0)
            tmp[subs,c("dn.LoF","dn.mis3")] <- c(0,1)
        #}else if(tmp[subs,"dn.LoF"] ==0 & tmp[subs,"dn.mis3"] == 0 & tmp[subs,"dn.mis"] > 0){
        #    tmp[subs,c("dn.LoF","dn.mis3","dn.mis")] <- c(0,0,1)
        }
        write.csv(tmp,file=paste("ASD/leaveone4/rand2_",i,".csv",sep=""),row.names=FALSE)
    }
    
    source("Network_analysis.R")
    source("Multi_net.R")
    dirstr <- "result/leaveone4/"
    for(i in 1:52){
            filename=paste("ASD/TADAresult/leaveone4/TADAdenovo_rand2_",i,".csv",sep="")
            strname <- paste("rand2_",i,sep="")
            TADAinput1(filename,strname,dirstr,genelist0=1)
    }
    
    
    ## the first case: different sample size
    n <- length(samples)
    for(j in 2:9){
        for(i in 1:20){
            subsamples <- sample(samples,round(j*n/10))
            ta <- Random_ASD_N(data,mutrate,subsamples,mutype)
            write.csv(ta,file=paste("ASD/randset_1/ASDrand1_",j,"_",i,".csv",sep=""),row.names=FALSE) 
        }
    }
    
    ## run tada
    
    ## tranform to tada input files
    source("Network_analysis.R")
    source("Multi_net.R")
    dirstr <- "result/randset_1/"
    for(j in 2:9){
        for(i in 1:20){
            filename=paste("ASD/TADAresult/randset_1/TADAdenovo_ASDrand1_",j,"_",i,".csv",sep="")
            strname <- paste("rand1",j,"_",i,sep="")
            TADAinput1(filename,strname,dirstr,genelist0=1)
        }
    }
    
}

ASD_mis <- function(){
    
    ## nature 13772, 1445 trios: 2270-825 trios
    deset1 <- read.csv("ASD/nature13772/nature13772-s4.csv")
    lofset1 <- c("CODON_CHANGE_PLUS_CODON_DELETION","CODON_DELETION","frameshift","nonsense","splice","START_GAINED")
    misset1 <- c("missense")
    #deset1 <- deset1[deset1[,"Category"] %in% union(lofset1,misset1),]
    deset1[,"AAchange"] <- paste(deset1[,"Ref"],">",deset1[,"Alt"],sep="")
    deset1 <- deset1[,c("Chr","Pos","AAchange","Child_ID","Gene_name","Category")]
    deset1[deset1[,"Category"] %in% lofset1, "MutaType"] <- "LOF"
    deset1[deset1[,"Category"] %in% misset1, "MutaType"] <- "MIS"
    deset1[,"inChild"] <- 0
    deset1 <- deset1[,c(1,2,3,c(1,5,3,2,4)+3)]
    colnames(deset1) <- c("Chr","Pos","AAchange","quadId","inChild","effectType","effectGenes","MutaType")
    deset1[,"From"] <- "nature 13772"

    
    ## nature 13908, 
    deset2 <- read.csv("ASD/nature13908/Supplementary Table 2.csv")
    #deset2 <- deset2[deset2[,1] %in% fam13908F,] ## 2517 trios
    inchild <- c("pF","pFsF","pFsM","pM","pMsF","pMsM")
    deset2 <- deset2[deset2[,"inChild"] %in% inchild,]
    lofset2 <- c("frame-shift","nonsense","splice-site") # LGD mutations: likely gene disrupting: nonsense, frameshift and splice site
    misset2 <- c("missense","3'UTR","3'UTR-intron","5'UTR","5'UTR-intron","intron","missense","no-frame-shift","no-frame-shift-newStop","noEnd","non-coding","non-coding-intron","noStart")
    #deset2 <- deset2[deset2[,"effectType"] %in% union(lofset2,misset2),]
    deset2[,"Chr"] <- sapply(1:dim(deset2)[1],function(i) unlist(strsplit(deset2[i,"location"],":"))[1])
    deset2[,"Pos"] <- sapply(1:dim(deset2)[1],function(i) unlist(strsplit(deset2[i,"location"],":"))[2])
    deset2[,"AAchange"] <- deset2[,"variant"] 
    
    deset2 <- deset2[,c("Chr","Pos","AAchange","familyId","inChild","effectType","effectGene")]
    deset2[deset2[,"effectType"] %in% lofset2, "MutaType"] <- "LOF"
    deset2[deset2[,"effectType"] %in% misset2, "MutaType"] <- "MIS"
    colnames(deset2) <- c("Chr","Pos","AAchange","quadId","inChild","effectType","effectGenes","MutaType")
    deset2[,"From"] <- "nature 13908"
    
    ## 2508 + 1445 = 3953 trios !!!!1751 + 1445 + 343 + 238 + 189 = 3966 trios
    data <- rbind(deset1,deset2)
    
    source("Network_analysis.R")
    data[,"effectGenes"] <- mapping_to(data[,"effectGenes"])
    write.csv(data,file="ASD/ASDmutationlists_all.csv",row.names=FALSE)
    
    mutrate <- read.csv("ASD/MutationRatem.csv")
    
    ## TADA input format and random samples
    samples <- unique(data[,"quadId"])
    n <- length(samples)
    ta <- Random_ASD(data,mutrate,samples)
    write.csv(ta,file="ASD/randset/ASDalldenovo_Mis.csv",row.names=FALSE)
 
    ### wannovar
    deset1 <- read.csv("ASD/nature13772/nature13772-s4.csv")
    deset1a <- deset1[,c("Chr","Pos","Pos","Ref","Alt","Gene_name","Category","Child_ID")]
    colnames(deset1a) <- c("Chr","Start","End","Ref","Alt","Gene","Category","sampleID")
    tmp <- sapply(1:dim(deset1)[1], function(i){
        subtmp <- unlist(strsplit(deset1[i,26],";"))
        if(any(grepl("SNPEFF_CODON_CHANGE",subtmp))){
            a1 <- subtmp[grepl("SNPEFF_CODON_CHANGE",subtmp)]
            subs <- gsub("SNPEFF_CODON_CHANGE=","",a1)
        }else{
            subs=""
        }
        subs
        })
    deset1a[,"AAchanges"] <- tmp
    deset1a[,"Disease"] <- "ASD"
    deset1a[,"From"] <- "nature 13772"
    
    ## nature 13908, 1751 trios,
    deset2 <- read.csv("ASD/nature13908/Supplementary Table 2.csv",stringsAsFactors=FALSE)
    inchild <- c("pF","pFsF","pFsM","pM","pMsF","pMsM")
    deset2 <- deset2[deset2[,"inChild"] %in% inchild,]
    

    tmp <- sapply(1:dim(deset2)[1],function(i) unlist(strsplit(deset2[i,"vcfVariant"],":"))[1])
    tmp1 <- sapply(1:dim(deset2)[1],function(i) unlist(strsplit(deset2[i,"vcfVariant"],":"))[2])
    deset2a <- cbind(tmp,tmp1)
    deset2a <- cbind(deset2a,tmp1)
    deset2a <- as.data.frame(deset2a)
#     subs <- !grepl("sub",deset2[,"variant"])
#     deset2[subs,"variant"] <- "sub(-->-)"  ##!!!!!!
#     tmp <- gsub("sub\\(","",deset2[,"variant"])
#     tmp <- gsub("\\)","",tmp)
    deset2a[,"Ref"] <- sapply(1:dim(deset2)[1],function(i) unlist(strsplit(deset2[i,"vcfVariant"],":"))[3])
    deset2a[,"Alt"] <- sapply(1:dim(deset2)[1],function(i) unlist(strsplit(deset2[i,"vcfVariant"],":"))[4])
    deset2a <- cbind(deset2a,deset2[,"effectGene"])
    deset2a <- cbind(deset2a,deset2[,"effectType"])
    deset2a <- cbind(deset2a,deset2[,"familyId"])
    colnames(deset2a) <- c("Chr","Start","End","Ref","Alt","Gene","Category","sampleID")
    deset2a[,"AAchanges"] <- ""
    deset2a[,"Disease"] <- "ASD"
    deset2a[,"From"] <- "nature 13908"
    ## 2508 + 1445 = 3953 trios !!!!1751 + 1445 + 343 + 238 + 189 = 3966 trios
    data <- rbind(deset1a,deset2a)
    
    source("Network_analysis.R")
    data[,"Gene"] <- mapping_to(data[,"Gene"])
    write.csv(data,file="ASD/ASDmutationlists_annovar.csv",row.names=FALSE)
    
}

TADAinput_6_7 <- function(){
    
    source("Network_analysis.R")
    source("Multi_net.R")
    dirstr <- "ASD/TADAresult/"
    filenames <- c("ASD/TADAresult/TADA_resultsall.csv","ASD/TADAresult/TADA_results13772.csv","ASD/TADAresult/TADA_results13908.csv","ASD/TADAresult/TADA_results932.csv")
    strname <- c("all","13772","13908","932")
    TADAinput1(filenames,strname,dirstr)
    
    source("Network_analysis.R")
    source("Multi_net.R")
    dirstr <- "result/PCGC/"
    filename <- c("ASD/TADAresult/TADAdenovo_resultsPCGC.csv")
    strname <- "PCGC"
    TADAinput1(filename,strname,dirstr,genelist0=1)
    
    ### control
    source("Network_analysis.R")
    source("Multi_net.R")
    dirstr <- "data/PCGCdata/"
    filename <- c("data/PCGCdata/TADAdenovo_control.csv")
    strname <- "control"
    TADAinput1(filename,strname,dirstr,genelist0=1)

    ### case
    filename <- c("data/PCGCdata/TADAdenovo_case.csv")
    strname <- "case"
    TADAinput1(filename,strname,dirstr,genelist0=1)
    
    #### meta 
    filename <- c("result/Nat1377_meta/TADAdenovo_Meta_dmis.csv")
    strname <- "meta_dmis"
    TADAinput1(filename,strname,dirstr,genelist0=1)
    filename <- c("result/Nat1377_meta/TADAdenovo_Meta_mis.csv")
    strname <- "meta_mis"
    TADAinput1(filename,strname,dirstr,genelist0=1)
    
    
    ### ASD all
    source("Network_analysis.R")
    source("Multi_net.R")
    dirstr <- "result/"
    filename <- c("ASD/TADAresult/TADAdenovo_ASD4_16.csv")
    strname <- "ASD4_16"
    TADAinput1(filename,strname,dirstr,genelist0=1)


    
    
    ### Nat13772 all mutations and variants 
    filename <- c("ASD/TADAresult/TADA_results.csv")
    strname <- "nat13772_0"
    TADAinput1(filename,strname,dirstr,genelist0="")
    source("Network_analysis.R")
    a <- read.csv(filename)
    a[,1] <- mapping_to(a[,1])
    write.csv(a,file="ASD/TADAresult/TADA_resultsm.csv",row.names=FALSE)
    

    ## nature 13772 with only de novo mutations
    source("Network_analysis.R")
    source("Multi_net.R")
    dirstr <- "result/"
    filename <- c("ASD/TADAresult/TADAdenovo_nat13772.csv")
    strname <- "nat13772"
    TADAinput1(filename,strname,dirstr,genelist0=1)   
    
    
}

DAWNresult <- function(){
    filenames <- c("result/DAWN13772.csv","result/DAWN13908.csv","result/DAWN932.csv")
    strname <- c("13772","13908","932")
    for(i in 1:length(filenames)){
        result <- read.csv(filenames[i])
        result <- result[,c("Gene","original_pvalue","Stage1_posterior","Stage2_posterior","FDR")]
        write.table(result,file=paste("result/DAWNresult",strname[i],".txt",sep=""),quote=FALSE,col.names=FALSE,row.names=FALSE)
    }
    
    filenames <- c("result/DAWNrand1result.csv","result/DAWNrand2result.csv","result/DAWNrand3result.csv","result/DAWNallPCGCresult.csv","result/DAWNPCGCresult.csv")
    strname <- c("rand1","rand2","rand3","allPCGC","PCGC")
    for(i in 1:length(filenames)){
        result <- read.csv(filenames[i])
        result <- result[,c("Gene","pval.TADA","DAWN.Stage1_posterior","DAWN.Stage2_posterior","DAWN.FDR")]
        write.table(result,file=paste("result/DAWNresult",strname[i],".txt",sep=""),quote=FALSE,col.names=FALSE,row.names=FALSE)
    }
    
}

Random_ASD_N <- function(data,mutrate,samples,mutype="mutation_type"){
    
    data <- data[data[,"sampleID"] %in% samples,]
    counttable <- mutrate
    
    dlist <- list()
    types <- c("LOF","dMIS","MIS","SYNONMOUS")
    numbs <- c("dn.LoF","dn.mis3","dn.mis","dn.syn")
    
    for(i in 1:length(types)){
        dlist <- data[data[,mutype]==types[i],]
        counttable[,numbs[i]] <- 0
        
        tmp <- table(dlist[,"Gene"])
        genes <- intersect(names(tmp),counttable[,"Gene"])
        counttable[match(genes,counttable[,"Gene"]),numbs[i]] <- tmp[match(genes,names(tmp))]
    }
    
    counttable
}

Random_ASD <- function(data,mutrate,samples){
    
    data <- data[data[,"quadId"] %in% samples,]
    
    dlist <- list()
    dlist[[1]] <- data[data[,"MutaType"]=="LOF",]
    dlist[[2]] <- data[data[,"MutaType"]=="MIS",]
    counttable <- mutrate
    counttable[,"dn.LoF"] <- 0
    counttable[,"dn.mis3"] <- 0
 
    tmp <- table(dlist[[1]][,"effectGenes"])
    genes <- intersect(names(tmp),counttable[,"Gene"])
    counttable[match(genes,counttable[,"Gene"]),"dn.LoF"] <- tmp[match(genes,names(tmp))]
    
    tmp <- table(dlist[[2]][,"effectGenes"])
    genes <- intersect(names(tmp),counttable[,"Gene"])
    counttable[match(genes,counttable[,"Gene"]),"dn.mis3"] <- tmp[match(genes,names(tmp))]
    
    counttable
}

randset_TADAinput <- function(){
    source("Network_analysis.R")
    source("Multi_net.R")
    dirstr <- "result/randset3/"
    
    strname <- "ASDall"
    filename <- "ASD/TADAresult/TADAdenovo_ASDall.csv"
    TADAinput1(filename,strname,dirstr,genelist0=1)
    
    ## different sample size and exclude sample sets
    for(j in 2:9){
        for(i in 1:10){
            filename=paste("ASD/TADAresult/randset3/TADAdenovo_part",j,"_",i,".csv",sep="")
            strname <- paste("part",j,"_",i,sep="")
            TADAinput1(filename,strname,dirstr,genelist0=1)
            
#             filename=paste("ASD/TADAresult/TADAdenovo_rest",j,"_",i,".csv",sep="")
#             strname <- paste("rest",j,"_",i,sep="")
#             TADAinput1(filename,strname,dirstr,genelist0=1)
        }
    }

    dirstr <- "result/randset2/"
    for(j in 1:3){
        for(i in 1:5){
            filename=paste("ASD/TADAresult/randset2/TADAdenovo_ASD1sub",j,"depart",i,".csv",sep="")
            strname <- paste("part",j,"_",i,sep="")
            TADAinput1(filename,strname,dirstr,genelist0=1)
            
            filename=paste("ASD/TADAresult/randset2/TADAdenovo_ASD2sub",j,"denrest",i,".csv",sep="")
            strname <- paste("rest",j,"_",i,sep="")
            TADAinput1(filename,strname,dirstr,genelist0=1)
        }
    }
}

randset2_TADAinput <- function(){
    source("Network_analysis.R")
    source("Multi_net.R")
    dirstr <- "result/leaveone/"
   
    ## leave one gene data set with FDR <= 0.2
    for(i in 1:162){
        filename=paste("ASD/TADAresult/leaveone/TADAdenovo_rand2_",i,".csv",sep="")
        strname <- paste("rand2_",i,sep="")
        TADAinput1(filename,strname,dirstr,genelist0=1)
    }
    
    source("Network_analysis.R")
    source("Multi_net.R")
    dirstr <- "result/leaveone3/"
    ## leave one missense mutation data set with FDR <= 0.2
    for(i in 1:251){
        filename=paste("ASD/TADAresult/leaveone3/TADAdenovo_rand2_de1_",i,".csv",sep="")
        strname <- paste("rand2_de1_",i,sep="")
        TADAinput1(filename,strname,dirstr,genelist0=1)
    }
}

leaveone_ASD <- function(){

    TADAFile <- "ASD/TADAresult/TADAdenovo_ASDall.csv"
    TADAf <- read.csv(TADAFile)
    genes <- TADAf[TADAf[,"qvalue.dn"] <= 0.3,"Gene"]
    ASDalldenovo <- read.csv("ASD/randset/ASDalldenovo.csv")
    
    for(i in 1:length(genes)){
        tmp <- ASDalldenovo
        subs <- which(tmp[,"Gene"]==genes[i])
        if(tmp[subs,3] > 0){
            tmp[subs,3:4] <- c(1,0) ## leaveone1: c(0,1) ## leaveone: 0 ## leaveone2: c(1,0) ## leaveone3: c()
        }else{
            tmp[subs,3:4] <- c(0,1)
        }
        write.csv(tmp,file=paste("ASD/leaveone3/rand2_",i,".csv",sep=""),row.names=FALSE)
    }

}

mapping_DAWNexp <- function(){

    source("Network_analysis.R")
    module <- read.csv("module.csv")
    
    module[,1] <- mapping_to(module[,1])
    subs <- duplicated(module[,1])
    module <- module[!subs,]
    
    save(module,file="DAWNmodule")
    write.csv(module,file="DAWNmodule.csv",row.names=FALSE)
    
    load("geneexp.Rdata")
    genes <- rownames(data)
    genes <- mapping_to(genes)
    
    rownames(data) <- genes
    save(data,file="DAWNgeneexp.Rdata")
    
    
}

net_DAWNexp <- function(){
    
    options(stringsAsFactors=FALSE)
    library(WGCNA)
    source("Network_analysis.R")
    
    module <- read.csv("data/module.csv")
    load("geneexp.Rdata")
    data <- data[rownames(data) %in% module[,1],]
    Labels <- unique(module[,2])
    allnet<- c()
    for(i in Labels){
        modgenes <- module[module[,2]==i,1]
        modexp <- data[match(modgenes,rownames(data)),]
        corm=abs(cor(t(modexp),use='pair',nThreads=3))
        diag(corm) <- 0
        edges <- which(corm > 0.7,arr.ind=TRUE)
        net <- cbind(modgenes[edges[,1]],modgenes[edges[,2]])
        net <- cbind(net,corm[edges])
        allnet <- rbind(allnet,net)
        print(i)
    }
    
    genes <- union(allnet[,1],allnet[,2])
    genes1 <- mapping_to(genes)
    
    allnet[,1] <- genes1[match(allnet[,1],genes)]
    allnet[,2] <- genes1[match(allnet[,2],genes)]
    write.table(allnet,file="data/network_inference/ASDexp_net.txt",quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
  
    
    ###===================
    source("ASD_data_set.R")
    source("Network_analysis.R")
    module <- read.csv("data/module.csv")
    load("geneexp.Rdata")
    data <- data[rownames(data) %in% module[,1],]
    Labels <- unique(module[,2])
    allnet<- c()    
    for(i in Labels){
        modgenes <- module[module[,2]==i,1]
        modexp <- data[match(modgenes,rownames(data)),]
        corm=abs(cor(t(modexp),use='pair',nThreads=3))
        diag(corm) <- 0
        
        corm1 <- select_para(corm,d=5)$corm1
        edges <- which(corm1==1,arr.ind=TRUE)
        net <- cbind(modgenes[edges[,1]],modgenes[edges[,2]])
        net <- cbind(net,corm[edges])
        allnet <- rbind(allnet,net)
        print(i)
    }
    
    genes <- union(allnet[,1],allnet[,2])
    genes1 <- mapping_to(genes)
    
    allnet[,1] <- genes1[match(allnet[,1],genes)]
    allnet[,2] <- genes1[match(allnet[,2],genes)]
    write.table(allnet,file="data/network_inference/ASDexp_net_top5.txt",quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
}

select_para <- function(corm,d=5){
    n <- dim(corm)[1]
    corm1 <- matrix(0,n,n)
    x <- rep(0,n)
    for(i in 1:n){
        tmp <- corm[i,]
        x[i] <- sort(tmp,decreasing = TRUE)[d]
        corm1[i,tmp >= x[i]] <- 1
    }
    corm1 <- corm1 + t(corm1)
    corm1[corm1>1] <- 1
    
    list(corm1=corm1,x=x)
}

select_para0 <- function(corm,d=5){
    library(optimx)
    x <- rowMeans(corm)
    n <- dim(corm)[1]
    re <- optim(x, fn, gr=NULL, corm=corm,method="L-BFGS-B",lower=rep(0,n),upper=rep(1,n))
    
    x <- re$par
    X <- matrix(x,n,n)
    corm1 <- corm - X
    corm1[corm1<=0] <- 0
    corm1[corm1>0] <- 1
    corm1 <- corm1 + t(corm1)
    corm1[corm1 >1] <- 1

    list(corm1=corm1,x=x)
}

fn <- function(x,corm){
    library(PCIT)
    ## mean/sd for edge weights + for all node clustering coefficient
    n <- dim(corm)[1]
    X <- matrix(x,n,n)
    corm1 <- corm - X
    corm1[corm1<=0] <- 0
    corm1[corm1>0] <- 1
    corm1 <- corm1 + t(corm1)
    corm1[corm1 >1] <- 1
    
    cc <- localClusteringCoefficient(corm1)
    cc[is.na(cc)] <- 0
    f1 <- sum(cc)
    
    corm2 <- matrix(0,n,n)
    corm2[corm1==1] <- corm[corm1==1]
    f2 <- sapply(1:n,function(i) mean(corm2[i,corm2[i,]>0]))
    f2[is.na(f2)] <- 0
    f3 <- sapply(1:n,function(i) sd(corm2[i,corm2[i,]>0]))
    f3[f3==0] <- 1
    
    #f <- -log(sum(f2)) - log(f1/2)
    f <- -sum(cc*f2)
    
    f
}

plot_hist <- function(){
    a <- read.csv("ASD/ASDmutationlists.csv")
    pdf(file="plot/denovocounts_samples.pdf")
    b <- table(a[,1])
    bv <- 1:max(b)
    for(i in 1:max(b)){
        bv[i] <- sum(b==i)
    }
    barplot(bv,width=0.1,ylab="Frequency",xlab="De novo count",main="",names.arg=1:max(b),col="Red",cex.axis=2)
    axis(1,cex.axis=2)
    dev.off()
    
    pdf(file="plot/denovocounts_genes.pdf",height=8,width=14)
    par(mai=c(2,2,2,2))
    b <- table(a[,"effectGenes"])
    bv <- 1:max(b)
    for(i in 1:max(b)){
        bv[i] <- sum(b==i)
    }
    barplot(bv,width=0.1,space=0.2,ylab="The number of genes",xlab="De novo count",main="",names.arg=1:max(b),col="Red",cex.axis=1.5,xaxt="n",cex.main=1.8,cex.lab=2,cex.names=1.8)
    axis(1,at=seq(0.07,1.75,length.out=15),labels=1:max(b),cex.axis=1.5,cex.lab=2)
    dev.off()
}

build_module <- function(){
    source("Multi_net.R")
    source("CRF_build.R")
    library(WGCNA)
    dirstr <- "data/expressiondata/"
    filenames <- c("data/network_inference/OverlapPPI.txt","data/network_inference/Overlap_O4_top5.txt","data/network_inference/Overlap_O4_top5_PPI.txt","data/network_inference/Joint_O4_top5.txt","data/network_inference/Joint_O4_top5_PPI.txt")
    strname <- c("data/network_inference/OverlapPPImodule","data/network_inference/Overlap_O4_top5module","data/network_inference/Overlap_O4_top5_PPImodule","data/network_inference/Joint_O4_top5module","data/network_inference/Joint_O4_top5_PPImodule")
    Kv <- c(19,15,17,16,18)
   
    for(i in c(1,3,5)){
        filename <- filenames[i]
        distgene = build_net(Kv[i],"")$matrix
        softPower <- 1
        minModuleSize = 100; # set the minimum module size relatively high:
        adjacency = adjacency.fromSimilarity(distgene, type = "unsigned", power=softPower)
        TOM = TOMsimilarity(adjacency);
        dissTOM = 1-TOM
        geneTree = flashClust(as.dist(dissTOM), method = "average"); # Call the hierarchical clustering function
        dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize); # Module identification using dynamic tree cut
        dynamicColors = labels2colors(dynamicMods) # Convert numeric lables into colors
        moduleColors = dynamicColors
        colorOrder = c("grey", standardColors());
        Labels = match(moduleColors, colorOrder)-1;
        module <- cbind(rownames(distgene),Labels)
        
        module[module[,2]==0,2] <- max(as.numeric(module[,2])) + 1
        save(module,file=strname[i])
    }
    
    
}

control_data <- function(){
    source("Network_analysis.R")
    source("DDD.R")
    
     data <- read.csv("data/controls_denovo_info0505.csv")
     data <- data[,c("proband.ID.GT.AD.DP.GQ.PL.","GeneName","VarClass")]
     colnames(data) <- c("Family_ID","Gene","Type")

     lof <- c("frameshiftdeletion","frameshiftinsertion","stopgainSNV")
     mis <- "nonsynonymousSNV"
     data[data[,"Type"] %in% lof,"MutaType"] <- "LOF"
     data[data[,"Type"] %in% mis,"MutaType"] <- "MIS"
     data[,"Disease"] <- "Control"
     data[,"From"] <- "SSC900control"
     data <- data[,c(1,5,3,2,4,6)]
     data[,4] <- mapping_to(data[,4])
     
     mutrate <- read.csv("ASD/MutationRatem.csv")
     t <- Random_DDD(data,mutrate)
    write.csv(t,file="ASD/controlcount.csv",row.names=FALSE)
    
}

MAGI_r <- function(){

    pathdir <- "/Users/qh2159/Downloads/DAWN_package/result/randresult3/MAGI/"
    filenames <- c()
   
    k=1
    #filename <- paste(pathdir,"rset1.o2863169.",k,sep="")
    filename <- paste(pathdir,"RandomGeneList.",k,sep="")
    tmp <- read.table(filename)
   
 
    k <- k + 1
    for(j in 2:9){
        for(i in 1:10){
            filename <- paste(pathdir,"RandomGeneList.",k,sep="")
           
            k=k+1
        }
    }
    
   

}

one_MAGI_module <- function(filename){

    #tmp <- read.delim(filename,stringsAsFactors=FALSE,sep="\t",header=FALSE)
    tmp <- readLines(con <- file(filename,"r"))
    close(con)
    k=1
    
    Kmodule <- 0
    Nmodule <- 1:10000
    Mogenes <- c()
    m <- 1
    while(k<=length(tmp)){
        Nmodule[m] <- as.numeric(tmp[k])
        Mogenes <- c(Mogenes,tmp[(k+1):(k+Nmodule[m])])
        Kmodule <- Kmodule + 1
        k <- k + Nmodule[m] + 2
        m <- m+1
    }
    Nmodule <- Nmodule[1:(m-1)]
    
    
}