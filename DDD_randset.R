DDD_randset <- function(){
    source("ASD_data_set.R")
    mutrate<- read.csv("data/mutarate/PCGCSCI.csv")
    mutrate <- mutrate[!duplicated(mutrate[,2]),-1]
    colnames(mutrate) <- c("Gene","LOF","dmis","mis","syn")
    
    ddd <- read.csv("DDD_mutations/datasheet/anno_DDDmutationlist_4_15.csv")
    asd <- read.csv("DDD_mutations/datasheet/anno_mutationlist_4_15.csv")
   
    data <- rbind(ddd,asd)
    mutype <- "mutation_type_metaSVM"
    
    samples <- unique(data[,"sampleID"])
    ta <- Random_ASD_N(data,mutrate,samples,mutype)
    write.csv(ta,file="DDD_mutations/datasheet/MutaTable5542.csv",row.names=FALSE)
    
    ## after running TADA
    source("Network_analysis.R")
    source("Multi_net.R")
    filename="../TADA_DAWN/result/TADAdenovo_meta_dmis.csv"
    strname <- "meta"
    dirstr <- "result/"
    TADAinput1(filename,strname,dirstr,genelist0=1)
    TADAinput2(filename,"meta2",dirstr,genelist0=1)
    
    n <- length(samples)
    #### the thrid case: independent data sets
    for(j in 3){
        for(i in 1:20){
            subsamples <- sample(samples,round(j*n/10))
            ta <- Random_ASD_N(data,mutrate,subsamples,mutype)
            write.csv(ta,file=paste("DDD_mutations/randset4/ASD1sub",j,"depart",i,".csv",sep=""),row.names=FALSE)
            
            restsamples <- setdiff(samples,subsamples)
            ta <- Random_ASD_N(data,mutrate,restsamples,mutype)
            write.csv(ta,file=paste("DDD_mutations/randset4/ASD2sub",j,"derest",i,".csv",sep=""),row.names=FALSE)  
        }
    }
    
    ### the second case: leave one mutation  ## changed by 8_12
    filename <- "../TADA_DAWN/result/TADAdenovo_meta_dmis.csv"
    asddata <- read.csv(filename)
    asddata[is.na(asddata[,"qvalue.dn"]),"qvalue.dn"] <- 1
    genes <- asddata[asddata[,"qvalue.dn"] < 0.1,"Gene"]
    
    for(i in 1:length(genes)){
        tmp <- asddata
        subs <- which(tmp[,"Gene"]==genes[i])
        if(length(subs) >1 ) print(i)
        if(tmp[subs,"dn.LoF"] > 0){
            tmp[subs,c("dn.LoF","dn.mis3")] <- c(1,0)
        }else if(tmp[subs,"dn.LoF"] ==0 & tmp[subs,"dn.mis3"] > 0){
            tmp[subs,c("dn.LoF","dn.mis3")] <- c(0,1)
        }
        write.csv(tmp,file=paste("DDD_mutations/leaveone4/rand2_",i,".csv",sep=""),row.names=FALSE)
    }
    
    ## the first case: different sample size
    n <- length(samples)
    for(j in 2:9){
        for(i in 1:20){
            subsamples <- sample(samples,round(j*n/10))
            ta <- Random_ASD_N(data,mutrate,subsamples,mutype)
            write.csv(ta,file=paste("DDD_mutations/randset_1/ASDrand1_",j,"_",i,".csv",sep=""),row.names=FALSE) 
        }
    }
        
}

TADAinputT <- function(){
    ##### after running TADA: randsetDDD.R
    source("Network_analysis.R")
    source("Multi_net.R")
    source("DDD_randset.R")
    #dirstr <- "result/randset4/"
    dirstr <- "result/randset4_2/"
    for(j in 3){
        for(i in 1:20){
            filename=paste("../TADA_DAWN/result/randset4/TADAdenovo_ASD1sub",j,"depart",i,".csv",sep="")
            strname <- paste("part",j,"_",i,sep="")
            #TADAinput1(filename,strname,dirstr,genelist0=1)
            TADAinput2(filename,strname,dirstr,genelist0=1)
        }
    }
    
    ## after running TADA: randsetDDD.R
    #dirstr <- "result/leaveone4/"
    dirstr <- "result/leaveone4_2/"
    filename <- "../TADA_DAWN/result/TADAdenovo_meta_dmis.csv"
    asddata <- read.csv(filename)
    asddata[is.na(asddata[,"qvalue.dn"]),"qvalue.dn"] <- 1
    genes <- asddata[asddata[,"qvalue.dn"] < 0.1,"Gene"]
    
    for(i in 1:length(genes)){
        filename=paste("../TADA_DAWN/result/leaveone4/TADAdenovo_rand2_",i,".csv",sep="")
        strname <- paste("rand2_",i,sep="")
        #TADAinput1(filename,strname,dirstr,genelist0=1)
        TADAinput2(filename,strname,dirstr,genelist0=1)
    }
    
    ## after running TADA: randsetDDD.R
    #dirstr <- "result/randset_1/"
    dirstr <- "result/randset_1_2/"
    for(j in 2:9){
        for(i in 1:20){
            filename=paste("../TADA_DAWN/result/randset_1/TADAdenovo_ASDrand1_",j,"_",i,".csv",sep="")
            strname <- paste("rand1",j,"_",i,sep="")
            #TADAinput1(filename,strname,dirstr,genelist0=1)
            TADAinput2(filename,strname,dirstr,genelist0=1)
        }
    }
    
}

control_data <- function(){

    source("ASD_data_set.R")
    mutrate<- read.csv("data/mutarate/PCGCSCI.csv")
    mutrate <- mutrate[!duplicated(mutrate[,2]),-1]
    colnames(mutrate) <- c("Gene","LOF","dmis","mis","syn")
    
    data <- read.csv("DDD_mutations/datasheet/anno_mutationlist_control.csv")
    mutype <- "mutation_type_metaSVM"
    
    samples <- unique(data[,"sampleID"])
    ta <- Random_ASD_N(data,mutrate,samples,mutype)
    write.csv(ta,file="DDD_mutations/datasheet/control_8_16.csv",row.names=FALSE)
    
    ## after running TADA
    source("Network_analysis.R")
    source("Multi_net.R")
    dirstr <- "result/control/"
    filename <- "../TADA_DAWN/result/TADAdenovo_control.csv"
    strname <- "control"
    TADAinput1(filename,strname,dirstr,genelist0=1)

    ### Nat 13772 TADA input
    ### Nat13772 all mutations and variants 
    source("Multi_net.R")
    dirstr <- "result/"
    
    filename <- c("../TADA_DAWN/result/TADAdenovo_nat8_20.csv")
    strname <- "nat13772_0"
    TADAinput1(filename,strname,dirstr,genelist0="")

    ## nature 13772 with only de novo mutations    
    strname <- "nat13772"
    TADAinput1(filename,strname,dirstr,genelist0=1)  
    
    ###### repeat DDD data set
    source("ASD_data_set.R")
    mutrate<- read.csv("data/mutarate/PCGCSCI.csv")
    mutrate <- mutrate[!duplicated(mutrate[,2]),-1]
    colnames(mutrate) <- c("Gene","LOF","dmis","mis","syn")
    
    ddd <- read.table("DDD_mutations/datasheet/rep-NDD.counts-and-Poisson.txt",header=TRUE)[,1:4]
    colnames(ddd) <- c("Gene","dn.mis3","dn.LoF","dn.mis")
    
    data <- mutrate
    ig <- intersect(mutrate[,1],ddd[,1])
    data[,c("dn.mis3","dn.LoF","dn.mis")] <- 0
    data[match(ig,data[,1]),c("dn.mis3","dn.LoF","dn.mis")] <- ddd[match(ig,ddd[,1]),c("dn.mis3","dn.LoF","dn.mis")]
    write.csv(data,file="DDD_mutations/datasheet/repDDD1984.csv",row.names=FALSE)
    
    ### after running TADA
    source("Multi_net.R")
    dirstr <- "result/"
    
    filename <- c("../TADA_DAWN/result/TADAdenovo_repDDD_8_24.csv")
    strname <- "repDDD"
    TADAinput1(filename,strname,dirstr,genelist0=1)  
}

TADAinput2 <- function(filenames,strname,dirstr="random_samples/",pi0=0.94,genelist0=""){
    
    if(genelist0==""){
        bfs <- "BF"
    }else{
        bfs <- "BF.dn"
    }    
    mapT <- as.matrix(read.delim("Fmap0121.txt",header=FALSE,sep="\t"))
    for(i in 1:length(filenames)){
        TADAFile <- filenames[i]
        geneinfo <- read.csv(TADAFile,stringsAsFactors=FALSE)
        geneinfo <- geneinfo[!is.na(geneinfo[,bfs]),]
        
        posP <- cbind(log10(geneinfo[,bfs]),-log10(geneinfo[,bfs])) # for matlab trainging and R inferring 
        
        infofile <- paste(dirstr,"hotnet_input",strname[i],".txt",sep="")
        write.table(cbind(geneinfo[,1],posP),file=infofile,quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")
        
        genelist <- geneinfo[,1]
        genelist <- mapT[match(genelist,mapT[,2]),1]
        nodesim <- posP
        rownames(nodesim) <- genelist
        subs <- !is.na(genelist)
        nodesim <- nodesim[subs,]
        infofile <- paste(dirstr,"CRF_input",strname[i],".txt",sep="")
        write.table(cbind(rownames(nodesim),nodesim),file=infofile,quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")
    }
    
}
