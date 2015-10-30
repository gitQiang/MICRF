source("Network_analysis.R")
TADAFile <- "TADA_lofmis1202.csv"
geneinfo <- as.matrix(read.csv(TADAFile))
tadagenes <- geneinfo[,1]
n.g <- length(tadagenes)
tadagenes <- mapping_to(tadagenes)
tadagenes294 <- tadagenes[as.numeric(geneinfo[,"pval.TADA"]) < 0.01]
Ourgenes <- as.matrix(read.table("Pheno_score.txt"))
Ourgenes294 <- Ourgenes[1:294]
length(intersect(Ourgenes294,tadagenes294))

genes <- intersect(Ourgenes294,tadagenes294)
x <- match(genes,Ourgenes)
y <- match(genes,tadagenes)
wilcox.test(x,y)
t.test(x,y)

genes <- intersect(Ourgenes,tadagenes)
x <- match(genes,Ourgenes)
y <- match(genes,tadagenes)
wilcox.test(x,y)
flag <- 2
if(flag ==1){
    # NDD topgenes
    nddfile <- "NDD6topgenes.txt"
    nddgene <- as.matrix(read.table(nddfile))
    # non-NDD topgenes
    nonddfile <- "nonNDD6topgenes.txt"
    nonddgene <- as.matrix(read.table(nonddfile))
}else if(flag==2){
    k = 200
    nddfile <- "NDD6LBPresult.txt"
    nddgene <- as.matrix(read.table(nddfile))[1:k,1]
    nonddfile <- "nonNDD6LBPresult.txt"
    nonddgene <- as.matrix(read.table(nonddfile))[1:k,1]
}

###  section 3.1
tadafin <- read.csv("TADA_resultsPCGC.csv")
rASDgene <- as.vector(tadafin[tadafin[,"qvalue"] < 0.5,1])
Our <- as.matrix(read.table("NewPheno_score.txt"))
Topgenes <- Our[1:length(rASDgene),1]
length(intersect(Topgenes,rASDgene))
genes <- intersect(Topgenes,rASDgene)
x <- match(genes,Topgenes)
y <- match(genes,rASDgene)
wilcox.test(x,y)
t.test(x,y)

## FDR < 0.5, compare with DAWN
dawnfin <- read.csv("PCGCall.csv")
rASDgene <- as.vector(dawnfin[,1])
fileflag <- 1
testname <- "DAWNFinal.txt"
mapT <- as.matrix(read.delim("Fmap.txt",header=FALSE,sep="\t"))
result <- read.csv("DAWN2result_4_1_0.65_0.5_0.5_1_12_2.csv")
result[,1] <- result[,3]
flag <- length(rASDgene)

dawnfin <- read.csv("PCGCresult.csv")
#dawnfin[is.na(dawnfin[,"DAWN.FDR"]),"DAWN.FDR"] <- 1
#rASDgene <- as.vector(dawnfin[as.numeric(dawnfin[,"DAWN.FDR"])<0.5,1])
dawnfin[is.na(dawnfin[,"DAWN.Stage2_posterior"]),"DAWN.Stage2_posterior"] <- 1
rASDgene <- as.vector(dawnfin[as.numeric(dawnfin[,"DAWN.Stage2_posterior"])<0.5,1])
fileflag <- 1
testname <- "DAWNFinal0126.txt"
mapT <- as.matrix(read.delim("Fmap.txt",header=FALSE,sep="\t"))
result <- read.csv("PCGCresult.csv")
flag <- length(rASDgene)

tadafin <- read.csv("TADA_resultsPCGC.csv")
tadafin[is.na(tadafin[,"qvalue"]),"qvalue"] <- 1
rASDgene <- as.vector(tadafin[as.numeric(tadafin[,"qvalue"])<0.4,1])
fileflag <- 1
testname <- "tadaFinal0126.txt"
mapT <- as.matrix(read.delim("Fmap.txt",header=FALSE,sep="\t"))
result <- tadafin

resultS <- as.matrix(read.table("NewPheno_score.txt"))
resultS[is.na(resultS[,2]),2] <- 1
rASDgene <- as.vector(resultS[as.numeric(resultS[,2])>0.7,1])
fileflag <- 1
testname <- "ScoreFinal0126.txt"
mapT <- as.matrix(read.delim("Fmap.txt",header=FALSE,sep="\t"))
result <- resultS

resultO <- as.matrix(read.table("result/Newscore_6_resultall.txt"))
rASDgene <- as.vector(resultO[as.numeric(resultO[,6]) > 0.7,1])
fileflag <- 1
testname <- "result/OurscoreFinal0126.txt"
mapT <- as.matrix(read.delim("Fmap.txt",header=FALSE,sep="\t"))
result <- resultO

resultS <- as.matrix(read.table("GammaPheno_score.txt"))
resultS[is.na(resultS[,2]),2] <- 1
rASDgene <- as.vector(resultS[as.numeric(resultS[,2])>0.9,1])
fileflag <- 1
testname <- "GammaFinal0126.txt"
mapT <- as.matrix(read.delim("Fmap.txt",header=FALSE,sep="\t"))
result <- resultS

resultO <- as.matrix(read.table("result/Newgamma_6_resultall.txt"))
rASDgene <- as.vector(resultO[as.numeric(resultO[,6]) > 0.9,1])
fileflag <- 1
testname <- "result/OursgammaFinal0126.txt"
mapT <- as.matrix(read.delim("Fmap.txt",header=FALSE,sep="\t"))
result <- resultO








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


mapT <- as.matrix(read.delim("Fmap0121.txt",header=FALSE,sep="\t"))
TADAFile <- "TADA_resultsPCGC.csv"
geneinfo <- as.matrix(read.csv(TADAFile))
genelist <- geneinfo[,1]
genelist <- mapping_to(genelist)
genelist <- mapT[match(genelist,mapT[,2]),1]
nodesim <- as.matrix(as.numeric(geneinfo[,"BF"])/(as.numeric(geneinfo[,"BF"]) + 1))
#nodesim <- as.matrix(1- as.numeric(geneinfo[,"pval.TADA"]))
names(nodesim) <- genelist
subs <- !is.na(genelist)
nodesim <- nodesim[subs]
write.table(cbind(names(nodesim),nodesim),file="TADAinfo0126.txt",quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")





