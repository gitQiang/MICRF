source("CRF_build.R")
source("Network_analysis.R")
source("Multi_net.R")
options(stringsAsFactors=FALSE)
fileexp <- "PCGCall.txt" 

netstr <- c("STRING/","iRef/","Infer/","coexp/","GNAT/")
k = 4
modulefiles <- c("module3","iRefmodule","DAWNmodule","DAWNmodule","data/GNATnet/module_brain")

netflag=7
    
if(netflag==3){ instr <- "CRF_input";idm=TRUE; }else{ instr <- "hotnet_input"; idm=FALSE;}
modulefile <- modulefiles[k]

dirstr <- "result/"
strname <- "ASD4_16"
filename <- "ASD/TADAresult/CRF_inputMeta.txt"
strn <- "result/meta"
cutn=20000
idm=FALSE;
beta=0.2;
net_e=4

allnet <- build_net(netflag,fileexp,net_e)
load(modulefile)
module <- module[module[,1] %in% allnet$node,]

nodeinfo <- build_node(module,1,filename)

library(CRF)
library(Corbi)

modLab <- setdiff(unique(module[,2]),0)
n.module <- length(modLab)
n.genes <- length(module[module[,2]!=0,1])

POSpro <- matrix(0,n.genes,3)
POSpro[,1] <- as.vector(module[module[,2]!=0,1]) # gene, probability; FDR, node score, edge score
for (i in 1:n.module){
    modgenes <- as.vector(module[module[,2]==modLab[i],1])
    if(length(netflag)==1){
        net <- discrete_net(allnet,modgenes,netflag,cutn,nodeinfo,net_e)
    }else{
        net <- combine_net(allnet,modgenes,netflag)
    }
    if(sum(net$matrix)/2==0){break;}
    print(sum(net$matrix)/2)
    
    modnodeinfo <- nodeinfo[match(modgenes,rownames(nodeinfo)),]
    ### change 6_3
    model <- build_model_0(modnodeinfo,net,net_e,beta)
    
    crfresult <- solve_crf(model, query.type=4)
    #result <- decode.lbp(model)
    POSpro[match(modgenes,POSpro[,1]),2] <- crfresult
}

write.table(POSpro,file=paste(strn,"LBP_",netflag,".txt",sep=""),quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")

