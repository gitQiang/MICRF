Network_inference <- function(){

    ## References: a multi-method approach for proteomic network inference in 11 human cancers.
    ## Top 6 methods: RIDGENET, ARACNE-M, ARACNE-A, LASSONET, CLR AND SEPEARMANCOR
    library(stats) #cor(x, y = NULL, method = "spearman")
    library(parcor) #ridge.net, adalasso.net
    library(parmigene) #aracne.a, aracne.m,clr
    options(stringsAsFactors=FALSE)
    
    module <- read.csv("module.csv")
    load("geneexp.Rdata")
    dexp <- data
    labels <- unique(module[,2])
    k <- 1
    genes <- rownames(dexp)
    allnet <- c()
    for(i in labels){
        modgenes <- module[module[,2]==i,1]
        modexp <- dexp[genes %in% modgenes,]
        M <- list() ## edge weights
        R <- list() ## edge ranks 
        ## cor spearman
        M[[1]] <- abs(cor(t(modexp), method = "spearman"))
        ## ridgenet
        M[[2]] <- abs(ridge.net(t(modexp))$pcor)
        ## adalasso.net
        M[[3]] <- abs(adalasso.net(t(modexp),both=FALSE)$pcor.lasso)
        ## aracne.a
        mi  <- knnmi.all(modexp)
        M[[4]] <- aracne.a(mi)
        ## aracne.m
        M[[5]] <- aracne.m(mi)
        ## clr
        M[[6]] <- clr(mi)
        # caluated the shared percentage for six methods,choose parameter 
        for(j in 1:6){
            R[[j]] <- trans_rank(M[[j]])
        }
        edges <- select_para(R)
        net <- cbind(modgenes[edges[,1]],modgenes[edges[,2]])
        allnet <- rbind(allnet,net)
        print(i)
    }
    
    allnet <- cbind(allnet,1)
    write.table(allnet,file="Infer_net.txt",quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
}

trans_rank <- function(M,n=20000){
    
    M[lower.tri(M,diag = TRUE)] <- 0
    edges <- which(M>0,arr.ind=TRUE)
    edges <- cbind(edges,M[edges])
    edges <- edges[order(edges[,3],decreasing=TRUE),]
    if(dim(edges)[1] > n){
        edges <- edges[1:n,]
    }
    M[upper.tri(M)] <- 0
    M[edges[,1:2]] <- 1:dim(edges)[1]
    
    M
}

select_para <- function(R,n=20000){
    y <- 20000 ## the maximum number of edges
    for(j in 1:6){
        y <- min(y,sum(R[[j]]>0))
    }   
    y <- min(y,n)
    ymin=dim(R[[1]])[1] - 1 ## the minimum number of edges
    
    for(i in 1:6){
        R[[i]][R[[i]]==0] <- Inf
    }
    
    p  <- rep(0,y-ymin+1)# fraction of shared edges
    for(i in ymin:y){
        M <- (R[[1]] <= i) + (R[[2]] <= i) + (R[[3]] <= i) + (R[[4]] <= i) + (R[[5]] <= i) + (R[[6]] <= i)
        p[i-ymin+1] <- sum(M>=2)/i
    }

    cutn <- which(p==min(p)) + ymin - 1
    plot(p)
    i <- cutn
    M <- (R[[1]] <= i) + (R[[2]] <= i) + (R[[3]] <= i) + (R[[4]] <= i) + (R[[5]] <= i) + (R[[6]] <= i)
    edges <- which(M>=2,arr.ind=TRUE)
    
    edges
}

infer_netmap <- function(){
    
    source("Network_analysis.R")
    net0 <- read.table("data/network_inference/Infer_net_O4.txt",sep="\t",header=FALSE)
    net0[,1] <- mapping_to(net0[,1])
    net0[,2] <- mapping_to(net0[,2])
    
    write.table(net0,file="data/network_inference/Infer_net_O4m.txt",quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")

}

