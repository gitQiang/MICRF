Bayes_build <- function(nodesim,module,datExpr,distgene){
    library(bnlearn)
    library(gRain)
    modLab <- as.numeric(setdiff(unique(module[,2]),0))
    n.module <- length(modLab)

    POSpro <- module[as.numeric(module[,2])!=0,] # gene, probability; FDR, node score, edge score
    POSpro[,2] <- 0
    for (i in 1:n.module){
        print(i)
        modgenes <- module[as.numeric(module[,2])==modLab[i],1]
        modExpr <- datExpr[,modgenes]
        modExpr <- as.data.frame(modExpr)
        modnodesim <- nodesim[match(modgenes,rownames(nodesim)),]
        modnodesim[is.na(modnodesim[,1]),] <- 0.5 
        rownames(modnodesim) <- modgenes
        modnodesim <- modnodesim[,1]
        
        ## == build Bayesian network
        #fitted <- build_net(modExpr)
        fitted <- use_net_hq(modExpr)
        
        ## == build conditional probability tables
        ## node feature: genes * num_features
        ## edge feature: one edge feature
        ## CPTable:  1/Z * exp(\sum node feature /num_features + edge feature)
        cplist <- list()
        n.node <- length(fitted)
        for(nn in 1:n.node){
            node <- fitted[[nn]]$node
            parents <- fitted[[nn]]$parents
            indeg <- length(parents)
            children <- fitted[[nn]]$children
            outdeg <- length(children)
            thet <- ifelse(indeg+outdeg==0,0.5,1/(indeg+outdeg))  # the parameter theta for edge feature (state transmission)
            cptone <- cptable_hq(node,parents,children,modnodesim,thet,nstate=2)
            cplist[[nn]] <- cptone
        }
        cplist <- compileCPT(cplist)
        
        ## == inference
        net <-  grain(cplist)
        #jtree = compile(net)
        for(nn in 1:n.node){
            POSpro[POSpro[,1]==fitted[[nn]]$node,2] <- querygrain(net, nodes = fitted[[nn]]$node)[[1]][1]
        }
        ## all probability and FDR cpmpute
    }
    
    POSpro <- POSpro[order(as.numeric(POSpro[,2]),decreasing = TRUE),]
    write.table(POSpro,file="Bayes_Posterior.txt",quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")
    
    mapT <- as.matrix(read.delim("Fmap.txt",header=FALSE,sep="\t"))
    pros1 <- mapT[match(POSpro[,1],mapT[,1]),2]
    pros1 <- cbind(pros1,POSpro[,2])
    write.table(pros1,file="Bayes_Posterior_1.txt",quote=FALSE,col.names=FALSE,row.names=FALSE,sep="\t")
    
    #POSpro
}

use_net_hq <- function(modExpr){
    modgenes <- colnames(modExpr)
    
    GRNFile <- "GRN_net_map.txt"
    GRNTable <- as.matrix(read.delim(GRNFile,sep="\t",header=FALSE))
    GRNTable <- GRNTable[GRNTable[,5]=="all",]
    mapT <- as.matrix(read.delim("Fmap.txt",header=FALSE,sep="\t"))
    GRNTable[,1] <- mapT[match(GRNTable[,1],mapT[,2]),1]
    GRNTable[,2] <- mapT[match(GRNTable[,2],mapT[,2]),1]
    GRNTable <- GRNTable[!is.na(GRNTable[,1]) & !is.na(GRNTable[,2]),]
    subs <- GRNTable[,1] %in% modgenes & GRNTable[,2] %in% modgenes
    subnet <- GRNTable[subs,]
    
    # delete the undirected edges or add the direction
#     nodes <- union(subnet[,1],subnet[,2])
#     deg <- rep(0,length(nodes))
#     names(deg) <- nodes
#     deg[names(table(subnet[,2]))] <- table(subnet[,2]) # only the TF
    subs <- sapply(1:dim(subnet)[1],function(i){
        a <- subnet[i,1]
        b <- subnet[i,2]
        b.a <- which(subnet[,1] %in% b & subnet[,2] %in% a)
        if(length(b.a)>0) {FALSE; }else{TRUE;}
        })
    
    subnet1 <- subnet[subs,]
    node1 <- union(subnet1[,1],subnet1[,2])
    # format to fitted
    fitted <- list()
    for(i in 1:length(node1)){
        fitted[[i]] <- list()
        fitted[[i]]$node <- node1[i]
        fitted[[i]]$parents <- unique(subnet1[subnet1[,1]==node1[i],2])
        fitted[[i]]$children <- unique(subnet1[subnet1[,2]==node1[i],1])
    }

    fitted
}

build_net_hq <- function(modExpr,distgene){
    modgenes <- colnames(modExpr)
    #moddist <- distgene[modgenes,modgenes]
    #cutf <- quantile(moddist,0.95)
    mim <- build.mim(modExpr,estimator="spearman")
    net2 <- minet( modExpr, estimator="pearson" ,method="aracne")
    adjw <- aracne(mim,eps=0)
    adj <- adjw
    adj[adj<=0] <- 0
    adj[adj>0] <- 1
    diag(adj) <- 0
    e = empty.graph(modgenes)
    amat(e) = adj
    net <- model2network(modelstring(e))
    dag = cextend(net)
    
}

build_net <- function(modExpr){

    #net <- rsmax2(modExpr, restrict = "si.hiton.pc", maximize = "hc", test="mi-g",alpha = 0.01, score = "bic-g")
    net <- hc(modExpr)
#     dsachs = discretize(modExpr, method = "hartemink", breaks = 2, ibreaks = 2, idisc = "quantile")
#     dnet = deal::network(dsachs[, bnlearn::node.ordering(net)])
#     dnet = deal::as.network(bnlearn::modelstring(net), dnet)
#     net = bnlearn::model2network(deal::modelstring(dnet))
    fitted = bn.fit(net, modExpr, method = "mle")
    fitted

}

read_net <- function(net.text){
    
    net.node <- unique(union(net.text[,1],net.text[,2]))
    net.node <- net.node[net.node != ""]
    net.size <- length(net.node)
    net.edge <- cbind(as.character(net.text[,1]), as.character(net.text[,2]))
    net.edge <- net.edge[net.edge[,2] != "", ]
    net.edge <- net.edge[net.edge[,1] != "", ]
    net.matrix <- matrix(0, net.size, net.size, dimnames=list(net.node, net.node))
    net.matrix[net.edge] <- as.numeric(net.text[,3])
    list(size=net.size, node=net.node, matrix=net.matrix)
    
}

cptable_hq <- function(node,parents, children, modnodesim,thet,nstate=2){
    library(gRain)
    yn <- c("1","2")
    
    n.ef <- 1   
    edgeM <- matrix(thet,nstate,nstate)
    edgeM[1,2] <- 1 - thet
    edgeM[2,1] <- 0
    edgeM[2,2] <- 1
    
    modnodesim <- as.matrix(modnodesim)
    n.nf <- dim(modnodesim)[2]
    if(n.nf ==1 ){
        n.pot <- modnodesim
    }else{
        n.pot <- rowSums(modnodesim)/n.nf
    }
    n.pot <- cbind(n.pot,1-n.pot) #!!!
    n.pa <- length(parents)
    
    if(n.pa>0){
        m <- 0 : (nstate^n.pa -1 )
        stateM <- as.binary(m)
        stateM <- stateM + 1
      
        cpt <- matrix(0,nstate,max(m)+1,)
        for(i in 1:dim(cpt)[2]){
            cpt[1,i] <- n.pot[node,1] + sum(n.pot[parents,stateM[,i]] * edgeM[stateM[,i],1])/n.pa
            cpt[2,i] <- n.pot[node,2] + sum(n.pot[parents,stateM[,i]] * edgeM[stateM[,i],2])/n.pa
        }
        #cpt <- exp(cpt)
        cpt <- cpt/sum(cpt)
        #print(cpt)
        cptarr <- cptable(c(node,parents),values=as.vector(cpt),levels=yn)
    }else{
        cpt <- n.pot[node,]
        cptarr <- cptable(node, values=cpt,levels=yn)
    }
    cptarr
    
}

as.binary <- function (x) { 
    #http://196sigma.wordpress.com/2012/02/09/converting-decimal-to-binary-in-r/
    base = 2 
    ndigits = 1+ floor(log(max(x), base)) 
    m = length(x)
    r = matrix(0,ndigits,m)
    for (i in 1:ndigits) { 
        r[i, ] <- x%%base 
        x <- x%/%base 
    } 
    #class(r) <- "binaryInt" 
    #attr(r, "base") <- base 
    return(r) 
}