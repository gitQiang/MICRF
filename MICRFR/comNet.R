comNet <- function(){
    source("misc.R")
    
    ## combined the two best performance networks: iRefIndex and co-expression networks
    filename <- "data/hotnet/iRefIndexm.txt";
    netiRef <- net_matrix(filename)
    netiRef$matrix[netiRef$matrix > 0] <- 1
    
    filename <- "data/network_inference/brainspan_net_top5.txt"; ## brain coexp 5
    netCo <- net_matrix(filename)
    netCo$matrix[netCo$matrix > 0] <- 1
    
    ## read the edge betweenness 
    filename <- "data/Network_betweenness/Betweenness_edge_iRef.txt"
    netiRef <- net_between(filename,netiRef)
    filename <- "data/Network_betweenness/Betweenness_edge_coexp5.txt"
    netCo <- net_between(filename,netCo)
    
    genes <- union(netiRef$node,netCo$node)
    ig <- intersect(netiRef$node,netCo$node)
    tmp <- netiRef$matrix[ig,ig] + netCo$matrix[ig,ig]
    tmp[tmp>0] <- 1
    uniRef <- setdiff(netiRef$node,ig)
    uniCo  <- setdiff(netCo$node,ig)
    
    net <- list()
    net$size <- length(genes)
    net$nodes <- genes
    net$matrix <- matrix(0,net$size,net$size,dimnames=list(genes,genes))
    net$matrix[uniRef,uniRef] <- netiRef$matrix[uniRef,uniRef]
    net$matrix[uniCo,uniCo] <- netCo$matrix[uniCo,uniCo]
    net$matrix[ig,ig] <- tmp
    
    net$we <- matrix(0,net$size,net$size,dimnames=list(genes,genes))
    net$we[uniRef,uniRef] <- netiRef$we[uniRef,uniRef]
    net$we[uniCo,uniCo] <- netCo$we[uniCo,uniCo]
    net$we[ig,ig] <- pmax(netiRef$we[ig,ig],netiCo$we[ig,ig])## this the method to combine network betweenness
    ## maybe recalculte the network betweenness for combined networks
    
    #net
    
    ## write the net file
    edges <- which(net$matrix>0,arr.ind=TRUE)
    net.text <- cbind(genes[edges[,1]],genes[edges[,2]],1)
    qwt(net.text,file="data/network_inference/ComNet.txt")
    
    ## write the net betweenness
    net.we <- cbind(genes[edges[,1]],"--",genes[edges[,2]],net$we[edges])
    qwt(net.we,file="data/Network_betweenness/Betweenness_edge_ComNet.txt")
    
}

net_matrix <- function(filename){
    net.text <- as.matrix(read.table(filename,sep="\t",header=FALSE))
    net.text <- rbind(net.text,net.text[,c(2,1,3)])
    net.text[,3] <- as.numeric(net.text[,3])/max(as.numeric(net.text[,3]))
    net<- read_net(net.text)
    net
}

net_between <- function(filename,net){
    edgebe <- read.table(filename)
    pedge <- ecdf(edgebe[,4])
    
    edgebe <- edgebe[edgebe[,1] %in% net$node & edgebe[,3] %in% net$node,]
    we <- matrix(0,net$size,net$size, dimnames= list(net$node,net$node))
    subs1 <- match(edgebe[,1],net$node)
    subs2 <- match(edgebe[,3],net$node)
    we[cbind(subs1,subs2)] <- pedge(edgebe[,4])
    we[cbind(subs2,subs1)] <- pedge(edgebe[,4])
    net$we <- we
    
    net
}

