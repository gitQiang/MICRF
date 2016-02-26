betweenness_hq <- function(){
    
    library(igraph)
  
    if(netflag==3){ filename <- "data/STRINGnetmap.txt";}
    if(netflag==7){ filename <- "data/network_inference/ASDexp_net_top5.txt";} ## ASD co-expression network
    if(netflag==20){ filename <- "data/StringNew_HPRDnet.txt";} ## MAGI
    
    net.text <- as.matrix(read.table(filename,sep="\t",header=FALSE))
    net.text <- rbind(net.text,net.text[,c(2,1,3)])
    
   
    allnet <- build_net(netflag,"")
    allnet$matrix[ allnet$matrix>0 ] <- 1
    
    g1 <- graph.adjacency(allnet$matrix, mode="undirected", weighted=NULL, diag=FALSE) #, add.colnames=as.vector(allnet$node), add.rownames=as.vector(allnet$node))
    
    bnode <- betweenness(g1, v=V(g1), directed = FALSE, weights = NULL)
    bedge <- edge.betweenness(g1, e=E(g1), directed = FALSE, weights = NULL)
    
    save(bnode,file="betweenness_node_3")
    save(bedge,file="betweenness_edge_3")
    
    
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