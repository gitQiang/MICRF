batch_CRF_train <- function(){
    source("CRF_build.R")
    source("MixNet.R")
    source("Network_analysis.R")
    source("Multi_net.R")
    library(CRF)
    library(Corbi)
    source("CRF_train.R")
    options(stringsAsFactors=FALSE)
    fileexp <- "PCGCall.txt" 
    
    netstr <- c("STRING/","iRef/","Infer/","coexp/")
    k = 4 
    modulefiles <- c("module3","iRefmodule","DAWNmodule","DAWNmodule")
    
    netflag=7
    if(netflag==3){ instr <- "CRF_input";idm=TRUE; }else{ instr <- "hotnet_input"; idm=FALSE;}
    modulefile <- modulefiles[k]
       
    allnet <- build_net(netflag,fileexp)
    load(modulefile)
    module <- module[module[,1] %in% allnet$node,]
    
    modLab <- setdiff(unique(module[,2]),0)
    n.module <- length(modLab)
    n.genes <- length(module[module[,2]!=0,1])
    cutn=20000
    dirstr <- "result/randset3/"

    Parcrf <- matrix(,n.module,4)
    # random set 1: different sample size and exclude sample sets
    for (kk in 1:n.module){
        modgenes <- as.vector(module[module[,2]==modLab[kk],1])
        modnodesim <- c()
        Y <- c()
        for(j in 2){
            for(i in 1:10){
                strname <- paste("part",j,"_",i,sep="")
                filename <- paste(dirstr,instr,strname,".txt",sep="")
                nodeinfo <- build_node(module,1,filename)
                net <- discrete_net(allnet,modgenes,netflag,cutn,nodeinfo)
                modnodeinfo <- nodeinfo[match(modgenes,rownames(nodeinfo)),]
                modnodesim <- cbind(modnodesim,modnodeinfo)
                TADAFile <- paste("ASD/TADAresult/randset3/TADAdenovo_part",j,"_",i,".csv",sep="")
                TADAall <- read.csv(TADAFile)
                Tset <- TADAall[TADAall[,"qvalue.dn"]<= 0.2,1]
                y <- rep(2,length(modnodeinfo)) 
                y[match(Tset,names(modnodeinfo))] <- 1
                
                if(all(y==2)){
                    y[which.max(modnodeinfo)] <- 1
                }
                
                Y <- cbind(Y,y)
            }
        }
        
        modenodesim <- as.matrix(modnodesim)
        Y <- as.matrix(Y)
        crf <- CRF_train(modnodesim,net,Y)
        Parcrf[kk,] <- crf$par
    }
    
    
}

CRF_train <- function(modnodesim,net,Y){
    
    ### Y and features
    feas <- build_fea(modnodesim,net,Y)
    FeaN <- feas$FeaN
    n.nf <- 1
    FeaE <- feas$FeaE
    n.ef <- 1
    
    K=1
    nsample <- dim(Y)[2]
    resample <- sample.cross(nsample,K)


    for(i in 1:K){
        node.fea <- list()
        edge.fea <- list()
        crf0 <- build.crf(net$matrix,net$size)
        subtrain <- resample$train[[i]]
        ntrain <- length(subtrain)
        Strain <- t(Y[,subtrain,drop=FALSE])
        for(j in 1:ntrain){
            node.fea[[j]] <- FeaN[[subtrain[j]]]
            edge.fea[[j]] <- FeaE[[subtrain[j]]]
        }
        crf <- train.crf(crf0, Strain, node.fea, edge.fea)   
    }
    
    crf
}

build_fea <- function(modnodesim,net,Y){
    
    n.nf <- 1
    FeaN <- list()
    FeaE <- list()
    nsample <- dim(Y)[2]
    
    n.states <- 2
    crf <- make.crf(net$matrix, rep(n.states,net$size))
    deg <- colSums(net$matrix)
        
    for(i in 1:nsample){
        S <- modnodesim[,i]
        S <- cbind(S,1-S)
        colnames(S) <- 1:2
        FeaN[[i]] <- S[1:dim(S)[1],Y[,i]]

        tmp <- 1:crf$n.edges
        W <- matrix(1,n.states,n.states)
        for (e in 1:crf$n.edges){
            n1 <- crf$edges[e, 1]
            n2 <- crf$edges[e, 2]
            W[1,1] <- max(S[n1,2],S[n2,2])/sqrt(deg[n1]*deg[n2])
            W[2,2] <- 0
            W[1,2] <- 0
            W[2,1] <- 0
            tmp[e] <- W[Y[n1,i],Y[n2,i]]
        }
        FeaE[[i]] <- tmp   
    }
    list(FeaN=FeaN,FeaE=FeaE)
}

build_model <- function(modnodesim,net,bflag=1){   
    
    n.nf <- dim(as.matrix(modnodesim))[2]
    #S <- modnodesim
    S <- node_information(modnodesim,n.nf,bflag)
    # bflag==1 # different node features as different features
    # bflag==2 # integrate different node features first
    S <- cbind(S,1-S)
    colnames(S) <- 1:2
    
    n.states <- 2
    n.node <- net$size
    query.net <- net$matrix
    crf <- make.crf(net$matrix, rep(n.states,net$size))
    
    crf$state.map <- matrix(n.node, nrow=crf$n.nodes, ncol=crf$max.state)
    for (i in 1:crf$n.nodes)
    {
        crf$state.map[i, 1:crf$n.states[i]] <- 1:n.states
        crf$node.pot[i,] <- exp(S[i, crf$state.map[i,]])
    }   
    
    deg <- colSums(net$matrix)
    W <- matrix(1,n.states,n.states)
    W1 <- W
    W2 <- W
    for (e in 1:crf$n.edges){
        n1 <- crf$edges[e, 1]
        n2 <- crf$edges[e, 2]
        m1 <- 1:crf$n.states[n1]
        m2 <- 1:crf$n.states[n2]
        S1 <- matrix(log(crf$node.pot[n1, m1]), crf$n.states[n1], crf$n.states[n2])
        S2 <- matrix(log(crf$node.pot[n2, m2]), crf$n.states[n1], crf$n.states[n2], byrow=T)
        
        W[1,1] <- exp(-abs(S1[1,1]-S2[1,1])^2)
        W[2,2] <- exp(-abs(S1[2,2]-S2[2,2])^2)
        W[1,2] <- min(exp(-abs(S1[1,1]-S2[1,2])^2),W[1,1])
        W[2,1] <- min(exp(-abs(S1[2,1]-S2[1,1])^2),W[2,2])
        #         crf$edge.pot[[e]] <- W*exp(S1+S2)/2
        
        W1[1,1] <- (max(S1[1,1],S2[1,1])+0.3)/sqrt(deg[n1]*deg[n2])
        W1[2,2] <- 0
        W1[1,2] <- 0
        W1[2,1] <- 0
        crf$edge.pot[[e]] <- exp(W1)
    }
    crf
}

allFea <- function(chr,C,nCG){
    
    
    
    list(X=X,FeaN=FeaN,FeaE=FeaE)
}

sample.cross <- function(nsample,K){
    
    train_sample <- list()
    pred_sample <- list()
    
    if(nsample==1 | K==1){
        train_sample[[1]] <- 1:nsample
    }else{
    nk <- floor(nsample/K)
    sam <- sample(nsample)
    #subsK <- matrix(K-1,nk)
    #subK <- sam[((K-1)*nk+1):nsample]
    
    for(i in 1:(K-1)){
        pred_sample[[i]] <- sam[((i-1)*nk+1):(i*nk)]
        train_sample[[i]] <- setdiff(sam,pred_sample[[i]])
    }
    pred_sample[[K]] <- sam[((K-1)*nk+1):nsample]
    train_sample[[K]] <- setdiff(sam,pred_sample[[K]])
    }
    
    list(train=train_sample,pred=pred_sample)
}

build.crf <- function(adj,n.nodes,n.states=2,n.nf=1,n.ef=1){
    
    library(CRF)
    crf0 <- make.crf(adj, n.states)
    crf0 <- make.features(crf0,n.nf,n.ef)
    
    n.par <- 0
    for (i in 1:n.nf) {
        for (j in 1:(n.states-1)) {
            n.par <- n.par + 1
            crf0$node.par[,j,i] <- n.par
        }
    }
    
    for (i in 1:n.states) {
        for (j in 1:n.states) {
            n.par <- n.par + 1
            crf0$edge.par[[1]][i,j,1] <- n.par
        }
    }
    n.par <- n.par - 1
    crf0$edge.par[[1]][n.states,n.states,1] <- 0
    
    for (i in 2:crf0$n.edges) {
        crf0$edge.par[[i]] <- crf0$edge.par[[1]]
    }
    
    crf0 <- make.par(crf0, n.par)
    
    crf0
    
}

predict.crf <- function(crf,FeaN,FeaE,subpred){
    
    npred <- length(subpred)
    result <- list()
    for(i in 1:npred){
        crf <- crf.update(crf,FeaN[[subpred[i]]],FeaE[[subpred[i]]])
        result[[i]] <- decode.chain(crf)
    }
    
    result
}

edge.Feas <- function(chr,nsample,TCG,nCG){
    
    FeaE <- list()
    subs <-  seq(1,TCG,nCG);
    # for(i in 1:nsample){
    # FeaE[[i]] <- chr[(subs[i]+1):(subs[i]+nCG-1),2] - chr[subs[i]:(subs[i]+nCG-2),2]
    # FeaE[[i]] <- as.matrix(FeaE[[i]]/max(FeaE[[i]]))
    # }
    betap <- mean(chr[2:TCG,2] - chr[1:(TCG-1),2])
    for(i in 1:nsample){
        FeaE[[i]] <- as.matrix(1-exp((-2*(chr[(subs[i]+1):(subs[i]+nCG-1),2] - chr[subs[i]:(subs[i]+nCG-2),2]))/betap))
    }
    
    FeaE
    
}

node.Feas <- function(chr){
    
    N <- dim(chr)[2]
    # depart <- matrix(0,5,2)
    # depart[1,] <- c(5,23)
    # depart[1,] <- c(13,23)
    # depart[2,] <- 12+c(12,16)
    # depart[3,] <- 12+c(17,21)
    # depart[4,] <- 12+c(22,26)
    # depart[5,] <- 12+c(27,N-12)
    # F <- chr[,depart[1,1]:depart[1,2]]
    # for(i in 1:4){
    # F <- cbind(F,rowSums(chr[,depart[i+1,1]:depart[i+1,2]]))
    # }
    F <- chr[,5:13]
    #F <- chr[,c(5,6,7,13)]
    
    F
    
}

