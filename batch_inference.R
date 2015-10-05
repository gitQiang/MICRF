Network_inference <- function(kk){

filenames <- c("data/network_inference/GTExbatch1.txt","data/network_inference/GTExbatch2.txt","data/network_inference/PCGCbatch1.txt","data/network_inference/PCGCbatch2.txt","data/network_inference/PCGCnewborn.txt","data/network_inference/PCGCbatch1_truned.txt","data/network_inference/PCGCbatch2_truned.txt")
strname <- c("data/network_inference/GTExbatch1module","data/network_inference/GTExbatch2module","data/network_inference/PCGCbatch1module","data/network_inference/PCGCbatch2module","data/network_inference/PCGCnewbornmodule","data/network_inference/PCGCbatch1module_t","data/network_inference/PCGCbatch2module_t")
 infernets <-  paste("data/network_inference/",c("GTExbatch1infer.txt","GTExbatch2infer.txt","PCGCbatch1infer.txt","PCGCbatch2infer.txt","PCGCnewborninfer.txt","PCGCbatch1infer_t.txt","PCGCbatch2infer_t.txt") ,sep="")
 
    ## References: a multi-method approach for proteomic network inference in 11 human cancers.
    ## Top 6 methods: RIDGENET, ARACNE-M, ARACNE-A, LASSONET, CLR AND SEPEARMANCOR
    library(stats) #cor(x, y = NULL, method = "spearman")
    library(parcor) #ridge.net, adalasso.net
    library(parmigene) #aracne.a, aracne.m,clr
    options(stringsAsFactors=FALSE)
    
    load(strname[kk])
    dexp <- read.delim(filenames[kk],check.names=F,row.names=1)

    ### infinite and missing values, and all same values
    if(kk==6){
    	dexp <- dexp[-11593,]	
    }
    dexp[is.na(dexp)] <- 0
    dexp[dexp==Inf] <- max(dexp[dexp<Inf]) + 1
    tmpsd <- apply(dexp,1,sd)
    subs <- tmpsd > 0
    dexp <- dexp[subs,]
    
    module <- module[subs,]	
    
    
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
        tmps1 <- rowSums(modexp) >0 
        tmps2 <- colSums(modexp) >0
        modexp1 <- modexp[tmps1,tmps2]
        tmp <- abs(ridge.net_hq(t(modexp1))$pcor)
        tmpM2 <- matrix(0,dim(modexp)[1],dim(modexp)[1])
        tmpM2[tmps1,tmps1] <- tmp
         M[[2]] <- tmpM2
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
    write.table(allnet,file=infernets[kk],quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
    
}

ridge.net_hq <- function (X, lambda = NULL, plot.it = FALSE, scale = TRUE, k = 10,  verbose = FALSE){
    if (is.null(lambda) == TRUE) {
        ss <- seq(-10, -1, length = 1000)
        ss <- 10^ss
        n <- nrow(X)
        nn <- n - floor(n/k)
        lambda <- ss * nn * ncol(X)
    }
    n <- nrow(X)
    p <- ncol(X)
    X <- scale(X, scale = scale)
    X[is.na(X)] <- 0  ## add by Qiang
    B <- matrix(0, nrow = p, ncol = p)
    lambda.opt <- rep(0, p)
    cat(paste("Performing local ridge regressions\n"))
    cat(paste("Vertex no "))
    for (i in 1:p) {
        if ((i/10) == floor(i/10)) {
            cat(paste(i, "..."))
        }
        noti <- (1:p)[-i]
        yi <- X[, i]
        Xi <- X[, noti]
        r.cv = ridge.cv_hq(Xi, yi, lambda = lambda, scale = scale, 
                        plot.it = plot.it, k = k)
        B[i, -i] = r.cv$coefficients
        lambda.opt[i] = r.cv$lambda.opt
    }
    pcor <- Beta2parcor(B, verbose = verbose)
    return(list(pcor = pcor))
}

ridge.cv_hq <- function (X, y, lambda = NULL, scale = TRUE, k = 10, plot.it = FALSE) 
{
    if (is.vector(X) == TRUE) {
        X <- matrix(X, ncol = 1)
    }
    if (is.null(lambda) == TRUE) {
        ss <- seq(-10, -1, length = 1000)
        ss <- 10^ss
        n <- nrow(X)
        nn <- n - floor(n/k)
        lambda <- ss * nn * ncol(X)
    }
    cv <- rep(0, length(lambda))
    n <- nrow(X)
    all.folds <- split(sample(1:n), rep(1:k, length = n))
    for (i in seq(k)) {
        omit <- all.folds[[i]]
        Xtrain = X[-omit, , drop = FALSE]
        ytrain = y[-omit]
        Xtest = X[omit, , drop = FALSE]
        ytest = y[omit]
        if ((is.vector(X) == TRUE) | (ncol(X) == 1)) {
            xtrain <- as.vector(Xtrain)
            coef.ll <- lm.ridge.univariate(xtrain, ytrain, lambda = lambda, 
                                           scale = scale)
        }
        else {
            ll <- lm.ridge_hq(ytrain ~ Xtrain, scale = scale, lambda = lambda)
            coef.ll <- coef(ll)
        }
        res <- matrix(, length(ytest), length(lambda))
        pred <- t(matrix(coef.ll[, 1], nrow = length(lambda), 
                         ncol = length(ytest))) + Xtest %*% t(coef.ll[, -1])
        res <- pred - matrix(ytest, nrow = length(ytest), ncol = length(lambda))
        cv <- cv + apply(res^2, 2, sum)
    }
    cv <- cv/n
    lambda.opt <- lambda[which.min(cv)]
    if (plot.it == TRUE) {
        plot(lambda, cv, type = "l")
    }
    if ((is.vector(X) == TRUE) | (ncol(X) == 1)) {
        x <- as.vector(X)
        coefficients <- as.vector(lm.ridge.univariate(x, y, scale = scale, 
                                                      lambda = lambda.opt))
    }
    else {
        rr <- lm.ridge_hq(y ~ X, scale = scale, lambda = lambda.opt)
        coefficients <- coef(rr)
    }
    intercept <- coefficients[1]
    coefficients <- coefficients[-1]
    return(list(intercept = intercept, coefficients = coefficients, 
                lambda.opt = lambda.opt))
}

lm.ridge_hq <- function (formula, data, subset, na.action, lambda = 0, model = FALSE, 
    x = FALSE, y = FALSE, contrasts = NULL, ...) 
{
    m <- match.call(expand.dots = FALSE)
    m$model <- m$x <- m$y <- m$contrasts <- m$... <- m$lambda <- NULL
    m[[1L]] <- quote(stats::model.frame)
    m <- eval.parent(m)
    Terms <- attr(m, "terms")
    Y <- model.response(m)
    X <- model.matrix(Terms, m, contrasts)
    n <- nrow(X)
    p <- ncol(X)
    offset <- model.offset(m)
    if (!is.null(offset)) 
        Y <- Y - offset
    if (Inter <- attr(Terms, "intercept")) {
        Xm <- colMeans(X[, -Inter])
        Ym <- mean(Y)
        p <- p - 1
        X <- X[, -Inter] - rep(Xm, rep(n, p))
        Y <- Y - Ym
    }
    else Ym <- Xm <- NA
    Xscale <- drop(rep(1/n, n) %*% X^2)^0.5
    X <- X/rep(Xscale, rep(n, p))
    X[is.na(X)] <- 0
    Xs <- svd(X)
    rhs <- t(Xs$u) %*% Y
    d <- Xs$d
    lscoef <- Xs$v %*% (rhs/d)
    lsfit <- X %*% lscoef
    resid <- Y - lsfit
    s2 <- sum(resid^2)/(n - p - Inter)
    HKB <- (p - 2) * s2/sum(lscoef^2)
    LW <- (p - 2) * s2 * n/sum(lsfit^2)
    k <- length(lambda)
    dx <- length(d)
    div <- d^2 + rep(lambda, rep(dx, k))
    a <- drop(d * rhs)/div
    dim(a) <- c(dx, k)
    coef <- Xs$v %*% a
    dimnames(coef) <- list(names(Xscale), format(lambda))
    GCV <- colSums((Y - X %*% coef)^2)/(n - colSums(matrix(d^2/div, 
        dx)))^2
    res <- list(coef = drop(coef), scales = Xscale, Inter = Inter, 
        lambda = lambda, ym = Ym, xm = Xm, GCV = GCV, kHKB = HKB, 
        kLW = LW)
    class(res) <- "ridgelm"
    res
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


GMRF_inference <- function(kk){
	filenames <- c("data/network_inference/GTExbatch1.txt","data/network_inference/GTExbatch2.txt","data/network_inference/PCGCbatch1.txt","data/network_inference/PCGCbatch2.txt","data/network_inference/PCGCnewborn.txt","data/network_inference/PCGCbatch1_truned.txt","data/network_inference/PCGCbatch2_truned.txt")
strname <- c("data/network_inference/GTExbatch1module","data/network_inference/GTExbatch2module","data/network_inference/PCGCbatch1module","data/network_inference/PCGCbatch2module","data/network_inference/PCGCnewbornmodule","data/network_inference/PCGCbatch1module_t","data/network_inference/PCGCbatch2module_t")
 infernets <-  paste("data/network_inference/",c("GTExbatch1infer.txt","GTExbatch2infer.txt","PCGCbatch1infer.txt","PCGCbatch2infer.txt","PCGCnewborninfer.txt","PCGCbatch1infer_t.txt","PCGCbatch2infer_t.txt") ,sep="")
 
    options(stringsAsFactors=FALSE)
    
    ##load(strname[kk])
    dexp <- read.delim(filenames[kk],check.names=F,row.names=1)

    ### infinite and missing values, and all same values
    if(kk==6){
    	dexp <- dexp[-11593,]	
    }
    dexp[is.na(dexp)] <- 0
    dexp[dexp==Inf] <- max(dexp[dexp<Inf]) + 1
    tmpsd <- apply(dexp,1,sd)
    subs <- tmpsd > 0
    dexp <- dexp[subs,]
    
	rho <- 0.2 ##!!!!

	library(glasso)
	s <- cov(t(dexp))
	r1 <- glasso(s,rho=rho)$wi
	save(r1,file=paste("Heart_infer_",kk,sep=""))
	
}
