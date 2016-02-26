CRF_edgeF <- function(){
X <- seq(0,1,0.001)
Y <- sapply(X,function(x){
    y=0.001
a <- matrix(c(x,1-x),1,2)
b <- matrix(c(y,1-y),2,1)
T0 <- matrix(1,2,2)
T1 <- matrix(1,2,2)

T0[1,1] <- exp(-abs(a[1]-b[1])^2)
T0[2,2] <- exp(-abs(a[2]-b[2])^2)
T0[1,2] <- min(exp(-abs(a[1]-b[2])^2),T0[1,1])
T0[2,1] <- min(exp(-abs(a[2]-b[1])^2),T0[2,2])

T1[1,1] <- 1/4
T1[2,2] <- 0
T1[1,2] <- 0
T1[2,1] <- 0
T1 <- exp(T1)

T0 <- T0*T1

y0 <- t(a%*%T0) * b
y1 <- t(a%*%T1) * b

t0 <- y0[1]/sum(y0)
t1 <- y1[1]/sum(y1)

c(t0,t1)
})
plot(X,Y[1,])
plot(X,Y[2,])
lines(x=c(0,1),y=c(0.9,0.95),type="l")
abline(h=0.4)
abline(v=0.5)
}

CRF_edgeF1 <- function(){
x <- seq(0,1,0.01)
y <- sapply(x, function(k){
    S1<- matrix(c(0.8,0.2),2,2)
    S2 <- matrix(c(k,1-k),2,2,byrow=TRUE)
    deg <- c(30,50)
    n1=1
    n2=2
    S <- matrix(0,2,2)
    #c <- min(1/(1- sqrt(S1[1,1]*S2[1,1]))*plogis((sqrt(S1[1,1]*S2[1,1])-0.5) ,location=0, scale=ifelse((sqrt(S1[1,1]*S2[1,1])-0.5) >= 0, (1- sqrt(S1[1,1]*S2[1,1])), sqrt(S1[1,1]*S2[1,1]) ))/(sqrt(deg[n1]*deg[n2])), 0.5)
    #c <- min((1-sqrt(S1[1,1]*S2[1,1]))*1/(sqrt(deg[n1]*deg[n2]))*exp(S1[1,1]*S2[1,1]), 0.5)
    #c <- min(plogis(k*0.7,location=k*0.7-0.1,scale=1/(1-k*0.7)),0.5)
    #c <- min(exp(((S1[1,1]+S2[1,1]))/abs(S1[1,1]*S2[1,1]))/(sqrt(deg[n1]*deg[n2])*exp(S1[1,1]*S2[1,1]-1.1)),0.5)
    c <- min(exp(abs(min(S1[1,1],S2[1,1])-1)/abs(S1[1,1]*S2[1,1]-sqrt(S1[1,1]*S2[1,1])))/(sqrt(deg[n1]*deg[n2])*exp(S1[1,1]*S2[1,1]-1)),0.5)
    print(c)
    
    S[1,1] <- 0.5+c; S[2,1] <- 1-S[1,1];
    S[2,2] <- 0.5+c; S[1,2] <- 1-S[2,2];
    
    S= 1/(1-log(1-abs(S1-S2)))
    if(S[1,1] > 0.5){
        S[2,1] <- S[1,1]/10;
        S[1,2] <- S[2,2]/10;
    }else{
        S[2,1] <- 0;
        S[1,2] <- 0;
    }
    pos2 <- (S*(S1+S2)/2) %*% (t(S2)[,1])
    pos2[1]/sum(pos2)
    
}
)
plot(x,y)
lines(x,x,type="l")
abline(v=0.5)
abline(h=0.7)
abline(v=0.7)

}
