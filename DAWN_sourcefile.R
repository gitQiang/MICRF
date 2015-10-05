#####################################################################################################
##DAWN package source code, most recent revision: Nov 16, 2014
#####################################################################################################


##load useful library
library("WGCNA")

###function: clean up expression data:
cleanexp=function(filename,outputname){
  filel=length(filename)
  brain.genes=read.table("Brain_expression.txt",header=T)
  table(brain.genes$Brain_expressed)
  brain.genes=brain.genes[brain.genes$Brain_expressed != 'No',]
  for (i in 1:filel){
    ### read  the original data
    data.in=read.table(filename[i], header=TRUE, sep="\t")
    gene.names=data.in$Gene.Symbol
    
    #retain column 2 and up
    data=as.matrix(data.in[,2:ncol(data.in)])
    row.names(data)=gene.names
    
    # delete genes not in the brain.genes
    i.gene=which(rownames(data) %in% brain.genes$Gene)
    data=data[i.gene,]
    
    # now select observations that are in the frontal cortex region
    # not needed when using all the data
    FCplus=c("MFC","DFC","OFC","MSC","VFC","M1C","S1C")
    region=unlist(lapply(colnames(data),sampleNames,first=1,last=2))
    i.region=which(region %in% FCplus)
    table(region[i.region])
    data=data[,i.region]
    
    #### take averages for duplicated genes
    gene.names=as.vector(rownames(data))
    sample.names=as.vector(colnames(data))
    data.upd=sapply(unique(gene.names), RowColMeans, row.col='row', data=data, all.names=gene.names)
    data.upd=t(matrix(as.numeric(unlist(data.upd)),ncol(data),length(unique(gene.names))))
    rownames(data.upd)=unique(gene.names)
    colnames(data.upd)=sample.names
    dim(data.upd)
    
    
    data=data.upd
    
    #### save the data set to use for the analysis
    save(data,file=outputname[i])
  } 
}

cleanexp_v2=function(filename,outputname){
  filel=length(filename)
  brain.genes=read.table("/Users/liliu/Documents/jing/expression_data/Brain_expression.txt",header=T)
  table(brain.genes$Brain_expressed)
  brain.genes=brain.genes[brain.genes$Brain_expressed != 'No',]
  for (i in 1:filel){
    ### read  the original data
    data.in=read.table(filename[i], header=TRUE, sep="\t")
    gene.names=data.in$Gene.Symbol
    
    #retain column 2 and up
    data=as.matrix(data.in[,2:ncol(data.in)])
    row.names(data)=gene.names
    
    # delete genes not in the brain.genes
    i.gene=which(rownames(data) %in% brain.genes$Gene)
    data=data[i.gene,]
    
    # now select observations that are in the frontal cortex region
    # not needed when using all the data
    FCplus=c("MD","CBC")
    region=unlist(lapply(colnames(data),sampleNames,first=1,last=2))
    i.region=which(region %in% FCplus)
    table(region[i.region])
    data=data[,i.region]
    
    #### take averages for duplicated genes
    gene.names=as.vector(rownames(data))
    sample.names=as.vector(colnames(data))
    data.upd=sapply(unique(gene.names), RowColMeans, row.col='row', data=data, all.names=gene.names)
    data.upd=t(matrix(as.numeric(unlist(data.upd)),ncol(data),length(unique(gene.names))))
    rownames(data.upd)=unique(gene.names)
    colnames(data.upd)=sample.names
    dim(data.upd)
    
    
    data=data.upd
    
    #### save the data set to use for the analysis
    save(data,file=outputname[i])
  } 
}



sampleNames <- function(x,first,last){
  xx=unlist(strsplit(x,""))
  i=which(xx == '.')
  if(first == 0){
    return(substr(x,1,i[last]-1))
  }else{
    return(substr(x,i[first]+1,i[last]-1))
  }
}

RowColMeans <- function(name,row.col,data,all.names){  
  i=which(all.names == name)
  if(i[1]%%500 == 0)print(c(i[1],name));flush.console()
  if(length(i) > 1){
    if(row.col == 'row'){ 	
      means=colMeans(data[i,],na.rm=T)
    }else{
      means=rowMeans(data[,i],na.rm=T)
    }	
    return(means)
  }else{
    return(data[i[1],])
  }
}

regionMeans <- function(sample, x.sample, x.region, x.composite, x.data){
  i = which(x.sample == sample & x.region %in% x.composite) 
  if(length(i) == 0){
    return(rep(NA,nrow(x.data)))
  } else if(length(i) == 1){
    return(x.data[,i])
  } else {
    return(rowMeans(x.data[,i]))
  }
}

cor_matrix_f=function(filename,rowindex){
  con<-file(filename,"r")
  line<-readLines(con,1)
  L=1
  corm=rep(0,length(rowindex))
  for (i in 1:17253){
    line<-readLines(con,1)
    if (sum(i==rowindex)!=0){
      mvector<-strsplit(line, ",")
      vec=mvector[[1]]
      corr=as.numeric(vec[2:length(vec)])
      corm=rbind(corm,corr[rowindex])
    }
    if (max(c(i,rowindex))==i){
      close(con)
      return(corm[-1,])
    }
  }
  
}


minim=function(p,label){ ##input p-values and labels, output pvalues for each subgroup
  k=length(table(label)) ##get the number of classes
  fishp=rep(0,k)
  for (i in 1:k){
    puse=p[label==i]
    fishp[i]=min(puse)
  }
  return(fishp)
}


minim_j=function(p,label){ ##input p-values and labels, output pvalues for each subgroup
  k=length(table(label)) ##get the number of classes
  fishpj=rep(0,k)
  fishpk=rep(0,k)
  fishpl=rep(0,k)
  for (i in 1:k){
    puse=p[label==i]
    n=length(puse)
    padj = min(p.adjust(puse,"fdr"))
    pmin = min(puse)
    pt = 1-pnorm(sqrt(n)*mean(qnorm(puse,lower.tail=FALSE)))       
    pj = 2*min(padj,pt)
    pk = pnorm(.5*(qnorm(1-padj)+qnorm(1-pt)),lower.tail=FALSE)
    pl = pnorm(.5*(qnorm(1-pmin)+qnorm(1-pt)),lower.tail=FALSE)    
    fishpj[i]=pj
    fishpk[i]=pk
    fishpl[i]=pl
  }
  return(data.frame(fishpj,fishpk,fishpl))
}


combinecorr=function(cmatrix,label){
  k=length(table(label))
  corm=matrix(0,k,k)
  for (i in 1:k){
    for (j in 1:k){
      corm[i,j]=sum(cmatrix[label==i,label==j])/(sum(label==i)*sum(label==j))
    }
  } 
  return(corm)
}

plf3=function(b,c,graph,zupdate){ ##psuedo likelihoood 1
  fvalue=1
  for (i in 1:length(zupdate)){
    new1=exp(b*zupdate[i]+c*zupdate[i]*t(graph[i,])%*%zupdate) 
    new2=exp(b*(1-zupdate[i])+c*(1-zupdate[i])*t(graph[i,])%*%zupdate) 
    fvalue=fvalue+log(new1/(new1+new2))
  }  
  return(fvalue) 
}


optimbc=function(plf2,graph1,Iupdate,times){
  biter=0
  citer=0
  for (k in 1:times){
    biter1=optimize(plf2,c(-10,10),c=citer,graph=graph1,zupdate=Iupdate,maximum=T)$maximum
    citer1=optimize(plf2,c(0,10),b=biter1,graph=graph1,zupdate=Iupdate,maximum=T)$maximum
    if (abs(biter1-biter)<10^-5&abs(citer1-citer)<10^-5){
      break
    }
    biter=biter1
    citer=citer1
  }
  return(c(biter,citer))
}

emf=function(zscore){
  d1=dnorm(zscore) 
  d2=dnorm(zscore,2,1) ##initial mu is 3
  pi=mean(d2/(d1+d2)) ##initial value
  for (j in 1:1000){
    mu=sum(d2*zscore/(d1+d2))/sum(d2/(d1+d2)) ##update mean
    pi=mean(d2/(d1+d2)) ##update pi
    d1=(1-pi)*dnorm(zscore) ## calculate posterior prob
    d2=pi*dnorm(zscore,mu,1)
  }    
  I=(d2>d1)+0
  return(list(mu,pi,I))
}

emf2=function(zscore){ ##with not fixed sigma
  d1=dnorm(zscore) 
  d2=dnorm(zscore,2,1) ##initial mu is 3
  pi=mean(d2/(d1+d2)) ##initial value
  for (j in 1:1000){
    mu=sum(d2*zscore/(d1+d2))/sum(d2/(d1+d2)) ##update mean
    sigmas=sum(d2*(zscore-mu)^2/(d1+d2))/sum(d2/(d1+d2))
    pi=mean(d2/(d1+d2)) ##update pi
    d1=(1-pi)*dnorm(zscore) ## calculate posterior prob
    d2=pi*dnorm(zscore,mu,sqrt(sigmas))
  }    
  I=(d2>d1)+0
  return(list(mu,pi,I,sigmas))
}


emf3=function(zscore){ ##with both null and alternative not fixed sigma
  d1=dnorm(zscore) 
  d2=dnorm(zscore,2,1) ##initial mu is 3
  pi=mean(d2/(d1+d2)) ##initial value
  for (j in 1:1000){
    mu=sum(d2*zscore/(d1+d2))/sum(d2/(d1+d2)) ##update mean
    sigmas2=sum(d2*(zscore-mu)^2/(d1+d2))/sum(d2/(d1+d2))
    sigmas1=sum(d1*(zscore)^2/(d1+d2))/sum(d1/(d1+d2))
    pi=mean(d2/(d1+d2)) ##update pi
    d1=(1-pi)*dnorm(zscore,sd=sqrt(sigmas1)) ## calculate posterior prob
    d2=pi*dnorm(zscore,mu,sqrt(sigmas2))
  }    
  I=(d2>d1)+0
  return(list(mu,pi,I,sigmas1,sigmas2))
}

emf4=function(zscore){ ##with both null and alternative not fixed sigma, but they are the same
  d1=dnorm(zscore) 
  d2=dnorm(zscore,2,1) ##initial mu is 2
  pi=mean(d2/(d1+d2)) ##initial value
  for (j in 1:1000){
    mu=sum(d2*zscore/(d1+d2))/sum(d2/(d1+d2)) ##update mean
    sigmas2=sum(d2*(zscore-mu)^2/(d1+d2))/sum(d2/(d1+d2))
    sigmas1=sum(d1*(zscore)^2/(d1+d2))/sum(d1/(d1+d2))
    sigmas=(sum(d2)*sigmas2+sum(d1)*sigmas1)/(sum(d1+d2))
    sigmas2=sigmas
    sigmas1=sigmas
    pi=mean(d2/(d1+d2)) ##update pi
    d1=(1-pi)*dnorm(zscore,sd=sqrt(sigmas1)) ## calculate posterior prob
    d2=pi*dnorm(zscore,mu,sqrt(sigmas2))
  }    
  I=(d2>d1)+0
  return(list(mu,pi,I,sigmas1,sigmas2))
}

emf5=function(zscore){ ##with both null and alternative not fixed sigma, but they are the same for 1000 times, then 3 times fix sigma_1=1
  d1=dnorm(zscore) 
  d2=dnorm(zscore,2,1) ##initial mu is 2
  pi=mean(d2/(d1+d2)) ##initial value
  for (j in 1:1000){
    mu=sum(d2*zscore/(d1+d2))/sum(d2/(d1+d2)) ##update mean
    sigmas2=sum(d2*(zscore-mu)^2/(d1+d2))/sum(d2/(d1+d2))
    sigmas1=sum(d1*(zscore)^2/(d1+d2))/sum(d1/(d1+d2))
    sigmas=(sum(d2)*sigmas2+sum(d1)*sigmas1)/(sum(d1+d2))
    sigmas2=sigmas
    sigmas1=sigmas
    pi=mean(d2/(d1+d2)) ##update pi
    d1=(1-pi)*dnorm(zscore,sd=sqrt(sigmas1)) ## calculate posterior prob
    d2=pi*dnorm(zscore,mu,sqrt(sigmas2))
  }  
  for (j in 1:3){
    mu=sum(d2*zscore/(d1+d2))/sum(d2/(d1+d2)) ##update mean
    sigmas2=sum(d2*(zscore-mu)^2/(d1+d2))/sum(d2/(d1+d2))
    sigmas1=1
    pi=mean(d2/(d1+d2)) ##update pi
    d1=(1-pi)*dnorm(zscore,sd=sqrt(sigmas1)) ## calculate posterior prob
    d2=pi*dnorm(zscore,mu,sqrt(sigmas2))
  }
  I=(d2>d1)+0
  return(list(mu,pi,I,sigmas1,sigmas2))
}

modulef8=function(mgindex,genename,genename2,pv,corm,thres1,thres2){ ##use simple initial cut
  genelist=genename[mgindex]
  
  glindex=rep(0,length(genelist))
  for (gi in 1:length(genelist)){
    if (sum(genelist[gi]==genename2)!=0){
      glindex[gi]=c(1:length(genename2))[genename2==genelist[gi]]
    }  
  }
  genelist=genelist[glindex!=0]  
  pv=pv[glindex[glindex!=0]] 
  
  corm=corm[glindex!=0,glindex!=0]
  dissim=1-abs(corm)
  distance=as.dist(dissim)
  result1=hclust(distance,method="average")
  
  #  sublabel=cutree(result1, h=0.25) 
  sublabel=cutree(result1, h=thres1) 
  
  corm2=abs(corm) 
  supercorm=combinecorr(corm2,sublabel)  
  thres=thres2  
  usecorrm=supercorm
  graph1=(usecorrm>thres+0)-diag(1,dim(usecorrm)[1])   ##build the graph
  
  pmeta=minim(pv,sublabel) ##obtain the pvalue of super nodes
  #  pmetas=minim_j(pv,sublabel)
  #  pmetaj=pmetas[,1]
  #  pmetak=pmetas[,2]
  #  pmetal=pmetas[,3]
  
  #  pmeta=pmetak
  # Istart=(pmeta<0.1)+0 #### Qiang Huang changed
  if(min(pmeta) < 0.1){
    Istart=(pmeta<0.1)+0
  }else{
    Istart=(pmeta<=min(pmeta))+0
  }
  
  zscore=qnorm(1-pmeta)
  
  seedindex=(pmeta<=10^{-8})+0
  mustart=mean(zscore[Istart==1&seedindex==0])

  sigma2=sd(zscore[Istart==1&seedindex==0])
  sigma1=sd(zscore[Istart==0])

  # Qiang Huang changed
  if(is.na(sigma2)) sigma2 <- 1
  if(is.na(sigma1)) sigma1 <- 1
  
#  print(mustart)
#  print(table(Istart))
  
  
  mainresult=main_f9(graph1,zscore,Istart,100,mustart,seedindex,sigma1^2,sigma2^2)
   
  output=list()
  output$glindex=glindex
  output$graph=graph1
  output$sublabel=sublabel
  output$Istart=Istart
  output$Ifinal=mainresult$Ifinal
  output$posterior=mainresult$pfinal
  output$b=mainresult$b
  output$c=mainresult$c
  output$mu=mainresult$mu1
  output$sigmas1=mainresult$sigmas1
  output$sigmas2=mainresult$sigmas2
  output$pmeta=pmeta
  output$genelist=genelist
  output$seedindex=seedindex
  return(output)
}



 
main_f9=function(graph,zscore,Istart,iter,mui,seedindex,sigmas1,sigmas2){   
  d=length(Istart)
  posterior=rep(0,d)
  Iupdate=Istart
  c1=0
  b1=0
  mu1=mui
  for (ite in 1:iter){
    # Qiang Huang changed
    # print(ite)
    
    res=optimbc(plf3,graph,Iupdate,20)
    b=res[1]
    c=res[2]
    if (abs(c1-c)<0.001&abs(b1-b)<0.001) { 
      break
    }
    c1=c
    b1=b  
    for (i in 1:d){
      new1=exp(b*Iupdate[i]+c*Iupdate[i]*t(graph[i,])%*%Iupdate)
      new2=exp(b*(1-Iupdate[i])+c*(1-Iupdate[i])*t(graph[i,])%*%Iupdate)
      # Qiang Huang changed
      
      p1=dnorm(zscore[i],mu1*Iupdate[i], ifelse(sqrt(sigmas2*Iupdate[i]+sigmas1*(1-Iupdate[i]))==0,1,sqrt(sigmas2*Iupdate[i]+sigmas1*(1-Iupdate[i]))) )*new1/(new1+new2) ##sigmas2 is the alternative, sigmas1 is the null var
      p2=dnorm(zscore[i],mu1*(1-Iupdate[i]), ifelse(sqrt(sigmas2*(1-Iupdate[i])+sigmas1*Iupdate[i])==0,1,sqrt(sigmas2*(1-Iupdate[i])+sigmas1*Iupdate[i]))  )*new2/(new1+new2)
      
      #p1=dnorm(zscore[i],mu1*Iupdate[i],sqrt(sigmas2*Iupdate[i]+sigmas1*(1-Iupdate[i])))*new1/(new1+new2) ##sigmas2 is the alternative, sigmas1 is the null var
      #p2=dnorm(zscore[i],mu1*(1-Iupdate[i]),sqrt(sigmas2*(1-Iupdate[i])+sigmas1*Iupdate[i]))*new2/(new1+new2)
      if (Iupdate[i]==1)
        posterior[i]=p1/(p1+p2)
      if (Iupdate[i]==0)
        posterior[i]=p2/(p1+p2)
      if (p2>p1){
        Iupdate[i]=1-Iupdate[i]
      }
      
    }
    mu1=sum(posterior[seedindex==0]*zscore[seedindex==0])/sum(posterior[seedindex==0])
    sigmas2=sum(posterior[seedindex==0]*(zscore[seedindex==0]-mu1)^2)/sum(posterior[seedindex==0])
    sigmas1=sum((1-posterior[seedindex==0])*(zscore[seedindex==0])^2)/sum(1-posterior[seedindex==0])
    sigmas=(sigmas1*sum(posterior[seedindex==0])+sigmas2*sum(1-posterior[seedindex==0]))/length(posterior)
    sigmas2=sigmas
    sigmas1=sigmas
  }
  for (ite in 1:4){
    res=optimbc(plf3,graph,Iupdate,20)
    b=res[1]
    c=res[2]
    sigmas1=1
    #    sigmas2=1
    #    sigmas2=1 ##assume it's 1, or we can give initial value?
    for (i in 1:d){
      new1=exp(b*Iupdate[i]+c*Iupdate[i]*t(graph[i,])%*%Iupdate)
      new2=exp(b*(1-Iupdate[i])+c*(1-Iupdate[i])*t(graph[i,])%*%Iupdate)
      # Qiang Huang changed
      p1=dnorm(zscore[i],mu1*Iupdate[i],  ifelse(sqrt(sigmas2*Iupdate[i]+sigmas1*(1-Iupdate[i]))==0,1,sqrt(sigmas2*Iupdate[i]+sigmas1*(1-Iupdate[i])))   )*new1/(new1+new2) ##sigmas2 is the alternative, sigmas1 is the null var
      p2=dnorm(zscore[i],mu1*(1-Iupdate[i]),  ifelse(sqrt(sigmas2*(1-Iupdate[i])+sigmas1*Iupdate[i])==0,1,sqrt(sigmas2*(1-Iupdate[i])+sigmas1*Iupdate[i]))  )*new2/(new1+new2)
      
      # p1=dnorm(zscore[i],mu1*Iupdate[i],sqrt(sigmas2*Iupdate[i]+sigmas1*(1-Iupdate[i])))*new1/(new1+new2) ##sigmas2 is the alternative, sigmas1 is the null var
      # p2=dnorm(zscore[i],mu1*(1-Iupdate[i]),sqrt(sigmas2*(1-Iupdate[i])+sigmas1*Iupdate[i]))*new2/(new1+new2)
      if (Iupdate[i]==1)
        posterior[i]=p1/(p1+p2)
      if (Iupdate[i]==0)
        posterior[i]=p2/(p1+p2)
      if (p2>p1){
        Iupdate[i]=1-Iupdate[i]
      }
    }
    mu1=sum(posterior[seedindex==0]*zscore[seedindex==0])/sum(posterior[seedindex==0])
    sigmas2=sum(posterior[seedindex==0]*(zscore[seedindex==0]-mu1)^2)/sum(posterior[seedindex==0])
    #   sigmas2=1
    
  }
  out=list()
  out$Ifinal=Iupdate
  out$pfinal=posterior
  out$b=b
  out$c=c
  out$mu1=mu1
  out$sigmas1=sigmas1
  out$sigmas2=sigmas2
  return(out)
}

 
reportpf=function(graph,pfinal,newconel,pmeta,keepgene,keepgenep,Istart,Ifinal){
  pfnew=1-pfinal
  clusterr=rank(pfnew, ties.method= "first")
  
  pbe=rep(0,length(keepgene))
  po=rep(0,length(keepgene))
  paf=rep(0,length(keepgene))
  geneclass=rep(0,length(keepgene))
  Initial=rep(0,length(keepgene))
  Final=rep(0,length(keepgene))
  riskedge=rep(0,length(keepgene))
  totaledge=rep(0,length(keepgene))
  
  for (i in 1:length(pfnew)){
    geneclass[newconel==i]=clusterr[i]
  }
  
  for (i in 1:length(pfnew)){
    pbe[newconel==i]=pmeta[i]
  }
  
  for (i in 1:length(pfnew)){
    Initial[newconel==i]=Istart[i]
  }
  
  for (i in 1:length(pfnew)){
    Final[newconel==i]=Ifinal[i]
  }
  
  for (i in 1:length(pfnew)){
    paf[newconel==i]=pfnew[i]
  }
  
  for (i in 1:length(pfnew)){
    riskedge[newconel==i]=sum(Ifinal[graph[i,]==1])
  }
  
  for (i in 1:length(pfnew)){
    totaledge[newconel==i]=sum(graph[i,]==1)
  }
  
  save1=data.frame(keepgene,geneclass,paf,pbe,keepgenep,Initial,Final,riskedge,totaledge)
  return(save1)
}

 
 
stage2_mixturemodel=function(save1,thres){ 
  paf=save1$paf
  keepi=c(1:length(paf))[paf<thres]
  kgene=as.character(save1$keepgene)[keepi]
  kclass=save1$geneclass[keepi]
  kmodule=save1$modulelabel[keepi]
  kpvalue=save1$keepgenep[keepi]
  
  zscore=qnorm(1-kpvalue)
  countc=table(kclass)
  counti=rep(0,length(kclass))
  for (i in 1:length(counti)){
    counti[i]=c(1:length(countc))[names(countc)==as.character(kclass[i])]
  }
  ikeep=c(1:length(kclass))[countc[counti]>9]
  #  browser()
  
  zscore3=zscore[ikeep]
  kclassuse3=kclass[ikeep] ##only keep the genes within the big nodes
  
  
  allclass=unique(kclass) ##assign new class levels to the genes kept
  newclass=rep(0,length(zscore3))
  k=1
  for (i in 1:length(allclass)){
    if (sum(kclassuse3==allclass[i])!=0){
      genei=c(1:length(kclassuse3))[kclassuse3==allclass[i]]
      newclass[genei]=k
      k=k+1
    }
  }
  
  theta=rep(0.2,length(unique(newclass))) ##begin to run the algorithm to estimate parameters
  out4=emf7_node(zscore3,newclass,theta,0,1.5,1,1)
  
#  print(out4$sigma1)
#  print(out4$sigma2)
#  print(out4$mu2)
# Qiang Huang changed
  sigma1=ifelse(out4$sigma1==0,1,out4$sigma1)
  sigma2=ifelse(out4$sigma2==0,1,out4$sigma2)

  #sigma1=out4$sigma1
  #sigma2=out4$sigma2
  mu2=out4$mu2
  post=out4$post
  theta=out4$theta
  
  mtheta=mean(theta) ##begin to 
  d1=dnorm(zscore,0,sigma1)*(1-mtheta)
  d2=dnorm(zscore,mu2,sigma2)*mtheta
  finalposter=d2/(d1+d2) ##step 1, calculate posterior
  nodeindex=c(1:length(zscore)) 
  nodefi=nodeindex[ikeep] ##find out the index of gene in the large nodes
  
  finaltheta=rep(mtheta,length(zscore))
  finaltheta[nodefi]=theta[newclass]
  finalposter[nodefi]=post
  finalgene=kgene 
  finalop=kpvalue 
  finalkclass=kclass
  
  usepost=1-finalposter
  rankpost=sort(usepost)
  localfdr=rep(0,length(usepost))
  for (i in 1:length(localfdr)){
    localfdr[i]=mean(rankpost[1:i])
  }
  
  flocalfdr=rep(0,length(localfdr))
  rankp=rank(usepost,ties.method="random")
  flocalfdr=localfdr[rankp]
  finalsave=data.frame(finalgene,finalkclass,finalop,finalposter,finaltheta,flocalfdr)
  return(finalsave)
  
}


 


##emf7 estimate sigma1 and sigma2
emf7_node=function(zscore,kclass,theta,mu1,mu2,sigma1,sigma2){ ##give four initial values
  poster=rep(0,length(zscore))
  #  browser()
  for (iter in 1:1000){
    # Qiang Huang changed
    d1=dnorm(zscore,0,ifelse(sigma1==0,1,sigma1))*(1-theta[kclass])
    d2=dnorm(zscore,mu2,ifelse(sigma2==0,1,sigma2))*(theta[kclass]) 
    
    
    #d1=dnorm(zscore,0,sigma1)*(1-theta[kclass])
    #d2=dnorm(zscore,mu2,sigma2)*(theta[kclass])
    poster=d2/(d1+d2) ##step 1, calculate posterior
    
    mu2=sum(poster*zscore)/sum(poster) ##step 2 ,update mu
    
    if (mu2<1.25){
      mu2=1.25
    }
    
    
    sigmas1=sum((1-poster)*((zscore-0)*(zscore-0)))/sum((1-poster))
    sigma1=sqrt(sigmas1)
    sigmas2=sum(poster*((zscore-mu2)*(zscore-mu2)))/sum(poster)
    sigma2=sqrt(sigmas2)
    
    for (j in 1:length(theta)){
      theta[j]= mean(poster[kclass==j])
    }   
  }
  
  #  browser()
  out=list()
  out$mu1=mu1
  out$mu2=mu2
  out$theta=theta
  out$posteroir=poster 
  out$sigma1=sigma1
  out$sigma2=sigma2
  return(out)
}


stage1_analysis=function(expression_data,module_data, pvalue_data,cor_thres1,cor_thres2,flag_q,distgene){
  ##module pre_process
  mlabel=module_data[,2] ##first column is the gene name, second column is the module membership.
  genename=as.character(module_data[,1])
  # num_module=length(table(mlabel)) ##number of module # Qiang Huang changed
  num_module=max(mlabel)
  
  ##expression data pre_process
  genenamecor=as.character(rownames(expression_data))
  
  ##pvalue
  xingene=as.character(pvalue_data[,1]) 
  pv=pvalue_data[,2]
  genename2=xingene

  mresult=list() 
  for (modi in 1:num_module){ #!!!!
    genembert=genename[mlabel==modi]
    mgindex=which(genenamecor %in% genembert)
    printdoc=paste("Analyzing Module", modi, ",  ", "Number of Genes:",length(genembert), "")
    print(printdoc)
    
    # changed by Qiang for different distance mersures
    # corm=abs(cor(t(data[mgindex,]),use='pair',nThreads=3))
    corm <- dist_q(t(expression_data[mgindex,]),flag_q,distgene,genenamecor[mgindex])
    
    mglist=genenamecor[mgindex]
    mresult[[modi]]=modulef8(mgindex,genenamecor,genename2,pv,corm,1-cor_thres1,cor_thres2)
  }
  return(mresult)
}

dist_q <- function(datExpr,flag_q,distgene,subgenes_q){
    # combine networks or combine batch samples all the same
	#if(flag_q==1){
		#corm <- abs(cor(datExpr,use='pair',nThreads=3))   
	#}else if(flag_q==2){
        corm <- distgene[subgenes_q,subgenes_q]
	#}
	
	corm
}


organize_stage1_result=function(mresult,pv,nASDthres,module_data){
  # Qiang Huang changed
  num_module=max(module_data[,2])
  # num_module=length(table(module_data[,2])) ##number of module
  
  
  rorganize=list()
  modi=1
  thres=1.01
  rorganize[[modi]]=reportpf(mresult[[modi]]$graph,mresult[[modi]]$post,mresult[[modi]]$sub,mresult[[modi]]$pmeta,mresult[[modi]]$genelist,pv[mresult[[modi]]$glindex[mresult[[modi]]$glindex!=0]],mresult[[modi]]$Istart,mresult[[modi]]$Ifinal)
  genenum=dim(rorganize[[modi]])[1]
  keepi=c(1:genenum)[rorganize[[modi]][,3]<thres&is.na(rorganize[[modi]][,3])==F]
  keepi=c(1:genenum)[rorganize[[modi]][,3]<thres]
  if (length(keepi)>0){
    modulelabel=rep(modi,length(keepi))
    save1=data.frame(modulelabel,rorganize[[modi]][keepi,])
  }
  
  for (modi in 2:num_module){
    rorganize[[modi]]=reportpf(mresult[[modi]]$graph,mresult[[modi]]$post,mresult[[modi]]$sub,mresult[[modi]]$pmeta,mresult[[modi]]$genelist,pv[mresult[[modi]]$glindex[mresult[[modi]]$glindex!=0]],mresult[[modi]]$Istart,mresult[[modi]]$Ifinal)
    genenum=dim(rorganize[[modi]])[1]
    keepi=c(1:genenum)[rorganize[[modi]][,3]<thres|is.na(rorganize[[modi]][,3])==T]
    #  keepi=c(1:genenum)[rorganize[[modi]][,3]<thres]
    if (length(keepi)>0){
      modulelabel=rep(modi,length(keepi))
      newsave1=data.frame(modulelabel,rorganize[[modi]][keepi,])
      save1=rbind(save1,newsave1)
    }
  }
  
  keepgene=as.character(save1[,2])
  fixindex=rep(0,length(keepgene))
#  fixindex[which(keepgene %in% fixgene)]=1 ##tell which genes are fixed
  nASD=(save1[,4]<nASDthres)+0 ##tell if this gene is nASD genes
  save1=data.frame(save1,nASD,fixindex)
  
  
  save2=rep(0,11)
  for (i in 1:num_module){
    mstudy=mresult[[i]]
    nnode=length(mstudy$glindex)
    nsupern=length(mstudy$Istart)
    npstart=sum(mstudy$Istart)
    #  npfinal=sum(mstudy$Ifinal)
    npfinal=sum(mstudy$posterior>nASDthres,na.rm=T)
    #  nchange=sum(mstudy$Istart==0&mstudy$Ifinal==1)
    nchange=sum(mstudy$Istart==0&mstudy$posterior>nASDthres)
    b=mstudy$b
    c=mstudy$c
    mu=mstudy$mu
    sigma1=mstudy$sigmas1
    sigma2=mstudy$sigmas2
    nedge=sum(mstudy$graph)/nsupern
    msave=c(nnode,nsupern,npstart,npfinal,nchange,b,c,nedge,mu,sigma1,sigma2)
    save2=rbind(save2,msave) 
  }
  
  save2=save2[-1,]
  
  output=list()
  output$save1=save1
  output$save2=save2
  return(output)
}

stage2_analysis=function(stage1_orag_output,nASDthres,fdrthres,filename){
  save1=stage1_orag_output$save1
  save2=stage1_orag_output$save2
  
  keepgene=as.character(save1[,2])
  ##revise the posterior prob by adjusting the estimated parameters
  numnode=save2[,2]
  b=save2[,6]
  c=save2[,7]
  mu=save2[,9]
  sigmas=save2[,11]
  
  meannode=numnode/sum(numnode)
  usemu=sum(mu*meannode ,na.rm=T)
  usesigmas=sum(sigmas*meannode ,na.rm=T)
  
  pos_adjust=rep(0,dim(save1)[1])
  
  zscore=qnorm(1-save1$pbe)
  for (i in 1:length(pos_adjust)){
    if (is.na(save1$riskedge[i])==T){
      save1$riskedge[i]=0
    }
    new1=exp(b[save1$module[i]]+c[save1$module[i]]*save1$riskedge[i])
    new2=1
    p1=dnorm(zscore[i],usemu,sqrt(usesigmas))*new1/(new1+new2)  
    #  p1=dnorm(zscore[i],usemu,1)*new1/(new1+new2)
    p2=dnorm(zscore[i],0,1)*new2/(new1+new2)
    pos_adjust[i]=p2/(p1+p2)
  }
  nASD_adjust=(pos_adjust<0.5)+0
  
  numnodes=save2[,2]
  newclass=save1$geneclass
  for (i in 1:length(newclass)){
    if (save1$modulelabel[i]>1)
      newclass[i]=save1$geneclass[i]+sum(numnodes[1:(save1$modulelabel[i]-1)])
  }
  save1$geneclass=newclass
  save1$paf=pos_adjust
  save1$nASD=nASD_adjust
  
  # Qiang Huang changed
  if(any(table(save1$geneclass[which(save1$paf<nASDthres)]) > 9)){
      finalsave=stage2_mixturemodel(save1,nASDthres)
  }else{
      print("There is no stage 2 analysis!!!!")
      finalsave <- data.frame(finalgene=save1$keepgene,finalkclass=save1$geneclass,finalop=save1$keepgenep,finalposter=0,finaltheta=0,flocalfdr=0)
  }
  #finalsave=stage2_mixturemodel(save1,nASDthres)
  keepi=c(1:dim(finalsave)[1])[finalsave[,6]<1.1]
  genekeep=finalsave[keepi,1]
  fdr=finalsave[keepi,6]
  posterior=1-finalsave[keepi,4]
  op=finalsave[keepi,3]
  save3=data.frame(genekeep,fdr,posterior,op)
  
  fdrall=rep(NA,length(keepgene))
  post_stage2=rep(NA,length(keepgene))
  rASD=rep(0,length(keepgene))
  
  reorderi=rep(0,length(genekeep))
  for (i in 1:length(reorderi)){
    reorderi[i]=c(1:length(keepgene))[genekeep[i]==keepgene]
  }
  fdrall[reorderi]=fdr
  post_stage2[reorderi]=posterior
  rASD[reorderi]=(fdr<fdrthres)+0
  
  save3=data.frame(save1,fdrall,post_stage2,rASD)

  
  print("Ouputting DAWN analysis results")
  finalsave=data.frame(save1$keepgene,save1$nASD,rASD,
                       save1$keepgenep,save1$modulelabel,save1$geneclass, fdrall,
                       save1$paf, post_stage2,  
                       save1$riskedge,save1$totaledge)
  
  names(finalsave)=c("Gene","nASD","rASD","original_pvalue",
                     "Module","Class","FDR",
                     "Stage1_posterior","Stage2_posterior",
                     "risk_neighbor","total_neighbor")
  sortindex=order(fdrall)
  newfinalsave=finalsave[sortindex,]
  
  row.names(newfinalsave)=c(1:length(rASD))
  # write.csv(newfinalsave,"DAWN_analysis_result.csv") # Qiang Huang changed
  write.csv(newfinalsave,filename)
  
}

DAWN_package=function(expression_data,module_data, pvalue_data,cor_thres1=0.75,cor_thres2=0.7,nASDthres=0.5,fdrthres=0.1,flag_q,filename,distgene){
  print("Starting Stage I Analysis")
  mresult=stage1_analysis(expression_data,module_data, pvalue_data,cor_thres1,cor_thres2,flag_q,distgene)
  pv=pvalue_data[,2]
  stage1_orag_output=organize_stage1_result(mresult,pv,nASDthres,module_data)
  print("Starting Stage II Analysis")
  stage2_analysis(stage1_orag_output,nASDthres,fdrthres,filename)
}

