args <- commandArgs(T)
kk <- as.numeric(args[1])
	
	print(kk)
	source("CRF_build.R")
    source("MixNet.R")
    source("Network_analysis.R")
    source("Multi_net.R")
    options(stringsAsFactors=FALSE,expressions=500000)
    fileexp <- "PCGCall.txt" 
    netstr <- c("STRING/","iRef/","Infer/","coexp/")
    modulefiles <- c("module3","iRefmodule","DAWNmodule","DAWNmodule")
    
    betaV <- c(0,0.1,0.2,0.3,0.4,0.5,1,2,5,10)
    
    beta <- betaV[kk]
    	    
    	k = 1 
    	netflag=3
    		
			if(netflag==3){ instr <- "CRF_input";idm=TRUE; }else{ instr <- "hotnet_input"; idm=FALSE;}
    		modulefile <- modulefiles[k]
    
    		dirstr <- "result/randset4/"
    		strname <- "ASD4_16"
    		filename <- paste(dirstr,instr,strname,".txt",sep="")
    		strn <- paste("result/randresult4/",netstr[k],"CRFresult_",beta,strname,sep="")
    		Multi_net(netflag,strn,filename,fileexp,modulefile,cutn=100000,idm,beta=beta);
 


