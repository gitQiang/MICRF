control_case_meta <- function(kk){

	source("CRF_build.R")
    source("MixNet.R")
    source("Network_analysis.R")
    source("Multi_net.R")
    options(stringsAsFactors=FALSE)
    fileexp <- "PCGCall.txt" 
    netstr <- c("STRING/","iRef/","coexp/","HPRD/","coexp1/","corr/","coexp5/","coexp7/")
    modulefiles <- c("module3","iRefmodule","DAWNmodule","HPRDmodule","DAWNmodule","module_coexp","module_coexp","module_coexp")
     
    betaV <- c(0,0.2,0.4,0.5,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0)
    beta <- betaV[kk]
    
    for(net_e in 1:4){
    	dirfold <- paste("result/control/v",net_e,"/",sep="")
    	k = 1 
    	dirstr <- "result/"
    	for(netflag in c(3,6,7,20,21,26,27,28)){
    		
		 if(netflag==3){ instr <- "CRF_input";idm=TRUE; }else{ instr <- "hotnet_input"; idm=FALSE;}
    		modulefile <- modulefiles[k]

#     		strname <- "control"
#             filename <- paste(dirstr,"control/",instr,strname,".txt",sep="")
#            	strn <- paste(dirfold,netstr[k],"CRFresult_",beta,strname,sep="")
#             Multi_net(netflag,strn,filename,fileexp,modulefile,cutn=100000,idm,beta=beta,net_e=net_e);    
#             
#             strname <- "case"
#             filename <- paste(dirstr,instr,strname,".txt",sep="")
#            	strn <- paste(dirfold,netstr[k],"CRFresult_",beta,strname,sep="")
#             Multi_net(netflag,strn,filename,fileexp,modulefile,cutn=100000,idm,beta=beta,net_e=net_e);  
            
#             strname <- "meta"
#             filename <- paste(dirstr,instr,strname,".txt",sep="")
#            	strn <- paste(dirfold,netstr[k],"CRFresult_",beta,strname,sep="")
#             Multi_net(netflag,strn,filename,fileexp,modulefile,cutn=100000,idm,beta=beta,net_e=net_e); 

#             strname <- "meta_mis"
#             filename <- paste(dirstr,instr,strname,".txt",sep="")
#             strn <- paste(dirfold,netstr[k],"CRFresult_",beta,strname,sep="")
#             Multi_net(netflag,strn,filename,fileexp,modulefile,cutn=100000,idm,beta=beta,net_e=net_e);  
         
# 		    strname <- "nat13772"
# 		    filename <- paste(dirstr,instr,strname,".txt",sep="")
# 		    strn <- paste(dirfold,netstr[k],"CRFresult_",beta,strname,sep="")
# 		    Multi_net(netflag,strn,filename,fileexp,modulefile,cutn=100000,idm,beta=beta,net_e=net_e);  
#          
# 		    strname <- "nat13772_0"
# 		    filename <- paste(dirstr,instr,strname,".txt",sep="")
# 		    strn <- paste(dirfold,netstr[k],"CRFresult_",beta,strname,sep="")
# 		    Multi_net(netflag,strn,filename,fileexp,modulefile,cutn=100000,idm,beta=beta,net_e=net_e); 
            
            strname <- "repDDD"
            filename <- paste(dirstr,instr,strname,".txt",sep="")
            strn <- paste(dirfold,netstr[k],"CRFresult_",beta,strname,sep="")
            Multi_net(netflag,strn,filename,fileexp,modulefile,cutn=100000,idm,beta=beta,net_e=net_e); 

    		k <- k + 1
		}
	}
		
}

batch_randset4 <- function(i){
    
    source("CRF_build.R")
    source("MixNet.R")
    source("Network_analysis.R")
    source("Multi_net.R")
    options(stringsAsFactors=FALSE)
    fileexp <- "PCGCall.txt" 
    netstr <- c("STRING/","iRef/","coexp/","HPRD/","coexp1/","corr/","coexp5/","coexp7/")
    modulefiles <- c("module3","iRefmodule","DAWNmodule","HPRDmodule","DAWNmodule","module_coexp","module_coexp","module_coexp")
    
    betaM=matrix(c(7.1162,6.5513),8,2,byrow=TRUE)
    dirstr <- "result/randset4_2/"
    net_e=5
    dirfold <- paste("result/randresult4_",net_e,"/",sep="")
        
    k = 1
    for(netflag in c(3,6,7,20,21,26,27,28)){   
        if(netflag==3){ instr <- "CRF_input";idm=TRUE; }else{ instr <- "hotnet_input"; idm=FALSE;}
        modulefile <- modulefiles[k]
        
        beta <- betaM[k,]
        ## random set 4: different sample size and more samples
        j=3
        strname <- paste("part",j,"_",i,sep="")
        filename <- paste(dirstr,instr,strname,".txt",sep="")
        strn <- paste(dirfold,netstr[k],"CRFresult_",netflag,strname,sep="")
        Multi_net(netflag,strn,filename,fileexp,modulefile,cutn=100000,idm,beta=beta,net_e=net_e);       
        
        k <- k + 1
    } 

}

batch_leaveone4 <- function(i){
    
    source("CRF_build.R")
    source("MixNet.R")
    source("Network_analysis.R")
    source("Multi_net.R")
    options(stringsAsFactors=FALSE,expressions=500000)
    fileexp <- "PCGCall.txt" 
    netstr <- c("STRING/","iRef/","coexp/","HPRD/","coexp1/","corr/","coexp5/","coexp7/")
    modulefiles <- c("module3","iRefmodule","DAWNmodule","HPRDmodule","DAWNmodule","module_coexp","module_coexp","module_coexp")
    
    betaM=matrix(c(7.1162,6.5513),8,2,byrow=TRUE)
    dirstr <- "result/leaveone4_2/"
    net_e=5
    dirfold=paste("result/leaveone4result_",net_e,"/",sep="")

    k = 1
    for(netflag in c(3,6,7,20,21,26,27,28)){
        if(netflag==3){ instr <- "CRF_input";idm=TRUE; }else{ instr <- "hotnet_input"; idm=FALSE;}
        modulefile <- modulefiles[k]
        
        beta <- betaM[k,]
        strname <- paste("rand2_",i,sep="")
        filename <- paste(dirstr,instr,strname,".txt",sep="")
        strn <- paste(dirfold,netstr[k],"CRFresult_",netflag,strname,sep="")
        Multi_net(netflag,strn,filename,fileexp,modulefile,cutn=100000,idm,beta=beta,net_e=net_e);      
        
        k <- k + 1
    }
    
}

batch_randset_1_faster <- function(k){
    
    source("CRF_build.R")
    source("MixNet.R")
    source("Network_analysis.R")
    source("Multi_net.R")
    fileexp <- "PCGCall.txt" 
    netstr <- c("STRING/","iRef/","coexp/","HPRD/","coexp1/","corr/","coexp5/","coexp7/")
    modulefiles <- c("module3","iRefmodule","DAWNmodule","HPRDmodule","DAWNmodule","module_coexp","module_coexp","module_coexp")
    
    ## bestP
    betaM=matrix(c(7.1162,6.5513),8,2,byrow=TRUE)
    netfs <- c(3,6,7,20,21,26,27,28)
    dirfold <- "result/randresult_5/"
    dirstr <- "result/randset_1_2/"
    
    for(netk in 1:8){
        netflag <- netfs[netk]
        net_e <- 5
        beta <- BetaM[netk,]
        
        if(netflag==3){ instr <- "CRF_input";idm=TRUE; }else{ instr <- "hotnet_input"; idm=FALSE;}
        modulefile <- modulefiles[netk]
        
        ## random set 1
        i <- k %% 20
        if(i==0) i=20
        j <- floor((k-1)/20) + 2
        
        strname <- paste("rand1",j,"_",i,sep="")
        filename <- paste(dirstr,instr,strname,".txt",sep="")
        strn <- paste(dirfold,netstr[netk],"CRFresult_",netflag,strname,sep="")
        Multi_net(netflag,strn,filename,fileexp,modulefile,cutn=100000,idm,beta=beta,net_e=net_e);       
        
    }
    
}

batch_ASD_randset4 <- function(kk){
	
    source("CRF_build.R")
    source("MixNet.R")
    source("Network_analysis.R")
    source("Multi_net.R")
    options(stringsAsFactors=FALSE,expressions=500000)
    fileexp <- "PCGCall.txt" 
    netstr <- c("STRING/","iRef/","coexp/","HPRD/","coexp1/","corr/","coexp5/","coexp7/")
    modulefiles <- c("module3","iRefmodule","DAWNmodule","HPRDmodule","DAWNmodule","module_coexp","module_coexp","module_coexp")

    #betaV <- c(0.1,0.5,1,2,3,4,5,6,7,10,100,1000)
    
    betaV <- c(0,0.2,0.4,0.5,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0)
    beta <- betaV[kk]
    dirstr <- "result/randset4/"
    
    for(net_e in 1:4){
    #net_e=4
    dirfold <- paste("result/randresult4_",net_e,"/",sep="")
    
    k = 1
    for(netflag in c(3,6,7,20,21,26,27,28)){   
        if(netflag==3){ instr <- "CRF_input";idm=TRUE; }else{ instr <- "hotnet_input"; idm=FALSE;}
        modulefile <- modulefiles[k]
                
        ## random set 4: different sample size and more samples
        j=3
        for(i in 1:20){
            strname <- paste("part",j,"_",i,sep="")
            filename <- paste(dirstr,instr,strname,".txt",sep="")
            strn <- paste(dirfold,netstr[k],"CRFresult_",beta,strname,sep="")
            Multi_net(netflag,strn,filename,fileexp,modulefile,cutn=100000,idm,beta=beta,net_e=net_e);       
        }
        
        k <- k + 1
    } 
    
    }

}

batch_ASD_leaveone4 <- function(i){
	
    source("CRF_build.R")
    source("MixNet.R")
    source("Network_analysis.R")
    source("Multi_net.R")
    options(stringsAsFactors=FALSE,expressions=500000)
    fileexp <- "PCGCall.txt" 
    netstr <- c("STRING/","iRef/","coexp/","HPRD/","coexp1/","corr/","coexp5/","coexp7/")
    modulefiles <- c("module3","iRefmodule","DAWNmodule","HPRDmodule","DAWNmodule","module_coexp","module_coexp","module_coexp")
    
    betaV <- c(0,0.2,0.4,0.5,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0)
    
    dirstr <- "result/leaveone4/"
    
    for(net_e in 1:4){
    dirfold=paste("result/leaveone4result_",net_e,"/",sep="")
    
    for(kk in 1:length(betaV)){
    	beta <- betaV[kk]
    	
        k = 1
    	for(netflag in c(3,6,7,20,21,26,27,28)){
		    if(netflag==3){ instr <- "CRF_input";idm=TRUE; }else{ instr <- "hotnet_input"; idm=FALSE;}
    		modulefile <- modulefiles[k]
    		
            strname <- paste("rand2_",i,sep="")
            filename <- paste(dirstr,instr,strname,".txt",sep="")
           	strn <- paste(dirfold,netstr[k],"CRFresult_",beta,strname,sep="")
            Multi_net(netflag,strn,filename,fileexp,modulefile,cutn=100000,idm,beta=beta,net_e=net_e);      
    		
    		k <- k + 1
		}
	}
    
    }

}

batch_ASD_randset_1 <- function(j){
    
    source("CRF_build.R")
    source("MixNet.R")
    source("Network_analysis.R")
    source("Multi_net.R")
    fileexp <- "PCGCall.txt" 
    netstr <- c("STRING/","iRef/","coexp/","HPRD/","coexp1/","corr/","coexp5/","coexp7/")
    modulefiles <- c("module3","iRefmodule","DAWNmodule","HPRDmodule","DAWNmodule","module_coexp","module_coexp","module_coexp")
    
    ## bestP
    load("BestP_8_19_4") 
    BestP[1,1] <- 3
    dirfold <- "result/randresult_1/"
    dirstr <- "result/randset_1/"
    
    for(netk in 1:8){
        netflag <- BestP[netk,1]
        net_e <- BestP[netk,2]
        beta <- BestP[netk,3]
        
        if(netflag==3){ instr <- "CRF_input";idm=TRUE; }else{ instr <- "hotnet_input"; idm=FALSE;}
        modulefile <- modulefiles[netk]
                
        ## random set 1
        for(i in 1:20){
            strname <- paste("rand1",j,"_",i,sep="")
            filename <- paste(dirstr,instr,strname,".txt",sep="")
            strn <- paste(dirfold,netstr[netk],"CRFresult_",beta,strname,sep="")
            Multi_net(netflag,strn,filename,fileexp,modulefile,cutn=100000,idm,beta=beta,net_e=net_e);       
        }       
    }
    
}

batch_ASD_randset_1_faster <- function(k){
    
    source("CRF_build.R")
    source("MixNet.R")
    source("Network_analysis.R")
    source("Multi_net.R")
    fileexp <- "PCGCall.txt" 
    netstr <- c("STRING/","iRef/","coexp/","HPRD/","coexp1/","corr/","coexp5/","coexp7/")
    modulefiles <- c("module3","iRefmodule","DAWNmodule","HPRDmodule","DAWNmodule","module_coexp","module_coexp","module_coexp")
    
    ## bestP
    load("BestP_8_19_4") 
    BestP[1,1] <- 3
    dirfold <- "result/randresult_1/"
    dirstr <- "result/randset_1/"
    
    for(netk in 1:8){
        netflag <- BestP[netk,1]
        net_e <- BestP[netk,2]
        beta <- BestP[netk,3]
        
        if(netflag==3){ instr <- "CRF_input";idm=TRUE; }else{ instr <- "hotnet_input"; idm=FALSE;}
        modulefile <- modulefiles[netk]
        
        ## random set 1
        i <- k %% 20
        if(i==0) i=20
        j <- floor((k-1)/20) + 2
            
        strname <- paste("rand1",j,"_",i,sep="")
        filename <- paste(dirstr,instr,strname,".txt",sep="")
        strn <- paste(dirfold,netstr[netk],"CRFresult_",beta,strname,sep="")
        Multi_net(netflag,strn,filename,fileexp,modulefile,cutn=100000,idm,beta=beta,net_e=net_e);       
               
    }
    
}

hotnet2_control <- function(){
    
#     source("CRF_build.R")
#     source("HotNet2.R")
#     source("Network_analysis.R")
#     source("Multi_net.R")
#     options(stringsAsFactors=FALSE)
#     fileexp=""
#     netstr <- c("STRING/","iRef/","coexp/","HPRD/","coexp1/","corr/","coexp5/","coexp7/")
#     k = 1
#     for(netflag in c(3,6,7,20,21,26,27,28)){
#         
#         if(netflag==3){ instr <- "CRF_input";idm=TRUE; }else{ instr <- "hotnet_input"; idm=FALSE;}
#         dirin <- "result/control/"
#         dirstr <- paste("result/control/hotnet/",netstr[k],sep="")
#         
#         strname <- "control"
#         filename <- paste(dirin,instr,strname,".txt",sep="")
#         HotNet2R(netflag,strname,filename,fileexp,dirstr,cutn=0.50,idm)
#      
#         k <- k +1
#     }
    
#     ### nature 13772
#     source("CRF_build.R")
#     source("HotNet2.R")
#     source("Network_analysis.R")
#     source("Multi_net.R")
#     options(stringsAsFactors=FALSE)
#     fileexp=""
#     netstr <- c("STRING/","iRef/","coexp/","HPRD/","coexp1/","corr/","coexp5/","coexp7/")
#     k = 1
#     for(netflag in c(3,6,7,20,21,26,27,28)){
#         
#         if(netflag==3){ instr <- "CRF_input";idm=TRUE; }else{ instr <- "hotnet_input"; idm=FALSE;}
#         dirin <- "result/"
#         dirstr <- paste("result/control/hotnet/",netstr[k],sep="")
# 
#         strname <- "nat13772_0"
#         filename <- paste(dirin,instr,strname,".txt",sep="")
#         HotNet2R(netflag,strname,filename,fileexp,dirstr,cutn=0.50,idm)
#         
#         strname <- "nat13772"
#         filename <- paste(dirin,instr,strname,".txt",sep="")
#         HotNet2R(netflag,strname,filename,fileexp,dirstr,cutn=0.50,idm)
#         
#         k <- k +1
#     }
    
#     ### meta analysis 
#     source("CRF_build.R")
#     source("HotNet2.R")
#     source("Network_analysis.R")
#     source("Multi_net.R")
#     options(stringsAsFactors=FALSE)
#     fileexp=""
#     netstr <- c("STRING/","iRef/","coexp/","HPRD/","coexp1/","corr/","coexp5/","coexp7/")
#     k = 1
#     for(netflag in c(3,6,7,20,21,26,27,28)){
#         
#         if(netflag==3){ instr <- "CRF_input";idm=TRUE; }else{ instr <- "hotnet_input"; idm=FALSE;}
#         dirin <- "result/"
#         dirstr <- paste("result/control/hotnet/",netstr[k],sep="")
#         
#         strname <- "meta"
#         filename <- paste(dirin,instr,strname,".txt",sep="")
#         HotNet2R(netflag,strname,filename,fileexp,dirstr,cutn=0.50,idm)
#                 
#         k <- k +1
#     }
#     
#     
    
    ### repDDD 1984
    source("CRF_build.R")
    source("HotNet2.R")
    source("Network_analysis.R")
    source("Multi_net.R")
    options(stringsAsFactors=FALSE)
    fileexp=""
    netstr <- c("STRING/","iRef/","coexp/","HPRD/","coexp1/","corr/","coexp5/","coexp7/")
    k = 1
    for(netflag in c(3,6,7,20,21,26,27,28)){
        
        if(netflag==3){ instr <- "CRF_input";idm=TRUE; }else{ instr <- "hotnet_input"; idm=FALSE;}
        dirin <- "result/"
        dirstr <- paste("result/control/hotnet/",netstr[k],sep="")
        
        strname <- "repDDD"
        filename <- paste(dirin,instr,strname,".txt",sep="")
        HotNet2R(netflag,strname,filename,fileexp,dirstr,cutn=0.50,idm)
       
        k <- k +1
    }
    
}

hotnet2_randset_1 <- function(i){
    
    source("CRF_build.R")
    source("HotNet2.R")
    source("Network_analysis.R")
    source("Multi_net.R")
    options(stringsAsFactors=FALSE)
    fileexp=""
    netstr <- c("STRING/","iRef/","coexp/","HPRD/","coexp1/","corr/","coexp5/","coexp7/")
    #k =4
    #for(netflag in c(20,21)){
    k = 1
    for(netflag in c(3,6,7,20,21,26,27,28)){
        if(netflag==3){ instr <- "CRF_input";idm=TRUE; }else{ instr <- "hotnet_input"; idm=FALSE;}
                         
        dirin <- "result/randset_1/"
        dirstr <- paste("result/randresult_1/",netstr[k],sep="")
   
        ## random set 1: different sample size and exclude sample sets
        for(j in 2:9){
                strname <- paste("rand1",j,"_",i,sep="")
                filename <- paste(dirin,instr,strname,".txt",sep="")
                HotNet2R(netflag,strname,filename,fileexp,dirstr,cutn=0.50,idm)
        }
        
        k <- k + 1
     }
       
}

hotnet2_leaveone4 <- function(i){

    source("CRF_build.R")
    source("HotNet2.R")
    source("Network_analysis.R")
    source("Multi_net.R")
    options(stringsAsFactors=FALSE)
    fileexp=""
    netstr <- c("STRING/","iRef/","coexp/","HPRD/","coexp1/","corr/","coexp5/","coexp7/")
    #k=4
    #for(netflag in c(20,21)){
    k = 1
    for(netflag in c(3,6,7,20,21,26,27,28)){
        
        if(netflag==3){ instr <- "CRF_input";idm=TRUE; }else{ instr <- "hotnet_input"; idm=FALSE;}
                         
        ### random set 2: leaveone mutation
        dirin <- "result/leaveone4/"
        dirstr <- paste("result/leaveone4result/",netstr[k],sep="")

        strname <- paste("rand2_",i,sep="")
        filename <- paste(dirin,instr,strname,".txt",sep="")
        HotNet2R(netflag,strname,filename,fileexp,dirstr,cutn=0.50,idm)
	k <- k + 1
	} 

}

hotnet2_randset4 <- function(i){
    source("CRF_build.R")
    source("HotNet2.R")
    source("Network_analysis.R")
    source("Multi_net.R")
    options(stringsAsFactors=FALSE)
    fileexp=""
    netstr <- c("STRING/","iRef/","coexp/","HPRD/","coexp1/","corr/","coexp5/","coexp7/")
    #k=4
    #for(netflag in c(20,21)){
    k = 1
    for(netflag in c(3,6,7,20,21,26,27,28)){
        
        if(netflag==3){ instr <- "CRF_input";idm=TRUE; }else{ instr <- "hotnet_input"; idm=FALSE;}
      
        ## random set 3: rand set 4
        dirin <- "result/randset4/"
        dirstr <- paste("result/randresult4/",netstr[k],sep="")
        
        ## random set 3: different sample size and exclude sample sets
        for(j in 3){
                strname <- paste("part",j,"_",i,sep="")
                filename <- paste(dirin,instr,strname,".txt",sep="")
                HotNet2R(netflag,strname,filename,fileexp,dirstr,cutn=0.50,idm)
        }
        
        k <- k +1
        }
}
