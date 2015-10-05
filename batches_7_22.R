batch_ASD_beta <- function(kk){
	
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
    	for(netflag in c(3,6,8,7)){
    		
			if(netflag==3){ instr <- "CRF_input";idm=TRUE; }else{ instr <- "hotnet_input"; idm=FALSE;}
    		modulefile <- modulefiles[k]
    
    		dirstr <- "result/randset4/"
    		strname <- "ASD4_16"
    		filename <- paste(dirstr,instr,strname,".txt",sep="")
    		strn <- paste("result/randresult4/",netstr[k],"CRFresult_",beta,strname,sep="")
    		Multi_net(netflag,strn,filename,fileexp,modulefile,cutn=20000,idm,beta=beta);
 
    		## random set 4: different sample size and more samples
    		j=3
    		for(i in 1:20){
            	strname <- paste("part",j,"_",i,sep="")
            	filename <- paste(dirstr,instr,strname,".txt",sep="")
           		strn <- paste("result/randresult4/",netstr[k],"CRFresult_",beta,strname,sep="")
            	Multi_net(netflag,strn,filename,fileexp,modulefile,cutn=20000,idm,beta=beta);       
            }
            
    		k <- k + 1
		} 

}

batch_ASD_beta1 <- function(i){
	
	source("CRF_build.R")
    source("MixNet.R")
    source("Network_analysis.R")
    source("Multi_net.R")
    options(stringsAsFactors=FALSE,expressions=500000)
    fileexp <- "PCGCall.txt" 
    netstr <- c("STRING/","iRef/","Infer/","coexp/")
    modulefiles <- c("module3","iRefmodule","DAWNmodule","DAWNmodule")
    
    betaV <- c(0,0.1,0.2,0.3,0.4,0.5,1,2,5,10)
    dirstr <- "result/leaveone4/"
    
    for(kk in 1:10){
    	beta <- betaV[kk]
    	k = 1 
    	for(netflag in c(3,6,8,7)){
			if(netflag==3){ instr <- "CRF_input";idm=TRUE; }else{ instr <- "hotnet_input"; idm=FALSE;}
    		modulefile <- modulefiles[k]
    		
            	strname <- paste("rand2_",i,sep="")
            	filename <- paste(dirstr,instr,strname,".txt",sep="")
           		strn <- paste("result/leaveone4result/",netstr[k],"CRFresult_",beta,strname,sep="")
            	Multi_net(netflag,strn,filename,fileexp,modulefile,cutn=20000,idm,beta=beta);      
    			k <- k + 1
		}
	} 

}

control_case_meta <- function(kk){

	source("CRF_build.R")
    source("MixNet.R")
    source("Network_analysis.R")
    source("Multi_net.R")
    options(stringsAsFactors=FALSE)
    fileexp <- "PCGCall.txt" 
    netstr <- c("STRING/","iRef/","Infer/","coexp/","GNAT/")
    modulefiles <- c("module3","iRefmodule","DAWNmodule","DAWNmodule","data/GNATnet/module_brain")
     
    betaV <- c(0,0.2,0.4,0.5,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0)
    beta <- betaV[kk]
    
    for(net_e in 1:4){
    	dirfold <- paste("result/control/v",net_e,"/",sep="")
    	k = 1 
    	dirstr <- "result/control/"
    	for(netflag in c(3,6,8,7,21)){
    		
		 if(netflag==3){ instr <- "CRF_input";idm=TRUE; }else{ instr <- "hotnet_input"; idm=FALSE;}
    		modulefile <- modulefiles[k]

#     		strname <- "control"
#             filename <- paste(dirstr,instr,strname,".txt",sep="")
#            	strn <- paste(dirfold,netstr[k],"CRFresult_",beta,strname,sep="")
#             Multi_net(netflag,strn,filename,fileexp,modulefile,cutn=20000,idm,beta=beta,net_e=net_e);    
#             
#             strname <- "case"
#             filename <- paste(dirstr,instr,strname,".txt",sep="")
#            	strn <- paste(dirfold,netstr[k],"CRFresult_",beta,strname,sep="")
#             Multi_net(netflag,strn,filename,fileexp,modulefile,cutn=20000,idm,beta=beta,net_e=net_e);  
#             
#             strname <- "Meta"
#             filename <- paste(dirstr,instr,strname,".txt",sep="")
#            	strn <- paste(dirfold,netstr[k],"CRFresult_",beta,strname,sep="")
#             Multi_net(netflag,strn,filename,fileexp,modulefile,cutn=20000,idm,beta=beta,net_e=net_e);  
         
		    strname <- "ASDnat13772"
		    filename <- paste(dirstr,instr,strname,".txt",sep="")
		    strn <- paste(dirfold,netstr[k],"CRFresult_",beta,strname,sep="")
		    Multi_net(netflag,strn,filename,fileexp,modulefile,cutn=20000,idm,beta=beta,net_e=net_e);  
            
    		k <- k + 1
		}
	}
		
}

nat13772 <- function(){
    
    source("CRF_build.R")
    source("MixNet.R")
    source("Network_analysis.R")
    source("Multi_net.R")
    options(stringsAsFactors=FALSE)
    fileexp <- "PCGCall.txt" 
    netstr <- c("STRING/","iRef/","Infer/","coexp/","GNAT/")
    modulefiles <- c("module3","iRefmodule","DAWNmodule","DAWNmodule","data/GNATnet/module_brain")
    
    kk=2
    betaV <- c(0,0.2,0.4,0.5,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0)
    beta <- betaV[kk]
    
    for(net_e in 4){
        dirfold <- paste("result/control/v",net_e,"/",sep="")
        k = 4 
        dirstr <- "result/control/"
        for(netflag in 7){
            
            if(netflag==3){ instr <- "CRF_input";idm=TRUE; }else{ instr <- "hotnet_input"; idm=FALSE;}
            modulefile <- modulefiles[k]
                       
            strname <- "nat13772_0"
            filename <- paste(dirstr,instr,strname,".txt",sep="")
            strn <- paste(dirfold,netstr[k],"CRFresult_",beta,strname,sep="")
            Multi_net(netflag,strn,filename,fileexp,modulefile,cutn=20000,idm,beta=beta,net_e=net_e);  
            
            k <- k + 1
        }
    }
    
}

batch_ASD_randset4 <- function(kk){
	
    source("CRF_build.R")
    source("MixNet.R")
    source("Network_analysis.R")
    source("Multi_net.R")
    options(stringsAsFactors=FALSE,expressions=500000)
    fileexp <- "PCGCall.txt" 
    netstr <- c("STRING/","iRef/","coexp/","HPRD/","coexp1")
    modulefiles <- c("module3","iRefmodule","DAWNmodule","HPRDmodule","DAWNmodule")

    #betaV <- c(0.1,0.5,1,2,3,4,5,6,7,10,100,1000)
    
    betaV <- c(0,0.2,0.4,0.5,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0)
    beta <- betaV[kk]
    	    
    net_e=4
    dirfold <- "result/randresult4_4/"
    dirstr <- "result/randset4/"
    
    k=4
    for(netflag in c(20,21)){  
    #k = 1
    #for(netflag in c(3,6,7,20,21)){   
        if(netflag==3){ instr <- "CRF_input";idm=TRUE; }else{ instr <- "hotnet_input"; idm=FALSE;}
        modulefile <- modulefiles[k]
        
        strname <- "ASD4_16"
        filename <- paste(dirstr,instr,strname,".txt",sep="")
        strn <- paste(dirfold,netstr[k],"CRFresult_",beta,strname,sep="")
        Multi_net(netflag,strn,filename,fileexp,modulefile,cutn=20000,idm,beta=beta,net_e=net_e);
        
        ## random set 4: different sample size and more samples
        j=3
        for(i in 1:20){
            strname <- paste("part",j,"_",i,sep="")
            filename <- paste(dirstr,instr,strname,".txt",sep="")
            strn <- paste(dirfold,netstr[k],"CRFresult_",beta,strname,sep="")
            Multi_net(netflag,strn,filename,fileexp,modulefile,cutn=20000,idm,beta=beta,net_e=net_e);       
        }
        
        k <- k + 1
    } 

}

batch_ASD_leaveone4 <- function(i){
	
    source("CRF_build.R")
    source("MixNet.R")
    source("Network_analysis.R")
    source("Multi_net.R")
    options(stringsAsFactors=FALSE,expressions=500000)
    fileexp <- "PCGCall.txt" 
    netstr <- c("STRING/","iRef/","coexp/","HPRD/","coexp1/")
    modulefiles <- c("module3","iRefmodule","DAWNmodule","HPRDmodule","DAWNmodule")
    
    betaV <- c(0,0.2,0.4,0.5,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0)
    dirstr <- "result/leaveone4/"
    
    net_e=4
    dirfold="result/leaveone4result_4/"
    
    for(kk in 1:length(betaV)){
    	beta <- betaV[kk]
    	
        #k = 1
    	#for(netflag in c(3,6,7,20,21)){
    	k=4
        for(netflag in c(20,21)){
		    if(netflag==3){ instr <- "CRF_input";idm=TRUE; }else{ instr <- "hotnet_input"; idm=FALSE;}
    		modulefile <- modulefiles[k]
    		
            strname <- paste("rand2_",i,sep="")
            filename <- paste(dirstr,instr,strname,".txt",sep="")
           	strn <- paste(dirfold,netstr[k],"CRFresult_",beta,strname,sep="")
            Multi_net(netflag,strn,filename,fileexp,modulefile,cutn=20000,idm,beta=beta,net_e=net_e);      
    		
    		k <- k + 1
		}
	} 

}

test_ori <- function(i){
	
	source("CRF_build.R")
    source("MixNet.R")
    source("Network_analysis.R")
    source("Multi_net.R")
    options(stringsAsFactors=FALSE,expressions=500000)
    fileexp <- "PCGCall.txt" 
    netstr <- c("STRING/","iRef/","Infer/","coexp/")
    modulefiles <- c("module3","iRefmodule","DAWNmodule","DAWNmodule")
    
    betaV <- c(0,0.2,0.4,0.5,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0)
    dirstr <- "result/leaveone4/"
    
    net_e=FALSE
    dirfold="result/leaveone4result_1/"
    kk=4
    beta <- betaV[kk]
    
    	k = 1 
    	for(netflag in c(3,6,8,7)){
    		
			if(netflag==3){ instr <- "CRF_input";idm=TRUE; }else{ instr <- "hotnet_input"; idm=FALSE;}
    		modulefile <- modulefiles[k]

            strname <- paste("rand2_",i,sep="")
            filename <- paste(dirstr,instr,strname,".txt",sep="")
           	strn <- paste(dirfold,netstr[k],"CRFresult_",beta,strname,sep="")
            Multi_net(netflag,strn,filename,fileexp,modulefile,cutn=20000,idm,beta=beta,net_e=net_e);      
    		k <- k + 1
		}

}

batch1_ASD <- function(i){
    
    source("CRF_build.R")
    source("MixNet.R")
    source("Network_analysis.R")
    source("Multi_net.R")
    options(stringsAsFactors=FALSE)
    fileexp <- "PCGCall.txt" 
    
    netstr <- c("Infer/","coexp/")
    k = 1 
    modulefiles <- c("DAWNmodule","DAWNmodule")
    
for(netflag in c(8,7)){
    
   		if(netflag==3){ instr <- "CRF_input";idm=TRUE; }else{ instr <- "hotnet_input"; idm=FALSE;}
    	modulefile <- modulefiles[k]
    
    	dirstr <- "result/leaveone3/"
        strname <- paste("rand2_de1_",i,sep="")
        filename <- paste(dirstr,instr,"rand2_de1_",i,".txt",sep="")
        strn <- paste("result/leaveone3result/",netstr[k],"CRFresult_",strname,sep="")
        Multi_net(netflag,strn,filename,fileexp,modulefile,cutn=20000,idm); 
    
    	k <- k + 1
} 

}

batch2_ASD <- function(i){
    
   source("CRF_build.R")
    source("MixNet.R")
    source("Network_analysis.R")
    source("Multi_net.R")
    options(stringsAsFactors=FALSE)
    fileexp <- "PCGCall.txt" 
    
    netstr <- c("Overlapexp/","Jointexp/","Overlapexp_PPI/","Jointexp_PPI/","Overlap/")
    k = 1
    modulefiles <- c("DAWNmodule","DAWNmodule","data/network_inference/Overlap_O4_top5_PPImodule","data/network_inference/Joint_O4_top5_PPImodule","data/network_inference/OverlapPPImodule")
    dirstr <- "result/randset3/"
for(netflag in c(15:19)){
    
    if(netflag==3){ instr <- "CRF_input";idm=TRUE; }else{ instr <- "hotnet_input"; idm=FALSE;}
    modulefile <- modulefiles[k]
   
    ## random set 1: different sample size and exclude sample sets
    for(j in 2:9){
            strname <- paste("part",j,"_",i,sep="")
            filename <- paste(dirstr,instr,strname,".txt",sep="")
            strn <- paste("result/randresult3/",netstr[k],"CRFresult_",strname,sep="")
            Multi_net(netflag,strn,filename,fileexp,modulefile,cutn=20000,idm);     
       }
       
    k <- k + 1
} # end for netflag

}

batch_ASD_randset_1 <- function(j){
    
    source("CRF_build.R")
    source("MixNet.R")
    source("Network_analysis.R")
    source("Multi_net.R")
    fileexp <- "PCGCall.txt" 
    netstr <- c("STRING/","iRef/","coexp/","HPRD/","coexp1/")
    modulefiles <- c("module3","iRefmodule","DAWNmodule","DAWNmodule","data/GNATnet/module_brain")
    
    ## bestP
    load("BestP_6_10")
    BestP[1,1] <- 3	  
    
    dirfold <- "result/randresult_1/"
    dirstr <- "result/randset_1/"
    
    for(netk in 1:5){
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
            Multi_net(netflag,strn,filename,fileexp,modulefile,cutn=20000,idm,beta=beta,net_e=net_e);       
        }       
    }
    
}

hotnet2_control <- function(){
    
#     source("CRF_build.R")
#     source("HotNet2.R")
#     source("Network_analysis.R")
#     source("Multi_net.R")
#     options(stringsAsFactors=FALSE)
#     fileexp=""
#     netstr <- c("STRING/","iRef/","coexp/","HPRD/","coexp1/")
#     k = 1
#     for(netflag in c(3,6,7,20,21)){
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
    
    ### nature 13772
    source("CRF_build.R")
    source("HotNet2.R")
    source("Network_analysis.R")
    source("Multi_net.R")
    options(stringsAsFactors=FALSE)
    fileexp=""
    netstr <- c("STRING/","iRef/","coexp/","HPRD/","coexp1/")
    k = 1
    for(netflag in c(3,6,7,20,21)){
        
        if(netflag==3){ instr <- "CRF_input";idm=TRUE; }else{ instr <- "hotnet_input"; idm=FALSE;}
        dirin <- "result/control/"
        dirstr <- paste("result/control/hotnet/",netstr[k],sep="")

        strname <- "nat13772_0"
        filename <- paste(dirin,instr,strname,".txt",sep="")
        HotNet2R(netflag,strname,filename,fileexp,dirstr,cutn=0.50,idm)
        
        strname <- "nat13772"
        filename <- paste(dirin,instr,strname,".txt",sep="")
        HotNet2R(netflag,strname,filename,fileexp,dirstr,cutn=0.50,idm)
        
        k <- k +1
    }
    
#     ### meta analysis 
#     source("CRF_build.R")
#     source("HotNet2.R")
#     source("Network_analysis.R")
#     source("Multi_net.R")
#     options(stringsAsFactors=FALSE)
#     fileexp=""
#     netstr <- c("STRING/","iRef/","coexp/","HPRD/","coexp1/")
#     k = 1
#     for(netflag in c(3,6,7,20,21)){
#         
#         if(netflag==3){ instr <- "CRF_input";idm=TRUE; }else{ instr <- "hotnet_input"; idm=FALSE;}
#         dirin <- "result/control/"
#         dirstr <- paste("result/control/hotnet/",netstr[k],sep="")
#         
#         strname <- "Meta_7_3"
#         filename <- paste(dirin,instr,strname,".txt",sep="")
#         HotNet2R(netflag,strname,filename,fileexp,dirstr,cutn=0.50,idm)
#         
#         k <- k +1
#     }
#     
#     
#     #### all ASD sample
#     source("CRF_build.R")
#     source("HotNet2.R")
#     source("Network_analysis.R")
#     source("Multi_net.R")
#     options(stringsAsFactors=FALSE)
#     fileexp=""
#     netstr <- c("STRING/","iRef/","coexp/","HPRD/","coexp1/")
#     k = 1
#     for(netflag in c(3,6,7,20,21)){
#         
#         if(netflag==3){ instr <- "CRF_input";idm=TRUE; }else{ instr <- "hotnet_input"; idm=FALSE;}
#         
#         dirin <- "result/randset_1/"
#         dirstr <- paste("result/control/hotnet/",netstr[k],sep="")
#         
#         strname <- "ASD4_16"
#         filename <- paste(dirin,instr,strname,".txt",sep="")
#         HotNet2R(netflag,strname,filename,fileexp,dirstr,cutn=0.50,idm)
#         
#         k <- k +1
#     }
    
}

hotnet2_randset_1 <- function(i){
    
    source("CRF_build.R")
    source("HotNet2.R")
    source("Network_analysis.R")
    source("Multi_net.R")
    options(stringsAsFactors=FALSE)
    fileexp=""
    netstr <- c("STRING/","iRef/","coexp/","HPRD/","coexp1/")
    k =4
    for(netflag in c(20,21)){
    #k = 1
    #for(netflag in c(3,6,7,8,21)){
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
    netstr <- c("STRING/","iRef/","coexp/","HPRD/","coexp1/")
    k=4
    for(netflag in c(20,21)){
    #k = 1
    #for(netflag in c(3,6,7,8,21)){
        
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
    netstr <- c("STRING/","iRef/","coexp/","HPRD/","coexp1/")
    k=4
    for(netflag in c(20,21)){
    #k = 1
    #for(netflag in c(3,6,7,20,21)){
        
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
