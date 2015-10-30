DDD <- function(){
    
    ### c("Family_ID","Gene","Type","MutaType","Disease","From")
    ## DDD
    filename <- "DDD_mutations/nature14135-s2/Table s2.csv" ## SNV and indels
    deddd <- read.csv(filename)
    dddlof <- c("frameshift_variant","inframe_deletion","inframe_insertion","splice_acceptor_variant","splice_donor_variant","splice_region_variant","stop_gained","stop_lost")
    dddmis <- c("missense_variant")
    deddd <- deddd[deddd[,"Consequence"] %in% union(dddlof,dddmis),c("Gene","Consequence")]
    colnames(deddd) <- c("Gene","Type")
    
    deddd[deddd[,"Type"] %in% dddlof,"MutaType"] <- "LOF"
    deddd[deddd[,"Type"] %in% dddmis,"MutaType"] <- "MIS"
    deddd[,"Family_ID"] <- 0
    deddd[,"Disease"] <- "DDD"
    deddd[,"From"] <- "nature14135_DDD"
    deddd <- deddd[,c(4,1,2,3,5,6)]
    
    ### ID
    filename <- "DDD_mutations/ID/plosone/journal.pgen.1004772.s002.csv"
    deidplos <- read.csv(filename,skip=1)
    idplof <- c("canonical splice site","Canonical splice site","Consensus splice site","deletion","frameshift_deletion","frameshift_insertion","stopgain_SNV","nonframeshift_deletion","nonframeshift_insertion")  
    idpmis <- c("nonsynonymous_SNV")
    #!!!
    deidplos[78,"Gene.symbol"] <- "PCDHGA1"
    deidplos[is.na(deidplos[,"Polyphen.2"]),"Polyphen.2"] <- 0
    subs1 <- deidplos[,"variant.type"] %in% idplof
    subs2 <- deidplos[,"variant.type"] %in% idpmis & deidplos[,"Polyphen.2"]>=0.957
    deidplos <- deidplos[subs1 | subs2,c("Family.ID","Gene.symbol","variant.type")]
    colnames(deidplos) <- c("Family_ID","Gene","Type")
    deidplos[deidplos[,"Type"] %in% idplof,"MutaType"] <- "LOF"
    deidplos[deidplos[,"Type"] %in% idpmis,"MutaType"] <- "MIS"
    deidplos[,"Disease"] <- "ID"
    deidplos[,"From"] <- "Hamdan et al., PLOS genetics"
    
    filename <- "DDD_mutations/ID/100/all.csv"
    deid100 <- read.csv(filename)
    deid100[is.na(deid100[,8]),8] <- ""
    subs <- !(deid100[,7]=="" & deid100[,8]=="" & deid100[,9]=="" & deid100[,10]=="" & deid100[,11]=="")
    deid100 <- deid100[subs,]
    deid100[deid100[,9]=="D","MutaType"] <- "LOF"
    deid100[deid100[,9]!="D","MutaType"] <- "MIS"
    deid100 <- deid100[,c("V1","V2","V3","MutaType")]
    colnames(deid100) <- c("Family_ID","Gene","Type","MutaType")
    deid100[,"Type"] <- 0
    deid100[,"Disease"] <- "ID"
    deid100[,"From"] <- "de Ligt et al., NEJM"
    
    filename <- "DDD_mutations/ID/nature13394/all.csv"
    deidna <- read.csv(filename,header=FALSE)
    deidna[,2] <- gsub(" $", "", deidna[,2])
    deidna[,6] <- gsub(" $", "", deidna[,6])
    idnalof <- c("Frameshift","Insertion","Nonsense","Splice site")
    idnamis <- "Missense"
    deidna <- deidna[deidna[,6] %in% union(idnalof,idnamis),c("V1","V2","V6")]
    colnames(deidna) <- c("Family_ID","Gene","Type")
    deidna[deidna[,"Type"] %in% idnalof,"MutaType"] <- "LOF"
    deidna[deidna[,"Type"] %in% idnamis,"MutaType"] <- "MIS"
    deidna[,"Disease"] <- "ID"
    deidna[,"From"] <- "nature13394"
    
    deID <- rbind(deidplos,deid100,deidna)
   
    ### EE
    filename <- "DDD_mutations/epileptic/nature12439/all.csv" 
    deeena <- read.csv(filename,header=FALSE)
    deeena[,5] <- gsub(" $", "", deeena[,5])
    deeena[,7] <- gsub(" $", "", deeena[,7])
    eelof <- c("frameshift","inframe deletion","splice acceptor","splice donor","stop gained")
    eemis <- c("missense") # 
    subs1 <- deeena[,7] %in% eemis & grepl("probably_damaging",deeena[,8])
    subs2 <- deeena[,7] %in% eelof
    deeena <- deeena[subs1 | subs2,c("V1","V5","V7")]
    colnames(deeena) <- c("Family_ID","Gene","Type")
    deeena[deeena[,"Type"] %in% eelof,"MutaType"] <- "LOF"
    deeena[deeena[,"Type"] %in% eemis,"MutaType"] <- "MIS"
    deeena[,"Disease"] <- "EE"
    deeena[,"From"] <- "nature12439"
    
    ## ASD
    asddata <- read.csv("ASD/ASDmutationlists.csv")
    colnames(asddata) <- c("Family_ID","Disease","Type","Gene","MutaType","From")
    
    source("Network_analysis.R")
    deddd[,"Gene"] <- mapping_to(deddd[,"Gene"])
    deID[,"Gene"] <- mapping_to(deID[,"Gene"])
    deeena[,"Gene"] <- mapping_to(deeena[,"Gene"])
    write.csv(deddd,file="DDD_mutations/datasheet/DDDmutationlist.csv",row.names=FALSE)
    write.csv(deID,file="DDD_mutations/datasheet/IDmutationlist.csv",row.names=FALSE)
    write.csv(deeena,file="DDD_mutations/datasheet/EEmutationlist.csv",row.names=FALSE)
    
    ALLde <- rbind(deddd,deID,deeena,asddata)
    write.csv(ALLde,file="DDD_mutations/datasheet/ALLdemutationlist.csv",row.names=FALSE)
    
    mutrate <- read.csv("ASD/MutationRatem.csv")
    
    ta <- Random_DDD(deddd,mutrate)
    write.csv(ta,file="DDD_mutations/datasheet/DDDcounttable.csv",row.names=FALSE)
    
    ta <- Random_DDD(deID,mutrate)
    write.csv(ta,file="DDD_mutations/datasheet/IDcounttable.csv",row.names=FALSE)
    
    ta <- Random_DDD(deeena,mutrate)
    write.csv(ta,file="DDD_mutations/datasheet/EEcounttable.csv",row.names=FALSE)
    
    ta <- Random_DDD(ALLde,mutrate)
    write.csv(ta,file="DDD_mutations/datasheet/ALLdecounttable.csv",row.names=FALSE)
    
    ta <- Random_DDD(asddata,mutrate)
    write.csv(ta,file="DDD_mutations/datasheet/ASDcounttable.csv",row.names=FALSE)
}         

DDD_mis <- function(){

    ### c("Family_ID","Gene","Type","MutaType","Disease","From")
    ## DDD
    filename <- "DDD_mutations/nature14135-s2/Table s2.csv" ## SNV and indels
    deddd <- read.csv(filename)
    dddlof <- c("frameshift_variant","inframe_deletion","inframe_insertion","splice_acceptor_variant","splice_donor_variant","splice_region_variant","stop_gained","stop_lost")
    #dddmis <- c("3_prime_UTR_variant","downstream_gene_variant","initiator_codon_variant","intron_variant","missense_variant","regulatory_region_variant","upstream_gene_variant")
    dddmis <- c("missense_variant")
    #deddd <- deddd[deddd[,"Consequence"] %in% union(dddlof,dddmis),c("Chr","Pos","AAchange","Gene","Consequence")]
    deddd <- deddd[,c("Chr","Pos","AAchange","Gene","Consequence")]
    colnames(deddd) <- c("Chr","Pos","AAchange","Gene","Type")
    
    deddd[deddd[,"Type"] %in% dddlof,"MutaType"] <- "LOF"
    deddd[deddd[,"Type"] %in% dddmis,"MutaType"] <- "MIS"
    deddd[,"Family_ID"] <- 0
    deddd[,"Disease"] <- "DDD"
    deddd[,"From"] <- "nature14135_DDD"
    deddd <- deddd[,c(1,2,3,c(4,1,2,3,5,6)+3)]
    
    ### ID
    filename <- "DDD_mutations/ID/plosone/journal.pgen.1004772.s002.csv"
    deidplos <- read.csv(filename,skip=1)
    idplof <- c("canonical splice site","Canonical splice site","Consensus splice site","deletion","frameshift_deletion","frameshift_insertion","stopgain_SNV","nonframeshift_deletion","nonframeshift_insertion")  
    idpmis <- c("nonsynonymous_SNV")
    #!!!
    deidplos[78,"Gene.symbol"] <- "PCDHGA1"
    #deidplos <- deidplos[deidplos[,"variant.type"] %in% union(idplof,idpmis),c("Chr","Position","Detailed.annotation.of.the.variant","Family.ID","Gene.symbol","variant.type")]
    deidplos <- deidplos[,c("Chr","Position","Detailed.annotation.of.the.variant","Family.ID","Gene.symbol","variant.type")]
    colnames(deidplos) <- c("Chr","Pos","AAchange","Family_ID","Gene","Type")
    deidplos[deidplos[,"Type"] %in% idplof,"MutaType"] <- "LOF"
    deidplos[deidplos[,"Type"] %in% idpmis,"MutaType"] <- "MIS"
    deidplos[,"Disease"] <- "ID"
    deidplos[,"From"] <- "Hamdan et al., PLOS genetics"
    
    filename <- "DDD_mutations/ID/100/all.csv"
    deid100 <- read.csv(filename)
    deid100[is.na(deid100[,8]),8] <- ""
    #subs <- !(deid100[,7]=="" & deid100[,8]=="" & deid100[,9]=="" & deid100[,10]=="" & deid100[,11]=="")
    #deid100 <- deid100[subs,]
    deid100[,"Chr"] <- sapply(1:dim(deid100)[1],function(i) unlist(strsplit(deid100[i,"V3"],"\\("))[1])
    deid100[,"Pos"] <- sapply(1:dim(deid100)[1],function(i) unlist(strsplit(deid100[i,"V3"],"g."))[2])
    deid100[,"AAchange"] <- deid100[,"V5"]
    deid100[deid100[,9]=="D","MutaType"] <- "LOF"
    deid100[deid100[,9]!="D","MutaType"] <- "MIS"
    deid100 <- deid100[,c("Chr","Pos","AAchange","V1","V2","V3","MutaType")]
    colnames(deid100) <- c("Chr","Pos","AAchange","Family_ID","Gene","Type","MutaType")
    deid100[,"Type"] <- 0
    deid100[,"Disease"] <- "ID"
    deid100[,"From"] <- "de Ligt et al., NEJM"
    
    filename <- "DDD_mutations/ID/nature13394/all.csv"
    deidna <- read.csv(filename,header=FALSE)
    deidna[,2] <- gsub(" $", "", deidna[,2])
    deidna[,6] <- gsub(" $", "", deidna[,6])
    idnalof <- c("Frameshift","Insertion","Nonsense","Splice site")
    idnamis <- "Missense"
    #deidna <- deidna[deidna[,6] %in% union(idnalof,idnamis),c("V1","V2","V6")]
    
    deidna[,"Chr"] <- sapply(1:dim(deidna)[1],function(i) unlist(strsplit(deidna[i,"V3"],"\\("))[1])
    deidna[,"Pos"] <- sapply(1:dim(deidna)[1],function(i) unlist(strsplit(deidna[i,"V3"],":"))[2])
    deidna[,"AAchange"] <- deidna[,"V4"]
    deidna <- deidna[,c("Chr","Pos","AAchange","V1","V2","V6")]
    colnames(deidna) <- c("Chr","Pos","AAchange","Family_ID","Gene","Type")
    deidna[deidna[,"Type"] %in% idnalof,"MutaType"] <- "LOF"
    deidna[deidna[,"Type"] %in% idnamis,"MutaType"] <- "MIS"
    deidna[,"Disease"] <- "ID"
    deidna[,"From"] <- "nature13394"
    
    deID <- rbind(deidplos,deid100,deidna)
    
    ### EE
    filename <- "DDD_mutations/epileptic/nature12439/all.csv" 
    deeena <- read.csv(filename,header=FALSE)
    deeena[,5] <- gsub(" $", "", deeena[,5])
    deeena[,7] <- gsub(" $", "", deeena[,7])
    eelof <- c("frameshift","inframe deletion","splice acceptor","splice donor","stop gained")
    eemis <- c("missense","3' UTR","5' UTR","downstream","upstream","intronic") # 
    subs1 <- deeena[,7] %in% eemis
    subs2 <- deeena[,7] %in% eelof
    #deeena <- deeena[subs1 | subs2,c("V1","V5","V7")]
    deeena[,"Chr"] <- sapply(1:dim(deeena)[1],function(i) unlist(strsplit(deeena[i,"V2"],":"))[1])
    deeena[,"Pos"] <- sapply(1:dim(deeena)[1],function(i) unlist(strsplit(deeena[i,"V2"],":"))[2])
    deeena[,"AAchange"] <- deeena[,"V3"]
    deeena <- deeena[,c("Chr","Pos","AAchange","V1","V5","V7")]
    colnames(deeena) <- c("Chr","Pos","AAchange","Family_ID","Gene","Type")
    deeena[deeena[,"Type"] %in% eelof,"MutaType"] <- "LOF"
    deeena[deeena[,"Type"] %in% eemis,"MutaType"] <- "MIS"
    deeena[,"Disease"] <- "EE"
    deeena[,"From"] <- "nature12439"
    
    ## ASD
    asddata <- read.csv("ASD/ASDmutationlists_all.csv")
    colnames(asddata) <- c("Chr","Pos","AAchange","Family_ID","Disease","Type","Gene","MutaType","From")
    
    asddata[,"Disease"] <- "ASD"
    
    
    
    source("Network_analysis.R")
    deddd[,"Gene"] <- mapping_to(deddd[,"Gene"])
    deID[,"Gene"] <- mapping_to(deID[,"Gene"])
    deeena[,"Gene"] <- mapping_to(deeena[,"Gene"])
    #write.csv(deddd,file="DDD_mutations/datasheet/DDDmutationlist_Mis.csv",row.names=FALSE)
    #write.csv(deID,file="DDD_mutations/datasheet/IDmutationlist_Mis.csv",row.names=FALSE)
    #write.csv(deeena,file="DDD_mutations/datasheet/EEmutationlist_Mis.csv",row.names=FALSE)
    
    ALLde <- rbind(deddd,deID,deeena,asddata)
    write.csv(ALLde,file="DDD_mutations/datasheet/ALLdemutationlist_all.csv",row.names=FALSE)
    
    #write.csv(ALLde,file="DDD_mutations/datasheet/ALLdemutationlist_Mis.csv",row.names=FALSE)
    
    mutrate <- read.csv("ASD/MutationRatem.csv")
    
    ta <- Random_DDD(deddd,mutrate)
    write.csv(ta,file="DDD_mutations/datasheet/DDDcounttable_Mis.csv",row.names=FALSE)
    
    ta <- Random_DDD(deID,mutrate)
    write.csv(ta,file="DDD_mutations/datasheet/IDcounttable_Mis.csv",row.names=FALSE)
    
    ta <- Random_DDD(deeena,mutrate)
    write.csv(ta,file="DDD_mutations/datasheet/EEcounttable_Mis.csv",row.names=FALSE)
    
    ta <- Random_DDD(ALLde,mutrate)
    write.csv(ta,file="DDD_mutations/datasheet/ALLdecounttable_Mis.csv",row.names=FALSE)
    
    ta <- Random_DDD(asddata,mutrate)
    write.csv(ta,file="DDD_mutations/datasheet/ASDcounttable_Mis.csv",row.names=FALSE)
    
    
    ###========================================================###    
    ### wannovar
    options(stringsAsFactors=FALSE)
    ## DDD
    filename <- "DDD_mutations/nature14135-s2/Table s2.csv" ## SNV and indels
    deddd <- read.csv(filename,stringsAsFactors=FALSE)
    deddda <- deddd[,c("Chr","Pos","Pos")]
    deddda[,"Ref"] <- sapply(1:dim(deddd)[1],function(i) unlist(strsplit(deddd[i,"Ref.Alt"],"/"))[1])
    deddda[,"Alt"] <- sapply(1:dim(deddd)[1],function(i) unlist(strsplit(deddd[i,"Ref.Alt"],"/"))[2])
    deddda[,"Gene"] <- deddd[,"Gene"]
    deddda[,"Category"] <- deddd[,"Consequence"]
    deddda[,"AAchanges"] <- deddd[,"AAchange"]
    colnames(deddda) <- c("Chr","Start","End","Ref","Alt","Gene","Category","AAchanges")
    deddda[,"sampleID"] <- ""
    deddda[,"Disease"] <- "DDD"
    deddda[,"From"] <- "nature14135_DDD"
    
    ### ID
    filename <- "DDD_mutations/ID/plosone/journal.pgen.1004772.s002.csv"
    deidplos <- read.csv(filename,skip=1)
    deidplosa <- deidplos[,c("Chr","Position","Position","Reference.Allele","Mutant.Allele","Gene.symbol","variant.type","Family.ID","Detailed.annotation.of.the.variant")]
    colnames(deidplosa) <- c("Chr","Start","End","Ref","Alt","Gene","Category","sampleID","AAchanges")
    deidplosa[,"Disease"] <- "ID"
    deidplosa[,"From"] <- "Hamdan et al., PLOS genetics"
    
    
    idsampleMap <- read.csv("DDD_mutations/ID/id100-idna13394.csv",skip=2,header=FALSE)
    filename <- "DDD_mutations/ID/100/all.csv"
    deid100 <- read.csv(filename)
    deid100 <- deid100[!(deid100[,"V1"] %in% idsampleMap[,2]),]
    
    deid100[is.na(deid100[,8]),8] <- ""
    tmp <- gsub("Chr","",deid100[,"V3"])
    deid100a <- sapply(1:dim(deid100)[1], function(i) {
        
        sub1 <- grepl("del",tmp[i])
        sub2 <- grepl("-",tmp[i])
        sub3 <- grepl("_",tmp[i])
        a1 <- unlist(strsplit(tmp[i],"\\("))[1]
        subtmp <- unlist(strsplit(tmp[i],"g\\."))[2]
        
        if(sub1 | sub2 | sub3){
            
                subtmp <-  gsub("\\D*$","",subtmp)
                if(sub2){
                    a2 <- as.numeric(unlist(strsplit(subtmp,"-"))[1])
                    a3 <- as.numeric(unlist(strsplit(subtmp,"-"))[2])
                }else if(sub3){
                    a2 <- as.numeric(unlist(strsplit(subtmp,"_"))[1])
                    a3 <- as.numeric(unlist(strsplit(subtmp,"_"))[2])
                }else{
                    a2 <- as.numeric(subtmp)
                    a3 <- a2
                }
            
        }else{
            a2 <- as.numeric(substr(subtmp,1,nchar(subtmp)-3))
            a3 <- a2        
        }
        
        if(grepl(">",subtmp)){
            a4 <- substr(subtmp,nchar(subtmp)-2,nchar(subtmp)-2)
            a5 <- substr(subtmp,nchar(subtmp),nchar(subtmp))
        }else{a4 <- "-";a5 <- "-";}
        
        c(a1,a2,a3,a4,a5)
        
        })
    deid100a <- t(deid100a)
    deid100a <- cbind(deid100a,deid100[,"V2"])
    deid100a <- cbind(deid100a,"synonmous")
    
    subs <- (deid100[,7]=="" & deid100[,8]=="" & deid100[,9]=="" & deid100[,10]=="" & deid100[,11]=="")
    deid100a[deid100[,9]=="D",7] <- "LOF"
    deid100a[deid100[,9]!="D" & !subs ,7] <- "MIS"
    deid100a <- cbind(deid100a,deid100[,"V5"])
    deid100a <- cbind(deid100a,deid100[,"V1"])
    deid100a <- cbind(deid100a,"ID")
    deid100a <- cbind(deid100a,"deLigt_et_al.,_NEJM")
    colnames(deid100a) <- c("Chr","Start","End","Ref","Alt","Gene","Category","AAchanges","sampleID","Disease","From")

    
    filename <- "DDD_mutations/ID/nature13394/all.csv"
    deidna <- read.csv(filename,header=FALSE)
    for(i in 1:16){
        deidna[,i] <- gsub(" $", "", deidna[,i])
    }
    tmp <- gsub("Chr","",deidna[,"V3"])
    deidnaa <- sapply(1:dim(deidna)[1], function(i) {
        
        sub1 <- grepl("del",tmp[i])
        sub2 <- grepl("-",tmp[i])
        sub3 <- grepl("_",tmp[i])
        sub4 <- grepl("dup",tmp[i])
        sub5 <- grepl("ins",tmp[i])
        a1 <- unlist(strsplit(tmp[i],"\\("))[1]
        subtmp <- unlist(strsplit(tmp[i],"g\\."))[2]
        
        if(sub1 | sub2 | sub3 | sub4 | sub5){
#             if(sub1){subtmp <-  unlist(strsplit(subtmp,"del"))[1];}
#             if(sub4){subtmp <-  unlist(strsplit(subtmp,"dup"))[1];}
#             if(sub5){subtmp <-  unlist(strsplit(subtmp,"ins"))[1];}
            subtmp <-  gsub("\\D*$","",subtmp)
            if(sub2){
                a2 <- as.numeric(unlist(strsplit(subtmp,"-"))[1])
                a3 <- as.numeric(unlist(strsplit(subtmp,"-"))[2])
            }else if(sub3){
                a2 <- as.numeric(unlist(strsplit(subtmp,"_"))[1])
                a3 <- as.numeric(unlist(strsplit(subtmp,"_"))[2])
            }else{
                a2 <- as.numeric(subtmp)
                a3 <- a2
            }
        }else{
            a2 <- as.numeric(substr(subtmp,1,nchar(subtmp)-3))
            a3 <- a2        
        }
        
        if(grepl(">",subtmp)){
            a4 <- substr(subtmp,nchar(subtmp)-2,nchar(subtmp)-2)
            a5 <- substr(subtmp,nchar(subtmp),nchar(subtmp))
        }else{a4 <- "-";a5 <- "-";}
        
        c(a1,a2,a3,a4,a5)
        
    })
    deidnaa <- t(deidnaa)
    deidnaa <- cbind(deidnaa,deidna[,"V2"])
    deidnaa <- cbind(deidnaa,deidna[,"V6"])
    deidnaa <- cbind(deidnaa,deidna[,"V1"])
    deidnaa <- cbind(deidnaa,deidna[,"V4"])
    deidnaa <- cbind(deidnaa,"ID")
    deidnaa <- cbind(deidnaa,"nature13394")
    colnames(deidnaa) <- c("Chr","Start","End","Ref","Alt","Gene","Category","sampleID","AAchanges","Disease","From")
    
    ## Lancet
    deidlan <- read.csv("DDD_mutations/ID/Lancet/all.csv",header=FALSE)
    tmp <- gsub("chr","",deidlan[,"V4"])
    tmp <- gsub(" $", "", tmp)
    deidlana <- sapply(1:dim(deidlan)[1], function(i) {
        sub1 <- grepl("del",tmp[i])
        sub2 <- grepl("-",tmp[i])
        sub3 <- grepl("_",tmp[i])
        a1 <- unlist(strsplit(tmp[i],":"))[1]
        subtmp <- unlist(strsplit(tmp[i],"g\\."))[2]
        
        if(sub1 | sub2 | sub3){
                subtmp <-  gsub("\\D*$","",subtmp)
                if(sub2){
                    a2 <- as.numeric(unlist(strsplit(subtmp,"-"))[1])
                    a3 <- as.numeric(unlist(strsplit(subtmp,"-"))[2])
                }else if(sub3){
                    a2 <- as.numeric(unlist(strsplit(subtmp,"_"))[1])
                    a3 <- as.numeric(unlist(strsplit(subtmp,"_"))[2])
                }else{
                    a2 <- as.numeric(subtmp)
                    a3 <- a2
                }
        }else{
            a2 <- as.numeric(substr(subtmp,1,nchar(subtmp)-3))
            a3 <- a2        
        }
        
        if(grepl(">",subtmp)){
            a4 <- substr(subtmp,nchar(subtmp)-2,nchar(subtmp)-2)
            a5 <- substr(subtmp,nchar(subtmp),nchar(subtmp))
        }else{a4 <- "-";a5 <- "-";}
        
        c(a1,a2,a3,a4,a5)
        
    })
    deidlana <- t(deidlana)
    deidlana <- cbind(deidlana,deidlan[,"V2"])
    deidlana <- cbind(deidlana,deidlan[,"V3"])
    deidlana <- cbind(deidlana,deidlan[,"V6"])
    deidlana <- cbind(deidlana,deidlan[,"V1"])
    deidlana <- cbind(deidlana,"ID")
    deidlana <- cbind(deidlana,"Rauch_et_al.,_Lancet")
    colnames(deidlana) <- c("Chr","Start","End","Ref","Alt","Gene","Category","AAchanges","sampleID","Disease","From")
    
#     ## Vissers
#     deidVis <- read.csv("DDD_mutations/ID/Vissers/all.csv",header=FALSE)
#     for(i in 1:17){
#         deidVis[,i] <- gsub(" $", "", deidVis[,i])
#     }
#     deidVisa <- deidVis[,c(1,2,3,4,5,10)]
#     deidVisa[,1] <- gsub("chr","",deidVisa[,1])
#     colnames(deidnaa) <- c("Chr","Start","End","Ref","Alt","Gene","Category","sampleID","AAchanges","Disease","From")

    deIDa <- rbind(deidplosa,deid100a,deidnaa,deidlana)


    
    ### EE
    filename <- "DDD_mutations/epileptic/nature12439/all.csv" 
    deeena <- read.csv(filename,header=FALSE)
    for(i in 1:9){
        deeena[,i] <- gsub(" $", "", deeena[,i])
    }
    deeenaa <- sapply(1:dim(deeena)[1],function(i) unlist(strsplit(deeena[i,"V2"],":"))[1])
    deeenaa <- as.data.frame(deeenaa)
    deeenaa[,"Start"] <- sapply(1:dim(deeena)[1],function(i) unlist(strsplit(deeena[i,"V2"],":"))[2])
    deeenaa[,"End"] <- deeenaa[,"Start"]
    deeenaa[,"Ref"] <- sapply(1:dim(deeena)[1],function(i) unlist(strsplit(deeena[i,"V3"],"/"))[1])
    deeenaa[,"Alt"] <- sapply(1:dim(deeena)[1],function(i) unlist(strsplit(deeena[i,"V3"],"/"))[2])
   # deeenaa[,"Alt"] <- gsub(" $", "", deeenaa[,"Alt"])
    deeenaa[,"Gene"] <- deeena[,"V5"]
    deeenaa[,"AAchanges"] <- ""
    deeenaa[,"Category"] <- deeena[,"V7"]
    deeenaa[,"sampleID"] <- deeena[,"V1"]
    deeenaa[,"Disease"] <- "EE"
    colnames(deeenaa) <- c("Chr","Start","End","Ref","Alt","Gene","AAchanges","Category","sampleID","Disease")
    deeenaa[,"From"] <- "nature12439"
   
    ## ASD
    asddataa <- read.csv("ASD/ASDmutationlists_annovar.csv")
    
    ALLde <- rbind(deIDa,deeenaa,asddataa)
    source("Network_analysis.R")
    ALLde[,"Gene"] <- mapping_to(ALLde[,"Gene"])
    write.csv(ALLde,file="DDD_mutations/datasheet/ALLdemutationlist_annovar.csv",row.names=FALSE)
    write.table(ALLde[,1:5],file="DDD_mutations/datasheet/ALLdemutationlist_annovar.txt",row.names=FALSE,sep="\t",quote=FALSE,col.names=FALSE)


    ALLde <- deddda
    source("Network_analysis.R")
    ALLde[,"Gene"] <- mapping_to(ALLde[,"Gene"])
    write.csv(ALLde,file="DDD_mutations/datasheet/DDDdemutationlist_annovar.csv",row.names=FALSE)
    write.table(ALLde[,1:5],file="DDD_mutations/datasheet/DDDdemutationlist_annovar.txt",row.names=FALSE,sep="\t",quote=FALSE,col.names=FALSE)

}

Random_DDD <- function(data,mutrate,samples=""){
    
    ### data <- data[data[,"Family_ID"] %in% samples,]
    
    dlist <- list()
    dlist[[1]] <- data[data[,"MutaType"]=="LOF",]
    dlist[[2]] <- data[data[,"MutaType"]=="MIS",]
    counttable <- mutrate
    counttable[,"dn.LoF"] <- 0
    counttable[,"dn.mis3"] <- 0
    
    tmp <- table(dlist[[1]][,"Gene"])
    genes <- intersect(names(tmp),counttable[,"Gene"])
    counttable[match(genes,counttable[,"Gene"]),"dn.LoF"] <- tmp[match(genes,names(tmp))]
    
    tmp <- table(dlist[[2]][,"Gene"])
    genes <- intersect(names(tmp),counttable[,"Gene"])
    counttable[match(genes,counttable[,"Gene"]),"dn.mis3"] <- tmp[match(genes,names(tmp))]
    
    counttable
}

TADA_result <- function(){
filenames <- c("DDD_mutations/datasheet/TADAdenovo_DDD.csv","DDD_mutations/datasheet/TADAdenovo_ID.csv","DDD_mutations/datasheet/TADAdenovo_EE.csv","DDD_mutations/datasheet/TADAdenovo_ASD.csv","DDD_mutations/datasheet/TADAdenovo_ALLde.csv")

TPcut=0.1
result <- list()
B <- rep(0,5)
C <- rep(0,5)
for(i in 1:5){
tmp <- read.csv(filenames[i])
result[[i]]  <- tmp[tmp[,"qvalue.dn"] <= TPcut,1]
B[i] <- length(result[[i]])
}

for(i in 1:5){
    C[i] <- length(intersect(result[[i]],result[[5]]))
}

pdf(file="DDD_mutations/bar.pdf")
legend=c("DDD","ID","EE","ASD","ALL")  
cols <- c("Black","Blue","Cyan","Brown","Coral","DarkViolet","Green","HotPink","Orange","Plum")
par(mar=c(4,4,1,1))
barplot(B,col=cols[1:5],names.arg=legend,ylim=c(0,max(B+15)),ylab="The number of genes with FDR < 0.1")
text(seq(0.8,12,1.2)-0.1,B,paste(B,"(",C,")",sep=""),pos=3,cex=2)
#text(rep(0,3),c(31,33,35),c("*    p < 0.01","**   p < 1e-30","***  p < 1e-40"),pos=4,cex=1.2)
dev.off()

a <- union(result[[1]],union(result[[2]],union(result[[3]],result[[4]])))
genes <- setdiff(result[[5]],a)
tmp <- read.csv(filenames[5])
subtab <- cbind(match(genes,tmp[,1]),tmp[match(genes,tmp[,1]),])
write.csv(subtab,"DDD_mutations/bargenes.csv",row.names=FALSE)


### TADA input
source("Network_analysis.R")
source("Multi_net.R")
dirstr <- "DDD_mutations/netinput/"
for(i in 1:5){
filename=filenames[i]
strname <- c("DDD","ID","EE","ASD","ALL")
TADAinput1(filename,strname[i],dirstr,genelist0=1)
}

}

exp_data <- function(){

### see "data/brain_expression/pre.R"
}

oddsratiof <- function(){
    source("Network_analysis.R")
    mutrate <- read.csv("ASD/GeneMutaRate.csv")
    mutrate[,1] <- mapping_to(mutrate[,1])
    
    library(fmsb)
    countT <- c("DDD_mutations/datasheet/DDDcounttable.csv","DDD_mutations/datasheet/IDcounttable.csv","DDD_mutations/datasheet/EEcounttable.csv","DDD_mutations/datasheet/ASDcounttable.csv","DDD_mutations/datasheet/ALLdecounttable.csv")
    strname <- c("DDD_mutations/datasheet/DDDoddstable.csv","DDD_mutations/datasheet/IDoddstable.csv","DDD_mutations/datasheet/EEoddstable.csv","DDD_mutations/datasheet/ASDoddstable.csv","DDD_mutations/datasheet/ALLdeoddstable.csv")
    
    countT <- c("DDD_mutations/datasheet/DDDcounttable_Mis.csv","DDD_mutations/datasheet/IDcounttable_Mis.csv","DDD_mutations/datasheet/EEcounttable_Mis.csv","DDD_mutations/datasheet/ASDcounttable_Mis.csv","DDD_mutations/datasheet/ALLdecounttable_Mis.csv")
    strname <- c("DDD_mutations/datasheet/DDDoddstable_Mis.csv","DDD_mutations/datasheet/IDoddstable_Mis.csv","DDD_mutations/datasheet/EEoddstable_Mis.csv","DDD_mutations/datasheet/ASDoddstable_Mis.csv","DDD_mutations/datasheet/ALLdeoddstable_Mis.csv")
    
    Nc <- c(1133,191,264,3953,5541)
    
    for(i in 1:length(countT)){
        count <- read.csv(countT[i])
        odds <-  sapply(1:dim(count)[1], function(k) {
            subs <- match(count[k,1],mutrate[,1])
            a <- count[k,3]+count[k,4]
            if(a==0){
                c(0,0,0)
                }else{
                b <- Nc[i] - a
                c <- 2 ## two allele
                d <- 1/mutrate[subs,5] - 2
                tmp <- oddsratio(a,b,c,d)
                c(tmp$estimate,tmp$conf.int[1],tmp$conf.int[2])
            }
            })
        odds <- t(odds)
        odds[is.na(odds)] <- 0
        count <- cbind(count,odds)
        colnames(count) <- c("Gene","mut.rate","dn.LoF","dn.mis3","oddsratio","Lower","Upper")
        write.csv(count,file=strname[i],row.names=FALSE)
    }
    
    
}

cytoscape_ana <- function(){

filename <- "DDD_mutations/datasheet/TADAdenovo_ALLde.csv"
a <- read.csv(filename)
b <- cbind(a[,1],a[,3]+a[,4])
colnames(b) <- c("Gene","SampleNumber")
b <- b[a[,"qvalue.dn"] <= 0.1,]
write.table(b,file="functions/DDD/genepair.txt",quote=FALSE,row.names=FALSE,sep="\t")
write.table(b[,1],file="functions/DDD/genelist.txt",quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
}

meta_analysis <- function(){
    
    source("DDD.R")
    source("Network_analysis.R")
    source("enrichana.R")
    
    TADAfile <- "DDD_mutations/datasheet/TADAdenovo_DDD.csv"
    MISfile <- "result/DDD/CRFresultDDD.txt"
    netflag <- 14
    idmp <- FALSE
    N=1133
    a <- infor_table(TADAfile,MISfile,netflag,N,idmap=FALSE)    
    
    tmpgenes <- c("POGZ","CHAMP1","BCL11A","NAA15","CUL3","TRIO")
    a[a[,1] %in% tmpgenes,]
    
#     TADAfile <- "DDD_mutations/datasheet/TADAdenovo_ALLde.csv"
#     MISmetafile <- "result/DDD/CRFresultALLde.txt"
#     b <- infor_table(TADAfile,MISmetafile,netflag,N=5541,idmap=FALSE)
#     
    genea <- a[a[,"MIS_FDR"]<=0.1,1]
    #geneb <- b[b[,"MIS_score"]>=0.6,1]
    
    #TADAg <-  read.csv("DDD_mutations/datasheet/TADAdenovo_DDD.csv")
    #TADAg <- TADAg[TADAg[,"qvalue.dn"]<=0.1,1]
    
#     tmp <- test_genesDDD(genea)
    #tmp <- test_genesDDD(TADAg)
    
#     subs <- c(1:3,5)
#     barplot(t(cbind(tmp[subs,1],tmp[subs,2]-tmp[subs,1])),col=c("red","blue"),names.arg=c("DDD knwon","EE/ID/Autism known","TADA FDR<=0.1","Total"),cex.names=2,cex.axis=2)
#     legend("topleft",col=c("red","blue"),fill=c("red","blue"),legend=c("Confirmed","Not-confirmed"))
    #test_genes(geneb)
    
    DDDgenes <- unlist(read.table("DDD_mutations/nature14135-s2/Genes.txt"))
    knowng <- DDDgenes[1:34]
    knowng <- mapping_to(knowng)
    
    pdf(file="plot/DDD1.pdf",width=12,height=12)
    pg <- intersect(knowng,a[,1])
    x <- -log10(apply(a[match(pg,a[,1]),c("poission_test_LOFMIS_pvalue","poission_test_LOF_pvalue","poission_test_MIS_pvalue")],1,min))
    #y <- -log10(3*apply(b[match(pg,b[,1]),c("poission_test_LOFMIS_pvalue","poission_test_LOF_pvalue","poission_test_MIS_pvalue")],1,min))
    y <- -log10(a[match(pg,a[,1]),c("MIS_FDR")])
    cuts <- -log10(0.05/(3*dim(a)[1])) 
    pg1 <- intersect(pg,genea)

    cols=rep(1,length(pg))
    cols[match(pg1,pg)] <- 2
    plot(x,y,col=cols,xlab="-log10(p) DDD-only",ylab="-log10(FDR)")
    abline(v=cuts,h=-log(0.1),lty=2)
    text(x,y,labels=pg,pos=3,col=cols,adj=c(1,1),cex=1)
    dev.off()
    
source("DDD.R")
source("Network_analysis.R")
source("enrichana.R")

TADAfile <- "DDD_mutations/datasheet/TADAdenovo_ALLde.csv"
MISfile <- "result/DDD/CRFresultALLde.txt"
netflag <- 14
idmp <- FALSE
N=5541
a <- infor_table(TADAfile,MISfile,netflag,N,idmap=FALSE)    
genea <- a[a[,"MIS_FDR"]<=0.1,1]

DDDgenes <- unlist(read.table("DDD_mutations/nature14135-s2/Genes.txt"))
knowng <- DDDgenes[1:34]
knowng <- mapping_to(knowng)

pdf(file="plot/DDD2.pdf",width=12,height=12)
pg <- intersect(knowng,a[,1])
x <- -log10(apply(a[match(pg,a[,1]),c("poission_test_LOFMIS_pvalue","poission_test_LOF_pvalue","poission_test_MIS_pvalue")],1,min))
y <- -log10(a[match(pg,a[,1]),c("MIS_FDR")])
cuts <- -log10(0.05/(3*dim(a)[1])) 
pg1 <- intersect(pg,genea)

cols=rep(1,length(pg))
cols[match(pg1,pg)] <- 2
plot(x,y,col=cols,xlab="-log10(p) DDD-only",ylab="-log10(FDR)")
abline(v=cuts,h=-log(0.1),lty=2)
text(x,y,labels=pg,pos=3,col=cols,adj=c(1,1),cex=1)
dev.off()


    subs <- which(y <= cuts)
    plot(x[subs],y[subs],col=cols[subs],ylim=c(0,max(y[subs])+3))
    abline(v=cuts,h=cuts,lty=2)
    text(x[subs],y[subs],labels=pg[subs],pos=3,col=cols[subs],adj=c(1,1),cex=1)
    
#     filename <- "DDD_mutations/nature14135-s2/Genes.txt"
#     genes <- unlist(read.table(filename))
#     
#     genes <- genes[1:34]
#     genes <- mapping_to(genes)
#     
#     genes1 <- c("ALG13","SCN1A","GABRB3","MECP2","POGZ","CHAMP1","FOXP1","DNM1","TRIO","CUL3","NAA15")
#     genes2 <- c("BCL11A","KMT2A","PURA","SMARCA2")
#     genes3 <- setdiff(genes,union(genes1,genes2))

}

infor_table <- function(TADAfile,MISfile,netflag,N,idmap=FALSE,miscut=0.1,tadacut=0.3){
    options(stringsAsFactors=FALSE)
    source("Multi_net.R")
    source("CRF_build.R")
    source("Network_analysis.R")
    source("enrichana.R")
    
    coexpnet <- build_net(netflag,"")
    if(idmp){
        mapT <- as.matrix(read.delim("Fmap0121.txt",header=FALSE,sep="\t"))
        coexpnet$node <- mapT[match(coexpnet$node,mapT[,1]),2]
        dimnames(coexpnet$matrix) <- list(coexpnet$node,coexpnet$node)
    }
    
    CRFr <- read.table(MISfile)
    TADAr <- read.csv(TADAfile)
    
    genes <- intersect(CRFr[,1],TADAr[,1])
    CRFr <- CRFr[CRFr[,1] %in% genes,]
    TADAr <- TADAr[TADAr[,1] %in% genes,]
    
    genes1 <- CRFr[CRFr[,5]<=miscut,1]
    genes2 <- TADAr[TADAr[,"qvalue.dn"]<=tadacut,1]
    
    subs <- match(CRFr[,1],TADAr[,1])
    resultT <- cbind(CRFr[,1],TADAr[subs,c("dn.LoF","dn.mis3")])
    #resultT[is.na(resultT)] <- 0
    
    source("enrichana.R")
    poipvalue <- poisson_testhq(resultT,N)
    
    resultT <- cbind(resultT,poipvalue,TADAr[subs,"qvalue.dn"])
    resultT <- cbind(resultT,CRFr[,c(2,5)])
    
    adj <- coexpnet$matrix[genes,genes]
    adj[adj>0] <- 1
    resultT <- cbind(resultT,rowSums(adj))
    
    resultT <- cbind(resultT,rowSums(adj[,genes1]))
    resultT <- cbind(resultT,rowSums(adj[,genes2]))
    
    colnames(resultT) <- c("Gene","dn.LoF","dn.mis3","poission_test_LOF_pvalue","poission_test_MIS_pvalue","poission_test_LOFMIS_pvalue","TADA_FDR","MIS_score","MIS_FDR","Degree_in_network","edges_in_FDR<=0.1_MISresult","edges_in_FDR<=0.3_TADAresult")
    
    #write.csv(resultT,file="result/PCGC/resultTable.csv",row.names=FALSE)
    resultT
}

test_genesDDD <- function(genea){

    ### known DDD genes
    DDDgenes <- unlist(read.table("DDD_mutations/nature14135-s2/Genes.txt"))
    knowng <- DDDgenes[1:34]
    zibrag <- DDDgenes[35:55]
    
    #DDDtmp <- unlist(read.table("ASD/Targets/DDDgenes.txt"))
    ### Other ID/Autism/EE genes (candidate genes)
    ID1 <- unlist(read.table("DDD_mutations/ID/100/Genes.txt"))
    ID2 <- unlist(read.table("DDD_mutations/ID/nature13394/Genes.txt"))
    ID3 <- unlist(read.table("DDD_mutations/ID/plosone/Genes.txt"))
    knownid <- union(ID1[1:14],union(ID2,ID3[1:12]))
    candi <- union(ID1[15:35],ID3[13:24])
    candi <- setdiff(candi,knownid)
    EEg <- unlist(read.table("DDD_mutations/epileptic/nature12439/Genes.txt"))
    ASDg <- read.csv("ASD/Targets/SFARIgene/gene-S-H-Strong.csv")[,1]
    
    #othdg <- union(union(knownid,candi),union(EEg,ASDg))
    othdg <- union(knownid,union(EEg,ASDg))
    ### confirmed with TADA 
    TADAg <-  read.csv("DDD_mutations/datasheet/TADAdenovo_DDD.csv")
    TADAg <- TADAg[TADAg[,"qvalue.dn"]<=0.1,1]
        
    ### Metag 
    metag <- read.csv("DDD_mutations/datasheet/TADAdenovo_ALLde.csv")
    metag <- metag[metag[,"qvalue.dn"]<=0.1,1]
    
    ### Others
    nog <- setdiff(genea,c(DDDgenes,othdg,TADAg))
    
    a1 <- length(intersect(genea,knowng))
    a2 <- length(intersect(genea,othdg))
    a3 <- length(intersect(genea,TADAg))
    a4 <- length(intersect(genea,metag))
    a5 <- length(setdiff(genea,nog))
    
    b1 <- length(unique(knowng))
    b2 <- length(unique(othdg))
    b3 <- length(unique(TADAg))
    b4 <- length(unique(metag))
    b5 <- length(unique(genea))
    
    cbind(c(a1,a2,a3,a4,a5),c(b1,b2,b3,b4,b5))
}

Clear_mis_annovar <- function(){

    reann <- read.csv("DDD_mutations/datasheet/wannovar.genome_summary.csv")
    reall <- read.csv("DDD_mutations/datasheet/ALLdemutationlist_annovar.csv")
    snps <- paste(reall[,1],reall[,2],reall[,3],reall[,4],reall[,5],sep="_")
    snpsann <- paste(reann[,1],reann[,2],reann[,3],reann[,4],reann[,5],sep="_")
 
    ## duplicatd snps
#     subs1 <- unique(snps[duplicated(snps)])
#     subs2 <- unique(snpsann[duplicated(snpsann)])
#     setdiff(subs1,subs2)
    
    subs1 <- snps %in% snpsann
    subs2 <- !subs1
    reall1 <- reall[subs2,]
    #table(reall1[,"Category"])

    reall1 <- cbind(reall1,"","","","","")
   #c("3' UTR","3'UTR","5'UTR","CODON_CHANGE_PLUS_CODON_DELETION","CODON_DELETION","deletion","downstream","frame-shift","frameshif","frameshift","Frameshift","frameshift_deletion","inframe deletion","Insertion","intergenic","intron","intronic","LOF","no-frame-shift","noEnd","non-coding","nonframeshift_deletion","noStart","splice","splice acceptor","splice-site")
    sublof <- c("CODON_CHANGE_PLUS_CODON_DELETION","CODON_DELETION","deletion","frame-shift","frameshif","frameshift","Frameshift","frameshift_deletion","inframe deletion","Insertion","LOF","splice","splice acceptor","splice-site")
    submis <- c("no-frame-shift","nonframeshift_deletion")
    suboth <- c("3' UTR","3'UTR","5'UTR","downstream","intergenic","intron","intronic","noEnd","non-coding","noStart")
    reall1[reall1[,"Category"] %in% sublof,"mutation_type"] <- "LOF"
    reall1[reall1[,"Category"] %in% submis,"mutation_type"] <- "MIS"
    reall1[reall1[,"Category"] %in% suboth,"mutation_type"] <- "Other"
    colnames(reall1) <- c("Chr","Start","End","Ref","Alt","Gene","Category","sampleID","AAchanges","Disease","From","Func.refgene","Gene.refgene","ExonicFunc.refgene","Polyphen2_HDIV_pred","AAChange.refgene","mutation_type")
        
    reall2 <- reall[subs1,]
    #subs <- match(snpsann,snps)
    reall2[,"Func.refgene"] <- reann[,"Func.refgene"]
    reall2[,"Gene.refgene"] <- reann[,"Gene.refgene"]
    reall2[,"ExonicFunc.refgene"] <- reann[,"ExonicFunc.refgene"]
    reall2[,"Polyphen2_HDIV_pred"] <- reann[,"Polyphen2_HDIV_pred"]
    
    reall2[,"AAChange.refgene"] <- reann[,"AAChange.refgene"]
    subs <- reann[,"AAChange.refgene"]==""
    reall2[subs,"AAChange.refgene"] <- reann[subs,"GeneDetail.refgene"]
    
    #subs1 <- reann[,"ExonicFunc.refgene"]=="unknown"
    #noann <- c("downstream","intergenic","intronic","ncRNA_exonic","ncRNA_exonic;ncRNA_exonic;ncRNA_exonic","ncRNA_exonic;ncRNA_exonic;ncRNA_exonic;ncRNA_exonic;ncRNA_exonic","ncRNA_exonic;splicing","ncRNA_intronic","splicing","splicing;splicing","upstream","UTR3","UTR3;UTR3;UTR3","UTR5")
    
    annlof <- c("frameshift substitution","stopgain","stopgain;stopgain","stoploss")
    subslof <- reann[,"ExonicFunc.refgene"] %in% annlof | (reann[,"ExonicFunc.refgene"]=="" & (reann[,"Func.refgene"] %in% c("splicing","splicing;splicing")))
    subsmisd <- (grepl("nonsynonymous SNV",reann[,"ExonicFunc.refgene"]) | grepl("nonframeshift substitution",reann[,"ExonicFunc.refgene"])) & reann[,"Polyphen2_HDIV_pred"]=="D"
    subsmis <- (grepl("nonsynonymous SNV",reann[,"ExonicFunc.refgene"]) | grepl("nonframeshift substitution",reann[,"ExonicFunc.refgene"])) & reann[,"Polyphen2_HDIV_pred"]!="D"
    subssyn <- grepl("synonymous SNV",reann[,"ExonicFunc.refgene"]) & !grepl("nonsynonymous SNV",reann[,"ExonicFunc.refgene"])
    subsintron <- reann[,"ExonicFunc.refgene"]=="" & !(reann[,"Func.refgene"] %in% c("splicing","splicing;splicing"))
    subsunkown <- reann[,"ExonicFunc.refgene"]=="unknown"
    
    reall2[subslof,"mutation_type"] <- "LOF"
    reall2[subsmisd,"mutation_type"] <- "dMIS"
    reall2[subsmis,"mutation_type"] <- "MIS"
    reall2[subssyn,"mutation_type"] <- "SYNONMOUS"
    reall2[subsintron,"mutation_type"] <- "Other"
    
    subsunkown <- reann[,"ExonicFunc.refgene"]=="unknown"
    subsunsyn <- reann[,"ExonicFunc.refgene"]=="unknown" & reall2[,"Category"] %in% c("synonmous","synonymous")
    subsunmis <- reann[,"ExonicFunc.refgene"]=="unknown" & reall2[,"Category"]== "missense" & reann[,"Polyphen2_HDIV_pred"]=="D"
    subsundmis <- reann[,"ExonicFunc.refgene"]=="unknown" & reall2[,"Category"]== "missense" & reann[,"Polyphen2_HDIV_pred"]!="D"
    reall2[subsunsyn,"mutation_type"] <- "SYNONMOUS"
    reall2[subsunmis,"mutation_type"] <- "MIS"
    reall2[subsundmis,"mutation_type"] <- "dMIS"
    colnames(reall2) <- c("Chr","Start","End","Ref","Alt","Gene","Category","sampleID","AAchanges","Disease","From","Func.refgene","Gene.refgene","ExonicFunc.refgene","Polyphen2_HDIV_pred","AAChange.refgene","mutation_type")

    reallnew <- rbind(reall1,reall2)
    reallnew <- reallnew[order(reallnew[,"Disease"]),]
table(reallnew[,"mutation_type"])
    write.csv(reallnew,file="DDD_mutations/datasheet/mutationlist_4_15_all.csv",row.names=FALSE)

    reallnew1 <- reallnew[,c("Chr","Start","Ref","Alt","sampleID","Disease","Gene","AAChange.refgene","mutation_type")]
    colnames(reallnew1) <- c("Chr","Pos","Ref","Alt","sampleID","Disease","Gene","AAChange","mutation_type")
    write.csv(reallnew1[,c("Gene","Chr","Pos","Ref","Alt","mutation_type","AAChange","sampleID","Disease")],file="DDD_mutations/datasheet/mutationlist_4_15.csv",row.names=FALSE)

    ### Table 1 (gene-centric):
    ## gene, number_of_LGD, number_of_D-mis, number_of_B-mis, number_of_silent (if available). The last four columns repeat for each disease. 
    
    ### Please list the exact number of probands profiled in each disease. 
    
    
    ### Table 2 (mutation centric):
    ###    gene, chr, pos, ref, alt, mutation_type,  aa-changes, meta-svm_prediction, sampleID, disease. 
    
}

Clear_misDDD_annovar <- function(){
    
    reann <- read.csv("DDD_mutations/datasheet/DDD.genome_summary.csv")
    reall <- read.csv("DDD_mutations/datasheet/DDDdemutationlist_annovar.csv")
    snps <- paste(reall[,1],reall[,2],reall[,3],reall[,4],reall[,5],sep="_")
    snpsann <- paste(reann[,1],reann[,2],reann[,3],reann[,4],reann[,5],sep="_")
    
    ## duplicatd snps
    #     subs1 <- unique(snps[duplicated(snps)])
    #     subs2 <- unique(snpsann[duplicated(snpsann)])
    #     setdiff(subs1,subs2)
    
    subs1 <- snps %in% snpsann
    subs2 <- !subs1
    reall1 <- reall[subs2,]
    #table(reall1[,"Category"])
    
    reall1 <- cbind(reall1,"","","","","")
    #c("frameshift_variant","inframe_deletion","inframe_insertion","intron_variant","splice_donor_variant","splice_region_variant","upstream_gene_variant")
    
    sublof <- c("frameshift_variant","inframe_deletion","inframe_insertion","splice_donor_variant","splice_region_variant")
    suboth <- c("intron_variant","upstream_gene_variant")
    reall1[reall1[,"Category"] %in% sublof,"mutation_type"] <- "LOF"
    reall1[reall1[,"Category"] %in% suboth,"mutation_type"] <- "Other"
    colnames(reall1) <- c("Chr","Start","End","Ref","Alt","Gene","Category","sampleID","AAchanges","Disease","From","Func.refgene","Gene.refgene","ExonicFunc.refgene","Polyphen2_HDIV_pred","AAChange.refgene","mutation_type")
    
    reall2 <- reall[subs1,]
    #subs <- match(snpsann,snps)
    reall2[,"Func.refgene"] <- reann[,"Func.refgene"]
    reall2[,"Gene.refgene"] <- reann[,"Gene.refgene"]
    reall2[,"ExonicFunc.refgene"] <- reann[,"ExonicFunc.refgene"]
    reall2[,"Polyphen2_HDIV_pred"] <- reann[,"Polyphen2_HDIV_pred"]
    
    reall2[,"AAChange.refgene"] <- reann[,"AAChange.refgene"]
    subs <- reann[,"AAChange.refgene"]==""
    reall2[subs,"AAChange.refgene"] <- reann[subs,"GeneDetail.refgene"]
    
    #table(reann[,"ExonicFunc.refgene"])
    #noann <- c("frameshift substitution","nonsynonymous SNV","nonsynonymous SNV;nonsynonymous SNV","nonsynonymous SNV;nonsynonymous SNV;nonsynonymous SNV","nonsynonymous SNV;nonsynonymous SNV;nonsynonymous SNV;nonsynonymous SNV","stopgain","stoploss","synonymous SNV","synonymous SNV;synonymous SNV","unknown")
    
    annlof <- c("frameshift substitution","stopgain","stoploss")
    subslof <- reann[,"ExonicFunc.refgene"] %in% annlof | (reann[,"ExonicFunc.refgene"]=="" & (reann[,"Func.refgene"] %in% c("splicing","splicing;splicing")))
    subsmisd <- (grepl("nonsynonymous SNV",reann[,"ExonicFunc.refgene"]) | grepl("nonframeshift substitution",reann[,"ExonicFunc.refgene"])) & reann[,"Polyphen2_HDIV_pred"]=="D"
    subsmis <- (grepl("nonsynonymous SNV",reann[,"ExonicFunc.refgene"]) | grepl("nonframeshift substitution",reann[,"ExonicFunc.refgene"])) & reann[,"Polyphen2_HDIV_pred"]!="D"
    subssyn <- grepl("synonymous SNV",reann[,"ExonicFunc.refgene"]) & !grepl("nonsynonymous SNV",reann[,"ExonicFunc.refgene"])
    subsintron <- reann[,"ExonicFunc.refgene"]=="" & !(reann[,"Func.refgene"] %in% c("splicing","splicing;splicing"))
    subsunkown <- reann[,"ExonicFunc.refgene"]=="unknown"
    
    reall2[subslof,"mutation_type"] <- "LOF"
    reall2[subsmisd,"mutation_type"] <- "dMIS"
    reall2[subsmis,"mutation_type"] <- "MIS"
    reall2[subssyn,"mutation_type"] <- "SYNONMOUS"
    reall2[subsintron,"mutation_type"] <- "Other"
    
    subsunkown <- reann[,"ExonicFunc.refgene"]=="unknown"
    subsunsyn <- reann[,"ExonicFunc.refgene"]=="unknown" & reall2[,"Category"]=="synonymous_variant"
    subsunmis <- reann[,"ExonicFunc.refgene"]=="unknown" & reall2[,"Category"]== "missense_variant" & reann[,"Polyphen2_HDIV_pred"]=="D"
    subsundmis <- reann[,"ExonicFunc.refgene"]=="unknown" & reall2[,"Category"]== "missense_variant" & reann[,"Polyphen2_HDIV_pred"]!="D"
    reall2[subsunsyn,"mutation_type"] <- "SYNONMOUS"
    reall2[subsunmis,"mutation_type"] <- "MIS"
    reall2[subsundmis,"mutation_type"] <- "dMIS"
    colnames(reall2) <- c("Chr","Start","End","Ref","Alt","Gene","Category","sampleID","AAchanges","Disease","From","Func.refgene","Gene.refgene","ExonicFunc.refgene","Polyphen2_HDIV_pred","AAChange.refgene","mutation_type")
    
    reallnew <- rbind(reall1,reall2)
    reallnew <- reallnew[order(reallnew[,"Disease"]),]
    table(reallnew[,"mutation_type"])
    write.csv(reallnew,file="DDD_mutations/datasheet/DDDmutationlist_4_15_all.csv",row.names=FALSE)
    reallnew1 <- reallnew[,c("Chr","Start","Ref","Alt","sampleID","Disease","Gene","AAChange.refgene","mutation_type")]
    colnames(reallnew1) <- c("Chr","Pos","Ref","Alt","sampleID","Disease","Gene","AAChange","mutation_type")
    write.csv(reallnew1[,c("Gene","Chr","Pos","Ref","Alt","mutation_type","AAChange","sampleID","Disease")],file="DDD_mutations/datasheet/DDDmutationlist_4_15.csv",row.names=FALSE)
    
    ### Table 1 (gene-centric):
    ## gene, number_of_LGD, number_of_D-mis, number_of_B-mis, number_of_silent (if available). The last four columns repeat for each disease. 
    
    ### Please list the exact number of probands profiled in each disease. 
    
    
    ### Table 2 (mutation centric):
    ###    gene, chr, pos, ref, alt, mutation_type,  aa-changes, meta-svm_prediction, sampleID, disease. 
    
}

Clear_mis_metaSVM <- function(){
    
    reann <- read.csv("DDD_mutations/datasheet/Annotated_mutationlist_4_15_metaSVM.csv")
    reall <- read.csv("DDD_mutations/datasheet/mutationlist_4_15_all.csv") 
    reall <- cbind(reall,reann[,c("Func.refgene","Gene.refgene","ExonicFunc.refgene","RadialSVM_pred","AAChange.refgene")])
    reall <- cbind(reall,"")
    colnames(reall) <- c("Chr","Start","End","Ref","Alt","Gene","Category","sampleID","AAchanges","Disease","From","Func.annovar","Gene.annovar","Category_annovar","Polyphen2_HDIV_pred_annovar","AAChange.annovar","mutation_type_annovar","Func.metaSVM","Gene.metaSVM","Category_metaSVM","RadialSVM_pred_metaSVM","AAChange.metaSVM","mutation_type_metaSVM")
    
    ## table(reann[,"ExonicFunc.refgene"])
    # c("","frameshift deletion","frameshift insertion","nonframeshift deletion","nonframeshift insertion","nonsynonymous SNV","nonsynonymous SNV;nonsynonymous SNV","nonsynonymous SNV;nonsynonymous SNV;nonsynonymous SNV","stopgain","stopgain;stopgain","stoploss","synonymous SNV","synonymous SNV;synonymous SNV","synonymous SNV;synonymous SNV;synonymous SNV","unknown")

    #subs <- reann[,"ExonicFunc.refgene"]==""
    #table(reann[subs,"Func.refgene"])
    ##c("","downstream","intergenic","intronic","ncRNA_exonic","ncRNA_exonic;ncRNA_exonic;ncRNA_exonic","ncRNA_exonic;ncRNA_exonic;ncRNA_exonic;ncRNA_exonic;ncRNA_exonic","ncRNA_exonic;splicing","ncRNA_intronic","splicing","splicing;splicing","upstream","UTR3","UTR3;UTR3;UTR3","UTR5") 
    annlof <- c("frameshift deletion","frameshift insertion","stopgain","stopgain;stopgain","stoploss")
    subslof <- reann[,"ExonicFunc.refgene"] %in% annlof | (reann[,"ExonicFunc.refgene"]=="" & (reann[,"Func.refgene"] %in% c("splicing","splicing;splicing","")))
    #subslof1  <- reann[,"ExonicFunc.refgene"]=="" & reann[,"Func.refgene"]=="" ## table(reall[subslof1,"Category"])
    subsmisd <- (grepl("nonsynonymous SNV",reann[,"ExonicFunc.refgene"]) | grepl("nonframeshift deletion",reann[,"ExonicFunc.refgene"]) | grepl("nonframeshift insertion",reann[,"ExonicFunc.refgene"])) & reann[,"RadialSVM_pred"]=="D"
    subsmis <- (grepl("nonsynonymous SNV",reann[,"ExonicFunc.refgene"]) | grepl("nonframeshift deletion",reann[,"ExonicFunc.refgene"]) | grepl("nonframeshift insertion",reann[,"ExonicFunc.refgene"])) & reann[,"RadialSVM_pred"]!="D"
    subssyn <- grepl("synonymous SNV",reann[,"ExonicFunc.refgene"]) & !grepl("nonsynonymous SNV",reann[,"ExonicFunc.refgene"])
    subsintron <- reann[,"ExonicFunc.refgene"]=="" & !(reann[,"Func.refgene"] %in% c("splicing","splicing;splicing",""))
    subsunkown <- reann[,"ExonicFunc.refgene"]=="unknown"
    
    reall[subslof,"mutation_type_metaSVM"] <- "LOF"
    reall[subsmisd,"mutation_type_metaSVM"] <- "dMIS"
    reall[subsmis,"mutation_type_metaSVM"] <- "MIS"
    reall[subssyn,"mutation_type_metaSVM"] <- "SYNONMOUS"
    reall[subsintron,"mutation_type_metaSVM"] <- "Other"
    
    subsunkown <- reann[,"ExonicFunc.refgene"]=="unknown"
    subsunlof <- reann[,"ExonicFunc.refgene"]=="unknown" & reall[,"Category"] %in% c("frame-shift")
    subsunsyn <- reann[,"ExonicFunc.refgene"]=="unknown" & reall[,"Category"] %in% c("synonmous","synonymous")
    subsundmis <- reann[,"ExonicFunc.refgene"]=="unknown" & reall[,"Category"]== "missense" & reann[,"RadialSVM_pred"]=="D"
    subsunmis <- reann[,"ExonicFunc.refgene"]=="unknown" & reall[,"Category"]== "missense" & reann[,"RadialSVM_pred"]!="D"
    
    reall[subsunlof,"mutation_type_metaSVM"] <- "LOF"
    reall[subsunsyn,"mutation_type_metaSVM"] <- "SYNONMOUS"
    reall[subsunmis,"mutation_type_metaSVM"] <- "MIS"
    reall[subsundmis,"mutation_type_metaSVM"] <- "dMIS"
    reall <- reall[,c("Chr","Start","End","Ref","Alt","Gene","Category","sampleID","AAchanges","Disease","From","Func.annovar","AAChange.annovar","Gene.annovar","Category_annovar","Polyphen2_HDIV_pred_annovar","mutation_type_annovar","Func.metaSVM","AAChange.metaSVM","Gene.metaSVM","Category_metaSVM","RadialSVM_pred_metaSVM","mutation_type_metaSVM")]
    
    write.csv(reall,file="DDD_mutations/datasheet/anno_mutationlist_4_15_all.csv",row.names=FALSE)
    
    reallnew1 <- reall[,c("Chr","Start","Ref","Alt","sampleID","Disease","Gene","AAChange.metaSVM","Category","Category_metaSVM","RadialSVM_pred_metaSVM","mutation_type_metaSVM")]
    write.csv(reallnew1,file="DDD_mutations/datasheet/anno_mutationlist_4_15.csv",row.names=FALSE)
        
}

Clear_misDDD_metaSVM <- function(){
    
    reann <- read.csv("DDD_mutations/datasheet/Annotated_DDDmutationlist_4_15_metaSVM.csv")
    reall <- read.csv("DDD_mutations/datasheet/DDDmutationlist_4_15_all.csv") 
    
    reall <- cbind(reall,reann[,c("Func.refgene","Gene.refgene","ExonicFunc.refgene","RadialSVM_pred","AAChange.refgene")])
    reall <- cbind(reall,"")
    colnames(reall) <- c("Chr","Start","End","Ref","Alt","Gene","Category","sampleID","AAchanges","Disease","From","Func.annovar","Gene.annovar","Category_annovar","Polyphen2_HDIV_pred_annovar","AAChange.annovar","mutation_type_annovar","Func.metaSVM","Gene.metaSVM","Category_metaSVM","RadialSVM_pred_metaSVM","AAChange.metaSVM","mutation_type_metaSVM")
    
    ## table(reann[,"ExonicFunc.refgene"])
    # c("","frameshift deletion","frameshift deletion;frameshift deletion","frameshift insertion","nonframeshift deletion","nonframeshift insertion","nonsynonymous SNV","nonsynonymous SNV;nonsynonymous SNV","nonsynonymous SNV;nonsynonymous SNV;nonsynonymous SNV","nonsynonymous SNV;nonsynonymous SNV;nonsynonymous SNV;nonsynonymous SNV","stopgain","stoploss","synonymous SNV","synonymous SNV;synonymous SNV","unknown")
    
    #subs <- reann[,"ExonicFunc.refgene"]==""
    #table(reann[subs,"Func.refgene"])
    ##c("","downstream","intergenic","intronic","ncRNA_exonic","ncRNA_intronic","splicing","UTR3") 
    
    annlof <- c("frameshift deletion","frameshift deletion;frameshift deletion","frameshift insertion","stopgain","stoploss")
    subslof <- reann[,"ExonicFunc.refgene"] %in% annlof | (reann[,"ExonicFunc.refgene"]=="" & (reann[,"Func.refgene"] %in% c("splicing")))
    #subslof1  <- reann[,"ExonicFunc.refgene"]=="" & reann[,"Func.refgene"]=="" ## table(reall0[subslof1,"Category"])
    subsmisd <- (grepl("nonsynonymous SNV",reann[,"ExonicFunc.refgene"]) | grepl("nonframeshift deletion",reann[,"ExonicFunc.refgene"]) | grepl("nonframeshift insertion",reann[,"ExonicFunc.refgene"])) & reann[,"RadialSVM_pred"]=="D"
    subsmis <- (grepl("nonsynonymous SNV",reann[,"ExonicFunc.refgene"]) | grepl("nonframeshift deletion",reann[,"ExonicFunc.refgene"]) | grepl("nonframeshift insertion",reann[,"ExonicFunc.refgene"])) & reann[,"RadialSVM_pred"]!="D"
    subssyn <- grepl("synonymous SNV",reann[,"ExonicFunc.refgene"]) & !grepl("nonsynonymous SNV",reann[,"ExonicFunc.refgene"])
    subsintron <- reann[,"ExonicFunc.refgene"]=="" & !(reann[,"Func.refgene"] %in% c("splicing"))
    subsunkown <- reann[,"ExonicFunc.refgene"]=="unknown"
    
    reall[subslof,"mutation_type_metaSVM"] <- "LOF"
    reall[subsmisd,"mutation_type_metaSVM"] <- "dMIS"
    reall[subsmis,"mutation_type_metaSVM"] <- "MIS"
    reall[subssyn,"mutation_type_metaSVM"] <- "SYNONMOUS"
    reall[subsintron,"mutation_type_metaSVM"] <- "Other"
    
    subsunkown <- reann[,"ExonicFunc.refgene"]=="unknown"
    subsunsyn <- reann[,"ExonicFunc.refgene"]=="unknown" & reall[,"Category"] %in% c("synonymous_variant")
    subsundmis <- reann[,"ExonicFunc.refgene"]=="unknown" & reall[,"Category"]== "missense_variant" & reann[,"RadialSVM_pred"]=="D"
    subsunmis <- reann[,"ExonicFunc.refgene"]=="unknown" & reall[,"Category"]== "missense_variant" & reann[,"RadialSVM_pred"]!="D"
    
    reall[subsunsyn,"mutation_type_metaSVM"] <- "SYNONMOUS"
    reall[subsunmis,"mutation_type_metaSVM"] <- "MIS"
    reall[subsundmis,"mutation_type_metaSVM"] <- "dMIS"
    reall <- reall[,c("Chr","Start","End","Ref","Alt","Gene","Category","sampleID","AAchanges","Disease","From","Func.annovar","AAChange.annovar","Gene.annovar","Category_annovar","Polyphen2_HDIV_pred_annovar","mutation_type_annovar","Func.metaSVM","AAChange.metaSVM","Gene.metaSVM","Category_metaSVM","RadialSVM_pred_metaSVM","mutation_type_metaSVM")]
    
    write.csv(reall,file="DDD_mutations/datasheet/anno_DDDmutationlist_4_15_all.csv",row.names=FALSE)
    
    reallnew1 <- reall[,c("Chr","Start","Ref","Alt","sampleID","Disease","Gene","AAChange.metaSVM","Category","Category_metaSVM","RadialSVM_pred_metaSVM","mutation_type_metaSVM")]
    write.csv(reallnew1,file="DDD_mutations/datasheet/anno_DDDmutationlist_4_15.csv",row.names=FALSE)
    
}

Clear_miscontrol_metaSVM <- function(){
    
    reann <- read.csv("DDD_mutations/datasheet/controls.csv")
    reall <- reann[,c("CHR","POS","POS","REF","ALT","Gene","ExonicFunc.refgene","Subject","AAChange.refgene","RadialSVM_pred")]
    reall[,"Disease"] <- "Control"
    reall[,"From"] <- "SSC"
    colnames(reall) <- c("Chr","Start","End","Ref","Alt","Gene","Category","sampleID","AAchanges","RadialSVM_pred_metaSVM","Disease","From")
    ## table(reann[,"ExonicFunc.refgene"])
    # c("","frameshift deletion","frameshift insertion","nonframeshift deletion","nonframeshift insertion","nonsynonymous SNV","stopgain","stopgain;stopgain","stoploss","synonymous SNV","unknown")
    #subs <- reann[,"ExonicFunc.refgene"]==""
    #table(reann[subs,"Func.refgene"])
    ##c("splicing") 
    annlof <- c("frameshift deletion","frameshift insertion","stopgain","stopgain;stopgain","stoploss")
    subslof <- reann[,"ExonicFunc.refgene"] %in% annlof | (reann[,"ExonicFunc.refgene"]=="" & (reann[,"Func.refgene"] %in% c("splicing")))
    subsmisd <- (grepl("nonsynonymous SNV",reann[,"ExonicFunc.refgene"]) | grepl("nonframeshift deletion",reann[,"ExonicFunc.refgene"]) | grepl("nonframeshift insertion",reann[,"ExonicFunc.refgene"])) & reann[,"RadialSVM_pred"]=="D"
    subsmis <- (grepl("nonsynonymous SNV",reann[,"ExonicFunc.refgene"]) | grepl("nonframeshift deletion",reann[,"ExonicFunc.refgene"]) | grepl("nonframeshift insertion",reann[,"ExonicFunc.refgene"])) & reann[,"RadialSVM_pred"]!="D"
    subssyn <- grepl("synonymous SNV",reann[,"ExonicFunc.refgene"]) & !grepl("nonsynonymous SNV",reann[,"ExonicFunc.refgene"])
    subsintron <- reann[,"ExonicFunc.refgene"]=="" & !(reann[,"Func.refgene"] %in% c("splicing","splicing;splicing",""))
    subsunkown <- reann[,"ExonicFunc.refgene"]=="unknown"
    
    reall[subslof,"mutation_type_metaSVM"] <- "LOF"
    reall[subsmisd,"mutation_type_metaSVM"] <- "dMIS"
    reall[subsmis,"mutation_type_metaSVM"] <- "MIS"
    reall[subssyn,"mutation_type_metaSVM"] <- "SYNONMOUS"
    reall[subsintron,"mutation_type_metaSVM"] <- "Other"
    reall[subsunkown,"mutation_type_metaSVM"] <- "unknown"
    
    write.csv(reall,file="DDD_mutations/datasheet/anno_mutationlist_control.csv",row.names=FALSE)
}

clear_ncounttable <- function(){
    source("DDD.R")
    tmp <- clear_number_metaSVM("DDD_mutations/datasheet/anno_mutationlist_4_15.csv")
    write.table(tmp,file="DDD_mutations/datasheet/mutation_4_15_table1.txt",quote=FALSE,sep="\t")
    
    tmp <- clear_number_metaSVM("DDD_mutations/datasheet/anno_DDDmutationlist_4_15.csv")
    write.table(tmp,file="DDD_mutations/datasheet/DDDmutation_4_15_table1.txt",quote=FALSE,sep="\t")
    
    tmp <- clear_number_metaSVM("DDD_mutations/datasheet/anno_mutationlist_control.csv")
    write.table(tmp,file="DDD_mutations/datasheet/control_mutation_table1.txt",quote=FALSE,sep="\t")
}

clear_number_metaSVM <- function(filename){
    tmp <- read.csv(filename)
    
    genes <- unique(tmp[,"Gene"])
    dis <- unique(tmp[,"Disease"])
    II <- length(dis)
    Mut <- c("LOF","dMIS","MIS","SYNONMOUS")
    JJ <- length(Mut)
    
    numberT <- matrix(0,length(genes),II*JJ)
    rownames(numberT) <- genes
    colnames(numberT) <- rep(c("number_of_LGD","number_of_D-mis","number_of_B-mis","number_of_silent"),II)
    for(i in 1:II){
        for(j in 1:JJ){
            tmpdata <- tmp[tmp[,"Disease"]==dis[i] & tmp[,"mutation_type_metaSVM"]==Mut[j],]
            tmpg <- names(table(tmpdata[,"Gene"]))
            numberT[match(tmpg,genes),((i-1)*JJ + j)] <- table(tmpdata[,"Gene"])
        }
    }
    
    numberT
    
}

clear_number <- function(filename){
    tmp <- read.csv(filename)
    
    genes <- unique(tmp[,"Gene"])
    dis <- unique(tmp[,"Disease"])
    II <- length(dis)
    Mut <- c("LOF","dMIS","MIS","SYNONMOUS")
    JJ <- length(Mut)
    
    numberT <- matrix(0,length(genes),II*JJ)
    rownames(numberT) <- genes
    colnames(numberT) <- rep(c("number_of_LGD","number_of_D-mis","number_of_B-mis","number_of_silent"),II)
    for(i in 1:II){
        for(j in 1:JJ){
            tmpdata <- tmp[tmp[,"Disease"]==dis[i] & tmp[,"mutation_type"]==Mut[j],]
            tmpg <- names(table(tmpdata[,"Gene"]))
            numberT[match(tmpg,genes),((i-1)*JJ + j)] <- table(tmpdata[,"Gene"])
        }
    }
    
    numberT

}

metana <- function(){

    source("ASD_data_set.R")
    mutrate0 <- read.csv("ASD/GeneMutaRatem.csv")
    mutrate1 <- read.csv("ASD/MutationRatem.csv")
    mutrate <- cbind(mutrate0,mutrate1[match(mutrate0[,1],mutrate1[,1]),"mut.rate"])
    colnames(mutrate)[6] <- "mut.rate"
    
    data <- read.csv("DDD_mutations/datasheet/mutation_4_15_table1.csv",skip=1)
    data <- data[,c(1,14:17)]
    subs <- data[,1] %in% mutrate[,1]
    data <- data[subs,]
    subs <- match(data[,1],mutrate[,1])
    data1 <- cbind(data[,1],mutrate[subs,2:6],data[,2:5])
    colnames(data1) <- c("Gene","dmis","mis","LOF","LOFDmis","mut.rate","dn.LoF","dn.mis3","dn.mis","dn.syn")
    write.csv(data1,file="DDD_mutations/datasheet/meta_ana.csv",row.names=FALSE)
    
    source("Poisson_test_hq.R")
    filename <- "DDD_mutations/datasheet/meta_ana.csv"
    ntrio <- 4409
    tmp <- Poisson_test_hq(filename,ntrio)
    tmp <- tmp[,c(1:13,16)]
    tmp[,"min_p"] <- apply(tmp[,11:13],1,min)
    tmp <- tmp[order(tmp[,"min_p"]),]
    write.csv(tmp,file="DDD_mutations/datasheet/Ptest_4_22.csv",row.names=FALSE)
    
    
    ### control =========================
    source("ASD_data_set.R")
    mutrate0 <- read.csv("ASD/GeneMutaRatem.csv")
    mutrate1 <- read.csv("ASD/MutationRatem.csv")
    mutrate <- cbind(mutrate0,mutrate1[match(mutrate0[,1],mutrate1[,1]),"mut.rate"])
    colnames(mutrate)[6] <- "mut.rate"
    
    data <- read.csv("DDD_mutations/datasheet/control_mutation_table1.csv",skip=1)
    subs <- data[,1] %in% mutrate[,1]
    data <- data[subs,]
    subs <- match(data[,1],mutrate[,1])
    data1 <- cbind(data[,1],mutrate[subs,2:6],data[,2:5])
    colnames(data1) <- c("Gene","dmis","mis","LOF","LOFDmis","mut.rate","dn.LoF","dn.mis3","dn.mis","dn.syn")
    write.csv(data1,file="DDD_mutations/datasheet/control_ana.csv",row.names=FALSE)
    
}