DDD_FDR <- function(){
    rat <- read.csv("D_rate.csv")
    mut <- read.csv("DDD_mutations/datasheet/Ptest_4_28.csv")

    ig <- intersect(mut[,1],rat[,2])
    
    mut <- mut[,1:10]
    mut[match(ig,mut[,1]),"LOF"] <- rat[match(ig,rat[,2]),"LOF"]
    mut[match(ig,mut[,1]),"dmis"] <- rat[match(ig,rat[,2]),"DMIS"]
    mut[match(ig,mut[,1]),"mis"] <- rat[match(ig,rat[,2]),"MIS"]
    
    
    write.csv(mut,file="DDD_dmis.csv",row.names=FALSE)
    
}