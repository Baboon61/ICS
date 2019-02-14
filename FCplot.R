#!/usr/bin/env Rscript
library(ggplot2)
FC <- read.csv("/media/sf_raid/Data/ICS/reads_mapped_regions_normalized.txt", header = FALSE, sep = "\t", col.names = c("Mark", "Locus", "NbReads"))

for (i in 5:56){
    if (as.vector(FC$Locus[i]) == "Sgamma3"){
        FC$NbReads[i] <- FC$NbReads[i] - FC$NbReads[1]
    }
    else if (as.vector(FC$Locus[i]) == "Sgamma1"){
        FC$NbReads[i] <- FC$NbReads[i] - FC$NbReads[2]
    }
    else if (as.vector(FC$Locus[i]) == "Smu"){
        FC$NbReads[i] <- FC$NbReads[i] - FC$NbReads[3]
    }
    else if (as.vector(FC$Locus[i]) == "3RR"){
        FC$NbReads[i] <- FC$NbReads[i] - FC$NbReads[4]
    }
    else{
        print("Error locus")
    }
}
for (i in 61:112){
    if (as.vector(FC$Locus[i]) == "Sgamma3"){
        FC$NbReads[i] <- FC$NbReads[i] - FC$NbReads[57]
    }
    else if (as.vector(FC$Locus[i]) == "Sgamma1"){
        FC$NbReads[i] <- FC$NbReads[i] - FC$NbReads[58]
    }
    else if (as.vector(FC$Locus[i]) == "Smu"){
        FC$NbReads[i] <- FC$NbReads[i] - FC$NbReads[59]
    }
    else if (as.vector(FC$Locus[i]) == "3RR"){
        FC$NbReads[i] <- FC$NbReads[i] - FC$NbReads[60]
    }
    else{
        print("Error locus")
    }
}

FC <- FC[-c(1,2,3,4,57,58,59,60), ]

FC[FC < 0] <- 0

FC_count <- FC

FC["FC"] <- 0

for (i in 1:52){
	if (as.vector(FC$NbReads[i]) >= as.vector(FC$NbReads[i+52])){
		if (as.vector(FC$NbReads[i+52]) != 0){
            FC[i,"FC"] <- log2(FC$NbReads[i])-log2(FC$NbReads[i+52])
        }
        else{
        	FC[i,"FC"] <- log2(FC$NbReads[i])
        }
    }
    else{
    	if (as.vector(FC$NbReads[i]) != 0){
            FC[i,"FC"] <- -(log2(FC$NbReads[i+52])-log2(FC$NbReads[i]))
        }
        else{
        	FC[i,"FC"] <- -(log2(FC$NbReads[i+52]))
        }
    }
}

FC <- FC[-c(53:104), ]
FC$Mark <- gsub("aB_wt_", "", as.character(FC$Mark))

for (locus in c("Sgamma3", "Sgamma1", "Smu", "3RR")){
	sub_FC <- FC[FC$Locus == locus,]
	print(sub_FC)
    ggplot(sub_FC) +
    geom_col(aes(reorder(sub_FC$Mark,-sub_FC$FC), sub_FC$FC, fill=sub_FC$FC))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x = element_text(angle = 70, hjust = 1))+
    xlab("Marks")+
    ylab(paste0(c("log2(Fold change in ",locus,")"), collapse=''))+
    geom_hline(yintercept=mean(sub_FC[sub_FC$FC != "-Inf","FC"]), linetype="dashed")+
    geom_hline(yintercept=-mean(sub_FC[sub_FC$FC != "-Inf","FC"]), linetype="dashed")+
    scale_fill_gradient2(midpoint=0, low="blue", mid="white", high="red", name="log2(Fold change)")
    ggsave(paste0(c("/media/sf_raid/Data/ICS/FC_plot/FC_",locus,".png"), collapse=''),plot = last_plot())   
}

FC_count[FC_count <= 0] <- 1
FC_count[,"NbReads"] <- log2(FC_count[,"NbReads"])
FC_count$Mark <- gsub("aB_wt_", "", as.character(FC_count$Mark))
FC_count[53:104,"NbReads"] <- -FC_count[53:104,"NbReads"]
FC_count$Mark <- gsub("rB_wt_", "", as.character(FC_count$Mark))

for (locus in c("Sgamma3", "Sgamma1", "Smu", "3RR")){
	sub_FC_count <- FC_count[FC_count$Locus == locus,]
	print(sub_FC_count)
    ggplot(sub_FC_count) +
    geom_col(aes(reorder(sub_FC_count$Mark,-sub_FC_count$NbReads), sub_FC_count$NbReads, fill=sub_FC_count$NbReads))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x = element_text(angle = 70, hjust = 1))+
    xlab("Marks")+
    ylab(paste0(c("log2(Number of reads in ",locus,")"), collapse=''))+
    geom_hline(yintercept=mean(sub_FC_count[(sub_FC_count$NbReads != "-Inf"),"NbReads"]), linetype="dashed")+
    geom_hline(yintercept=-mean(sub_FC_count[sub_FC_count$NbReads != "-Inf","NbReads"]), linetype="dashed")+
    scale_fill_gradient2(midpoint=0, low="blue", mid="white", high="red", name="log2(Number of reads)")
    ggsave(paste0(c("/media/sf_raid/Data/ICS/FC_plot/NumberOfReads_",locus,".png"), collapse=''),plot = last_plot())   
}