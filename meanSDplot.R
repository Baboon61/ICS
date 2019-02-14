#!/usr/bin/env Rscript
library("DESeq2")
library("vsn")

#Create matrix
data <- as.matrix(read.csv("/media/sf_raid/Data/ICS/bedgraph_normalized/merged.normalized.bdg", row.names=1, sep="\t", header=TRUE))

#Create design
data_coldata <- read.csv("/media/sf_raid/Data/ICS/bedgraph_normalized/design_MA_plot.txt", row.names=1, sep="\t", header=TRUE)


run_meanSDplot <- function(cts, coldata,aB,rB,aB_input,rB_input,name,threshold){

	#Apply input
    cts[, aB] <- (cts[, aB] - cts[, aB_input])
    cts[, rB] <- (cts[, rB] - cts[, rB_input])

    #Remove input columns
    cts <- cts[,-c(1,2)]

    #Round float
    cts <-round(cts,0)

    #Apply 1 to value < 1, due to input files
    cts[cts < 0] <- 0
    cts <- cts[(cts[,aB] > 0) | (cts[,rB] > 0),]

    print(all(rownames(coldata) == colnames(cts)))

    dds <- DESeqDataSetFromMatrix(countData = cts,colData = coldata,design = ~ condition)                       

    rld <- rlog(dds)

    png(filename=paste0(c("/media/sf_raid/Data/ICS/meanSDplot/",name,".png"), collapse=''), width = 4, height = 4, units = 'in', res = 300)
    meanSdPlot(assay(rld), ranks = FALSE)
    dev.off()

    best_bins <- assay(rld)[(assay(rld)[,aB]+assay(rld)[,rB])/2 >= threshold,]
    write.table(rownames(best_bins), paste0(c("/media/sf_raid/Data/ICS/bed_high_coverage/",name,".bed"), collapse=''), row.names = FALSE, col.names=FALSE, sep = "\t")
    

}

run_meanSDplot(data[,c("aB24h_wt_ChIP_input", "rB_wt_ChIP_input", "aB_wt_H2AZ", "rB_wt_H2AZ")], data_coldata[c("aB_wt_H2AZ", "rB_wt_H2AZ"),], "aB_wt_H2AZ", "rB_wt_H2AZ", "aB24h_wt_ChIP_input", "rB_wt_ChIP_input", "H2AZ", 3)

run_meanSDplot(data[,c("aB24h_wt_ChIP_input", "rB_wt_ChIP_input", "aB_wt_H3K27Ac", "rB_wt_H3K27Ac")], data_coldata[c("aB_wt_H3K27Ac", "rB_wt_H3K27Ac"),], "aB_wt_H3K27Ac", "rB_wt_H3K27Ac", "aB24h_wt_ChIP_input", "rB_wt_ChIP_input", "H3K27Ac", 3)

run_meanSDplot(data[,c("aB24h_wt_ChIP_input", "rB_wt_ChIP_input", "aB_wt_H3K36me1", "rB_wt_H3K36me1")], data_coldata[c("aB_wt_H3K36me1", "rB_wt_H3K36me1"),], "aB_wt_H3K36me1", "rB_wt_H3K36me1", "aB24h_wt_ChIP_input", "rB_wt_ChIP_input", "H3K36me1", 3)
run_meanSDplot(data[,c("aB24h_wt_ChIP_input", "rB_wt_ChIP_input", "aB_wt_H3K36me2", "rB_wt_H3K36me2")], data_coldata[c("aB_wt_H3K36me2", "rB_wt_H3K36me2"),], "aB_wt_H3K36me2", "rB_wt_H3K36me2", "aB24h_wt_ChIP_input", "rB_wt_ChIP_input", "H3K36me2", 3)
run_meanSDplot(data[,c("aB24h_wt_ChIP_input", "rB_wt_ChIP_input", "aB_wt_H3K36me3", "rB_wt_H3K36me3")], data_coldata[c("aB_wt_H3K36me3", "rB_wt_H3K36me3"),], "aB_wt_H3K36me3", "rB_wt_H3K36me3", "aB24h_wt_ChIP_input", "rB_wt_ChIP_input", "H3K36me3", 3)

run_meanSDplot(data[,c("aB24h_wt_ChIP_input", "rB_wt_ChIP_input", "aB_wt_H3K4me1", "rB_wt_H3K4me1")], data_coldata[c("aB_wt_H3K4me1", "rB_wt_H3K4me1"),], "aB_wt_H3K4me1", "rB_wt_H3K4me1", "aB24h_wt_ChIP_input", "rB_wt_ChIP_input", "H3K4me1", 3)
run_meanSDplot(data[,c("aB24h_wt_ChIP_input", "rB_wt_ChIP_input", "aB_wt_H3K4me3", "rB_wt_H3K4me3")], data_coldata[c("aB_wt_H3K4me3", "rB_wt_H3K4me3"),], "aB_wt_H3K4me3", "rB_wt_H3K4me3", "aB24h_wt_ChIP_input", "rB_wt_ChIP_input", "H3K4me3", 3)

run_meanSDplot(data[,c("aB24h_wt_ChIP_input", "rB_wt_ChIP_input", "aB_wt_H3K9me1", "rB_wt_H3K9me1")], data_coldata[c("aB_wt_H3K9me1", "rB_wt_H3K9me1"),], "aB_wt_H3K9me1", "rB_wt_H3K9me1", "aB24h_wt_ChIP_input", "rB_wt_ChIP_input", "H3K9me1", 3)
run_meanSDplot(data[,c("aB24h_wt_ChIP_input", "rB_wt_ChIP_input", "aB_wt_H3K9me2", "rB_wt_H3K9me2")], data_coldata[c("aB_wt_H3K9me2", "rB_wt_H3K9me2"),], "aB_wt_H3K9me2", "rB_wt_H3K9me2", "aB24h_wt_ChIP_input", "rB_wt_ChIP_input", "H3K9me2", 3)
run_meanSDplot(data[,c("aB24h_wt_ChIP_input", "rB_wt_ChIP_input", "aB_wt_H3K9me3", "rB_wt_H3K9me3")], data_coldata[c("aB_wt_H3K9me3", "rB_wt_H3K9me3"),], "aB_wt_H3K9me3", "rB_wt_H3K9me3", "aB24h_wt_ChIP_input", "rB_wt_ChIP_input", "H3K9me3", 3)

run_meanSDplot(data[,c("aB24h_wt_ChIP_input", "rB_wt_ChIP_input", "aB_wt_H4K16Ac", "rB_wt_H4K16Ac")], data_coldata[c("aB_wt_H4K16Ac", "rB_wt_H4K16Ac"),], "aB_wt_H4K16Ac", "rB_wt_H4K16Ac", "aB24h_wt_ChIP_input", "rB_wt_ChIP_input", "H4K16Ac", 3)

run_meanSDplot(data[,c("aB24h_wt_ChIP_input", "rB_wt_ChIP_input", "aB_wt_H4K20me1", "rB_wt_H4K20me1")], data_coldata[c("aB_wt_H4K20me1", "rB_wt_H4K20me1"),], "aB_wt_H4K20me1", "rB_wt_H4K20me1", "aB24h_wt_ChIP_input", "rB_wt_ChIP_input", "H4K20me1", 3)
run_meanSDplot(data[,c("aB24h_wt_ChIP_input", "rB_wt_ChIP_input", "aB_wt_H4K20me3", "rB_wt_H4K20me3")], data_coldata[c("aB_wt_H4K20me3", "rB_wt_H4K20me3"),], "aB_wt_H4K20me3", "rB_wt_H4K20me3", "aB24h_wt_ChIP_input", "rB_wt_ChIP_input", "H4K20me3", 3)
