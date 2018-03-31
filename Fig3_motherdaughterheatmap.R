###Code for a heatmap. Modeled after code written by Anna Seekatz on her GitHub.com account
###under Rcode_erinsubset/erinsubet_Fig3.heatmap.R 

library(RColorBrewer)
library(ggplot2)
library(plyr)
library(gplots)

setwd("~/Desktop/workingmotherdaughterscriptfiles/")
shared<-read.table(file = "motherdaughter.files.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.shared", header = TRUE, row.names = 2)
dim(shared)
otu<-subset(shared, select =-c(label, numOtus))
otu.filtered<-otu[,which(colSums(otu)>=200)]  
otu.filtered<- otu.filtered[ order(row.names(otu.filtered)), ]   
otu.matrix<-as.matrix(otu.filtered)
otu.rel<-otu.matrix/rowSums(otu.matrix)
otu.rel.max<-apply(otu.rel,2,max)
otu.rel.filtered<-otu.rel[,otu.rel.max>0.02]
colorscale <- c("#4d4d4d","#878787","#bababa","#e0e0e0","#fddbc7","#f4a582","#d6604d","#b2182b")
myBreaks2 <- c(0, 0.001, 0.003, 0.01, 0.05, 0.10, 0.50, 0.80, 1) 
heatmap.2(otu.rel.filtered, dendrogram='none', col=colorscale, breaks=myBreaks2, Rowv=TRUE, Colv=TRUE, trace='none')

meta1<-read.table(file="motherdaughterquestionairre_data_meta1.txt", header=TRUE)
subjects.and.samples<-read.table(file = "subject.sample.pair.txt", header = TRUE)
meta<-merge(meta1, subjects.and.samples, by.x=c("Subject"), by.y=c("Subject"))
#write.table(meta, "motherdaughter.meta.heatmap.txt", quote = FALSE, sep = "\t", col.names = NA) 

meta<-read.table(file="motherdaughter.meta.heatmap.txt", header = TRUE) #rewrites meta table
meta.otu.rel.filtered<-merge(otu.rel.filtered, meta, by.x=c("row.names"), by.y=c("Sample"))
meta3<-meta.otu.rel.filtered
colnames(meta3)[1]<- "Sample"
#write.table(meta.otu.rel.filtered, "motherdaughter.otu.meta.heatmap.txt", quote=FALSE, sep = "\t", col.names = NA)

###made motherdaughter.otu.heatmap.heatorder.txt in excel with 
###motherdaughter.otu.heatmap.ordered.wcolumns.txt as the base file, comprising of meta 
###data and relative abundances of the OTUs that had more than 200 sequences. 
###motherdaughter.otu.heatmap.ordered.wcolumns.txt was sorted first by Pairs, then by Subject1, 
###then by timepoint, then by Heatorder. The number for "Heatorder" was assigned to each sample based on 
###based on the daughter's birth mode, gestation time, and the number of samples received from the pair.
###0=vaginal birth, full term, 3 or more weeks collected,
###1=vaginal birth, full term, 1 sample collected, 
###2=vaginal birth, pre term, 3 or more weeks collected,
###3=C-section, pre term, 3 or more weeks collected,
###4=C-section, full term, 3 or more weeks collected,
###5=C-section, full term, 3 or more weeks collected,
###After, I deleted the Sample column, Timepoint column, Pair column, Subject1 column, Subject column
###numbered Rows column, and Heatorder column, leaving only the sample ID as the row and its relative
###abundance in the various OTU columns in motherdaughter.otu.heatmap.heatorder.txt.

meta4<-read.table(file = "motherdaughter.otu.heatmap.ordered.wcolumns.txt", header = TRUE, row.names = 1)
meta2<-read.table(file = "motherdaughter.otu.heatmap.heatorder.wcolumns.txt", header = TRUE, row.names = 1)
otu.ordered.rel.filtered<-read.table("motherdaughter.otu.heatmap.heatorder.txt", header=TRUE, row.names=1)
otu.ordered.rel.filtered.matrix<-as.matrix(otu.ordered.rel.filtered)

###Made motherdaughter.files.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.0.03.cons.taxonomy.subset.txt with Alyx Schubert's code and manipulation in Excel

###Makes column/taxonomy labels
tax<-read.table(file="motherdaughter.files.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.0.03.cons.taxonomy.subset.txt", header=TRUE)
tax<- tax[order(tax$GROUP),]
tax.col<-tax[,c("GROUP", "Phylum_order", "Size")]
tax.col$color<-mapvalues(tax.col$Phylum_order, from = c("01_Firmicutes","02_Firmicutes", "03_Actinobacteria", "04_Actinobacteria", "05_Bacteroidetes", "06_Bacteroidetes", "07_Fusobacteria", "08_Tenericutes"), to = c("#0c2c84","#225ea8","#1d91c0","#41b6c4","#7fcdbb","#c7e9b4","#edf8b1","#ffffd9")) 
tax.col<- tax.col[order(tax.col$Phylum_order, -tax.col$Size) , ]
rows<-tax.col[,1]
rownames(tax.col)<-rows
tax.col<-t(tax.col)
phyla.col<-tax.col[4,] 
col.order<-as.character(tax.col[1, ]) 

heatmap.2(otu.ordered.rel.filtered.matrix, col=colorscale, breaks=myBreaks2, cexRow=0.5, cexCol=0.5, trace="none", tracecol="black", rowsep=c(1, 9, 18, 26, 36, 46, 56, 64, 66, 68, 78, 88, 98, 108), density.info="none", Colv=F, Rowv=F, labCol = rows, dendrogram="none", ColSideColors = as.character(phyla.col))

###Makes row labels
meta2$color<-mapvalues(meta2$Pair, from = c("CONT_0", "PAIR_A", "PAIR_B", "PAIR_D", "PAIR_E", "PAIR_F", "PAIR_J", "PAIR_L", "PAIR_N", "PAIR_P", "PAIR_Q", "PAIR_S", "PAIR_T", "PAIR_U"), to = c("#F0DE14", "#D2E7D5", "#f768a1", "#e5f5e0", "#006d2c", "#c7e9c0", "#9e9ac8", "#ae017e", "#a1d99b", "#6a51a3", "#74c476", "#41ab5d", "#238b45", "#00441b"))
heatmap.2(otu.ordered.rel.filtered.matrix, col=colorscale, breaks=myBreaks2, cexRow=0.5, cexCol=0.5, trace="none", tracecol="black", rowsep=c(1, 9, 18, 26, 36, 46, 56, 64, 66, 68, 78, 88, 98, 108), density.info="none", Colv=FALSE, Rowv=FALSE, labCol = rows, RowSideColors=as.character(meta2$color), ColSideColors = as.character(phyla.col), dendrogram="none")

#heatmap.2(otu.ordered.rel.filtered.matrix, col=colorscale, breaks=myBreaks2, cexRow=0.5, cexCol=0.5, trace="none", dendrogram="none") #clears first row label 

#meta2$relation<-mapvalues(meta2$Subject1, from = c("Mother", "Daughter", "na"), to = c("darkred", "darksalmon", "cyan2"))
#heatmap.2(otu.ordered.rel.filtered.matrix, col=colorscale, breaks=myBreaks2, cexRow=0.5, cexCol=0.5, trace="none", tracecol="black", density.info="none", Colv=FALSE, Rowv=FALSE, labCol=rows, dendrogram="none", RowSideColors = as.character(meta2$relation), ColSideColors = as.character(phyla.col))

###To get rid of row labels (sample names):
#heatmap.2(otu.ordered.rel.filtered.matrix, col=colorscale, breaks=myBreaks2, cexRow=0.5, cexCol=0.5, trace="none", density.info="none", labRow=F, Colv=F, Rowv=F, dendrogram="none", RowSideColors = as.character(meta2$relation), ColSideColors = as.character(phyla.col))

###To get colors for key:
#display.brewer.pal(8, "Blues")
