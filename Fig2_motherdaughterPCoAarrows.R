###This code outlines the commands for a PCoA table by withing individual mother-daughter pairs and mothers versus daughters in general
###list.shared, and corr.axes file. This code was model after code, Rcode_erinsubset/erinsubset_Fig4.pcoa.R,
###written by Anna Seekatz, viewed on her GitHub account.

library(shape)
setwd("~/Desktop/workingmotherdaughterscriptfiles/")

###Subset to get significant OTUs 
tax<-read.table(file="motherdaughter.files.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.0.03.cons.taxonomy", header = TRUE)
corr.axes<-read.table(file="motherdaughter.files.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.shared_minus_control.pearson.corr.axes", header = TRUE)
corr.tax<-merge(corr.axes, tax, by="OTU", all.y=TRUE)
corr.tax.1<-subset(corr.tax, p.value < 0.001)
dim(corr.tax.1) #5 OTUs
corr.tax.2<-subset(corr.tax.1, p.value.1 < 0.001)
dim(corr.tax.2) #2 OTUs

###Make plot

#par(mfrow=c(1,2)) 
pcoa.axes<-read.table(file="motherdaughter.files.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.shared_minus_control.thetayc.0.03.lt.pcoa.axes", header = TRUE)
pcoa.col<-pcoa.col<-c(rep("#D2E7D5", 4), rep("#e5f5e0", 4), rep("#c7e9c0", 4), rep("#a1d99b", 5), rep("#6a51a3", 5), rep("#f768a1", 5), rep("#00441b", 1), rep("#74c476", 5), rep("#238b45", 4), rep("#41ab5d", 5), rep("#9e9ac8", 5), rep("#006d2c", 1), rep("#ae017e", 5), rep("#D2E7D5", 4), rep("#f768a1", 5), rep("#F0DE14", 1), rep("#e5f5e0", 5), rep("#006d2c", 1), rep("#c7e9c0", 4), rep("#9e9ac8", 5), rep("#ae017e", 5), rep("#a1d99b", 5), rep("#6a51a3", 5), rep("#74c476", 5), rep("#41ab5d", 5), rep("#238b45", 4), rep("#00441b", 1))
###Labeling just Mothers, Daughters, Control 
pcoa.pch<-pcoa.pch<-c(rep(16, 53), rep(17, 9), rep(18, 1), rep(17, 45))
plot(pcoa.axes$axis1, pcoa.axes$axis2, col=pcoa.col, pch=pcoa.pch, main= "Community Distances of Mother-Daughter Pairs", xlab="Axis1 (23%)", ylab="Axis 2 (13%)", ylim=c(-0.6, 0.6), xlim=c(-0.6, 0.6), cex=1.3)
legend(x=0.33, y=.5, legend=c("Mothers","Daughters", "Control"), pch=c(17, 16, 18), col=c("black","black","black"), bty=("n"), cex=0.6)

###Labeling Mothers, Daughters-Vaginal Full-Term, Daughters Vaginal Pre-term, Daughters C-section Full-Term, Daughters Pre-term C-section, Control
#pcoa.pch.all<-c(rep(16, 17), rep(0, 5), rep(1, 5), rep(16, 15), rep(15, 5), rep(16, 1), rep(0, 5), rep(17, 9), rep(18, 1), rep(17, 45))
#plot(pcoa.axes$axis1, pcoa.axes$axis2, col=pcoa.col, pch=pcoa.pch.all, main= "Community Distances of Individual Mother-Daughter Pairs", xlab="Axis1 (23%)", ylab="Axis 2 (13%)", ylim=c(-0.6, 0.6), xlim=c(-0.6, 0.6), cex=1.3)
#legend(x=0.33, y=.67, legend=c("Mothers","Daughters Full-Term Vaginal", "Daughters Pre-Term Vaginal", "Daughters Full-Term C-section", "Daughters Pre-Term C-section", "Control"), pch=c(17, 16, 1, 15, 0, 18), col=c("black","black","black"), bty=("n"), cex=0.6)

###**directions for arrows: "0, 0" is where you want the arrow to start/its origin, x1 is the number under axis1 in the corr.tax.1 or .2 table, x2 is the number under axis2 in the corr.tax.1 or .2  
###Arrows(0, 0, x1=0.930454, y1=-0.251720, lty = 1, arr.length = 0.1, arr.type = "triangle")

###Make OTU Arrows
###Arrows(0, 0, x1= , y1= , lty = 1, arr.length = 0.1, arr.type = "triangle") <- base code for new arrows

Arrows(0, 0, x1=-0.4, y1=-0.44, lty = 1, arr.length = 0.1, arr.type = "triangle")
text(-0.4,-0.5, label="OTU 2: \nL. iners", cex=0.6)
###NOTE: for Arrows(0, 0, x1=-0.4, y1=-0.44, lty = 1, arr.length = 0.1, arr.type = "triangle") -0.640645 is the TRUE x1 value and -0.716888 is the true y1 value
Arrows(0, 0, x1=-0.314334, y1=0.492207, lty = 1, arr.length = 0.1, arr.type = "triangle")
text(-0.23, 0.27, label="OTU 8: \nAtopobium", cex=0.6)
Arrows(0, 0, x1=0.54 , y1=-0.15 , lty = 1, arr.length = 0.1, arr.type = "triangle")
text(0.55,-0.20, label="OTU 1: \nL. crispatus", cex = 0.6)
###NOTE: for Arrows(0, 0, x1=0.656432…) and text(-0.5, 0.5, label="OTU 1: \nL.cripatus, cex = 0.6) 0.927861 is the TRUE x1 value and -0.256968 is the true y1 value
###0.656432
Arrows(0, 0, x1=0.318463 , y1=-0.158469 , lty = 1, arr.length = 0.1, arr.type = "triangle")
text(0.25, -0.188469, label="OTU 30: \nClostridia", cex = 0.6) #excluded from final Figure 2
Arrows(0, 0, x1=-0.338637, y1=0.215843, lty = 1, arr.length = 0.1, arr.type = "triangle")
text(-0.338637, 0.155843, label="OTU 11: \nPrevotella", cex = 0.6)

