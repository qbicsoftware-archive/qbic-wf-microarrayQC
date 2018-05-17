######################################
## Quality control analysis report of Microarray Data 


######################################
#load packages needed for oligo run
library("siggenes")
library("RColorBrewer")
library("multtest")
library("limma")
library("oligo")
library("genefilter")
library("gplots")
library("ggplot2")
library("dendextend")
library("statmod")
library("annotate")


######################################
#clean up
rm(list = ls(all = TRUE)) # clear all variables
graphics.off()
path <- getwd()
path


######################################
# create directories needed
dir.create("qc_results")
dir.create("qc_results/raw_data")
dir.create("qc_results/plots")
dir.create("qc_results/tables")
dir.create("qc_results/metadata")
dir.create("qc_results/final")
######################################


######################################
#check raw data
#put raw data (.cel files) raw data/
files = list.files(path,recursive = T, pattern = ".CEL")
for (i in files)
{
  cmd=paste("cp ", paste(i), " qc_results/raw_data/",sep="")
  system(cmd)
}

setwd(path)


######################################
#metadata table
m <- list.files(path, recursive = T, pattern = "sample_preparation")
for (i in m)
{
  cmd=paste("cp ", paste(i), " qc_results/metadata/",sep="")
  system(cmd)
}

m <- list.files("qc_results/metadata/")
m <- read.table(m,  header = T,sep = "\t",na.strings =c("","NaN"),quote=NULL,stringsAsFactors=F,dec=".",fill=TRUE)
# create a new column `x` with all the columns collapsed together
cols <- names(m)
m$x <- apply(m[ ,cols],1,paste, collapse = "-")
m$x <- gsub(" ","",m$x)

#use filename to create a vector of regexes
files <- list.files("qc_results/raw_data/",recursive = T, pattern = ".CEL")
files <- gsub(".CEL","",files)

m$xx <- 0
for (i in files) {
  ifelse(m[grepl(i, m$x), "xx"] <- i,0)
}
m <- subset(m, !xx == 0)
m$xx <- paste(m$xx,".CEL",sep="")
m$filenames <- m$xx
m$xx <- NULL
m$x <- NULL

write.table(m, "qc_results/metadata/metadata.txt",append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = F,  col.names = T, qmethod = c("escape", "double"))

setwd(path)


######################################
#Metadata preparation
pd = read.AnnotatedDataFrame("qc_results/metadata/metadata.txt", header = TRUE)
pd
row.names(pd) = pd[["filenames"]]
row.names(pd)


######################################
#load raw data with that additonal piece of information
setwd(path)
celFiles <- list.celfiles("qc_results/raw_data",full.names=T)
data <- read.celfiles(celFiles,phenoData = pd,verbose = T)  
#changed to filenames, more secure then sampleNames



######################################
#a check:
identical(sampleNames(data),row.names(pd))


#fit Probe Level Models (PLMs) with probe- level and sample-level parameters.
#The resulting object is an oligoPLM object, which stores parameter estimates, residuals and weights.
Pset <- fitProbeLevelModel(data)


##################
#Diagnostic plots and analysis
##################
####set-up color code####
grps <- as.character(pData(data)$treatment)
grps <- as.factor(grps)
grps
display.brewer.all()
col = brewer.pal(length(levels(grps)),"Set1")
col = c(c(col)[grps])
col


#1)Correlation analysis between the samples and their repeated measurements
cor <- cor(exprs(data), use = "everything",method = c("pearson"))
write.table(cor, "qc_results/tables/pearson_correlation_all_data.tsv", append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,  col.names = NA, qmethod = c("escape", "double"))
#use that for matrix eventually such as in the DESeq2 package for RNA-Seq experiments


#2)boxplot of raw log intensities
pdf("qc_results/plots/Boxplot_raw_intensities.pdf")
par(oma=c(12,3,3,3))
par(mfrow = c(1,2))
boxplot(data, which='all', col=col,xlab="", main="", ylab="log2 signal intensity (PM+bg)", cex.axis=0.5, las=2)
legend("topright",col=levels(factor(col)),lwd=1,cex=0.5, legend=levels(grps))
boxplot(data, which='pm', col=col,xlab="", main="", ylab="log2 signal intensity (PM only)", cex.axis=0.5, las=2)
#boxplot(data, which='mm', col=col, xlab="", main="", ylab="log2 signal intensity (MM only)", cex.axis=0.5, las=2)
legend("topright",col=levels(factor(col)),lwd=1,cex=0.5, legend=levels(grps))
mtext("Visualization of raw data using Boxplots of log2 transformed intensity values. The Legend depicts colouring based on sample groups.\nPM=Perfect Match, bg=background", outer = T,side=1,cex=0.5,adj = 0)
dev.off()


#3)RLE and NUSE plots on dataset
#RLE:relative log expression
#NUSE:normalized unscaled standard error 
#In the NUSE plot, low-quality arrays are those that are significantly elevated or more spread out, relative to the other arrays. NUSE values are useful for comparing arrays within one dataset, but their magnitudes are not comparable across different datasets.
#In the RLE plot (Figure 3.3, bottom), problematic arrays are indicated by larger spread, by a center location different from y = 0, or both.

pdf("qc_results/plots/NUSE_plot.pdf")
par(oma=c(12,3,3,3)) 
NUSE(Pset, main="NUSE",ylim=c(0.5,2),outline=FALSE,col=col,las=2,cex.axis=0.5,ylab="Normalized Unscaled Error (NUSE) values",whisklty="dashed",staplelty=1,cex.axis=0.75)
dev.off()
pdf("qc_results/plots/RLE_plot.pdf")
par(oma=c(12,3,3,3)) 
RLE(Pset, main="RLE", ylim = c(-8, 8), outline = FALSE, col=col,las=2, cex.axis=0.5,ylab="Relative Log Expression (RLE) values",whisklty="dashed", staplelty=1,cex.axis=0.75)
dev.off()


#4)Histogram to compare log2 intensities vs density between arrays
#density plots of log base 2 intensities (log2(PMij) for array i and probe j) of perfect match probes for comparison of probe intensity behavior between different arrays. If you see differences in shape or center of the distributions, it means that normalization is required.
pdf("qc_results/plots/Histogramm_log2_intensities_vs_density.pdf")
hist(data,col = col, lty = 1, xlab="log2 intensity", ylab="density", xlim = c(2, 12), type="l")
legend("topright",col = col, lwd=1, legend=sampleNames(data),cex=0.5)
dev.off()


#5)MA plots raw data
#The MAplot also allows summarization, so groups can be compared more easily:

pdf("qc_results/plots/MA_plot_before_normalization_groups.pdf")
MAplot(data, pairs=TRUE, groups=grps,na.rm=TRUE)
dev.off()


#6)PCA plot before normalization
#You want to see which genes that mean the most for the differences between the samples, and therefore your samples should be in the rows and your genes should be in the columns.
#there t():
pca_before <- prcomp(t(exprs(data)), scores=TRUE, scale. = TRUE, cor=TRUE)
summary(pca_before)
# sqrt of eigenvalues
pca_before$sdev
#loadings
head(pca_before$rotation)
#PCs (aka scores)
head(pca_before$x)

# create data frame with scores
scores_before = as.data.frame(pca_before$x)
# plot of observations

#reorder grps just for pca plot
grps_pca <- grps

pdf("qc_results/plots/PCA_before_normalization.pdf")
ggplot(data = scores_before, aes(x = PC1, y = PC2,colour=grps_pca)) +
  geom_hline(yintercept = 0, colour = "gray65") +
  geom_vline(xintercept = 0, colour = "gray65") +
  #geom_text(colour = "black",label=sampleNames(data), size = 2,angle=40) +
  #scale_fill_manual(values=c("#E41A1C", "#377EB8", "#4DAF4A"), breaks=c("Parabel", "Simbox", "Texus"), labels=c("Parabel", "Simbox", "Texus")) +
  geom_point(aes(shape = factor(data$treatment)),size=2) + 
  #scale_colour_manual(values = c("#E41A1C","#377EB8", "#4DAF4A"))
  #scale_shape_manual(values=1:nlevels(col)) +
  theme(legend.title=element_blank()) +  ## turn off legend title
  ggtitle("PCA plot before normalization")
dev.off()

##--end QC before normalization
setwd(path)



##################
#Data Normalization
##################
#as default do rma here....
eset <- rma(data)  #depending on chip such as ST chips, default used is target="core"
#it does
#Background correcting
#Normalizing
#Calculating Expression

#Currently the rma function implements RMA in the following manner
#1. Probe specific correction of the PM probes using a model based on observed intensity being the sum of signal and noise
#2. Normalization of corrected PM probes using quantile normalization (Bolstad et al., 2003)
#3. Calculation of Expression measure using median polish.


##########################
#Quality plots
#########################

#1)boxplot after normalization
pdf("qc_results/plots/Boxplot_after_normalization.pdf")
par(oma=c(10,2,2,2))
boxplot(exprs(eset), col=col,which='both', xlab="", main="", ylab="log2 signal intensity", cex.axis=0.4, las=2)
dev.off()

#2)Scatter matrix of arrays against one another
png("out/Scatter_plot_after_normalization.png")
scatter <- pairs(exprs(eset), pch=".",main="Scatter plots", cex=0.5)
dev.off()


#3)The MAplot also allows summarization, so groups can be compared more easily:
pdf("qc_results/plots/MA_plot_after_normalization_groups.pdf")
MAplot(exprs(eset), pairs=TRUE, groups=grps)
dev.off()


#4)clustering
#for arrays, problem: arrays are not row names but at column position, thus transpose is needed
d <- dist(t(exprs(eset))) # find distance matrix
hc <- hclust(d)               # apply hierarchical clustering
##
dend <- as.dendrogram(hc)
#remember groups, assign a new color code
labels_colors(dend) <- col[order.dendrogram(dend)]
#colorCodes = c("red", "blue", "green")
pdf("qc_results/plots/Cluster_Dendogram.pdf")
par(oma=c(10,2,2,2))
dend %>% set("labels_cex",0.5) %>% plot()
legend("topright",col=levels(factor(col)),lwd=1,cex=0.5, legend=levels(grps))
dev.off()


#5)PCA after normalization
pca <- prcomp(t(exprs(eset)), scores=TRUE, cor=TRUE)
summary(pca)
# sqrt of eigenvalues
pca$sdev
#loadings
head(pca$rotation)
#PCs (aka scores)
head(pca$x)

# create data frame with scores
scores_after = as.data.frame(pca$x)
# plot of observations
pdf("qc_results/plots/PCA_after_normalization.pdf")
ggplot(data = scores_after, aes(x = PC1, y = PC2,colour=grps_pca)) +
  geom_hline(yintercept = 0, colour = "gray65") +
  geom_vline(xintercept = 0, colour = "gray65") +
  #geom_text(colour = "black",label=sampleNames(data), size = 2,angle=40) +
  #scale_fill_manual(values=c("#E41A1C", "#377EB8", "#4DAF4A"), breaks=c("Parabel", "Simbox", "Texus"), labels=c("Parabel", "Simbox", "Texus")) +
  geom_point(aes(shape = factor(data$treatment)),size=2) + 
  #scale_colour_manual(values = c("#E41A1C","#377EB8", "#4DAF4A"))
  #scale_shape_manual(values=1:nlevels(col)) +
  theme(legend.title=element_blank()) +  ## turn off legend title
  ggtitle("PCA plot after normalization")
dev.off()



#6)shows scree plot to verify plotting of PC1 vs PC2
library("affycoretools")
pdf("qc_results/plots/PCs.pdf")
plotPCA(exprs(eset),main="Principal component analysis (PCA)", screeplot=TRUE, outside=TRUE)
dev.off()

#unload affy related packages again as analysis is focused on using oligo package function:
detach("package:affycoretools", unload=TRUE)
#detach("package:affy", unload=TRUE)


##---end QC after normalization



#write to file and save esets
write.exprs(eset, "qc_results/tables/RMAnorm_nonfiltered.txt")
#write.exprs(eset_filt_sds, "qc_results/tables/RMAnorm_sds.filtered.txt")
save(list=c("pd","eset"), file="qc_results/final/eset.Rdata")
setwd(path)


#end of script
####-------------save Sessioninfo
fn <- paste("qc_results/tables/sessionInfo_",format(Sys.Date(), "%d_%m_%Y"),".txt",sep="")
sink(fn)
sessionInfo()
sink()
####---------END----save Sessioninfo 








