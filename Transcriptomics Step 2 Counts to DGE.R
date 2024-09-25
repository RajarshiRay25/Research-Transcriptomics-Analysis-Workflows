# Summary of Workflow for RNA-seq Data:

# Data Normalization: limma and edgeR help in normalizing RNA-seq data.
# Differential Expression: Both limma and edgeR provide statistical tools for identifying differentially expressed genes.
# Visualization: Glimma, gplots, and RColorBrewer are useful for visualizing the results with MA plots, heatmaps, and interactive plots.
# Annotation: org.Mm.eg.db helps map gene IDs to annotations for mouse datasets.

if (!requireNamespace("BiocManager"))
  install.packages("BiocManager")
BiocManager::install(c("limma","edgeR","Glimma","org.Mm.eg.db","gplots","RColorBrewer","NMF","BiasedUrn"),force = TRUE)
BiocManager::install("GO.db")
library(edgeR)
library(limma)
library(Glimma)
library(org.Mm.eg.db)
library(gplots)
library(RColorBrewer)
library(NMF)
library(GO.db)

# Input the Feature counts data obtained from the shell script file

seqdata <- read.delim("GSE60450_LactationGenewiseCounts.txt",stringsAsFactors = FALSE)
head(seqdata)
View(seqdata) # View the data in another tab


# Reading the sample information corresponding to columns related to gene counts. This is the information of conditions names for each sample

sampleinfo <- read.delim("SampleInfo.txt",stringsAsFactors = TRUE)
head(sampleinfo)
View(sampleinfo)

# Only keep the gene counts and remove other values 

countdata <- seqdata[,-(1:2)]
head(countdata)
View(countdata)

# Row name set Gene IDs - setting the gene name as identifiers at rows

rownames(countdata) <- seqdata[,1]

View(countdata)

# Col name - Verifying the column names as samples for which genes have been sequenced 

colnames(countdata)

# Optional - shorten the column names to first 7 characters - Make it same as the sample info

colnames(countdata) <- substr(colnames(countdata), start=1,stop=7)
head(countdata)
View(countdata)

# Check for same order of colnames of count data and sample info - Checking the same order makes it easy to render information 

table(colnames(countdata) == sampleinfo$SampleName)

# We use edgeR to create DGEList object to take in gene count data 

y <- DGEList(countdata)
y

# Visualise names of various properties in the DGE  object

names(y)

# Visualise the groups by joining them together - OPTIONAL

group <- paste(sampleinfo$CellType,sampleinfo$Status,sep = ".")
group

# Convert to category factor  - Setting categories for DGE analysis with factor function

group <- factor(group)
group

# Adding groups from DGE Object to counts data

y$samples$group <- group
y$samples


# Adding Annotation to gene ids based on the provided IDs

columns(org.Mm.eg.db)
View(y$counts)

# Create annotation object to obtain gene names and information against the IDs in the counts data

ann <- select(org.Mm.eg.db,keys = rownames(y$counts),columns = c("ENTREZID","SYMBOL","GENENAME"))

View(ann)

# Checking the order of annotation and sequence data

table(ann$ENTREZID == rownames(y$counts))

# Include gene names into DGE object for proper annotation 

y$genes <- ann

# Read data

head(y)
View(y)

# Filtering lowly expressed genes -  counts per milion (cpm) method - ensures removal of genes with null or negligent expression based on counts 

# obtain cpm values for all the genes in the counts data

myCPM <- cpm(countdata)
head(myCPM)
View(myCPM)

# We set a threshold cpm value of 0.5 where genes with CPM < 0.5 are removed.

thresh <- myCPM > 0.5
head(thresh)

# How many TRUE values are in each rows (Samples)

table(rowSums(thresh))

# Screening counts cpm to only keep genes which have atleast 2 TRUE in each genes

keep <- rowSums(thresh) >= 2
summary(keep)

# Plot to assess cpm score against total count for understanding sequencing depth relationship

plot(myCPM[,1], countdata[,1])

# Keeping filtered genes DGE List

y <- y[keep, keep.lib.sizes=FALSE]
View(y)

# Quality Control

## Check read size

y$samples$lib.size

## Bar plot to analyse library gene counts for each sample conditions

barplot(y$samples$lib.size,names=colnames(y),las=2)
title("Barplot of library sizes")

## Perform adjustments to the labelling on y axis

barplot(y$samples$lib.size/1e06, names=colnames(y),las=2,ann=FALSE,cex.names=0.75)
mtext(side=1, text="Samples", line=4)
mtext(side=2, text="Library size in millions",line = 3)

## Creating log scaling of the counts to ensure proper distribution

# get log counts per million

logcounts <- cpm(y,log = TRUE)

# Analyse the logcounts with boxplot

boxplot(logcounts, xlab="",ylab="Log2 scaled counts per million",las=2)

# Add median line

abline(h=median(logcounts),col="blue")
title("Boxplots of logCPMS")

# Create Multi Dimensional Scaling Plots to analyse cluster.

plotMDS(y)


# Creating the same with more interactive manner

## Align two plots side by side

par(mfrow=c(1,2))

# Type of cells

levels(sampleinfo$CellType)

# Color code 

col.cell <- c("purple","orange")[sampleinfo$CellType]
data.frame(sampleinfo$CellType,col.cell)

# MDS

plotMDS(y,col=col.cell)
legend("topleft",fill = c("purple","orange"),legend=levels(sampleinfo$CellType))

title("Cell Type")

# Status analysis

levels(sampleinfo$Status)

col.status <- c("blue","red","black")[sampleinfo$Status]

col.status
plotMDS(y,col=col.status)

legend("topleft",fill = c("blue","red","black"),legend=levels(sampleinfo$Status),cex=0.8)
title("Status")

# Corrected Sampleinfo file

sampleinfo <- read.delim("SampleInfo_Corrected.txt",stringsAsFactors = TRUE)
sampleinfo

group <- factor(paste(sampleinfo$CellType,sampleinfo$Status,sep = "."))

par(mfrow=c(1,2))
col.cell <- c("purple","orange")[sampleinfo$CellType]
col.status <- c("blue","red","black")[sampleinfo$Status]
plotMDS(y,col=col.cell)
legend("topleft",fill=c("purple","orange"),legend=levels(sampleinfo$Status))
title("Cell Type")

plotMDS(y,col=col.status)
legend("topleft",fill=c("blue","red","black"),legend=levels(sampleinfo$Status),cex=0.8)
title("Status")

# Develop more variation of MDS plots for better analysis

plotMDS(y,dim=c(3,4),col=col.status,pch=13,cex=2)
legend("topright",legend=levels(sampleinfo$Status),col = cols,pch=16)
legend("bottomright",legend = levels(sampleinfo$CellType),pch=c(1,4))


# Plots using Glimma ---- USE THIS for high quality analysis

labels <- paste(sampleinfo$SampleName,sampleinfo$CellType,sampleinfo$Status)
glMDSPlot(y,labels = labels,groups = group,folder = "mds")

# Development of hierarchy heat maps clustering based on variance genes

## we select most 500 variable genes for heat map

var_genes <- apply(logcounts,1,var)
head(var_genes)
# Obtain gene names

select_var <- names(sort(var_genes,decreasing = TRUE))[1:500]
head(select_var)

# obtain log counts matrix for genes selected

highly_variable_lcpm <- logcounts[select_var,]
dim(highly_variable_lcpm)


head(highly_variable_lcpm)

# Plot the graph

mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)

col.cell <- c("purple","orange")[sampleinfo$CellType]

heatmap.2(highly_variable_lcpm,col = rev(morecols(50)),trace = "none",main = "Top 500 variable genes across samples",ColSideColors = col.cell,scale = "row")

# Save plots

png(file="High_var_image.png")
heatmap.2(highly_variable_lcpm,col = rev(morecols(50)),trace = "none",main = "Top 500 variable genes across samples",ColSideColors = col.cell,scale = "row")
dev.off()

# Normalisation technique to deal with count composition bias

y <- calcNormFactors(y)
y$samples

# Analyse un normalised data

par(mfrow=c(1,2))
plotMD(logcounts,column = 7)
abline(h=0,col="grey")
plotMD(logcounts,column = 11)
abline(h=0,col="grey")

View(logcounts)

# Analyse normalised data

par(mfrow=c(1,2))
plotMD(y,column = 7)
abline(h=0,col="grey")
plotMD(y,column = 11)
abline(h=0,col="grey")

# SAVE PROGRESS until now to preserve the progress

save(group,y,logcounts,sampleinfo,file = "day1objects.Rdata")


# Differential Gene Expression Analysis with Limma-Voom

## load data

load("day1objects.Rdata")
objects()

# Group variable

group

# Design matrix to confirm the sample category for proper DEG process. the format is [~intercept or control + groups]. 0 as we dont have control but only various groups 

design <- model.matrix(~ 0 + group)
design

# Set column names of design matrix in proper way

colnames(design) <- levels(group)
design

# Use VOOM to transform the design data into normalised format

par(mfrow=c(1,1))
v <- voom(y,design,plot = TRUE)

View(v)

names(v)

# E column in voom is the expression value

# see boxplots to analyse expression data on categorical levels

par(mfrow=c(1,2))
boxplot(logcounts,xlab="",ylab="Log2 CPM",las=2,main="Without normalisation")
abline(h=median(logcounts),col="blue")
boxplot(v$E,xlab="",ylab="Log2 CPM",las=2,main="VOOM Normalisation")
abline(h=median(v$E),col="blue")


# Performing Differential gene expression model using linear fit model with limma

# fit the model

fit <- lmFit(v)
names(fit)

# Contrast analysis between 2 samples. null hypothesis will say the difference between 2 groups is 0 hence same and vice versa
# Differences among biological states and samples based on expression levels.
# Single sample analysis

cont.matrix <- makeContrasts(B.PregVsLac=basal.pregnant - basal.lactate,levels=design)
cont.matrix

# Contrast analysis between 2 samples. null hypothesis will say the difference between 2 groups is 0 hence same and vice versa
# Multiple sample analysis 

cont.matrix <- makeContrasts(B.PregVsLac=basal.pregnant - basal.lactate,L.PregVsLac=luminal.pregnant - luminal.lactate,levels=design)
cont.matrix


# fit the contrasts

fit.cont <- contrasts.fit(fit,cont.matrix)

# Bayes shrinkage - Basically for reduction in noises in data
fit.cont <- eBayes(fit.cont)
dim(fit.cont)

# Generate summary of DEG

summa.fit <- decideTests(fit.cont)
summary(summa.fit)
summa.fit
View(summa.fit)
# Testing out DE plots

## Highlight the significant genes

par(mfrow=c(1,2))
plotMD(fit.cont,coef = 1,status = summa.fit[,"B.PregVsLac"],values=c(-1,1),hl.col=c("blue","red"))

View(fit.cont)

volcanoplot(fit.cont,cef=1,highlight = 100,names = fit.cont$genes$SYMBOL, main="B.PregVsLac")

# Analysing individual gene and their expression levels across conditions using 1D scatter stripcharts

par(mfrow=c(1,3))

stripchart(v$E["24117",]~group)

# Modification of the above stripcharts for better visualisations to assess expression of separate genes across biological conditions
stripchart(v$E["24117",]~group,vertical=TRUE,las=2,cex.axis=0.8,pch=16,col=1:6,method="jitter")

nice.col <- brewer.pal(6,name = "Dark2")

stripchart(v$E["24117",]~group,vertical=TRUE,las=2,cex.axis=0.8,pch=16,cex=1.3,col=nice.col,method="jitter",ylab="Normalised log2 Expression",main="Wif1")

# HTML Result page 

group2 <- group
levels(group2) <- c("basal.lactate","basal.preg","basal.virgin","lum.lactate","lum.preg","lum.virgin")
glXYPlot(x=fit.cont$coefficients[,1],y=fit.cont$lods[,1],
         xlab = "LogFC",ylab="B",main="B.PregVsLac",
         counts = v$E,groups=group2,status = summa.fit[,1],
         anno = fit.cont$genes,side.main = "ENTREZID",folder="volcano")

# Testing relative to a threshold TREAT method

## Assessing statistical screening of DEGs based on log2FC values threshold- filter with logFC and recalc p values
## logFC threshold = 1

fit.treat <- treat(fit.cont,lfc=1)
res.treat <- decideTests(fit.treat)
summary(res.treat)

# see the table

gene_info <- topTable(fit.treat,coef = 1,sort.by = "p")
df_gene<-data.frame(gene_info)
View(df_gene)

# Develop volcano plots

par(mfrow=c(1,2))
plotMD(fit.treat,coef=1,status=res.treat[,"B.PregVsLac"],values=c(-1,1),hl.col=c("blue","red"))
abline(h=0,col="grey")
plotMD(fit.treat,coef=2,status=res.treat[,"L.PregVsLac"],values=c(-1,1),hl.col=c("blue","red"))
abline(h=0,col="grey")

# GO Analysis

go <- goana(fit.cont, coef = "B.PregVsLac",species = "Mm")
topGO(go,n=10)


