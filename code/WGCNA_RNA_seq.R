library(WGCNA)
library(DESeq2)
library(tidyverse)
library(dplyr)


options(stringsAsFactors = FALSE)
enableWGCNAThreads()

#===============================================================================
#
#  Read the gene counts table and plot the sample tree
#
#===============================================================================

# read sample meta data
sample.md = read.csv("~/emmynoether/intermediate_data/rna-seq/sample_metadata.csv", header = T)
sample.md$Genotype = as.factor(sample.md$Genotype)
# read vst normalized RNA data
rna.data = read.csv("~/emmynoether/intermediate_data/rna-seq/filtered_counts_table_ref_v5.csv")
rownames(rna.data) = rna.data$X
rna.data = rna.data[, -1]
# filter samples for each compartment
rna.data.comp = rna.data[, sample.md[sample.md$Compartment == 'S', ]$sampleID]
# Filtering to remove lowly expressed genes 
thresh = rna.data.comp > 5
# we would like to keep genes that have at least 10 TRUES in each row of thresh
keep = rowSums(thresh) >= 10
# Subset the rows of count data to keep the more highly expressed genes
counts.keep = rna.data.comp[keep,]
dim(counts.keep)  
# prepare normalized gene expression for WGCNA analysis
dds = DESeqDataSetFromMatrix(countData = counts.keep,
                              colData = sample.md[sample.md$Compartment == "S", ],
                              design = ~ Genotype+Treatment)
# Variance Stablizing Transformation 
vsd = varianceStabilizingTransformation(dds, blind = F)  
saveRDS(vsd, '~/emmynoether/results/rna-seq/WGCNA/SR/normalized_gene_counts_vst_SR_WGCNA_input.RDS')
# cluster samples
datExpr.rna = t(assay(vsd))
sampleTree = hclust(dist(datExpr.rna), method = "average")
# plot sample tree
pdf(file = "./results/rna-seq/WGCNA/filtered/1-n-sampleClustering.pdf", width = 20, height = 9)
par(cex = 1.2)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="",
     cex.lab = 1.5,cex.axis = 1.5, cex.main = 2)
dev.off()

#===============================================================================
#
#    Choose soft threshold parameter
#
#===============================================================================

# Choose a set of soft threshold parameters
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr.rna, powerVector = powers, verbose = 5) 
# Scale-free topology fit index as a function of the soft-thresholding power
pdf(file = "./results/rna-seq/WGCNA/SR/2-n-sft.pdf", width = 9, height = 5)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
	xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
	main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
labels=powers,cex=cex1,col="red")
# this line corresponds to using an R^2 cut-off of h
abline(h=0.8,col="red") 
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
	xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
	main = paste("Mean connectivity")) 
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

#===============================================================================
#
#     Turn gene expression into topological overlap matrix
#
#===============================================================================

# Turn gene expression into topological overlap matrix
power=sft$powerEstimate
TOM.rna = TOMsimilarityFromExpr(datExpr.rna, power = power, TOMType="unsigned")
saveRDS(TOM.rna, file = "~/emmynoether/results/rna-seq/WGCNA/SR/TOM_power10.RDS")
dissTOM = 1-TOM.rna
# Plot gene tree
geneTree = hclust(as.dist(dissTOM), method = "average")
pdf(file = "./results/rna-seq/WGCNA/3-gene_cluster.pdf", width = 12, height = 9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
dev.off()

#===============================================================================
#
#   Construct modules
#
#===============================================================================

# Module identification using dynamic tree cut
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2, 
                            pamRespectsDendro = FALSE, minClusterSize = 30)
table(dynamicMods)
length(table(dynamicMods)) 
# Convert numeric labels into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
pdf(file = "4-module_tree.pdf", width = 8, height = 6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",dendroLabels = FALSE,
                    hang = 0.03,addGuide = TRUE, guideHang = 0.05,main = "Gene dendrogram and module colors")
dev.off()

#===============================================================================
#
#   Merge close modules
#
#===============================================================================

# Merge close modules
MEDissThres=0.25
merge = mergeCloseModules(datExpr.rna, dynamicColors, cutHeight = MEDissThres, verbose = 3) 
mergedColors.rna = merge$colors 
names(mergedColors.rna) = colnames(datExpr.rna)
table(mergedColors.rna)

# Plot merged module tree
pdf(file = "./results/rna-seq/WGCNA/5-merged_Module_Tree.pdf", width = 20, height = 9)  
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors.rna), 
                    c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels = FALSE, 
                    hang = 0.03, addGuide = TRUE, guideHang = 0.05)  
dev.off()

write.table(merge(merge$newMEs, sample.md, by.x="row.names", by.y = "sampleID"), 
            "~/emmynoether/results/rna-seq/WGCNA/SR/merged_MEs.txt", 
            sep = "\t", quote = F, row.names = F)



