library(WGCNA)
library(DESeq2)
library(tidyverse)
library(dplyr)


options(stringsAsFactors = FALSE)
enableWGCNAThreads()

#===============================================================================
#
#  Read bacteria data and plot the sample tree
#
#===============================================================================

bac.ps = readRDS("~/emmynoether/intermediate_data/16S/pre-analysis/ps.abund.RDS")
bac.comp = subset_samples(bac.ps, Compartment=='AR') 
# keep taxa with relative abundance > 0 at least in 10 samples
bac.RA = transform_sample_counts(bac.comp, function(x) x/sum(x))
bac.RA.abund = filter_taxa(bac.RA, 
                           function(x) sum(x > 0) >= 10, TRUE)  
bac.comp.abund = subset_taxa(bac.comp, taxa_names(bac.comp)%in%taxa_names(bac.RA.abund))

# DESeq2 Normalization 
bac.deseq = phyloseq_to_deseq2(bac.comp.abund, ~Genotype+Treatment)
bac.deseq = estimateSizeFactors(bac.deseq, type="poscount")
bac.vst = varianceStabilizingTransformation(bac.deseq, blind=F, fitType = 'local')   
bac.vst.norm.comp = bac.comp.abund
otu_table(bac.vst.norm.comp) = otu_table(assay(bac.vst), taxa_are_rows = T)
# cluster samples
datExpr.bac = t(as.data.frame(otu_table(bac.vst.norm.comp)))
sampleTree = hclust(dist(datExpr.bac), method = "average")
# plot sample tree
pdf(file = "./1-n-sampleClustering.pdf", width = 50, height = 9)
par(cex = 0.8)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="",
     cex.lab = 1.5,cex.axis = 1.5, cex.main = 2)
dev.off()

#===============================================================================
#
#       Choose soft threshold parameter
#
#===============================================================================

# Choose a set of soft threshold parameters
powers = c(c(1:10), seq(from = 12, to=30, by=2))
sft = pickSoftThreshold(datExpr.bac, powerVector = powers, verbose = 5, networkType = 'unsigned', 
                        corOptions = list(method='spearman')) 
# Scale-free topology fit index as a function of the soft-thresholding power
pdf(file = "./2-n-sft.pdf", width = 9, height = 5)
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
#     Turn data expression into topological overlap matrix
#
#===============================================================================

# Turn data expression into topological overlap matrix
power=sft$powerEstimate   # or selected based on above figures
TOM = TOMsimilarityFromExpr(datExpr.bac, power = power, TOMType="unsigned")
saveRDS(TOM, file = "./TOM_power5.RDS")
dissTOM = 1-TOM
# Plot gene tree
geneTree = hclust(as.dist(dissTOM), method = "average")
pdf(file = "./WGCNA/3-gene_cluster.pdf", width = 12, height = 9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
dev.off()

#===============================================================================
#
#      Construct modules
#
#===============================================================================

# Module identification using dynamic tree cut
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 4, # The higher the value, the more smaller clusters will be produced. 
                            pamRespectsDendro = FALSE, minClusterSize = 3)
table(dynamicMods)
length(table(dynamicMods)) 
# Convert numeric labels into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

#===============================================================================
#
#      Merge modules
#
#===============================================================================

# Merge close modules
MEDissThres=0.25
merge = mergeCloseModules(datExpr.bac, dynamicColors, cutHeight = MEDissThres, verbose = 3) 
mergedColors.otu = merge$colors  
table(mergedColors.otu)
# Plot merged module tree
pdf(file = "./WGCNA/5-merged_Module_Tree.pdf", width = 20, height = 9)  
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), 
                    c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels = FALSE, 
                    hang = 0.03, addGuide = TRUE, guideHang = 0.05)  
dev.off()

# Heatmap of new module eigen-OTUs and samples
pdf(file="~/emmynoether/results/16S/WGCNA/AR/newMEs.pdf",heigh=60,width=20)
row.names(merge$newMEs) = rownames(datExpr.bac)
pheatmap(merge$newMEs,cluster_col=T,cluster_row=T, breaks = c(seq(-0.5, 0.5, 0.01)),
         show_rownames=T,show_colnames=T,fontsize=20)
dev.off()
# save new eigen-OTUs to text file for integration analysis
write.table(merge(merge$newMEs, sample_data(bac.vst.norm.comp), by="row.names"), 
            "~/emmynoether/results/16S/WGCNA/AR/ME_merged.txt", 
            sep = "\t", quote = F, row.names = F)




