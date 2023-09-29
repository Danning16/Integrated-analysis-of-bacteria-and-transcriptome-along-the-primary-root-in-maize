library(Rsubread)
library(DESeq2)
library(pheatmap)
library(apeglm)
library(RColorBrewer)
library(ggplot2)
library(tidyverse)
library(limma)
library(edgeR)


#===============================================================================
#
#              RNA-seq data analysis
#
#===============================================================================

# read all bam files
all_files = list.files(pattern = ".bam", path = "/home/data/Precision-5820_4TB/emmynoether/rna/bam_files")
# generate gene reads counts table
setwd("/home/data/Precision-5820_4TB/emmynoether/rna/bam_files")
fc = featureCounts(all_files, annot.ext = "~/emmynoether/data/zm_ref_v5/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.53.chr.gtf", 
                    isGTFAnnotationFile = TRUE, isPairedEnd = TRUE, 
                    countMultiMappingReads = FALSE, useMetaFeatures = TRUE,
                    countChimericFragments = FALSE, nthreads = 18)
setwd("~/emmynoether")
write.csv(fc$counts, file = "./intermediate_data/rna-seq/raw_counts_table_ref_v5.csv", row.names = TRUE)
write.csv(fc$stat, file = "./intermediate_data/rna-seq/featurecounts_stat.csv", row.names = F)

# ------------------------------------------------------------------------------
#    subset expressed genes for each compartment
# ------------------------------------------------------------------------------

raw.counts = fc$counts
colnames(raw.counts) = substr(colnames(raw.counts), 1, 7)
# Obtain CPMs
#myCPM = cpm(raw.counts)
# Which values are greater than 10
thresh = raw.counts >= 10
table(rowSums(thresh))
# we would like to keep genes that have at least 4 TRUES in each row of thresh
keep = rowSums(thresh) >= 4
# Subset the rows of counts table to keep the more highly expressed genes
counts.keep = raw.counts[keep,]
dim(counts.keep)   # 29335  x 204
write.csv(counts.keep, file = "./intermediate_data/rna-seq/filtered_counts_table_ref_v5.csv", row.names = TRUE)
# read samples meta data file
sample.md = read.csv("~/emmynoether/intermediate_data/rna-seq/sample_metadata.csv", header = T)
table(sample.md$sampleID == colnames(counts.keep))
# subset expressed genes for each compartment
AR.genes = counts.keep[, sample.md[sample.md$Compartment == 'A', 1]]
AR.genes = AR.genes[rowSums(AR.genes) > 0, ]
write.table(rownames(AR.genes), '~/emmynoether/intermediate_data/rna-seq/expressed_genes/PR_genes.txt', sep = '\n', quote = F, row.names = F, col.names = F)

BR.genes = counts.keep[, sample.md[sample.md$Compartment == 'B', 1]]
BR.genes = BR.genes[rowSums(BR.genes) > 0, ]
write.table(rownames(BR.genes), '~/emmynoether/intermediate_data/rna-seq/expressed_genes/LR_genes.txt', sep = '\n', quote = F, row.names = F, col.names = F)

CR.genes = counts.keep[, sample.md[sample.md$Compartment == 'C', 1]]
CR.genes = CR.genes[rowSums(CR.genes) > 0, ]
write.table(rownames(CR.genes), '~/emmynoether/intermediate_data/rna-seq/expressed_genes/Cortex_genes.txt', sep = '\n', quote = F, row.names = F, col.names = F)

SR.genes = counts.keep[, sample.md[sample.md$Compartment == 'S', 1]]
SR.genes = SR.genes[rowSums(SR.genes) > 0, ]
write.table(rownames(SR.genes), '~/emmynoether/intermediate_data/rna-seq/expressed_genes/Stele_genes.txt', sep = '\n', quote = F, row.names = F, col.names = F)

#-------------------------------------------------------------------------------
#             Principle Component Analysis
#-------------------------------------------------------------------------------

dds = DESeqDataSetFromMatrix(countData = counts.keep,
                              colData = sample.md,
                              design = ~ Compartment+Treatment+Genotype)
# Variance Stablizing Transformation 
vsd = varianceStabilizingTransformation(dds, blind = F)  
# run PCA
pcaDat = prcomp(t(assay(vsd)))
screeplot(pcaDat, type = "l", npcs = 10, main = "Screeplot of the first 10 PCs")
percentage = round(pcaDat$sdev^2 / sum(pcaDat$sdev^2) * 100, 2)  
percentage = paste( colnames(pcaDat$x), "(", paste( as.character(percentage), "%", ")", sep=""), sep = "" )
# prepare data for PCA plot
pcs.df = as.data.frame(pcaDat$x[, 1:5])
table(colnames(assay(vsd)) == colData(vsd)$sampleID)
pcs.df$Treatment = colData(vsd)$Treatment
pcs.df$Compartment = colData(vsd)$Compartment
pcs.df$Genotype = colData(vsd)$Genotype

mytheme = theme_bw() +
  theme(plot.margin = unit(c(1,1,1,1), "cm"), plot.title = element_text(size = 10),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black", size = 8)) 

pdf(file = "./gene-PCA-Compartment.pdf", width = 6, height = 6)
ggplot(pcs.df,aes(x=PC1,y=PC2)) +
  geom_point(aes(color=Compartment), size = 3)  + 
  #scale_shape_manual(values=c(3, 16))+
  #facet_wrap(~Compartment) +
  stat_ellipse(aes(color=Compartment), linetype = 2,  type = "norm") +
    # compartment color
    scale_color_manual(values = c('seagreen', 'limegreen', 'royalblue', 'cornflowerblue'),
                       limits = c("A", "B", "C", "S")) +
    # genotype color
    #scale_color_brewer(palette='Dark2') +
    # treatment color
    # scale_color_manual(values = c('grey', 'orange', ' purple4'), 
    #                    limits = c("C", "N", "P")) +
    scale_x_continuous(limits=c(-140, 200), breaks = seq(-100, 150, 50)) +
    scale_y_continuous(limits=c(-110, 110), breaks = seq(-100, 100, 50)) +
  xlab(percentage[1]) + ylab(percentage[2]) + 
  coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on") +
  mytheme
dev.off()

#-------------------------------------------------------------------------------
#                  PERMANOVA test
#-------------------------------------------------------------------------------

library(vegan)

vsd = readRDS("~/emmynoether/intermediate_data/rna-seq/normalized_gene(29335)_counts_vst.RDS")
# Extract coldata
coldata1 = colData(vsd) 
# Extract normalized count data
countdata1 = assay(vsd)
# !! change genotype to factor !!
sample.md$Genotype = factor(sample.md$Genotype)
all.dist = dist(t(countdata1)) 
# PERMANOVA test marginal effects of each factor
adonis2(all.dist~Compartment+Genotype+Treatment, 
       data = sample.md, 
       permutations = 3999, parallel = 16, by="margin")

# PERMANOVA test for treatment and genotype under each compartment 
permanova.res = data.frame()
for (comp in c("A", "B", "C", "S")) {
    print(comp)
    countdata = countdata1[, coldata1$Compartment == comp]
    dist.comp = dist(t(countdata))
    perm1 = adonis2(dist.comp~Genotype+Treatment, 
            data = sample.md[sample.md$Compartment == comp, ], 
            permutations = 3999, parallel = 16, by="margin")
    print(perm1)
    permanova.res = rbind(permanova.res, cbind(data.frame(Compartment=comp), perm1[1:2, ], type = rownames(perm1)[1:2]))
    
  }
 
write.table(permanova.res, "~/emmynoether/results/rna-seq/permanova/permanova_results_4_compartments.txt", 
            row.names = F, sep = "\t", quote = F)

#  PERMANOVA results plot 
pdf('./rna-permanova-res-barplot.pdf', width = 6, height = 6)
ggplot(data = permanova.res, aes(x=type, y = R2)) +
    geom_bar(aes(fill = Compartment),stat="identity", position = position_dodge(),
             alpha = 0.8, color='black') +
    # compartment color
    scale_fill_manual(values = c('seagreen', 'limegreen', 'royalblue', 'cornflowerblue'),
                       limits = c("A", "B", "C", "S")) +
    ylim(c(0,1))+
    # genotype color
    #scale_color_brewer(palette='Dark2') +
    # treatment color
    # scale_color_manual(values = c('grey', 'orange', ' purple4'), 
    #                    limits = c("C", "N", "P")) +
    mytheme
dev.off()

#-------------------------------------------------------------------------------
#           Differential expression analysis
#-------------------------------------------------------------------------------

setwd("~/emmynoether/results/rna-seq/DEG/")

sample.md$Type = 'WT'
sample.md[sample.md$Genotype == '2', ]$Type = 'WT'
sample.md[sample.md$Genotype %in% c('3', '4'), ]$Type = 'lateralroot_mu'
sample.md[sample.md$Genotype %in% c('5', '6', '7'), ]$Type = 'roothair_mu'

# do the following for each compartment
counts.geno = counts.keep[, sample.md[sample.md$Compartment=='A', ]$sampleID]
print(dim(counts.geno))
# remove lowly expressed genes
thresh = counts.geno > 10
# we would like to keep genes that have at least 4 TRUES in each row of thresh
keep = rowSums(thresh) >= 4
# Subset the rows of count data to keep the more highly expressed genes
counts.geno1 = counts.geno[keep,]
print(dim(counts.geno1))
dds = DESeqDataSetFromMatrix(countData = counts.geno1,
                              colData = sample.md[sample.md$Compartment=='A', ],
                              design = ~ Type)
# DESeq differential expession analysis
dds = DESeq(dds)
print(resultsNames(dds))
for(g in resultsNames(dds)[-1]){
      print("......genotype.......")
      print(g)
      res = results(dds, name = g, alpha = 0.05, test="Wald")
      print(summary(res))
      # subset significant genes
      res.Sig = subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
      
     
      if(dim(res.Sig)[1] == 0){
        next
      }
      else{
          write.table(res.Sig, paste0(paste0("./DEgenes_AR_", g), '.txt'), sep = '\t', quote = F)
      }

}

    



sessionInfo()         # show the version of packages


