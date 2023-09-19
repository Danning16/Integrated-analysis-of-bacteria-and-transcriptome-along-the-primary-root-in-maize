#===============================================================================
#
#            Function prediction using PICRUST2 
#
#===============================================================================

library(phyloseq)
library(DESeq2)
library(ggplot2)
library(viridis)

# read bacteria data in phyloseq object
ps.abund = readRDS("~/emmynoether/intermediate_data/16S/pre-analysis/ps.abund.RDS")
# extract OTU ids 
ids = taxa_names(ps.abund)
write.table(ids, "~/emmynoether/intermediate_data/16S/picrust2/ps.abund.OTUid.txt",
            sep = "\n", quote = F, row.names = F, col.names = F)
# extract above OTU sequences from original fasta file
# seqtk subseq OTU_final.fasta ~/emmynoether/intermediate_data/16S/picrust2/ps.abund.OTUid.txt > ~/emmynoether/intermediate_data/16S/picrust2/ps.abund.fasta
# grep ">" ps.abund.fasta | wc -l

# !!!! The first column of this table needs to match the ids in the FASTA file. 
# $ grep '>' ps.abund.fasta > ps.ids.txt
ps.ids = read.table('./intermediate_data/16S/picrust2/ps.ids.txt')
ps.ids$V1 = gsub('>', '', ps.ids$V1)
table(ps.ids$V1 %in% ids)
# prepare data input for PICRUST2
write.table(cbind(id=ps.ids$V1, otu_table(ps.abund)[ps.ids$V1, ]), "~/emmynoether/intermediate_data/16S/picrust2/ps.abund.asvTAB.tsv",
            sep = "\t", row.names = F)
# $ conda activate picrust2
# $ picrust2_pipeline.py -s ps.abund.fasta -i ps.abund.asvTAB.tsv -o ~/emmynoether/results/16S/prediction_function/picrust2_out -p 16


# read KO table
KO.tab = read.csv("~/emmynoether/results/16S/prediction_function/picrust2_out/KO_metagenome_out/pred_metagenome_unstrat.tsv", header = T, sep = "\t")
rownames(KO.tab) = KO.tab$function.
KO.tab = KO.tab[, -1]
# remove lowly expressed KOs
KO.tab.ra = apply(KO.tab, 2, function(x){x/sum(x)})
KO.tab.ra.abund = KO.tab.ra[which(rowSums(KO.tab.ra > 0.0001) >= 0.05*ncol(KO.tab.ra)), ]
KO.tab.filter = KO.tab[rownames(KO.tab.ra.abund), ]
KO.tab.filter = round(KO.tab.filter, 0)
# normalization and PCA plot
dds = DESeqDataSetFromMatrix(countData = KO.tab.filter,
                              colData = sample_data(ps.abund),
                              design = ~ Compartment + Treatment)
# Variance Stablizing Transformation 
vsd = varianceStabilizingTransformation(dds, blind = F)  

# run PCA
pcaDat = prcomp(t(assay(vsd)))
screeplot(pcaDat, type = "l", npcs = 10, main = "Screeplot of the first 10 PCs")
percentage = round(pcaDat$sdev^2 / sum(pcaDat$sdev^2) * 100, 2)  
percentage = paste( colnames(pcaDat$x), "(", paste( as.character(percentage), "%", ")", sep=""), sep = "" )
# prepare data for PCA plot
pcs.df = as.data.frame(pcaDat$x[, 1:3])
pcs.df$Compartment = colData(vsd)$Compartment
pcs.df$Treatment = colData(vsd)$Treatment
pcs.df$Genotype = colData(vsd)$Genotype

mytheme = theme_bw() +
    theme(plot.margin = unit(c(1,1,1,1), "cm"), plot.title = element_text(size = 10),
          panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          axis.text = element_text(colour = "black", size = 8)) 

pdf(file = "./results/16S/prediction_function/predicted-KO-PCA-plot.pdf", width = 6, height = 6)
ggplot(pcs.df, aes(x=PC1,y=PC2)) +
    geom_point(aes(color=Compartment, shape = Treatment), size = 2, alpha = 0.6)  + 
    #stat_ellipse(aes(color=stage), linetype = 2,  type = "norm") +
    scale_color_brewer(palette='Dark2') +
    # treatment color
    # scale_color_manual(values = c('grey', 'orange', ' purple4'), 
    #                    limits = c("C", "N", "P")) +
    xlab(percentage[1]) + ylab(percentage[2]) + 
    coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on") +
    mytheme
dev.off()

# Read KEGG map which is downloaded from PICRUST2 wiki website
kegg_brite_map = read.table("~/emmynoether/data/OTU/picrust2_out/picrust1_KO_BRITE_map.tsv",
                        header=TRUE, sep="\t", quote = "", stringsAsFactors = FALSE, comment.char="", row.names=1)
# group KOs into KEGG pathway level 2 or 3
KO.tab.filter.L3 = categorize_by_function_l3(KO.tab.filter, kegg_brite_map)
# differential abundance analysis
md = sample_data(ps.abund)
table(colnames(KO.tab.filter.L3) == md$sampleID)
KO.tab.L3.part = KO.tab.filter.L3[, md$Compartment == 'AS']
# classify genotypes to wild type, rtcs, lateral root mutants and root hair mutants
md$Type = 'WT'
md[md$Genotype == '2', ]$Type = 'rtcs'
md[md$Genotype %in% c('3', '4'), ]$Type = 'lateralroot_mu'
md[md$Genotype %in% c('5', '6', '7'), ]$Type = 'roothair_mu'
# DESeq2 analysis
dds = DESeqDataSetFromMatrix(countData = KO.tab.L3.part,
                              colData = md[md$Compartment == 'AS', ],
                              design = ~ Type)
dds = DESeq(dds)
resultsNames(dds)
# set FDR level to 0.05
res = results(dds, contrast=c("Type", "roothair_mu", "WT"), alpha = 0.05, test="Wald")  
mcols(res, use.names = TRUE)  # show the meaning of the columns of res
summary(res)
# subset significant pathways
resSig = subset(res, padj < 0.05 & abs(log2FoldChange) > 0)
write.table(resSig, "~/emmynoether/results/16S/prediction_function/DEGs/AS/new/DEgenes_AS_lrtmu_vs_B73_log2FC_0.txt", sep = '\t', quote = F)
# prepare data for differential pathways bar plot
resSig$pathway = rownames(resSig)
resSig = data.frame(resSig)
resSig$lab = cut(resSig$padj, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), 
                  labels = c("***", "**", "*", "n.s."), right = FALSE)

pdf(file = "~/emmynoether/results/16S/prediction_function/DEGs/AS/new/Differential_pathways_between_lrtmu_vs_B73_in_PR.pdf", 
    width = 8, height = 15)
ggplot(resSig, aes(x = pathway, y = log2FoldChange)) +
    geom_bar(stat="identity", color = 'black', fill='orange', width = 0.7) + 
    geom_text(aes(y = log2FoldChange*1.2, label = lab)) + 
    coord_flip() +
    mytheme +
    ggtitle('Differential pathways between lateral root mutants and B73')
dev.off()

### The following code is downloaded from PICRUST website
### Reproducing the categorize by function (level 3) functionality in plain-text tables.
### Doing this because adding a column of KEGG Pathways to a table and then converting
### that table to BIOM is difficult.

categorize_by_function_l3 <- function(in_ko, kegg_brite_mapping) {
    # Function to create identical output as categorize_by_function.py script,
    # but with R objects instead of BIOM objects in Python.
    # Input KO table is assumed to have rownames as KOs and sample names as columns.
    
    out_pathway <- data.frame(matrix(NA, nrow=0, ncol=(ncol(in_ko) + 1)))
    
    colnames(out_pathway) <- c("pathway", colnames(in_ko))
    
    for(ko in rownames(in_ko)) {
        
        # Skip KO if not in KEGG BRITE mapping df
        # (this occurs with newer KOs that weren't present in PICRUSt1).
        if(! ko %in% rownames(kegg_brite_mapping)) {
            next
        }
        
        pathway_list <- strsplit(kegg_brite_mapping[ko, "metadata_KEGG_Pathways"], "\\|")[[1]]
        
        for(pathway in pathway_list) {
            # level 3 then choose [3], level 2 change to [2]
            pathway <- strsplit(pathway, ";")[[1]][3]
            
            new_row <- data.frame(matrix(c(NA, as.numeric(in_ko[ko,])), nrow=1, ncol=ncol(out_pathway)))
            colnames(new_row) <- colnames(out_pathway)
            new_row$pathway <- pathway
            out_pathway <- rbind(out_pathway, new_row)
        }
        
    }
    
    out_pathway = data.frame(aggregate(. ~ pathway, data = out_pathway, FUN=sum))
    
    rownames(out_pathway) <- out_pathway$pathway
    
    out_pathway <- out_pathway[, -which(colnames(out_pathway) == "pathway")]
    
    if(length(which(rowSums(out_pathway) == 0)) > 0) {
        out_pathway <- out_pathway[-which(rowSums(out_pathway) == 0), ]
    }
    
    return(out_pathway)
    
}

