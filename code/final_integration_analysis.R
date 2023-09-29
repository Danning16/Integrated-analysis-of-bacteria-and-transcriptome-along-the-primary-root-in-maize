library(WGCNA)
library(psych)
library(igraph)


options(stringsAsFactors = FALSE)

### do the following integration analysis for each compartment and for each phenotypic trait ###

#------------------------------------------------------------------------------#
#    phenotype correlation with bacterial modules                              #
#------------------------------------------------------------------------------#
# read phenotypic traits file
biomass.mean = read.table('~/emmynoether/intermediate_data/biomass_mean.txt', header = T, sep = '\t')
# read eigenOTUs of baterial WGCNA modules
bac.ME = read.table('~/emmynoether/results/16S/WGCNA/AR/ME_merged.txt', header = T)
bac.mods = 5
bac.ME$comb = paste0(bac.ME$Treatment, bac.ME$Genotype)
table(bac.ME$comb)
# calculate mean value across 4 replicates
bac.ME.mean = aggregate( . ~ comb, bac.ME[, c(2:(bac.mods+1), ncol(bac.ME))], mean) 
# Calculate Spearman correlation coefficients between module eigen-OTUs and traits
moduleTraitCor = cor(bac.ME.mean[, 2:ncol(bac.ME.mean)], 
                     (biomass.mean[match(bac.ME.mean$comb, biomass.mean$comb), 2:3]), 
                     use = "p", method = 'spearman')
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(bac.ME.mean))

# heatmap display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
rownames(textMatrix) = colnames(bac.ME.mean[, 2:ncol(bac.ME.mean)])
colnames(textMatrix) = colnames(biomass.mean[, 2:3])
pdf("~/emmynoether/results/16S/WGCNA/AR/biomass_OTU_AR_module_corr.pdf", width = 6, height = 10)
par(mar = c(15, 12, 5, 5))
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(textMatrix),
               yLabels = rownames(textMatrix),
               ySymbols = rownames(textMatrix),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("OTU-Module-Biomass relationships"))
dev.off()

#------------------------------------------------------------------------------#
#     find hub OTUs in each significant module with phenotype                  # 
#------------------------------------------------------------------------------#

# read normalized bacteria reads table
datExpr.bac = readRDS('~/emmynoether/intermediate_data/16S/wgcna/wgcna_input_AR_vst_norm.RDS')
datExpr.bac = t(as.data.frame(otu_table(datExpr.bac)))
datExpr.bac = as.data.frame(datExpr.bac)
# calculate mean value across 4 replicates
datExpr.bac$Treatment = substr(rownames(datExpr.bac), 2, 2)
datExpr.bac$Genotype = substr(rownames(datExpr.bac), 3, 3)
datExpr.bac$Compartment = substr(rownames(datExpr.bac), 4, 4)
datExpr.bac$comb = paste0(datExpr.bac$Treatment, datExpr.bac$Genotype)
datExpr.bac.mean = aggregate( . ~ comb, datExpr.bac[, c(1:ncol(datExpr.bac), ncol(datExpr.bac))], mean)
table(datExpr.bac.mean$comb)
# Spearman correlation between normalized OTU reads and phenotypic traits
otuTraitSignificance = as.data.frame(cor(datExpr.bac.mean[, 2:ncol(datExpr.bac.mean)], 
                                         biomass.mean[match(datExpr.bac.mean$comb, biomass.mean$comb), 2:3], 
                                         method = 'spearman'))
GSPvalue.otu = as.data.frame(corPvalueStudent(as.matrix(otuTraitSignificance), nrow(datExpr.bac.mean)))
colnames(otuTraitSignificance) = paste("OTUsig.", names(otuTraitSignificance), sep="")
colnames(GSPvalue.otu) = paste("p.GS.otu.", names(GSPvalue.otu), sep="")
otuSignificance= cbind.data.frame(otuTraitSignificance, GSPvalue.otu)

# calculate OTU module membership value
datKME.otu=signedKME(datExpr.bac, 
                 bac.ME[match(rownames(datExpr.bac), bac.ME$Row.names), 2:(bac.mods+1)], 
                 outputColumnName="MM.", corOptions = "method='spearman'")
# there are some OTUs with std =0, so set it to 0
datKME.otu[unique(which(is.na(datKME.otu), arr.ind = T)[,1]), ] = 0
# select all OTU modules which significantly correlate with phenotype
selected.otu.mods = gsub('ME', '', rownames(moduleTraitPvalue)[moduleTraitPvalue[,2]<0.05])
# find hub OTUs within above selected modules
hubOTUs.tot = c()
for (m in selected.otu.mods) {
    # hub OTU is defined as the OTU whose |correlation with phenotype| >0.7 and p value <0.05 and OTU module membership >0.7
    FilterOTUs= abs(otuSignificance$OTUsig.Dryweight) > 0.7 & 
        otuSignificance$p.GS.otu.Dryweight < 0.05 & 
        abs(datKME.otu[, paste0('MM.', m)]) > 0.7
    hubOTUs = colnames(datExpr.bac)[FilterOTUs]
    print(hubOTUs)
    hubOTUs.tot = c(hubOTUs.tot, hubOTUs)
}
# save to file
hubOTUs.tot = unique(hubOTUs.tot)
bac.ps = readRDS('~/emmynoether/intermediate_data/16S/pre-analysis/ps.abund.RDS')
write.table(cbind(hubOTUs.tot, tax_table(bac.ps)[hubOTUs.tot,]),
            '~/emmynoether/results/integration/AR/hubOTUs_taxa_AR.txt', 
            sep = '\t', quote = F, row.names = F)

#------------------------------------------------------------------------------#
#    phenotype correlation with gene modules                                   #
#------------------------------------------------------------------------------#
# read gene module eigen-gene file
gene.ME = read.table('~/emmynoether/results/rna-seq/WGCNA/AR/merged_MEs.txt', header = T, sep = '\t')
rownames(gene.ME) = gene.ME$Row.names
gene.ME = gene.ME[, -1]
# this needs to be changed accordingly
gene.mods = 25  
# calculate mean value across 4 replicates
gene.ME$comb = paste0(gene.ME$Treatment, gene.ME$Genotype)
gene.ME.mean = aggregate( . ~ comb, gene.ME[, c(1:gene.mods, ncol(gene.ME))], mean) 
# !!!! sample names should be consistent in eigen genes and traits 
# Calculate Spearman correlation coefficients between module eigen-genes and traits
moduleTraitCor = cor(gene.ME.mean[, 2:ncol(gene.ME.mean)], 
                     biomass.mean[match(gene.ME.mean$comb, biomass.mean$comb), 2:3], 
                     use = "p", method = 'spearman')
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(gene.ME.mean))
# heatmap display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
rownames(textMatrix) = colnames(gene.ME.mean[,2:ncol(gene.ME.mean)])
colnames(textMatrix) = colnames(biomass.mean[, 2:3])

pdf("~/emmynoether/results/rna-seq/WGCNA/AR/biomass_gene_AR_module_corr.pdf", width = 6, height = 10)
par(mar = c(15, 12, 5, 5))
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(biomass.mean[, 2:3]),
               yLabels = colnames(gene.ME.mean[,2:ncol(gene.ME.mean)]),
               ySymbols = colnames(gene.ME.mean[,2:ncol(gene.ME.mean)]),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Gene-AR-Module-Biomass relationships"))
dev.off()

#------------------------------------------------------------------------------#
#     find hub genes in each significant module with phenotype                 # 
#------------------------------------------------------------------------------#
# read normalized gene reads table
datExpr.rna = readRDS('~/emmynoether/intermediate_data/rna-seq/WGCNA/normalized_gene(20065)_counts_vst_AB_WGCNA_input.RDS')
sample.md = read.csv('~/emmynoether/intermediate_data/rna-seq/sample_metadata.csv')
datExpr.rna1 = merge(datExpr.rna, sample.md, by.x = 'row.names', by.y = 'sampleID')
# calculate mean gene reads across 4 replicates
datExpr.rna1$comb = paste0(datExpr.rna1$Treatment, datExpr.rna1$Genotype)
datExpr.rna.mean = aggregate( . ~ comb, datExpr.rna1[, c(2:(ncol(datExpr.rna)+1), ncol(datExpr.rna1))], mean)
# calculate gene-trait correlation and p value
geneTraitSignificance = as.data.frame(cor(datExpr.rna.mean[, 2:ncol(datExpr.rna.mean)], 
                                          biomass.mean[match(datExpr.rna.mean$comb, biomass.mean$comb), 2:3], 
                                          method = 'spearman'))
GSPvalue.gene = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nrow(datExpr.rna.mean)))
colnames(geneTraitSignificance) = paste("Genesig.", names(geneTraitSignificance), sep="")
colnames(GSPvalue.gene) = paste("p.GS.gene.", names(GSPvalue.gene), sep="")
geneSignificance= cbind.data.frame(geneTraitSignificance, GSPvalue.gene)
# calculate gene module membership
datKME=signedKME(datExpr.rna, gene.ME[, 1:(gene.mods)], 
                 outputColumnName="MM.", corOptions = "method='spearman'")
# select all gene modules which significantly correlate with phenotype
selected.mods = gsub('ME', '', rownames(moduleTraitPvalue[moduleTraitPvalue[,2] <0.05, ]))
# find hub genes within each selected module
hubGenes.tot = c()
for (m in selected.mods) {
    Filtergenes= abs(geneSignificance$Genesig.Dryweight) > 0.7 & 
        geneSignificance$p.GS.gene.Dryweight < 0.05 & 
        abs(datKME[, paste0('MM.', m)]) > 0.7
    hubGenes = dimnames(data.frame(datExpr.rna))[[2]][Filtergenes]
    hubGenes.tot = c(hubGenes.tot, hubGenes)
}
hubGenes.tot = unique(hubGenes.tot)

#------------------------------------------------------------------------------#
#   Spearman correlation betweeen hub genes and hub OTUs                       #
#------------------------------------------------------------------------------#
# Spearman correlation
cor.res.hub = corr.test(datExpr.rna[, hubGenes.tot], 
                    datExpr.bac[match(substr(rownames(datExpr.rna), 1,4), gsub('R', '', rownames(datExpr.bac))), hubOTUs.tot], 
                    method = "spearman", adjust = "fdr", ci=F, normal = F)
cor.mat = as.data.frame(as.table(cor.res.hub$r))
cor.padj = as.data.frame(as.table(cor.res.hub$p.adj))
corr.results.hub = cbind.data.frame(cor.mat, padj=cor.padj$Freq)
colnames(corr.results.hub)[1:3] = c("row", "col", "rho")
head(corr.results.hub)
range(corr.results.hub$rho)
# save significant results to file
corr.res.sig.hub = corr.results.hub[abs(corr.results.hub$rho) >0.7 & corr.results.hub$padj <0.05, ]
write.table(corr.results.hub, "~/emmynoether/results/integration/AR/corr_results_hubgene_hubOTU_dryweight_AR.txt", sep = "\t", row.names = F)

#------------------------------------------------------------------------------#
#       create edges list and nodes list  to build networks in Cytoscape       #
#------------------------------------------------------------------------------#
# !! the column 'row' is factor, so must change to vector, otherwise it is in factor level order, not the list order
# get gene-trait correlation
corr.gene.dryWT = cbind(row = as.vector(unique(corr.res.sig.hub$row)), col = 'DryWT',
                        geneSignificance[as.vector(unique(corr.res.sig.hub$row)), c(2, 4)])
colnames(corr.gene.dryWT)[3:4] = c('rho', 'padj')
# get OTU-trait correlation
corr.OTU.dryWT = cbind(row = as.vector(unique(corr.res.sig.hub$col)), col = 'DryWT',
                        otuSignificance[as.vector(unique(corr.res.sig.hub$col)), c(2, 4)])
colnames(corr.OTU.dryWT)[3:4] = c('rho', 'padj')
# combine gene-OTU correlation with above gene-trait and OTU-trait correlations
edges.list = rbind(corr.res.sig.hub, corr.gene.dryWT, corr.OTU.dryWT)
write.table(edges.list, "~/emmynoether/results/integration/AR/AR_net_edges_list.txt", sep = "\t", row.names = F, quote = F)
# create nodes list
AR.net = graph_from_data_frame(edges.list, directed = F)
nodes.AR.dryWT.net = data.frame(ID=V(AR.net)$name, 
                          hubScore=hub_score(AR.net)$vector, 
                          Degree=igraph::degree(AR.net))
nodes.AR.dryWT.net$type = substr(nodes.AR.dryWT.net$ID, 1, 3)
mean_vst_Expr = c(apply(datExpr.rna[, nodes.AR.dryWT.net[nodes.AR.dryWT.net$type == 'Zm0', 'ID']], 2, mean), 
                  apply(datExpr.bac[, nodes.AR.dryWT.net[nodes.AR.dryWT.net$type == 'OTU', 'ID']], 2, mean), 
                  DryWT=1)
nodes.AR.dryWT.net$mean_vst = mean_vst_Expr[match(nodes.AR.dryWT.net$ID, names(mean_vst_Expr))]
write.csv(nodes.AR.dryWT.net, '~/emmynoether/results/integration/AR/AR_net_nodelist.csv', quote = F, row.names = F)



