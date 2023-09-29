library(phyloseq)
library(biomformat)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(DESeq2)
library(vegan)
library(FSA)
library(lme4)
library(rcompanion)
library(viridis)


mytheme = theme_bw() +
  theme(plot.margin = unit(c(1,1,1,1), "cm"), plot.title = element_text(size = 10),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black", size = 8)) 


setwd("~/emmynoether/")

#===============================================================================
#
#      import to phyloseq
#
#===============================================================================

# read OTU table in biom format
biom_table = read_biom("./data/OTU/OTU_table.biom")
feature.table = as.matrix(biom_data(biom_table))

# read taxonomy file
taxonomy = read.table("./data/OTU/OTU_tax_assignments.txt", sep = "\t", header = F)
colnames(taxonomy) = c("OTU_ID", "Taxon", "Confidence")
taxa.table = taxonomy %>%
  separate(Taxon, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";") 
table(rownames(feature.table) %in% taxa.table$OTU_ID)
rownames(taxa.table) = taxa.table$OTU_ID
taxa.table = taxa.table[, -1]
taxa.table = as.matrix(taxa.table)

# read metadata file
md = read.table("./intermediate_data/16S/metadata_OTU_final.txt", header = T)

# read phylogetic tree 
tree = read_tree("./data/OTU/OTU_final_phylogeny_tree.txt")

# create phyloseq object
OTU.tab = otu_table(feature.table, taxa_are_rows = TRUE)
TAX.tab = tax_table(taxa.table[rownames(feature.table), ])
rownames(md) = md$sampleID
md = md[colnames(OTU.tab), ]
sample.tab = sample_data(md)
ps = phyloseq(OTU.tab, TAX.tab, sample.tab, tree)

# keep only bacteria OTUs
ps.tree.prefiltered = subset_taxa(ps, Phylum!='Not_Available' & Kingdom=="Bacteria")

# check if every sample has >1000 reads
table(sample_sums(ps.tree.prefiltered) >5000 )    
range(rowSums(otu_table(ps.tree.prefiltered)) )

# only keep taxa with RA > 0.1% in at least one sample
ps.tree.prefiltered.ra = transform_sample_counts(ps.tree.prefiltered, function(x){x/sum(x)})
ps.tree.filtered.ra.abund = filter_taxa(ps.tree.prefiltered.ra, function(x) sum(x > 0.001) >= 1, TRUE)  
ps.tree.filtered.abund = subset_taxa(ps.tree.prefiltered, taxa_names(ps.tree.prefiltered)%in%taxa_names(ps.tree.filtered.ra.abund))

#===============================================================================
#
#     overlap OTUs between each compartment  
#
#===============================================================================

# filter OTUs for different compartments
ps.abund.comp = subset_samples(ps.tree.filtered.abund, Compartment=="SR")
ps.abund.comp.ra = transform_sample_counts(ps.abund.comp, function(x){x/sum(x)})
ps.abund.comp.ra.abund = filter_taxa(ps.abund.comp.ra, function(x) sum(x > 0.001) >= 0.05*nsamples(ps.abund.comp.ra), TRUE)  
ps.abund.comp.abund = subset_taxa(ps.abund.comp, taxa_names(ps.abund.comp)%in%taxa_names(ps.abund.comp.ra.abund))
# do above filtering for each compartment
ps.AS = ps.abund.comp.abund
ps.BS = ps.abund.comp.abund
ps.AR = ps.abund.comp.abund
ps.BR = ps.abund.comp.abund
ps.CR = ps.abund.comp.abund
ps.SR = ps.abund.comp.abund

write.table(taxa_names(ps.SR), './intermediate_data/16S/pre-analysis/OTU_names/SR.OTUs.names.txt', 
            sep = '\n', quote = F, row.names = F, col.names = F)

#===============================================================================
#
#       alpha diversity
#
#===============================================================================

pdf("./results/16S/alpha_div/rarecurve.pdf", height = 6, width = 6)
rarecurve(t(as.data.frame(otu_table(ps.tree.filtered.abund))), step=100, label = F)
dev.off()
# the minimum sample reads is 245(our cortex and stele samples have less reads), so choose 200 as rarefaction level
# rarefaction
ps.rare = rarefy_even_depth(ps.tree.filtered.abund, 200, rngseed = 1625458, replace = FALSE)

# generate a data.frame with alpha diversity measures
adiv = data.frame(
  "Observed" = phyloseq::estimate_richness(ps.rare, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(ps.rare, measures = "Shannon"))
adiv = cbind(adiv, sample_data(ps.rare))
write.table(adiv, file = "./results/16S/alpha_div/alpha_diversity_table_200.txt", 
            sep = "\t", quote = F, row.names = F)

## BOX plot of alpha diversity
setwd("./results/16S/alpha_div/")

# do Dunn test 
dunn.res = dunnTest(Shannon~Compartment, data = adiv, method = "bh")
dunn.res$res$Compartment = gsub("_.*$", "", dunn.res$res$Comparison)
label.dunn = cldList(P.adj ~ Comparison, data = dunn.res$res, threshold = 0.05)[, 1:2]

temp = ggplot(data= adiv, aes(x = Compartment, y = Shannon, color=Compartment)) +
    geom_boxplot(outlier.colour = NA) +
    geom_point(position=position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75), 
               aes(color=Compartment), alpha=0.6, size=1, stroke=0) +
    scale_color_manual(values = c('chocolate4', 'orange', 'chocolate1', 'seagreen', 
                                  'limegreen', 'royalblue', 'cornflowerblue'), 
                        limits = c("SoilS", "AS", "BS", "AR", "BR", "CR", "SR")) +
    geom_text(data = label.dunn, aes(x = Group, y=rep(6, nrow(label.dunn)), 
                                     label=Letter), 
              position = position_dodge(width = 1), inherit.aes = F) +
    scale_y_continuous(limits=c(-0.5, 6.5), breaks = seq(0, 6, 2)) +
    mytheme

ggsave(temp, file = "bac_alpha_div_rarefy_200.pdf", width = 15, height = 9, units = "cm")


#===============================================================================
#
#              beta diversity
#
#===============================================================================

# DESeq2 Normalization 
bac.deseq = phyloseq_to_deseq2(ps.tree.filtered.abund, ~Compartment+Treatment)
bac.deseq = estimateSizeFactors(bac.deseq, type="poscount")
bac.vst = varianceStabilizingTransformation(bac.deseq, blind=F)   
bac.vst.norm = ps.tree.filtered.abund
otu_table(bac.vst.norm) = otu_table(assay(bac.vst), taxa_are_rows = T)

# Round negative values up to zeroes, to enable Bray-Curtis calculations
bac.vst.norm.comp = transformSampleCounts(bac.vst.norm.comp,function(x) ifelse(x<0,0,x)) 
saveRDS(bac.vst.norm.comp, "./intermediate_data/16S/beta_div/bac.vst.norm.RDS")

# PERMANOVA test 
bac.perm = data.frame()
# do the following for each compartment
bac.vst.norm.comp = subset_samples(bac.vst.norm, Compartment == 'BS')
bac.vst.norm.comp = subset_taxa(bac.vst.norm.comp, taxa_sums(bac.vst.norm.comp) > 0)
bray.dist = phyloseq::distance(bac.vst.norm.comp, method = "bray") 
perm.res = adonis2(bray.dist~Genotype+Treatment, 
        data = as(sample_data(bac.vst.norm.comp), "data.frame"), 
        permutations = 3999, parallel = 16, by="margin")
# save results in data frame 
bac.perm = rbind(bac.perm, cbind(data.frame(Compartment='BS'), perm.res[1:2, ], type = rownames(perm.res)[1:2]))
write.table(bac.perm, "./results/16S/PERMANOVA/permanova_results_7_compartments.txt", 
            row.names = F, sep = "\t", quote = F)

#  PERMANOVA results plot 
bac.perm$Compartment = factor(bac.perm$Compartment, 
                              levels = c("SoilS", "AS", "BS", "AR", "BR", "CR", "SR"))
pdf('./bac-permanova-res-barplot.pdf', width = 6, height = 6)
ggplot(data = bac.perm, aes(x=type, y = R2)) +
    geom_bar(aes(fill = Compartment),stat="identity", position = position_dodge(),
             alpha = 0.8, color='black') +
    # compartment color
    scale_fill_manual(values = c('chocolate4', 'orange', 'chocolate1', 'seagreen', 
                                 'limegreen', 'royalblue', 'cornflowerblue'), 
                       limits = c("SoilS", "AS", "BS", "AR", "BR", "CR", "SR")) +
    ylim(c(0,1))+
    mytheme
dev.off()

# PCoA plot
# do the following for each compartment
bac.vst.norm.comp = subset_samples(bac.vst.norm, Compartment == 'AS')
bac.vst.norm.comp = subset_taxa(bac.vst.norm.comp, taxa_sums(bac.vst.norm.comp) > 0)
# PCoA ordination analysis
PCoA.bray.vst = ordinate(bac.vst.norm.comp, method = "PCoA", distance = "bray")
pcoa.df = cbind.data.frame(PCoA.bray.vst$vectors[, 1:5], sample_data(bac.vst.norm.comp))
pct = PCoA.bray.vst$values$Relative_eig[1:2]
percentage = paste( c("PCoA1", "PCoA2"), "(", paste( as.character(round(pct*100,2)), "%", ")", sep=""), sep = "" )
pcoa.df$Type = 'WT'
pcoa.df[pcoa.df$Genotype == '2', ]$Type = 'rtcs'
pcoa.df[pcoa.df$Genotype %in% c('3', '4'), ]$Type = 'lateralroot_mu'
pcoa.df[pcoa.df$Genotype %in% c('5', '6', '7'), ]$Type = 'roothair_mu'
# PCoA plot
pcoa1.fig = ggplot(data= pcoa.df, aes(x = Axis.1, y = Axis.2, color = Type)) +
    geom_point(alpha=0.8, size=2) +
    #facet_grid(cols = vars(Treatment)) +
    #stat_ellipse(linetype = 2) +
    scale_color_manual(values = c('brown', 'royalblue', 'black')) +
    ggtitle("PCoA plot for primary root rhizosphere") + 
    mytheme +
    xlab(percentage[1]) + ylab(percentage[2]) +
    coord_fixed() 
ggsave(pcoa.fig, file = "./PCOA-rhizosphere-PR.pdf", width = 15, height = 15, units = "cm")

#===============================================================================
#
#          top 10 family relative abundance bar plot
# 
#===============================================================================
# group taxa into family level and calculate relative abundance
ps.abund.fa = tax_glom(ps.tree.filtered.abund, taxrank = "Family")
ps.abund.fa.RA = transform_sample_counts(ps.abund.fa, function(x){x/sum(x)})
fa.otuTab = otu_table(ps.abund.fa.RA)
rownames(fa.otuTab) = tax_table(ps.abund.fa.RA)[rownames(fa.otuTab), "Family"]
fa.plotdf = cbind(t(fa.otuTab), sample_data(ps.abund.fa.RA))
fa.plotdf = fa.plotdf[, colnames(fa.plotdf)!='Not_Available']
write.table(fa.plotdf, "~/emmynoether/results/16S/family_RA_comp/family_otu_table.txt", 
            sep = "\t", row.names = F, quote = F)
# do Dunn test between compartments for each family
dunnTest.res.fa = data.frame()
dunnTest.label.fa = data.frame()
for (g in colnames(fa.plotdf)[1:118]) {
    formula1 = as.formula(paste(g, "Compartment", sep = "~"))
    krus.res = kruskal.test(formula1, data = fa.plotdf) 
    if(krus.res$p.value < 0.05 ){
        dunn.res = dunnTest(formula1, data = fa.plotdf, method = "bh")
        dunnTest.res.fa = rbind(dunnTest.res.fa, cbind(Family=g, dunn.res$res))
        # create significance label with letters
        label.dunn = cldList(P.adj ~ Comparison, data = dunn.res$res, threshold = 0.05)[, 1:2]
        dunnTest.label.fa = rbind(dunnTest.label.fa, cbind(Family=g, label.dunn))
    }
}
write.csv(dunnTest.res.fa, './diff_family_RA_compartment_dunntest_results.csv', row.names = F, quote = F)
write.csv(dunnTest.label.fa, './diff_family_RA_compartment_dunntest_label.csv', row.names = F, quote = F)

# find top 10 families
# remove family with 'not_available'
ps.fa.sub = subset_taxa(ps.tree.filtered.abund, Family != 'Not_Available')
fa.ra = ps.fa.sub %>% tax_glom(., taxrank = "Family") %>% otu_table(.) %>%
    as.data.frame(.) %>% mutate("RA"=rowSums(.)/sum(rowSums(.)))
top10.fa.names = rownames(fa.ra[order(fa.ra$RA, decreasing = T),])[1:10]

# prepare data for box plot
p2.df = pivot_longer(cbind(t(fa.otuTab[top10.fa.names, ]), sample_data(ps.abund.fa.RA)), 
                     cols = top10.fa.names, 
                     names_to = 'Family', values_to = 'RA')
p2.df$Compartment = factor(p2.df$Compartment, levels = c("SoilS", "AS", "BS", "AR", "BR", "CR", "SR"))

temp1 = ggplot(data= p2.df, aes(x = Family, y = RA, fill = Compartment)) +
    geom_boxplot() +
    facet_wrap(~Family, scale="free", nrow = 2, ncol = 5) +
    scale_fill_manual(values = c( 'brown', 'orange', 'chocolate1', 'seagreen', 'limegreen', 'royalblue', 'cornflowerblue'), 
                      limits = c("SoilS", "AS", "BS", "AR", "BR", "CR", "SR")) +
    mytheme

ggsave(temp1, file = "./results/16S/family_RA_comp/top_10_family_RA_compartment.pdf", 
       width = 30, height = 15, units = "cm")

#===============================================================================
#
#          Differential abundance analysis
#
#===============================================================================

setwd("~/emmynoether/results/16S/DA/")
ps.abund = readRDS("~/emmynoether/intermediate_data/16S/pre-analysis/ps.abund.RDS")

# differential abundance analysis between genotypes 
for (treat in c("C", "N", "P")) {
  if(treat=="C"){
    ps.abund.treat = subset_samples(ps.abund, Treatment=="C")
  }else if(treat=="N"){
    ps.abund.treat = subset_samples(ps.abund, Treatment=="N")
  }else{
    ps.abund.treat = subset_samples(ps.abund, Treatment=="P")
  }

  for (comp in c('AS', 'BS', 'CR', 'SR', "AR", "BR")) {
    if(comp=="AS"){
      ps.abund.treat.comp = subset_samples(ps.abund.treat, Compartment=="AS")
    }else if (comp=="BS"){
      ps.abund.treat.comp = subset_samples(ps.abund.treat, Compartment=="BS")
    }else if(comp=="CR"){
      ps.abund.treat.comp = subset_samples(ps.abund.treat, Compartment=="CR")
    }else if(comp=="SR"){
      ps.abund.treat.comp = subset_samples(ps.abund.treat, Compartment=="SR")
    }else if(comp=="AR"){
      ps.abund.treat.comp = subset_samples(ps.abund.treat, Compartment=="AR")
    }else{
      ps.abund.treat.comp = subset_samples(ps.abund.treat, Compartment=="BR")
    }
      # calculate relative abundance and remove rare OTUs
    ps.abund.treat.comp.ra = transform_sample_counts(ps.abund.treat.comp, function(x){x/sum(x)})
    abund.taxa = filter_taxa(ps.abund.treat.comp.ra, function(x) sum(x > 0.001) >= 1, TRUE)
    ps.input = subset_taxa(ps.abund.treat.comp, taxa_names(ps.abund.treat.comp)%in%taxa_names(abund.taxa))
    print(ps.input)
    deseq.abund = phyloseq_to_deseq2(ps.input, ~Genotype)
    dds = DESeq(deseq.abund, parallel = T)
    print(resultsNames(dds))
    # for each comparison, save results
    for(g in resultsNames(dds)[-1]){
      print("......genotype.......")
      print(g)
      res = results(dds, name = g, alpha = 0.05, test="Wald")
      print(summary(res))
      # subset significant genes
      res.Sig = subset(res, padj < 0.05 & abs(log2FoldChange) > 0)
      if(dim(res.Sig)[1] == 0){
        next
      }
      else{
        res.Sig.taxa = cbind(res.Sig, tax_table(ps.input)[rownames(res.Sig), ])
        write.table(res.Sig.taxa, paste0( paste(paste(paste0('./genotype_under_specific_treatment/significant_DEGs_betw_', g), comp, sep="_"), treat, sep="_"), 'treatment.txt'), 
                    quote = F, row.names = T, sep = "\t")
      }
      
    }
    
  }
  
}


#===============================================================================
#
#         indicator species analysis
#
#===============================================================================

library(indicspecies)

# group taxa into family level and calculate relative abundance
ps.abund.fa = tax_glom(ps.tree.filtered.abund, taxrank = "Family")
ps.abund.fa.RA = transform_sample_counts(ps.abund.fa, function(x){x/sum(x)})
ps.abund.fa.RA = subset_taxa(ps.abund.fa.RA, !Family %in% c('Not_Available') )
# only keep families with RA > 1% in at least 5% samples
ps.fa = filter_taxa(ps.abund.fa.RA, function(x) sum(x > 0.01) >= 0.05*nsamples(ps.abund.fe.RA) , TRUE)
fa.otutab = data.frame(otu_table(ps.fa))
rownames(fa.otutab) = tax_table(ps.fa)[, 'Family']
# indicator species analysis
#comb = combinespecies(as.matrix(t(fa.otutab)), max.order = 2)$XC
indvalspcomb = multipatt(as.matrix(t(fa.otutab)), unlist(as.vector(sample_data(ps.fa)[, 'Compartment'])), duleg = TRUE, 
                         control = how(nperm=1999))
summary(indvalspcomb, indvalcomp = TRUE)




