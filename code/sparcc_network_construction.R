library(phyloseq)
library(SpiecEasi)
library(igraph)

#===============================================================================
#
#             SparCC network construction         
#
#===============================================================================

# read bacteria data 
bac.ps = readRDS("~/emmynoether/intermediate_data/16S/pre-analysis/ps.abund.RDS")
# do the following for each compartment
bac.comp = subset_samples(bac.ps, Compartment=='AR') 
# keep taxa with relative abundance > 0.1% at least in 10% samples
bac.RA = transform_sample_counts(bac.comp, function(x) x/sum(x))
bac.RA.abund = filter_taxa(bac.RA, 
                           function(x) sum(x > 0.001) >= 0.1*nsamples(bac.RA), TRUE)  
bac.comp.abund = subset_taxa(bac.comp, taxa_names(bac.comp)%in%taxa_names(bac.RA.abund))

# prepare OTU table and taxonomy table 
bac.comp.otuTab = as.data.frame(otu_table(bac.comp.abund))
bac.comp.taxTab = as.data.frame(tax_table(bac.comp.abund))
# SparCC analysis and do 100 bootstraps to get p value 
sparcc.net = sparcc(as.matrix(t(bac.comp.otuTab)), iter = 100, inner_iter = 50)
sparcc.bt = sparccboot(as.matrix(t(bac.comp.otuTab)), ncpus = 16, R=100)
sparcc.pv = pval.sparccboot(sparcc.bt, sided = "both")
# get correlation data frame and p value 
pvals = sparcc.pv$pvals
sparcc.cor = sparcc.net$Cor
rownames(sparcc.cor) = colnames(sparcc.cor) = rownames(bac.comp.otuTab)
ind = which(upper.tri(sparcc.cor, diag = F), arr.ind = TRUE)
dn = dimnames(sparcc.cor)
cor.mat = data.frame(row = dn[[1]][ind[, 1]],
                     col = dn[[2]][ind[, 2]],
                     corr.sp = sparcc.cor[ind])
cor.mat$pval = pvals
write.csv(cor.mat, './AR/all_sparcc_cor_AR.csv', row.names = F, quote = F)

# Define arbitrary threshold for SparCC correlation matrix for the graph
cor.sig = cor.mat[abs(cor.mat$corr.sp) > 0.4 & cor.mat$pval < 0.05, ]
# Create igraph objects and output edges list
ig.sparcc = graph_from_data_frame(cor.sig, directed = F)
edges.AR = as_edgelist(ig.sparcc)
edges.AR = cbind(edges.AR, cor.sig$corr.sp)
colnames(edges.AR) = c('Source', 'Target', 'Weight')
edges.AR = as.data.frame(edges.AR)
write.csv(edges.AR, './AR/OTU-net-AR-edgelist.csv', quote = F, row.names = F)
# output nodes list
vsize = apply(otu_table(bac.RA.abund), 1, mean) 
names(vsize) = taxa_names(bac.RA.abund)
nodes.AR = data.frame(ID=V(ig.sparcc)$name, 
                      meanRA=vsize[V(ig.sparcc)$name], 
                      hubScore=hub_score(ig.sparcc)$vector, 
                      Degree=igraph::degree(ig.sparcc))
nodes.AR = cbind(nodes.AR, bac.comp.taxTab[nodes.AR$ID, ])
write.csv(nodes.AR, './AR/OTU-net-AR-nodelist.csv', quote = F, row.names = F)




