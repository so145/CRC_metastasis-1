# The pipeline of a basic NicheNet analysis consist mainly of the following steps:
# 1. Define a “sender/niche” cell population and a “receiver/target” cell 
# population present in your expression data and determine which genes are 
# expressed in both populations
# 2. Define a gene set of interest: these are the genes in the “receiver/target”
# cell population that are potentially affected by ligands expressed by interacting 
# cells (e.g. genes differentially expressed upon cell-cell interaction)
# 3. Define a set of potential ligands: these are ligands that are expressed by 
# the “sender/niche” cell population and bind a (putative) receptor expressed
# by the “receiver/target” population
# 4. Perform NicheNet ligand activity analysis: rank the potential ligands based
# on the presence of their target genes in the gene set of interest (compared to 
# the background set of genes)
# 5. Infer top-predicted target genes of ligands that are top-ranked in the 
# ligand activity analysis

library(nichenetr)
library(Seurat) 
library(tidyverse)

microenv2consider = "all_celltypes" #"0" # to change depending on what you want # 0_with_Stem_NOTUM
NICHENET_DIR <- "/data/BCI-CRC/nasrine/data/CRC/spatial/CRC_LM_VISIUM/CRC_LM_VISIUM_04_08_09_11/nichenet/concat_withWu2022/"
DIR2LOAD <- paste0(NICHENET_DIR, "prepareInput/")
DIR2SAVE <- paste0(NICHENET_DIR, "nichenet_microenv",microenv2consider,"/", "intersect_cellphonedb/")
dir.create(DIR2SAVE)

# Load seurat object
seo <- readRDS(paste0(DIR2LOAD, "counts_microenv", microenv2consider,".rds"))
seo@meta.data %>% head()

seo@meta.data$Annotation_scVI_detailed %>% table()

# Maybe we have to normalise the data? yes
# Seurat object should contain normalized expression data (numeric matrix)
# according to this: https://rdrr.io/github/browaeysrobin/nichenetr/src/R/application_prediction.R
# get_expressed_genes requires normalised data
# We will also normalize the data to ensure that the count depth is equalized for
# all cells, given that we need the gene expression values to be comparable across the cell types.
seo <- NormalizeData(seo, verbose = TRUE)

# 1. Ligand-target model:
# This model denotes the prior potential that a particular ligand might regulate 
# the expression of a specific target gene. 
# Read in NicheNet’s ligand-target prior model, ligand-receptor network 
# and weighted integrated networks:
options(timeout = 600)
#lr_network = readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
#ligand_target_matrix = readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"))
#weighted_networks = readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final.rds"))

lr_network = readRDS("/data/BCI-CRC/nasrine/data/CRC/spatial/CRC_LM_VISIUM/CRC_LM_VISIUM_04_08_09_11/nichenet/databases/lr_network_human_21122021.rds")
ligand_target_matrix = readRDS("/data/BCI-CRC/nasrine/data/CRC/spatial/CRC_LM_VISIUM/CRC_LM_VISIUM_04_08_09_11/nichenet/databases/ligand_target_matrix_nsga2r_final.rds")
weighted_networks = readRDS("/data/BCI-CRC/nasrine/data/CRC/spatial/CRC_LM_VISIUM/CRC_LM_VISIUM_04_08_09_11/nichenet/databases/weighted_networks_nsga2r_final.rds")

# keep only unique/distinct rows from a data frame.
lr_network = lr_network %>% distinct(from, to)
head(lr_network)

ligand_target_matrix[1:5,1:5] # target genes in rows, ligands in columns

weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network, by = c("from","to"))
head(weighted_networks$lr_sig) 
# interactions and their weights in the ligand-receptor + signaling network

head(weighted_networks$gr) # interactions and their weights in the gene regulatory network

# If your expression data has the older gene symbols, you may want to use 
# our alias conversion function to avoid the loss of gene names
seo = alias_to_symbol_seurat(seo, "human")

# 2. Define a “sender/niche” cell population and a “receiver/target” cell 
# population present in expression data and define which genes are expressed 
# in sender and receiver populations
# We will consider a gene to be expressed when it is expressed in at least 10% 
# of cells in one cluster.

# need to change indent so it knows what is the "cluster annot")
seo <- SetIdent(seo, value = seo@meta.data$Annotation_scVI_detailed)

## receiver
receiver = "ipEMT" #"ipEMT" Stem (NOTUM high)
expressed_genes_receiver = get_expressed_genes(receiver, seo, pct = 0.10) #assay_oi = "RNA"
background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

## sender cell types : TO MANUALLY CHANGE depending on what you are using as signature
sender_celltypes = c("SPP1 Mac", "Neutrophil", "IL1B Mac", "NLRP3 Mac", "ECM CAF", "Myofibroblast", "CD8 Tex", "Treg")
# for AP1 signature: 
# c("SPP1 Mac", "Neutrophil", "IL1B Mac", "NLRP3 Mac", "ECM CAF", "Myofibroblast", "Pericyte")
# for AP1 regulon:
# c("SPP1 Mac", "Neutrophil", "IL1B Mac", "NLRP3 Mac", "ECM CAF", "Myofibroblast", "Pericyte")
# for NFKB regulon:
# c("SPP1 Mac", "Neutrophil", "IL1B Mac", "NLRP3 Mac", "ECM CAF", "Myofibroblast", "CD8 Tex", "Treg")
# for ipEMT signature:
# c("SPP1 Mac", "Neutrophil", "IL1B Mac", "NLRP3 Mac","ECM CAF", "Myofibroblast", "Pericyte", "CD8 Tex", "Treg")
# for MP17 interferon:
# c("SPP1 Mac", "Neutrophil", "IL1B Mac", "NLRP3 Mac","ECM CAF", "Myofibroblast", "Pericyte", "CD8 Tex", "Treg")

# lapply to get the expressed genes of every sender cell type separately here
list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, seo, 0.10) 
expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()

# Check the number of expressed genes: should be a 'reasonable' number of total 
# expressed genes in a cell type, e.g. between 5000-10000 (and not 500 or 20000)
length(expressed_genes_sender)
length(expressed_genes_receiver)

# 3. Define gene set of interest
sign2consider = "NFKB_regulon_combined" #"AP1_regulon_updated" # "AP1", "ipEMT" "AP1_regulon" "pEMT" "NFKB_regulon", "AP1_regulon_updated"

if (sign2consider %in% c("ipEMT", "pEMT")) {
  geneset_oi = readr::read_tsv(paste0(DIR2LOAD, sign2consider,"_signatures.txt"), col_names = FALSE) # pEMT_signatures.txt
  head(geneset_oi)
  geneset_oi = geneset_oi$X1
}

if (sign2consider %in% c("AP1")) {
  geneset_oi = read.csv(paste0("/data/BCI-CRC/nasrine/data/gene_sets/", sign2consider,"_signature.csv"))
  geneset_oi = geneset_oi$Shared
}

if (sign2consider %in% c("AP1_regulon")) {
  geneset_oi = readr::read_tsv(paste0(DIR2LOAD,"scenicplus_AP1_regulon.txt"), col_names = FALSE)
  head(geneset_oi)
  geneset_oi = geneset_oi$X1
}

if (sign2consider %in% c("NFKB_regulon")) {
  geneset_oi = readr::read_tsv(paste0(DIR2LOAD,"scenicplus_NFKB_regulon.txt"), col_names = FALSE)
  head(geneset_oi)
  geneset_oi = geneset_oi$X1
}

if (sign2consider %in% c("NFKB_regulon_combined")) {
  geneset_oi = readr::read_tsv(paste0(DIR2LOAD,"scenicplus_NFKB_regulon_combined.txt"), col_names = FALSE)
  head(geneset_oi)
  geneset_oi = geneset_oi$X1
}

if (sign2consider %in% c("MP17_interferon")) {
  library("readxl")
  geneset_oi = read_excel("/data/BCI-CRC/nasrine/data/gene_sets/gavishHallmarksTranscriptio2023_41586_2023_6130_MOESM6_ESM.xlsx", 
                          sheet = "Cancer MPs",
                          col_names = TRUE)
  head(geneset_oi)
  geneset_oi = geneset_oi$`MP17 Interferon/MHC-II (I)`
}

if (sign2consider %in% c("AP1_regulon_updated")) {
  geneset_oi = readr::read_tsv(paste0(DIR2LOAD,"AP1targetGenes_activators.txt"), col_names = FALSE)
  head(geneset_oi)
  geneset_oi = geneset_oi$X1
}


# 4.  Define a set of potential ligands: these are ligands that are expressed 
# by the “sender/niche” cell population and bind a (putative) receptor expressed 
# by the “receiver/target” population
# As potentially active ligands, we will use ligands that are 
# 1) expressed by sender and 2) can bind a (putative) receptor expressed by receiver. 
# Putative ligand-receptor links were gathered from NicheNet’s ligand-receptor data sources.

ligands = lr_network %>% pull(from) %>% unique()
expressed_ligands = intersect(ligands,expressed_genes_sender)

receptors = lr_network %>% pull(to) %>% unique()
expressed_receptors = intersect(receptors,expressed_genes_receiver)

lr_network_expressed = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) 
head(lr_network_expressed)

# This ligand-receptor network contains the expressed ligand-receptor interactions.
# As potentially active ligands for the NicheNet analysis, we will consider the
# ligands from this network.
potential_ligands = lr_network_expressed %>% pull(from)  %>% unique()
head(potential_ligands)

# 5. Perform NicheNet’s ligand activity analysis on the gene set of interest
# we will calculate the ligand activity of each ligand, i.e. we will assess how
# well each sender-ligand can predict the signature gene set compared to the 
# background of expressed genes (predict whether a gene belongs to the signature program or not)
ligand_activities = predict_ligand_activities(geneset = geneset_oi, 
                                              background_expressed_genes = background_expressed_genes,
                                              ligand_target_matrix = ligand_target_matrix, 
                                              potential_ligands = potential_ligands)

# Now, we want to rank the ligands based on their ligand activity. 
# In our validation study, we showed that the area under the precision-recall 
# curve (AUPR) between a ligand’s target predictions and the observed transcriptional
# response was the most informative measure to define ligand activity 
# (this was the Pearson correlation for v1). Therefore, we will rank the ligands 
# based on their AUPR. This allows us to prioritize p-EMT-regulating ligands.
ligand_activities = ligand_activities %>% arrange(-aupr_corrected) %>% mutate(rank = rank(desc(aupr_corrected)))
ligand_activities

# NEW: select only ligands that are significant in cellphoneDB
if (sign2consider %in% c("AP1_regulon")) {
  common_ligands = readr::read_tsv(paste0(DIR2SAVE,sign2consider,receiver,"_common_ligands_4nichenet.txt"), # _common_ligands.txt
                                   col_names = FALSE)
}

if (sign2consider %in% c("AP1_regulon_updated")) {
  common_ligands = readr::read_tsv(paste0(DIR2SAVE,sign2consider,receiver,"_common_ligands_4nichenet.txt"), # _common_ligands.txt
                                   col_names = FALSE)
}

if (sign2consider %in% c("NFKB_regulon")) {
  common_ligands = readr::read_tsv(paste0(DIR2SAVE,sign2consider,receiver,"_common_ligands_4nichenet.txt"), # _common_ligands.txt
                                   col_names = FALSE)
}

if (sign2consider %in% c("NFKB_regulon_combined")) {
  common_ligands = readr::read_tsv(paste0(DIR2SAVE,sign2consider,receiver,"_common_ligands_4nichenet.txt"), # _common_ligands.txt
                                   col_names = FALSE)
}

# only select common ligands 
ligand_activities = ligand_activities[ligand_activities$test_ligand %in% common_ligands$X1, ]

# The number of top-ranked ligands that are further used to predict active target genes 
# and construct an active ligand-receptor network is here 30.
nligands = dim(ligand_activities)[1] # should be set to ligand_activities.shape[0]
best_upstream_ligands = ligand_activities %>% top_n(nligands, aupr_corrected) %>% arrange(-aupr_corrected) %>% pull(test_ligand) %>% unique()

# These ligands are expressed by one or more of the input sender cells. 
# To see which cell population expresses which of these top-ranked ligands, 
# you can run the following
pdf(paste0(DIR2SAVE,sign2consider,receiver, "_nligands", nligands,"_dotplot_ligand_expression_sender_cells.pdf"), 
    height=4.5, width=20)
lig_dotplot <- DotPlot(seo, features = best_upstream_ligands %>% rev(), 
                       cols = "RdBu", idents=c( receiver, sender_celltypes)) + RotatedAxis() 
print(lig_dotplot)
dev.off()


# 6.  Infer target genes of top-ranked ligands and visualize in a heatmap
# look at the regulatory potential scores between ligands and target genes of interest. 
# In the ligand-target heatmaps, we show here regulatory potential scores for 
# interactions between the 30 top-ranked ligands and following target genes: 
# genes that belong to the gene set of interest and to the 250 most strongly 
# predicted targets of at least one of the 30 top-ranked ligands 
active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,
                                                                 geneset = geneset_oi,
                                                                 ligand_target_matrix = ligand_target_matrix, 
                                                                 n = 250) %>% bind_rows() %>% drop_na()
nrow(active_ligand_target_links_df)
head(active_ligand_target_links_df)

# For visualization purposes, we adapted the ligand-target regulatory potential 
# matrix as follows. Regulatory potential scores were set as 0 if their score
# was below a predefined threshold, which was here the 0.25 quantile of scores 
# of interactions between the 30 top-ranked ligands and each of their respective top targets 
active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df,
                                                                 ligand_target_matrix = ligand_target_matrix, 
                                                                 cutoff = 0.25)
# returns df with rownames = targets, colnames = ligands

# order the ligands according to the ranking according to the ligand activity prediction.
order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23

vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()

# save to file the different target links
write.table(active_ligand_target_links_df, file=paste0(DIR2SAVE,sign2consider,receiver, "_nligands", nligands,"_active_ligand_target_links_df.csv"), 
            row.names = FALSE, sep=",",  col.names=TRUE)
write.table(active_ligand_target_links, file=paste0(DIR2SAVE,sign2consider,receiver, "_nligands", nligands,"_active_ligand_target_links_matrix.csv"), 
            row.names = TRUE, sep=",",  col.names=TRUE)
write.table(vis_ligand_target, file=paste0(DIR2SAVE,sign2consider,receiver, "_nligands", nligands,"_vis_ligand_target.csv"), 
            row.names = TRUE, sep=",",  col.names=TRUE)

p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized ligands",
                                                                    "Predicted target genes", 
                                                                    color = "hotpink",
                                                                    legend_position = "top", 
                                                                    x_axis_position = "top",
                                                                    legend_title = "Regulatory potential")  + theme(axis.text.x = element_text(face = "italic")) + scale_fill_gradient2(low = "whitesmoke",  high = "hotpink", breaks = c(0,0.0045,0.0090))
p_ligand_target_network

pdf(paste0(DIR2SAVE,sign2consider,receiver,"_nligands", nligands,"_dotplot_ligand_target_network.pdf"), height=12, width=50)
print(p_ligand_target_network)
dev.off()

# save ligands and target genees in txt file
#write.table(order_ligands, file=paste0(DIR2SAVE,sign2consider,receiver, "_nligands", nligands,"_ligands.csv"), 
            #row.names = FALSE, sep=",",  col.names=FALSE)
#write.table(order_targets, file=paste0(DIR2SAVE,sign2consider,receiver, "_nligands", nligands,"_targets.csv"), 
            #row.names = FALSE, sep=",",  col.names=FALSE)

# 7. Ligand-receptor network inference for top-ranked ligands
# looking at which receptors of the receiver cell population can potentially bind 
# to the prioritized ligands from the sender cell population 

# get the ligand-receptor network of the top-ranked ligands
lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()

# save ligand-receptor network to csv 
write.table(lr_network_top, file=paste0(DIR2SAVE,sign2consider,receiver, "_nligands", nligands,"_LRnetwork.csv"), 
            row.names = FALSE, sep=",",  col.names=TRUE)

# get the weights of the ligand-receptor interactions as used in the NicheNet model
lr_network_top_df_large = weighted_networks_lr %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)
lr_network_top_df = lr_network_top_df_large %>% spread("from","weight",fill = 0)
# convert to a matrix
lr_network_top_matrix = lr_network_top_df %>% dplyr::select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)

# perform hierarchical clustering to order the ligands and receptors
dist_receptors = dist(lr_network_top_matrix, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]

dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix))
order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix))

vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
rownames(vis_ligand_receptor_network) = order_receptors %>% make.names()
colnames(vis_ligand_receptor_network) = order_ligands_receptor %>% make.names()

# maybe save this matrix as csv to be able to sort interactions by strength?
write.table(vis_ligand_receptor_network, file=paste0(DIR2SAVE,sign2consider,receiver, "_nligands", nligands,"_vis_ligand_receptor_network.csv"), 
            row.names = TRUE, sep=",",  col.names=TRUE)

# reorder ligands for LR in function of order of ligands in the ligand - target plot
if (setequal(order_ligands, order_ligands_receptor)) {
  vis_ligand_receptor_network = vis_ligand_receptor_network[order_receptors, order_ligands]
}


p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "mediumorchid1", x_axis_position = "top",legend_title = "Prior interaction potential")
p_ligand_receptor_network

pdf(paste0(DIR2SAVE,sign2consider,receiver, "_nligands", nligands,"_dotplot_ligand_receptor_network.pdf"), height=7, width=10)
print(p_ligand_receptor_network)
dev.off()

# write receptors to file
#write.table(order_receptors, file=paste0(DIR2SAVE,sign2consider,receiver, "_nligands", nligands,"_receptors.csv"), 
            #row.names = FALSE, sep=",",  col.names=FALSE)

# if we want too do a plot together order_ligands and order_ligands_receptor need 
# to have same order i think. we have to first make sure they are the same vectors
# i.e contain same elements even if not in same ordr  (done see above with setequal)

### Additional visualisation 
library(RColorBrewer)
library(cowplot)
library(ggpubr)

# ligand activity
ligand_aupr_matrix = ligand_activities %>% dplyr::select(aupr_corrected) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities$test_ligand)

vis_ligand_aupr = ligand_aupr_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("AUPR")
p_ligand_aupr = vis_ligand_aupr %>% make_heatmap_ggplot("Prioritized ligands","Ligand activity", color = "sienna1",legend_position = "top", x_axis_position = "top", legend_title = "AUPR\n(target gene prediction ability)")
p_ligand_aupr
pdf(paste0(DIR2SAVE,sign2consider,receiver, "_nligands", nligands,"_aupr_ligand.pdf"), height=10, width=4.5)
print(p_ligand_aupr)
dev.off()
