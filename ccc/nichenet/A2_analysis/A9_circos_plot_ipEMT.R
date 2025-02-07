library(nichenetr)
library(Seurat) 
library(tidyverse)
library(circlize)

microenv2consider = "all_celltypes" 
NICHENET_DIR <- "/data/BCI-CRC/nasrine/data/CRC/spatial/CRC_LM_VISIUM/CRC_LM_VISIUM_04_08_09_11/nichenet/concat_withWu2022/"
DIR2LOAD <- paste0(NICHENET_DIR, "prepareInput/")
DIR2SAVE <- paste0(NICHENET_DIR, "nichenet_microenv",microenv2consider,"/", "intersect_cellphonedb/")
dir.create(DIR2SAVE)

sign2consider = "AP1_regulon_updated"
nligands = 54 
receiver = "ipEMT"

## receiver
receiver = "ipEMT" 

## sender cell types : TO MANUALLY CHANGE depending on what you are using as signature
sender_celltypes = c("SPP1 Mac", "Neutrophil", "IL1B Mac", "NLRP3 Mac", 
                     "ECM CAF", "Myofibroblast", "Pericyte")

# load target gene links
active_ligand_target_links_df = read.csv(paste0(DIR2SAVE,sign2consider,receiver, "_nligands", nligands,"_active_ligand_target_links_df.csv"))

# Load seurat object
#seo <- readRDS(paste0(DIR2LOAD, "counts_microenv", microenv2consider,".rds"))
#seo@meta.data %>% head()

#seo@meta.data$Annotation_scVI_detailed %>% table()

# Maybe we have to normalise the data? yes
# Seurat object should contain normalized expression data (numeric matrix)
# according to this: https://rdrr.io/github/browaeysrobin/nichenetr/src/R/application_prediction.R
# get_expressed_genes requires normalised data
# We will also normalize the data to ensure that the count depth is equalized for
# all cells, given that we need the gene expression values to be comparable across the cell types.
#seo <- NormalizeData(seo, verbose = TRUE)

# If your expression data has the older gene symbols, you may want to use 
# our alias conversion function to avoid the loss of gene names
#seo = alias_to_symbol_seurat(seo, "human")

# need to change indent so it knows what is the "cluster annot")
#seo <- SetIdent(seo, value = seo@meta.data$Annotation_scVI_detailed)

### Calculate average ligand expression in sender cells
#avg_expression_ligands = AverageExpression(seo %>% subset(subset = Annotation_scVI_detailed %in% sender_celltypes), 
                                           #features = active_ligand_target_links_df$ligand %>% unique(),
                                           #slot="data"
                                           #)
#avg_expression_ligands = AverageExpression(seo, features = active_ligand_target_links_df$ligand %>% unique())
avg_expression_ligands = read.csv(paste0(DIR2SAVE,sign2consider,receiver, "_final_LR_1Ligand_average_geneexp.csv"),
                                    row.names = 1, check.names=FALSE) 
# check.names=FALSE because want to import column names with spaces and keep them

### Assign ligands to sender cells
# To assign ligands to sender cell type, we can e.g. look for which sender cell 
# types show an expression that is higher than the average + SD.
#sender_ligand_assignment = avg_expression_ligands$RNA %>% apply(1, function(ligand_expression){
 # ligand_expression > (ligand_expression %>% mean() + ligand_expression %>% sd())
#}) %>% t()

sender_ligand_assignment = avg_expression_ligands %>% apply(1, function(ligand_expression){
  ligand_expression > (ligand_expression %>% mean() + ligand_expression %>% sd())
}) %>% t()


sender_ligand_assignment = sender_ligand_assignment %>% apply(2, function(x){x[x == TRUE]}) %>% purrr::keep(function(x){length(x) > 0})
names(sender_ligand_assignment)

# We will know also look at which ligands are common across multiple cell types
# (= those that are specific to > 1 cell type, or those that were not assigned to 
# a cell type in the previous block of code)
# Determine now which prioritized ligands are expressed by sender cells
all_assigned_ligands = sender_ligand_assignment %>% lapply(function(x){names(x)}) %>% unlist()
unique_ligands = all_assigned_ligands %>% table() %>% .[. == 1] %>% names()
general_ligands = active_ligand_target_links_df$ligand %>% unique() %>% setdiff(unique_ligands)

SPP1Mac_specific_ligands = sender_ligand_assignment$"SPP1 Mac" %>% names() %>% setdiff(general_ligands)
Neutrophil_specific_ligands = sender_ligand_assignment$"Neutrophil" %>% names() %>% setdiff(general_ligands)
IL1BMac_specific_ligands = sender_ligand_assignment$"IL1B Mac" %>% names() %>% setdiff(general_ligands)
NLRP3Mac_specific_ligands = sender_ligand_assignment$"NLRP3 Mac" %>% names() %>% setdiff(general_ligands)
ECMCAF_specific_ligands = sender_ligand_assignment$"ECM CAF" %>% names() %>% setdiff(general_ligands)
Myofibroblast_specific_ligands = sender_ligand_assignment$"Myofibroblast" %>% names() %>% setdiff(general_ligands)
Pericyte_specific_ligands = sender_ligand_assignment$"Pericyte" %>% names() %>% setdiff(general_ligands)

ligand_type_indication_df = tibble(
  ligand_type = c(rep("SPP1 Mac-specific", times = SPP1Mac_specific_ligands %>% length()),
                  rep("Neutrophil-specific", times = Neutrophil_specific_ligands %>% length()),
                  rep("IL1B Mac-specific", times = IL1BMac_specific_ligands %>% length()),
                  rep("NLRP3 Mac-specific", times = NLRP3Mac_specific_ligands %>% length()),
                  rep("ECM CAF-specific", times = ECMCAF_specific_ligands %>% length()),
                  rep("Myofibroblast-specific", times = Myofibroblast_specific_ligands %>% length()),
                  rep("Pericyte-specific", times = Pericyte_specific_ligands %>% length()),
                  rep("Common", times = general_ligands %>% length())),
  ligand = c(SPP1Mac_specific_ligands, Neutrophil_specific_ligands, IL1BMac_specific_ligands, 
             NLRP3Mac_specific_ligands, ECMCAF_specific_ligands, 
             Myofibroblast_specific_ligands, Pericyte_specific_ligands,
             general_ligands))

write.table(ligand_type_indication_df, file=paste0(DIR2SAVE,sign2consider,receiver, "_ligand_specificity_senders.csv"), 
            row.names = FALSE, sep=",",  col.names=TRUE)

# Define the ligand-target links of interest
# some filtering on targets based on some of weights from ligands (THIS IS OPTIONAL)
# so we don't display all target genes
# there are two ways to do this option1: 
# todo change based on which option you want to have
opt2chose = "option1" # option1 is myway, option2 is nichenetway

if (opt2chose == "option1") {
  my_thresh = 0.8
  target_filter_df = active_ligand_target_links_df %>% group_by(target) %>% 
    summarise(sum_potential=sum(weight), .groups="drop")
  target_filter_df = target_filter_df %>% filter(sum_potential > my_thresh)
}

# To avoid making a circos plots with too many ligand-target links, 
# we will show only links with a weight higher than a predefined cutoff: 
# links belonging to the 40% of lowest scores were removed.
active_ligand_target_links_df = active_ligand_target_links_df %>% mutate(target_type = "AP1 regulon") %>% inner_join(ligand_type_indication_df)

cutoff_include_all_ligands = active_ligand_target_links_df$weight %>% quantile(0.4) # 0.66

if (opt2chose == "option2") {
  active_ligand_target_links_df_circos = active_ligand_target_links_df %>% filter(weight > cutoff_include_all_ligands)
}
if (opt2chose == "option1") {
  active_ligand_target_links_df_circos = active_ligand_target_links_df %>% filter(target %in% target_filter_df$target)
}

ligands_to_remove = setdiff(active_ligand_target_links_df$ligand %>% unique(), active_ligand_target_links_df_circos$ligand %>% unique())
targets_to_remove = setdiff(active_ligand_target_links_df$target %>% unique(), active_ligand_target_links_df_circos$target %>% unique())

circos_links = active_ligand_target_links_df %>% filter(!target %in% targets_to_remove &!ligand %in% ligands_to_remove)

# Prepare the circos visualization: give each segment of ligands and targets a specific color and order

grid_col_ligand =c("Common" = "lavender", # "lavenderblush2" "darkolivegreen2", #"darkolivegreen1",#"#ffd8e5",
                   "SPP1 Mac-specific" = "#8dd3c7",
                   "Neutrophil-specific" = "#1f78b4",
                   "IL1B Mac-specific" = "#fb8072",
                   "NLRP3 Mac-specific" = "#fdb462",
                   "ECM CAF-specific" = "hotpink",  #"#fe4fe1"
                   "Myofibroblast-specific" = "#4fe1fe",
                   "Pericyte-specific" = "#4f8afe"
                   )
grid_col_target =c(
  "AP1 regulon" = "gray")

grid_col_tbl_ligand = tibble(ligand_type = grid_col_ligand %>% names(), color_ligand_type = grid_col_ligand)
grid_col_tbl_target = tibble(target_type = grid_col_target %>% names(), color_target_type = grid_col_target)

circos_links = circos_links %>% mutate(ligand = paste(ligand," ")) # extra space: make a difference between a gene as ligand and a gene as target!

# this step is not working
circos_links = circos_links %>% inner_join(grid_col_tbl_ligand) %>% inner_join(grid_col_tbl_target)
links_circle = circos_links %>% select(ligand,target, weight)

ligand_color = circos_links %>% distinct(ligand,color_ligand_type)
grid_ligand_color = ligand_color$color_ligand_type %>% set_names(ligand_color$ligand)
target_color = circos_links %>% distinct(target,color_target_type)
grid_target_color = target_color$color_target_type %>% set_names(target_color$target)

grid_col =c(grid_ligand_color,grid_target_color)

# give the option that links in the circos plot will be transparant ~ ligand-target potential score
transparency = circos_links %>% mutate(weight =(weight-min(weight))/(max(weight)-min(weight))) %>% mutate(transparency = 1-weight) %>% .$transparency 

# Prepare the circos visualization: order ligands and targets
target_order = circos_links$target %>% unique()
ligand_order = c(SPP1Mac_specific_ligands, Neutrophil_specific_ligands,
                 IL1BMac_specific_ligands, NLRP3Mac_specific_ligands,
                 ECMCAF_specific_ligands, Myofibroblast_specific_ligands,
                 Pericyte_specific_ligands,general_ligands) %>% c(paste(.," ")) %>% intersect(circos_links$ligand)
order = c(ligand_order,target_order)

# Prepare the circos visualization: define the gaps between the different segments
width_same_cell_same_ligand_type = 1.5 #0.5
width_different_cell = 5 #6 distance between cell types
width_ligand_target = 15
width_same_cell_same_target_type = 0.99 #0.5

gaps = c(
  # width_ligand_target,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "SPP1 Mac-specific") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "Neutrophil-specific") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "IL1B Mac-specific") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "NLRP3 Mac-specific") %>% distinct(ligand) %>% nrow() -1)),
  #width_different_cell,
  #rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "ECM CAF-specific") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "Myofibroblast-specific") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "Pericyte-specific") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "Common") %>% distinct(ligand) %>% nrow() -1)),
  width_ligand_target,
  rep(width_same_cell_same_target_type, times = (circos_links %>% filter(target_type == "AP1 regulon") %>% distinct(target) %>% nrow() -1)),
  width_ligand_target
)

# Render the circos plot (all links same transparancy).
# Only the widths of the blocks that indicate each target gene is proportional 
# the ligand-target regulatory potential (~prior knowledge supporting the regulatory interaction).
FIG2SAVE <- paste0(DIR2SAVE, "figures/final/", sign2consider, "_", receiver, "/")

svg(paste0(FIG2SAVE,opt2chose,"ligand_target_circos.svg"), width = 10, height = 10)
circos.par(gap.degree = gaps)
chordDiagram(links_circle, directional = 1,
             order=order,link.sort = TRUE, 
             link.decreasing = FALSE, 
             grid.col = grid_col,
             transparency = 0, 
             diffHeight = 0.005, 
             direction.type = c("diffHeight", "arrows"),
             link.arr.type = "big.arrow", 
             link.visible = links_circle$weight >= cutoff_include_all_ligands,
             annotationTrack = "grid",
             preAllocateTracks = list(track.height = 0.075))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 1)
}, bg.border = NA) #
circos.clear()
dev.off()

svg(paste0(FIG2SAVE,opt2chose,"ligand_target_circos_with_transparency.svg"), width = 10, height = 10)
circos.par(gap.degree = gaps)
chordDiagram(links_circle, directional = 1,
             order=order,link.sort = TRUE, 
             link.decreasing = FALSE, 
             grid.col = grid_col,
             transparency = transparency, 
             diffHeight = 0.005, 
             direction.type = c("diffHeight", "arrows"),
             link.arr.type = "big.arrow", 
             link.visible = links_circle$weight >= cutoff_include_all_ligands,
             annotationTrack = "grid",
             preAllocateTracks = list(track.height = 0.075))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 1)
}, bg.border = NA) #
circos.clear()
dev.off()

library(ComplexHeatmap)

# create legend
# only select cell types to displayy that have specific ligands or common ligands
# i.e. ECM CAF has no specific ligand 
grid_col_tbl_ligand2plot = grid_col_tbl_ligand %>% subset(subset=ligand_type %in% (ligand_type_indication_df$ligand_type %>% unique()))

pdf(paste0(FIG2SAVE,opt2chose,"ligand_target_circos_legend.pdf"), width = 10, height = 10)
L1 = Legend(labels = grid_col_tbl_ligand2plot$ligand_type,
            type = "points", pch = as.character(""),
            legend_gp = gpar(col = "white", cex = 0.7),
            background = grid_col_tbl_ligand2plot$color_ligand_type)
draw(L1)
dev.off()

