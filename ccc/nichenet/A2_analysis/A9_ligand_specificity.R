library(Seurat) 
library(tidyverse)
library(circlize)

microenv2consider = "all_celltypes" 
NICHENET_DIR <- "/data/BCI-CRC/nasrine/data/CRC/spatial/CRC_LM_VISIUM/CRC_LM_VISIUM_04_08_09_11/nichenet/concat_withWu2022/"
DIR2LOAD <- paste0(NICHENET_DIR, "prepareInput/")
DIR2SAVE <- paste0(NICHENET_DIR, "nichenet_microenv",microenv2consider,"/", "intersect_cellphonedb/")
dir.create(DIR2SAVE)

sign2consider = "NFKB_regulon_combined"
nligands = 63 
receiver = "ipEMT"

## receiver
receiver = "ipEMT" 

## sender cell types : TO MANUALLY CHANGE depending on what you are using as signature
sender_celltypes = c("SPP1 Mac", "Neutrophil", "IL1B Mac", "NLRP3 Mac", 
                     "ECM CAF", "Myofibroblast", "CD8 Tex", "Treg")
# to change depending on regulon

# load target gene links
active_ligand_target_links_df = read.csv(paste0(DIR2SAVE,sign2consider,receiver, "_nligands", nligands,"_active_ligand_target_links_df.csv"))

# Maybe we have to normalise the data? yes

# We will also normalize the data to ensure that the count depth is equalized for
# all cells, given that we need the gene expression values to be comparable across the cell types.


### Calculate average ligand expression in sender cells

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
CD8Tex_specific_ligands = sender_ligand_assignment$"CD8 Tex" %>% names() %>% setdiff(general_ligands)
Treg_specific_ligands = sender_ligand_assignment$"Treg" %>% names() %>% setdiff(general_ligands)

ligand_type_indication_df = tibble(
  ligand_type = c(rep("SPP1 Mac-specific", times = SPP1Mac_specific_ligands %>% length()),
                  rep("Neutrophil-specific", times = Neutrophil_specific_ligands %>% length()),
                  rep("IL1B Mac-specific", times = IL1BMac_specific_ligands %>% length()),
                  rep("NLRP3 Mac-specific", times = NLRP3Mac_specific_ligands %>% length()),
                  rep("ECM CAF-specific", times = ECMCAF_specific_ligands %>% length()),
                  rep("Myofibroblast-specific", times = Myofibroblast_specific_ligands %>% length()),
                  rep("CD8Tex-specific", times = CD8Tex_specific_ligands %>% length()),
                  rep("Treg-specific", times = Treg_specific_ligands %>% length()),
                  rep("Common", times = general_ligands %>% length())),
  ligand = c(SPP1Mac_specific_ligands, Neutrophil_specific_ligands, IL1BMac_specific_ligands, 
             NLRP3Mac_specific_ligands, ECMCAF_specific_ligands, 
             Myofibroblast_specific_ligands, CD8Tex_specific_ligands, Treg_specific_ligands,
             general_ligands))


write.table(ligand_type_indication_df, file=paste0(DIR2SAVE,sign2consider,receiver, "_ligand_specificity_senders.csv"), 
            row.names = FALSE, sep=",",  col.names=TRUE)
