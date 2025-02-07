# .libPaths("/data/home/hfx941/R/x86_64-pc-linux-gnu-library/4.2")
library(Matrix)
library(Seurat)
source("~/code/R_utils/create_seurat_object.R")

microenv = "all_celltypes" # 0_with_Stem_NOTUM, 0
name = "all_celltypes" # Stem_NOTUM_ipEMT_microenv, "ipEMT_microenv"
 
# DIR2LOAD mtx data from
DIR2LOAD <- "/data/BCI-CRC/nasrine/data/CRC/spatial/CRC_LM_VISIUM/CRC_LM_VISIUM_04_08_09_11/nichenet/concat_withWu2022/prepareInput/"
counts_mtx <- paste0("counts_microenv", microenv, ".mtx")
obs_df <- paste0("obs_microenv", microenv, ".csv")
var_df <- paste0("var_microenv", microenv, ".csv")
DIR2SAVE <- "/data/BCI-CRC/nasrine/data/CRC/spatial/CRC_LM_VISIUM/CRC_LM_VISIUM_04_08_09_11/nichenet/concat_withWu2022/prepareInput/"

seo <- Anndata2Seurat(
  input_dir = DIR2LOAD,
  counts_mtx = counts_mtx, # matrix needs to be cellsxgenes here (not already transposed!)
  obs_df = obs_df,
  var_df = var_df,
  project_name = name #"ipEMT_microenv"
  
)
head(seo@meta.data)
rownames(seo@meta.data)

print("Saving the Seurat object as .rds file:")
saveRDS(seo, paste0(DIR2SAVE,"counts_microenv", microenv, ".rds"))
