#!/bin/bash
#$ -j y
#$ -pe smp 10 # 8 cores 
#$ -l h_vmem=32G
#$ -l h_rt=1:0:0 # to start request an hour
#$ -l highmem

module load R/4.2.0

# plot heatmaps for microenvs 

Rscript /data/home/hfx941/code/CRC/st/CRC_LM_VISIUM/A10_cellphonedb/A3_plots/cpdb_plot_microenvs.R --microenvs /data/BCI-CRC/nasrine/data/CRC/spatial/CRC_LM_VISIUM/CRC_LM_VISIUM_04_08_09_11/cellphonedb3/concat_withWu2022/microenviroments_cell2loc_spatialde2.tsv --pvalues /data/BCI-CRC/nasrine/data/CRC/spatial/CRC_LM_VISIUM/CRC_LM_VISIUM_04_08_09_11/cellphonedb3/concat_withWu2022/cpdb3_output/pvalues.txt --out /data/BCI-CRC/nasrine/data/CRC/spatial/CRC_LM_VISIUM/CRC_LM_VISIUM_04_08_09_11/cellphonedb3/concat_withWu2022/heatmap/ --pvalue 0.05 

