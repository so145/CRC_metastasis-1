#!/bin/bash
#$ -j y
#$ -pe smp 10 # 8 cores 
#$ -l h_vmem=32G
#$ -l h_rt=240:0:0 # to start request an hour
#$ -l highmem

module load anaconda3  R/4.2.0 hdf5/1.10.2
conda activate cpdb3_numpy1203

date
echo "got $NSLOTS cpu slots"

# run cellphonedb with statistical analysis and microenvironments file 
# CellphoneDB uses empirical shuffling to calculate which ligand–receptor pairs display significant cell-state specificity. Specifically, it estimates a null distribution of the mean of the average ligand and receptor expression in the interacting clusters by randomly permuting the cluster labels of all cells.
cellphonedb method statistical_analysis \
/data/BCI-CRC/nasrine/data/CRC/spatial/CRC_LM_VISIUM/CRC_LM_VISIUM_04_08_09_11/cellphonedb3/concat_withWu2022/meta.tsv \
/data/BCI-CRC/nasrine/data/CRC/spatial/CRC_LM_VISIUM/CRC_LM_VISIUM_04_08_09_11/cellphonedb3/concat_withWu2022/Multiome_Che_Wu_CRC_LM_counts_normalised_compatible_anndata.h5ad \
--microenvs /data/BCI-CRC/nasrine/data/CRC/spatial/CRC_LM_VISIUM/CRC_LM_VISIUM_04_08_09_11/cellphonedb3/concat_withWu2022/microenviroments_cell2loc_spatialde2.tsv \
--counts-data gene_name \
--threshold 0.1 \
--output-path /data/BCI-CRC/nasrine/data/CRC/spatial/CRC_LM_VISIUM/CRC_LM_VISIUM_04_08_09_11/cellphonedb3/concat_withWu2022/cpdb3_output/ \
--threads ${NSLOTS}

echo "finished job"
date 
