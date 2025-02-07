library(GSVA)

# Human annotation package we'll use for gene identifier conversion
library(org.Hs.eg.db)

genes2consider = 10 #100 40 # nb genes for each signature in geneset
DIR2LOAD <- "/data/BCI-CRC/nasrine/data/CRC/Primary_CRC_dataset/final_object/20mt/GSVA/filtered_signatures/" #iREC_AP1regulon/"

geneset = readr::read_tsv(paste0(DIR2LOAD, "gene_signatures_", genes2consider, "_DE_genes.tsv"))
# log normalised data
#tcga_expression = readr::read_tsv(paste0(DIR2LOAD,"TCGA_COADREAD_TPonly_uniquePatients_zscore.tsv"))
tcga_expression = readr::read_tsv(paste0(DIR2LOAD,"TCGA_COADREAD_TPonly_VST2.tsv"))
tcga_expression <- tcga_expression %>%
# Here we are going to store the gene IDs as row names so that we can have a numeric matrix to perform calculations on later
tibble::column_to_rownames("index")

# set gene name as row id?
#tcga_expression <- tcga_expression %>%
  # Here we are going to store the gene IDs as row names so that we can have a numeric matrix to perform calculations on later
  #tibble::column_to_rownames("index")

# The function that we will use to run GSVA wants the gene sets to be in a list, 
# where each entry in the list is a vector of genes that comprise the pathway the element is named for.
# To make this into the list format we need, we can use the split() function. 
# We want a list where each element of the list is a vector that contains 
# the Entrez gene IDs/gene names that are in a particular pathway set.
hallmarks_list <- split(
  geneset$member, # The genes we want split into pathways
  geneset$name # The pathways made as the higher levels of the list
)

keytypes(org.Hs.eg.db)
# lets see how many ensembls id map to a single gene name
sum(duplicated(tcga_expression$Gene))

# Handling duplicate gene identifiers
# As we mentioned earlier, we will not want any duplicate gene identifiers in our data 
# frame when we convert it into a matrix in preparation for running GSVA.
# For RNA-seq processing in refine.bio, transcripts were quantified (Ensembl transcript IDs) 
# and aggregated to the gene-level (Ensembl gene IDs). 

# For a single Entrez ID that maps to multiple Ensembl gene IDs, we will use the
# values associated with the Ensembl gene ID that seems to be most highly expressed. 
# Specifically, we’re going retain the Ensembl gene ID with maximum mean expression value.
# We expect that this approach may be a better reflection of the reads that were quantified
# than taking the mean or median of the values for multiple Ensembl gene IDs would be.

# First, we first need to calculate the gene means, but we’ll need to move our 
# non-numeric variables (the gene ID columns) out of the way for that calculation.
# First let's determine the gene means
gene_means <- rowMeans(tcga_expression %>% dplyr::select(-Gene))

# Let's add this as a column in our `mapped_df`.
tcga_expression <- tcga_expression %>%
  # Add gene_means as a column called gene_means
  dplyr::mutate(gene_means) %>%
  # Reorder the columns so `gene_means` column is upfront
  dplyr::select(Gene, gene_means, dplyr::everything())

# Now we can filter out the duplicate gene identifiers using the gene mean values. 
# First, we’ll use dplyr::arrange() by gene_means such that the the rows will be in order
# of highest gene mean to lowest gene mean. For the duplicate values of gene name,
#the row with the lower index will be the one that’s kept by dplyr::distinct(). 
# In practice, this means that we’ll keep the instance of the Entrez ID with the highest 
# gene mean value as intended.
                               
filtered_tcga_expression <- tcga_expression %>%
  # Sort so that the highest mean expression values are at the top
  dplyr::arrange(dplyr::desc(gene_means)) %>%
  # Filter out the duplicated rows using `dplyr::distinct()`
  dplyr::distinct(Gene, .keep_all = TRUE)

# let’s do our check again to see if we still have duplicates.
sum(duplicated(filtered_tcga_expression$Gene))

# Now we should prep this data so GSVA can use it.
#   remove_rownames  %>% 
filtered_tcga_expression <- filtered_tcga_expression %>% remove_rownames()
filtered_tcga_matrix <- filtered_tcga_expression %>%
  # GSVA can't the Ensembl IDs so we should drop this column as well as the means
  dplyr::select(-gene_means) %>%
  # We need to store our gene identifiers as row names
  tibble::column_to_rownames("Gene") %>%
  # Now we can convert our object into a matrix
  as.matrix()
# Note that if we had duplicate gene identifiers here, we would not be able to set them as row names.

# Gene Set Variation Analysis
# GSVA fits a model and ranks genes based on their expression level relative to the sample distribution
# (Hänzelmann et al. 2013a). The pathway-level score calculated is a way of asking how genes within 
# a gene set vary as compared to genes that are outside of that gene set (Malhotra 2018).

# The idea here is that we will get pathway-level scores for each sample that indicate if genes 
# in a pathway vary concordantly in one direction (over-expressed or under-expressed relative to 
# the overall population) (Hänzelmann et al. 2013a). This means that GSVA scores will depend on 
# the samples included in the dataset when you run GSVA; if you added more samples and ran GSVA again,
# you would expect the scores to change (Hänzelmann et al. 2013a).
# The output is a gene set by sample matrix of GSVA scores.
gsva_results <- gsva(
  filtered_tcga_matrix,
  hallmarks_list,
  method = "gsva",
  # Appropriate for our vst transformed data
  kcdf = "Gaussian",
  # Minimum gene set size
  min.sz = 10,
  # Maximum gene set size
  max.sz = 500,
  # Compute Gaussian-distributed scores
  mx.diff = TRUE,
  # Don't print out the progress bar
  verbose = FALSE
)

# Note that the gsva() function documentation says we can use kcdf = "Gaussian"
# if we have expression values that are continuous such as log-CPMs, log-RPKMs or
# log-TPMs, but we would use kcdf = "Poisson" on integer counts. Our vst() transformed
# data is on a log2-like scale, so Gaussian works for us.
# Print 6 rows,
head(gsva_results[, 1:10])
# let’s write all of our GSVA results to file.
gsva_results %>%
  as.data.frame() %>%
  tibble::rownames_to_column("pathway") %>%
  readr::write_tsv(paste0(DIR2LOAD, "TCGA_gsva_results", genes2consider,".tsv")
  )

pathway_heatmap <- pheatmap::pheatmap(gsva_results,
                                      #annotation_col = annot_df, # Add metadata labels!
                                      show_colnames = FALSE, # Don't show sample labels
                                      fontsize_row = 6 # Shrink the pathway labels a tad
)

# Print out heatmap here
pathway_heatmap

