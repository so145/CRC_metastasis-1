source('/data/home/hfx941/code/CRC/st/Visium_Ozato_2023/A8_cellphonedb/A3_plots/cpdb_heatmap_microenvs.R')
library("optparse")

option_list = list(
  make_option(c("-m", "--microenvs"), type="character", default=NULL, 
              help="microenv file (full path)", metavar="file"),
  make_option(c("-f", "--pvalues"), type="character", default=NULL, 
              help="pvalues file (full path)", metavar="file"),
  make_option(c("-o", "--out"), type="character", default=NULL, 
              help="output directory path", metavar="directory"),
  make_option(c("-p", "--pvalue"), type="double", default=0.05, 
              help="pvalue threshold", metavar="double")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

print(opt)

heatmaps_plot(microenvs_file = opt$microenvs,
              pvalues_file= opt$pvalues,
              output_dir =  opt$out,
              #count_filename='cellphonedb_heatmap_count_test_fct.pdf', 
              #log_filename='cellphonedb_heatmap_logcount_test_fct.pdf', 
              #count_network_filename='count_network.txt', 
              #interaction_count_filename='interactions_count.txt', 
              #count_network_separator='\t', 
              #interaction_count_separator='\t', 
              show_rownames = T, 
              show_colnames = T, scale="none", cluster_cols = T,
              border_color='white', cluster_rows = T, fontsize_row=11,
              fontsize_col = 11, main = '', treeheight_row=0, 
              family='Arial', treeheight_col = 0, col1 = "dodgerblue4", 
              col2 = 'peachpuff', col3 = 'deeppink4', meta_sep='\t', 
              pvalues_sep='\t', pvalue=opt$pvalue)

