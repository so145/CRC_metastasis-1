library(pheatmap)

heatmaps_plot = function(microenvs_file, 
                         pvalues_file,
                         output_dir,
                         show_rownames = T, 
                         show_colnames = T, scale="none", cluster_cols = T,
                         border_color='white', cluster_rows = T, fontsize_row=11,
                         fontsize_col = 11, main = '', treeheight_row=0, 
                         family='Arial', treeheight_col = 0, col1 = "dodgerblue4", 
                         col2 = 'peachpuff', col3 = 'deeppink4', meta_sep='\t', 
                         pvalues_sep='\t', pvalue=0.05) {
  
  ### Network 
  
  # set to microenvironment we want to look at 
  #micro2look = '8'
  
  # read pvalues file
  all_intr = read.table(pvalues_file, header=T, stringsAsFactors = F, 
                        sep=pvalues_sep, comment.char = '', check.names = F)
  # select interaction LR 
  intr_pairs = all_intr$interacting_pair
  all_intr = all_intr[,-c(1:11)] # only select columns that have cell-type interactions
  
  # read microenvironment file 
  microenv = read.table(microenvs_file, header=T, stringsAsFactors = F, 
                        sep=pvalues_sep, comment.char = '', check.names = F)
  
  
  # generate data and heatmaps for each microenvironment
  for (micro2look in unique(microenv[,2])) {
  
  # select cell types that are only present in microenv 'micro2look'
  microenv2plot = microenv[microenv$microenviroment==micro2look,]
  
  # start constructing network: nb interactions between pairs of cell types 
  pairs1_all = unique(microenv2plot[,1])
  
  split_sep = '\\|'
  join_sep = '|'
  
  pairs1 = c()
  for (i in 1:length(pairs1_all))
    for (j in 1:length(pairs1_all))
      pairs1 = c(pairs1,paste(pairs1_all[i],pairs1_all[j],sep=join_sep))
  
  all_count = matrix(ncol=3)
  colnames(all_count) = c('SOURCE','TARGET','count')
  count1 = c()
  for(i in 1:length(pairs1))
  {
    p1 = strsplit(pairs1[i], split_sep)[[1]][1]
    p2 = strsplit(pairs1[i], split_sep)[[1]][2]
    
    n1 = intr_pairs[which(all_intr[,pairs1[i]]<=pvalue)]
    pairs_rev = paste(p2, p1, sep=join_sep)
    n2 = intr_pairs[which(all_intr[,pairs_rev]<=pvalue)]
    
    if(p1!=p2)
      count1 = length(unique(n1))+length(unique(n2))
    else
      count1 = length(unique(n1))
    
    new_count = c(p1,p2,count1)
    names(new_count) = c('SOURCE','TARGET','count')
    all_count = rbind(all_count, new_count)
  }
  
  all_count = all_count[-1,]
  count_network_filename = paste('microenv',micro2look,'count_network.txt',sep = '_')
  count_network_filename = paste0(output_dir, count_network_filename)
  count_network_separator = '\t'
  write.table(all_count, count_network_filename, sep=count_network_separator, 
              quote=F, row.names = F)
  
  #######   count interactions
  
  count1 = c()
  for(i in 1:length(pairs1))
  {
    p1 = strsplit(pairs1[i], split_sep)[[1]][1]
    p2 = strsplit(pairs1[i], split_sep)[[1]][2]
    
    n1 = intr_pairs[which(all_intr[,pairs1[i]]<=pvalue)]
    
    pairs_rev = paste(p2, p1, sep=join_sep)
    n2 = intr_pairs[which(all_intr[,pairs_rev]<=pvalue)]
    if(p1!=p2)
      count1 = c(count1,length(unique(n1))+length(unique(n2)))
    else
      count1 = c(count1,length(unique(n1)))
    
  }
  
  if (any(count1)>0)
  {
    count_matrix = matrix(count1, nrow=length(unique(microenv2plot[,1])), 
                          ncol=length(unique(microenv2plot[,1])))
    rownames(count_matrix)= unique(microenv2plot[,1])
    colnames(count_matrix)= unique(microenv2plot[,1])
    
    all_sum = rowSums(count_matrix)
    all_sum = cbind(names(all_sum), all_sum)
    interaction_count_filename = paste('microenv',micro2look,'interactions_count.txt',sep = '_')
    interaction_count_filename = paste0(output_dir, interaction_count_filename)
    count_network_separator = '\t'
    write.table(all_sum, file=interaction_count_filename, quote=F, 
                sep=count_network_separator, row.names=F)
    
    col.heatmap <- colorRampPalette(c(col1,col2,col3 ))( 1000 )
    
    count_filename = paste('microenv',micro2look,'heatmap_count.pdf',sep = '_')
    count_filename = paste0(output_dir, count_filename)
    pheatmap(count_matrix, show_rownames = show_rownames, show_colnames = show_colnames, scale=scale, cluster_cols = cluster_cols,
             border_color=border_color, cluster_rows = cluster_rows, fontsize_row = fontsize_row, fontsize_col = fontsize_col,
             main = main, treeheight_row = treeheight_row, family = family,color = col.heatmap, treeheight_col = treeheight_col, 
             filename = count_filename)
    
    log_filename = paste('microenv',micro2look,'heatmap_log_count.pdf',sep = '_')
    log_filename = paste0(output_dir, log_filename)
    pheatmap(log(count_matrix+1), show_rownames = show_rownames, show_colnames = show_colnames, scale=scale, cluster_cols = cluster_cols,
             border_color=border_color, cluster_rows = cluster_rows, fontsize_row = fontsize_row, fontsize_col = fontsize_col,
             main = main, treeheight_row = treeheight_row, family = family,color = col.heatmap, treeheight_col = treeheight_col, 
             filename = log_filename)
  } else {
    stop("There are no significant results using p-value of: ", pvalue, call.=FALSE)
  }
  
  }
}