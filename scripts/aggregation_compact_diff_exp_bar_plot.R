suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(stringi))
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(forcats))


#snakemake = readRDS("snakemake_aggregation_compact_diff_exp_bar_plot.rds")

set.seed(snakemake@params$parameters_props$set_seed)

print("aggregation_compact_diff_exp_bar_plot")


#---------------------------------- Functions ----------------------------------
get_reg_rnas_local = function(reg_rnas) {
  # first, combine all tissues as before, then search for features that have an occurrence of == 1
  # this means, they are only occurring in 1 tissue
  # then, we use %in% to get TRUE/FALSE values if the features in the unique list are in this tissue and sum() over them
  # to do this for all tissues, we use sapply
  unique_reg_rnas = names(table(unlist(reg_rnas))[table(unlist(reg_rnas)) == 1])
  reg_rnas_local = list()
  for (sub_name in names(reg_rnas)) {
    # go through sub lists
    sub_list = reg_rnas[[sub_name]]
    num_reg_local = sum(unique_reg_rnas %in% sub_list)
    reg_rnas_local[[sub_name]] = num_reg_local
  }
  return(unlist(reg_rnas_local))
}

get_reg_rnas_global = function(reg_rnas, num_dereg_tissues) {
  # we combine the list of all tissues knowing that the features are unique in each list with unlist()
  # table() counts the occurence of each feature in the union of all tissues
  # then, we filter with num_dereg_tissues and sum() over all TRUE/FALSE
  return(sum(table(unlist(reg_rnas)) >= num_dereg_tissues))
}

get_reg_rnas_local_names = function(reg_rnas, num_dereg_tissues) {
  # we combine the list of all tissues knowing that the features are unique in each list with unlist()
  # table() counts the occurence of each feature in the union of all tissues
  # then, we filter with num_dereg_tissues return the list of names
  mask = (table(unlist(reg_rnas)) == num_dereg_tissues)
  feature_names = names(table(unlist(reg_rnas)))[mask]
  return(feature_names)
}

get_reg_rnas_global_names = function(reg_rnas, num_dereg_tissues) {
  # we combine the list of all tissues knowing that the features are unique in each list with unlist()
  # table() counts the occurence of each feature in the union of all tissues
  # then, we filter with num_dereg_tissues return the list of names
  mask = table(unlist(reg_rnas)) >= num_dereg_tissues
  feature_names = names(table(unlist(reg_rnas)))[mask]
  return(feature_names)
}

plot_barplot = function(up_reg_rnas_global, up_reg_rnas_local, down_reg_rnas_global, down_reg_rnas_local, tissue, data_input, xticks_names, colors, plot_props, output_folder_fig, output_folder_tab, file_name) {
  
  colors[[tissue]]$all = "#D5A021"
  
  up_reg_rnas_local_df = data.frame(count = up_reg_rnas_local)
  colnames(up_reg_rnas_local_df) = "count"
  up_reg_rnas_local_df$tissue = rownames(up_reg_rnas_local_df)
  rownames(up_reg_rnas_local_df) = c()
  tmp = order(up_reg_rnas_local_df$count)
  up_reg_rnas_local_df = up_reg_rnas_local_df[tmp,]
  up_reg_rnas_local_df$stack = "Unique"
  up_reg_rnas_local_df = rbind(up_reg_rnas_local_df, c(up_reg_rnas_global, "all", "Multiple"))
  up_reg_rnas_local_df$bar_group = 1
  up_reg_rnas_local_df$ordering = 1:length(up_reg_rnas_local_df$count)
  
  down_reg_rnas_local_df = data.frame(count = down_reg_rnas_local)
  colnames(down_reg_rnas_local_df) = "count"
  down_reg_rnas_local_df$tissue = rownames(down_reg_rnas_local_df)
  rownames(down_reg_rnas_local_df) = c()
  tmp = order(down_reg_rnas_local_df$count)
  down_reg_rnas_local_df = down_reg_rnas_local_df[tmp,]
  down_reg_rnas_local_df$stack = "Unique"
  down_reg_rnas_local_df = rbind(down_reg_rnas_local_df, c(down_reg_rnas_global, "all", "Multiple"))
  down_reg_rnas_local_df$bar_group = 2
  down_reg_rnas_local_df$ordering = 1:length(down_reg_rnas_local_df$count)
  
  
  plot_df <- rbind(up_reg_rnas_local_df, down_reg_rnas_local_df) 
  plot_df$bar_group = factor(plot_df$bar_group, levels = c(1,2), labels = c("Upregulated", "Downregulated"))
  plot_df$count = as.numeric(plot_df$count)

  plot_df$stack = factor(plot_df$stack, levels = c("Unique", "Multiple"))
  
  p = ggplot(plot_df, aes(x=stack, y=count, fill=tissue, group = ordering)) + 
    geom_bar(position="stack", stat="identity") +
    #geom_text(aes(label = count), size = 3, position = position_stack(vjust = 0.5)) + 
    facet_grid(~bar_group) +
    labs(x="", y=sprintf("Number of %ss", data_input$rna_class)) +
    scale_fill_manual(values=colors[[tissue]]) +
    theme_classic() +
    theme(legend.position="none",
          legend.title=element_blank(),
          text = element_text(family = plot_props$font_family, size = plot_props$font_size),
          axis.text = element_text(family = plot_props$font_family, size = plot_props$font_size),
          axis.title = element_text(family = plot_props$font_family, size = plot_props$font_size),
          plot.title = element_text(family = plot_props$font_family, size = plot_props$font_size_header),
          legend.text = element_text(family = plot_props$font_family, size = plot_props$font_size),
          # make x invisible
          #axis.text.x=element_blank(),
          panel.spacing = unit(2, "lines"),  # spacing between groups
          strip.text = element_text(family = plot_props$font_family, size = plot_props$font_size),
    )
  #p
  
  # Save results 
  p_adj = p + theme(axis.ticks.x = element_blank())
  ggsave(sprintf("%s/%s.png", output_folder_fig, file_name), p_adj, dpi = plot_props$dpi, width = plot_props$image_width, height = plot_props$image_height, units = plot_props$image_units)
  ggsave(sprintf("%s/%s.svg", output_folder_fig, file_name), p_adj, width = plot_props$image_width, height = plot_props$image_height, units = plot_props$image_units, limitsize = FALSE, scale = 1)
  
  p_flip_right = p + coord_flip() + theme(axis.text.x = element_text(vjust = 0.5, hjust = 1), #angle = -90, 
                                          axis.ticks.y = element_blank())
  height = plot_props$image_height
  ggsave(sprintf("%s/%s_%sx%s_flip_right.png", output_folder_fig, file_name, plot_props$image_width, height), p_flip_right, dpi = plot_props$dpi, width = plot_props$image_width, height = height, units = plot_props$image_units)
  ggsave(sprintf("%s/%s_%sx%s_flip_right.svg", output_folder_fig, file_name, plot_props$image_width, height), p_flip_right, width = plot_props$image_width, height = height, units = plot_props$image_units, limitsize = FALSE, scale = 1)
  height = plot_props$image_height / 2
  ggsave(sprintf("%s/%s_%sx%s_flip_right.png", output_folder_fig, file_name, plot_props$image_width, height), p_flip_right, dpi = plot_props$dpi, width = plot_props$image_width, height = height, units = plot_props$image_units)
  ggsave(sprintf("%s/%s_%sx%s_flip_right.svg", output_folder_fig, file_name, plot_props$image_width, height), p_flip_right, width = plot_props$image_width, height = height, units = plot_props$image_units, limitsize = FALSE, scale = 1)
  
  plot_df$bar_group = fct_rev(plot_df$bar_group)

  p = ggplot(plot_df, aes(x=stack, y=count, fill=tissue, group = ordering)) + 
    geom_bar(position="stack", stat="identity") +
    #geom_text(aes(label = count), size = 3, position = position_stack(vjust = 0.5)) + 
    facet_grid(~bar_group) +
    labs(x="", y=sprintf("Number of %ss", data_input$rna_class)) +
    scale_fill_manual(values=colors[[tissue]]) +
    theme_classic() +
    theme(legend.position="none",
          legend.title=element_blank(),
          text = element_text(family = plot_props$font_family, size = plot_props$font_size),
          axis.text = element_text(family = plot_props$font_family, size = plot_props$font_size),
          axis.title = element_text(family = plot_props$font_family, size = plot_props$font_size),
          plot.title = element_text(family = plot_props$font_family, size = plot_props$font_size_header),
          legend.text = element_text(family = plot_props$font_family, size = plot_props$font_size),
          # make x invisible
          #axis.text.x=element_blank(),
          panel.spacing = unit(2, "lines"),  # spacing between groups
          strip.text = element_text(family = plot_props$font_family, size = plot_props$font_size),
    )
  
  p_flip_left = p + scale_y_reverse () +
    scale_x_discrete(name = "", position = "top") +
    coord_flip () +
    theme(axis.text.x = element_text(vjust = 0.5, hjust = 1), #angle = 90,
          axis.ticks.y = element_blank())
  
  height = plot_props$image_height                              
  ggsave(sprintf("%s/%s_%sx%s_flip_left.png", output_folder_fig, file_name, plot_props$image_width, height), p_flip_left, dpi = plot_props$dpi, width = plot_props$image_width, height = height, units = plot_props$image_units)
  ggsave(sprintf("%s/%s_%sx%s_flip_left.svg", output_folder_fig, file_name, plot_props$image_width, height), p_flip_left, width = plot_props$image_width, height = height, units = plot_props$image_units, limitsize = FALSE, scale = 1)
  height = plot_props$image_height / 2
  ggsave(sprintf("%s/%s_%sx%s_flip_left.png", output_folder_fig, file_name, plot_props$image_width, height), p_flip_left, dpi = plot_props$dpi, width = plot_props$image_width, height = height, units = plot_props$image_units)
  ggsave(sprintf("%s/%s_%sx%s_flip_left.svg", output_folder_fig, file_name, plot_props$image_width, height), p_flip_left, width = plot_props$image_width, height = height, units = plot_props$image_units, limitsize = FALSE, scale = 1)
  
  fwrite(plot_df, sprintf("%s/%s.csv", output_folder_tab, file_name), sep = "\t", row.names = TRUE)
  write.xlsx(plot_df, sprintf("%s/%s.xlsx", output_folder_tab, file_name), colNames = TRUE, rowNames = TRUE, append = FALSE)
}


#---------------------------------- Load data ---------------------------------- 
data_input = snakemake@params$data_input
diff_log = snakemake@params$diff_exp_log
log_fc_thres = log2(snakemake@params$thresholds$fc)
sig_lvl = snakemake@params$thresholds$adj_p_value
plot_props = snakemake@params$plots_props
colors = snakemake@params$colors
xticks_names = snakemake@params$xticks_names
results_folder = snakemake@params$results_folder
number_of_deregulated_combinations_thresholds = snakemake@params$thresholds$num_dereg_combinations_th
number_of_deregulated_tissues_thresholds = snakemake@params$thresholds$num_dereg_tissues_th
props = snakemake@params$props

# save rdata
#saveRDS(snakemake, file = "snakemake_aggregation_compact_diff_exp_bar_plot.rds")

annot = fread(sprintf("%s/annotation_%s_%s.csv", data_input$annotation_data_folder, data_input$data_sub_set, data_input$feature_filtering), sep='\t', colClasses=c(ID="character")) 
diff_exp = fread(sprintf("%s/%s_%s/results_%s/matrices/diff_exp/diff_exp%s.csv", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, diff_log), sep='\t') 

output_folder_fig = sprintf("%s/%s_%s/results_%s/figures/diff_exp/bar_plots", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set)
dir.create(output_folder_fig, recursive=TRUE)
output_folder_tab = sprintf("%s/%s_%s/results_%s/matrices/aggregation_diff_exp", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set)
dir.create(output_folder_tab, recursive=TRUE)


#------------------------------------ Script ----------------------------------- 
for (prop in props) {
  
  tissue = prop[1]
  time = prop[2]

  ####
  up_reg = c()
  down_reg = c()
  sig_up_reg = c()
  sig_down_reg = c()
  ####
  
  up_reg_rnas = c()
  down_reg_rnas = c()
  sig_up_reg_rnas = c()
  sig_down_reg_rnas = c()
  for (num_dereg_combinations in number_of_deregulated_combinations_thresholds){
    for (num_dereg_tissues in number_of_deregulated_tissues_thresholds) {
      #stop()
      for (ti in unique(annot[[tissue]])) {
        up_reg_tissue = c()
        down_reg_tissue = c()
        sig_up_reg_tissue = c()
        sig_down_reg_tissue = c()
        for (time_case in sort(as.numeric(unique(annot[[time]])))[-1]) {
          info_table = data.table(feature=diff_exp[[data_input$rna_class]],
                                  ttest_adj_pval=diff_exp[[sprintf("ttest_adjp_%s__%s_vs_%s_%s=%s", time, time_case, sort(as.numeric(unique(annot[[time]])))[1], tissue, ti)]],  #ttest_adjp_age__12_vs_3_brain_region_fixed=cc
                                  fc=diff_exp[[sprintf("fc_%s__%s_vs_%s_%s=%s", time, time_case, sort(as.numeric(unique(annot[[time]])))[1], tissue, ti)]]) #fc_age__12_vs_3_brain_region_fixed=cc
          info_table[,logfc := log2(fc)]
          
          info_table$fc_cat_up = rep(FALSE, dim(info_table)[1])
          info_table$fc_cat_down = rep(FALSE, dim(info_table)[1])
          info_table$ttest_adj_cat_up = rep(FALSE, dim(info_table)[1])
          info_table$ttest_adj_cat_down = rep(FALSE, dim(info_table)[1])
          
          #info_table[ttest_adj_pval >= sig_lvl, ttest_adj_cat:="not_sig"]
          #info_table[ttest_adj_pval < sig_lvl, ttest_adj_cat:="sig"]
          info_table[logfc >= log_fc_thres, fc_cat_up:=TRUE]
          info_table[logfc <= -log_fc_thres, fc_cat_down:=TRUE]
          info_table[logfc >= log_fc_thres & ttest_adj_pval < sig_lvl, ttest_adj_cat_up:=TRUE]
          info_table[logfc <= -log_fc_thres & ttest_adj_pval < sig_lvl, ttest_adj_cat_down:=TRUE]
          
          up_reg_tissue[[sprintf("%s_case=%s", time, time_case)]] = data.frame(feature = info_table[info_table$fc_cat_up,]$feature)
          down_reg_tissue[[sprintf("%s_case=%s", time, time_case)]] = data.frame(feature = info_table[info_table$fc_cat_down,]$feature)
          sig_up_reg_tissue[[sprintf("%s_case=%s", time, time_case)]] = data.frame(feature = info_table[info_table$ttest_adj_cat_up,]$feature)
          sig_down_reg_tissue[[sprintf("%s_case=%s", time, time_case)]] = data.frame(feature = info_table[info_table$ttest_adj_cat_down,]$feature)
        }
        
        up_reg_tissue_contain = data.table(feature = unique(unlist(up_reg_tissue, use.names = FALSE)))
        down_reg_tissue_contain = data.table(feature = unique(unlist(down_reg_tissue, use.names = FALSE)))
        sig_up_reg_tissue_contain = data.table(feature = unique(unlist(sig_up_reg_tissue, use.names = FALSE)))
        sig_down_reg_tissue_contain = data.table(feature = unique(unlist(sig_down_reg_tissue, use.names = FALSE)))
        
        for (time_case in sort(unique(annot[[time]]))[-1]) {
          up_reg_tissue_contain[[sprintf("%s_case=%s", time, time_case)]] = rep(FALSE, length(up_reg_tissue_contain$feature))
          down_reg_tissue_contain[[sprintf("%s_case=%s", time, time_case)]] = rep(FALSE, length(down_reg_tissue_contain$feature))
          sig_up_reg_tissue_contain[[sprintf("%s_case=%s", time, time_case)]] = rep(FALSE, length(sig_up_reg_tissue_contain$feature))
          sig_down_reg_tissue_contain[[sprintf("%s_case=%s", time, time_case)]] = rep(FALSE, length(sig_down_reg_tissue_contain$feature))
          
          up_reg_tissue_contain[feature %in% up_reg_tissue[[sprintf("%s_case=%s", time, time_case)]]$feature, sprintf("%s_case=%s", time, time_case):=TRUE]
          down_reg_tissue_contain[feature %in% down_reg_tissue[[sprintf("%s_case=%s", time, time_case)]]$feature, sprintf("%s_case=%s", time, time_case):=TRUE]
          sig_up_reg_tissue_contain[feature %in% sig_up_reg_tissue[[sprintf("%s_case=%s", time, time_case)]]$feature, sprintf("%s_case=%s", time, time_case):=TRUE]
          sig_down_reg_tissue_contain[feature %in% sig_down_reg_tissue[[sprintf("%s_case=%s", time, time_case)]]$feature, sprintf("%s_case=%s", time, time_case):=TRUE]
        } 
        
        up_reg_tissue_number = data.table("feature" = up_reg_tissue_contain$feature, "sum" = rowSums(up_reg_tissue_contain[,-"feature"]))
        down_reg_tissue_number = data.table("feature" = down_reg_tissue_contain$feature, "sum" = rowSums(down_reg_tissue_contain[,-"feature"]))
        sig_up_reg_tissue_number = data.table("feature" = sig_up_reg_tissue_contain$feature, "sum" = rowSums(sig_up_reg_tissue_contain[,-"feature"]))
        sig_down_reg_tissue_number = data.table("feature" = sig_down_reg_tissue_contain$feature, "sum" = rowSums(sig_down_reg_tissue_contain[,-"feature"]))
      
        up_reg_tissue_number[sum >= num_dereg_combinations, keep:=TRUE]
        up_reg_tissue_number[sum < num_dereg_combinations, keep:=FALSE]
        down_reg_tissue_number[sum >= num_dereg_combinations, keep:=TRUE]
        down_reg_tissue_number[sum < num_dereg_combinations, keep:=FALSE]
        sig_up_reg_tissue_number[sum >= num_dereg_combinations, keep:=TRUE]
        sig_up_reg_tissue_number[sum < num_dereg_combinations, keep:=FALSE]
        sig_down_reg_tissue_number[sum >= num_dereg_combinations, keep:=TRUE]
        sig_down_reg_tissue_number[sum < num_dereg_combinations, keep:=FALSE]
        
        ####
        up_reg[[ti]] = length(up_reg_tissue_number[up_reg_tissue_number$keep,]$feature)
        down_reg[[ti]] = length(down_reg_tissue_number[down_reg_tissue_number$keep,]$feature)
        sig_up_reg[[ti]] = length(sig_up_reg_tissue_number[sig_up_reg_tissue_number$keep,]$feature)
        sig_down_reg[[ti]] = length(sig_down_reg_tissue_number[sig_down_reg_tissue_number$keep,]$feature)
        ####
        
        up_reg_rnas[[ti]] = up_reg_tissue_number[up_reg_tissue_number$keep,]$feature
        down_reg_rnas[[ti]] = down_reg_tissue_number[down_reg_tissue_number$keep,]$feature
        sig_up_reg_rnas[[ti]] = sig_up_reg_tissue_number[sig_up_reg_tissue_number$keep,]$feature
        sig_down_reg_rnas[[ti]] = sig_down_reg_tissue_number[sig_down_reg_tissue_number$keep,]$feature
      }
      
      # occurrence of feature in the list >= num_dereg_tissues -> number
      up_reg_rnas_global = get_reg_rnas_global(up_reg_rnas, num_dereg_tissues)
      down_reg_rnas_global = get_reg_rnas_global(down_reg_rnas, num_dereg_tissues)
      sig_up_reg_rnas_global = get_reg_rnas_global(sig_up_reg_rnas, num_dereg_tissues)
      sig_down_reg_rnas_global = get_reg_rnas_global(sig_down_reg_rnas, num_dereg_tissues)
      
      # occurrence of feature in the list == 1 -> feature list per tissue
      up_reg_rnas_local = get_reg_rnas_local(up_reg_rnas)
      down_reg_rnas_local = get_reg_rnas_local(down_reg_rnas)
      sig_up_reg_rnas_local = get_reg_rnas_local(sig_up_reg_rnas)
      sig_down_reg_rnas_local = get_reg_rnas_local(sig_down_reg_rnas)
      
      # all features with deregulation
      file_name = sprintf("diff_exp_aggregation_at_least_in_%s_%s_comps_for_%s_%s_complete", num_dereg_combinations, time, num_dereg_tissues, tissue)
      plot_barplot(up_reg_rnas_global, up_reg_rnas_local, down_reg_rnas_global, down_reg_rnas_local, tissue, data_input, xticks_names, colors, plot_props, output_folder_fig, output_folder_tab, file_name)
      # features with a significant deregulation
      file_name = sprintf("diff_exp_sig_aggregation_at_least_in_%s_%s_comps_for_%s_%s_complete", num_dereg_combinations, time, num_dereg_tissues, tissue)
      plot_barplot(sig_up_reg_rnas_global, sig_up_reg_rnas_local, sig_down_reg_rnas_global, sig_down_reg_rnas_local, tissue, data_input, xticks_names, colors, plot_props, output_folder_fig, output_folder_tab, file_name)
      ################
      
      ################
      # calculate a table with features as rows and tissues as columns
      # the values are TRUE/FALSE if the feature-tissue combination is (sig) up/down regulated in at least num_dereg_combinations comparisons
      up_reg_rnas_table = as.data.frame.matrix(table(melt(up_reg_rnas)))
      down_reg_rnas_table = as.data.frame.matrix(table(melt(down_reg_rnas)))
      sig_up_reg_rnas_table = as.data.frame.matrix(table(melt(sig_up_reg_rnas)))
      sig_down_reg_rnas_table = as.data.frame.matrix(table(melt(sig_down_reg_rnas)))

      # filter such that only globas features survive
      # num_dereg_tissues determines if the feature is global
      # get_reg_rnas_global_names
      up_reg_rnas_global_table = up_reg_rnas_table[get_reg_rnas_global_names(up_reg_rnas, num_dereg_tissues), , drop=FALSE]
      down_reg_rnas_global_table = down_reg_rnas_table[get_reg_rnas_global_names(down_reg_rnas, num_dereg_tissues), , drop=FALSE]
      sig_up_reg_rnas_global_table = sig_up_reg_rnas_table[get_reg_rnas_global_names(sig_up_reg_rnas, num_dereg_tissues), , drop=FALSE]
      sig_down_reg_rnas_global_table = sig_down_reg_rnas_table[get_reg_rnas_global_names(sig_down_reg_rnas, num_dereg_tissues), , drop=FALSE]
      
      # TODO not working
      # filter such that only local features survive
      # a feature is local if it occur only in one brain region
      # get_reg_rnas_local_names
      #up_reg_rnas_local_table = up_reg_rnas_table[get_reg_rnas_local_names(up_reg_rnas, 1), , drop=FALSE]
      #down_reg_rnas_local_table = down_reg_rnas_table[get_reg_rnas_local_names(down_reg_rnas, 1), , drop=FALSE]
      #sig_up_reg_rnas_local_table = sig_up_reg_rnas_table[get_reg_rnas_local_names(sig_up_reg_rnas, 1), , drop=FALSE]
      #sig_down_reg_rnas_local_table = sig_down_reg_rnas_table[get_reg_rnas_local_names(sig_down_reg_rnas, 1), , drop=FALSE]
      
      #save tables
      # global up
      file_name = sprintf("diff_exp_global_up_%s_at_least_in_%s_%s_comps_for_%s_%s", data_input$rna_class, num_dereg_combinations, time, num_dereg_tissues, tissue)
      fwrite(up_reg_rnas_global_table, sprintf("%s/%s.csv", output_folder_tab, file_name), sep = "\t", row.names = TRUE)
      write.xlsx(up_reg_rnas_global_table, sprintf("%s/%s.xlsx", output_folder_tab, file_name), colNames = TRUE, rowNames = TRUE, append = FALSE)
      
      # global down
      file_name = sprintf("diff_exp_global_down_%s_at_least_in_%s_%s_comps_for_%s_%s", data_input$rna_class, num_dereg_combinations, time, num_dereg_tissues, tissue)
      fwrite(down_reg_rnas_global_table, sprintf("%s/%s.csv", output_folder_tab, file_name), sep = "\t", row.names = TRUE)
      write.xlsx(down_reg_rnas_global_table, sprintf("%s/%s.xlsx", output_folder_tab, file_name), colNames = TRUE, rowNames = TRUE, append = FALSE)
      
      # global sig up
      file_name = sprintf("diff_exp_global_sig_up_%s_at_least_in_%s_%s_comps_for_%s_%s", data_input$rna_class, num_dereg_combinations, time, num_dereg_tissues, tissue)
      fwrite(sig_up_reg_rnas_global_table, sprintf("%s/%s.csv", output_folder_tab, file_name), sep = "\t", row.names = TRUE)
      write.xlsx(sig_up_reg_rnas_global_table, sprintf("%s/%s.xlsx", output_folder_tab, file_name), colNames = TRUE, rowNames = TRUE, append = FALSE)
      
      # global sig down
      file_name = sprintf("diff_exp_global_sig_down_%s_at_least_in_%s_%s_comps_for_%s_%s", data_input$rna_class, num_dereg_combinations, time, num_dereg_tissues, tissue)
      fwrite(sig_down_reg_rnas_global_table, sprintf("%s/%s.csv", output_folder_tab, file_name), sep = "\t", row.names = TRUE)
      write.xlsx(sig_down_reg_rnas_global_table, sprintf("%s/%s.xlsx", output_folder_tab, file_name), colNames = TRUE, rowNames = TRUE, append = FALSE)
    
      # up
      file_name = sprintf("diff_exp_up_%s_at_least_in_%s_%s_comps_for_%s_%s", data_input$rna_class, num_dereg_combinations, time, num_dereg_tissues, tissue)
      fwrite(up_reg_rnas_table, sprintf("%s/%s.csv", output_folder_tab, file_name), sep = "\t", row.names = TRUE)
      write.xlsx(up_reg_rnas_table, sprintf("%s/%s.xlsx", output_folder_tab, file_name), colNames = TRUE, rowNames = TRUE, append = FALSE)
      
      # down
      file_name = sprintf("diff_exp_down_%s_at_least_in_%s_%s_comps_for_%s_%s", data_input$rna_class, num_dereg_combinations, time, num_dereg_tissues, tissue)
      fwrite(down_reg_rnas_table, sprintf("%s/%s.csv", output_folder_tab, file_name), sep = "\t", row.names = TRUE)
      write.xlsx(down_reg_rnas_table, sprintf("%s/%s.xlsx", output_folder_tab, file_name), colNames = TRUE, rowNames = TRUE, append = FALSE)
      
      # sig up
      file_name = sprintf("diff_exp_sig_up_%s_at_least_in_%s_%s_comps_for_%s_%s", data_input$rna_class, num_dereg_combinations, time, num_dereg_tissues, tissue)
      fwrite(sig_up_reg_rnas_table, sprintf("%s/%s.csv", output_folder_tab, file_name), sep = "\t", row.names = TRUE)
      write.xlsx(sig_up_reg_rnas_table, sprintf("%s/%s.xlsx", output_folder_tab, file_name), colNames = TRUE, rowNames = TRUE, append = FALSE)
      
      # sig down
      file_name = sprintf("diff_exp_sig_down_%s_at_least_in_%s_%s_comps_for_%s_%s", data_input$rna_class, num_dereg_combinations, time, num_dereg_tissues, tissue)
      fwrite(sig_down_reg_rnas_table, sprintf("%s/%s.csv", output_folder_tab, file_name), sep = "\t", row.names = TRUE)
      write.xlsx(sig_down_reg_rnas_table, sprintf("%s/%s.xlsx", output_folder_tab, file_name), colNames = TRUE, rowNames = TRUE, append = FALSE)
    }
  }
}


#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])
