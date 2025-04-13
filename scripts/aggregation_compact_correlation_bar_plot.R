suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(stringi))
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(forcats))

#snakemake = readRDS("snakemake_aggregation_compact_correlation_bar_plot.rds")

set.seed(snakemake@params$parameters_props$set_seed)

print("aggregation_compact_correlation_bar_plot")


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
  
  plot_df$bar_group = factor(plot_df$bar_group, levels = c(1,2), labels = c("Pos. correlated", "Neg. correlated"))
  plot_df$count = as.numeric(plot_df$count)
  
  plot_df$stack = factor(plot_df$stack, levels = c("Unique", "Multiple"))
  
  p = ggplot(plot_df, aes(x=stack, y=count, fill=tissue, group = ordering)) + 
    geom_bar(position="stack", stat="identity") +
    #geom_text(aes(label = count), size = 3, position = position_stack(vjust = 0.5)) +
    facet_grid(~bar_group) +
    labs(x="", y=sprintf("Number of %s", data_input$rna_class)) +
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
  
  # Save results 
  p_adj = p + theme(axis.ticks.x = element_blank())
  ggsave(sprintf("%s/%s.png", output_folder_fig, file_name), p_adj, dpi = plot_props$dpi, width = plot_props$image_width, height = plot_props$image_height, units = plot_props$image_units)
  ggsave(sprintf("%s/%s.svg", output_folder_fig, file_name), p_adj, width = plot_props$image_width, height = plot_props$image_height, units = plot_props$image_units, limitsize = FALSE, scale = 1)
  
  p_flip_right = p + coord_flip() + theme(axis.text.x = element_text(vjust = 0.5, hjust = 1), 
                                          axis.ticks.y = element_blank()) #angle = -90, 
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
    labs(x="", y=sprintf("Number of %s", data_input$rna_class)) +
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
                    theme(axis.text.x = element_text(vjust = 0.5, hjust = 1),
                          axis.ticks.y = element_blank()) #angle = 90,
         
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
#tissues = snakemake@params$tissues 
plot_props = snakemake@params$plots_props
colors = snakemake@params$colors
xticks_names = snakemake@params$xticks_names
corr_list = snakemake@params$thresholds$corr_th
number_of_corr_tissues_thresholds = snakemake@params$thresholds$num_corr_tissues_th
sig_lvl = snakemake@params$thresholds$adj_p_value
method = snakemake@params$method
props = snakemake@params$props
results_folder = snakemake@params$results_folder

# save rdata
#saveRDS(snakemake, file = "snakemake_aggregation_compact_correlation_bar_plot.rds")
#stop()

output_folder_fig = sprintf("%s/%s_%s/results_%s/figures/corr_plots/bar_plots",results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set)
dir.create(output_folder_fig, recursive=TRUE)
output_folder_tab = sprintf("%s/%s_%s/results_%s/matrices/aggregation_corr", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set)
dir.create(output_folder_tab, recursive=TRUE)
output_folder_tab_mieaa = sprintf("%s/%s_%s/results_%s/matrices/mieaa_input", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set)
dir.create(output_folder_tab_mieaa, recursive=TRUE)

annot = fread(sprintf("%s/annotation_%s_%s.csv", data_input$annotation_data_folder, data_input$data_sub_set, data_input$feature_filtering), sep='\t', colClasses=c(ID="character")) 


#------------------------------------ Script ----------------------------------- 
for (m in method) {
  for (prop in props) {
    
    tissue = prop[1]
    time = prop[2]
     
    # load data
    corr = fread(sprintf("%s/corr_method=%s_%ss_with_%s_per_%s.csv", output_folder_tab, m, data_input$rna_class, time, tissue), sep='\t')
    corr$V1 = c()
    pvalue = fread(sprintf("%s/corr_method=%s_adj_pvalues_%ss_with_%s_per_%s.csv", output_folder_tab, m, data_input$rna_class, time, tissue), sep='\t')
    pvalue$V1 = c()  
    #pvalue_raw = fread(sprintf("%s/corr_method=%s_pvalues_%ss_with_%s_per_%s.csv", output_folder_tab, m, data_input$rna_class, time, tissue), sep='\t')
    #pvalue_raw$V1 = c()
    
    for (th in corr_list){
      for (num_corr_tissues in number_of_corr_tissues_thresholds) {
        
        ####
        pos_corr_number = c()
        neg_corr_number = c()
        sig_pos_corr_number = c()
        sig_neg_corr_number = c()
        ####
        
        pos_corr_rnas = c()
        neg_corr_rnas = c()
        sig_pos_corr_rnas = c()
        sig_neg_corr_rnas = c()
        #info_table_tissue_tmp = c()
        for (ti in unique(annot[[tissue]])) {
          info_table_tissue = c()
          for (miR in corr[[data_input$rna_class]]) {
            info_table = data.table(feature=miR,
                                    tissue=ti,
                                    ttest_adj_pval=pvalue[pvalue[[data_input$rna_class]] == miR,][[ti]], 
                                    #ttest_pval=pvalue_raw[pvalue_raw[[data_input$rna_class]] == miR,][[ti]], 
                                    corr=corr[corr[[data_input$rna_class]] == miR,][[ti]]
            )
            
            info_table$corr_cat_pos = rep(FALSE, dim(info_table)[1])
            info_table$corr_cat_neg = rep(FALSE, dim(info_table)[1])
            info_table$ttest_adj_cat_corr_pos = rep(FALSE, dim(info_table)[1])
            info_table$ttest_adj_cat_corr_neg = rep(FALSE, dim(info_table)[1])
            
            #info_table[ttest_adj_pval >= sig_lvl, ttest_adj_cat:="not_sig"]
            #info_table[ttest_adj_pval < sig_lvl, ttest_adj_cat:="sig"]
            info_table[corr >= th, corr_cat_pos:=TRUE]
            info_table[corr <= -th, corr_cat_neg:=TRUE]
            info_table[corr >= th & ttest_adj_pval < sig_lvl, ttest_adj_cat_corr_pos:=TRUE]
            info_table[corr <= -th & ttest_adj_pval < sig_lvl, ttest_adj_cat_corr_neg:=TRUE]
            
            info_table_tissue = rbind(info_table_tissue, info_table)
          }
          
          #tmp = log10(info_table_tissue$ttest_pval) * sign(info_table_tissue$corr)
          #info_table_tissue_mieaa_input = info_table_tissue[order(tmp)]$feature
          
          #info_table_tissue_mieaa_input_df = data.frame(info_table_tissue_mieaa_input)
          #fwrite(info_table_tissue_mieaa_input_df, sprintf("%s/corr_method=%s_corr_th=%s_pvalues_%ss_with_%s_per_%s=%s.csv", output_folder_tab_mieaa, m, th, data_input$rna_class, time, tissue, ti), sep = "\t", row.names = FALSE, col.names = FALSE)
          #write.xlsx(info_table_tissue_mieaa_input_df, sprintf("%s/corr_method=%s_corr_th=%s_pvalues_%ss_with_%s_per_%s=%s.xlsx", output_folder_tab_mieaa, m, th, data_input$rna_class, time, tissue, ti), colNames = FALSE, rowNames = FALSE, append = FALSE)
          
          #info_table_tissue_tmp = rbind(info_table_tissue_tmp, info_table_tissue[info_table_tissue$corr == max(info_table_tissue$corr),])
          
          pos_corr = data.frame("feature" = info_table_tissue[info_table_tissue$corr_cat_pos,]$feature)
          neg_corr = data.frame("feature" = info_table_tissue[info_table_tissue$corr_cat_neg,]$feature)
          sig_pos_corr = data.frame("feature" = info_table_tissue[info_table_tissue$ttest_adj_cat_corr_pos,]$feature)
          sig_neg_corr = data.frame("feature" = info_table_tissue[info_table_tissue$ttest_adj_cat_corr_neg,]$feature)
          
          pos_corr_contain = data.table("feature" = unique(unlist(pos_corr, use.names = FALSE)))
          neg_corr_contain = data.table("feature" = unique(unlist(neg_corr, use.names = FALSE)))
          sig_pos_corr_contain = data.table("feature" = unique(unlist(sig_pos_corr, use.names = FALSE)))
          sig_neg_corr_contain = data.table("feature" = unique(unlist(sig_neg_corr, use.names = FALSE)))
          
          pos_corr_number[[ti]] = length(pos_corr_contain$feature)
          neg_corr_number[[ti]] = length(neg_corr_contain$feature)
          sig_pos_corr_number[[ti]] = length(sig_pos_corr_contain$feature)
          sig_neg_corr_number[[ti]] = length(sig_neg_corr_contain$feature)
          
          pos_corr_rnas[[ti]] = pos_corr_contain$feature
          neg_corr_rnas[[ti]] = neg_corr_contain$feature
          sig_pos_corr_rnas[[ti]] = sig_pos_corr_contain$feature
          sig_neg_corr_rnas[[ti]] = sig_neg_corr_contain$feature
        }
        
        # occurrence of feature in the list >= num_corr_tissues -> number
        pos_corr_rnas_global = get_reg_rnas_global(pos_corr_rnas, num_corr_tissues)
        neg_corr_rnas_global = get_reg_rnas_global(neg_corr_rnas, num_corr_tissues)
        sig_pos_corr_rnas_global = get_reg_rnas_global(sig_pos_corr_rnas, num_corr_tissues)
        sig_neg_corr_rnas_global = get_reg_rnas_global(sig_neg_corr_rnas, num_corr_tissues)
        
        # occurrence of feature in the list == 1 -> feature list per tissue
        pos_corr_rnas_local = get_reg_rnas_local(pos_corr_rnas)
        neg_corr_rnas_local = get_reg_rnas_local(neg_corr_rnas)
        sig_pos_corr_rnas_local = get_reg_rnas_local(sig_pos_corr_rnas)
        sig_neg_corr_rnas_local = get_reg_rnas_local(sig_neg_corr_rnas)
        
        # all features with deregulation
        file_name = sprintf("corr_method=%s_th=%s_aggregation_for_%s_in_%s_complete", m, th, num_corr_tissues, tissue)
        plot_barplot(pos_corr_rnas_global, pos_corr_rnas_local, neg_corr_rnas_global, neg_corr_rnas_local, tissue, data_input, xticks_names, colors, plot_props, output_folder_fig, output_folder_tab, file_name)
        # features with a significant deregulation
        file_name = sprintf("corr_method=%s_th=%s_sig_aggregation_for_%s_in_%s_complete", m, th, num_corr_tissues, tissue)
        plot_barplot(sig_pos_corr_rnas_global, sig_pos_corr_rnas_local, sig_neg_corr_rnas_global, sig_neg_corr_rnas_local, tissue, data_input, xticks_names, colors, plot_props, output_folder_fig, output_folder_tab, file_name)
        ################
        
        ################
        # calculate a table with features as rows and tissues as columns
        # the values are TRUE/FALSE if the feature-tissue combination is pos/neg correlatad in at least num_dereg_combinations comparisons
        pos_rnas_table = as.data.frame.matrix(table(melt(pos_corr_rnas)))
        neg_rnas_table = as.data.frame.matrix(table(melt(neg_corr_rnas)))
        sig_pos_rnas_table = as.data.frame.matrix(table(melt(sig_pos_corr_rnas)))
        sig_neg_rnas_table = as.data.frame.matrix(table(melt(sig_neg_corr_rnas)))

        # filter such that only global features survive
        # num_dereg_tissues determines if the feature is global
        # get_reg_rnas_global_names gets the names of the globas rna
        pos_rnas_global_table = pos_rnas_table[get_reg_rnas_global_names(pos_corr_rnas, num_corr_tissues), ]
        neg_rnas_global_table = neg_rnas_table[get_reg_rnas_global_names(neg_corr_rnas, num_corr_tissues), ]
        sig_pos_rnas_global_table = sig_pos_rnas_table[get_reg_rnas_global_names(sig_pos_corr_rnas, num_corr_tissues), ]
        sig_neg_rnas_global_table = sig_neg_rnas_table[get_reg_rnas_global_names(sig_neg_corr_rnas, num_corr_tissues), ]
        
        # TODO not working
        # filter such that only local features survive
        # a feature is local if it occur only in one tissue
        # get_reg_rnas_local_names
        #pos_rnas_local_table = pos_rnas_table[get_reg_rnas_local_names(pos_corr_rnas, 1), ]
        #neg_rnas_local_table = neg_rnas_table[get_reg_rnas_local_names(neg_corr_rnas, 1), ]
        #sig_pos_rnas_local_table = pos_rnas_table[get_reg_rnas_local_names(sig_pos_corr_rnas, 1), ]
        #sig_neg_rnas_local_table = neg_rnas_table[get_reg_rnas_local_names(sig_neg_corr_rnas, 1), ]
        
        # save tables
        # global pos
        file_name = sprintf("corr_method=%s_global_pos_for_%s_%s_%s_corr_th=%s", m, num_corr_tissues, tissue, data_input$rna_class, th)
        if (class(pos_rnas_global_table) == "data.frame") {
          fwrite(pos_rnas_global_table, sprintf("%s/%s.csv", output_folder_tab, file_name), sep = "\t", row.names = TRUE)
          write.xlsx(pos_rnas_global_table, sprintf("%s/%s.xlsx", output_folder_tab, file_name), colNames = TRUE, rowNames = TRUE, append = FALSE)
        } else {
          print(sprintf("%s is empty", file_name))
        }
                
        # global neg
        file_name = sprintf("corr_method=%s_global_neg_for_%s_%s_%s_corr_th=%s", m, num_corr_tissues, tissue, data_input$rna_class, th)
        if (class(neg_rnas_global_table) == "data.frame") {
          fwrite(neg_rnas_global_table, sprintf("%s/%s.csv", output_folder_tab, file_name), sep = "\t", row.names = TRUE)
          write.xlsx(neg_rnas_global_table, sprintf("%s/%s.xlsx", output_folder_tab, file_name), colNames = TRUE, rowNames = TRUE, append = FALSE)
        } else {
          print(sprintf("%s is empty", file_name))
        }
          
        # global sig pos
        file_name = sprintf("corr_method=%s_global_sig_pos_for_%s_%s_%s_corr_th=%s", m, num_corr_tissues, tissue, data_input$rna_class, th)
        if (class(sig_pos_rnas_global_table) == "data.frame") {
          fwrite(sig_pos_rnas_global_table, sprintf("%s/%s.csv", output_folder_tab, file_name), sep = "\t", row.names = TRUE)
          write.xlsx(sig_pos_rnas_global_table, sprintf("%s/%s.xlsx", output_folder_tab, file_name), colNames = TRUE, rowNames = TRUE, append = FALSE)
        } else {
          print(sprintf("%s is empty", file_name))
        }
        
        # global sig neg
        file_name = sprintf("corr_method=%s_global_sig_neg_for_%s_%s_%s_corr_th=%s", m, num_corr_tissues, tissue, data_input$rna_class, th)
        if (class(sig_neg_rnas_global_table) == "data.frame") {
          fwrite(sig_neg_rnas_global_table, sprintf("%s/%s.csv", output_folder_tab, file_name), sep = "\t", row.names = TRUE)
          write.xlsx(sig_neg_rnas_global_table, sprintf("%s/%s.xlsx", output_folder_tab, file_name), colNames = TRUE, rowNames = TRUE, append = FALSE)
        } else {
          print(sprintf("%s is empty", file_name))
        }
        
        # pos
        file_name = sprintf("corr_method=%s_pos_for_%s_%s_%s_corr_th=%s", m, num_corr_tissues, tissue, data_input$rna_class, th)
        if (class(pos_rnas_table) == "data.frame") {
          fwrite(pos_rnas_table, sprintf("%s/%s.csv", output_folder_tab, file_name), sep = "\t", row.names = TRUE)
          write.xlsx(pos_rnas_table, sprintf("%s/%s.xlsx", output_folder_tab, file_name), colNames = TRUE, rowNames = TRUE, append = FALSE)
        } else {
          print(sprintf("%s is empty", file_name))
        }
        
        # neg
        file_name = sprintf("corr_method=%s_neg_for_%s_%s_%s_corr_th=%s", m, num_corr_tissues, tissue, data_input$rna_class, th)
        if (class(neg_rnas_table) == "data.frame") {
          fwrite(neg_rnas_table, sprintf("%s/%s.csv", output_folder_tab, file_name), sep = "\t", row.names = TRUE)
          write.xlsx(neg_rnas_table, sprintf("%s/%s.xlsx", output_folder_tab, file_name), colNames = TRUE, rowNames = TRUE, append = FALSE)
        } else {
          print(sprintf("%s is empty", file_name))
        }
        
        # sig pos
        file_name = sprintf("corr_method=%s_sig_pos_for_%s_%s_%s_corr_th=%s", m, num_corr_tissues, tissue, data_input$rna_class, th)
        if (class(sig_pos_rnas_table) == "data.frame") {
          fwrite(sig_pos_rnas_table, sprintf("%s/%s.csv", output_folder_tab, file_name), sep = "\t", row.names = TRUE)
          write.xlsx(sig_pos_rnas_table, sprintf("%s/%s.xlsx", output_folder_tab, file_name), colNames = TRUE, rowNames = TRUE, append = FALSE)
        } else {
          print(sprintf("%s is empty", file_name))
        }
        
        # sig neg
        file_name = sprintf("corr_method=%s_sig_neg_for_%s_%s_%s_corr_th=%s", m, num_corr_tissues, tissue, data_input$rna_class, th)
        if (class(sig_neg_rnas_table) == "data.frame") {
          fwrite(sig_neg_rnas_table, sprintf("%s/%s.csv", output_folder_tab, file_name), sep = "\t", row.names = TRUE)
          write.xlsx(sig_neg_rnas_table, sprintf("%s/%s.xlsx", output_folder_tab, file_name), colNames = TRUE, rowNames = TRUE, append = FALSE)
        } else {
          print(sprintf("%s is empty", file_name))
        }
      }
    }
  }
}


#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])
