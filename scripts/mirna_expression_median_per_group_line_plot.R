suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(svglite))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(openxlsx))

#snakemake = readRDS("snakemake_mirna_expression_median_per_group_line_plot.rds")

data.table::setDTthreads(1)
set.seed(snakemake@params$parameters_props$set_seed)

print("mirna_expression_median_per_group_line_plot")


#---------------------------------- Functions ----------------------------------


line_plot = function(plot_df, colour_set, prop, data_input, xticks_names, plot_props, output_folder, file_name) {
  
  if (grepl(" ", xticks_names$categories[[data_input$detection_group]]) == TRUE) {
    detection_group = strsplit(xticks_names$categories[[data_input$detection_group]], " ")[[1]][1]
  }
  #y_axis_name = sprintf("Expr (%s,\ndetection rate = %s%%\nper %s)", gsub("_norm", "", data_input$norm), gsub("p", "", data_input$detection_rate), detection_group)
  y_axis_name = ""
  
  plot_df$feature_visible = gsub("mmu-", "", plot_df$feature)
  mask = (plot_df$prop == "young")
  plot_df[mask,]$feature_visible = ""
  
  if (unique(unique(plot_df$prop) %in% c("young", "old")) == TRUE) {
    plot_df$prop = factor(plot_df$prop, levels = c("young","old"),ordered = TRUE)
  }
  
  font_size_correction = 6/17.07
  
  # no log
  p = ggplot(plot_df, aes(x = prop, y = value, group = feature)) + 
    geom_line(aes(color=feature), lwd=1.5) +
    geom_text_repel(aes(label = feature_visible), color = "black", size = plot_props$font_size*font_size_correction, direction="y", nudge_x = 0.5, hjust = 0, force = 0.5, label.padding = 0.3, box.padding = 0.3, min.segment.length = 0) + #segment.color = NA, segment.size = 0
    scale_x_discrete(expand = c(0, 0)) +
    scale_color_manual(values=colour_set) +
    xlab("") + 
    ylab(y_axis_name) +
    theme_classic() +
    theme(plot.margin = margin(-0.2, 0, -0.35, 0, "cm"), #top, right, bottom, and left left=-1.95
          legend.position="none", 
          #legend.margin = margin(t = -10, r = 5, b = -5, l = -5),
          line = element_blank(),
          text = element_text(family = plot_props$font_family, size = plot_props$font_size),
          axis.text = element_text(family = plot_props$font_family, size = plot_props$font_size),
          axis.title = element_text(family = plot_props$font_family, size = plot_props$font_size),
          #axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.title = element_text(family = plot_props$font_family, size = plot_props$font_size_header),
          legend.text = element_text(family = plot_props$font_family, size = plot_props$font_size),
          legend.title = element_text(family = plot_props$font_family, size = plot_props$font_size)
    )
  p
  
  # save bar plot
  ggsave(sprintf("%s/%s_quad.png", output_folder, file_name), p, dpi=plot_props$dpi, width=2/3*plot_props$image_width, height=plot_props$image_height, units = plot_props$image_units)
  ggsave(sprintf("%s/%s_quad.svg", output_folder, file_name), p, width=2/3*plot_props$image_width, height=plot_props$image_height, units =plot_props$image_units, limitsize = FALSE, scale = 1)
  
  ggsave(sprintf("%s/%s.png", output_folder, file_name), p, dpi=plot_props$dpi, width=plot_props$image_width, height=plot_props$image_height, units = plot_props$image_units)
  ggsave(sprintf("%s/%s.svg", output_folder, file_name), p, width=plot_props$image_width, height=plot_props$image_height, units =plot_props$image_units, limitsize = FALSE, scale = 1)
  
  # log10
  p = ggplot(plot_df, aes(x = prop, y = value_log, group = feature)) + 
    geom_line(aes(color=feature), lwd=1.5) +
    geom_text_repel(aes(label = feature_visible), color = "black", size = plot_props$font_size*font_size_correction, direction="y", nudge_x = 0.5, hjust = 0, force = 0.5, label.padding = 0.3, box.padding = 0.3, min.segment.length = 0) + #segment.color = NA, segment.size = 0
    scale_x_discrete(expand = c(0, 0)) +
    scale_color_manual(values=colour_set) +
    xlab("") + 
    ylab(y_axis_name) +
    theme_classic() +
    theme(plot.margin = margin(-0.2, 0, -0.35, 0, "cm"), #top, right, bottom, and left left=-1.95
          legend.position="none", 
          #legend.margin = margin(t = -10, r = 5, b = -5, l = -5),
          line = element_blank(),
          text = element_text(family = plot_props$font_family, size = plot_props$font_size),
          axis.text = element_text(family = plot_props$font_family, size = plot_props$font_size),
          axis.title = element_text(family = plot_props$font_family, size = plot_props$font_size),
          #axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.title = element_text(family = plot_props$font_family, size = plot_props$font_size_header),
          legend.text = element_text(family = plot_props$font_family, size = plot_props$font_size),
          legend.title = element_text(family = plot_props$font_family, size = plot_props$font_size)
    )
  p
  
  # save bar plot
  ggsave(sprintf("%s/%s_log10_quad.png", output_folder, file_name), p, dpi=plot_props$dpi, width=2/3*plot_props$image_width, height=plot_props$image_height, units = plot_props$image_units)
  ggsave(sprintf("%s/%s_log10_quad.svg", output_folder, file_name), p, width=2/3*plot_props$image_width, height=plot_props$image_height, units =plot_props$image_units, limitsize = FALSE, scale = 1)
  
  ggsave(sprintf("%s/%s_log10.png", output_folder, file_name), p, dpi=plot_props$dpi, width=plot_props$image_width, height=plot_props$image_height, units = plot_props$image_units)
  ggsave(sprintf("%s/%s_log10.svg", output_folder, file_name), p, width=plot_props$image_width, height=plot_props$image_height, units =plot_props$image_units, limitsize = FALSE, scale = 1)
}


#---------------------------------- Load data ----------------------------------
data_input = snakemake@params$data_input
diff_exp_log = snakemake@params$diff_log
properties = snakemake@params$properties
features_black = snakemake@params$features_black
features_grey = snakemake@params$features_grey
colours = snakemake@params$colors
plot_props = snakemake@params$plots_props
results_folder = snakemake@params$results_folder
xticks_names = snakemake@params$xticks_names

log_fc_thres = log2(1.5)
sig_lvl = 0.05
top_dereg_features = 5

# save rdata
#saveRDS(snakemake, file = "snakemake_mirna_expression_median_per_group_line_plot.rds")

#expr = fread(sprintf("%s_%s/%s_%s_quantification_%s.detection_rate_%s_per_%s.csv", data_input$raw_data_folder, data_input$data_sub_set, data_input$rna_class, data_input$feature_filtering, data_input$norm, data_input$detection_rate, data_input$detection_group), sep='\t', header=T)
annot = fread(sprintf("%s/annotation_%s_%s.csv", data_input$annotation_data_folder, data_input$data_sub_set, data_input$feature_filtering), sep='\t', colClasses=c(ID="character")) 
diff_exp = fread(sprintf("%s/%s_%s/results_%s/matrices/diff_exp/diff_exp%s.csv", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, diff_exp_log), sep='\t') 

output_folder_fig = sprintf("%s/%s_%s/results_%s/figures/expr/line_plots", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set)
dir.create(output_folder_fig, recursive=TRUE)
output_folder_tab = sprintf("%s/%s_%s/results_%s/matrices/diff_exp/sig_deregulated", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set)
dir.create(output_folder_tab, recursive=TRUE)

  
#------------------------------------ Script ----------------------------------- 
for (prop in properties) {
  # given features
  feature_set = c(features_black, features_grey)
  colour_set = c(rep("black", length(features_black)), rep("grey", length(features_grey)))
  names(colour_set) = feature_set
  
  diff_exp_filtered = diff_exp[diff_exp[[data_input$rna_class]] %in% feature_set,]
  
  median_table_df = data.frame()
  for (sub_group in unique(annot[[prop]])){
    median_table = data.table(feature_prop=sprintf("%s_%s", diff_exp_filtered[[data_input$rna_class]], sub_group))
    median_table$prop = sub_group
    median_table$feature = diff_exp_filtered[[data_input$rna_class]]
    median_table$value= 10^(as.numeric(diff_exp_filtered[[sprintf("median_%s__%s", prop, sub_group)]]))
    median_table$value_log= as.numeric(diff_exp_filtered[[sprintf("median_%s__%s", prop, sub_group)]])
    
    median_table_df = rbind(median_table_df, data.frame(median_table, check.names = FALSE))
  }
  
  file_name = sprintf("prop_comparison_for_%s_%ss", length(feature_set), data_input$rna_class)
  line_plot(median_table_df, colour_set, prop, data_input, xticks_names, plot_props, output_folder_fig, file_name)
  
  
  # sig dregulated features
  info_table = data.table(feature=diff_exp[[data_input$rna_class]],
                          ttest_rawp_pval=diff_exp[[sprintf("ttest_rawp_%s__%s_vs_%s", prop, "old", "young")]],  
                          fc=diff_exp[[sprintf("fc_%s__%s_vs_%s", prop, "old", "young")]]) 
  info_table[,logfc := log2(fc)]
  info_table[,logp := log10(ttest_rawp_pval)]
  
  info_table$fc_cat_up = rep(FALSE, dim(info_table)[1])
  info_table$fc_cat_down = rep(FALSE, dim(info_table)[1])
  info_table$ttest_rawp_cat_up = rep(FALSE, dim(info_table)[1])
  info_table$ttest_rawp_cat_down = rep(FALSE, dim(info_table)[1])
  
  info_table[logfc >= log_fc_thres, fc_cat_up:=TRUE]
  info_table[logfc <= -log_fc_thres, fc_cat_down:=TRUE]
  info_table[logfc >= log_fc_thres & (ttest_rawp_pval < sig_lvl), ttest_rawp_cat_up:=TRUE]
  info_table[logfc <= -log_fc_thres & (ttest_rawp_pval < sig_lvl), ttest_rawp_cat_down:=TRUE]

  sig_dereg_table = data.frame(feature = c(info_table[info_table$ttest_rawp_cat_up == TRUE]$feature, info_table[info_table$ttest_rawp_cat_down == TRUE]$feature), 
                            direction = c(rep("up", length(info_table[info_table$ttest_rawp_cat_up == TRUE]$feature)), rep("down", length(info_table[info_table$ttest_rawp_cat_down == TRUE]$feature))))
  
  fwrite(sig_dereg_table, sprintf("%s/sig_deregulated_%ss.csv", output_folder_tab, data_input$rna_class), sep = "\t", row.names = FALSE)
  write.xlsx(sig_dereg_table, sprintf("%s/sig_deregulated_%ss.xlsx", output_folder_tab, data_input$rna_class), colNames = TRUE, rowNames = FALSE, append = FALSE)
  
  info_table_sig_de_reg = info_table[((info_table$ttest_rawp_cat_up == TRUE) | (info_table$ttest_rawp_cat_down == TRUE)),]

  info_table_sig_de_reg$sort_col = (sign(info_table_sig_de_reg$logfc) * log10(info_table_sig_de_reg$ttest_rawp_pval))
  info_table_sig_de_reg = info_table_sig_de_reg[order(info_table_sig_de_reg$sort_col, decreasing = TRUE),]  
  
  info_table_cat = rbind(tail(info_table_sig_de_reg, top_dereg_features), head(info_table_sig_de_reg, top_dereg_features))
  info_table_cat = info_table_cat[order(info_table_cat$sort_col, decreasing = TRUE),]  
  
  
  colour_set = c(rep("#488468", dim(info_table_cat[info_table_cat$sort_col > 0,])[1]), rep("#E09900", dim(info_table_cat[info_table_cat$sort_col < 0,])[1]))
  names(colour_set) = info_table_cat$feature
  

  diff_exp_filtered = diff_exp[diff_exp[[data_input$rna_class]] %in% info_table_cat$feature,]
  
  median_table_df = data.frame()
  for (sub_group in unique(annot[[prop]])){
    median_table = data.table(feature_prop=sprintf("%s_%s", diff_exp_filtered[[data_input$rna_class]], sub_group))
    median_table$prop = sub_group
    median_table$feature = diff_exp_filtered[[data_input$rna_class]]
    median_table$value= 10^(as.numeric(diff_exp_filtered[[sprintf("median_%s__%s", prop, sub_group)]]))
    median_table$value_log= as.numeric(diff_exp_filtered[[sprintf("median_%s__%s", prop, sub_group)]])
    
    median_table_df = rbind(median_table_df, data.frame(median_table, check.names = FALSE))
  }
  
  file_name = sprintf("prop_comparison_for_top%s_sig_dereg_%ss", top_dereg_features, data_input$rna_class)
  line_plot(median_table_df, colour_set, prop, data_input, xticks_names, plot_props, output_folder_fig, file_name)

  
  
}


#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])
