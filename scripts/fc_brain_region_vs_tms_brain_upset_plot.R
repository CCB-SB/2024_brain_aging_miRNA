library(data.table)
library(UpSetR)
library(ComplexHeatmap)
#library(ggplot2)
#library(stringr)
#library(ggVennDiagram)
#library(ggvenn)
library(gridExtra)
library(svglite)

#snakemake = readRDS("snakemake_fc_brain_region_vs_tms_brain_upset_plot.rds")

print("fc_brain_region_vs_tms_brain_upset_plot")

set.seed(snakemake@params$parameters_props$set_seed)

upset_plot = function(df, mode, save_path, save_path_file_name, plots_props, th){
  plot_table = make_comb_mat(df, mode=mode)
  
  # only show N intersections
  plot_table = plot_table[1:20]

  #plot_table_thresholded = plot_table_thresholded[1:(dim(plot_table_thresholded)[1]-1),]
  cs = comb_size(plot_table)
  nc = ncol(plot_table)
  
  upset_plot = ComplexHeatmap::UpSet(plot_table, comb_order = order(cs, decreasing = TRUE),
                                     pt_size = unit(1, "mm"), lwd = 1,
                                     top_annotation = HeatmapAnnotation("Intersecting\nfeatures" = anno_barplot(cs, 
                                                                                                             ylim = c(0, max(cs)*1.1), 
                                                                                                             border = FALSE, gp = gpar(fill = "black"), 
                                                                                                             labels_gp = gpar(col = "black", family = plots_props$font_family, fontsize = plots_props$font_size), 
                                                                                                             height = unit(1, "cm")),
                                                                        annotation_name_side = "left", annotation_name_rot = 90,
                                                                        annotation_name_gp = gpar(fontsize = plots_props$font_size, fontfamily = plots_props$font_family)),
                             #top_annotation = upset_top_annotation(plot_table_thresholded[comb_degree(plot_table_thresholded) > 1], add_numbers = TRUE),
                             right_annotation = rowAnnotation("Number of\nupregulated\nfeatures" = anno_barplot(set_size(plot_table[comb_degree(plot_table) > 1]), 
                                                                                                             border = FALSE, width = unit(3, "cm")),
                                                              annotation_name_gp = gpar(fontsize = plots_props$font_size, fontfamily = plots_props$font_family)),
                             row_names_gp = gpar(fontsize = 4),  # size of brain_region labels
                             #right_annotation = NULL,
                             column_title_gp = gpar(fontsize = plots_props$font_size_header, fontfamily=plots_props$font_family),
                             column_title = sprintf("%s", save_path_file_name))
                             # gp = gpar(fill = sample_colors, col = sample_colors), 
  
  #decorate_annotation("Intersection\nsize", {grid.text(cs[co], x = 1:nc, y = unit(cs[co], "native") + unit(1, "mm"), gp = gpar(fontsize = 5), just = "bottom", default.units = "native")})
  png(sprintf("%s/%s.png", save_path, save_path_file_name), res=plots_props$dpi, width=plots_props$image_width, height=plots_props$image_height, units = plots_props$image_units)
  upset_plot = draw(upset_plot); co = column_order(upset_plot)
  dev.off()
  svglite(sprintf("%s/%s.svg", save_path, save_path_file_name), width=plots_props$image_width/2.54, height=plots_props$image_height/2.54)
  upset_plot = draw(upset_plot); co = column_order(upset_plot)
  dev.off()
}


#---------------------------------- Load data ---------------------------------- 
data_input = snakemake@params$data_input
log_fc_thres_list = as.numeric(snakemake@params$params$log_fc_thres_list)
sig_lvl = snakemake@params$params$sig_lvl
plots_props = snakemake@params$plots_props
colors = snakemake@params$colors
tissues_property = snakemake@params$params$tissues
tissue_control = snakemake@params$params$tissues_control
diff_log = snakemake@params$params$diff_log
id = data_input$identifier_column
mode = snakemake@params$params$mode
th = snakemake@params$params$plot_th
results_folder = snakemake@params$results_folder

# save rdata
#saveRDS(snakemake, file = "snakemake_fc_brain_region_vs_tms_brain_upset_plot.rds")

diff_exp = fread(sprintf("%s/%s_%s/results_%s/matrices/diff_exp/diff_exp%s.csv", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, diff_log), sep='\t') 
annot = fread(sprintf("%s/annotation_%s_%s.csv", data_input$annotation_data_folder, data_input$data_sub_set, data_input$feature_filtering), sep='\t', colClasses=c(ID="character")) 

output_folder = sprintf("%s/%s_%s/results_%s/figures/fc_comparisons/upset_plots", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set)
dir.create(output_folder, recursive=TRUE)


#------------------------------------ Script -----------------------------------
for (log_fc_thres in log_fc_thres_list) {
  
  # save folder
  output_file_name = sprintf("fc_th=%s", log_fc_thres)
  output_file_name_sig = sprintf("sig_fc_th=%s", log_fc_thres)

  up_reg_tissue = c()
  down_reg_tissue = c()
  sig_up_reg_tissue = c()
  sig_down_reg_tissue = c()
  for (tissue in unique(annot[[tissues_property]])) {
    if (tissue == tissue_control) {
      # skip if this is the control tissue
      next
    }
    info_table = data.table(feature=diff_exp[[data_input$rna_class]],
                            ttest_adj_pval=diff_exp[[sprintf("ttest_adjp_%s__%s_vs_%s", tissues_property, tissue, tissue_control)]],  #ttest_adjp_age__12_vs_3_brain_region_fixed=cc
                            fc=diff_exp[[sprintf("fc_%s__%s_vs_%s", tissues_property, tissue, tissue_control)]]) #fc_age__12_vs_3_brain_region_fixed=cc
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
    
    up_reg_tissue[[tissue]] = data.frame("feature" = info_table[info_table$fc_cat_up,]$feature)
    down_reg_tissue[[tissue]] = data.frame("feature" = info_table[info_table$fc_cat_down,]$feature)
    sig_up_reg_tissue[[tissue]] = data.frame("feature" = info_table[info_table$ttest_adj_cat_up,]$feature)
    sig_down_reg_tissue[[tissue]] = data.frame("feature" = info_table[info_table$ttest_adj_cat_down,]$feature)
  }
  up_reg_tissue_contain = data.table(feature = unique(unlist(up_reg_tissue, use.names = FALSE)))
  down_reg_tissue_contain = data.table(feature = unique(unlist(down_reg_tissue, use.names = FALSE)))
  sig_up_reg_tissue_contain = data.table(feature = unique(unlist(sig_up_reg_tissue, use.names = FALSE)))
  sig_down_reg_tissue_contain = data.table(feature = unique(unlist(sig_down_reg_tissue, use.names = FALSE)))
  
  for (tissue in unique(annot[[tissues_property]])) {
    up_reg_tissue_contain[[tissue]] = rep(FALSE, length(up_reg_tissue_contain$feature))
    down_reg_tissue_contain[[tissue]] = rep(FALSE, length(down_reg_tissue_contain$feature))
    sig_up_reg_tissue_contain[[tissue]] = rep(FALSE, length(sig_up_reg_tissue_contain$feature))
    sig_down_reg_tissue_contain[[tissue]] = rep(FALSE, length(sig_down_reg_tissue_contain$feature))
    
    up_reg_tissue_contain[up_reg_tissue_contain$feature %in% up_reg_tissue[[tissue]]$feature, sprintf("%s", tissue):=TRUE]
    down_reg_tissue_contain[down_reg_tissue_contain$feature %in% down_reg_tissue[[tissue]]$feature, sprintf("%s", tissue):=TRUE]
    sig_up_reg_tissue_contain[sig_up_reg_tissue_contain$feature %in% sig_up_reg_tissue[[tissue]]$feature, sprintf("%s", tissue):=TRUE]
    sig_down_reg_tissue_contain[sig_down_reg_tissue_contain$feature %in% sig_down_reg_tissue[[tissue]]$feature, sprintf("%s", tissue):=TRUE]
  } 
  
  # Upset plot 
  upset_plot(up_reg_tissue_contain, mode, output_folder, output_file_name, plots_props, th)
  upset_plot(sig_up_reg_tissue_contain, mode, output_folder, output_file_name_sig, plots_props, th)
}


#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])
