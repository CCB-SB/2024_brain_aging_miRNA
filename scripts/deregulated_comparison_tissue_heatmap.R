suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(svglite))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(openxlsx))

#snakemake = readRDS("snakemake_deregulated_comparison_tissue_heatmap.rds")

set.seed(snakemake@params$parameters_props$set_seed)

print("deregulated_comparison_tissue_heatmap")


#---------------------------------- Functions ----------------------------------
#plot_heatmap_number_dereg = function(df, row_names, xticks_names, plots_props, prop, annotation_row, colors, color_bar_name, color_bar_colors, color_bar_max) {
#  
#  # convert to data.frame beacause of rownames
#  df = as.data.frame(df)
#  rownames(df) = row_names
#  
#  col_fun = colorRamp2(c(-color_bar_max, - (color_bar_max/10), 0,  color_bar_max/10, color_bar_max), color_bar_colors)
#  
#  #cell_fun = function(j, i, x, y, width, height, fill) {
#  #  grid.text(sprintf("%.0f", as.integer(df[i, j])), x, y, gp = gpar(fontsize = 4, fontfamily=plots_props$font_family))
#  #}
#  
#  annotation_colors = list()
#  ucolors <- unlist(colors[[prop]])
#  if(is.vector(ucolors) && all(annotation_row %in% names(ucolors))) {
#    annotation_colors[[prop]] = ucolors
#  }
#  
#  # sort the annotation according to the dataframe
#  #annotation_row = annotation_row[rownames(df), , drop=FALSE]
#  
#  # change rownames for df and annotation
#  rownames(df) = unlist(xticks_names[[tissue]][rownames(df)])
#  rownames(annotation_row) = unlist(xticks_names[[tissue]][rownames(annotation_row)])
#  
#  column_ha = rowAnnotation(df = annotation_row, col = annotation_colors, gp = gpar(row = "white"), simple_anno_size = unit(0.2, "cm"),
#                            annotation_name_side = "top", annotation_name_gp = gpar(fontsize = plots_props$font_size, fontface = "bold", fontfamily=plots_props$font_family),
#                            show_legend=FALSE, show_annotation_name=FALSE)
#  p = ComplexHeatmap::Heatmap(df, col = col_fun, #cell_fun = cell_fun,
#                              #rect_gp = gpar(col = "white", lwd = .5),
#                              # height = unit(0.3, "cm") * nrow(t(df)), width = unit(0.3, "cm") * ncol(t(df)),
#                              #height = unit(plots_props$image_height  - 2, "cm"), width = unit(plots_props$image_width - 4, "cm"),
#                              height = unit(5, "cm"), width = unit(2.5, "cm"),
#                              show_column_names = TRUE,
#                              show_row_names = TRUE,
#                              row_names_side = "left", row_names_gp = gpar(fontsize = plots_props$font_size, fontfamily=plots_props$font_family), 
#                              column_names_side = "bottom", column_names_gp = gpar(fontsize = plots_props$font_size, fontfamily=plots_props$font_family),
#                              column_names_rot = 0,
#                              column_names_centered = TRUE,
#                              # cluster_columns = function(m) hclust(dist(m), method = "complete"), show_column_dend = TRUE, clustering_distance_columns = "euclidean", 
#                              cluster_columns = FALSE,
#                              #cluster_rows = function(m) hclust(dist(m), method = "complete"), show_row_dend = TRUE, clustering_distance_rows = "euclidean", 
#                              column_dend_reorder = FALSE, #cluster_columns = TRUE,
#                              row_dend_reorder = FALSE, cluster_rows = FALSE,
#                              show_heatmap_legend = FALSE,
#                              left_annotation = column_ha
#  )
#  legend_ticks = seq(from = -color_bar_max, to = color_bar_max, length.out = 5)
#  legend_tick_labels = round(legend_ticks, digits = 2)
#  
#  legend = ComplexHeatmap::Legend(title = color_bar_name, col_fun = col_fun,
#                                  at = legend_tick_labels, #direction="horizontal",
#                                  title_gp = gpar(fontsize = plots_props$font_size, fontface = "bold", fontfamily=plots_props$font_family), 
#                                  labels_gp = gpar(fontsize = plots_props$font_size, fontfamily=plots_props$font_family),
#                                  legend_width=unit(3, "cm"), legend_height=unit(3, "cm"))
#  
#  return(list("heatmap" = p, "legend" = legend))
#}

plot_heatmap_number_dereg_adj = function(df, sig_tag, th, row_names, xticks_names, plots_props, prop, props_col, annotation_row, colors, color_bar_name, color_bar_colors, color_bar_max) {
 
  # convert to data.frame beacause of rownames
  df = as.data.frame(df)
  rownames(df) = row_names
  
  #print(sig_tag)
  if (sig_tag == TRUE) {
    scaling_factor = 40
    color_bar_max = color_bar_max
  } else {
    scaling_factor = 10
    color_bar_max = th
  }
  
  col_fun = colorRamp2(c(-color_bar_max, - (color_bar_max/scaling_factor), 0,  color_bar_max/scaling_factor, color_bar_max), color_bar_colors)
  
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (abs(df[i, j]) >= th){
      grid.text(sprintf("%.0f", as.integer(abs(df[i, j]))), x, y, gp = gpar(fontsize = 6, fontfamily=plots_props$font_family, col = "white"))
    } else {
      grid.text("", x, y, gp = gpar(fontsize = 6, fontfamily=plots_props$font_family))
    }
  }
  
  mask = df > th
  df_th = df
  # commenting in this line destroys the colors of the upregulated
  #df_th[mask] = th

  annotation_colors = list()
  ucolors = unlist(colors[[prop]])
  if(is.vector(ucolors) && all(annotation_row %in% names(ucolors))) {
    annotation_colors[[prop]] = ucolors
  }
  
  # sort the annotation according to the dataframe
  #annotation_row = annotation_row[rownames(df_th), , drop=FALSE]
  
  # change rownames for df and annotation
  rownames(df_th) = unlist(xticks_names[[tissue]][rownames(df_th)])
  rownames(annotation_row) = unlist(xticks_names[[tissue]][rownames(annotation_row)])
  
  annotation_row_df = data.frame(col=annotation_row)
  colnames(annotation_row_df) = c(prop)
  
  column_ha = rowAnnotation(df = annotation_row_df, col = annotation_colors, gp = gpar(row = "white"), simple_anno_size = unit(0.2, "cm"),
                            annotation_name_side = "top", annotation_name_gp = gpar(fontsize = plots_props$font_size, fontface = "bold", fontfamily=plots_props$font_family),
                            show_legend=FALSE, show_annotation_name=FALSE)
  
  p = ComplexHeatmap::Heatmap(df_th, col = col_fun, cell_fun = cell_fun,
                              #rect_gp = gpar(col = "white", lwd = .5),
                              # height = unit(0.3, "cm") * nrow(t(df)), width = unit(0.3, "cm") * ncol(t(df)),
                              #height = unit(plots_props$image_height  - 2, "cm"), width = unit(plots_props$image_width - 4, "cm"),
                              height = unit(4.5, "cm"), width = unit(2.5, "cm"),
                              show_column_names = TRUE,
                              column_title = xticks_names$categories[[sprintf("%s_capital", props_col)]], column_title_side = "bottom", column_title_gp = gpar(fontsize = plots_props$font_size, fontfamily=plots_props$font_family),
                              show_row_names = TRUE,
                              row_title = "",
                              row_names_side = "left", row_names_gp = gpar(fontsize = plots_props$font_size, fontfamily=plots_props$font_family), 
                              column_names_side = "bottom", column_names_gp = gpar(fontsize = plots_props$font_size, fontfamily=plots_props$font_family),
                              column_names_rot = 0,
                              column_names_centered = TRUE,
                              # cluster_columns = function(m) hclust(dist(m), method = "complete"), show_column_dend = TRUE, clustering_distance_columns = "euclidean", 
                              cluster_columns = FALSE,
                              #cluster_rows = function(m) hclust(dist(m), method = "complete"), show_row_dend = TRUE, clustering_distance_rows = "euclidean", 
                              column_dend_reorder = FALSE, #cluster_columns = TRUE,
                              row_dend_reorder = FALSE, cluster_rows = FALSE,
                              show_heatmap_legend = FALSE,
                              #left_annotation = column_ha
  )

  if (sig_tag == TRUE) {
    legend_ticks = seq(from = -th, to = th, length.out = 5)
  } else {
    legend_ticks = seq(from = -color_bar_max, to = color_bar_max, length.out = 5)
  }
  
  legend_tick_labels = round(legend_ticks, digits = 2)

  legend = ComplexHeatmap::Legend(title = color_bar_name, col_fun = col_fun,
                                  at = legend_tick_labels, #direction="horizontal",
                                  title_gp = gpar(fontsize = plots_props$font_size, fontface = "bold", fontfamily=plots_props$font_family), 
                                  labels_gp = gpar(fontsize = plots_props$font_size, fontfamily=plots_props$font_family),
                                  legend_width=unit(3, "cm"), legend_height=unit(3, "cm"))
  return(list("heatmap" = p, "legend" = legend))
}

row_name_squares = function(tissue_list, xticks_names, tissue, colors, plots_props) {
  # Create a sample matrix
  data =- matrix(1:length(tissue_list), nrow = length(tissue_list), ncol = 1)
  rownames(data) = tissue_list
  colnames(data) = "Square"
  
  # Create a function to draw text inside squares
  draw_text = function(j, i, x, y, width, height, fill) { 
    wording = xticks_names[[tissue]][[rownames(data)[i]]]
    text_colour = unlist(colors[[sprintf("%s_text", tissue)]])[[rownames(data)[i]]]
    rect_colour = unlist(colors[[tissue]])[[rownames(data)[i]]]
    grid.roundrect(x=x, y=y, width=width, height=height, r = unit(0.1, "cm"), gp = gpar(fill = rect_colour, col = "white", lwd = 2))
    grid.text(wording, x = x, y = y, gp = gpar(col = text_colour, fontsize = 6)) #fontfamily=plots_props$font_family
  }
  
  #col_fun  = structure(unname(unlist(colors[[tissue]])), names = rownames(data))
  
  # Create the heatmap
  p = Heatmap(data,
              #col = col_fun,
              height = unit(4.5, "cm"), width = unit(1, "cm"),
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              show_row_names = FALSE,
              show_column_names = FALSE,
              show_heatmap_legend = FALSE,
              cell_fun = draw_text,
              rect_gp = gpar(type = "none")
  )
  return(p)
}


#---------------------------------- Load data ---------------------------------- 
data_input = snakemake@params$data_input
diff_log = snakemake@params$diff_exp_log
log_fc_thres = log2(snakemake@params$thresholds$fc)
sig_lvl = snakemake@params$thresholds$adj_p_value
props = snakemake@params$props
xticks_names = snakemake@params$xticks_names
plots_props = snakemake@params$plots_props
colors = snakemake@params$colors
results_folder = snakemake@params$results_folder
number_of_deregulated_tissues_thresholds = snakemake@params$thresholds$num_dereg_tissues_th
number_of_deregulated_comparison_thresholds = snakemake@params$thresholds$num_dereg_combinations_th
ID = data_input$identifier_column

# save rdata
#saveRDS(snakemake, file = "snakemake_deregulated_comparison_tissue_heatmap.rds")

annot = fread(sprintf("%s/annotation_%s_%s.csv", data_input$annotation_data_folder, data_input$data_sub_set, data_input$feature_filtering), sep='\t', colClasses=c(ID="character")) 
diff_exp = fread(sprintf("%s/%s_%s/results_%s/matrices/diff_exp/diff_exp%s.csv", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, diff_log), sep='\t') 


#------------------------------------ Script ----------------------------------- 
for (prop in props) {
  
  tissue = prop[1]
  time = prop[2]
  
  output_folder_fig = sprintf("%s/%s_%s/results_%s/figures/diff_exp/heatmap_complex/num_fc_th_%s_per_%s", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, snakemake@params$thresholds$fc, tissue)
  dir.create(output_folder_fig, recursive=TRUE)
  output_folder_tab = sprintf("%s/%s_%s/results_%s/matrices/aggregation_diff_exp", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set)
  dir.create(output_folder_tab, recursive=TRUE)
  
  up_reg = c()
  down_reg = c()
  up_reg_contain = c()
  down_reg_contain = c()
  sig_up_reg = c()
  sig_down_reg = c()
  sig_up_reg_contain = c()
  sig_down_reg_contain = c()
  for (ti in unique(annot[[tissue]])) {
    up_reg_tissue = c()
    down_reg_tissue = c()
    sig_up_reg_tissue = c()
    sig_down_reg_tissue = c()
    for (time_case in sort(as.numeric(unique(annot[[time]])))[-1]) {
      info_table = data.table(feature=diff_exp[[data_input$rna_class]],
                              ttest_adj_pval=diff_exp[[sprintf("ttest_adjp_%s__%s_vs_%s_%s=%s", time, time_case, sort(as.numeric(unique(annot[[time]])))[1], tissue, ti)]],  #ttest_adjp_age__12_vs_3_brain_region=cc
                              fc=diff_exp[[sprintf("fc_%s__%s_vs_%s_%s=%s", time, time_case, sort(as.numeric(unique(annot[[time]])))[1], tissue, ti)]]) #fc_age__12_vs_3_brain_region=cc
      info_table[,logfc := log2(fc)]
      
      info_table$fc_cat_up = rep(FALSE, dim(info_table)[1])
      info_table$fc_cat_down = rep(FALSE, dim(info_table)[1])
      info_table$ttest_adj_cat_up = rep(FALSE, dim(info_table)[1])
      info_table$ttest_adj_cat_down = rep(FALSE, dim(info_table)[1])
      
      #info_table[ttest_adj_pval >= sig_lvl, ttest_adj_cat:="not_sig"]
      #info_table[ttest_adj_pval < sig_lvl, ttest_adj_cat:="sig"]
      info_table[logfc >= log_fc_thres, fc_cat_up:=TRUE]
      info_table[logfc <= -log_fc_thres, fc_cat_down:=TRUE]
      info_table[logfc >= log_fc_thres & (ttest_adj_pval < sig_lvl), ttest_adj_cat_up:=TRUE]
      info_table[logfc <= -log_fc_thres & (ttest_adj_pval < sig_lvl), ttest_adj_cat_down:=TRUE]
      
      up_reg_tissue[[sprintf("%s", time_case)]] = data.frame("feature" = info_table[info_table$fc_cat_up,]$feature)
      down_reg_tissue[[sprintf("%s", time_case)]] = data.frame("feature" = info_table[info_table$fc_cat_down,]$feature)
      
      sig_up_reg_tissue[[sprintf("%s", time_case)]] = data.frame("feature" = info_table[info_table$ttest_adj_cat_up,]$feature)
      sig_down_reg_tissue[[sprintf("%s", time_case)]] = data.frame("feature" = info_table[info_table$ttest_adj_cat_down,]$feature)
    }
    
    up_reg_tissue_contain = data.table(feature = unique(unlist(up_reg_tissue, use.names = FALSE)))
    down_reg_tissue_contain = data.table(feature = unique(unlist(down_reg_tissue, use.names = FALSE)))
    
    sig_up_reg_tissue_contain = data.table(feature = unique(unlist(sig_up_reg_tissue, use.names = FALSE)))
    sig_down_reg_tissue_contain = data.table(feature = unique(unlist(sig_down_reg_tissue, use.names = FALSE)))
    
    # add tissue column
    up_reg_tissue_contain$tissue = ti
    down_reg_tissue_contain$tissue = ti
    
    sig_up_reg_tissue_contain$tissue = ti
    sig_down_reg_tissue_contain$tissue = ti
    
    for (time_case in sort(as.numeric(unique(annot[[time]])))[-1]) {
      up_reg_tissue_contain[[sprintf("%s", time_case)]] = rep(FALSE, length(up_reg_tissue_contain$feature))
      down_reg_tissue_contain[[sprintf("%s", time_case)]] = rep(FALSE, length(down_reg_tissue_contain$feature))
      
      up_reg_tissue_contain[feature %in% up_reg_tissue[[sprintf("%s", time_case)]]$feature, sprintf("%s",time_case):=TRUE]
      down_reg_tissue_contain[feature %in% down_reg_tissue[[sprintf("%s", time_case)]]$feature, sprintf("%s",time_case):=TRUE]
      
      
      sig_up_reg_tissue_contain[[sprintf("%s", time_case)]] = rep(FALSE, length(sig_up_reg_tissue_contain$feature))
      sig_down_reg_tissue_contain[[sprintf("%s", time_case)]] = rep(FALSE, length(sig_down_reg_tissue_contain$feature))
      
      sig_up_reg_tissue_contain[feature %in% sig_up_reg_tissue[[sprintf("%s", time_case)]]$feature, sprintf("%s",time_case):=TRUE]
      sig_down_reg_tissue_contain[feature %in% sig_down_reg_tissue[[sprintf("%s", time_case)]]$feature, sprintf("%s",time_case):=TRUE]
    }
    
    up_reg_contain = rbind(up_reg_contain, up_reg_tissue_contain)
    down_reg_contain = rbind(down_reg_contain, down_reg_tissue_contain)
    
    sig_up_reg_contain = rbind(sig_up_reg_contain, sig_up_reg_tissue_contain)
    sig_down_reg_contain = rbind(sig_down_reg_contain, sig_down_reg_tissue_contain)
  }
  
  #######################
  # up-down
  # make dataframe per tissue, per time_case -> number of TRUE
  up_reg_heatmap = c()
  down_reg_heatmap = c()
  for (ti in unique(annot[[tissue]])) {
    up_reg_count_tissue_time = data.table("tissue" = ti)
    down_reg_count_tissue_time = data.table("tissue" = ti)
    
    # make colSums for every tissue, do not use first 2 columns (feature and tissue)
    up_reg_count_per_tissue = cbind(up_reg_count_tissue_time, t(colSums(up_reg_contain[up_reg_contain$tissue == ti, -1:-2])))
    down_reg_count_per_tissue = cbind(down_reg_count_tissue_time, t(colSums(down_reg_contain[down_reg_contain$tissue == ti, -1:-2])))
    
    up_reg_heatmap = rbind(up_reg_heatmap, up_reg_count_per_tissue)
    down_reg_heatmap = rbind(down_reg_heatmap, down_reg_count_per_tissue)
  }
  
  # plot and save 2 heatmaps
  # calculate max for the colorbar
  color_bar_max = round(max(abs(up_reg_heatmap[, -1]), abs(down_reg_heatmap[, -1])),-2)
  if (color_bar_max == 0) {
    color_bar_max = max(abs(up_reg_heatmap[, -1]), abs(down_reg_heatmap[, -1]))
  }
  if (color_bar_max == 0) {
    color_bar_max = 1
  }
  annot_row = data.frame(Tissue=unique(annot[[tissue]]))
  # rename to brain_region in this case
  colnames(annot_row) = c(tissue)
  rownames(annot_row) = unique(annot[[tissue]])
  row_names = up_reg_heatmap$tissue
  #colorbar_name = sprintf("num. of\ndereg. %ss\n(DEGs relative\nto 3m)", data_input$rna_class)
  colorbar_name = sprintf("Num. of\ndereg.\n%ss", data_input$rna_class)
  # colors for 0 (darkgrey) and a gradient from red4 to gold for values > 0
  #colorbar_colors = c("#2B5D93", "#7DC8DC", "#F4F5F5", "#F5AD29", "#ED2E26")
  #colorbar_colors = c("#488468", "#91C2AB", "#F4F5F5", "#E9BD5B", "#E09900") # down to up
  colorbar_colors = c(colors$direction$down, colors$direction$down_light, colors$direction$light, colors$direction$up_light, colors$direction$up)
    
  #return_up = plot_heatmap_number_dereg(up_reg_heatmap[, -1], row_names, xticks_names, plots_props, tissue, annot_row, colors, colorbar_name, colorbar_colors, color_bar_max)
  #
  #row_names = down_reg_heatmap$tissue
  #return_down = plot_heatmap_number_dereg(-down_reg_heatmap[, -1], row_names, xticks_names, plots_props, tissue, annot_row, colors, colorbar_name, colorbar_colors, color_bar_max)
  #
  #p_list = return_down$heatmap + return_up$heatmap
  #
  # save
  # figures
  #file_name = sprintf("%s/diff_exp_global_de_reg_%s_per_%s", output_folder_fig, time, data_input$rna_class)
  #png(sprintf("%s.png", file_name), width = plots_props$image_width, height = plots_props$image_height, units = plots_props$image_units, res = plots_props$dpi)
  #draw(p_list, padding = unit(c(0.25, -1.2, 0.25, 0.25), "cm"))  # bottom, left, top and right margins
  #draw(return_up$legend, x = unit((plots_props$image_width - 0.25), "cm"), y = unit((plots_props$image_height / 2), "cm"), just = c("right", "center"))
  #dev.off()
  #svglite(sprintf("%s.svg", file_name), width = (plots_props$image_width / 2.54), height = plots_props$image_height / 2.54)
  #draw(p_list, padding = unit(c(0.25, -1.2, 0.25, 0.25), "cm"))  # bottom, left, top and right margins
  #draw(return_up$legend, x = unit((plots_props$image_width - 0.25), "cm"), y = unit((plots_props$image_height / 2), "cm"), just = c("right", "center"))
  #dev.off()
  
  # order the brain regions in an alphabetical order given by the orer of the colour list
  target_order = names(colors$brain_region)[names(colors$brain_region) %in% up_reg_heatmap$tissue]
  order_indices = match(target_order, up_reg_heatmap$tissue)
  up_reg_heatmap = up_reg_heatmap[order_indices,]
  order_indices = match(target_order, down_reg_heatmap$tissue)
  down_reg_heatmap = down_reg_heatmap[order_indices,]
  
  sig_tag = FALSE
  th = 200
  return_up_adj = plot_heatmap_number_dereg_adj(up_reg_heatmap[, -1], sig_tag, th, row_names, xticks_names, plots_props, tissue, time, target_order, colors, colorbar_name, colorbar_colors, th)
                                              
  #row_names = down_reg_heatmap$tissue
  row_names = target_order
  return_down_adj = plot_heatmap_number_dereg_adj(-down_reg_heatmap[, -1], sig_tag, th, row_names, xticks_names, plots_props, tissue, time, target_order, colors, colorbar_name, colorbar_colors, th)
  
  row_name_square_plot = row_name_squares(row_names, xticks_names, tissue, colors, plots_props)
  
  p_list_adj = row_name_square_plot + return_up_adj$heatmap + return_down_adj$heatmap
  
  # save
  # figures
  file_name = sprintf("%s/diff_exp_global_de_reg_%s_per_%s_adj", output_folder_fig, time, data_input$rna_class)
  png(sprintf("%s.png", file_name), width = plots_props$image_width, height = plots_props$image_height, units = plots_props$image_units, res = plots_props$dpi)
  draw(p_list_adj, padding = unit(c(0.25, -1.2, 0.25, 0.25), "cm"))  # bottom, left, top and right margins
  draw(return_up_adj$legend, x = unit((plots_props$image_width - 0.25), "cm"), y = unit((plots_props$image_height / 2), "cm"), just = c("right", "center"))
  dev.off()
  svglite(sprintf("%s.svg", file_name), width = (plots_props$image_width / 2.54), height = plots_props$image_height / 2.54)
  draw(p_list_adj, padding = unit(c(0.25, -1.2, 0.25, 0.25), "cm"))  # bottom, left, top and right margins
  draw(return_up_adj$legend, x = unit((plots_props$image_width - 0.25), "cm"), y = unit((plots_props$image_height / 2), "cm"), just = c("right", "center"))
  dev.off()
  
  file_name = sprintf("%s/diff_exp_global_de_reg_%s_per_%s_adj_without_legend", output_folder_fig, time, data_input$rna_class)
  png(sprintf("%s.png", file_name), width = plots_props$image_width, height = plots_props$image_height, units = plots_props$image_units, res = plots_props$dpi)
  draw(p_list_adj, padding = unit(c(0.25, 0, 0.25, 0.25), "cm"))  # bottom, left, top and right margins
  #draw(return_up_adj$legend, x = unit((4 * plots_props$image_width / 3), "cm"), y = unit((plots_props$image_height / 2), "cm"), just = c("right", "center"))
  dev.off()
  svglite(sprintf("%s.svg", file_name), width = (plots_props$image_width/ 2.54), height = plots_props$image_height / 2.54)
  draw(p_list_adj, padding = unit(c(0.25, 0, 0.25, 0.25), "cm"))  # bottom, left, top and right margins
  #draw(return_up_adj$legend, x = unit((4 * plots_props$image_width / 3), "cm"), y = unit((plots_props$image_height / 2), "cm"), just = c("right", "center"))
  dev.off()
  
  # table
  fwrite(up_reg_heatmap, sprintf("%s/diff_exp_global_up_reg_%s_per_%s.csv", output_folder_tab, time, data_input$rna_class), sep = "\t", row.names = TRUE)
  write.xlsx(up_reg_heatmap, sprintf("%s/diff_exp_global_up_reg_%s_per_%s.xlsx", output_folder_tab, time, data_input$rna_class), colNames = TRUE, rowNames = TRUE, append = FALSE)
  
  fwrite(down_reg_heatmap, sprintf("%s/diff_exp_global_down_reg_%s_per_%s.csv", output_folder_tab, time, data_input$rna_class), sep = "\t", row.names = TRUE)
  write.xlsx(down_reg_heatmap, sprintf("%s/diff_exp_global_down_reg_%s_per_%s.xlsx", output_folder_tab, time, data_input$rna_class), colNames = TRUE, rowNames = TRUE, append = FALSE)
  #######################
  
  #######################
  # sig up-down
  # make dataframe per tissue, per time_case -> number of TRUE
  sig_up_reg_heatmap = c()
  sig_down_reg_heatmap = c()
  for (ti in unique(annot[[tissue]])) {
    sig_up_reg_count_tissue_time = data.table("tissue" = ti)
    sig_down_reg_count_tissue_time = data.table("tissue" = ti)
    
    # make colSums for every tissue, do not use first 2 columns (feature and tissue)
    sig_up_reg_count_per_tissue = cbind(sig_up_reg_count_tissue_time, t(colSums(sig_up_reg_contain[sig_up_reg_contain$tissue == ti, -1:-2])))
    sig_down_reg_count_per_tissue = cbind(sig_down_reg_count_tissue_time, t(colSums(sig_down_reg_contain[sig_down_reg_contain$tissue == ti, -1:-2])))
    
    sig_up_reg_heatmap = rbind(sig_up_reg_heatmap, sig_up_reg_count_per_tissue)
    sig_down_reg_heatmap = rbind(sig_down_reg_heatmap, sig_down_reg_count_per_tissue)
  }
  
  # plot and save 2 heatmaps
  # calculate max for the colorbar
  color_bar_max = round(max(abs(sig_up_reg_heatmap[, -1]), abs(sig_down_reg_heatmap[, -1])),-2)
  if (color_bar_max == 0) {
    color_bar_max = max(abs(sig_up_reg_heatmap[, -1]), abs(sig_down_reg_heatmap[, -1]))
  }
  if (color_bar_max == 0) {
    color_bar_max = 1
  }
  annot_row = data.frame(Tissue=unique(annot[[tissue]]))
  # rename to brain_region in this case
  colnames(annot_row) = c(tissue)
  rownames(annot_row) = unique(annot[[tissue]])
  row_names = sig_up_reg_heatmap$tissue
  #colorbar_name = sprintf("num. of sig.\ndereg. %ss\n(DEGs relative\nto 3m)", data_input$rna_class)
  colorbar_name = sprintf("num. of\nsig.\ndereg.\n%ss", data_input$rna_class)
  # colors for 0 (darkgrey) and a gradient from red4 to gold for values > 0
  
  #return_up = plot_heatmap_number_dereg(sig_up_reg_heatmap[, -1], row_names, xticks_names, plots_props, tissue, annot_row, colors, colorbar_name, colorbar_colors, color_bar_max)
  #  
  #row_names = sig_down_reg_heatmap$tissue
  #return_down = plot_heatmap_number_dereg(-sig_down_reg_heatmap[, -1], row_names, xticks_names, plots_props, tissue, annot_row, colors, colorbar_name, colorbar_colors, color_bar_max)
  #
  #p_list = return_down$heatmap + return_up$heatmap
  #
  # save
  # figures
  #file_name = sprintf("%s/diff_exp_global_sig_de_reg_%s_per_%s", output_folder_fig, time, data_input$rna_class)
  #png(sprintf("%s.png", file_name), width = plots_props$image_width, height = plots_props$image_height, units = plots_props$image_units, res = plots_props$dpi)
  #draw(p_list, padding = unit(c(0.25, -1.2, 0.25, 0.25), "cm"))  # bottom, left, top and right margins
  #draw(return_up$legend, x = unit((plots_props$image_width - 0.25), "cm"), y = unit((plots_props$image_height / 2), "cm"), just = c("right", "center"))
  #dev.off()
  #svglite(sprintf("%s.svg", file_name), width = (plots_props$image_width / 2.54), height = plots_props$image_height / 2.54)
  #draw(p_list, padding = unit(c(0.25, -1.2, 0.25, 0.25), "cm"))  # bottom, left, top and right margins
  #draw(return_up$legend, x = unit((plots_props$image_width - 0.25), "cm"), y = unit((plots_props$image_height / 2), "cm"), just = c("right", "center"))
  #dev.off()
  
  # order the brain regions in an alphabetical order given by the orer of the colour list
  target_order = names(colors$brain_region)[names(colors$brain_region) %in% sig_up_reg_heatmap$tissue]
  order_indices = match(target_order, sig_up_reg_heatmap$tissue)
  sig_up_reg_heatmap = sig_up_reg_heatmap[order_indices,]
  order_indices = match(target_order, sig_down_reg_heatmap$tissue)
  sig_down_reg_heatmap = sig_down_reg_heatmap[order_indices,]
  
  sig_tag = TRUE
  th = 10
  #print(color_bar_max)
  return_up_adj = plot_heatmap_number_dereg_adj(sig_up_reg_heatmap[, -1], sig_tag, th, row_names, xticks_names, plots_props, tissue, time, target_order, colors, colorbar_name, colorbar_colors, color_bar_max)
  #row_names = sig_down_reg_heatmap$tissue
  row_names = target_order
  return_down_adj = plot_heatmap_number_dereg_adj(-sig_down_reg_heatmap[, -1], sig_tag, th, row_names, xticks_names, plots_props, tissue, time, target_order, colors, colorbar_name, colorbar_colors, color_bar_max)
  
  p_list_adj = row_name_square_plot + return_up_adj$heatmap + return_down_adj$heatmap
  
  # save
  # figures
  file_name = sprintf("%s/diff_exp_global_sig_de_reg_%s_per_%s_adj", output_folder_fig, time, data_input$rna_class)
  png(sprintf("%s.png", file_name), width = plots_props$image_width, height = plots_props$image_height, units = plots_props$image_units, res = plots_props$dpi)
  draw(p_list_adj, padding = unit(c(0.25, -1.2, 0.25, 0.25), "cm"))  # bottom, left, top and right margins
  draw(return_up_adj$legend, x = unit((plots_props$image_width - 0.25), "cm"), y = unit((plots_props$image_height / 2), "cm"), just = c("right", "center"))
  dev.off()
  svglite(sprintf("%s.svg", file_name), width = (plots_props$image_width / 2.54), height = plots_props$image_height / 2.54)
  draw(p_list_adj, padding = unit(c(0.25, -1.2, 0.25, 0.25), "cm"))  # bottom, left, top and right margins
  draw(return_up_adj$legend, x = unit((plots_props$image_width - 0.25), "cm"), y = unit((plots_props$image_height / 2), "cm"), just = c("right", "center"))
  dev.off()
  
  # figures
  file_name = sprintf("%s/diff_exp_global_sig_de_reg_%s_per_%s_adj_without_legend", output_folder_fig, time, data_input$rna_class)
  png(sprintf("%s.png", file_name), width = plots_props$image_width, height = plots_props$image_height, units = plots_props$image_units, res = plots_props$dpi)
  draw(p_list_adj, padding = unit(c(0.25, 0, 0.25, 0.25), "cm"))  # bottom, left, top and right margins
  #draw(return_up_adj$legend, x = unit((4 * plots_props$image_width / 3), "cm"), y = unit((plots_props$image_height / 2), "cm"), just = c("right", "center"))
  dev.off()
  svglite(sprintf("%s.svg", file_name), width = (plots_props$image_width / 2.54), height = plots_props$image_height / 2.54)
  draw(p_list_adj, padding = unit(c(0.25, 0, 0.25, 0.25), "cm"))  # bottom, left, top and right margins
  #draw(return_up_adj$legend, x = unit((4 * plots_props$image_width / 3), "cm"), y = unit((plots_props$image_height / 2), "cm"), just = c("right", "center"))
  dev.off()
  
  #tables
  fwrite(sig_up_reg_heatmap, sprintf("%s/diff_exp_global_sig_up_reg_%s_per_%s.csv", output_folder_tab, time, data_input$rna_class), sep = "\t", row.names = TRUE)
  write.xlsx(sig_up_reg_heatmap, sprintf("%s/diff_exp_global_sig_up_reg_%s_per_%s.xlsx", output_folder_tab, time, data_input$rna_class), colNames = TRUE, rowNames = TRUE, append = FALSE)
  
  fwrite(sig_down_reg_heatmap, sprintf("%s/diff_exp_global_sig_down_reg_%s_per_%s.csv", output_folder_tab, time, data_input$rna_class), sep = "\t", row.names = TRUE)
  write.xlsx(sig_down_reg_heatmap, sprintf("%s/diff_exp_global_sig_down_reg_%s_per_%s.xlsx", output_folder_tab, time, data_input$rna_class), colNames = TRUE, rowNames = TRUE, append = FALSE)
  #######################
}


#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])
