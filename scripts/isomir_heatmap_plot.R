suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(svglite))
#suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(matrixStats))
suppressPackageStartupMessages(library(openxlsx))


#snakemake = readRDS("snakemake_isomir_heatmap_plot.rds")

set.seed(snakemake@params$parameters_props$set_seed)

print("isomir_heatmap_plot")


#---------------------------------- Functions ----------------------------------
scale = function (x, rows, columns) {
  if(rows){
    m = apply(x, 1, mean, na.rm = T)
    s = apply(x, 1, sd, na.rm = T)
    x = (x - m)/s
  }
  if(columns){
    x = t(x)
    m = apply(x, 1, mean, na.rm = T)
    s = apply(x, 1, sd, na.rm = T)
    x = (x - m)/s
    x = t(x)
  }
  return(x)
}

create_heatmap = function(expr_feature_plot_group_medians, zscored, keep_canonical, plot_group, color_bar_name_bottom, color_bar_name_side, data_input, cutoff, feature_colours, colors, plots_props, file_name, file_name_tab) {
  # if there are too high values for one or more cells we can select a cutoff for which we set allways the same colour but display the exact values in the cells
  if (cutoff != 0) {
    col_fun = colorRamp2(c(-cutoff, - (cutoff/10), 0,  cutoff/10, cutoff), c("#DFDFB9","#525200"))
    
    cell_fun = function(j, i, x, y, width, height, fill) {
      if (abs(expr_feature_plot_group_medians[i, j]) >= cutoff){
        grid.text(sprintf("%.0f", as.integer(abs(expr_feature_plot_group_medians[i, j]))), x, y, gp = gpar(fontsize = 6, fontfamily=plots_props$font_family, col = "white"))
      } else {
        grid.text("", x, y, gp = gpar(fontsize = 6, fontfamily=plots_props$font_family))
      }
    }
    
    mask = expr_feature_plot_group_medians > cutoff
    df_th = expr_feature_plot_group_medians
    df_th[mask] = cutoff
  } else {
    # if no cutoff is need we take the values as there are given
    df_th = expr_feature_plot_group_medians
    maxi = max(df_th[,2:ncol(df_th)])
    mini = min(df_th[,2:ncol(df_th)])
    if (zscored) {
      col_fun = colorRamp2(c(mini, maxi), c("#FFF6D6","#D5A021"))
    } else {
      col_fun = colorRamp2(c(mini, maxi), c("#DFDFB9","#525200"))
    }
    cell_fun = NULL
  }
  
  # set rownames for expr_feature_plot_group_medians 
  rownames(df_th) = df_th$iso_type
  df_th$iso_type = c()
  
  
  # create annotation for the columns
  annotation_colors = list()
  ucolors = unlist(colors[[plot_group]])
  if(is.vector(ucolors) && all(colnames(df_th) %in% names(ucolors))) {
    annotation_colors[[plot_group]] = ucolors
  }
  
  annotation_col_df = data.frame(col=colnames(df_th))
  colnames(annotation_col_df) = c(plot_group)
  
  column_annotation = HeatmapAnnotation(df = annotation_col_df, col = annotation_colors, gp = gpar(row = "white"), simple_anno_size = unit(0.2, "cm"),
                                        annotation_name_side = "left", annotation_name_gp = gpar(fontsize = plots_props$font_size, fontface = "plain", fontfamily=plots_props$font_family),
                                        show_legend=FALSE, show_annotation_name=FALSE
  )
  
  # create annotation for the rows if needed
  if (keep_canonical) {
    row_colours = c()
    for (name in rownames(df_th)) {
      if (name == "0F_0T") {
        row_colours = append(row_colours, feature_colours$canonical)
      } else {
        row_colours = append(row_colours, feature_colours$rest)
      }
    }
    names(row_colours) = rownames(df_th)
    
    annotation_row_df = data.frame(iso_type=rownames(df_th))
    
    row_annotation = rowAnnotation(df = annotation_row_df, col = list("iso_type" = row_colours), gp = gpar(row = "white"), simple_anno_size = unit(0.2, "cm"),
                                   annotation_name_side = "top", annotation_name_gp = gpar(fontsize = plots_props$font_size, fontface = "plain", fontfamily=plots_props$font_family),
                                   show_legend=FALSE, show_annotation_name=FALSE
    )
  } else {
    row_annotation = NULL
  }
  
  # create plot
  heatmap_plot = ComplexHeatmap::Heatmap(df_th,
                                         col = col_fun,
                                         cell_fun = cell_fun,
                                         #rect_gp = gpar(col = "white", lwd = .5),
                                         # 9x6 with horizontal legend at the bottom
                                         #height = unit(2, "cm"), width = unit(6.25, "cm"),
                                         height = unit(4, "cm"), width = unit(6.25, "cm"),
                                         show_row_names = FALSE,
                                         show_column_names = FALSE,
                                         cluster_columns = function(m) hclust(dist(m), method = "complete"), show_column_dend = TRUE, clustering_distance_columns = "euclidean", 
                                         column_dend_reorder = TRUE, #column_split = 3, #column_km = 3,
                                         #cluster_rows = function(m) hclust(dist(m), method = "complete"), show_row_dend = TRUE, clustering_distance_rows = "euclidean", 
                                         #row_dend_reorder = TRUE, 
                                         row_dend_reorder = FALSE, cluster_rows = FALSE,
                                         show_heatmap_legend = FALSE,
                                         left_annotation = row_annotation,
                                         top_annotation = column_annotation,
  )
  
  # # build legend
  if (cutoff != 0) {
    legend_tick_locations = seq(from = -color_bar_max, to = color_bar_max, length.out = 5)
    legend_tick_labels = legend_tick_locations
  } else {
    if (maxi -  mini < 1) {
      length_out = 4
      round_digits = 2
      legend_tick_locations = round(seq(from = mini, to = maxi, length.out = length_out), digits = round_digits)
      legend_tick_labels = round(legend_tick_locations, digits = round_digits)
    } else {
      max_ticks = 4 
      legend_tick_locations = round(seq(from = floor(mini), to = ceiling(maxi)), digits = 0)
      if (length(legend_tick_locations) > max_ticks) {
        # if too many, limit to max_ticks
        legend_tick_locations = round(seq(from = floor(mini), to = ceiling(maxi), length.out=max_ticks), digits = 0)
      }
      legend_tick_labels = round(legend_tick_locations, digits = 0)
    }
  }
  
  #print(legend_tick_locations)
  #print(legend_tick_labels)
  legend_bottom = ComplexHeatmap::Legend(title = color_bar_name_bottom, col_fun = col_fun,
                                         at = legend_tick_locations, labels = legend_tick_labels, direction="horizontal", 
                                         legend_width = unit(3, plots_props$image_units), legend_height=unit(0.5, "cm"),
                                         title_gp = gpar(fontsize = plots_props$font_size, fontface = "bold", fontfamily=plots_props$font_family), labels_gp = gpar(fontsize = plots_props$font_size, fontfamily=plots_props$font_family))
  legend_side = ComplexHeatmap::Legend(title = color_bar_name_side, col_fun = col_fun,
                                       at = legend_tick_locations, labels = legend_tick_labels, 
                                       legend_width = unit(0.75, plots_props$image_units), legend_height=unit(3, "cm"),
                                       title_gp = gpar(fontsize = plots_props$font_size, fontface = "bold", fontfamily=plots_props$font_family), labels_gp = gpar(fontsize = plots_props$font_size, fontfamily=plots_props$font_family))
  
  # saving of the plot
  
  # 9x6 with horizontal legend at the bottom
  #plot_width = plots_props$image_width
  #plot_height = plots_props$image_height
  #png(sprintf("%s_%sx%s.png", file_name, plot_width, plot_height), width = plot_width, height = plot_height, units = plots_props$image_units, res = plots_props$dpi)
  #draw(heatmap_plot, padding = unit(c(1, 0.1, -1, 0.1), "cm"))  # bottom, left, top and right margins
  #draw(legend_bottom, x = unit((6.5), "cm"), y = unit((1), "cm"), just = c("right"))
  #dev.off()
  #svglite(sprintf("%s_%sx%s.svg", file_name, plot_width, plot_height), width = (plot_width / 2.54), height = plot_height / 2.54)
  #draw(heatmap_plot, padding = unit(c(1, 0.1, -1, 0.1), "cm"))  # bottom, left, top and right margins
  #draw(legend_bottom, x = unit((6.5), "cm"), y = unit((1), "cm"), just = c("right"))
  #dev.off()
  
  plot_width = plots_props$image_width
  plot_height = plots_props$image_height
  png(sprintf("%s_%sx%s.png", file_name, plot_width, plot_height), width = plot_width, height = plot_height, units = plots_props$image_units, res = plots_props$dpi)
  draw(heatmap_plot, padding = unit(c(0.1, -1, 0.1, 1), "cm"))  # bottom, left, top and right margins
  draw(legend_side, x = unit((plot_width - 0.1), "cm"), y = unit((plot_height / 2), "cm"), just = c("right", "center"))
  dev.off()
  svglite(sprintf("%s_%sx%s.svg", file_name, plot_width, plot_height), width = (plot_width / 2.54), height = plot_height / 2.54)
  draw(heatmap_plot, padding = unit(c(0.1, -1, 0.1, 1), "cm"))  # bottom, left, top and right margins
  draw(legend_side, x = unit((plot_width - 0.1), "cm"), y = unit((plot_height / 2), "cm"), just = c("right", "center"))
  dev.off()
  
  plot_width = plots_props$image_width
  plot_height = 4
  png(sprintf("%s_%sx%s.png", file_name, plot_width, plot_height), width = plot_width, height = plot_height, units = plots_props$image_units, res = plots_props$dpi)
  draw(heatmap_plot, padding = unit(c(0.1, -1, 0.1, 1), "cm"))  # bottom, left, top and right margins
  draw(legend_side, x = unit((plot_width - 0.1), "cm"), y = unit((plot_height / 2), "cm"), just = c("right", "center"))
  dev.off()
  svglite(sprintf("%s_%sx%s.svg", file_name, plot_width, plot_height), width = (plot_width / 2.54), height = plot_height / 2.54)
  draw(heatmap_plot, padding = unit(c(0.1, -1, 0.1, 1), "cm"))  # bottom, left, top and right margins
  draw(legend_side, x = unit((plot_width - 0.1), "cm"), y = unit((plot_height / 2), "cm"), just = c("right", "center"))
  dev.off()
  
  # table plotting table
  fwrite(df_th, sprintf("%s.csv", file_name_tab), sep = "\t", row.names = TRUE)
  write.xlsx(df_th, sprintf("%s.xlsx", file_name_tab), colNames = TRUE, rowNames = TRUE, append = FALSE)
}


#---------------------------------- Load data ----------------------------------
data_input = snakemake@params$data_input
plot_group = snakemake@params$plot_group
feature = snakemake@params$feature
top_list = snakemake@params$top_list
feature_colours = snakemake@params$colours
feature_col = snakemake@params$feature_col
plots_props = snakemake@params$plots_props
xticks_names = snakemake@params$xticks_names
colors = snakemake@params$colors
results_folder = snakemake@params$results_folder

cutoff = 0

# save rdata
#saveRDS(snakemake, file = "snakemake_isomir_heatmap_plot.rds")

expr = fread(sprintf("%s_%s/%s_%s_quantification_%s.detection_rate_%s_per_%s.csv", data_input$raw_data_folder, data_input$data_sub_set, data_input$rna_class, data_input$feature_filtering, data_input$norm, data_input$detection_rate, data_input$detection_group), sep='\t', header=T)
annot = fread(sprintf("%s/annotation_%s_%s.csv", data_input$annotation_data_folder, data_input$data_sub_set, data_input$feature_filtering), sep='\t', colClasses=c(ID="character")) 


#------------------------------------ Script ----------------------------------- 
# create output folders
output_folder_expr = sprintf("%s/%s_%s/results_%s/figures/%s/heatmap_complex/expr", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, data_input$rna_class)
dir.create(output_folder_expr, recursive=TRUE)
output_folder_zscore_row = sprintf("%s/%s_%s/results_%s/figures/%s/heatmap_complex/zscore=row", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, data_input$rna_class)
dir.create(output_folder_zscore_row, recursive=TRUE)
output_folder_tab = sprintf("%s/%s_%s/results_%s/matrices/%s/heatmap_complex", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, data_input$rna_class)
dir.create(output_folder_tab, recursive=TRUE)

expr_feature = expr[expr[[feature_col]] == feature,]

print(sprintf("There are %s different iso types for %s", dim(expr_feature)[1] - 1, feature ))

if (length(expr_feature$iso_type) != length(unique(expr_feature$iso_type))) {
  for (duplicate in expr_feature$iso_type[!is.na(expr_feature$iso_type)]) {
    if (sum(expr_feature$iso_type[!is.na(expr_feature$iso_type)] == duplicate) > 1) {
      print(sprintf("There are more than one occurances of the same iso type (%s) for %s present", duplicate, feature))
      #print(sprintf("Skip %s", type))
    }
  }
}

for (top in top_list) {
  print(sprintf("top: %s %smiRs", top, data_input$rna_class))
  order_indices = order(rowSums(as.matrix(expr_feature[,3:ncol(expr_feature)])), decreasing = TRUE)
  expr_feature_filtered_top = expr_feature[order_indices][1:top]

  expr_feature_plot_group_medians = data.frame(iso_type=expr_feature_filtered_top$iso_type)
  for (group_1 in unique(annot[[plot_group]])) {
    group_1 = as.character(group_1)
    group_ID = annot[annot[[plot_group]] == group_1,][[data_input$identifier_column]]
    expr_feature_plot_group = expr_feature_filtered_top[, ..group_ID]
    expr_feature_plot_group_medians = merge(expr_feature_plot_group_medians, data.frame(iso_type=expr_feature_filtered_top$iso_type, median=rowMedians(as.matrix(expr_feature_plot_group))), by = "iso_type", sort = FALSE)
  }
  colnames(expr_feature_plot_group_medians) = c("iso_type", unique(annot[[plot_group]]))
  
  #####
  # rpmm 
  zscored = FALSE
  color_bar_name_bottom = sprintf("Expr. (%s, median \ndetection rate = %s%% per %s)", gsub("_norm", "", data_input$norm), gsub("p", "", data_input$detection_rate), gsub("_", " ", data_input$detection_group))
  color_bar_name_side = sprintf("Expr. (%s,\nmedian)", gsub("_norm", "", data_input$norm), gsub("p", "", data_input$detection_rate), gsub("_", " ", data_input$detection_group))
  
  # with the canonical isomiR
  file_name = sprintf("%s/top%s_%smiRs_%s", output_folder_expr, top, data_input$rna_class, feature)
  file_name_tab = sprintf("%s/top%s_%smiRs_expr_%s", output_folder_tab, top, data_input$rna_class, feature)
  keep_canonical = TRUE
  create_heatmap(expr_feature_plot_group_medians, zscored, keep_canonical, plot_group, color_bar_name_bottom, color_bar_name_side, data_input, cutoff, feature_colours, colors, plots_props, file_name, file_name_tab)
    
  # without the canonical isomiR
  file_name = sprintf("%s/top%s_%smiRs_without_canonical_%s", output_folder_expr, top, data_input$rna_class, feature)
  file_name_tab = sprintf("%s/top%s_%smiRs_without_canonical_expr_%s", output_folder_tab, top, data_input$rna_class, feature)
  keep_canonical = FALSE
  expr_feature_plot_group_medians_filtered = expr_feature_plot_group_medians[!(expr_feature_plot_group_medians$iso_type == "0F_0T"),]
  create_heatmap(expr_feature_plot_group_medians_filtered, zscored, keep_canonical, plot_group, color_bar_name_bottom, color_bar_name_side, data_input, cutoff, feature_colours, colors, plots_props, file_name, file_name_tab)
  
  #zsocres for the rows
  zscored = TRUE
  color_bar_name_bottom = sprintf("Stand. expr. (%s, median\ndetection rate = %s%% per %s)", gsub("_norm", "", data_input$norm), gsub("p", "", data_input$detection_rate), gsub("_", " ", data_input$detection_group))
  color_bar_name_side = sprintf("Stand. expr.\n(%s,\nmedian)", gsub("_norm", "", data_input$norm), gsub("p", "", data_input$detection_rate), gsub("_", " ", data_input$detection_group))
  
  # with the canonical isomiR
  expr_feature_plot_group_medians_zscores = expr_feature_plot_group_medians
  expr_feature_plot_group_medians_zscores[2:ncol(expr_feature_plot_group_medians)] = scale(as.matrix(expr_feature_plot_group_medians[2:ncol(expr_feature_plot_group_medians)]), TRUE, FALSE)
  file_name = sprintf("%s/top%s_%smiRs_zscores=row_%s", output_folder_zscore_row, top, data_input$rna_class, feature)
  file_name_tab = sprintf("%s/top%s_%smiRs_zscores=row_%s", output_folder_tab, top, data_input$rna_class, feature)
  keep_canonical = TRUE
  create_heatmap(expr_feature_plot_group_medians_zscores, zscored, keep_canonical, plot_group, color_bar_name_bottom, color_bar_name_side, data_input, cutoff, feature_colours, colors, plots_props, file_name, file_name_tab)

  # without the canonical isomiR
  expr_feature_plot_group_medians_zscores_filtered = expr_feature_plot_group_medians_zscores[!(expr_feature_plot_group_medians_zscores$iso_type == "0F_0T"),]
  file_name = sprintf("%s/top%s_%smiRs_without_canonical_zscores=row_%s", output_folder_zscore_row, top, data_input$rna_class, feature)
  file_name_tab = sprintf("%s/top%s_%smiRs_without_canonical_zscores=row_%s", output_folder_tab, top, data_input$rna_class, feature)
  keep_canonical = FALSE
  create_heatmap(expr_feature_plot_group_medians_zscores_filtered, zscored, keep_canonical, plot_group, color_bar_name_bottom, color_bar_name_side, data_input, cutoff, feature_colours, colors, plots_props, file_name, file_name_tab)
  
  #####
  # log10 rpmm
  expr_feature_plot_group_medians_log10 = expr_feature_plot_group_medians
  expr_feature_plot_group_medians_log10[2:ncol(expr_feature_plot_group_medians_log10)] = log10(as.matrix(expr_feature_plot_group_medians_log10[2:ncol(expr_feature_plot_group_medians_log10)]) + 1)
  
  zscored = FALSE
  color_bar_name_bottom = sprintf("Expr. (%s, median, log10,\ndetection rate = %s%% per %s)", gsub("_norm", "", data_input$norm), gsub("p", "", data_input$detection_rate), gsub("_", " ", data_input$detection_group))
  color_bar_name_side = sprintf("Expr. (%s,\nmedian,\nlog10)", gsub("_norm", "", data_input$norm), gsub("p", "", data_input$detection_rate), gsub("_", " ", data_input$detection_group))
  
  # with the canonical isomiR
  file_name = sprintf("%s/top%s_%smiRs_log10_%s", output_folder_expr, top, data_input$rna_class, feature)
  file_name_tab = sprintf("%s/top%s_%smiRs_expr_log10_%s", output_folder_tab, top, data_input$rna_class, feature)
  keep_canonical = TRUE
  create_heatmap(expr_feature_plot_group_medians_log10, zscored, keep_canonical, plot_group, color_bar_name_bottom, color_bar_name_side, data_input, cutoff, feature_colours, colors, plots_props, file_name, file_name_tab)
  
  # without the canonical isomiR
  file_name = sprintf("%s/top%s_%smiRs_without_canonical_log10_%s", output_folder_expr, top, data_input$rna_class, feature)
  file_name_tab = sprintf("%s/top%s_%smiRs_without_canonical_expr_log10_%s", output_folder_tab, top, data_input$rna_class, feature)
  keep_canonical = FALSE
  expr_feature_plot_group_medians_log10_filtered = expr_feature_plot_group_medians_log10[!(expr_feature_plot_group_medians_log10$iso_type == "0F_0T"),]
  create_heatmap(expr_feature_plot_group_medians_log10_filtered, zscored, keep_canonical, plot_group, color_bar_name_bottom, color_bar_name_side, data_input, cutoff, feature_colours, colors, plots_props, file_name, file_name_tab)
  
  # zscores for the rows
  zscored = TRUE
  color_bar_name_bottom = sprintf("Stand. expr. (%s, median, log10,\ndetection rate = %s%% per %s)", gsub("_norm", "", data_input$norm), gsub("p", "", data_input$detection_rate), gsub("_", " ", data_input$detection_group))
  color_bar_name_side = sprintf("Stand. expr.\n(%s,\nmedian,\nlog10)", gsub("_norm", "", data_input$norm), gsub("p", "", data_input$detection_rate), gsub("_", " ", data_input$detection_group))
  
  # with the canonical isomiR
  expr_feature_plot_group_medians_log10_zscores = expr_feature_plot_group_medians_log10
  expr_feature_plot_group_medians_log10_zscores[2:ncol(expr_feature_plot_group_medians_log10)] = scale(as.matrix(expr_feature_plot_group_medians_log10[2:ncol(expr_feature_plot_group_medians_log10)]), TRUE, FALSE)
  file_name = sprintf("%s/top%s_%smiRs_log10_zscores=row_%s", output_folder_zscore_row, top, data_input$rna_class, feature)
  file_name_tab = sprintf("%s/top%s_%smiRs_log10_zscores=row_%s", output_folder_tab, top, data_input$rna_class, feature)
  keep_canonical = TRUE
  create_heatmap(expr_feature_plot_group_medians_log10_zscores, zscored, keep_canonical, plot_group, color_bar_name_bottom, color_bar_name_side, data_input, cutoff, feature_colours, colors, plots_props, file_name, file_name_tab)
  
  # without the canonical isomiR
  expr_feature_plot_group_medians_log10_zscores_filtered = expr_feature_plot_group_medians_log10_zscores[!(expr_feature_plot_group_medians_log10_zscores$iso_type == "0F_0T"),]
  file_name = sprintf("%s/top%s_%smiRs_without_canonical_log10_zscores=row_%s", output_folder_zscore_row, top, data_input$rna_class, feature)
  file_name_tab = sprintf("%s/top%s_%smiRs_without_canonical_log10_zscores=row_%s", output_folder_tab, top, data_input$rna_class, feature)
  keep_canonical = FALSE
  create_heatmap(expr_feature_plot_group_medians_log10_zscores_filtered, zscored, keep_canonical, plot_group, color_bar_name_bottom, color_bar_name_side, data_input, cutoff, feature_colours, colors, plots_props, file_name, file_name_tab)
}


#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])
