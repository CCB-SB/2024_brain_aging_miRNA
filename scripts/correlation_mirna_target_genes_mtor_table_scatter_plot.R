suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(Hmisc))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(svglite))
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(foreach))

#snakemake = readRDS("snakemake_correlation_mirna_target_genes_mtor_table_scatter_plot.rds")

set.seed(snakemake@params$parameters_props$set_seed)

print("correlation_mirna_target_genes_mtor_table_scatter_plot")


#---------------------------------- Functions ----------------------------------
plot_heatmap_annot = function(df, pvalues_df, row_names, data_input, xticks_names, plots_props, prop, time, m, annotation_row, colors, th, sig_lvl) {
  
  df = t(df)
  #df[abs(df) < th] = 0
  #df[df <= -th] = -0.5
  #df[df >= th] = 0.5
  
  pvalues_df = t(pvalues_df)
  #pvalues_df = pvalues_df[, colSums(df != 0) > 0]
  #row_names = row_names[colSums(df != 0) > 0]
  #df = df[, colSums(df != 0) > 0]
  #row_names = row_names[colSums(is.na(df)) != nrow(df)]
  #df = df[, colSums(is.na(df)) != nrow(df)]
  
  #tmp = !(row_names %in% feature_list[[data_input$feature_column]])
  #row_names[tmp] = ""
  #df_label = matrix("|", nrow(df), ncol(df))
  #df_label[,tmp] = ""
  
  if (m == "spearman") {
    color_bar_name = sprintf("Spearman correlation (%s)", time) 
  } else if (m == "pearson") {
    color_bar_name = sprintf("Pearson correlation (%s)", time) 
  } else {
    print("Wrong corr method selected")
    stop()
  }
  
  #cmap = colorRamp2(c(-1, 0, 1), hcl_palette = "Green-Brown", reverse = TRUE)
  #col_val = 0.75
  #col_fun = c(cmap(-col_val), cmap(0), cmap(col_val))
  #col_legend = c(cmap(-col_val), cmap(col_val))
  #col_fun = c(colors$direction$neg, colors$direction$light, colors$direction$pos)
  color_bar_max = 1
  scaling_factor = 5
  colorbar_colors = c(colors$direction$down, colors$direction$down_light, colors$direction$light, colors$direction$up_light, colors$direction$up)
  col_fun = colorRamp2(c(-color_bar_max, -(color_bar_max/scaling_factor), 0, color_bar_max/scaling_factor, color_bar_max), colorbar_colors)
  col_legend = c(colors$direction$neg, colors$direction$pos)
  #names(col_legend) = c(sprintf("negative (R ≤ -%s)", th), sprintf("positive (R ≥ %s)", th))
  
  #df_label = matrix("", nrow(df), ncol(df))
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (abs(df[i, j]) >= th) {
      grid.rect(x, y, width, height, gp = gpar(fill=fill, col = "black", lwd = 1))  # Black cell border
    }
    if (is.na(pvalues_df[i, j])) {
    }
    else if ((pvalues_df)[i, j] < 0.05 & (df[i,j] != 0)){ 
      grid.text("*", x, y, gp = gpar(fontsize = 2, col = "black"))
    #} else if (pvalues_df[i, j] < 0.01){ 
    #    grid.text("**", x, y, gp = gpar(fontsize = 4))
    #} else if (pvalues_df[i, j] < 0.001){ 
    #    grid.text("***", x, y, gp = gpar(fontsize = 4))
    }
  }
  
  #annotation_row = data.frame(dummy=1:nrow(df))
  #annotation_row[[prop]] = as.character(annot[[prop]][match(rownames(df), annot[[prop]])])
  #colnames(annotation_row) = rownames(df)
  #annotation_row$dummy = NULL
  
  annotation_colors = list()
  ucolors = unlist(colors[[prop]])
  if(is.vector(ucolors) && all(annotation_row %in% names(ucolors))) {
    annotation_colors[[prop]] = ucolors
  }
  
  # if(!is.null(additional_props)){
  #   if(is.vector(additional_props)){
  #     for(p in additional_props){
  #       annotation_col[[p]] = as.character(annot[[p]][match(colnames(df), annot[[data_input$identifier_column]])])
  #       ucolors <- unlist(colors[[p]])
  #       if(is.vector(ucolors) && all(annotation_col[[p]] %in% names(ucolors))) {
  #         annotation_colors[[p]] = ucolors
  #       }
  #     }
  #   } else {
  #     annotation_col[[additional_props]] = as.character(annot[[additional_props]][match(colnames(df), annot[[data_input$identifier_column]])])
  #     annotation_colors[[additional_props]] = unlist(colors[[additional_props]])
  #   }
  # }
  
  # annotation_row = row_names
  # annotation_colors = colors[[prop]]
  
  # sort the annotation according to the dataframe
  annotation_row = annotation_row[rownames(df), drop=FALSE]
  
  colnames(df) = row_names
  rownames(df) = unlist(xticks_names[[prop]][rownames(df)])
  rownames(annotation_row) = unlist(xticks_names[[prop]][rownames(annotation_row)])
  
  # print(annotation_row)
  # print(annotation_colors)
  column_ha = rowAnnotation(df = annotation_row, col = annotation_colors, gp = gpar(row = "white"), simple_anno_size = unit(0.2, "cm"), 
                            annotation_name_side = "top", annotation_name_gp = gpar(fontsize = plots_props$font_size, fontface = "bold", fontfamily=plots_props$font_family),
                            show_legend=FALSE, show_annotation_name=FALSE)
  #cell_fun = function(j, i, x, y, width, height, fill) {grid.text(df_label[i, j], x, y, gp = gpar(fontsize = 6))}
  
  p = ComplexHeatmap::Heatmap(df, col = col_fun,
                              cell_fun = cell_fun,
                              #rect_gp = gpar(col = "white", lwd = .5),
                              #height = unit(0.3, "cm") * nrow(t(df)), width = unit(0.3, "cm") * ncol(t(df)),
                              #height = unit(plots_props$image_height  - 2, "cm"), width = unit(plots_props$image_width - 4, "cm"),
                              height = unit(2, "cm"), width = unit(6.25, "cm"),
                              #show_column_names = FALSE,
                              row_names_side = "left", row_names_gp = gpar(fontsize = 6, fontfamily=plots_props$font_family), 
                              column_names_side = "bottom", column_names_gp = gpar(fontsize = 6, fontfamily=plots_props$font_family),
                              cluster_columns = function(m) hclust(dist(m), method = "complete"), show_column_dend = TRUE, clustering_distance_columns = "euclidean", 
                              #cluster_rows = function(m) hclust(dist(m), method = "complete"), show_row_dend = TRUE, clustering_distance_rows = "euclidean", 
                              column_dend_reorder = TRUE, #cluster_columns = FALSE,
                              row_dend_reorder = FALSE, cluster_rows = FALSE,
                              show_heatmap_legend = FALSE,
                              #left_annotation = column_ha
  )
  #legend_ticks = seq(from = min(t(df)), to = max(t(df)), length.out = 5)
  #legend_tick_labels = round(legend_ticks, digits = 2)
  
  #legend = ComplexHeatmap::Legend(title = color_bar_name, col_fun = col_fun,
  #                                at = legend_tick_labels, direction="horizontal", #col_fun = circlize::colorRamp2(legend_ticks, viridis(5))
  #                                title_gp = gpar(fontsize = plots_props$font_size, fontface = "bold", fontfamily=plots_props$font_family), labels_gp = gpar(fontsize = plots_props$font_size, fontfamily=plots_props$font_family),
  #                                legend_width=unit(4, "cm"))
  
  #legend = ComplexHeatmap::Legend(title = color_bar_name, labels = names(col_legend), legend_gp = gpar(fill = col_legend),
  #                                direction="horizontal", ncol = 2,
  #                                title_gp = gpar(fontsize = plots_props$font_size, fontface = "bold", fontfamily=plots_props$font_family), labels_gp = gpar(fontsize = plots_props$font_size, fontfamily=plots_props$font_family),
  #                                legend_width=unit(4, "cm"))
  
  legend_ticks = seq(from = -1, to = 1, length.out = 5)
  legend_tick_labels = round(legend_ticks, digits = 2)
  
  legend = ComplexHeatmap::Legend(title = color_bar_name, col_fun = col_fun,
                                  at = legend_tick_labels, 
                                  direction="horizontal",
                                  title_gp = gpar(fontsize = plots_props$font_size, fontface = "bold", fontfamily=plots_props$font_family), 
                                  labels_gp = gpar(fontsize = plots_props$font_size, fontfamily=plots_props$font_family),
                                  legend_width=unit(3, "cm"), legend_height=unit(3, "cm")) #title_position = "topcenter")
  
  #p = pheatmap(t(df), color = viridis(100), border_color = "white", show_colnames = TRUE, show_rownames = TRUE,
  #             fontsize = plots_props$font_size, fontsize_row = 6, fontsize_col = plots_props$font_size, 
  #             cellwidth = 8, cellheight = 8,
  #             cluster_cols = FALSE, clustering_distance_cols = "euclidean", cluster_rows = TRUE, clustering_distance_rows = "euclidean", treeheight_row = 0
  #)
  
  return(list("heatmap" = p, "legend" = legend))
}

plot_heatmap_annot_flipped = function(df, pvalues_df, row_names, data_input, xticks_names, plots_props, prop, time, m, annotation_row, colors, th, sig_lvl) {
  
  #df = t(df)
  #df[abs(df) < th] = 0
  #df[df <= -th] = -0.5
  #df[df >= th] = 0.5
  
  #pvalues_df = t(pvalues_df)
  #pvalues_df = pvalues_df[, colSums(df != 0) > 0]
  #row_names = row_names[colSums(df != 0) > 0]
  #df = df[, colSums(df != 0) > 0]
  #row_names = row_names[colSums(is.na(df)) != nrow(df)]
  #df = df[, colSums(is.na(df)) != nrow(df)]
  
  #tmp = !(row_names %in% feature_list[[data_input$feature_column]])
  #row_names[tmp] = ""
  #df_label = matrix("|", nrow(df), ncol(df))
  #df_label[,tmp] = ""
  
  if (!time == "") {
    if (m == "spearman") {
      color_bar_name = sprintf("Spearman\ncorr. (%s)", time) 
    } else if (m == "pearson") {
      color_bar_name = sprintf("Pearson\ncorr. (%s)", time) 
    } else {
      print("Wrong corr method selected")
      stop()
    }
  } else {
    if (m == "spearman") {
      color_bar_name = "Spearman\ncorr."
    } else if (m == "pearson") {
      color_bar_name = "Pearson\ncorr."
    } else {
      print("Wrong corr method selected")
      stop()
    }
  }
  
  #cmap = colorRamp2(c(-1, 0, 1), hcl_palette = "Green-Brown", reverse = TRUE)
  #col_val = 0.75
  #col_fun = c(cmap(-col_val), cmap(0), cmap(col_val))
  #col_legend = c(cmap(-col_val), cmap(col_val))
  #col_fun = c(colors$direction$neg, colors$direction$light, colors$direction$pos)
  color_bar_max = 1
  scaling_factor = 5
  colorbar_colors = c(colors$direction$down, colors$direction$down_light, colors$direction$light, colors$direction$up_light, colors$direction$up)
  col_fun = colorRamp2(c(-color_bar_max, -(color_bar_max/scaling_factor), 0, color_bar_max/scaling_factor, color_bar_max), colorbar_colors)
  col_legend = c(colors$direction$neg, colors$direction$pos)
  #names(col_legend) = c(sprintf("negative (R ≤ -%s)", th), sprintf("positive (R ≥ %s)", th))
  
  #df_label = matrix("", nrow(df), ncol(df))
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (abs(df[i, j]) >= th) {
      grid.rect(x, y, width, height, gp = gpar(fill=fill, col = "black", lwd = 1))  # Black cell border
    }
    if (is.na(pvalues_df[i, j])) {
    }
    else if ((pvalues_df)[i, j] < 0.05 & (df[i,j] != 0)){ 
      grid.text("✱", x, y, gp = gpar(fontsize = 5, col = "black"), just = "center")
      #} else if (pvalues_df[i, j] < 0.01){ 
      #    grid.text("**", x, y, gp = gpar(fontsize = 4))
      #} else if (pvalues_df[i, j] < 0.001){ 
      #    grid.text("***", x, y, gp = gpar(fontsize = 4))
    }
  }
  
  #annotation_row = data.frame(dummy=1:nrow(df))
  #annotation_row[[prop]] = as.character(annot[[prop]][match(rownames(df), annot[[prop]])])
  #colnames(annotation_row) = rownames(df)
  #annotation_row$dummy = NULL
  
  annotation_colors = list()
  ucolors = unlist(colors[[prop]])
  if(is.vector(ucolors) && all(annotation_row %in% names(ucolors))) {
    annotation_colors[[prop]] = ucolors
  }
  
  # if(!is.null(additional_props)){
  #   if(is.vector(additional_props)){
  #     for(p in additional_props){
  #       annotation_col[[p]] = as.character(annot[[p]][match(colnames(df), annot[[data_input$identifier_column]])])
  #       ucolors <- unlist(colors[[p]])
  #       if(is.vector(ucolors) && all(annotation_col[[p]] %in% names(ucolors))) {
  #         annotation_colors[[p]] = ucolors
  #       }
  #     }
  #   } else {
  #     annotation_col[[additional_props]] = as.character(annot[[additional_props]][match(colnames(df), annot[[data_input$identifier_column]])])
  #     annotation_colors[[additional_props]] = unlist(colors[[additional_props]])
  #   }
  # }
  
  # annotation_row = row_names
  # annotation_colors = colors[[prop]]
  
  # sort the annotation according to the dataframe
  annotation_row = annotation_row[order(match(unlist(annotation_row), colnames(df)))]
  #annotation_row = annotation_row[colnames(df), drop=FALSE]
  names(annotation_row) = unname(unlist(xticks_names[[prop]][annotation_row]))
  
  rownames(df) = row_names
  colnames(df) = unname(unlist(xticks_names[[prop]][colnames(df)]))
  
  # print(annotation_row)
  # print(annotation_colors)
  column_ha = HeatmapAnnotation(df = annotation_row, col = annotation_colors, gp = gpar(row = "white"), simple_anno_size = unit(0.2, "cm"),
                                annotation_name_side = "left", annotation_name_gp = gpar(fontsize = plots_props$font_size, fontface = "bold", fontfamily=plots_props$font_family),
                                show_legend=FALSE, show_annotation_name=FALSE)
  #cell_fun = function(j, i, x, y, width, height, fill) {grid.text(df_label[i, j], x, y, gp = gpar(fontsize = 6))}
  
  p = ComplexHeatmap::Heatmap(df, col = col_fun,
                              cell_fun = cell_fun,
                              #rect_gp = gpar(col = "white", lwd = .5),
                              #height = unit(0.3, "cm") * nrow(t(df)), width = unit(0.3, "cm") * ncol(t(df)),
                              #height = unit(plots_props$image_height  - 2, "cm"), width = unit(plots_props$image_width - 4, "cm"),
                              height = unit(11, "cm"), width = unit(12, "cm"),
                              show_column_names = FALSE,
                              row_names_side = "left", row_names_gp = gpar(fontsize = 6, fontfamily=plots_props$font_family), 
                              #column_names_side = "bottom", column_names_gp = gpar(fontsize = 6, fontfamily=plots_props$font_family),
                              #cluster_columns = function(m) hclust(dist(m), method = "complete"), show_column_dend = TRUE, clustering_distance_columns = "euclidean", 
                              cluster_rows = function(m) hclust(dist(m), method = "complete"), show_row_dend = FALSE, clustering_distance_rows = "euclidean", 
                              column_dend_reorder = FALSE, cluster_columns = FALSE,
                              row_dend_reorder = TRUE, #cluster_rows = FALSE,
                              show_heatmap_legend = FALSE,
                              #left_annotation = column_ha
  )
  #legend_ticks = seq(from = min(t(df)), to = max(t(df)), length.out = 5)
  #legend_tick_labels = round(legend_ticks, digits = 2)
  
  #legend = ComplexHeatmap::Legend(title = color_bar_name, col_fun = col_fun,
  #                                at = legend_tick_labels, direction="horizontal", #col_fun = circlize::colorRamp2(legend_ticks, viridis(5))
  #                                title_gp = gpar(fontsize = plots_props$font_size, fontface = "bold", fontfamily=plots_props$font_family), labels_gp = gpar(fontsize = plots_props$font_size, fontfamily=plots_props$font_family),
  #                                legend_width=unit(4, "cm"))
  
  #legend = ComplexHeatmap::Legend(title = color_bar_name, labels = names(col_legend), legend_gp = gpar(fill = col_legend),
  #                                direction="horizontal", ncol = 2,
  #                                title_gp = gpar(fontsize = plots_props$font_size, fontface = "bold", fontfamily=plots_props$font_family), labels_gp = gpar(fontsize = plots_props$font_size, fontfamily=plots_props$font_family),
  #                                legend_width=unit(4, "cm"))
  
  legend_ticks = seq(from = -1, to = 1, length.out = 5)
  legend_tick_labels = round(legend_ticks, digits = 2)
  
  legend = ComplexHeatmap::Legend(title = color_bar_name, col_fun = col_fun,
                                  at = legend_tick_labels, 
                                  #direction="horizontal",
                                  title_gp = gpar(fontsize = plots_props$font_size, fontface = "bold", fontfamily=plots_props$font_family), 
                                  labels_gp = gpar(fontsize = plots_props$font_size, fontfamily=plots_props$font_family),
                                  legend_width=unit(3, "cm"), legend_height=unit(3, "cm")) #title_position = "topcenter")
  
  #p = pheatmap(t(df), color = viridis(100), border_color = "white", show_colnames = TRUE, show_rownames = TRUE,
  #             fontsize = plots_props$font_size, fontsize_row = 6, fontsize_col = plots_props$font_size, 
  #             cellwidth = 8, cellheight = 8,
  #             cluster_cols = FALSE, clustering_distance_cols = "euclidean", cluster_rows = TRUE, clustering_distance_rows = "euclidean", treeheight_row = 0
  #)
  
  return(list("heatmap" = p, "legend" = legend))
}

plot_heatmap_annot_flipped_tissue_flipped = function(df, pvalues_df, row_names, data_input, xticks_names, plots_props, prop, time, m, annotation_row, colors, th, sig_lvl) {
  
  #df = t(df)
  #df[abs(df) < th] = 0
  #df[df <= -th] = -0.5
  #df[df >= th] = 0.5
  
  #pvalues_df = t(pvalues_df)
  #pvalues_df = pvalues_df[, colSums(df != 0) > 0]
  #row_names = row_names[colSums(df != 0) > 0]
  #df = df[, colSums(df != 0) > 0]
  #row_names = row_names[colSums(is.na(df)) != nrow(df)]
  #df = df[, colSums(is.na(df)) != nrow(df)]
  
  #tmp = !(row_names %in% feature_list[[data_input$feature_column]])
  #row_names[tmp] = ""
  #df_label = matrix("|", nrow(df), ncol(df))
  #df_label[,tmp] = ""
  
  if (!time == "") {
    if (m == "spearman") {
      color_bar_name = sprintf("Spearman\ncorr. (%s)", time) 
    } else if (m == "pearson") {
      color_bar_name = sprintf("Pearson\ncorr. (%s)", time) 
    } else {
      print("Wrong corr method selected")
      stop()
    }
  } else {
    if (m == "spearman") {
      color_bar_name = "Spearman\ncorr."
    } else if (m == "pearson") {
      color_bar_name = "Pearson\ncorr."
    } else {
      print("Wrong corr method selected")
      stop()
    }
  }
  
  #cmap = colorRamp2(c(-1, 0, 1), hcl_palette = "Green-Brown", reverse = TRUE)
  #col_val = 0.75
  #col_fun = c(cmap(-col_val), cmap(0), cmap(col_val))
  #col_legend = c(cmap(-col_val), cmap(col_val))
  #col_fun = c(colors$direction$neg, colors$direction$light, colors$direction$pos)
  color_bar_max = 1
  scaling_factor = 5
  colorbar_colors = c(colors$direction$down, colors$direction$down_light, colors$direction$light, colors$direction$up_light, colors$direction$up)
  col_fun = colorRamp2(c(-color_bar_max, -(color_bar_max/scaling_factor), 0, color_bar_max/scaling_factor, color_bar_max), colorbar_colors)
  col_legend = c(colors$direction$neg, colors$direction$pos)
  #names(col_legend) = c(sprintf("negative (R ≤ -%s)", th), sprintf("positive (R ≥ %s)", th))
  
  #df_label = matrix("", nrow(df), ncol(df))
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (abs(df[i, j]) >= th) {
      grid.rect(x, y, width, height, gp = gpar(fill=fill, col = "black", lwd = 1))  # Black cell border
    }
    if (is.na(pvalues_df[i, j])) {
    }
    else if ((pvalues_df)[i, j] < 0.05 & (df[i,j] != 0)){ 
      grid.text("✱", x, y, gp = gpar(fontsize = 5, col = "black"), just = "center")
      #} else if (pvalues_df[i, j] < 0.01){ 
      #    grid.text("**", x, y, gp = gpar(fontsize = 4))
      #} else if (pvalues_df[i, j] < 0.001){ 
      #    grid.text("***", x, y, gp = gpar(fontsize = 4))
    }
  }
  
  #annotation_row = data.frame(dummy=1:nrow(df))
  #annotation_row[[prop]] = as.character(annot[[prop]][match(rownames(df), annot[[prop]])])
  #colnames(annotation_row) = rownames(df)
  #annotation_row$dummy = NULL
  
  annotation_colors = list()
  ucolors = unlist(colors[[prop]])
  if(is.vector(ucolors) && all(annotation_row %in% names(ucolors))) {
    annotation_colors[[prop]] = ucolors
  }
  
  # if(!is.null(additional_props)){
  #   if(is.vector(additional_props)){
  #     for(p in additional_props){
  #       annotation_col[[p]] = as.character(annot[[p]][match(colnames(df), annot[[data_input$identifier_column]])])
  #       ucolors <- unlist(colors[[p]])
  #       if(is.vector(ucolors) && all(annotation_col[[p]] %in% names(ucolors))) {
  #         annotation_colors[[p]] = ucolors
  #       }
  #     }
  #   } else {
  #     annotation_col[[additional_props]] = as.character(annot[[additional_props]][match(colnames(df), annot[[data_input$identifier_column]])])
  #     annotation_colors[[additional_props]] = unlist(colors[[additional_props]])
  #   }
  # }
  
  # annotation_row = row_names
  # annotation_colors = colors[[prop]]
  
  # sort the annotation according to the dataframe
  annotation_row = annotation_row[order(match(unlist(annotation_row), colnames(df)))]
  #annotation_row = annotation_row[colnames(df), drop=FALSE]
  names(annotation_row) = unname(unlist(xticks_names[[prop]][annotation_row]))
  
  rownames(df) = row_names
  colnames(df) = unname(unlist(xticks_names[[prop]][colnames(df)]))
  
  # print(annotation_row)
  # print(annotation_colors)
  column_ha = HeatmapAnnotation(df = annotation_row, col = annotation_colors, gp = gpar(row = "white"), simple_anno_size = unit(0.2, "cm"),
                                annotation_name_side = "left", annotation_name_gp = gpar(fontsize = plots_props$font_size, fontface = "bold", fontfamily=plots_props$font_family),
                                show_legend=FALSE, show_annotation_name=FALSE)
  #cell_fun = function(j, i, x, y, width, height, fill) {grid.text(df_label[i, j], x, y, gp = gpar(fontsize = 6))}
  
  p = ComplexHeatmap::Heatmap(df, col = col_fun,
                              cell_fun = cell_fun,
                              #rect_gp = gpar(col = "white", lwd = .5),
                              #height = unit(0.3, "cm") * nrow(t(df)), width = unit(0.3, "cm") * ncol(t(df)),
                              #height = unit(plots_props$image_height  - 2, "cm"), width = unit(plots_props$image_width - 4, "cm"),
                              height = unit(11, "cm"), width = unit(7, "cm"),
                              show_column_names = FALSE,
                              row_names_side = "left", row_names_gp = gpar(fontsize = 6, fontfamily=plots_props$font_family), 
                              #column_names_side = "bottom", column_names_gp = gpar(fontsize = 6, fontfamily=plots_props$font_family),
                              #cluster_columns = function(m) hclust(dist(m), method = "complete"), show_column_dend = TRUE, clustering_distance_columns = "euclidean", 
                              cluster_rows = function(m) hclust(dist(m), method = "complete"), show_row_dend = FALSE, clustering_distance_rows = "euclidean", 
                              column_dend_reorder = FALSE, cluster_columns = FALSE,
                              row_dend_reorder = TRUE, #cluster_rows = FALSE,
                              show_heatmap_legend = FALSE,
                              #left_annotation = column_ha
  )
  #legend_ticks = seq(from = min(t(df)), to = max(t(df)), length.out = 5)
  #legend_tick_labels = round(legend_ticks, digits = 2)
  
  #legend = ComplexHeatmap::Legend(title = color_bar_name, col_fun = col_fun,
  #                                at = legend_tick_labels, direction="horizontal", #col_fun = circlize::colorRamp2(legend_ticks, viridis(5))
  #                                title_gp = gpar(fontsize = plots_props$font_size, fontface = "bold", fontfamily=plots_props$font_family), labels_gp = gpar(fontsize = plots_props$font_size, fontfamily=plots_props$font_family),
  #                                legend_width=unit(4, "cm"))
  
  #legend = ComplexHeatmap::Legend(title = color_bar_name, labels = names(col_legend), legend_gp = gpar(fill = col_legend),
  #                                direction="horizontal", ncol = 2,
  #                                title_gp = gpar(fontsize = plots_props$font_size, fontface = "bold", fontfamily=plots_props$font_family), labels_gp = gpar(fontsize = plots_props$font_size, fontfamily=plots_props$font_family),
  #                                legend_width=unit(4, "cm"))
  
  legend_ticks = seq(from = -1, to = 1, length.out = 5)
  legend_tick_labels = round(legend_ticks, digits = 2)
  
  legend = ComplexHeatmap::Legend(title = color_bar_name, col_fun = col_fun,
                                  at = legend_tick_labels, 
                                  #direction="horizontal",
                                  title_gp = gpar(fontsize = plots_props$font_size, fontface = "bold", fontfamily=plots_props$font_family), 
                                  labels_gp = gpar(fontsize = plots_props$font_size, fontfamily=plots_props$font_family),
                                  legend_width=unit(3, "cm"), legend_height=unit(3, "cm")) #title_position = "topcenter")
  
  #p = pheatmap(t(df), color = viridis(100), border_color = "white", show_colnames = TRUE, show_rownames = TRUE,
  #             fontsize = plots_props$font_size, fontsize_row = 6, fontsize_col = plots_props$font_size, 
  #             cellwidth = 8, cellheight = 8,
  #             cluster_cols = FALSE, clustering_distance_cols = "euclidean", cluster_rows = TRUE, clustering_distance_rows = "euclidean", treeheight_row = 0
  #)
  
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
              height = unit(4, "cm"), width = unit(1.5, "cm"),
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

row_name_squares_top = function(tissue_list, xticks_names, tissue, colors, plots_props) {
  # Create a sample matrix
  data =- matrix(1:length(tissue_list), nrow = 1, ncol = length(tissue_list))
  colnames(data) = tissue_list
  rownames(data) = "Square"
  
  # Create a function to draw text inside squares
  draw_text = function(j, i, x, y, width, height, fill) { 
    wording = xticks_names[[tissue]][[colnames(data)[j]]]
    text_colour = unlist(colors[[sprintf("%s_text", tissue)]])[[colnames(data)[j]]]
    rect_colour = unlist(colors[[tissue]])[[colnames(data)[j]]]
    grid.roundrect(x=x, y=y, width=width, height=height, r = unit(0.1, "cm"), gp = gpar(fill = rect_colour, col = "white", lwd = 2))
    grid.text(wording, x = x, y = y, gp = gpar(col = text_colour, fontsize = 6)) #fontfamily=plots_props$font_family
  }
  
  #col_fun  = structure(unname(unlist(colors[[tissue]])), names = rownames(data))
  
  # Create the heatmap
  p = Heatmap(data,
              #col = col_fun,
              height = unit(0.3, "cm"), width = unit(12, "cm"),
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

row_name_squares_top_flipped = function(tissue_list, xticks_names, tissue, colors, plots_props) {
  # Create a sample matrix
  data =- matrix(1:length(tissue_list), nrow = 1, ncol = length(tissue_list))
  colnames(data) = tissue_list
  rownames(data) = "Square"
  
  # Create a function to draw text inside squares
  draw_text = function(j, i, x, y, width, height, fill) { 
    wording = xticks_names[[tissue]][[colnames(data)[j]]]
    text_colour = unlist(colors[[sprintf("%s_text", tissue)]])[[colnames(data)[j]]]
    rect_colour = unlist(colors[[tissue]])[[colnames(data)[j]]]
    grid.roundrect(x=x, y=y, width=width, height=height, r = unit(0.1, "cm"), gp = gpar(fill = rect_colour, col = "white", lwd = 2))
    grid.text(wording, x = x, y = y, gp = gpar(col = text_colour, fontsize = 6), rot = 90) #fontfamily=plots_props$font_family
  }
  
  #col_fun  = structure(unname(unlist(colors[[tissue]])), names = rownames(data))
  
  # Create the heatmap
  p = Heatmap(data,
              #col = col_fun,
              height = unit(1.5, "cm"), width = unit(7, "cm"),
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


# ---------------------- Load data ----------------------
data_input = snakemake@params$data_input
feature = snakemake@params$feature
gene_list_cas_folder_path = snakemake@params$gene_list_cas_folder_path
gene_list_target_genes_martin_file = snakemake@params$gene_list_target_genes
gene_list_target_genes_folder_path = snakemake@params$gene_list_target_genes_folder_path
gene_list_mtor = snakemake@params$gene_list_mtor
split_prop = snakemake@params$split_prop
corr_prop = snakemake@params$corr_prop
corr_params = snakemake@params$corr_params
manual_plots_parameters = snakemake@params$manual_plots_parameters
plots_props = snakemake@params$plots_props
xticks_names = snakemake@params$xticks_names
colors = snakemake@params$colors
results_folder = snakemake@params$results_folder

# save rdata
#saveRDS(snakemake, file = "snakemake_correlation_mirna_target_genes_mtor_table_scatter_plot.rds")
#stop()

# load annot
annot = fread(sprintf("%s/annotation_%s_%s.csv", data_input$annotation_data_folder, data_input$data_sub_set, data_input$feature_filtering), sep='\t', colClasses=c(ID="character"))

# load mrna expr
tbl = read.csv(sprintf("%s_%s/mRNA_normalized_counts.csv", data_input$raw_data_folder, data_input$data_sub_set), sep='\t', header=TRUE, check.names = FALSE)
names_expressed_mrna = rownames(tbl)
rownames(tbl) = c()
names_mrna = toupper(names_expressed_mrna)
expr_mrna = as.data.table(tbl)
print(sprintf("Shape of gene expression matrix: %s", paste(dim(expr_mrna), collapse = "x")))

# load mirna expr
tbl = fread(sprintf("%s_%s/%s_%s_quantification_%s.detection_rate_%s_per_%s.csv", data_input$raw_data_folder, data_input$data_sub_set, data_input$rna_class, data_input$feature_filtering, data_input$norm, data_input$detection_rate, data_input$detection_group), sep='\t', header=T)
expr_mirna = tbl[, sapply(tbl, is.numeric), with=F]
names_mirna = tbl[[data_input$rna_class]]
print(sprintf("Shape of miRNA expression matrix: %s", paste(dim(expr_mirna), collapse = "x")))

# load cas list
gene_list_cas = fread(gene_list_cas_folder_path, sep='\t', header=T)
gene_list_cas$Gene = toupper(gene_list_cas$Gene)
print(sprintf("Length of cas gene list: %s", length(gene_list_cas$Gene)))

# load target gene list from mirtargetlink
gene_list_target_genes_all = fread(gene_list_target_genes_folder_path, sep='\t', header=T)
gene_list_target_genes_all$`Target / Target pathway` = toupper(gene_list_target_genes_all$`Target / Target pathway`)
print(sprintf("Length of unique target gene list: %s", length(unique(gene_list_target_genes_all$`Target / Target pathway`))))
gene_list_target_genes_functional = gene_list_target_genes_all[gene_list_target_genes_all$Support == "Functional MTI",]
print(sprintf("Length of unique target gene list without weak: %s", length(unique(gene_list_target_genes_functional$`Target / Target pathway`))))

# load supplementary table from paper 10.1080/15476286.2025.2449775 Hart et al.
tmp = read.xlsx(gene_list_target_genes_martin_file$path, sheet = gene_list_target_genes_martin_file$sheet, startRow = 2, colNames = TRUE, rowNames = FALSE)$Hart.et.al.
gene_list_target_genes_martin = as.vector(na.omit(tmp))
print(sprintf("Length of unique target gene from paper: %s", length(unique(gene_list_target_genes_martin))))

# load genes from mtor pathway
gene_list_target_genes_mtor = fread(gene_list_mtor, sep='\t', header=T)$mRNA
print(sprintf("Length of unique target gene from mtor: %s", length(unique(gene_list_target_genes_all))))

# create output folders
output_folder = sprintf("%s/%s_%s/results_%s/matrices/%s_target_genes_mtor", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, feature)
dir.create(output_folder, recursive=TRUE)
output_folder_fig = sprintf("%s/%s_%s/results_%s/figures/%s_target_genes_mtor", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, feature)
dir.create(output_folder_fig, recursive=TRUE)

#------------------------------------ Script ----------------------------------- 
# take eather all target genes or only non weak ones
#i = 1
#for (gene_list_target_genes in list(unique(gene_list_target_genes_all$`Target / Target pathway`), unique(gene_list_target_genes_functional$`Target / Target pathway`))) {
#  if (i == 1) {
#    target_genes_source = "target_genes_all"
#  } else {
#    target_genes_source = "target_genes_without_weak"
#  }
gene_list_target_genes = unique(gene_list_target_genes_all$`Target / Target pathway`)
target_genes_source = "target_genes_mtor"

# is there an overlap between the cas genes and the mtor genes
common_genes = intersect(gene_list_cas$Gene, gene_list_target_genes_mtor)
print(sprintf("Between CAS and the mtor genes, there is an overlap of %s genes (%s)", length(common_genes), paste(common_genes, collapse = ", ")))

# is there an overlap between the mirtargetlink results for mir-155-5p and the mtor genes
common_genes = intersect(gene_list_target_genes, gene_list_target_genes_mtor)
print(sprintf("Between all targets from mir-155-5p based on mirtargetlink and the mtor genes, there is an overlap of %s genes (%s)", length(common_genes), paste(common_genes, collapse = ", ")))

# is there an overlap between the targets from martins paper and the mtor genes
common_genes = intersect(gene_list_target_genes_martin, gene_list_target_genes_mtor)
print(sprintf("Between target genes from mir-155-5p based on martins paper and the mtor genes, there is an overlap of %s genes (%s)", length(common_genes), paste(common_genes, collapse = ", ")))

# filter samples
# get intersection of column names
common_samples = intersect(colnames(expr_mrna), colnames(expr_mirna))

# select columns from intersection
expr_mrna = expr_mrna[, ..common_samples]
expr_mirna = expr_mirna[, ..common_samples]
print(sprintf("Shape after intersection filtering of gene expression matrix: %s", paste(dim(expr_mrna), collapse = "x")))
print(sprintf("Shape after intersection filtering of miRNA expression matrix: %s", paste(dim(expr_mirna), collapse = "x")))

# filter annot also to only keep the samples which are in both expression matrices
print(sprintf("Number of samples in annot before filtering: %s", dim(annot)[1]))
tmp = (annot[[data_input$identifier_column]] %in% common_samples)
annot_intersection = annot[tmp,]
print(sprintf("Shape of annot after filtering: %s", dim(annot_intersection)[1]))

# filter features
# filter genes by target gene list
tmp = (names_mrna %in% gene_list_target_genes_mtor)
expr_mrna_genes = expr_mrna[tmp,]
names_mrna_genes = names_mrna[tmp]

# if there are some genes not included in the expression file
if (length(setdiff(gene_list_target_genes_mtor, names_mrna)) > 0) {
  print(sprintf("There are %s (%s) genes missing in the expression which a contained in the inserted target gene file", length(setdiff(gene_list_target_genes_mtor, names_mrna)), paste(setdiff(gene_list_target_genes_mtor, names_mrna), collapse = ", ")))
}
print(sprintf("Shape of gene expression matrix: %s", paste(dim(expr_mrna_genes), collapse = "x")))

# filter miRNA expression by given feature from the config
expr_mirna_feature = expr_mirna[names_mirna == feature,]

correlation_ti = c()
padj_ti = c()
for (ti in unique(annot_intersection[[split_prop]])) {
  # filter expression files for tissue ti
  group_ids = annot_intersection[annot_intersection[[split_prop]] == ti,][[data_input$identifier_column]]
  expr_mrna_genes_ti = expr_mrna_genes[,..group_ids]
  expr_mirna_feature_ti = expr_mirna_feature[,..group_ids]
  if (all(colnames(expr_mrna_genes_ti) == colnames(expr_mirna_feature_ti))) {
  } else {
    print(ti)
  }
  
  # correlation calculation   
  for (method in corr_params$method) {
    #print(sprintf("%s: %s", method, ti))
    # create folder
    dir.create(sprintf("%s/%s/%s", output_folder, target_genes_source, method), recursive=TRUE)
    
    correlation = matrix(NA, nrow(expr_mirna_feature_ti), nrow(expr_mrna_genes_ti))
    pval = matrix(NA, nrow(expr_mirna_feature_ti), nrow(expr_mrna_genes_ti))
    
    for (i in 1:nrow(expr_mirna_feature_ti)) {
      x = as.numeric(expr_mirna_feature_ti[i, ])
      correlation[i, ] = foreach (j=1:nrow(expr_mrna_genes_ti), .combine = "rbind") %do% {
        # calculate correlation
        y = as.numeric(expr_mrna_genes_ti[j, ])
        cor(x, y, method=method)
      }
      pval[i, ] = foreach (j=1:nrow(expr_mrna_genes_ti), .combine = "rbind") %do% {
        # calculate correlation
        y = as.numeric(expr_mrna_genes_ti[j, ])
        cor.test(x, y, method=method)$p.value
      }
      if ((i %% 10) == 0) {
        print(i)
      }
    }
    
    # adjustment
    padj = matrix(p.adjust(pval, corr_params$adjustment), nrow(pval), ncol(pval))
    
    # correlation
    # add colnames and mirna column
    correlation = data.frame(correlation)
    colnames(correlation) = names_mrna_genes
    correlation$rna = feature
    # reverse column order such that mirna are in front
    correlation = correlation[, rev(colnames(correlation))]
    
    # pvalue
    pval = data.frame(pval)
    colnames(pval) = names_mrna_genes
    pval$rna = feature
    # reverse column order such that mirna are in front
    pval = pval[, rev(colnames(pval))]
    
    # adjusted
    padj = data.frame(padj)
    colnames(padj) = names_mrna_genes
    padj$rna = feature
    padj = padj[, rev(colnames(padj))]

    # save corr data files
    write.table(correlation, sprintf("%s/%s/%s/correlation_%s.csv", output_folder, target_genes_source, method, ti), sep='\t', row.names = F)
    write.table(pval, sprintf("%s/%s/%s/pval_%s.csv", output_folder, target_genes_source, method, ti), sep='\t', row.names = F)
    write.table(padj, sprintf("%s/%s/%s/padj_%s.csv", output_folder, target_genes_source, method, ti), sep='\t', row.names = F)
    
    write.xlsx(correlation, sprintf("%s/%s/%s/correlation_%s.xlsx", output_folder, target_genes_source, method, ti), colNames = TRUE, rowNames = FALSE, append = FALSE)
    write.xlsx(pval, sprintf("%s/%s/%s/pval_%s.xlsx", output_folder, target_genes_source, method, ti), colNames = TRUE, rowNames = FALSE, append = FALSE)
    write.xlsx(padj, sprintf("%s/%s/%s/padj_%s.xlsx", output_folder, target_genes_source, method, ti), colNames = TRUE, rowNames = FALSE, append = FALSE)
    
    correlation_ti[[ti]] = correlation
    padj_ti[[ti]] = padj
  }
}
#  i = i + 1
#}
for (method in corr_params$method) {
  for (corr_th in c(0.5, 0.3)) {
    sig_neg_correlated_genes = c()
    info_table = c()
    for (ti in unique(annot_intersection[[split_prop]])){
      # if padj is not sorted like correlation we do so
      padj = padj_ti[[ti]][colnames(correlation_ti[[ti]])]
      
      # filter for sig. neg. correlated genes
      keep_genes = ((correlation_ti[[ti]] <= -corr_th) & (padj < corr_params$sig_lvl))
      sig_neg_correlated_genes[[ti]] = colnames(correlation_ti[[ti]][,keep_genes, drop = FALSE])
      
      for (rna in sig_neg_correlated_genes[[ti]]) {
        info_table = rbind(info_table, c(rna, ti, correlation_ti[[ti]][[rna]], padj[[rna]]))
      }
    }
    colnames(info_table) = c("mRNA", split_prop, "corr_value", "adj_p_value")
    
    write.table(info_table, sprintf("%s/%s/%s/sig_neg_correlation_corr_the=%s.csv", output_folder, target_genes_source, method, corr_th), sep='\t', row.names = F)
    write.xlsx(info_table, sprintf("%s/%s/%s/sig_neg_correlation_corr_the=%s.xlsx", output_folder, target_genes_source, method, corr_th), colNames = TRUE, rowNames = FALSE, append = FALSE)
  }
}

# calculate correlation with age from the mRNA under consideration
for (method in corr_params$method) {
  
  correlation_age = c()
  corr_pvalues_age = c()
  for(group in unique(annot_intersection[[split_prop]])){
    
    groups_ID = annot_intersection[annot_intersection[[split_prop]] == group,][[data_input$identifier_column]]
    tmp = annot_intersection[[data_input$identifier_column]] %in% groups_ID
    timepoints = annot_intersection[tmp,][[corr_prop]]
    #tmp = colnames(expr_mrna_genes) %in% groups_ID
    #expr_sort = data.frame(feature=names_mrna_genes, expr_mrna_genes[,..tmp], check.names = FALSE)
    expr_sort = data.frame(feature=names_mrna_genes, expr_mrna_genes[,..groups_ID], check.names = FALSE)
    
    if (!is.numeric(timepoints)) {
      unique_strings = unique(timepoints)
      string_to_number = setNames(seq_along(unique_strings), unique_strings)
      
      numeric_timepoints = string_to_number[timepoints]
    } else {
      numeric_timepoints = timepoints
    }
    
    if (length(timepoints) <= 4) {
      print(sprintf("For group %s, the number of samples is %s (must at least be 5)", group, length(timepoints)))
      print("skip")
    } else {
      correlation_age[[group]] = as.data.frame(cor(t(expr_sort[,2:dim(expr_sort)[2]]), numeric_timepoints, method = method))
      # if all expr for every timepoint for a sample is 0 for a feature
      # then the correlation is NA
      # we then force the correlation to be 0
      na_cor = is.na(correlation_age[[group]])
      # uncomment this to show the corresponding lines in the expr_matrix
      # print(expr_sort[,2:dim(expr_sort)[2]][na_cor, ])
      # replace NA with 0
      correlation_age[[group]][na_cor, ] = 0
      # get pvalues from rcorr
      # $P gets p-values
      pvalues_tmp = as.data.frame(rcorr(t(expr_sort[,2:dim(expr_sort)[2]]), numeric_timepoints, type=method)$P)
      # select y because x=matrix, y=timepoints
      # and remove last value because it is correlation of timepoints with timepoints = NA
      corr_pvalues_age[[group]] = head(pvalues_tmp$y, -1)
      # adjustment
      corr_pvalues_age[[group]] = as.data.frame(p.adjust(corr_pvalues_age[[group]], method=corr_params$adjustment))
      # remove all pvalues == NA and replace them with 1 (then it is not significant)
      na_pvalue = is.na(corr_pvalues_age[[group]])
      corr_pvalues_age[[group]][na_pvalue] = 1
      # check if anything is still NA
      #print(sprintf("NA in correlation_age: %s", any(is.na(correlation_age[[group]]))))
      #print(sprintf("NA in pvalues: %s", any(is.na(corr_pvalues_age[[group]]))))
    }
  }
  
  correlation_age_df = as.data.frame(correlation_age)
  colnames(correlation_age_df) = unique(annot[[split_prop]])
  correlation_age_df$RNA = expr_sort$feature
  correlation_age_df = correlation_age_df[, rev(colnames(correlation_age_df))]
  
  corr_pvalues_age_df = as.data.frame(corr_pvalues_age)
  colnames(corr_pvalues_age_df) = unique(annot[[split_prop]])
  corr_pvalues_age_df$RNA = expr_sort$feature
  corr_pvalues_age_df = corr_pvalues_age_df[, rev(colnames(corr_pvalues_age_df))]
  
  # save corr data files
  write.table(correlation_age_df, sprintf("%s/%s/%s/correlation_with_%s_%s.csv", output_folder, target_genes_source, method, corr_prop, ti), sep='\t', row.names = F)
  write.table(corr_pvalues_age_df, sprintf("%s/%s/%s/padj_with_%s_%s.csv", output_folder, target_genes_source, method, corr_prop, ti), sep='\t', row.names = F)
  
  write.xlsx(correlation_age_df, sprintf("%s/%s/%s/correlation_with_%s_%s.xlsx", output_folder, target_genes_source, method, corr_prop, ti), colNames = TRUE, rowNames = FALSE, append = FALSE)
  write.xlsx(corr_pvalues_age_df, sprintf("%s/%s/%s/padj_with_%s_%s.xlsx", output_folder, target_genes_source, method, corr_prop, ti), colNames = TRUE, rowNames = FALSE, append = FALSE)
  
  # create scatter plot per brain region
  correlation_age_df_melted = melt(correlation_age_df)
  
  for (corr_th in corr_params$corr_th) {
    correlation_age_scatter = ggplot(correlation_age_df_melted, aes(x = RNA, y = value, color = variable)) +
      geom_point(size = 3, alpha = 0.8) +  # Scatter points
      scale_color_manual(values = unlist(colors[[split_prop]][unique(correlation_age_df_melted$variable)])) +  # Apply predefined colors
      geom_hline(yintercept = c(corr_th, -corr_th), linetype = "dashed", color = "grey") +  # Add threshold lines
      xlab("RNA") +
      ylab("Value") +
      ylim (-1, 1) +
      theme_classic() +
      theme(#plot.margin = unit(c(0.1,-1.9,0.25,-1.75), plots_props$image_units), #top, right, bottom, left
        plot.margin = unit(c(0,0.1,0,0.1), plots_props$image_units), #top, right, bottom, left
        #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 9), 
        plot.background = element_rect(fill='transparent', color=NA),
        legend.position="none", 
        text = element_text(family = plots_props$font_family, size = plots_props$font_size),
        axis.text = element_text(family = plots_props$font_family, size = plots_props$font_size),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title = element_text(family = plots_props$font_family, size = plots_props$font_size),
        plot.title = element_text(family = plots_props$font_family, size = plots_props$font_size_header),
        legend.text = element_text(family = plots_props$font_family, size = plots_props$font_size),
        legend.title = element_text(family = plots_props$font_family, size = plots_props$font_size)
      )
    #correlation_age_scatter
    
    correlation_age_df_sorted = correlation_age_df[append("RNA", names(unlist(colors[[split_prop]][unique(correlation_age_df_melted$variable)])))]
    corr_pvalues_age_df_sorted = corr_pvalues_age_df[append("RNA", names(unlist(colors[[split_prop]][unique(correlation_age_df_melted$variable)])))]
    
    row_names = correlation_age_df_sorted$RNA
    correlation_age_df_sorted$RNA = c()
    
    annotation_row = names(unlist(colors[[split_prop]][unique(correlation_age_df_melted$variable)]))
  
    #heatmap_plot = plot_heatmap_annot(correlation_age_df_sorted, corr_pvalues_age_df_sorted, row_names, data_input, xticks_names, plots_props, split_prop, corr_prop, method, annotation_row, colors, corr_th, corr_params$sig_lvl)
    heatmap_plot = plot_heatmap_annot_flipped(correlation_age_df_sorted, corr_pvalues_age_df_sorted, row_names, data_input, xticks_names, plots_props, split_prop, corr_prop, method, annotation_row, colors, corr_th, corr_params$sig_lvl)
    #df = correlation_age_df_sorted
    #pvalues_df = corr_pvalues_age_df_sorted
    #prop = split_prop
    #time = corr_prop
    #m = method
    #th = corr_th
    #sig_lvl = corr_params$sig_lvl
    
    #row_name_square_plot = row_name_squares(names(unlist(colors[[split_prop]][unique(correlation_age_df_melted$variable)])), xticks_names, split_prop, colors, plots_props)
    row_name_square_plot = row_name_squares_top(names(unlist(colors[[split_prop]][unique(correlation_age_df_melted$variable)])), xticks_names, split_prop, colors, plots_props)
    
    heatmap_plot_comb = row_name_square_plot %v% heatmap_plot$heatmap
    
    # save figures
    height = 2*plots_props$image_height
    width = 15
    filename = sprintf("%s/corr_method=%s_%s_over_%s_per_%s", output_folder_fig, method, data_input$rna_class, corr_prop, split_prop)
    png(sprintf("%s_corr_th=%s_%sx%s.png", filename, corr_th, width, height), width = width, height = height, units = plots_props$image_units, res = plots_props$dpi)
    draw(heatmap_plot_comb, gap = unit(0, plots_props$image_units), padding = unit(c(-1, 0.25, -1, 1.75), "cm"))  # bottom, left, top and right margins
    draw(heatmap_plot$legend, x = unit(1.5*plots_props$image_width - 0.1, "cm"), y = unit(0.4, "cm"), just = c("left", "bottom"))
    dev.off()
    svglite(sprintf("%s_corr_th=%s_%sx%s.svg", filename, corr_th, width, height), width = width / 2.54, height = height / 2.54)
    draw(heatmap_plot_comb, gap = unit(0, plots_props$image_units), padding = unit(c(-1, 0.25, -1, 1.75), "cm"))  # bottom, left, top and right margins
    draw(heatmap_plot$legend, x = unit(1.5*plots_props$image_width - 0.1, "cm"), y = unit(0.4, "cm"), just = c("left", "bottom"))
    dev.off()
    
    for(gr in unique(annot_intersection[[split_prop]])){
      correlation_ti[[gr]]$group = gr
      padj_ti[[gr]]$group = gr
    }
    correlation_ti_df = do.call(rbind, correlation_ti)
    correlation_ti_df$rna = c()
    padj_ti_df = do.call(rbind, padj_ti)
    padj_ti_df$rna = c()
    
    correlation_ti_df_sorted = t(correlation_ti_df[order(match(correlation_ti_df$group, names(colors[[split_prop]]))),])
    padj_ti_df_sorted = t(padj_ti_df[order(match(padj_ti_df$group, names(colors[[split_prop]]))),])
    
    # remove last row
    correlation_ti_df_sorted = data.frame(correlation_ti_df_sorted[1:nrow(correlation_ti_df_sorted) - 1, ])
    padj_ti_df_sorted = data.frame(padj_ti_df_sorted[1:nrow(padj_ti_df_sorted) - 1, ])
    
    row_names = rownames(correlation_ti_df_sorted)
    rownames(correlation_ti_df_sorted) = c()
    rownames(padj_ti_df_sorted) = c()
    
    # force numeric
    correlation_ti_df_sorted = data.frame(lapply(correlation_ti_df_sorted, as.numeric))
    padj_ti_df_sorted = data.frame(lapply(padj_ti_df_sorted, as.numeric))
    
    annotation_row = names(unlist(colors[[split_prop]][colnames(correlation_ti_df_sorted)]))
    
    #heatmap_plot = plot_heatmap_annot(correlation_age_df_sorted, corr_pvalues_age_df_sorted, row_names, data_input, xticks_names, plots_props, split_prop, corr_prop, method, annotation_row, colors, corr_th, corr_params$sig_lvl)
    heatmap_plot = plot_heatmap_annot_flipped(correlation_ti_df_sorted, padj_ti_df_sorted, row_names, data_input, xticks_names, plots_props, split_prop, "", method, annotation_row, colors, corr_th, corr_params$sig_lvl)
    #df = correlation_ti_df_sorted
    #pvalues_df = padj_ti_df_sorted
    #prop = split_prop
    #time = corr_prop
    #m = method
    #th = corr_th
    #sig_lvl = corr_params$sig_lvl
    
    #row_name_square_plot = row_name_squares(names(unlist(colors[[split_prop]][unique(correlation_age_df_melted$variable)])), xticks_names, split_prop, colors, plots_props)
    row_name_square_plot = row_name_squares_top(names(unlist(colors[[split_prop]][unique(correlation_age_df_melted$variable)])), xticks_names, split_prop, colors, plots_props)
    
    heatmap_plot_comb = row_name_square_plot %v% heatmap_plot$heatmap
    
    # save figures
    height = 2*plots_props$image_height
    width = 15
    filename = sprintf("%s/corr_method=%s_%s_per_%s", output_folder_fig, method, feature, split_prop)
    png(sprintf("%s_corr_th=%s_%sx%s.png", filename, corr_th, width, height), width = width, height = height, units = plots_props$image_units, res = plots_props$dpi)
    draw(heatmap_plot_comb, gap = unit(0.2, plots_props$image_units), padding = unit(c(-1, 0.25, -1, 1.75), "cm"))  # bottom, left, top and right margins
    draw(heatmap_plot$legend, x = unit(1.5*plots_props$image_width - 0.1, "cm"), y = unit(0.3, "cm"), just = c("left", "bottom"))
    dev.off()
    svglite(sprintf("%s_corr_th=%s_%sx%s.svg", filename, corr_th, width, height), width = width / 2.54, height = height / 2.54)
    draw(heatmap_plot_comb, gap = unit(0, plots_props$image_units), padding = unit(c(-1, 0.25, -1, 1.75), "cm"))  # bottom, left, top and right margins
    draw(heatmap_plot$legend, x = unit(1.5*plots_props$image_width - 0.1, "cm"), y = unit(0.3, "cm"), just = c("left", "bottom"))
    dev.off()
    
    #heatmap_plot = plot_heatmap_annot(correlation_age_df_sorted, corr_pvalues_age_df_sorted, row_names, data_input, xticks_names, plots_props, split_prop, corr_prop, method, annotation_row, colors, corr_th, corr_params$sig_lvl)
    heatmap_plot = plot_heatmap_annot_flipped_tissue_flipped(correlation_ti_df_sorted, padj_ti_df_sorted, row_names, data_input, xticks_names, plots_props, split_prop, "", method, annotation_row, colors, corr_th, corr_params$sig_lvl)
    #df = correlation_ti_df_sorted
    #pvalues_df = padj_ti_df_sorted
    #prop = split_prop
    #time = corr_prop
    #m = method
    #th = corr_th
    #sig_lvl = corr_params$sig_lvl
    
    #row_name_square_plot = row_name_squares(names(unlist(colors[[split_prop]][unique(correlation_age_df_melted$variable)])), xticks_names, split_prop, colors, plots_props)
    row_name_square_plot = row_name_squares_top_flipped(names(unlist(colors[[split_prop]][unique(correlation_age_df_melted$variable)])), xticks_names, split_prop, colors, plots_props)
    
    heatmap_plot_comb = row_name_square_plot %v% heatmap_plot$heatmap
    
    # save figures
    height = 14
    width = 10
    filename = sprintf("%s/corr_method=%s_%s_per_%s", output_folder_fig, method, feature, split_prop)
    png(sprintf("%s_corr_th=%s_tissue_flipped_%sx%s.png", filename, corr_th, width, height), width = width, height = height, units = plots_props$image_units, res = plots_props$dpi)
    draw(heatmap_plot_comb, gap = unit(0.2, plots_props$image_units), padding = unit(c(-1, 0.25, -1, 1.75), "cm"))  # bottom, left, top and right margins
    draw(heatmap_plot$legend, x = unit(1*width - 1.6, "cm"), y = unit(0.7, "cm"), just = c("left", "bottom"))
    dev.off()
    svglite(sprintf("%s_corr_th=%s_tissue_flipped_%sx%s.svg", filename, corr_th, width, height), width = width / 2.54, height = height / 2.54)
    draw(heatmap_plot_comb, gap = unit(0, plots_props$image_units), padding = unit(c(-1, 0.25, -1, 1.75), "cm"))  # bottom, left, top and right margins
    draw(heatmap_plot$legend, x = unit(1*width - 1.6, "cm"), y = unit(0.7, "cm"), just = c("left", "bottom"))
    dev.off()
  }
}


#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])
