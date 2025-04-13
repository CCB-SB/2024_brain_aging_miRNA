suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(viridisLite))
suppressPackageStartupMessages(library(svglite))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(stringi))

#snakemake = readRDS("snakemake_hclust_expression_brain_region_for_mirna_list_zscores_heatmap_plot.rds")

data.table::setDTthreads(1)
set.seed(snakemake@params$parameters_props$set_seed)

print("hclust_expression_brain_region_for_mirna_list_zscore_heatmap_plot")


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

plot_heatmap_row_feature_portrait = function(df, row_names, prop, params, legend_title, cluster_colours_list, plot_props, width, height, filename_fig, filename_tab){
  df = as.matrix(df)
  rownames(df) = row_names
  
  annotation_col = data.frame(dummy=1:ncol(df))
  annotation_col[[prop]] = c(as.character(annot[[prop]][match(colnames(df), annot[[prop]])]))
  rownames(annotation_col) = colnames(df)
  annotation_col$dummy = NULL
  
  annotation_colors = list()
  ucolors <- unlist(colors[[prop]])
  if(is.vector(ucolors) && all(annotation_col[[prop]] %in% names(ucolors))) {
    annotation_colors[[prop]] = ucolors
  }
  
  #show_rownames = nrow(df) < 50
  show_rownames = TRUE
  # print(annotation_col)
  column_ha = HeatmapAnnotation(df = annotation_col, col = annotation_colors, gp = gpar(col = "white"), simple_anno_size = unit(0.15, "cm"), show_annotation_name = FALSE, 
                                annotation_name_side = "left", annotation_name_gp = gpar(fontsize = plot_props$font_size, fontface = "bold", fontfamily=plot_props$font_family))
  
  col_hc = hclust(dist(t(df)), method = params$method)
  # number of cluster is params$number_of_cluster
  # if heatmap has only one unique values, use k=1 (else there is an error)
  all_values = as.numeric(df)
  if (length(unique(all_values)) == 1) {
    k = 1
  } else {
    k = params$number_of_cluster
  }
  col_clusters = cutree(col_hc, k = k)
  # make roman numerals
  # also make them strings so that we can set names
  col_clusters_roman = as.character(as.roman(col_clusters))
  names(col_clusters_roman) = names(col_clusters)
  # make a dataframe to save
  col_clusters_df = data.table(cluster=col_clusters, cluster_roman=col_clusters_roman)
  col_clusters_df[[prop]] = names(col_clusters)

  # save cluster info values for all tissues 
  fwrite(col_clusters_df, sprintf("%s.csv", filename_tab), sep = "\t", row.names = FALSE)
  write.xlsx(col_clusters_df, sprintf("%s.xlsx", filename_tab), colNames = TRUE, rowNames = FALSE, append = FALSE)

  if (length(unique(as.numeric(df))) <= 2) {
    #legend_ticks = seq(from = 0, to = max(t(df)), length.out = 5)
    maxi = max(df)
    #col_fun = colorRamp2(c(0, maxi), c(viridis(100)[1],viridis(100)[100]))
    if (maxi == 0) {
      maxi = 1 
    }
    col_fun = colorRamp2(c(0, maxi), c("#FFF6D6","#D5A021"))
    #col_fun = colorRamp2(c(0,maxi), hcl_palette = "Greens", reverse = TRUE)
    col_fun_legend = col_fun
  } else {
    maxi = max(df)
    #col_fun = colorRamp2(c(0,maxi), hcl_palette = "Greens", reverse = TRUE)
    col_fun = colorRamp2(c(0,maxi), c("#FFF6D6","#D5A021"), reverse = TRUE)
    col_fun_legend = col_fun
  }

  # this makes the cluster names more centered in color block
  ht_opt$TITLE_PADDING = unit(c(4, 4), "points")
  # features rows
  p = ComplexHeatmap::Heatmap(df, col = col_fun,
                              rect_gp = gpar(col = "white", lwd = .05),
                              #height = unit(0.005, "cm") * nrow(df), width = unit(0.3, "cm") * ncol(df),
                              height = unit(height + 2, "cm"), width = unit(width - 3, "cm"),
                              row_names_side = "left", row_names_gp = gpar(fontsize = plot_props$font_size, fontfamily=plot_props$font_family),
                              column_names_side = "bottom", column_names_gp = gpar(fontsize = plot_props$font_size, fontfamily=plot_props$font_family),
                              cluster_columns = function(m) hclust(dist(m), method = params$method), clustering_distance_columns = params$dist_method, show_column_dend = FALSE,
                              cluster_rows = function(m) hclust(dist(m), method = params$method), clustering_distance_rows = params$dist_method, show_row_dend = FALSE,
                              column_split = col_clusters_roman, cluster_column_slices = FALSE,  # last argument keeps the clusters ordered from I to IV
                              column_title_gp = gpar(col="white", fill = cluster_colours_list),
                              column_dend_reorder = FALSE, #cluster_columns = FALSE,
                              row_dend_reorder = FALSE, #cluster_rows = TRUE,
                              show_column_names = FALSE,
                              show_row_names = show_rownames,
                              show_heatmap_legend = FALSE,
                              top_annotation = column_ha
  )
  
  if (length(unique(as.numeric(df))) <= 2) {
    # discrete legend
    # legend_ticks = c(0, maxi)
    legend_tick_locations = c(0, 1)
    legend_tick_labels = c("0" = sprintf("|zscore| < %s" , zscore_th), "1" = sprintf("|zscore| ≥ %s", zscore_th))
    legend = ComplexHeatmap::Legend(title = legend_title, at = legend_tick_locations, labels=legend_tick_labels, direction="horizontal",
                                    legend_gp = gpar(fill = col_fun_legend(legend_tick_locations)),
                                    title_gp = gpar(fontsize = plot_props$font_size, fontface = "bold", fontfamily=plot_props$font_family), labels_gp = gpar(fontsize = plot_props$font_size, fontfamily=plot_props$font_family))
    
  } else {
    # continuous (colorbar) legend
    legend_ticks = round(seq(from = min(t(df)), to = maxi, length.out = 5))
    legend_tick_locations = round(legend_ticks, digits = 2)
    legend = ComplexHeatmap::Legend(title = legend_title, col_fun = col_fun_legend, at = legend_tick_locations, direction="horizontal",
                                    title_gp = gpar(fontsize = plot_props$font_size, fontface = "bold", fontfamily=plot_props$font_family), labels_gp = gpar(fontsize = plot_props$font_size, fontfamily=plot_props$font_family))
    
  }
  
  png(sprintf("%s_row=%s_portrait.png", filename_fig, data_input$rna_class), width = width, height = 2 * height, units = plot_props$image_units, res = plot_props$dpi)
  draw(p, show_annotation_legend = FALSE, padding = unit(c(1.5, 0.25, 0, 0.25), "cm")) 
  draw(legend, x = unit(1.3, "cm"), y = unit(0.02, "cm"), just = c("left", "bottom"))
  dev.off()
  svglite(sprintf("%s_row=%s_portrait.svg", filename_fig, data_input$rna_class), width = width / 2.54, height = 2 * height / 2.54)
  draw(p, show_annotation_legend = FALSE, padding = unit(c(1.5, 0.25, 0, 0.25), "cm"))
  draw(legend, x = unit(1.3, "cm"), y = unit(0.02, "cm"), just = c("left", "bottom"))
  dev.off()
}

plot_heatmap_row_feature_portrait_mod = function(df, df_bin, row_names, prop, params, legend_title, cluster_colours_list, plot_props, rownames_pos, show_legend, legend_position, width, height, filename_fig, filename_tab){
  df = as.matrix(df)
  rownames(df) = row_names
  
  df_bin = as.matrix(df_bin)
  rownames(df_bin) = row_names
  
  annotation_col = data.frame(dummy=1:ncol(df))
  annotation_col[[prop]] = c(as.character(annot[[prop]][match(colnames(df), annot[[prop]])]))
  rownames(annotation_col) = colnames(df)
  annotation_col$dummy = NULL
  
  annotation_colors = list()
  ucolors <- unlist(colors[[prop]])
  if(is.vector(ucolors) && all(annotation_col[[prop]] %in% names(ucolors))) {
    annotation_colors[[prop]] = ucolors
  }
  
  #show_rownames = nrow(df) < 50
  show_rownames = TRUE
  # print(annotation_col)
  column_ha = HeatmapAnnotation(df = annotation_col, col = annotation_colors, gp = gpar(col = "white"), simple_anno_size = unit(0.15, "cm"), show_annotation_name = FALSE, 
                                annotation_name_side = "left", annotation_name_gp = gpar(fontsize = plot_props$font_size, fontface = "bold", fontfamily=plot_props$font_family))
  
  col_hc = hclust(dist(t(df_bin)), method = params$method)
  col_clusters = cutree(col_hc, k = params$number_of_cluster)
  # make roman numerals
  # also make them strings so that we can set names
  col_clusters_roman = as.character(as.roman(col_clusters))
  names(col_clusters_roman) = names(col_clusters)
  # make a dataframe to save
  col_clusters_df = data.table(cluster=col_clusters, cluster_roman=col_clusters_roman)
  col_clusters_df[[prop]] = names(col_clusters)
  
  # save cluster info values for all tissues 
  fwrite(col_clusters_df, sprintf("%s.csv", filename_tab), sep = "\t", row.names = FALSE)
  write.xlsx(col_clusters_df, sprintf("%s.xlsx", filename_tab), colNames = TRUE, rowNames = FALSE, append = FALSE)
  
  if (length(unique(as.numeric(df))) <= 2) {
    #legend_ticks = seq(from = 0, to = max(t(df)), length.out = 5)
    maxi = max(df)
    #col_fun = colorRamp2(c(0, maxi), c(viridis(100)[1],viridis(100)[100]))
    col_fun = colorRamp2(c(0, maxi), c("#FFF6D6","#D5A021"))
    #col_fun = colorRamp2(c(0,maxi), hcl_palette = "Greens", reverse = TRUE)
    col_fun_legend = col_fun
  } else {
    maxi = max(df)
    #col_fun = colorRamp2(c(0,maxi), hcl_palette = "Greens", reverse = TRUE)
    col_fun = colorRamp2(c(0,maxi), c("#FFF6D6","#D5A021"), reverse = TRUE)
    col_fun_legend = col_fun
  }
  
  # this makes the cluster names more centered in color block
  ht_opt$TITLE_PADDING = unit(c(4, 4), "points")
  
  # cell_fun styles the cells where df_bin == 1 as a circle
  cell_fun = function(i, j, x, y, width, height, fill) {
    # linewidth outer border (white)
    lwd_outer = unit(1, "mm")
    # marker linewidth if df_bin == 1
    lwd_inner = unit(1, "mm")
    # color for marker
    col_marker = "black"
    # fill color
    color = col_fun(df[j, i])
    # rectangle, fill with fill color, line in white, lwd to separate cells
    grid.rect(x = x, y = y, width = width, height = height, gp = gpar(fill=color, col = "white", lwd = lwd_outer))
    if (df_bin[j, i] == 1) {
      # draw circle over rectangle
      # grid.circle(x = x, y = y, r = unit(.5, "mm"), gp = gpar(fill = "red", col = NA))
      # draw outer rectangle as border
      grid.rect(x = x, y = y, width = (width - lwd_outer / 2), height = (height - lwd_outer / 2), gp = gpar(fill=color, col = col_marker, lwd = lwd_inner))
    }
  }
  
  # save plot_df
  #print(df_bin)
  #print(sprintf("%s_plot_df.csv", filename_tab))
  write.table(df_bin, sprintf("%s_plot_df.csv", filename_tab), sep = "\t", row.names = TRUE, col.names = TRUE)
  write.xlsx(df_bin, sprintf("%s_plot_df.xlsx", filename_tab), colNames = TRUE, rowNames = TRUE, append = FALSE)
  
  if (show_legend == FALSE) {
    height_heatmap = height + 2
    width_heatmap = width
  } else {
    if (legend_position == "right") {
      height_heatmap = height - 1 
      width_heatmap = width 
    } else if (legend_position == "bottom")
    height_heatmap = height
    width_heatmap = width
  }
  
  # remove mmu-miR
  rownames(df_bin) = gsub("mmu-miR-", "", rownames(df_bin))

  # features rows
  p = ComplexHeatmap::Heatmap(df_bin, col = col_fun,
                              rect_gp = gpar(type = "none"), 
                              #height = unit(0.005, "cm") * nrow(df), width = unit(0.3, "cm") * ncol(df),
                              height = unit(height_heatmap + 2.5, "cm"), width = unit(width_heatmap - 3, "cm"),
                              row_names_side = rownames_pos, row_names_gp = gpar(fontsize = 6, fontfamily=plot_props$font_family),
                              column_names_side = "bottom", column_names_gp = gpar(fontsize = plot_props$font_size, fontfamily=plot_props$font_family),
                              cluster_columns = function(m) hclust(dist(m), method = params$method), clustering_distance_columns = params$dist_method, show_column_dend = FALSE,
                              cluster_rows = function(m) hclust(dist(m), method = params$method), clustering_distance_rows = params$dist_method, show_row_dend = FALSE,
                              column_split = col_clusters_roman, cluster_column_slices = FALSE,  # last argument keeps the clusters ordered from I to IV
                              column_title_gp = gpar(col="white", fill = cluster_colours_list),
                              column_dend_reorder = FALSE, #cluster_columns = FALSE,
                              row_dend_reorder = FALSE, #cluster_rows = TRUE,
                              show_column_names = FALSE,
                              show_row_names = show_rownames,
                              show_heatmap_legend = FALSE,
                              top_annotation = column_ha,
                              cell_fun=cell_fun
  )
  
  if (length(unique(as.numeric(df))) <= 2) {
    # discrete legend
    legend_ticks = c(0, maxi)
    legend_tick_locations = c(0, 1)
    legend_tick_labels = c("0" = sprintf("|zscore| < %s", zscore_th), "1" = sprintf("|zscore| ≥ %s", zscore_th))
    if (legend_position == "bottom") {
      legend = ComplexHeatmap::Legend(title = legend_title, at = legend_tick_locations, labels=legend_tick_labels, direction="horizontal",
                                      legend_gp = gpar(fill = col_fun_legend(legend_tick_locations)),
                                      title_gp = gpar(fontsize = plot_props$font_size, fontface = "bold", fontfamily=plot_props$font_family), labels_gp = gpar(fontsize = plot_props$font_size, fontfamily=plot_props$font_family))
    } else if (legend_position == "right") {
      legend = ComplexHeatmap::Legend(title = legend_title, at = legend_tick_locations, labels=legend_tick_labels, #direction="horizontal",
                                      legend_gp = gpar(fill = col_fun_legend(legend_tick_locations)),
                                      title_gp = gpar(fontsize = plot_props$font_size, fontface = "bold", fontfamily=plot_props$font_family), labels_gp = gpar(fontsize = plot_props$font_size, fontfamily=plot_props$font_family))
    }
  } else {
    # continuous (colorbar) legend
    if (legend_position == "bottom") {
      legend_ticks = round(seq(from = min(t(df)), to = max(t(df))))#, length.out = 5))
      legend_tick_locations = round(legend_ticks, digits = 2)
      legend = ComplexHeatmap::Legend(title = legend_title, col_fun = col_fun_legend, at = legend_tick_locations, direction="horizontal",
                                      title_gp = gpar(fontsize = plot_props$font_size, fontface = "bold", fontfamily=plot_props$font_family), labels_gp = gpar(fontsize = plot_props$font_size, fontfamily=plot_props$font_family))
    } else if (legend_position == "right") {
      legend_ticks = round(seq(from = min(t(df)), to = max(t(df))))#, length.out = 5))
      legend_tick_locations = round(legend_ticks, digits = 2)
      legend = ComplexHeatmap::Legend(title = legend_title, col_fun = col_fun_legend, at = legend_tick_locations, #direction="horizontal",
                                      title_gp = gpar(fontsize = plot_props$font_size, fontface = "bold", fontfamily=plot_props$font_family), labels_gp = gpar(fontsize = plot_props$font_size, fontfamily=plot_props$font_family))
    }
  }
  
  if (show_legend == TRUE) {
    if (legend_position == "bottom") {
      png(sprintf("%s_row=%s_portrait_mod.png", filename_fig, data_input$rna_class), width = width, height = 2 * height, units = plot_props$image_units, res = plot_props$dpi)
      draw(p, show_annotation_legend = FALSE, padding = unit(c(1.5, 0.25, 0, 0.25), "cm")) #top, right, bottom, left
      draw(legend, x = unit(1.3, "cm"), y = unit(0.02, "cm"), just = c("left", "bottom"))
    } else if (legend_position == "right") {
      png(sprintf("%s_row=%s_portrait_mod.png", filename_fig, data_input$rna_class), width = width, height = 3/2 * height, units = plot_props$image_units, res = plot_props$dpi)
      draw(p, show_annotation_legend = FALSE, padding = unit(c(0.25, 0.25, 0, 2), "cm")) #top, right, bottom, left
      draw(legend, x = unit(width, "cm"), y = unit(height / 2, "cm"), just = c("right", "center"))
    }
  } else {
    png(sprintf("%s_row=%s_portrait_mod.png", filename_fig, data_input$rna_class), width = width, height = 3/2 * height, units = plot_props$image_units, res = plot_props$dpi)
    draw(p, show_annotation_legend = FALSE, padding = unit(c(0.25, 0.25, 0, 0.25), "cm"))
  }
  dev.off()
  if (show_legend == TRUE) {
    if (legend_position == "bottom") {
      svglite(sprintf("%s_row=%s_portrait_mod.svg", filename_fig, data_input$rna_class), width = width / 2.54, height = 2 * height / 2.54)
      draw(p, show_annotation_legend = FALSE, padding = unit(c(1.5, 0.25, 0, 0.25), "cm"))
      draw(legend, x = unit(1.3, "cm"), y = unit(0.02, "cm"), just = c("left", "bottom"))
    } else if (legend_position == "right") {
      svglite(sprintf("%s_row=%s_portrait_mod.svg", filename_fig, data_input$rna_class), width = width / 2.54, height = 3/2 * height / 2.54)
      draw(p, show_annotation_legend = FALSE, padding = unit(c(0.25, 0.25, 0, 2), "cm"))
      draw(legend, x = unit(width, "cm"), y = unit(height / 2, "cm"), just = c("right", "center"))
    }
  } else {
    svglite(sprintf("%s_row=%s_portrait_mod.svg", filename_fig, data_input$rna_class), width = width / 2.54, height = 3/2 * height / 2.54)
    draw(p, show_annotation_legend = FALSE, padding = unit(c(0.25, 0.25, 0, 0.25), "cm"))
  }
  
  dev.off()
}


#---------------------------------- Load data ----------------------------------
data_input = snakemake@params$data_input
method = snakemake@params$method
top_list = snakemake@params$top_list
diff_log = snakemake@params$diff_exp_log
zscore_th_list = snakemake@params$zscore_th
prop = snakemake@params$prop
params = snakemake@params$params
colors = snakemake@params$colors
cluster_colours = snakemake@params$colors
plot_props = snakemake@params$plot_props
results_folder = snakemake@params$results_folder

# save rdata
#saveRDS(snakemake, file = "snakemake_hclust_expression_brain_region_for_mirna_list_zscores_heatmap_plot.rds")

# cluster colours
cluster_colours_list = unlist(cluster_colours$clusters)
#tmp = names(cluster_colours_list) %in% names(features_cluster)
#cluster_colours_list = cluster_colours_list[tmp]

expr = fread(sprintf("%s_%s/%s_%s_quantification_%s.detection_rate_%s_per_%s.csv", data_input$raw_data_folder, data_input$data_sub_set, data_input$rna_class, data_input$feature_filtering, data_input$norm, data_input$detection_rate, data_input$detection_group), sep='\t', header=T)
annot = fread(sprintf("%s/annotation_%s_%s.csv", data_input$annotation_data_folder, data_input$data_sub_set, data_input$feature_filtering), sep='\t', colClasses=c(ID="character")) 
diff_exp = fread(sprintf("%s/%s_%s/results_%s/matrices/diff_exp/diff_exp%s.csv", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, diff_log), sep='\t') 


#------------------------------------ Script ----------------------------------- 
method_split = data.frame(stri_split_fixed(method, "_"))[1,]
for (m in method_split) {
  print(m)
  output_folder_fig = sprintf("%s/%s_%s/results_%s/figures/%s_clustering/heatmap_complex", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, m)
  dir.create(output_folder_fig, recursive=TRUE)
  output_folder_tab = sprintf("%s/%s_%s/results_%s/matrices/%s_clustering", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, m)
  dir.create(output_folder_tab, recursive=TRUE)
  
  for(i in top_list){
    print(i)
    
    feature_list = fread(sprintf("%s/%s_%s/results_%s/matrices/%s/%s.csv", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, m, i))
    feature_list = feature_list[[data_input$rna_class]]
    
    # filter expr with feature_list
    diff_exp_filtered = diff_exp[diff_exp[[data_input$rna_class]] %in% feature_list,]
    
    median_table = data.table(feature=diff_exp_filtered[[data_input$rna_class]])
    for (ti in unique(annot[[prop]])) {
      #print(colnames(diff_exp_filtered))
      #print(sprintf("median_%s__%s", prop, ti))
      median_table[[ti]] = 10^(as.numeric(diff_exp_filtered[[sprintf("median_%s__%s", prop, ti)]]))
    }
    
    # save rnas and create row names for plotting
    rna_ids = median_table[[colnames(median_table)[1]]]
    #if(i <= 50) {
    row_names = rna_ids
    #} else {
    #  row_names = NULL
    #}
    
    
    # keep only expression
    median_table_adj = median_table[, sapply(median_table, is.numeric), with=F]
    
    for (zscore_th in zscore_th_list) {
      print(zscore_th)
      
      # calculate zscores
      params$scale = "row"
      if(params$scale != "none"){
        scale_row = params$scale %in% c("row", "both")
        scale_col = params$scale %in% c("column", "both")
        median_table_adj_zscores = scale(median_table_adj, scale_row, scale_col) #base::scale(t(median_table_adj[i,]))
        median_table_adj_zscores[is.na(median_table_adj_zscores)] = 0
        median_table_adj_zscores = abs(median_table_adj_zscores)
      }
      
      tmp = median_table_adj_zscores
      rownames(tmp) = row_names
      write.table(tmp, sprintf("%s/%s_zscore_plot_df.csv", output_folder_tab, i), sep = "\t", row.names = TRUE, col.names = TRUE)
      write.xlsx(tmp, sprintf("%s/%s_zscore_plot_df.xlsx", output_folder_tab, i), colNames = TRUE, rowNames = TRUE, append = FALSE)
      
      #print(median_table_adj)
      median_table_adj_zscores_th = ifelse(abs(median_table_adj_zscores) >= zscore_th, 1, 0)

      #print(median_table_adj_zscores_th)
      legend_title = sprintf("Bin. stand. expr. value\n(median per %s)", prop)
      plot_heatmap_row_feature_portrait(median_table_adj_zscores_th, row_names, prop, params, legend_title, cluster_colours_list, plot_props, width = plot_props$image_width, height = plot_props$image_height, filename_fig = sprintf("%s/%s_zscore_th=%s_num_clusters=%s", output_folder_fig, i, zscore_th, params$number_of_cluster), filename_tab = sprintf("%s/%s_zscore_th=%s_clustered_by_%s", output_folder_tab, i, zscore_th, m))
      
      rownames(median_table_adj_zscores_th) = row_names
      # save plot_df
      median_table_adj_zscores_th_tmp = as.data.frame(median_table_adj_zscores_th)
      median_table_adj_zscores_th_tmp[[data_input$rna_class]] = rownames(median_table_adj_zscores_th)
      rownames(median_table_adj_zscores_th_tmp) = c()
      # re-order cols
      median_table_adj_zscores_th_tmp = median_table_adj_zscores_th_tmp[, c(data_input$rna_class, colnames(median_table_adj_zscores_th))]
      write.table(median_table_adj_zscores_th_tmp, sprintf("%s/%s_zscore_th=%s_plot_df.csv", output_folder_tab, i, zscore_th), sep = "\t", row.names = FALSE, col.names = TRUE)
      write.xlsx(median_table_adj_zscores_th_tmp, sprintf("%s/%s_zscore_th=%s_plot_df.xlsx", output_folder_tab, i, zscore_th), colNames = TRUE, rowNames = FALSE, append = FALSE)
      
      # remove rows containing only 0s
      mask = rowSums(median_table_adj_zscores_th) > 0
      median_table_adj_th_removed_zeros = median_table_adj_zscores_th[mask, , drop=FALSE]
      row_names_removed_zeros = row_names[mask]
      plot_heatmap_row_feature_portrait(median_table_adj_th_removed_zeros, row_names_removed_zeros, prop, params, legend_title, cluster_colours_list, plot_props, width = plot_props$image_width, height = plot_props$image_height, filename_fig = sprintf("%s/%s_zscore_th=%s_num_clusters=%s_removed_zero_rows", output_folder_fig, i, zscore_th, params$number_of_cluster), filename_tab = sprintf("%s/%s_zscore_th=%s_clustered_by_%s_removed_zero_rows", output_folder_tab, i, zscore_th, m))
      
      # removed rows containing only zeros after binarisation
      # red border around expression values higher than the threhold
      median_table_adj_removed_zeros = median_table_adj_zscores[mask, , drop=FALSE]
      legend_title = sprintf("Abs. stand.\nexpr.\n(%s,\n median)", gsub("_norm", "", data_input$norm))
      if (nrow(median_table_adj_removed_zeros) > 0) {
        show_legend = TRUE
        rownames_pos = "left"
        legend_position = "right"
        plot_heatmap_row_feature_portrait_mod(median_table_adj_removed_zeros, median_table_adj_th_removed_zeros, row_names_removed_zeros, prop, params, legend_title, cluster_colours_list, plot_props, rownames_pos, show_legend, legend_position, width = plot_props$image_width, height = plot_props$image_height, filename_fig = sprintf("%s/%s_zscore_th=%s_num_clusters=%s_removed_zero_rows_rownames_pos=%s_show_legend=%s", output_folder_fig, i, zscore_th, params$number_of_cluster, rownames_pos, show_legend), filename_tab = sprintf("%s/%s_zscore_th=%s_clustered_by_%s_removed_zero_rows", output_folder_tab, i, zscore_th, m))
        rownames_pos = "right"
        legend_position = "bottom"
        plot_heatmap_row_feature_portrait_mod(median_table_adj_removed_zeros, median_table_adj_th_removed_zeros, row_names_removed_zeros, prop, params, legend_title, cluster_colours_list, plot_props, rownames_pos, show_legend, legend_position, width = plot_props$image_width, height = plot_props$image_height, filename_fig = sprintf("%s/%s_zscore_th=%s_num_clusters=%s_removed_zero_rows_rownames_pos=%s_show_legend=%s", output_folder_fig, i, zscore_th, params$number_of_cluster, rownames_pos, show_legend), filename_tab = sprintf("%s/%s_zscore_th=%s_clustered_by_%s_removed_zero_rows", output_folder_tab, i, zscore_th, m))
        
        show_legend = FALSE
        rownames_pos = "left"
        legend_position = "right"
        plot_heatmap_row_feature_portrait_mod(median_table_adj_removed_zeros, median_table_adj_th_removed_zeros, row_names_removed_zeros, prop, params, legend_title, cluster_colours_list, plot_props, rownames_pos, show_legend, legend_position, width = plot_props$image_width, height = plot_props$image_height, filename_fig = sprintf("%s/%s_zscore_th=%s_num_clusters=%s_removed_zero_rows_rownames_pos=%s_show_legend=%s", output_folder_fig, i, zscore_th, params$number_of_cluster, rownames_pos, show_legend), filename_tab = sprintf("%s/%s_zscore_th=%s_clustered_by_%s_removed_zero_rows", output_folder_tab, i, zscore_th, m))
        rownames_pos = "right"
        legend_position = "bottom"
        plot_heatmap_row_feature_portrait_mod(median_table_adj_removed_zeros, median_table_adj_th_removed_zeros, row_names_removed_zeros, prop, params, legend_title, cluster_colours_list, plot_props, rownames_pos, show_legend, legend_position, width = plot_props$image_width, height = plot_props$image_height, filename_fig = sprintf("%s/%s_zscore_th=%s_num_clusters=%s_removed_zero_rows_rownames_pos=%s_show_legend=%s", output_folder_fig, i, zscore_th, params$number_of_cluster, rownames_pos, show_legend), filename_tab = sprintf("%s/%s_zscore_th=%s_clustered_by_%s_removed_zero_rows", output_folder_tab, i, zscore_th, m))
      }
    }
  }
}


#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])

