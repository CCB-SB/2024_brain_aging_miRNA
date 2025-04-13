suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(viridisLite))
suppressPackageStartupMessages(library(svglite))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(stringi))
suppressPackageStartupMessages(library(cowplot))


#snakemake = readRDS("snakemake_hclust_expression_brain_region_for_mirna_list_zscores_heatmap_plot_male_female.rds")

data.table::setDTthreads(1)
set.seed(snakemake@params$parameters_props$set_seed)

print("hclust_expression_brain_region_for_mirna_list_zscores_heatmap_plot_male_female")


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

plot_heatmap_row_feature_portrait = function(df, row_names, annot, prop, params, legend_title, cluster_colours_list, plot_props, width, height, filename_fig, filename_tab){
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

plot_heatmap_row_feature_portrait_mod = function(df, df_bin, annot, row_names, max_value, prop, params, legend_title, cluster_colours_list, plot_props, rownames_pos, show_legend, legend_position, width, height, filename_fig, filename_tab){
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
  column_ha = HeatmapAnnotation(df = annotation_col, col = annotation_colors, gp = gpar(col = "white"), simple_anno_size = unit(0.15, "cm"), show_annotation_name = FALSE, show_legend = FALSE,
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
    #maxi = max(df)
    #col_fun = colorRamp2(c(0, max_value), c(viridis(100)[1],viridis(100)[100]))
    col_fun = colorRamp2(c(0, max_value), c("#FFF6D6","#D5A021"))
    #col_fun = colorRamp2(c(0,max_value), hcl_palette = "Greens", reverse = TRUE)
    col_fun_legend = col_fun
  } else {
    #maxi = max(df)
    #col_fun = colorRamp2(c(0,max_value), hcl_palette = "Greens", reverse = TRUE)
    col_fun = colorRamp2(c(0,max_value), c("#FFF6D6","#D5A021"), reverse = TRUE)
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
  write.table(df_bin, sprintf("%s_plot_df.csv", filename_tab), sep = "\t", row.names = TRUE, col.names = TRUE)
  write.xlsx(df_bin, sprintf("%s_plot_df.xlsx", filename_tab), colNames = TRUE, rowNames = TRUE, append = FALSE)
  
  if ((show_legend == FALSE) && (legend_position != "right")) {
    height_heatmap = height + 2
  } else {
    height_heatmap = height
  }
  
  # remove mmu-miR
  rownames(df_bin) = gsub("mmu-miR-", "", rownames(df_bin))

  # features rows
  p = ComplexHeatmap::Heatmap(df_bin, col = col_fun,
                              rect_gp = gpar(type = "none"), 
                              #height = unit(0.005, "cm") * nrow(df), width = unit(0.3, "cm") * ncol(df),
                              height = unit(height_heatmap + 4.5, "cm"), width = unit(width - 4.5, "cm"),
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
  if (show_legend == TRUE) { 
    if (length(unique(as.numeric(df))) <= 2) {
      # discrete legend
      legend_ticks = c(0, max_value)
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
  }
  #png(sprintf("%s_row=%s_portrait_mod_single.png", filename_fig, data_input$rna_class), width = width, height = 2 * height, units = plot_props$image_units, res = plot_props$dpi)
  #if (show_legend == TRUE) {
  #  draw(p, show_annotation_legend = FALSE, padding = unit(c(1.5, 0.25, 0, 0.25), "cm")) #top, right, bottom, left
  #  if (legend_position == "bottom") {
  #    draw(legend, x = unit(1.3, "cm"), y = unit(0.02, "cm"), just = c("left", "bottom"))
  #  } else if (legend_position == "right") {
  #    print("legend: right")
  #    draw(legend, x = unit(plot_props$image_width - 1.25, "cm"), y = unit(plot_props$image_height, "cm"))
  #  }
  #} else {
  #  draw(p, show_annotation_legend = FALSE, padding = unit(c(0.25, 0.25, 0, 0.25), "cm"))
  #}
  #dev.off()
  #svglite(sprintf("%s_row=%s_portrait_mod.svg", filename_fig, data_input$rna_class), width = width / 2.54, height = 2 * height / 2.54)
  #if (show_legend == TRUE) {
  #  draw(p, show_annotation_legend = FALSE, padding = unit(c(1.5, 0.25, 0, 0.25), "cm"))
  #  if (legend_position == "bottom") {
  #    draw(legend, x = unit(1.3, "cm"), y = unit(0.02, "cm"), just = c("left", "bottom"))
  #  } else if (legend_position == "right") {
  #    draw(legend, x = unit(1.3, "cm"), y = unit(0.02, "cm"), just = c("right"))
  #  }
  #} else {
  #  draw(p, show_annotation_legend = FALSE, padding = unit(c(0.25, 0.25, 0, 0.25), "cm"))
  #}
  #dev.off()
  
  return(list("heatmap" = p, "legend" = legend, "row_names" = rownames(df_bin)))
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
xticks_names = snakemake@params$xticks_names
cluster_colours = snakemake@params$colors
plot_props = snakemake@params$plot_props
results_folder = snakemake@params$results_folder

# save rdata
#saveRDS(snakemake, file = "snakemake_hclust_expression_brain_region_for_mirna_list_zscores_heatmap_plot_male_female.rds")

# cluster colours
cluster_colours_list = unlist(cluster_colours$clusters)
#tmp = names(cluster_colours_list) %in% names(features_cluster)
#cluster_colours_list = cluster_colours_list[tmp]
annot = c()
print(data_input)
expr_male = fread(sprintf("%s_%s/%s_%s_quantification_%s.detection_rate_%s_per_%s.csv", data_input$raw_data_folder, data_input$data_sub_set_male, data_input$rna_class, data_input$feature_filtering, data_input$norm, data_input$detection_rate, data_input$detection_group), sep='\t', header=T)
annot$male = fread(sprintf("%s/annotation_%s_%s.csv", data_input$annotation_data_folder, data_input$data_sub_set_male, data_input$feature_filtering), sep='\t', colClasses=c(ID="character")) 
diff_exp_male = fread(sprintf("%s/%s_%s/results_%s/matrices/diff_exp/diff_exp%s.csv", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set_male, diff_log), sep='\t') 

expr_female = fread(sprintf("%s_%s/%s_%s_quantification_%s.detection_rate_%s_per_%s.csv", data_input$raw_data_folder, data_input$data_sub_set_female, data_input$rna_class, data_input$feature_filtering, data_input$norm, data_input$detection_rate, data_input$detection_group), sep='\t', header=T)
annot$female = fread(sprintf("%s/annotation_%s_%s.csv", data_input$annotation_data_folder, data_input$data_sub_set_female, data_input$feature_filtering), sep='\t', colClasses=c(ID="character")) 
diff_exp_female = fread(sprintf("%s/%s_%s/results_%s/matrices/diff_exp/diff_exp%s.csv", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set_female, diff_log), sep='\t') 


#------------------------------------ Script ----------------------------------- 
method_split = data.frame(stri_split_fixed(method, "_"))[1,]
for (m in method_split) {
  print(m)
  output_folder_fig = sprintf("%s/%s_%s/results_%s_%s/figures/%s_clustering/heatmap_complex", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set_male, data_input$data_sub_set_female, m)
  dir.create(output_folder_fig, recursive=TRUE)
  output_folder_tab = sprintf("%s/%s_%s/results_%s_%s/matrices/%s_clustering", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set_male, data_input$data_sub_set_female, m)
  dir.create(output_folder_tab, recursive=TRUE)
  
  for(i in top_list){
    print(i)
    feature_list = c()
    feature_list$male = fread(sprintf("%s/%s_%s/results_%s/matrices/%s/%s.csv", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set_male, m, i))
    feature_list$male = feature_list$male[[data_input$rna_class]]
    
    feature_list$female = fread(sprintf("%s/%s_%s/results_%s/matrices/%s/%s.csv", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set_female, m, i))
    feature_list$female = feature_list$female[[data_input$rna_class]]
    
    # filter expr with feature_list
    diff_exp_filtered = c()
    diff_exp_filtered$male = diff_exp_male[diff_exp_male[[data_input$rna_class]] %in% feature_list$male,]
    diff_exp_filtered$female = diff_exp_female[diff_exp_female[[data_input$rna_class]] %in% feature_list$female,]
    
    median_table = c()
    median_table$male = data.table(feature=diff_exp_filtered[["male"]][[data_input$rna_class]])
    for (ti in unique(annot[["male"]][[prop]])) {
      median_table[["male"]][[ti]] = 10^(as.numeric(diff_exp_filtered[["male"]][[sprintf("median_%s__%s", prop, ti)]]))
    }
    median_table$female = data.table(feature=diff_exp_filtered[["female"]][[data_input$rna_class]])
    for (ti in unique(annot[["female"]][[prop]])) {
      median_table[["female"]][[ti]] = 10^(as.numeric(diff_exp_filtered[["female"]][[sprintf("median_%s__%s", prop, ti)]]))
    }
    
    # save rnas and create row names for plotting
    rna_ids = c()
    rna_ids$male = median_table$male[[colnames(median_table$male)[1]]]
    rna_ids$female = median_table$female[[colnames(median_table$female)[1]]]
    #if(i <= 50) {
    row_names = c()
    row_names$male = rna_ids$male
    row_names$female = rna_ids$female
    #} else {
    #  row_names = NULL
    #}
    
    # keep only expression
    median_table_adj = c()
    median_table_adj$male = median_table[["male"]][, sapply(median_table[["male"]], is.numeric), with=F]
    median_table_adj$female = median_table[["female"]][, sapply(median_table[["female"]], is.numeric), with=F]
    
    p_left = list()
    p_right = list()
    p_legend_left = list()
    p_legend_right = list()
    
    
    for (zscore_th in zscore_th_list) {
      
      median_table_adj_th_removed_zeros = list()
      median_table_adj_removed_zeros = list()
      row_names_removed_zeros = list()
      
      for (sex_key in c("male", "female")) {
        #print(zscore_th)
        
        # calculate zscores
        params$scale = "row"
        if(params$scale != "none"){
          scale_row = params$scale %in% c("row", "both")
          scale_col = params$scale %in% c("column", "both")
          median_table_adj_zscores = scale(median_table_adj[[sex_key]], scale_row, scale_col) #base::scale(t(median_table_adj[i,]))
          median_table_adj_zscores[is.na(median_table_adj_zscores)] = 0
          median_table_adj_zscores = abs(median_table_adj_zscores)
        }
        
        #print(median_table_adj)
        median_table_adj_zscores_th = ifelse(abs(median_table_adj_zscores) >= zscore_th, 1, 0)
  
        #print(median_table_adj_zscores_th)
        legend_title = sprintf("Bin. stand. expr. value\n(median per %s)", prop)
        #plot_heatmap_row_feature_portrait(median_table_adj_zscores_th, row_names[[sex_key]], annot[[sex_key]], prop, params, legend_title, cluster_colours_list, plot_props, width = plot_props$image_width, height = plot_props$image_height, filename_fig = sprintf("%s/%s_zscore_th=%s_num_clusters=%s", output_folder_fig, i, zscore_th, params$number_of_cluster), filename_tab = sprintf("%s/%s_zscore_th=%s_clustered_by_%s", output_folder_tab, i, zscore_th, m))
        
        rownames(median_table_adj_zscores_th) = row_names[[sex_key]]
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
        median_table_adj_th_removed_zeros[[sex_key]] = median_table_adj_zscores_th[mask, , drop=FALSE]
        row_names_removed_zeros[[sex_key]] = row_names[[sex_key]][mask]
        #plot_heatmap_row_feature_portrait(median_table_adj_th_removed_zeros[[sex_key]], row_names_removed_zeros[[sex_key]], annot[[sex_key]], prop, params, legend_title, cluster_colours_list, plot_props, width = plot_props$image_width, height = plot_props$image_height, filename_fig = sprintf("%s/%s_zscore_th=%s_num_clusters=%s_removed_zero_rows", output_folder_fig, i, zscore_th, params$number_of_cluster), filename_tab = sprintf("%s/%s_zscore_th=%s_clustered_by_%s_removed_zero_rows", output_folder_tab, i, zscore_th, m))
        
        # removed rows containing only zeros after binarisation
        # red border around expression values higher than the threhold
        median_table_adj_removed_zeros[[sex_key]] = median_table_adj_zscores[mask, , drop=FALSE]
        
      }
      
      max_value = ceiling(max(max(median_table_adj_removed_zeros$male), max(median_table_adj_removed_zeros$female)))
      
      for (sex_key in c("male", "female")) {
        
        print(sex_key)
        
        legend_title = sprintf("Abs. stand.\nexpr. (%s,\nmedian)", gsub("_norm", "", data_input$norm))
        if (nrow(median_table_adj_removed_zeros[[sex_key]]) > 0) {
          show_legend = TRUE
          legend_position = "right"
          rownames_pos = "left"
          p_legend_left[[sex_key]] = plot_heatmap_row_feature_portrait_mod(median_table_adj_removed_zeros[[sex_key]], median_table_adj_th_removed_zeros[[sex_key]], annot[[sex_key]], row_names_removed_zeros[[sex_key]], max_value, prop, params, legend_title, cluster_colours_list, plot_props, rownames_pos, show_legend, legend_position, width = plot_props$image_width, height = plot_props$image_height, filename_fig = sprintf("%s/%s_zscore_th=%s_num_clusters=%s_removed_zero_rows_rownames_pos=%s_show_legend=%s", output_folder_fig, i, zscore_th, params$number_of_cluster, rownames_pos, show_legend), filename_tab = sprintf("%s/%s_zscore_th=%s_clustered_by_%s_removed_zero_rows", output_folder_tab, i, zscore_th, m))
          rownames_pos = "right"
          p_legend_right[[sex_key]] = plot_heatmap_row_feature_portrait_mod(median_table_adj_removed_zeros[[sex_key]], median_table_adj_th_removed_zeros[[sex_key]], annot[[sex_key]], row_names_removed_zeros[[sex_key]], max_value, prop, params, legend_title, cluster_colours_list, plot_props, rownames_pos, show_legend, legend_position, width = plot_props$image_width, height = plot_props$image_height, filename_fig = sprintf("%s/%s_zscore_th=%s_num_clusters=%s_removed_zero_rows_rownames_pos=%s_show_legend=%s", output_folder_fig, i, zscore_th, params$number_of_cluster, rownames_pos, show_legend), filename_tab = sprintf("%s/%s_zscore_th=%s_clustered_by_%s_removed_zero_rows", output_folder_tab, i, zscore_th, m))
          
          show_legend = FALSE
          rownames_pos = "left"
          p_left[[sex_key]] = plot_heatmap_row_feature_portrait_mod(median_table_adj_removed_zeros[[sex_key]], median_table_adj_th_removed_zeros[[sex_key]], annot[[sex_key]], row_names_removed_zeros[[sex_key]], max_value, prop, params, legend_title, cluster_colours_list, plot_props, rownames_pos, show_legend, legend_position, width = plot_props$image_width, height = plot_props$image_height, filename_fig = sprintf("%s/%s_zscore_th=%s_num_clusters=%s_removed_zero_rows_rownames_pos=%s_show_legend=%s", output_folder_fig, i, zscore_th, params$number_of_cluster, rownames_pos, show_legend), filename_tab = sprintf("%s/%s_zscore_th=%s_clustered_by_%s_removed_zero_rows", output_folder_tab, i, zscore_th, m))
          rownames_pos = "right"
          p_right[[sex_key]] = plot_heatmap_row_feature_portrait_mod(median_table_adj_removed_zeros[[sex_key]], median_table_adj_th_removed_zeros[[sex_key]], annot[[sex_key]], row_names_removed_zeros[[sex_key]], max_value, prop, params, legend_title, cluster_colours_list, plot_props, rownames_pos, show_legend, legend_position, width = plot_props$image_width, height = plot_props$image_height, filename_fig = sprintf("%s/%s_zscore_th=%s_num_clusters=%s_removed_zero_rows_rownames_pos=%s_show_legend=%s", output_folder_fig, i, zscore_th, params$number_of_cluster, rownames_pos, show_legend), filename_tab = sprintf("%s/%s_zscore_th=%s_clustered_by_%s_removed_zero_rows", output_folder_tab, i, zscore_th, m))
        }
      }
      
      #p = p_right[["male"]]$heatmap + p_left[["female"]]$heatmap
      
      # get both plots as grob objects
      grob_male = grid.grabExpr(draw(p_right[["male"]]$heatmap))
      grob_female = grid.grabExpr(draw(p_left[["female"]]$heatmap))
      legend_obj = grid.grabExpr(draw(p_legend_left[["female"]]$legend))
      # put both side by side (with cowplot) (see https://jokergoo.github.io/2023/04/03/align-heatmaps/)
      #p = plot_grid(plot_grid(grob_male, grob_female, ncol=2, align='h'), plot_grid(NULL, p_legend_left[["female"]]$legend, ncol=1)) #, rel_widths=c(1, 0.2))
      p_18x12 = plot_grid(NULL, grob_male, grob_female, legend_obj, ncol=4, align='h', rel_widths=c(0.3, 1, 1, 0.3))
      p_15x12 = plot_grid(grob_male, grob_female, legend_obj, ncol=3, align='h', rel_widths=c(1, 1, 0.3))
      
      #p = plot_grid(grob_male, grob_female, ncol=2, align='h')
      #draw(legend, x = unit(plot_props$image_width - 1.25, "cm"), y = unit(plot_props$image_height, "cm"))
      
      # creating concatinated plot
      # save
      # figures
      # calculate x shift 
      #tmp = unlist(xticks_names$brain_region[unique(annot[[sex_key]][[prop]])])
      #x_shift = (max(nchar(tmp)) * 0.17) + 0.2
      
      filename_fig = sprintf("%s/%s_zscore_th=%s_num_clusters=%s_removed_zero_show_legend=TRUE", output_folder_fig, i, zscore_th, params$number_of_cluster)
      file_name = sprintf("%s_row=%s_portrait_mod", filename_fig, data_input$rna_class)
      #png(sprintf("%s.png", file_name), width = 2*plot_props$image_width, height = 2*plot_props$image_height, units = plot_props$image_units, res = plot_props$dpi)
      #ggdraw() + draw_plot(p)
      #dev.off()
      #svglite(sprintf("%s.svg", file_name), width = 2*(plot_props$image_width / 2.54), height = (plot_props$image_height / 2.54))
      #ggdraw() + draw_plot(p)
      #dev.off()
      
      save_plot(plot=p_18x12, filename=sprintf("%s_18x12.png", file_name), nrow=1, ncol=1, base_width = 2*plot_props$image_width, base_height = 2*plot_props$image_height, units = plot_props$image_units, dpi = plot_props$dpi)
      save_plot(plot=p_18x12, filename=sprintf("%s_18x12.svg", file_name), nrow=1, ncol=1, base_width = 2*plot_props$image_width, base_height = 2*plot_props$image_height, units = plot_props$image_units)
      
      save_plot(plot=p_15x12, filename=sprintf("%s_15x12.png", file_name), nrow=1, ncol=1, base_width = 2*plot_props$image_width-3, base_height = 2*plot_props$image_height, units = plot_props$image_units, dpi = plot_props$dpi)
      save_plot(plot=p_15x12, filename=sprintf("%s_15x12.svg", file_name), nrow=1, ncol=1, base_width = 2*plot_props$image_width-3, base_height = 2*plot_props$image_height, units = plot_props$image_units)
      
      #filename_fig = sprintf("%s/%s_zscore_th=%s_num_clusters=%s_removed_zero_show_legend=FALSE", output_folder_fig, i, zscore_th, params$number_of_cluster)
      #file_name = sprintf("%s_row=%s_portrait_mod", filename_fig, data_input$rna_class)
      #png(sprintf("%s.png", file_name), width = 2*plot_props$image_width, height = plot_props$image_height, units = plot_props$image_units, res = plot_props$dpi)
      #draw(p, padding = unit(c(0.25, 0.25-x_shift, 0.25, 0.25), "cm"), ht_gap = unit(c(0, 15, 0, 0), "mm"))  # bottom, left, top and right margins
      #draw(p_legend_left[["female"]]$legend, x = unit((2*plot_props$image_width - 1), "cm"), y = unit((plot_props$image_height / 2), "cm"), just = c("right", "center"))
      #dev.off()
      #svglite(sprintf("%s.svg", file_name), width = 2*(plot_props$image_width / 2.54), height = (plot_props$image_height / 2.54))
      #draw(p, padding = unit(c(0.25, 0.25-x_shift, 0.25, 0.25), "cm"), ht_gap = unit(c(0, 15, 0, 0), "mm"))  # bottom, left, top and right margins
      #draw(p_legend_left[["female"]]$legend, x = unit((2*plot_props$image_width - 1), "cm"), y = unit((plot_props$image_height / 2), "cm"), just = c("right", "center"))
      #dev.off()
      
      filename_tab = sprintf("%s/%s_zscore_th=%s_removed_zero", output_folder_tab, i, zscore_th)
      file_name_tab = sprintf("%s_overlapping_%ss", filename_tab, data_input$rna_class)
      
      overlapping_rnas = data.frame(features = intersect(p_left[["male"]]$row_names, p_left[["female"]]$row_names))
      fwrite(overlapping_rnas, sprintf("%s.csv", file_name_tab), sep = "\t", row.names = FALSE)
      write.xlsx(overlapping_rnas, sprintf("%s.xlsx", file_name_tab), colNames = TRUE, rowNames = FALSE, append = FALSE)
    }
  }
}


#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])

