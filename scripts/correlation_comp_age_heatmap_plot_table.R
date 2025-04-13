suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(svglite))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(openxlsx))
#library(corrplot)
suppressPackageStartupMessages(library(Hmisc))

# save rdata
#saveRDS(snakemake, file = "snakemake_correlation_comp_age_heatmap_plot_table.rds")
#stop()

#snakemake = readRDS("snakemake_correlation_comp_age_heatmap_plot_table.rds")

set.seed(snakemake@params$parameters_props$set_seed)

print("correlation_comp_age_heatmap_plot_table")


#---------------------------------- Functions ----------------------------------
plot_heatmap_annot = function(df, pvalues_df,  row_names, xticks_names, plots_props, prop, time, m, annotation_row, data_input, colors, th) {
  df = t(df)
  df[abs(df) <= th] = 0
  df[df < -th] = -0.5
  df[df > th] = 0.5
  
  color_bar_name = sprintf("%s correlation (%s)", time, m) 
  
  #cmap = colorRamp2(c(-1, 0, 1), hcl_palette = "Green-Brown", reverse = TRUE)
  #col_val = 0.75
  #col_fun = c(cmap(-col_val), cmap(0), cmap(col_val))
  #col_legend = c(cmap(-col_val), cmap(col_val))
  col_fun = c(colors$direction$neg, colors$direction$light, colors$direction$pos)
  col_legend = c(colors$direction$neg, colors$direction$pos)
  names(col_legend) = c(sprintf("negative (R ≤ -%s)", th), sprintf("positive (R ≥ %s)", th))
  
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (is.na(t(pvalues_df)[i, j])) {
    }
    else if ((t(pvalues_df)[i, j] < 0.05) & (df[i,j] != 0)){ 
      grid.text("*", x, y, gp = gpar(fontsize = 2, col = "white"))
    #} else if (t(pvalues_df)[i, j] < 0.01){ 
    #  grid.text("**", x, y, gp = gpar(fontsize = 4))
    #} else if (t(pvalues_df)[i, j] < 0.001){ 
    #  grid.text("***", x, y, gp = gpar(fontsize = 4))
    }
  }
  
  #annotation_row = data.frame(dummy=1:nrow(df))
  #annotation_row[[prop]] = as.character(annot[[prop]][match(rownames(df), annot[[prop]])])
  #colnames(annotation_row) = rownames(df)
  #annotation_row$dummy = NULL
  
  annotation_colors = list()
  ucolors <- unlist(colors[[prop]])
  if(is.vector(ucolors) && all(annotation_row[[prop]] %in% names(ucolors))) {
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
  
  annotation_row = annotation_row[rownames(df), , drop=FALSE]
  
  colnames(df) = row_names
  rownames(df) = unlist(xticks_names[[tissue]][rownames(df)])
  rownames(annotation_row) = unlist(xticks_names[[tissue]][rownames(annotation_row)])
  #print(annotation_row)
  #print(annotation_colors)
  column_ha = rowAnnotation(df = annotation_row, col = annotation_colors, gp = gpar(row = "white"), simple_anno_size = unit(0.2, "cm"), 
                            annotation_name_side = "top", annotation_name_gp = gpar(fontsize = plots_props$font_size, fontface = "bold", fontfamily=plots_props$font_family),
                            show_legend=FALSE)
  
  p = ComplexHeatmap::Heatmap(df, col = col_fun, cell_fun = cell_fun,
                              #rect_gp = gpar(col = "white", lwd = .5),
                              # height = unit(0.3, "cm") * nrow(t(df)), width = unit(0.3, "cm") * ncol(t(df)),
                              #height = unit(plots_props$image_height  - 2, "cm"), width = unit(plots_props$image_width - 4, "cm"),
                              height = unit(2, "cm"), width = unit(6.25*3, "cm"),
                              show_column_names = FALSE,
                              row_names_side = "left", row_names_gp = gpar(fontsize = 6, fontfamily=plots_props$font_family), 
                              #column_names_side = "bottom", column_names_gp = gpar(fontsize = plots_props$font_size, fontfamily=plots_props$font_family),
                              cluster_columns = function(m) hclust(dist(m), method = "complete"), show_column_dend = TRUE, clustering_distance_columns = "euclidean", 
                              #cluster_rows = function(m) hclust(dist(m), method = "complete"), show_row_dend = TRUE, clustering_distance_rows = "euclidean", 
                              column_dend_reorder = TRUE, #cluster_columns = TRUE,
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
  
  legend = ComplexHeatmap::Legend(title = color_bar_name, labels = names(col_legend), legend_gp = gpar(fill = col_legend),
                                  direction="horizontal", ncol = 2,
                                  title_gp = gpar(fontsize = plots_props$font_size, fontface = "bold", fontfamily=plots_props$font_family), labels_gp = gpar(fontsize = plots_props$font_size, fontfamily=plots_props$font_family),
                                  legend_width=unit(4, "cm"))
  
  #p = pheatmap(t(df), color = viridis(100), border_color = "white", show_colnames = TRUE, show_rownames = TRUE,
  #             fontsize = plots_props$font_size, fontsize_row = 6, fontsize_col = plots_props$font_size, 
  #             cellwidth = 8, cellheight = 8,
  #             cluster_cols = FALSE, clustering_distance_cols = "euclidean", cluster_rows = TRUE, clustering_distance_rows = "euclidean", treeheight_row = 0
  #)
  
  return(list("heatmap" = p, "legend" = legend))
}

plot_heatmap_annot_filtered = function(df, pvalues_df, row_names, data_input, xticks_names, plots_props, prop, time, m, annotation_row, feature_list, colors, th) {

  df = t(df)
  df[abs(df) < th] = 0
  df[df <= -th] = -0.5
  df[df >= th] = 0.5
  
  pvalues_df = t(pvalues_df)
  pvalues_df = pvalues_df[, colSums(df != 0) > 0]
  row_names = row_names[colSums(df != 0) > 0]
  df = df[, colSums(df != 0) > 0]
  row_names = row_names[colSums(is.na(df)) != nrow(df)]
  df = df[, colSums(is.na(df)) != nrow(df)]
  
  tmp = !(row_names %in% feature_list[[data_input$feature_column]])
  row_names[tmp] = ""
  df_label = matrix("|", nrow(df), ncol(df))
  df_label[,tmp] = ""

  color_bar_name = sprintf("%s correlation (%s)", m, time) 
  
  #cmap = colorRamp2(c(-1, 0, 1), hcl_palette = "Green-Brown", reverse = TRUE)
  #col_val = 0.75
  #col_fun = c(cmap(-col_val), cmap(0), cmap(col_val))
  #col_legend = c(cmap(-col_val), cmap(col_val))
  col_fun = c(colors$direction$neg, colors$direction$light, colors$direction$pos)
  col_legend = c(colors$direction$neg, colors$direction$pos)
  names(col_legend) = c(sprintf("negative (R ≤ -%s)", th), sprintf("positive (R ≥ %s)", th))
  
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (is.na(pvalues_df[i, j])) {
    }
    else if ((pvalues_df)[i, j] < 0.05 & (df[i,j] != 0)){ 
      grid.text(sprintf("*%s", df_label[i, j]), x, y, gp = gpar(fontsize = 2, col = "white"))
      #} else if (pvalues_df[i, j] < 0.01){ 
      #  grid.text("**", x, y, gp = gpar(fontsize = 4))
      #} else if (pvalues_df[i, j] < 0.001){ 
      #  grid.text("***", x, y, gp = gpar(fontsize = 4))
    }
  }
  
  #annotation_row = data.frame(dummy=1:nrow(df))
  #annotation_row[[prop]] = as.character(annot[[prop]][match(rownames(df), annot[[prop]])])
  #colnames(annotation_row) = rownames(df)
  #annotation_row$dummy = NULL
  
  annotation_colors = list()
  ucolors <- unlist(colors[[prop]])
  if(is.vector(ucolors) && all(annotation_row[[prop]] %in% names(ucolors))) {
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
  annotation_row = annotation_row[rownames(df), , drop=FALSE]
  
  colnames(df) = row_names
  rownames(df) = unlist(xticks_names[[tissue]][rownames(df)])
  rownames(annotation_row) = unlist(xticks_names[[tissue]][rownames(annotation_row)])
  
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
  
  legend = ComplexHeatmap::Legend(title = color_bar_name, labels = names(col_legend), legend_gp = gpar(fill = col_legend),
                                  direction="horizontal", ncol = 2,
                                  title_gp = gpar(fontsize = plots_props$font_size, fontface = "bold", fontfamily=plots_props$font_family), labels_gp = gpar(fontsize = plots_props$font_size, fontfamily=plots_props$font_family),
                                  legend_width=unit(4, "cm"))
  
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


#---------------------------------- Load data ---------------------------------- 
data_input = snakemake@params$data_input
method = snakemake@params$method
props = snakemake@params$props
file = snakemake@params$file
xticks_names = snakemake@params$xticks_names
plots_props = snakemake@params$plots_props
colors = snakemake@params$colors
results_folder = snakemake@params$results_folder
adjustment = snakemake@params$adjustment
ID = data_input$identifier_column
corr_ths = snakemake@params$corr_th

expr = fread(sprintf("%s_%s/%s_%s_quantification_%s.detection_rate_%s_per_%s.csv", data_input$raw_data_folder, data_input$data_sub_set, data_input$rna_class, data_input$feature_filtering, data_input$norm, data_input$detection_rate, data_input$detection_group), sep='\t', header=T)
annot = fread(sprintf("%s/annotation_%s_%s.csv", data_input$annotation_data_folder, data_input$data_sub_set, data_input$feature_filtering), sep='\t', colClasses=c(ID="character")) 
fixed_cols = c(data_input$feature_column, annot[[data_input$identifier_column]])
expr = expr[, ..fixed_cols]

feature_list = fread(file, sep='\t')  # , colClasses=c(ID="character")
feature_list = feature_list[FALSE,]


#------------------------------------ Script ----------------------------------- 
for (m in method) {
  for (prop in props) {
    tissue = prop[1]
    time = prop[2]
    
    output_folder_fig = sprintf("%s/%s_%s/results_%s/figures/corr_plots/heatmap_complex/all_%s_with_%s", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, data_input$rna_class, time)
    dir.create(output_folder_fig, recursive=TRUE)
    output_folder_tab = sprintf("%s/%s_%s/results_%s/matrices/aggregation_corr", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set)
    dir.create(output_folder_tab, recursive=TRUE)
    
    correlation = c()
    corr_pvalues = c()
    corr_pvalues_raw = c()
    for(group in unique(annot[[tissue]])){
      groups_ID = annot[annot[[tissue]] == group,][[ID]]
      tmp = colnames(expr) %in% groups_ID
      expr_sort = data.frame(feature=expr[[data_input$feature_column]], expr[,..tmp], check.names = FALSE)
      tmp = annot[[ID]] %in% groups_ID
      
      timepoints = annot[tmp,][[time]]
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
        correlation[[group]] = as.data.frame(cor(t(expr_sort[,2:dim(expr_sort)[2]]), numeric_timepoints, method = m))
        # if all expr for every timepoint for a sample is 0 for a feature
        # then the correlation is NA
        # we then force the correlation to be 0
        na_cor = is.na(correlation[[group]])
        # uncomment this to show the corresponding lines in the expr_matrix
        # print(expr_sort[,2:dim(expr_sort)[2]][na_cor, ])
        # replace NA with 0
        correlation[[group]][na_cor, ] = 0
        # get pvalues from rcorr
        # $P gets p-values
        pvalues_tmp = as.data.frame(rcorr(t(expr_sort[,2:dim(expr_sort)[2]]), numeric_timepoints, type=m)$P)
        # select y because x=matrix, y=timepoints
        # and remove last value because it is correlation of timepoints with timepoints = NA
        corr_pvalues[[group]] = head(pvalues_tmp$y, -1)
        corr_pvalues_raw[[group]] = corr_pvalues[[group]]
        # adjustment
        corr_pvalues[[group]] = as.data.frame(p.adjust(corr_pvalues[[group]], method=adjustment))
        # remove all pvalues == NA and replace them with 1 (then it is not significant)
        na_pvalue = is.na(corr_pvalues[[group]])
        corr_pvalues[[group]][na_pvalue] = 1
        # check if anything is still NA
        #print(sprintf("NA in correlation: %s", any(is.na(correlation[[group]]))))
        #print(sprintf("NA in pvalues: %s", any(is.na(corr_pvalues[[group]]))))
      }
    }
    
    correlation_df = as.data.frame(correlation)
    colnames(correlation_df) = unique(annot[[tissue]])
    
    corr_pvalues_raw_df = as.data.frame(corr_pvalues_raw)
    colnames(corr_pvalues_raw_df) = colnames(correlation_df)
    
    corr_pvalues_df = as.data.frame(corr_pvalues)
    colnames(corr_pvalues_df) = colnames(correlation_df)
    
    annot_row = data.frame(Tissue=unique(annot[[tissue]]))
    # rename to brain_region this case
    colnames(annot_row) = c(tissue)
    rownames(annot_row) = unique(annot[[tissue]])
    
    # order cols according to color config
    col_order_fixed = names(colors[[tissue]])
    # find the ones actually present
    col_order_fixed = col_order_fixed[col_order_fixed %in% colnames(correlation_df)]
    
    for (th in corr_ths) {
      # all features
      heatmap_plot = plot_heatmap_annot(correlation_df[, col_order_fixed], corr_pvalues_df[, col_order_fixed], expr_sort[["feature"]], xticks_names, plots_props, tissue, time,  m, annot_row, data_input, colors, th)
      # only features for which one group fulfill the threshold
      heatmap_plot_filtered = plot_heatmap_annot_filtered(correlation_df[, col_order_fixed], corr_pvalues_df[, col_order_fixed], expr_sort[["feature"]], data_input, xticks_names, plots_props, tissue, time, m, annot_row, feature_list, colors, th)
    }
 
    row_name_square_plot = row_name_squares(col_order_fixed, xticks_names, tissue, colors, plots_props)
    
    heatmap_plot_comb = row_name_square_plot + heatmap_plot$heatmap
    heatmap_plot_filtered_comb = row_name_square_plot + heatmap_plot_filtered$heatmap
    
    # save figures
    filename = sprintf("%s/corr_method=%s_%s_over_%s_per_%s", output_folder_fig, m, data_input$rna_class, time, tissue)
    png(sprintf("%s_corr_th=%s.png", filename, th), width = plots_props$image_width*2.5, height = plots_props$image_height*1.25, units = plots_props$image_units, res = plots_props$dpi)
    draw(heatmap_plot_comb, padding = unit(c(0.25, 0.25, -1, 0.25), "cm"))  # bottom, left, top and right margins
    draw(heatmap_plot$legend, x = unit(.5, "cm"), y = unit(.1, "cm"), just = c("left", "bottom"))
    dev.off()
    svglite(sprintf("%s_corr_th=%s.svg", filename, th), width = plots_props$image_width*2.5 / 2.54, height = plots_props$image_height*1.25 / 2.54)
    draw(heatmap_plot_comb, padding = unit(c(0.25, 0.25, -1, 0.25), "cm"))  # bottom, left, top and right margins
    draw(heatmap_plot$legend, x = unit(.5, "cm"), y = unit(.1, "cm"), just = c("left", "bottom"))
    dev.off()
    
    filename = sprintf("%s/corr_method=%s_%s_over_%s_per_%s_filtered", output_folder_fig, m, data_input$rna_class, time, tissue)
    png(sprintf("%s_corr_th=%s.png", filename, th), width = plots_props$image_width, height = plots_props$image_height, units = plots_props$image_units, res = plots_props$dpi)
    draw(heatmap_plot_filtered_comb, padding = unit(c(-0.5, 0.25, -1.5, 0.25), "cm"))  # bottom, left, top and right margins
    draw(heatmap_plot_filtered$legend, x = unit(.5, "cm"), y = unit(.1, "cm"), just = c("left", "bottom"))
    dev.off()
    svglite(sprintf("%s_corr_th=%s.svg", filename, th), width = plots_props$image_width / 2.54, height = plots_props$image_height / 2.54)
    draw(heatmap_plot_filtered_comb, padding = unit(c(-0.5, 0.25, -1.5, 0.25), "cm"))  # bottom, left, top and right margins
    draw(heatmap_plot_filtered$legend, x = unit(.5, "cm"), y = unit(.1, "cm"), just = c("left", "bottom"))
    dev.off()
    
    # save tables
    correlation_df = data.frame(feature = expr[[data_input$feature_column]], correlation_df)
    colnames(correlation_df) = gsub("feature", paste(data_input$rna_class), colnames(correlation_df))
    fwrite(correlation_df, sprintf("%s/corr_method=%s_%ss_with_%s_per_%s.csv", output_folder_tab, m, data_input$rna_class, time, tissue), sep = "\t", row.names = TRUE)
    write.xlsx(correlation_df, sprintf("%s/corr_method=%s_%ss_with_%s_per_%s.xlsx", output_folder_tab, m, data_input$rna_class, time, tissue), colNames = TRUE, rowNames = TRUE, append = FALSE)
    
    corr_pvalues_raw_df = data.frame(feature = expr[[data_input$feature_column]], corr_pvalues_raw_df)
    colnames(corr_pvalues_raw_df) = gsub("feature", paste(data_input$rna_class), colnames(corr_pvalues_raw_df))
    fwrite(corr_pvalues_raw_df, sprintf("%s/corr_method=%s_pvalues_%ss_with_%s_per_%s.csv", output_folder_tab, m, data_input$rna_class, time, tissue), sep = "\t", row.names = TRUE)
    write.xlsx(corr_pvalues_raw_df, sprintf("%s/corr_method=%s_pvalues_%ss_with_%s_per_%s.xlsx", output_folder_tab, m, data_input$rna_class, time, tissue), colNames = TRUE, rowNames = TRUE, append = FALSE)
    
    corr_pvalues_df = data.frame(feature = expr[[data_input$feature_column]], corr_pvalues_df)
    colnames(corr_pvalues_df) = gsub("feature", paste(data_input$rna_class), colnames(corr_pvalues_df))
    fwrite(corr_pvalues_df, sprintf("%s/corr_method=%s_adj_pvalues_%ss_with_%s_per_%s.csv", output_folder_tab, m, data_input$rna_class, time, tissue), sep = "\t", row.names = TRUE)
    write.xlsx(corr_pvalues_df, sprintf("%s/corr_method=%s_adj_pvalues_%ss_with_%s_per_%s.xlsx", output_folder_tab, m, data_input$rna_class, time, tissue), colNames = TRUE, rowNames = TRUE, append = FALSE)
  }
  
  print(sprintf("Done: %s", m))
}


#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])
