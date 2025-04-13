suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(Hmisc))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(svglite))
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(foreach))

#snakemake = readRDS("snakemake_correlation_mirna_human_with_age_heatmap_plot.rds")

set.seed(snakemake@params$parameters_props$set_seed)

print("correlation_mirna_human_with_age_heatmap_plot")


#---------------------------------- Functions ----------------------------------
plot_heatmap_annot = function(df, pvalues_df, row_names, group_size_df, data_input, xticks_names, plots_props, time, m, annotation_row, colors, th, sig_lvl) {
  
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
    color_bar_name = sprintf("Spearman correlation (age)") 
  } else if (m == "pearson") {
    color_bar_name = sprintf("Pearson correlation (age)") 
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
      grid.rect(x, y, width, height, gp = gpar(fill=fill, col = "black", lwd = 0.3))  # Black cell border
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
  
  #annotation_colors = list()
  #ucolors = unlist(colors[[prop]])
  #if(is.vector(ucolors) && all(annotation_row %in% names(ucolors))) {
  #  annotation_colors[[prop]] = ucolors
  #}
  
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
  #annotation_row = annotation_row[rownames(df), drop=FALSE]
  
  colnames(df) = row_names
  #rownames(df) = unlist(xticks_names[[prop]][rownames(df)])
  #rownames(annotation_row) = unlist(xticks_names[[prop]][rownames(annotation_row)])
  
  # print(annotation_row)
  # print(annotation_colors)
  #column_ha = rowAnnotation(df = annotation_row, col = annotation_colors, gp = gpar(row = "white"), simple_anno_size = unit(0.2, "cm"), 
  #                          annotation_name_side = "top", annotation_name_gp = gpar(fontsize = plots_props$font_size, fontface = "bold", fontfamily=plots_props$font_family),
  #                          show_legend=FALSE, show_annotation_name=FALSE)
  #cell_fun = function(j, i, x, y, width, height, fill) {grid.text(df_label[i, j], x, y, gp = gpar(fontsize = 6))}
  
  # Create row annotation
  row_anno = rowAnnotation(
    extra_labels = anno_text(group_size_df, 
                             location = 0, # Aligns text with the heatmap rows
                             just = "left",
                             gp = gpar(fontsize = plots_props$font_size, fontfamily = plots_props$font_family)) # Positions text on the left
  )
  
  rownames(df) = annotation_row
  
  p = ComplexHeatmap::Heatmap(df, col = col_fun,
                              cell_fun = cell_fun,
                              #rect_gp = gpar(col = "white", lwd = .5),
                              #height = unit(0.3, "cm") * nrow(t(df)), width = unit(0.3, "cm") * ncol(t(df)),
                              #height = unit(plots_props$image_height  - 2, "cm"), width = unit(plots_props$image_width - 4, "cm"),
                              height = unit(1, "cm"), width = unit(5.5, "cm"),
                              show_column_names = FALSE,
                              row_names_side = "left", row_names_gp = gpar(fontsize = plots_props$font_size, fontfamily=plots_props$font_family), 
                              column_names_side = "bottom", column_names_gp = gpar(fontsize = 6, fontfamily=plots_props$font_family),
                              cluster_columns = function(m) hclust(dist(m), method = "complete"), show_column_dend = TRUE, clustering_distance_columns = "euclidean", 
                              #cluster_rows = function(m) hclust(dist(m), method = "complete"), show_row_dend = TRUE, clustering_distance_rows = "euclidean", 
                              column_dend_reorder = TRUE, #cluster_columns = FALSE,
                              row_dend_reorder = FALSE, cluster_rows = FALSE,
                              show_heatmap_legend = FALSE,
                              right_annotation = row_anno,
                              #left_annotation = column_ha
  )
  #p
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


# ---------------------- Load data ----------------------
data_input = snakemake@params$data_input
feature = snakemake@params$feature
split_prop = snakemake@params$split_prop
corr_prop = snakemake@params$corr_prop
corr_params = snakemake@params$corr_params
plots_props = snakemake@params$plots_props
xticks_names = snakemake@params$xticks_names
colors = snakemake@params$colors
results_folder = snakemake@params$results_folder

# save rdata
#saveRDS(snakemake, file = "snakemake_correlation_mirna_human_with_age_heatmap_plot.rds")
#stop()

# load annot
annot = fread(sprintf("%s/annotation_%s_%s.csv", data_input$annotation_data_folder, data_input$data_sub_set, data_input$feature_filtering), sep='\t', colClasses=c(ID="character")) 

# load mirna expr
tbl = fread(sprintf("%s/data_%s/%s_%s_quantification_%s.detection_rate_%s_per_%s.csv", data_input$raw_data_folder, data_input$data_sub_set, data_input$rna_class, data_input$feature_filtering, data_input$norm, data_input$detection_rate, data_input$detection_group), sep='\t', header=T)
expr_mirna = tbl[, sapply(tbl, is.numeric), with=F]
names_mirna = tbl[[data_input$rna_class]]
print(sprintf("Shape of miRNA expression matrix: %s", paste(dim(expr_mirna), collapse = "x")))
#sort expression according to annot
fixed_expr_order = annot[[data_input$identifier_column]]
expr_mirna = expr_mirna[, ..fixed_expr_order]

# create output folders
output_folder = sprintf("%s/%s_%s/results_%s/matrices/corr_%s_human_with_%s", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, data_input$rna_class, corr_prop)
dir.create(output_folder, recursive=TRUE)
output_folder_fig = sprintf("%s/%s_%s/results_%s/figures/corr_%s_human_with_%s", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, data_input$rna_class, corr_prop)
dir.create(output_folder_fig, recursive=TRUE)

#------------------------------------ Script ----------------------------------- 
# calculate correlation with age from the mRNA under consideration
for (method in corr_params$method) {
  
  correlation_age = c()
  corr_pvalues_age = c()
  group_names = c()
  group_size = c()
  for(group in split_prop){
    #filter for the investigating group
    if (group$g2[2] == "all") {
      groups_ID = annot[annot[[group$g1[1]]] == group$g1[2],][[data_input$identifier_column]]
    } else if (group$g2[2] == "male") {
      groups_ID = annot[((annot[[group$g1[1]]] == group$g1[2]) & (annot$msex == 1)),][[data_input$identifier_column]]
    } else if (group$g2[2] == "female") {
      groups_ID = annot[((annot[[group$g1[1]]] == group$g1[2]) & (annot$msex == 0)),][[data_input$identifier_column]]
    }
    tmp = colnames(expr_mirna) %in% groups_ID
    expr_sort = data.frame(feature=names_mirna, expr_mirna[,..tmp], check.names = FALSE)
    tmp = annot[[data_input$identifier_column]] %in% groups_ID
    timepoints = annot[tmp,][[corr_prop]]
    
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
      group_name = sprintf("%s=%s_%s=%s", group$g1[1], group$g1[2], group$g2[1], group$g2[2])
      group_size[[group_name]] = length(timepoints)
      correlation_age[[group_name]] = as.data.frame(cor(t(expr_sort[,2:dim(expr_sort)[2]]), numeric_timepoints, method = method))
      # if all expr for every timepoint for a sample is 0 for a feature
      # then the correlation is NA
      # we then force the correlation to be 0
      na_cor = is.na(correlation_age[[group_name]])
      # uncomment this to show the corresponding lines in the expr_matrix
      # print(expr_sort[,2:dim(expr_sort)[2]][na_cor, ])
      # replace NA with 0
      correlation_age[[group_name]][na_cor, ] = 0
      # get pvalues from rcorr
      # $P gets p-values
      pvalues_tmp = as.data.frame(rcorr(t(expr_sort[,2:dim(expr_sort)[2]]), numeric_timepoints, type=method)$P)
      # select y because x=matrix, y=timepoints
      # and remove last value because it is correlation of timepoints with timepoints = NA
      corr_pvalues_age[[group_name]] = head(pvalues_tmp$y, -1)
      # adjustment
      corr_pvalues_age[[group_name]] = as.data.frame(p.adjust(corr_pvalues_age[[group_name]], method=corr_params$adjustment))
      # remove all pvalues == NA and replace them with 1 (then it is not significant)
      na_pvalue = is.na(corr_pvalues_age[[group_name]])
      corr_pvalues_age[[group_name]][na_pvalue] = 1
      # check if anything is still NA
      #print(sprintf("NA in correlation_age: %s", any(is.na(correlation_age[[group_name]]))))
      #print(sprintf("NA in pvalues: %s", any(is.na(corr_pvalues_age[[group_name]]))))
      group_names = append(group_names, group_name)
    }
  }
  group_size_df = as.data.frame(group_size)
  colnames(group_size_df) = names(group_size)
  
  correlation_age_df = as.data.frame(correlation_age)
  colnames(correlation_age_df) = group_names
  correlation_age_df$RNA = expr_sort$feature
  correlation_age_df = correlation_age_df[, rev(colnames(correlation_age_df))]
  
  corr_pvalues_age_df = as.data.frame(corr_pvalues_age)
  colnames(corr_pvalues_age_df) = group_names
  corr_pvalues_age_df$RNA = expr_sort$feature
  corr_pvalues_age_df = corr_pvalues_age_df[, rev(colnames(corr_pvalues_age_df))]
  
  # save corr data files
  output_folder_tab = sprintf("%s/%s", output_folder, method)
  dir.create(output_folder_tab, recursive=TRUE)
  write.table(correlation_age_df, sprintf("%s/correlation_with_%s.csv", output_folder_tab, corr_prop), sep='\t', row.names = F)
  write.table(corr_pvalues_age_df, sprintf("%s/padj_with_%s.csv", output_folder_tab, corr_prop), sep='\t', row.names = F)
  
  write.xlsx(correlation_age_df, sprintf("%s/correlation_with_%s.xlsx", output_folder_tab, corr_prop), colNames = TRUE, rowNames = FALSE, append = FALSE)
  write.xlsx(corr_pvalues_age_df, sprintf("%s/padj_with_%s.xlsx", output_folder_tab, corr_prop), colNames = TRUE, rowNames = FALSE, append = FALSE)
  
  for (corr_th in corr_params$corr_th) {
    
    row_names = correlation_age_df$RNA
    correlation_age_df$RNA = c()
    
    names(group_names) = group_names
    group_names = group_names[colnames(correlation_age_df)]
    annotation_row = names(group_names)
    if (length(split_prop) == 2) {
      if (all(annotation_row %in% list("disease_status=0_samples=female", "disease_status=0_samples=male")))
        annotation_row = ifelse(annotation_row == "disease_status=0_samples=male", "Healthy\nmale", "Healthy\nfemale")
    }
  
    group_size_df_sorted = group_size_df[colnames(correlation_age_df)]
    group_size_df_sorted_rename = sprintf("Number of\nsamples: %s", group_size_df_sorted)
    
    if (any(correlation_age_df$`disease_status=0_samples=male` <= -corr_th)) {
      print(method)
      print(-corr_th)
      print(row_names[correlation_age_df$`disease_status=0_samples=male` <= -corr_th]) 
    } 
    if (any(correlation_age_df$`disease_status=0_samples=male` >= corr_th)) {
      print(method)
      print(corr_th)
      print(row_names[correlation_age_df$`disease_status=0_samples=male` >= corr_th]) 
    }
    if (any(correlation_age_df$`disease_status=0_samples=female` <= -corr_th)) {
      print(method)
      print(-corr_th)
      print(row_names[correlation_age_df$`disease_status=0_samples=female` <= -corr_th]) 
    }
    if (any(correlation_age_df$`disease_status=0_samples=female` >= corr_th)) {
      print(method)
      print(corr_th)
      print(row_names[correlation_age_df$`disease_status=0_samples=female` >= corr_th]) 
    }
    
    heatmap_plot = plot_heatmap_annot(correlation_age_df, corr_pvalues_age_df, row_names, group_size_df_sorted_rename, data_input, xticks_names, plots_props, corr_prop, method, annotation_row, colors, corr_th, corr_params$sig_lvl)
    #heatmap_plot = plot_heatmap_annot_flipped(correlation_age_df, corr_pvalues_age_df, row_names, data_input, xticks_names, plots_props, corr_prop, method, annotation_row, colors, corr_th, corr_params$sig_lvl)
    #df = correlation_age_df
    #pvalues_df = corr_pvalues_age_df
    #time = corr_prop
    #m = method
    #th = corr_th
    #sig_lvl = corr_params$sig_lvl
    
    #row_name_square_plot = row_name_squares(names(unlist(colors[[split_prop]][unique(correlation_age_df_melted$variable)])), xticks_names, split_prop, colors, plots_props)
    
    #row_name_square_plot = row_name_squares_bottom(names(unlist(colors[[split_prop]][unique(correlation_age_df_melted$variable)])), xticks_names, split_prop, colors, plots_props)
    
    #heatmap_plot_comb = row_name_square_plot %v% heatmap_plot$heatmap
    
    # save figures
    height = 3.5
    width = plots_props$image_width
    filename = sprintf("%s/corr_method=%s_%s_over_%s_per_%s", output_folder_fig, method, data_input$rna_class, corr_prop, paste(group_names, collapse = "_"))
    png(sprintf("%s_corr_th=%s_%sx%s.png", filename, corr_th, width, height), width = width, height = height, units = plots_props$image_units, res = plots_props$dpi)
    draw(heatmap_plot$heatmap, gap = unit(0, plots_props$image_units), padding = unit(c(0.25, 1.5, -1, 1.6), "cm"))  # bottom, left, top and right margins
    draw(heatmap_plot$legend, x = unit(1.55, "cm"), y = unit(0, "cm"), just = c("left", "bottom"))
    dev.off()
    svglite(sprintf("%s_corr_th=%s_%sx%s.svg", filename, corr_th, width, height), width = width / 2.54, height = height / 2.54)
    draw(heatmap_plot$heatmap, gap = unit(0, plots_props$image_units), padding = unit(c(0.25, 1.5, -1, 1.75), "cm"))  # bottom, left, top and right margins
    draw(heatmap_plot$legend, x = unit(1.55, "cm"), y = unit(0, "cm"), just = c("left", "bottom"))
    dev.off()
  }
}


#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])
