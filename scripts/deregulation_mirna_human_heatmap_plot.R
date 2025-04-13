suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(svglite))
suppressPackageStartupMessages(library(openxlsx))

# save rdata
#saveRDS(snakemake, file = "snakemake_deregulation_mirna_human_heatmap_plot.rds")
#stop()

#snakemake = readRDS("snakemake_deregulation_mirna_human_heatmap_plot.rds")

set.seed(snakemake@params$parameters_props$set_seed)

print("deregulation_mirna_human_heatmap_plot")


#---------------------------------- Functions ----------------------------------
split_into_three_intervals = function(values, number_of_bins) {
  min_val = min(values)
  max_val = max(values)
  
  # Calculate interval size
  interval_size = (max_val - min_val) / number_of_bins
  
  # Define interval boundaries
  threshold1 = min_val + interval_size
  threshold2 = min_val + 2 * interval_size
  
  return(list(
    min_val = floor(min_val),
    threshold1 = round(threshold1),
    threshold2 = round(threshold2),
    max_val = ceiling(max_val) + 1
  ))
}

plot_heatmap_annot = function(df, pvalues_df, row_names, group_size_df, data_input, xticks_names, plots_props, time, log, annotation_row, colors, th, sig_lvl) {
  
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
  
  #if (m == "spearman") {
  #  color_bar_name = sprintf("Spearman correlation (age)") 
  #} else if (m == "pearson") {
  #  color_bar_name = sprintf("Pearson correlation (age)") 
  #} else {
  #  print("Wrong corr method selected")
  #  stop()
  #}
  
  color_bar_name = sprintf("Fold change (log2)") 
  
  
  #cmap = colorRamp2(c(-1, 0, 1), hcl_palette = "Green-Brown", reverse = TRUE)
  #col_val = 0.75
  #col_fun = c(cmap(-col_val), cmap(0), cmap(col_val))
  #col_legend = c(cmap(-col_val), cmap(col_val))
  #col_fun = c(colors$direction$neg, colors$direction$light, colors$direction$pos)
  color_bar_max = max(max(log_fc_values_df), abs(min(log_fc_values_df)))
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
    else if ((pvalues_df[i, j] < 0.05) & (df[i,j] != 0)){ 
      print(pvalues_df[i, j])
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
  
  legend_ticks = seq(from = -color_bar_max, to = color_bar_max, length.out = 5)
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
property_to_be_binned = snakemake@params$bin_params$property_to_be_binned
number_of_bins = snakemake@params$bin_params$number_of_bins
comps = snakemake@params$comparisons
diff_log = snakemake@params$params$diff_log
fc_th = snakemake@params$params$fc_th
sig_lvl = snakemake@params$params$sig_lvl
plots_props = snakemake@params$plots_props
xticks_names = snakemake@params$xticks_names
colors = snakemake@params$colors
results_folder = snakemake@params$results_folder

# load annot
annot = fread(sprintf("%s/annotation_%s_%s.csv", data_input$annotation_data_folder, data_input$data_sub_set, data_input$feature_filtering), sep='\t', colClasses=c(ID="character")) 

# load mirna expr
#tbl = fread(sprintf("%s_%s/%s_%s_quantification_%s.detection_rate_%s_per_%s.csv", data_input$raw_data_folder, data_input$data_sub_set, data_input$rna_class, data_input$feature_filtering, data_input$norm, data_input$detection_rate, data_input$detection_group), sep='\t', header=T)
#expr_mirna = tbl[, sapply(tbl, is.numeric), with=F]
#names_mirna = tbl[[data_input$rna_class]]
#print(sprintf("Shape of miRNA expression matrix: %s", paste(dim(expr_mirna), collapse = "x")))
#sort expression according to annot
#fixed_expr_order = annot[[data_input$identifier_column]]
#expr_mirna = expr_mirna[, ..fixed_expr_order]

diff_exp = fread(sprintf("%s/%s_%s/results_%s/matrices/diff_exp/diff_exp%s.csv", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, diff_log), sep='\t') 

# create output folders
output_folder_fig = sprintf("%s/%s_%s/results_%s/figures/deregulated_%ss_human_%s", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, data_input$rna_class, names(comps[1]))
dir.create(output_folder_fig, recursive=TRUE)

#------------------------------------ Script ----------------------------------- 
age_intervalls = unlist(unname(split_into_three_intervals(floor(annot[[property_to_be_binned]]), number_of_bins)))

log_fc_values = c()
ttest_adj_pval_values = c()
wilcox_adj_pval_values = c()
cohend_values = c()
group_names = c()
group_size = c()
for(i in 1:length(comps)){
  c = comps[[i]]
  
  # this is only needed if comparisons are loaded from the volcano part of config.yaml
  c$g1 = c$g[1]
  c$g2 = c$g[2]
  
  prop = names(comps)[i]
  
  paired_string_measures = ifelse(!is.null(c$paired), sprintf("_paired_%s_%s_vs_%s", c$paired, c$g1, c$g2), "")
  paired_string = ifelse(!is.null(c$paired), sprintf("_paired_%s", c$paired), "")
  sub_string = ifelse(!is.null(c$subset), sprintf("_%s=%s", c$subset[1], c$subset[2]), "")
  
  group_name = sprintf("%s__%s_vs_%s%s", prop, c$g1, c$g2, sub_string)
  
  plot_df_comp = data.table(feature=diff_exp[[data_input$rna_class]],
                       ttest_raw_pval=diff_exp[[sprintf("ttest_rawp_%s__%s_vs_%s%s%s", prop, c$g1, c$g2, paired_string, sub_string)]],
                       ttest_adj_pval=diff_exp[[sprintf("ttest_adjp_%s__%s_vs_%s%s%s", prop, c$g1, c$g2, paired_string, sub_string)]],
                       wilcox_raw_pval=diff_exp[[sprintf("wilcox_rawp_%s__%s_vs_%s%s%s", prop, c$g1, c$g2, paired_string, sub_string)]],
                       wilcox_adj_pval=diff_exp[[sprintf("wilcox_adjp_%s__%s_vs_%s%s%s", prop, c$g1, c$g2, paired_string, sub_string)]],
                       fc=diff_exp[[sprintf("fc_%s__%s_vs_%s%s%s", prop, c$g1, c$g2, paired_string, sub_string)]],
                       cohend=diff_exp[[sprintf("cohend_estimate_%s__%s_vs_%s%s%s", prop, c$g1, c$g2, paired_string, sub_string)]])
  plot_df_comp[,ttest_raw_pval_log10 := -log10(ttest_raw_pval)]
  plot_df_comp[,ttest_adj_pval_log10 := -log10(ttest_adj_pval)]
  plot_df_comp[,wilcox_raw_pval_log10 := -log10(wilcox_raw_pval)]
  plot_df_comp[,wilcox_adj_pval_log10 := -log10(wilcox_adj_pval)]
  plot_df_comp[,logfc := log2(fc)]
  plot_df_comp[,cohend_abs := abs(cohend)]
  
  log_fc_values[[group_name]] = plot_df_comp$logfc
  ttest_adj_pval_values[[group_name]] = plot_df_comp$ttest_adj_pval
  wilcox_adj_pval_values[[group_name]] = plot_df_comp$wilcox_adj_pval
  cohend_values[[group_name]] = plot_df_comp$cohend
  
  #group_size[[group_name]] = sprintf("(%s samples) vs. (%s)", diff_exp[[sprintf("num_%s__%s%s", prop, c$g1, sub_string)]][1], diff_exp[[sprintf("num_%s__%s%s", prop, c$g2, sub_string)]][1])
  group_size[[group_name]] = sprintf("%s to %s years (%s samples)\nvs. %s to %s (%s)", floor(age_intervalls[3]), floor(age_intervalls[4]-1), diff_exp[[sprintf("num_%s__%s%s", prop, c$g1, sub_string)]][1], floor(age_intervalls[1]), floor(age_intervalls[2]-1), diff_exp[[sprintf("num_%s__%s%s", prop, c$g2, sub_string)]][1])
  group_names = append(group_names, group_name)
}

group_size_df = as.data.frame(group_size)
colnames(group_size_df) = names(group_size)

log_fc_values_df = as.data.frame(log_fc_values)
colnames(log_fc_values_df) = group_names
log_fc_values_df$RNA = diff_exp[[data_input$rna_class]]
log_fc_values_df = log_fc_values_df[, rev(colnames(log_fc_values_df))]

ttest_adj_pval_values_df = as.data.frame(ttest_adj_pval_values)
colnames(ttest_adj_pval_values_df) = group_names
ttest_adj_pval_values_df$RNA = diff_exp[[data_input$rna_class]]
ttest_adj_pval_values_df = ttest_adj_pval_values_df[, rev(colnames(ttest_adj_pval_values_df))]

wilcox_adj_pval_values_df = as.data.frame(wilcox_adj_pval_values)
colnames(wilcox_adj_pval_values_df) = group_names
wilcox_adj_pval_values_df$RNA = diff_exp[[data_input$rna_class]]
wilcox_adj_pval_values_df = wilcox_adj_pval_values_df[, rev(colnames(wilcox_adj_pval_values_df))]

cohend_values_df = as.data.frame(cohend_values)
colnames(cohend_values_df) = group_names
cohend_values_df$RNA = diff_exp[[data_input$rna_class]]
cohend_values_df = cohend_values_df[, rev(colnames(cohend_values_df))]

for (th in log2(fc_th)) {
  row_names = log_fc_values_df$RNA
  log_fc_values_df$RNA = c()

  annotation_row = 
    c("age_bin_3_groups__92_vs_71_msex=1" = "Healthy\nmale",
      "age_bin_3_groups__92_vs_71_msex=0" = "Healthy\nfemale")

  annotation_row = annotation_row[colnames(log_fc_values_df)]
  annotation_row = unname(annotation_row)
  
  group_size_df_sorted = group_size_df[colnames(log_fc_values_df)]
  group_size_df_sorted_rename = unname(group_size_df_sorted)
  
  if (any(log_fc_values_df$`age_bin_3_groups__92_vs_71_msex=0` >= th)) {
    print("female")
    print("up")
    print(2^th)
    print(row_names[log_fc_values_df$`age_bin_3_groups__92_vs_71_msex=0` >= th]) 
    print(length(row_names[log_fc_values_df$`age_bin_3_groups__92_vs_71_msex=0` >= th]))
  } 
  if (any(log_fc_values_df$`age_bin_3_groups__92_vs_71_msex=0` <= -th)) {
    print("female")
    print("down")
    print(-2^th)
    print(row_names[log_fc_values_df$`age_bin_3_groups__92_vs_71_msex=0` <= -th])
    print(length(row_names[log_fc_values_df$`age_bin_3_groups__92_vs_71_msex=0` <= -th])) 
  } 
  if (any(log_fc_values_df$`age_bin_3_groups__92_vs_71_msex=1` >= th)) {
    print("male")
    print("up")
    print(2^th)
    print(row_names[log_fc_values_df$`age_bin_3_groups__92_vs_71_msex=1` >= th]) 
    print(length(row_names[log_fc_values_df$`age_bin_3_groups__92_vs_71_msex=1` >= th]))
  }
  if (any(log_fc_values_df$`age_bin_3_groups__92_vs_71_msex=1` <= -th)) {
    print("male")
    print("down")
    print(-2^th)
    print(row_names[log_fc_values_df$`age_bin_3_groups__92_vs_71_msex=1` <= -th]) 
    print(length(row_names[log_fc_values_df$`age_bin_3_groups__92_vs_71_msex=1` <= -th]))
  }
  
  if (any((log_fc_values_df$`age_bin_3_groups__92_vs_71_msex=0` >= th) & (abs(cohend_values_df$`age_bin_3_groups__92_vs_71_msex=0`) > 0.5))) {
    print("female")
    print("up + cohens")
    print(2^th)
    print(row_names[((log_fc_values_df$`age_bin_3_groups__92_vs_71_msex=0` >= th) & (abs(cohend_values_df$`age_bin_3_groups__92_vs_71_msex=0`) > 0.5))]) 
    print(length(row_names[((log_fc_values_df$`age_bin_3_groups__92_vs_71_msex=0` >= th) & (abs(cohend_values_df$`age_bin_3_groups__92_vs_71_msex=0`) > 0.5))]))
  } 
  if (any((log_fc_values_df$`age_bin_3_groups__92_vs_71_msex=0` <= -th) & (abs(cohend_values_df$`age_bin_3_groups__92_vs_71_msex=0`) > 0.5))) {
    print("female")
    print("down + cohens")
    print(-2^th)
    print(row_names[((log_fc_values_df$`age_bin_3_groups__92_vs_71_msex=0` <= -th) & (abs(cohend_values_df$`age_bin_3_groups__92_vs_71_msex=0`) > 0.5))])
    print(length(row_names[((log_fc_values_df$`age_bin_3_groups__92_vs_71_msex=0` <= -th) & (abs(cohend_values_df$`age_bin_3_groups__92_vs_71_msex=0`) > 0.5))])) 
  } 
  if (any((log_fc_values_df$`age_bin_3_groups__92_vs_71_msex=1` >= th) & (abs(cohend_values_df$`age_bin_3_groups__92_vs_71_msex=1`) > 0.5))) {
    print("male")
    print("up + cohens")
    print(2^th)
    print(row_names[((log_fc_values_df$`age_bin_3_groups__92_vs_71_msex=1` >= th) & (abs(cohend_values_df$`age_bin_3_groups__92_vs_71_msex=1`) > 0.5))]) 
    print(length(row_names[((log_fc_values_df$`age_bin_3_groups__92_vs_71_msex=1` >= th) & (abs(cohend_values_df$`age_bin_3_groups__92_vs_71_msex=1`) > 0.5))]))
  }
  if (any((log_fc_values_df$`age_bin_3_groups__92_vs_71_msex=1` <= -th) & (abs(cohend_values_df$`age_bin_3_groups__92_vs_71_msex=1`) > 0.5))) {
    print("male")
    print("down + cohens")
    print(-2^th)
    print(row_names[((log_fc_values_df$`age_bin_3_groups__92_vs_71_msex=1` <= -th) & (abs(cohend_values_df$`age_bin_3_groups__92_vs_71_msex=1`) > 0.5))]) 
    print(length(row_names[((log_fc_values_df$`age_bin_3_groups__92_vs_71_msex=1` <= -th) & (abs(cohend_values_df$`age_bin_3_groups__92_vs_71_msex=1`) > 0.5))]))
  }
  
  ttest_adj_pval_values_df$RNA = c()
  
  heatmap_plot = plot_heatmap_annot(log_fc_values_df, ttest_adj_pval_values_df, row_names, group_size_df_sorted_rename, data_input, xticks_names, plots_props, corr_prop, diff_log, annotation_row, colors, th, corr_params$sig_lvl)
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
  width = 12
  filename = sprintf("%s/deregulation_fc_th=%s_per_%s", output_folder_fig, (2^th), paste(group_names, collapse = "_"))
  png(sprintf("%s_%sx%s.png", filename, width, height), width = width, height = height, units = plots_props$image_units, res = plots_props$dpi)
  draw(heatmap_plot$heatmap, gap = unit(0, plots_props$image_units), padding = unit(c(0.25, 1.5, -1, 1.6), "cm"))  # bottom, left, top and right margins
  draw(heatmap_plot$legend, x = unit(1.55, "cm"), y = unit(0, "cm"), just = c("left", "bottom"))
  dev.off()
  svglite(sprintf("%s_%sx%s.svg", filename, width, height), width = width / 2.54, height = height / 2.54)
  draw(heatmap_plot$heatmap, gap = unit(0, plots_props$image_units), padding = unit(c(0.25, 1.5, -1, 1.75), "cm"))  # bottom, left, top and right margins
  draw(heatmap_plot$legend, x = unit(1.55, "cm"), y = unit(0, "cm"), just = c("left", "bottom"))
  dev.off()
}



#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])
