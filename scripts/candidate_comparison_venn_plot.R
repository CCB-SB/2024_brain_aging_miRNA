suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(viridisLite))
suppressPackageStartupMessages(library(svglite))
suppressPackageStartupMessages(library(ggnetwork))
suppressPackageStartupMessages(library(ggVennDiagram))
suppressPackageStartupMessages(library(eulerr))
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(grid))
suppressMessages(library(stringi))
suppressPackageStartupMessages(library(ComplexHeatmap))


#snakemake = readRDS("snakemake_candidate_comparison_venn_plot.rds")

data.table::setDTthreads(1)
set.seed(snakemake@params$parameters_props$set_seed)

print("candidate_comparison_venn_plot")


#---------------------------------- Functions ----------------------------------
replace_null_with_empty_vector = function(input) {
  # because of some strange behaviour, the venn function thinks c() != NULL != character(0)
  # it seems that it needs the last
  if (length(input) == 0) {
    return(character(0))
  } else {
    return(input)
  }
}

reduce_landscape_spacing = function() {
  return(theme(plot.margin=unit(c(t=-0.15, r=-0.25, b=-0.15, l=-0.25), "cm")))
}

# write the tissue name on the left side of a venn
tissue_name_side = function(plot, tissue_group, tissue, xticks_names) {
  plot_with_label = plot + annotation_custom(
    grob = textGrob(sprintf("%s", xticks_names[[tissue_group]][tissue]), rot = 90, vjust = 0.5, hjust = 0.5),
    xmin = 50, xmax = -Inf, ymin = -Inf, ymax = Inf
    )
  return(plot_with_label)
}

# write the tissue name on the top side of a venn
tissue_name_top = function(plot, tissue_group, tissue, xticks_names) {
  plot_with_label = plot + annotation_custom(
    grob = textGrob(sprintf("%s", xticks_names[[tissue_group]][tissue]), rot = 0, vjust = 0.5, hjust = 0.5),
    xmin = -Inf, xmax = Inf, ymin = 800, ymax = Inf
  )
  return(plot_with_label)
}

venn_plot = function(diff_exp_ts_rnas, corr_ts_rnas, hclust_ts_rnas, maximal_value, plot_props, file_name, file_name_tab, show_set_name = TRUE) {
  
  font_size_factor = 2.857
  
  plot_df = list(diff_exp=diff_exp_ts_rnas, corr=corr_ts_rnas, hclust=hclust_ts_rnas)
  # manually calculate the stuff ggVennDiagramm should do automatically
  # this way we can color the edges manually
  venn = Venn(plot_df)
  features_cluster_manual_venn = process_data(venn)

  items = features_cluster_manual_venn@region$item
  item_combined = c()
  for (item in items) {
    item_combined = append(item_combined, do.call(paste, c(list(item), collapse = "\n")))
  }
  # remove mmu-miR
  item_combined = gsub("mmu-miR-", "", item_combined)
  features_cluster_manual_venn@region$item_combined = item_combined

  
  cluster_colours_list = c("diff_exp"="grey", "corr"="red", "hclust"="blue")
  # use colorRampPalette to create function that interpolates colors 
  #colfunc <- colorRampPalette(cluster_colours_list)
  # call function and create vector of 15 colors
  #col <- colfunc(15)
  
  # cmap = colorRamp2(c(0, maximal_value), hcl_palette = "YlOrRd")
  cmap = colorRamp2(c(0, maximal_value), colors = c("#ffffb2", "#bd0026"))
  # find maximal color for this specific plot
  max_value_here = (max(venn_setedge(features_cluster_manual_venn)$count))
  min_colour = cmap(0)
  max_colour = cmap(max_value_here)
  
  if (show_set_name) {
    set_size = plot_props$font_size_header /3
  } else {
    set_size = 0
  }
  p_labels = ggVennDiagram(plot_df,
                           set_size = set_size,
                           label="none", label_color = "black", label_size = 6/font_size_factor, label_alpha=0,
                           edge_size = 0.5,
                           ) +
    # text with boxes
    # geom_sf_label(aes(label = item_combined), data = venn_region(features_cluster_manual_venn), alpha=0.5) + 
    # no boxes
    geom_sf_text(aes(label = item_combined), data = venn_region(features_cluster_manual_venn), size=plot_props$font_size /3) + 
    geom_sf(aes(color = name), data = venn_setedge(features_cluster_manual_venn), show.legend = F) +
    #scale_color_manual(values = cluster_colours_list) +
    # expand the x axis
    scale_x_continuous(expand = expansion(mult = 0.2)) +
    scale_y_continuous(expand = expansion(mult = 0.3)) +
    #scale_fill_distiller(palette = "YlOrRd", direction = -1, values=c(0, 5, 10)) +
    # geom_sf(aes(fill = name), data = venn_region(features_cluster_manual_venn), show.legend = F) +
    # scale_fill_manual(values =  alpha(cluster_colours_list, .2)) +
    #scale_fill_manual(values = c("grey","grey", "grey", "grey")) +
    # scale_fill_gradient(low = "#D5FBFF", high = "#0092BF") +
    #scale_fill_gradient(low = "#F7D748", high = "#d62728") +
    #scale_fill_gradient(low = "#ffffb2", high = "#bd0026") +
    scale_fill_gradient(low = min_colour, high = max_colour) + # pre-calculated
    scale_color_manual(values = rep("black", 2*length(plot_df))) + 
    #scale_fill_viridis_c() + 
    #scale_y_continuous(expand = expansion(mult = .1)) +
    theme(text = element_text(size = plot_props$font_size, family = plot_props$font_family),
          plot.title = element_text(family = plot_props$font_family, size = plot_props$font_size_header),
          legend.text = element_text(family = plot_props$font_family, size = plot_props$font_size),
          legend.title = element_text(family = plot_props$font_family, size = plot_props$font_size, face = "bold"),
          legend.position = "none",# legend.box = "horizontal",
          plot.margin = unit(c(-0.25,-0.25,-0.25,-0.25), "cm"),
          #panel.spacing = margin(t = 0, r = 0, b = -3, l = 0, unit = "cm"),
    ) #+ 
    #guides(fill = guide_colourbar(title.position = "top", title.hjust = 0), size = guide_legend(title.position = "top", title.hjust = 0))
  p_labels
  
  ggsave(sprintf("%s_labels.png", file_name), p_labels, dpi=plot_props$dpi, width=plot_props$image_width, height=plot_props$image_height, units = plot_props$image_units)
  ggsave(sprintf("%s_labels.svg", file_name), p_labels, width=plot_props$image_width, height=plot_props$image_height, units =plot_props$image_units, limitsize = FALSE, scale = 1)

  if (show_set_name) {
    set_size = plot_props$font_size_header /3
  } else {
    set_size = 0
  }
  p = ggVennDiagram(plot_df, 
                    set_size = set_size,
                    label="count", label_color = "black", label_size = 6/font_size_factor, label_alpha=0,
                    edge_size = 0.5,
  ) +
    # text with boxes
    # geom_sf_label(aes(label = item_combined), data = venn_region(features_cluster_manual_venn), alpha=0.5) + 
    # no boxes
    #geom_sf_text(aes(label = item_combined), data = venn_region(features_cluster_manual_venn), size=6 * 25.4/72) + 
    geom_sf(aes(color = name), data = venn_setedge(features_cluster_manual_venn), show.legend = F) +
    #scale_color_manual(values = cluster_colours_list) +
    # expand the x axis
    scale_x_continuous(expand = expansion(mult = 0.2)) +
    scale_y_continuous(expand = expansion(mult = 0.3)) +
    # geom_sf(aes(fill = name), data = venn_region(features_cluster_manual_venn), show.legend = F) +
    # scale_fill_manual(values =  alpha(cluster_colours_list, .2)) +
    # scale_fill_manual(values = c("grey","grey","grey","grey")) +
    # scale_fill_gradient(low = "#D5FBFF", high = "#0092BF") +
    #scale_fill_gradient(low = "#F7D748", high = "#d62728") +
    # scale_fill_gradient(low = "#ffffb2", high = "#bd0026") +
    scale_fill_gradient(low = min_colour, high = max_colour) + # pre-calculated
    scale_color_manual(values = rep("black", 2*length(plot_df))) + 
    #scale_fill_gradient(low = "white", high = "white", guide = "none") + 
    #scale_fill_viridis_c() + 
    theme(text = element_text(size = plot_props$font_size, family = plot_props$font_family),
          plot.title = element_text(family = plot_props$font_family, size = plot_props$font_size_header),
          legend.text = element_text(family = plot_props$font_family, size = plot_props$font_size),
          legend.title = element_text(family = plot_props$font_family, size = plot_props$font_size, face = "bold"),
          legend.position = "none",# legend.box = "horizontal",
          plot.margin = unit(c(-0.25,-0.25,-0.25,-0.25), "cm"),
          #panel.spacing = margin(t = 0, r = 0, b = -3, l = 0, unit = "cm"),
    ) #+ 
  #guides(fill = guide_colourbar(title.position = "top", title.hjust = 0), size = guide_legend(title.position = "top", title.hjust = 0))
  #p

  ggsave(sprintf("%s.png", file_name), p, dpi=plot_props$dpi, width=plot_props$image_width, height=plot_props$image_height, units = plot_props$image_units)
  ggsave(sprintf("%s.svg", file_name), p, width=plot_props$image_width, height=plot_props$image_height, units =plot_props$image_units, limitsize = FALSE, scale = 1)
  
  # pad both lists
  max_length = max(length(plot_df$diff_exp), length(plot_df$corr))
  if (max_length == 0) {
    # make a least one Element if both are empty
    max_length = 1
  }
  list_diff_exp_padded = c(plot_df$diff_exp, rep(NA, max_length - length(plot_df$diff_exp)))
  list_corr_padded = c(plot_df$corr, rep(NA, max_length - length(plot_df$corr)))
  
  plot_df_padded = data.frame(diff_exp=list_diff_exp_padded, corr=list_corr_padded)

  fwrite(plot_df_padded, sprintf("%s.csv", file_name_tab), sep = "\t", row.names = FALSE)
  write.xlsx(plot_df_padded, sprintf("%s.xlsx", file_name_tab), colNames = TRUE, rowNames = FALSE, append = FALSE)
  
  return(list(p_labels = p_labels, p = p))
}

upset_plot = function(feature_candidates_tissues, tissue, xticks_names, colours, plot_props, file_name) {

  # upset plot
  mode = "distinct"
  intersection_th_list = c(1,2,5,10)
  
  for (th in intersection_th_list) {
  
    plot_table = make_comb_mat(feature_candidates_tissues, mode=mode)
    plot_table_thresholded = plot_table[comb_size(plot_table) >= th]
    
    cs = comb_size(plot_table_thresholded)
    ss = set_size(plot_table_thresholded)
    nc = ncol(plot_table_thresholded)
    
    # sort color vector according to plot_df
    colour_list = unlist(colours[[tissue]])[rownames(plot_table_thresholded)]
    rownames(plot_table_thresholded) =xticks_names[[tissue]][rownames(plot_table_thresholded)]
    
    upset_plot = ComplexHeatmap::UpSet(plot_table_thresholded, pt_size = unit(1, "mm"), lwd = 1, comb_order = order(comb_size(plot_table_thresholded), decreasing = TRUE), #pt_size = unit(1, "mm"), lwd = 1, bg_col = c("#ECECEC", "#CBCBCB")
                                       height = unit(3.5, "cm"), width = unit(4.5, "cm"),
                                       row_names_side = "left", row_names_gp = gpar(fontsize = plot_props$font_size, fontfamily=plot_props$font_family), 
                                       top_annotation = HeatmapAnnotation("miRNA\noverlap" = anno_barplot(comb_size(plot_table_thresholded),
                                                                                                          #add_numbers=TRUE,
                                                                                                          ylim = c(0, max(comb_size(plot_table_thresholded))*1.1), 
                                                                                                          border = FALSE, 
                                                                                                          gp = gpar(fill = "black", fontsize = plot_props$font_size, fontfamily=plot_props$font_family),
                                                                                                          labels_gp = gpar(col = "black", fontsize = 6), 
                                                                                                          height = unit(1, "cm")
                                       ), 
                                       annotation_name_side = "left", annotation_name_rot = 90, annotation_name_gp = gpar(fontsize = plot_props$font_size, fontfamily=plot_props$font_family)
                                       ),
                                       
                                       right_annotation = rowAnnotation("Num. of candidates\nper brain regions" = anno_barplot(set_size(plot_table_thresholded), 
                                                                                                                                 #add_numbers=TRUE,
                                                                                                                                 border = FALSE, 
                                                                                                                                 gp = gpar(fill = colour_list, col = "white", fontsize = plot_props$font_size, fontfamily=plot_props$font_family),
                                                                                                                                 labels_gp = gpar(col = "black", fontsize = 6), 
                                                                                                                                 width = unit(2.5, "cm"),
                                                                                                                                 axis_param = list(labels_rot = 0),
                                       ),
                                       annotation_name_side = "bottom", annotation_name_gp = gpar(fontsize = plot_props$font_size, fontfamily=plot_props$font_family)
                                       ),
                                       #column_title = sprintf("%s", plot_title)
                                       )
    co = column_order(upset_plot)
    ro = row_order(upset_plot)
    png(sprintf("%s_intersection_th=%s.png", file_name, th), width = plot_props$image_width, height = plot_props$image_height, units = plot_props$image_units, res = plot_props$dpi)
    upset_plot = draw(upset_plot, padding = unit(c(0.25, -0.25, 0.6, 0.25), "cm"));  # bottom, left, top and right margins
    decorate_annotation("miRNA\noverlap", {grid.text(cs[co], x = 1:nc, y = unit(cs[co], "native") + unit(1, "mm"), gp = gpar(fontsize = 6), rot = 90, vjust = 0.35, hjust = 0.15, just = "bottom", default.units = "native")})
    decorate_annotation("Num. of candidates\nper brain regions", {grid.text(ss[rev(ro)], x = unit(ss[rev(ro)], "native") + unit(1, "mm"), y = 1:length(ro), gp = gpar(fontsize = 6), vjust = 0.35, hjust = 0.15, just = "left", default.units = "native")})
    dev.off()
    svglite(sprintf("%s_intersection_th=%s.svg", file_name, th), width = plot_props$image_width / 2.54, height = plot_props$image_height / 2.54)
    upset_plot = draw(upset_plot, padding = unit(c(0.25, -0.25, 0.6, 0.25), "cm"));  # bottom, left, top and right margins
    decorate_annotation("miRNA\noverlap", {grid.text(cs[co], x = 1:nc, y = unit(cs[co], "native") + unit(1, "mm"), gp = gpar(fontsize = 6), rot = 90, vjust = 0.35, hjust = 0.15, just = "bottom", default.units = "native")})
    decorate_annotation("Num. of candidates\nper brain regions", {grid.text(ss[rev(ro)], x = unit(ss[rev(ro)], "native") + unit(1, "mm"), y = 1:length(ro), gp = gpar(fontsize = 6), vjust = 0.35, hjust = 0.15, just = "left", default.units = "native")})
    dev.off()
  }
}


#---------------------------------- Load data ----------------------------------
data_input = snakemake@params$data_input
diff_exp_props = snakemake@params$diff_exp_props
corr_props = snakemake@params$corr_props
hclust_props = snakemake@params$hclust_props
props = snakemake@params$props
tissue_specific_fig_size= snakemake@params$tissue_specific_fig_size
all_features_fig_size = snakemake@params$all_features_fig_size
xticks_names = snakemake@params$xticks_names
colours = snakemake@params$colors
plot_props = snakemake@params$plots_props
results_folder = snakemake@params$results_folder

# save rdata
#saveRDS(snakemake, file = "snakemake_candidate_comparison_venn_plot.rds")

output_folder_venn = sprintf("%s/%s_%s/results_%s/figures/candidate_comparison/venn_plots", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set)
dir.create(output_folder_venn, recursive=TRUE)
output_folder_venn_tab = sprintf("%s/%s_%s/results_%s/matrices/candidate_comparison/venn_plots", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set)
dir.create(output_folder_venn_tab, recursive=TRUE)
output_folder_upset = sprintf("%s/%s_%s/results_%s/figures/candidate_comparison/upset_plots", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set)
dir.create(output_folder_upset, recursive=TRUE)


#------------------------------------ Script ----------------------------------- 
for (prop in props) {
  
  tissue = prop[1]
  time = prop[2]
  for (num_dereg_combinations in diff_exp_props$num_dereg_combinations_th) {
    for (num_dereg_tissues in diff_exp_props$num_dereg_tissues_th) {
      for (m in corr_props$method) {
        for (th in corr_props$corr_th) {
          for (num_corr_tissues in corr_props$num_corr_tissues_th) {
            for (m_hclust in hclust_props$method) {
              # only first part (for example cv from cv_zscore)
              m_hclust_split = data.frame(stri_split_fixed(m_hclust, "_"))[1,]
              for (top in hclust_props$top) {
                for (th_hclust in hclust_props$zscores_th) {
    
                  # diff exp
                  file_to_load_up = sprintf("%s/%s_%s/results_%s/matrices/aggregation_diff_exp/diff_exp_sig_up_%s_at_least_in_%s_%s_comps_for_%s_%s.csv", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, data_input$rna_class, num_dereg_combinations, time, num_dereg_tissues, tissue)
                  file_to_load_down = sprintf("%s/%s_%s/results_%s/matrices/aggregation_diff_exp/diff_exp_sig_down_%s_at_least_in_%s_%s_comps_for_%s_%s.csv", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, data_input$rna_class, num_dereg_combinations, time, num_dereg_tissues, tissue)
                  
                  diff_exp_list_up = fread(file_to_load_up)
                  # remove columns containing only 0 
                  tmp = colSums(diff_exp_list_up != 0) > 0
                  diff_exp_list_up = diff_exp_list_up[, ..tmp]
                  diff_exp_list_down = fread(file_to_load_down)
                  tmp = colSums(diff_exp_list_down != 0) > 0
                  diff_exp_list_down = diff_exp_list_down[, ..tmp]

                  diff_exp_list_up_rownames = diff_exp_list_up$V1
                  diff_exp_list_up$V1 = c()
                  tmp = (rowSums(diff_exp_list_up) == 1)
                  diff_exp_list_up_ts = diff_exp_list_up[tmp,]
                  diff_exp_list_up_rownames_ts = diff_exp_list_up_rownames[tmp]
                  tmp = colSums(diff_exp_list_up_ts != 0) > 0
                  diff_exp_list_up_ts = diff_exp_list_up_ts[, ..tmp]
                  
                  diff_exp_list_down_rownames = diff_exp_list_down$V1
                  diff_exp_list_down$V1 = c()
                  tmp = (rowSums(diff_exp_list_down) == 1)
                  diff_exp_list_down_ts = diff_exp_list_down[tmp,]
                  diff_exp_list_down_rownames_ts = diff_exp_list_down_rownames[tmp]
                  tmp = colSums(diff_exp_list_down_ts != 0) > 0
                  diff_exp_list_down_ts = diff_exp_list_down_ts[, ..tmp]
                  
                  
                  # corr
                  file_to_load_pos = sprintf("%s/%s_%s/results_%s/matrices/aggregation_corr/corr_method=%s_sig_pos_for_%s_%s_%s_corr_th=%s.csv", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, m, num_corr_tissues, tissue, data_input$rna_class, th)
                  file_to_load_neg = sprintf("%s/%s_%s/results_%s/matrices/aggregation_corr/corr_method=%s_sig_neg_for_%s_%s_%s_corr_th=%s.csv", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, m, num_corr_tissues, tissue, data_input$rna_class, th)
                  
                  corr_list_pos = fread(file_to_load_pos)
                  tmp = colSums(corr_list_pos != 0) > 0
                  corr_list_pos = corr_list_pos[, ..tmp]
                  corr_list_neg = fread(file_to_load_neg)
                  tmp = colSums(corr_list_neg != 0) > 0
                  corr_list_neg = corr_list_neg[, ..tmp]
                  
                  corr_list_pos_rownames = corr_list_pos$V1
                  corr_list_pos$V1 = c()
                  tmp = (rowSums(corr_list_pos) == 1)
                  corr_list_pos_ts = corr_list_pos[tmp,]
                  corr_list_pos_rownames_ts = corr_list_pos_rownames[tmp]
                  tmp = colSums(corr_list_pos_ts != 0) > 0
                  corr_list_pos_ts = corr_list_pos_ts[, ..tmp]
                  
                  corr_list_neg_rownames = corr_list_neg$V1
                  corr_list_neg$V1 = c()
                  tmp = (rowSums(corr_list_neg) == 1)
                  corr_list_neg_ts = corr_list_neg[tmp,]
                  corr_list_neg_rownames_ts = corr_list_neg_rownames[tmp]
                  tmp = colSums(corr_list_neg_ts != 0) > 0
                  corr_list_neg_ts = corr_list_neg_ts[, ..tmp]
                  
                  
                  # hclust
                  file_to_load = sprintf("%s/%s_%s/results_%s/matrices/%s_clustering/%s_zscore_th=%s_plot_df.csv", results_folder, data_input$rna_class, data_input$detection_rate, hclust_props$subset, m_hclust_split, top, th_hclust)

                  hclust_table = fread(file_to_load)
                  tmp = colSums(hclust_table != 0) > 0
                  hclust_table = hclust_table[, ..tmp]
                  
                  hclust_table_rownames = hclust_table[[data_input$rna_class]]
                  hclust_table[[data_input$rna_class]] = c()
                  tmp = (rowSums(hclust_table) == 1)
                  hclust_table_ts = hclust_table[tmp]
                  hclust_table_rownames_ts = hclust_table_rownames[tmp]
                  tmp = colSums(hclust_table_ts != 0) > 0
                  hclust_table_ts = hclust_table_ts[, ..tmp]
                  
                  
                  
                  
                  # venn plot
                  # only if not empty
                  
                  # merge
                  if (dim(diff_exp_list_up_ts)[1] == 0) {
                    diff_exp_list_up_ts = c()
                  }
                  if (dim(diff_exp_list_down_ts)[1] == 0) {
                    diff_exp_list_down_ts = c()
                  }
                  if (dim(corr_list_pos_ts)[1] == 0) {
                    corr_list_pos_ts = c()
                  } 
                  if (dim(corr_list_neg_ts)[1] == 0) {
                    corr_list_neg_ts = c()
                  }
                  if (dim(hclust_table_ts)[1] == 0) {
                    hclust_table_ts = c()
                  }
                  
                  maximal_value_ts = 0
                  for (sub_list in list(diff_exp_list_up_ts, diff_exp_list_down_ts, corr_list_pos_ts, corr_list_neg_ts, hclust_table_ts)) {
                    # only go into if clause if it is really a matrix
                    if (!is.null(nrow(sub_list)) && !is.null(ncol(sub_list))) {
                      maximal_value_ts = max(maximal_value_ts, colSums(sub_list))
                    }
                  }
                  # does not work if one is NULL
                  #maximal_value_ts = max(colSums(diff_exp_list_up_ts), colSums(diff_exp_list_down_ts), colSums(corr_list_pos_ts), colSums(corr_list_neg_ts), colSums(hclust_table_ts))
                  
                  tissue_list_up_pos_ts = union(colnames(diff_exp_list_up_ts), colnames(corr_list_pos_ts))
                  tissue_list_down_neg_ts = union(colnames(diff_exp_list_down_ts), colnames(corr_list_neg_ts))
                  tissue_list_ts_merged = unique(c(tissue_list_up_pos_ts, tissue_list_down_neg_ts, colnames(hclust_table_ts)))
                  
                  # tissue specific
                  if (length(tissue_list_ts_merged) > 0) {
                    p_both_labels = list()
                    p_both = list()
                    p_both_landscape_labels = list()
                    p_both_landscape = list()
                    i = 1
                    feature_candidates_tissues_ts = c()
                    feature_candidates_tissues_ts_without_hclust = c()
                    for (ti in tissue_list_ts_merged) {
                      
                      diff_exp_list_up_ts_rnas = diff_exp_list_up_rownames_ts[diff_exp_list_up_ts[[ti]] == 1] 
                      corr_list_pos_ts_rnas = corr_list_pos_rownames_ts[corr_list_pos_ts[[ti]] == 1]
                      
                      diff_exp_list_down_ts_rnas = diff_exp_list_down_rownames_ts[diff_exp_list_down_ts[[ti]] == 1] 
                      corr_list_neg_ts_rnas = corr_list_neg_rownames_ts[corr_list_neg_ts[[ti]] == 1]
                      
                      hclust_table_ts_rnas = hclust_table_rownames_ts[hclust_table_ts[[ti]] == 1]
                      
                      if (length(intersect(diff_exp_list_up_ts_rnas, corr_list_neg_ts_rnas)) > 0) {
                        print(sprintf("!Same %s in sig. up regulated and sig negative correlated in tissue %s for given parameters!", data_input$rna_class, ti))
                        opt = options(show.error.messages = FALSE)
                        on.exit(options(opt))
                        stop()
                      } else if (length(intersect(diff_exp_list_down_ts_rnas, corr_list_pos_ts_rnas)) > 0) {
                        print(sprintf("!Same %s in sig. down regulated and sig postive correlated in tissue %s for given parameters!", data_input$rna_class, ti))
                        opt = options(show.error.messages = FALSE)
                        on.exit(options(opt))
                        stop()
                      }
                      
                      diff_exp_list_ts_rnas = unique(c(diff_exp_list_up_ts_rnas, diff_exp_list_down_ts_rnas))
                      corr_list_ts_rnas = unique(c(corr_list_pos_ts_rnas, corr_list_neg_ts_rnas))
                      
                      # remove NULL
                      diff_exp_list_ts_rnas = replace_null_with_empty_vector(diff_exp_list_ts_rnas)
                      corr_list_ts_rnas = replace_null_with_empty_vector(corr_list_ts_rnas)
                      hclust_table_ts_rnas = replace_null_with_empty_vector(hclust_table_ts_rnas)
                      
                      # if all three sets are empty
                      if ((length(diff_exp_list_ts_rnas) == 0) && (length(corr_list_ts_rnas) == 0) && (length(hclust_table_ts_rnas) == 0)) {
                        print(sprintf("no sig corr, dregulated or hclust noticeable %ss", data_input$rna_class))
                        next
                      }
                      
                      file_name = sprintf("candidate_comp_%s_specific_%s_sig_diff_exp_both_%s_at_least_in_%s_%s_comps_for_%s_%s_sig_corr_method=%s_both_for_%s_%s_%s_corr_th=%s_hclust_method=%s_th=%s_top=%s_%ss", tissue, ti, data_input$rna_class, num_dereg_combinations, time, num_dereg_tissues, tissue, m, num_corr_tissues, tissue, data_input$rna_class, th, m_hclust, th_hclust, top, data_input$rna_class)
                      file_path_venn_plot = sprintf("%s/%s", output_folder_venn, file_name)
                      file_path_venn_tab = sprintf("%s/%s", output_folder_venn_tab, file_name)
        
                      # first, we run the plot function to fill without set name for the variable p_output
                      # then, we plot it again with show_set_name = TRUE to save the right plot
                      # and fill the variable p_output_with_set_name
                      p_output = venn_plot(diff_exp_list_ts_rnas, corr_list_ts_rnas, hclust_table_ts_rnas, maximal_value_ts, plot_props, file_path_venn_plot, file_path_venn_tab, show_set_name = FALSE)
                      p_output_with_set_name = venn_plot(diff_exp_list_ts_rnas, corr_list_ts_rnas, hclust_table_ts_rnas, maximal_value_ts, plot_props, file_path_venn_plot, file_path_venn_tab, show_set_name = TRUE)
        
                      #p_both_labels[[i]] = p_output$p_labels + theme(plot.title = element_text(angle = 90, vjust = 0.5, hjust = 0.5)) + labs(title = sprintf("%s", ti))
                      #p_both[[i]] = p_output$p + theme(plot.title = element_text(angle = 90, vjust = 0.5, hjust = 0.5)) + labs(title = sprintf("%s", ti))
                      
                      # write the tissue names on the left side (vertical)
                      #if (i == 1) {
                      #  # first plot with set names
                      #  p_both_labels[[i]] = tissue_name_side(p_output_with_set_name$p_labels, ti)
                      #  p_both[[i]] = tissue_name_side(p_output_with_set_name$p, ti)
                      #} else {
                      # the rest without
                      p_both_labels[[i]] = tissue_name_side(p_output$p_labels, tissue, ti, xticks_names)
                      p_both[[i]] = tissue_name_side(p_output$p, tissue, ti, xticks_names)
                      #}
                      # on top (landscape)
                      # never print set name
                      p_both_landscape_labels[[i]] = tissue_name_top(p_output$p_labels, tissue, ti, xticks_names) + reduce_landscape_spacing()
                      p_both_landscape[[i]] = tissue_name_top(p_output$p, tissue, ti, xticks_names) + reduce_landscape_spacing()
                      
                      # euler plot
                      #file_name_euler = sprintf("%s/sig_tissue_specific_comparison_diff_exp_%s_%s_at_least_in_%s_%s_comps_for_%s_%s_corr_method=%s_%s_for_%s_%s_%s_corr_th=%s", output_folder_euler, diff_exp_dir, data_input$rna_class, num_dereg_combinations, time, num_dereg_tissues, tissue, m, corr_dir, num_corr_tissues, tissue, data_input$rna_class, th)
                      #euler_plot(diff_exp_ts_rnas, corr_ts_rnas, plot_props, file_name_euler)
                      
                      feature_candidates_tissues_ts[[ti]] = unique(c(diff_exp_list_ts_rnas, corr_list_ts_rnas, hclust_table_ts_rnas))
                      feature_candidates_tissues_ts_without_hclust[[ti]] = unique(c(diff_exp_list_ts_rnas, corr_list_ts_rnas))
                      
                      i = i + 1
                    }
                    
                    # merged plots
                    number_of_plots = length(tissue_list_ts_merged)
                    file_name_venn_merged = sprintf("candidate_comp_%s_specific_sig_diff_exp_both_%s_at_least_in_%s_%s_comps_for_%s_%s_sig_corr_method=%s_both_for_%s_%s_%s_corr_th=%s_hclust_method=%s_th=%s_top=%s_%ss", tissue, data_input$rna_class, num_dereg_combinations, time, num_dereg_tissues, tissue, m, num_corr_tissues, tissue, data_input$rna_class, th, m_hclust, th_hclust, top, data_input$rna_class)
                    file_path_venn_merged = sprintf("%s/%s", output_folder_venn, file_name_venn_merged)
                    
                    # save plots
                    # vertical
                    # number_of_plots x 1
                    all_p_labels = arrangeGrob(grobs=p_both_labels, ncol=1, nrow=number_of_plots, gp=gpar(fontsize=plot_props$font_size, fontfamily=plot_props$font_family))
                    all_p = arrangeGrob(grobs=p_both, ncol=1, nrow=number_of_plots, gp=gpar(fontsize=plot_props$font_size, fontfamily=plot_props$font_family))
                    
                    ggsave(sprintf("%s_vertical_vector_num=%s_labels.png", file_path_venn_merged, number_of_plots), all_p_labels, dpi=plot_props$dpi, width=0.25*plot_props$image_width, height=0.25*number_of_plots*plot_props$image_height, units = plot_props$image_units)
                    ggsave(sprintf("%s_vertical_vector_num=%s_labels.svg", file_path_venn_merged, number_of_plots), all_p_labels, width=0.25*plot_props$image_width, height=0.25*number_of_plots*plot_props$image_height, units =plot_props$image_units, limitsize = FALSE, scale = 1)
                    
                    ggsave(sprintf("%s_vertical_vector_num=%s.png", file_path_venn_merged, number_of_plots), all_p, dpi=plot_props$dpi, width=0.25*plot_props$image_width, height=0.25*number_of_plots*plot_props$image_height, units = plot_props$image_units)
                    ggsave(sprintf("%s_vertical_vector_num=%s.svg", file_path_venn_merged, number_of_plots), all_p, width=0.25*plot_props$image_width, height=0.25*number_of_plots*plot_props$image_height, units =plot_props$image_units, limitsize = FALSE, scale = 1)
                    
                    # horizontal
                    # 1 x number_of_plots
                    all_p_labels = arrangeGrob(grobs=p_both_landscape_labels, ncol=number_of_plots, nrow=1, gp=gpar(fontsize=plot_props$font_size, fontfamily=plot_props$font_family))
                    all_p = arrangeGrob(grobs=p_both_landscape, ncol=number_of_plots, nrow=1, gp=gpar(fontsize=plot_props$font_size, fontfamily=plot_props$font_family))
                    
                    ggsave(sprintf("%s_landscape_vector_num=%s_labels.png", file_path_venn_merged, number_of_plots), all_p_labels, dpi=plot_props$dpi, width=0.25*number_of_plots*plot_props$image_width, height=0.5*plot_props$image_height, units = plot_props$image_units)
                    ggsave(sprintf("%s_landscape_vector_num=%s_labels.svg", file_path_venn_merged, number_of_plots), all_p_labels, width=0.25*number_of_plots*plot_props$image_width, height=0.5*plot_props$image_height, units =plot_props$image_units, limitsize = FALSE, scale = 1)
                    
                    ggsave(sprintf("%s_landscape_vector_num=%s.png", file_path_venn_merged, number_of_plots), all_p, dpi=plot_props$dpi, width=0.25*number_of_plots*plot_props$image_width, height=0.5*plot_props$image_height, units = plot_props$image_units)
                    ggsave(sprintf("%s_landscape_vector_num=%s.svg", file_path_venn_merged, number_of_plots), all_p, width=0.25*number_of_plots*plot_props$image_width, height=0.5*plot_props$image_height, units =plot_props$image_units, limitsize = FALSE, scale = 1)
                    
                    # 3 x 5
                    all_p_labels = arrangeGrob(grobs=p_both_landscape_labels, ncol=tissue_specific_fig_size$plot_ncols, nrow=tissue_specific_fig_size$plot_nrows, gp=gpar(fontsize=plot_props$font_size, fontfamily=plot_props$font_family))
                    all_p = arrangeGrob(grobs=p_both_landscape, ncol=tissue_specific_fig_size$plot_ncols, nrow=tissue_specific_fig_size$plot_nrows, gp=gpar(fontsize=plot_props$font_size, fontfamily=plot_props$font_family))
                    
                    ggsave(sprintf("%s_landscape_num=%s_labels.png", file_path_venn_merged, number_of_plots), all_p_labels, dpi=plot_props$dpi, width=plot_props$image_width, height=plot_props$image_height, units = plot_props$image_units)
                    ggsave(sprintf("%s_landscape_num=%s_labels.svg", file_path_venn_merged, number_of_plots), all_p_labels, width=plot_props$image_width, height=plot_props$image_height, units =plot_props$image_units, limitsize = FALSE, scale = 1)
                    
                    ggsave(sprintf("%s_landscape_num=%s.png", file_path_venn_merged, number_of_plots), all_p, dpi=plot_props$dpi, width=plot_props$image_width, height=plot_props$image_height, units = plot_props$image_units)
                    ggsave(sprintf("%s_landscape_num=%s.svg", file_path_venn_merged, number_of_plots), all_p, width=plot_props$image_width, height=plot_props$image_height, units =plot_props$image_units, limitsize = FALSE, scale = 1)
                    
                    # upset plot
                    file_name_upset = sprintf("candidate_comp_%s_specific_sig_diff_exp_both_%s_at_least_in_%s_%s_comps_for_%s_%s_sig_corr_method=%s_both_for_%s_%s_%s_corr_th=%s_hclust_method=%s_th=%s_top=%s_%ss", tissue, data_input$rna_class, num_dereg_combinations, time, num_dereg_tissues, tissue, m, num_corr_tissues, tissue, data_input$rna_class, th, m_hclust, th_hclust, top, data_input$rna_class)
                    file_name_upset = sprintf("%s/%s", output_folder_upset, file_name_upset)
                    upset_plot(feature_candidates_tissues_ts, tissue, xticks_names, colours, plot_props, file_name_upset)
                    
                    file_name_upset = sprintf("candidate_comp_%s_specific_sig_diff_exp_both_%s_at_least_in_%s_%s_comps_for_%s_%s_sig_corr_method=%s_both_for_%s_%s_%s_corr_th=%s", tissue, data_input$rna_class, num_dereg_combinations, time, num_dereg_tissues, tissue, m, num_corr_tissues, tissue, data_input$rna_class, th)
                    file_name_upset = sprintf("%s/%s", output_folder_upset, file_name_upset)
                    upset_plot(feature_candidates_tissues_ts_without_hclust, tissue, xticks_names, colours, plot_props, file_name_upset)
                  }
                
                  
                  # all
                  # merge          
                  if (dim(diff_exp_list_up)[1] == 0) {
                    diff_exp_list_up = c()
                  }
                  if (dim(diff_exp_list_down)[1] == 0) {
                    diff_exp_list_down = c()
                  }
                  if (dim(corr_list_pos)[1] == 0) {
                    corr_list_pos = c()
                  } 
                  if (dim(corr_list_neg)[1] == 0) {
                    corr_list_neg = c()
                  }
                  if (dim(hclust_table)[1] == 0) {
                    hclust_table = c()
                  }
                  
                  maximal_value = max(colSums(diff_exp_list_up), colSums(diff_exp_list_down), colSums(corr_list_pos), colSums(corr_list_neg), colSums(hclust_table))
                  
                  tissue_list_up_pos = union(colnames(diff_exp_list_up), colnames(corr_list_pos))
                  tissue_list_down_neg = union(colnames(diff_exp_list_down), colnames(corr_list_neg))
                  tissue_list_merged = unique(c(tissue_list_up_pos, tissue_list_down_neg, colnames(hclust_table)))
                  
                  if (length(tissue_list_merged) > 0) {
                    p_both_labels = list()
                    p_both = list()
                    p_both_landscape_labels = list()
                    p_both_landscape = list()
                    i = 1
                    feature_candidates_tissues = c()
                    feature_candidates_tissues_without_hclust = c()
                    for (ti in tissue_list_merged) {
                      diff_exp_list_up_rnas = diff_exp_list_up_rownames[diff_exp_list_up[[ti]] == 1] 
                      corr_list_pos_rnas = corr_list_pos_rownames[corr_list_pos[[ti]] == 1]
                      
                      diff_exp_list_down_rnas = diff_exp_list_down_rownames[diff_exp_list_down[[ti]] == 1] 
                      corr_list_neg_rnas = corr_list_neg_rownames[corr_list_neg[[ti]] == 1]
                      
                      hclust_table_rnas = hclust_table_rownames[hclust_table[[ti]] == 1]
                      
                      if (length(intersect(diff_exp_list_up_rnas, corr_list_neg_rnas)) > 0) {
                        print(sprintf("Same %s in sig. up regulated and sig negative correlated in tissue %s for given parameters!", data_input$rna_class, ti))
                        opt = options(show.error.messages = FALSE)
                        on.exit(options(opt))
                        stop()
                      } else if (length(intersect(diff_exp_list_down_rnas, corr_list_pos_rnas)) > 0) {
                        print(sprintf("Same %s in sig. down regulated and sig postive correlated in tissue %s for given parameters!", data_input$rna_class, ti))
                        opt = options(show.error.messages = FALSE)
                        on.exit(options(opt))
                        stop()
                      }
                      
                      diff_exp_list_rnas = unique(c(diff_exp_list_up_rnas, diff_exp_list_down_rnas))
                      corr_list_rnas = unique(c(corr_list_pos_rnas, corr_list_neg_rnas))
                      
                      # remove NULL
                      diff_exp_list_rnas = replace_null_with_empty_vector(diff_exp_list_rnas)
                      corr_list_rnas = replace_null_with_empty_vector(corr_list_rnas)
                      hclust_table_rnas = replace_null_with_empty_vector(hclust_table_rnas)
                      
                      # if all three sets are empty
                      if ((length(diff_exp_list_rnas) == 0) && (length(corr_list_rnas) == 0) && (length(hclust_table_rnas) == 0)) {
                        print(sprintf("no sig corr, dregulated or hclust noticeable %ss", data_input$rna_class))
                        next
                      }
                      
                      file_name = sprintf("candidate_comp_%s_sig_diff_exp_both_%s_at_least_in_%s_%s_comps_for_%s_%s_sig_corr_method=%s_both_for_%s_%s_%s_corr_th=%s_hclust_method=%s_th=%s_top=%s_%ss", ti, data_input$rna_class, num_dereg_combinations, time, num_dereg_tissues, tissue, m, num_corr_tissues, tissue, data_input$rna_class, th, m_hclust, th_hclust, top, data_input$rna_class)
                      file_path_venn_plot = sprintf("%s/%s", output_folder_venn, file_name)
                      file_path_venn_tab = sprintf("%s/%s", output_folder_venn_tab, file_name)
                      
                      # first, we run the plot function to fill without set name for the variable p_output
                      # then, we plot it again with show_set_name = TRUE to save the right plot
                      # and fill the variable p_output_with_set_name
                      p_output = venn_plot(diff_exp_list_rnas, corr_list_rnas, hclust_table_rnas, maximal_value, plot_props, file_path_venn_plot, file_path_venn_tab, show_set_name = FALSE)
                      p_output_with_set_name = venn_plot(diff_exp_list_rnas, corr_list_rnas, hclust_table_rnas, maximal_value, plot_props, file_path_venn_plot, file_path_venn_tab, show_set_name = TRUE)
                      
                      #p_both_labels[[i]] = p_output$p_labels + theme(plot.title = element_text(angle = 90, vjust = 0.5, hjust = 0.5)) + labs(title = sprintf("%s", ti))
                      #p_both[[i]] = p_output$p + theme(plot.title = element_text(angle = 90, vjust = 0.5, hjust = 0.5)) + labs(title = sprintf("%s", ti))
                      
                      # write the tissue names on the left side (vertical)
                      #if (i == 1) {
                      #  # first plot with set names
                      #  p_both_labels[[i]] = tissue_name_side(p_output_with_set_name$p_labels, ti)
                      #  p_both[[i]] = tissue_name_side(p_output_with_set_name$p, ti)
                      #} else {
                      # the rest without
                      p_both_labels[[i]] = tissue_name_side(p_output$p_labels, tissue, ti, xticks_names)
                      p_both[[i]] = tissue_name_side(p_output$p, tissue, ti, xticks_names)
                      #}
                      # on top (landscape)
                      # never print set name
                      p_both_landscape_labels[[i]] = tissue_name_top(p_output$p_labels, tissue, ti, xticks_names) + reduce_landscape_spacing()
                      p_both_landscape[[i]] = tissue_name_top(p_output$p, tissue, ti, xticks_names) + reduce_landscape_spacing()
                      
                      # euler plot
                      #file_name_euler = sprintf("%s/sig_tissue_specific_comparison_diff_exp_%s_%s_at_least_in_%s_%s_comps_for_%s_%s_corr_method=%s_%s_for_%s_%s_%s_corr_th=%s", output_folder_euler, diff_exp_dir, data_input$rna_class, num_dereg_combinations, time, num_dereg_tissues, tissue, m, corr_dir, num_corr_tissues, tissue, data_input$rna_class, th)
                      #euler_plot(diff_exp_ts_rnas, corr_ts_rnas, plot_props, file_name_euler)
                      
                      feature_candidates_tissues[[ti]] = unique(c(diff_exp_list_rnas, corr_list_rnas, hclust_table_rnas))
                      feature_candidates_tissues_without_hclust[[ti]] = unique(c(diff_exp_list_rnas, corr_list_rnas))
                      
                      i = i + 1
                    }
                    
                    # merged plots
                    number_of_plots = length(tissue_list_merged)
                    file_name_venn_merged = sprintf("candidate_comp_sig_diff_exp_both_%s_at_least_in_%s_%s_comps_for_%s_%s_sig_corr_method=%s_both_for_%s_%s_%s_corr_th=%s_hclust_method=%s_th=%s_top=%s_%ss", data_input$rna_class, num_dereg_combinations, time, num_dereg_tissues, tissue, m, num_corr_tissues, tissue, data_input$rna_class, th, m_hclust, th_hclust, top, data_input$rna_class)
                    file_path_venn_merged = sprintf("%s/%s", output_folder_venn, file_name_venn_merged)
                    
                    # save plots
                    # vertical
                    # number_of_plots x 1
                    all_p_labels = arrangeGrob(grobs=p_both_labels, ncol=1, nrow=number_of_plots, gp=gpar(fontsize=plot_props$font_size, fontfamily=plot_props$font_family))
                    all_p = arrangeGrob(grobs=p_both, ncol=1, nrow=number_of_plots, gp=gpar(fontsize=plot_props$font_size, fontfamily=plot_props$font_family))
                     
                    ggsave(sprintf("%s_vertical_vector_num=%s_labels.png", file_path_venn_merged, number_of_plots), all_p_labels, dpi=plot_props$dpi, width=0.25*plot_props$image_width, height=0.25*number_of_plots*plot_props$image_height, units = plot_props$image_units)
                    ggsave(sprintf("%s_vertical_vector_num=%s_labels.svg", file_path_venn_merged, number_of_plots), all_p_labels, width=0.25*plot_props$image_width, height=0.25*number_of_plots*plot_props$image_height, units =plot_props$image_units, limitsize = FALSE, scale = 1)
                    
                    ggsave(sprintf("%s_vertical_vector_num=%s.png", file_path_venn_merged, number_of_plots), all_p, dpi=plot_props$dpi, width=0.25*plot_props$image_width, height=0.25*number_of_plots*plot_props$image_height, units = plot_props$image_units)
                    ggsave(sprintf("%s_vertical_vector_num=%s.svg", file_path_venn_merged, number_of_plots), all_p, width=0.25*plot_props$image_width, height=0.25*number_of_plots*plot_props$image_height, units =plot_props$image_units, limitsize = FALSE, scale = 1)
                    
                    # horizontal
                    # 1 x number_of_plots
                    all_p_labels = arrangeGrob(grobs=p_both_landscape_labels, ncol=number_of_plots, nrow=1, gp=gpar(fontsize=plot_props$font_size, fontfamily=plot_props$font_family))
                    all_p = arrangeGrob(grobs=p_both_landscape, ncol=number_of_plots, nrow=1, gp=gpar(fontsize=plot_props$font_size, fontfamily=plot_props$font_family))
                    
                    ggsave(sprintf("%s_landscape_vector_num=%s_labels.png", file_path_venn_merged, number_of_plots), all_p_labels, dpi=plot_props$dpi, width=0.25*number_of_plots*plot_props$image_width, height=0.5*plot_props$image_height, units = plot_props$image_units)
                    ggsave(sprintf("%s_landscape_vector_num=%s_labels.svg", file_path_venn_merged, number_of_plots), all_p_labels, width=0.25*number_of_plots*plot_props$image_width, height=0.5*plot_props$image_height, units =plot_props$image_units, limitsize = FALSE, scale = 1)
                    
                    ggsave(sprintf("%s_landscape_vector_num=%s.png", file_path_venn_merged, number_of_plots), all_p, dpi=plot_props$dpi, width=0.25*number_of_plots*plot_props$image_width, height=0.5*plot_props$image_height, units = plot_props$image_units)
                    ggsave(sprintf("%s_landscape_vector_num=%s.svg", file_path_venn_merged, number_of_plots), all_p, width=0.25*number_of_plots*plot_props$image_width, height=0.5*plot_props$image_height, units =plot_props$image_units, limitsize = FALSE, scale = 1)
                    
                    # 3 x 5
                    all_p_labels = arrangeGrob(grobs=p_both_landscape_labels, ncol=all_features_fig_size$plot_ncols, nrow=all_features_fig_size$plot_nrows, gp=gpar(fontsize=plot_props$font_size, fontfamily=plot_props$font_family))
                    all_p = arrangeGrob(grobs=p_both_landscape, ncol=all_features_fig_size$plot_ncols, nrow=all_features_fig_size$plot_nrows, gp=gpar(fontsize=plot_props$font_size, fontfamily=plot_props$font_family))
                    
                    ggsave(sprintf("%s_landscape_num=%s_labels.png", file_path_venn_merged, number_of_plots), all_p_labels, dpi=plot_props$dpi, width=plot_props$image_width, height=plot_props$image_height, units = plot_props$image_units)
                    ggsave(sprintf("%s_landscape_num=%s_labels.svg", file_path_venn_merged, number_of_plots), all_p_labels, width=plot_props$image_width, height=plot_props$image_height, units =plot_props$image_units, limitsize = FALSE, scale = 1)
                    
                    ggsave(sprintf("%s_landscape_num=%s.png", file_path_venn_merged, number_of_plots), all_p, dpi=plot_props$dpi, width=plot_props$image_width, height=plot_props$image_height, units = plot_props$image_units)
                    ggsave(sprintf("%s_landscape_num=%s.svg", file_path_venn_merged, number_of_plots), all_p, width=plot_props$image_width, height=plot_props$image_height, units =plot_props$image_units, limitsize = FALSE, scale = 1)
                    
                    
                    # upset plot
                    file_name_upset = sprintf("candidate_comp_sig_diff_exp_both_%s_at_least_in_%s_%s_comps_for_%s_%s_sig_corr_method=%s_both_for_%s_%s_%s_corr_th=%s_hclust_method=%s_th=%s_top=%s_%ss", data_input$rna_class, num_dereg_combinations, time, num_dereg_tissues, tissue, m, num_corr_tissues, tissue, data_input$rna_class, th, m_hclust, th_hclust, top, data_input$rna_class)
                    file_name_upset = sprintf("%s/%s", output_folder_upset, file_name_upset)
                    upset_plot(feature_candidates_tissues, tissue, xticks_names, colours, plot_props, file_name_upset)
                    
                    file_name_upset = sprintf("candidate_comp_sig_diff_exp_both_%s_at_least_in_%s_%s_comps_for_%s_%s_sig_corr_method=%s_both_for_%s_%s_%s_corr_th=%s", data_input$rna_class, num_dereg_combinations, time, num_dereg_tissues, tissue, m, num_corr_tissues, tissue, data_input$rna_class, th)
                    file_name_upset = sprintf("%s/%s", output_folder_upset, file_name_upset)
                    upset_plot(feature_candidates_tissues_without_hclust, tissue, xticks_names, colours, plot_props, file_name_upset)
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}


#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])
