suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(svglite))
suppressPackageStartupMessages(library(ggVennDiagram))
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(openxlsx))

# save rdata
#saveRDS(snakemake, file = "snakemale_files/snakemake_correlation_mirna_target_genes_cas_upset_plot.rds")
#stop()

#snakemake = readRDS("snakemale_files/snakemake_correlation_mirna_target_genes_cas_upset_plot.rds")

set.seed(snakemake@params$parameters_props$set_seed)

print("correlation_mirna_target_genes_cas_upset_plot")


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

row_name_squares = function(tissue_list, object_height, object_width, xticks_names, tissue, colors, plots_props) {
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
    grid.text(wording, x = x, y = y, gp = gpar(col = text_colour, fontsize = plots_props$font_size)) #fontfamily=plots_props$font_family
  }
  
  #col_fun  = structure(unname(unlist(colors[[tissue]])), names = rownames(data))
  
  # Create the heatmap
  p = Heatmap(data,
              #col = col_fun,
              height = unit(object_height, "cm"), width = unit(object_width, "cm"),
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

row_name_squares_adj = function(tissue_list, object_height, object_width, xticks_names, genes, colors, plots_props) {
  # Create a sample matrix
  data =- matrix(1:length(tissue_list), nrow = length(tissue_list), ncol = 1)
  rownames(data) = tissue_list
  colnames(data) = "Square"
  
  #print(setdiff(tissue_list, gene_colours$genes))
  #print(colors)
  
  gene_colours_list = colors$colour
  names(gene_colours_list) = colors[[genes]]
  gene_colours_text_list = colors$text_colour
  names(gene_colours_text_list) = colors[[genes]]
  
  # Create a function to draw text inside squares
  draw_text = function(j, i, x, y, width, height, fill) {
    wording = rownames(data)[i]
    text_colour = gene_colours_text_list[[rownames(data)[i]]]
    rect_colour = gene_colours_list[[rownames(data)[i]]]
    #print(sprintf("draw color: %s %s", wording, rect_colour))
    grid.roundrect(x=x, y=y, width=width, height=height, r = unit(0.1, "cm"), gp = gpar(fill = rect_colour, col = "white", lwd = 0.5))
    grid.text(wording, x = x, y = y, gp = gpar(col = text_colour, fontsize = plots_props$font_size)) #fontfamily=plots_props$font_family
  }
  
  #col_fun  = structure(unname(unlist(colors[[tissue]])), names = rownames(data))
  
  # Create the heatmap
  p = Heatmap(data,
              #col = col_fun,
              height = unit(object_height, "cm"), width = unit(object_width, "cm"),
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

col_name_squares = function(tissue_list, object_height, object_width, xticks_names, tissue, colors, plots_props) {
  # Create a sample matrix
  data =- matrix(1:length(tissue_list), nrow = 1, ncol = length(tissue_list))
  rownames(data) = "Square"
  colnames(data) = tissue_list
  
  # Create a function to draw text inside squares
  draw_text = function(j, i, x, y, width, height, fill) { 
    wording = xticks_names[[tissue]][[colnames(data)[j]]]
    text_colour = unlist(colors[[sprintf("%s_text", tissue)]])[[colnames(data)[j]]]
    rect_colour = unlist(colors[[tissue]])[[colnames(data)[j]]]
    grid.roundrect(x=x, y=y, width=width, height=height, r = unit(0.1, "cm"), gp = gpar(fill = rect_colour, col = "white", lwd = 2))
    grid.text(wording, x = x, y = y, gp = gpar(col = text_colour, fontsize = plots_props$font_size)) #fontfamily=plots_props$font_family
  }
  
  #col_fun  = structure(unname(unlist(colors[[tissue]])), names = rownames(data))
  
  # Create the heatmap
  p = Heatmap(data,
              #col = col_fun,
              height = unit(object_height, "cm"), width = unit(object_width, "cm"),
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

# create_venn_plot = function(age_specific, tissue_specific, maximal_value, ti, colours, plot_props, file_name, file_name_tab, show_set_name = TRUE) {
#   
#   font_size_factor = 2.857
#   
#   plot_df = list(Age=age_specific, Region=tissue_specific)
#   # manually calculate the stuff ggVennDiagramm should do automatically
#   # this way we can color the edges manually
#   venn = Venn(plot_df)
#   features_tissue_manual_venn = process_data(venn)
# 
#   items = features_tissue_manual_venn@region$item
#   item_combined = c()
#   for (item in items) {
#     item_combined = append(item_combined, do.call(paste, c(list(item), collapse = "\n")))
#   }
#   # remove mmu-miR
#   item_combined = gsub("mmu-miR-", "", item_combined)
#   features_tissue_manual_venn@region$item_combined = item_combined
# 
#   
#   set_colours_list = c("Age"="#333990", "Region"="#E09900")
#   tissue_colours_list = unlist(colours[[tissue]])
#   tissue_colour = tissue_colours_list[[ti]]
#   # use colorRampPalette to create function that interpolates colors 
#   #colfunc <- colorRampPalette(tissue_colours_list)
#   # call function and create vector of 15 colors
#   #col <- colfunc(15)
#   
#   # cmap = colorRamp2(c(0, maximal_value), hcl_palette = "YlOrRd")
#   cmap = colorRamp2(c(0, maximal_value), colors = c("#ffffb2", "#bd0026"))
#   # find maximal color for this specific plot
#   max_value_here = (max(venn_setedge(features_tissue_manual_venn)$count))
#   min_colour = cmap(0)
#   max_colour = cmap(max_value_here)
#   
#   #if (show_set_name) {
#   #  set_size = plot_props$font_size
#   #} else {
#   set_size = 0
#   #}
#   p_labels = ggVennDiagram(plot_df,
#                            set_size = set_size/font_size_factor, set_color = rep(tissue_colour, length(plot_df)),
#                            label="none", label_color = "black", label_size = plot_props$font_size/font_size_factor, label_alpha=0,
#                            edge_size = 0.5,
#                            ) +
#     # text with boxes
#     # geom_sf_label(aes(label = item_combined), data = venn_region(features_tissue_manual_venn), alpha=0.5) + 
#     # set names
#     geom_sf_text(aes(label = name), data = venn_setlabel(features_tissue_manual_venn), nudge_y = 80, size=plot_props$font_size * 25.4/72, colour = tissue_colour) +
#     # no boxes
#     geom_sf_text(aes(label = item_combined), data = venn_region(features_tissue_manual_venn), size=plot_props$font_size /3) + 
#     geom_sf(aes(color = name), data = venn_setedge(features_tissue_manual_venn), show.legend = F) +
#     # expand the x axis
#     scale_x_continuous(expand = expansion(mult = 0.2)) +
#     scale_y_continuous(expand = expansion(mult = c(0.1, 0.2))) +
#     #scale_fill_distiller(palette = "YlOrRd", direction = -1, values=c(0, 5, 10)) +
#     # geom_sf(aes(fill = name), data = venn_region(features_tissue_manual_venn), show.legend = F) +
#     # scale_fill_manual(values =  alpha(tissue_colours_list, .2)) +
#     #scale_fill_manual(values = c("grey","grey", "grey", "grey")) +
#     # scale_fill_gradient(low = "#D5FBFF", high = "#0092BF") +
#     #scale_fill_gradient(low = "#F7D748", high = "#d62728") +
#     #scale_fill_gradient(low = "#ffffb2", high = "#bd0026") +
#     scale_fill_gradient(low = min_colour, high = max_colour) + # pre-calculated
#     scale_color_manual(values = rep(tissue_colour, 2*length(plot_df))) + 
#     #scale_fill_viridis_c() + 
#     #scale_y_continuous(expand = expansion(mult = .1)) +
#     theme(text = element_text(size = plot_props$font_size, family = plot_props$font_family),
#           plot.title = element_text(family = plot_props$font_family, size = plot_props$font_size_header, color = tissue_colour),
#           legend.text = element_text(family = plot_props$font_family, size = plot_props$font_size),
#           legend.title = element_text(family = plot_props$font_family, size = plot_props$font_size, face = "bold"),
#           legend.position = "none",# legend.box = "horizontal",
#           plot.margin = unit(c(-1,0,-1,0), "cm"),
#           #plot.margin = unit(c(-0.25,-0.25,-0.25,-0.25), "cm"),
#           #panel.spacing = margin(t = 0, r = 0, b = -3, l = 0, unit = "cm"),
#     ) #+ 
#     #guides(fill = guide_colourbar(title.position = "top", title.hjust = 0), size = guide_legend(title.position = "top", title.hjust = 0))
#   #p_labels
#   
#   #ggsave(sprintf("%s_labels.png", file_name), p_labels, dpi=plot_props$dpi, width=plot_props$image_width, height=plot_props$image_height, units = plot_props$image_units)
#   #ggsave(sprintf("%s_labels.svg", file_name), p_labels, width=plot_props$image_width, height=plot_props$image_height, units =plot_props$image_units, limitsize = FALSE, scale = 1)
# 
#   #if (show_set_name) {
#   #  set_size = plot_props$font_size 
#   #} else {
#   set_size = 0
#   #}
#   p = ggVennDiagram(plot_df, 
#                     set_size = set_size/font_size_factor, set_color = rep(tissue_colour, length(plot_df)),
#                     label="count", label_color = "black", label_size = plot_props$font_size/font_size_factor, label_alpha=0,
#                     edge_size = 0.5,
#   ) +
#     # text with boxes
#     # geom_sf_label(aes(label = item_combined), data = venn_region(features_tissue_manual_venn), alpha=0.5) + 
#     # set names
#     geom_sf_text(aes(label = name), data = venn_setlabel(features_tissue_manual_venn), nudge_y = 80, size=plot_props$font_size * 25.4/72, colour = tissue_colour) +
#     # no boxes
#     #geom_sf_text(aes(label = item_combined), data = venn_region(features_tissue_manual_venn), size=6 * 25.4/72) + 
#     geom_sf(aes(color = name), data = venn_setedge(features_tissue_manual_venn), show.legend = F) +
#     # expand the x axis
#     scale_x_continuous(expand = expansion(mult = 0.2)) +
#     scale_y_continuous(expand = expansion(mult = c(0.1, 0.2))) +
#     # geom_sf(aes(fill = name), data = venn_region(features_tissue_manual_venn), show.legend = F) +
#     # scale_fill_manual(values =  alpha(tissue_colours_list, .2)) +
#     # scale_fill_manual(values = c("grey","grey","grey","grey")) +
#     # scale_fill_gradient(low = "#D5FBFF", high = "#0092BF") +
#     #scale_fill_gradient(low = "#F7D748", high = "#d62728") +
#     # scale_fill_gradient(low = "#ffffb2", high = "#bd0026") +
#     scale_fill_gradient(low = min_colour, high = max_colour) + # pre-calculated
#     scale_color_manual(values = rep(tissue_colour, 2*length(plot_df))) + 
#     #scale_fill_gradient(low = "white", high = "white", guide = "none") + 
#     #scale_fill_viridis_c() + 
#     theme(text = element_text(size = plot_props$font_size, family = plot_props$font_family),
#           plot.title = element_text(family = plot_props$font_family, size = plot_props$font_size_header, color = tissue_colour),
#           legend.text = element_text(family = plot_props$font_family, size = plot_props$font_size),
#           legend.title = element_text(family = plot_props$font_family, size = plot_props$font_size, face = "bold"),
#           legend.position = "none",# legend.box = "horizontal",
#           plot.margin = unit(c(-1,0,-1,0), "cm"),
#           #plot.margin = unit(c(-0.25,-0.25,-0.25,-0.25), "cm"),
#           #panel.spacing = margin(t = 0, r = 0, b = -3, l = 0, unit = "cm"),
#     ) #+ 
#   #guides(fill = guide_colourbar(title.position = "top", title.hjust = 0), size = guide_legend(title.position = "top", title.hjust = 0))
#   #p
# 
#   #ggsave(sprintf("%s.png", file_name), p, dpi=plot_props$dpi, width=plot_props$image_width, height=plot_props$image_height, units = plot_props$image_units)
#   #ggsave(sprintf("%s.svg", file_name), p, width=plot_props$image_width, height=plot_props$image_height, units =plot_props$image_units, limitsize = FALSE, scale = 1)
#   
#   # pad both lists
#   max_length = max(length(plot_df$Age), length(plot_df$Region))
#   if (max_length == 0) {
#     # make a least one Element if both are empty
#     max_length = 1
#   }
# 
#   list_age_padded = c(plot_df$Age, rep(NA, max_length - length(plot_df$Age)))
#   list_tissue_padded = c(plot_df$Region, rep(NA, max_length - length(plot_df$Region)))
#   plot_df_padded = data.frame(Age=list_age_padded, Region=list_tissue_padded)
# 
#   fwrite(plot_df_padded, sprintf("%s.csv", file_name_tab), sep = "\t", row.names = FALSE)
#   write.xlsx(plot_df_padded, sprintf("%s.xlsx", file_name_tab), colNames = TRUE, rowNames = FALSE, append = FALSE)
#   
#   return(list(p_labels = p_labels, p = p))
# }

create_upset_plot_rot = function(sig_neg_correlated_genes_tissue, split_prop, xticks_names, colours, plot_props, file_name_upset, file_name_tab) {
  # save mat
  write.table(list_to_matrix(sig_neg_correlated_genes_tissue), sprintf("%s.csv", file_name_tab), sep='\t', row.names = T)
  write.xlsx(list_to_matrix(sig_neg_correlated_genes_tissue), sprintf("%s.xlsx", file_name_tab), colNames = TRUE, rowNames = TRUE, append = FALSE)
  
  # upset plot
  mode = "distinct"
  intersection_th_list = c(1,2,5,10)
  
  for (th in intersection_th_list) {
    plot_table = make_comb_mat(sig_neg_correlated_genes_tissue, mode=mode)
    #print(dim(plot_table))
    # skip if all combinattions are < th
    if (all(comb_size(plot_table) < th)) {
      #print(sprintf("skipping intersection_th=%s", th))
      next
    }
    plot_table_thresholded = plot_table[, comb_size(plot_table) >= th]
    #print(dim(plot_table_thresholded))
    
    #if (th == 5) {
    #  interesting_regions = c("pon", "plx", "svz", "olf", "vis", "med", "th", "cor", "cc")
    #  tmp = rownames(plot_table_thresholded)[rownames(plot_table_thresholded) %in% interesting_regions]
    #  plot_table_thresholded = plot_table_thresholded[tmp, ]
    #  # now filter the comb_size again!
    #  plot_table_thresholded = plot_table_thresholded[, comb_size(plot_table_thresholded) >= th]
    #  
    #  print(rownames(plot_table_thresholded))
    #}
    #print(dim(plot_table_thresholded))
    #print(rownames(plot_table_thresholded))
    
    if (length(set_size(plot_table_thresholded)) >= 15) {
      gene_font_size = 6
    } else {
      gene_font_size = plot_props$font_size
    }
    
    # sort color vector according to plot_df
    #colour_list = unlist(colours[[split_prop]])[rownames(plot_table_thresholded)]
    
    object_height = 3
    object_width = 6
    
    plot_table_thresholded = t(plot_table_thresholded)
    
    
    # find which col belongs to which brain region only if the size of all combinations is exactly 1, then change the order of rows for the plot
    if ((length(unique(comb_size(plot_table_thresholded))) == 1) && (unique(comb_size(plot_table_thresholded)) == 1)) {
      # transform combination list into matrix
      mat = do.call(rbind, type.convert(strsplit(comb_name(plot_table_thresholded, readable = FALSE), ""), as.is = TRUE))
      # set colnames according to upset plot object. order is by implementation eqaul
      colnames(mat) = set_name(plot_table_thresholded)
      
      # sort according to image
      tmp = names(sort(colSums(mat), decreasing = TRUE))
      mat = mat[,tmp]
      
      mat_orig = list_to_matrix(sig_neg_correlated_genes_tissue)
      mat_orig = mat_orig[, colnames(mat)]
      
      # make list of row vectors
      target = split(mat, row(mat))
      orig = split(mat_orig, row(mat_orig))
      
      # find the ordering
      target_order = match(target, orig)
      #mat_orig[target_order, ]
      
      sorted_rownames = rownames(mat_orig[target_order, ])
      
      # sort plot rows and rownames alphabetically as given in the colours list
      ordered_rownames = names(colours[[split_prop]])[names(colours[[split_prop]]) %in% sorted_rownames]
      
      # extract the keys (number for reordering)
      target_row_order = match(ordered_rownames, sorted_rownames)
    } else {
      target_row_order = order(comb_size(plot_table_thresholded), decreasing = TRUE)
    }
    
    upset_plot = ComplexHeatmap::UpSet(plot_table_thresholded, pt_size = unit(1, "mm"), comb_order = target_row_order, lwd = 1, #pt_size = unit(1, "mm"), lwd = 1, bg_col = c("#ECECEC", "#CBCBCB")
                                       height = unit(object_height, "cm"), width = unit(object_width, "cm"),
                                       column_names_side = "bottom", column_names_gp = gpar(fontsize = gene_font_size, fontfamily=plot_props$font_family),
                                       #top_annotation = NULL,
                                       top_annotation = HeatmapAnnotation("Num. of brain reg.\n(with sig. neg. corr.)" = anno_barplot(set_size(plot_table_thresholded), 
                                                                                                                                       add_numbers=TRUE,
                                                                                                                                       numbers_rot = 0,                                                                                                                                       border = FALSE, 
                                                                                                                                       #gp = gpar(fill = colour_list, col = "white", fontsize = plot_props$font_size, fontfamily=plot_props$font_family),
                                                                                                                                       gp = gpar(fill = "black", col = "white", fontsize = 4, fontfamily=plot_props$font_family),
                                                                                                                                       labels_gp = gpar(col = "black", fontsize = 4),
                                                                                                                                       numbers_gp = gpar(col = "black", fontsize = 4), 
                                                                                                                                       height = unit(1.25, "cm"),
                                                                                                                                       axis_param = list(labels_rot = 0),
                                                                                                                                       ),
                                                                          annotation_name_side = "left", annotation_name_rot = 0, annotation_name_gp = gpar(fontsize = plot_props$font_size, fontfamily=plot_props$font_family)
                                                                          ),
                                       
                                       right_annotation = rowAnnotation("Brain region\noverlap" = anno_barplot(comb_size(plot_table_thresholded),
                                                                                                                 #add_numbers=TRUE,
                                                                                                                 #numbers_rot=0,
                                                                                                                 ylim = c(0, max(comb_size(plot_table_thresholded))*1.1),
                                                                                                                 border = FALSE,
                                                                                                                 gp = gpar(fill = "black", fontsize = plot_props$font_size, fontfamily=plot_props$font_family),
                                                                                                                 labels_gp = gpar(col = "black", fontsize = 4), 
                                                                                                                 width = unit(1.25, "cm"),
                                                                                                                 axis_param = list(labels_rot = 0)
                                                                                                                 ),
                                                                          annotation_name_side = "bottom", annotation_name_gp = gpar(fontsize = plot_props$font_size, fontfamily=plot_props$font_family)
                                                                          ),
                                       #column_title = sprintf("%s", plot_title)
                                       )
    #co = column_order(upset_plot)
    #ro = row_order(upset_plot)
    
    #upset_plot
    
    # if the size of all combinations is exactly 1, then display the brain regions and remove right bar plot
    if ((length(unique(comb_size(plot_table_thresholded))) == 1) && (unique(comb_size(plot_table_thresholded)) == 1)) {
      # for some strange reason, we have to create the names in the wrong order
      row_name_square_plot = row_name_squares(sorted_rownames, object_height, 1.5, xticks_names, split_prop, colours, plot_props)
      
      # then, we have to change the order of the heatmap object because it must be the same as the order of the upset_plot object (which we alreadz sorted)
      row_name_square_plot@row_order = upset_plot@row_order
      
      # remove right annotaion
      upset_plot@right_annotation = NULL
      
      upset_plot = row_name_square_plot + upset_plot
      
    }
    
    plot_height = plot_props$image_height
    plot_width = plot_props$image_width
    png(sprintf("%s_intersection_th=%s_%sx%s.png", file_name_upset, th, plot_height, plot_width), width = plot_width, height = plot_height, units = plot_props$image_units, res = plot_props$dpi)
    upset_plot = draw(upset_plot, padding = unit(c(0.1, 0.5, 0.1, -0.5), "cm"));  # bottom, left, top and right margins
    #decorate_annotation("Brain region\noverlap", {grid.text(cs[co], x = 1:nc, y = unit(cs[co], "native") + unit(1, "mm"), gp = gpar(fontsize = 6), rot = 90, vjust = 0.35, hjust = 0.15, just = "bottom", default.units = "native")})
    #decorate_annotation("Num. of brain regions", {grid.text(ss[rev(ro)], x = unit(ss[rev(ro)], "native") + unit(1, "mm"), y = 1:length(ro), gp = gpar(fontsize = 6), vjust = 0.35, hjust = 0.15, just = "left", default.units = "native")})
    dev.off()
    svglite(sprintf("%s_intersection_th=%s_%sx%s.svg", file_name_upset, th, plot_height, plot_width), width = plot_width / 2.54, height = plot_height / 2.54)
    upset_plot = draw(upset_plot, padding = unit(c(0.1, 0.5, 0.1, -0.5), "cm"));  # bottom, left, top and right margins
    #decorate_annotation("Brain region\noverlap", {grid.text(cs[co], x = 1:nc, y = unit(cs[co], "native") + unit(1, "mm"), gp = gpar(fontsize = 6), rot = 90, vjust = 0.35, hjust = 0.15, just = "bottom", default.units = "native")})
    #decorate_annotation("Num. of brain regions", {grid.text(ss[rev(ro)], x = unit(ss[rev(ro)], "native") + unit(1, "mm"), y = 1:length(ro), gp = gpar(fontsize = 6), vjust = 0.35, hjust = 0.15, just = "left", default.units = "native")})
    dev.off()
    
  }
}

create_upset_plot = function(sig_neg_correlated_genes_tissue, split_prop, xticks_names, gene_colours, colours, plot_props, file_name_upset) {
  # save mat
  #write.table(list_to_matrix(sig_neg_correlated_genes_tissue), sprintf("%s.csv", file_name_tab), sep='\t', row.names = T)
  #write.xlsx(list_to_matrix(sig_neg_correlated_genes_tissue), sprintf("%s.xlsx", file_name_tab), colNames = TRUE, rowNames = TRUE, append = FALSE)
  
  # upset plot
  mode = "distinct"
  intersection_th_list = c(1,2,5,10)
  
  for (th in intersection_th_list) {
    plot_table = make_comb_mat(sig_neg_correlated_genes_tissue, mode=mode)
    #print(dim(plot_table))
    # skip if all combinattions are < th
    if (all(comb_size(plot_table) < th)) {
      #print(sprintf("skipping intersection_th=%s", th))
      next
    }
    plot_table_thresholded = plot_table[, comb_size(plot_table) >= th]
    #print(dim(plot_table_thresholded))
    
    #if (th == 5) {
    #  interesting_regions = c("pon", "plx", "svz", "olf", "vis", "med", "th", "cor", "cc")
    #  tmp = rownames(plot_table_thresholded)[rownames(plot_table_thresholded) %in% interesting_regions]
    #  plot_table_thresholded = plot_table_thresholded[tmp, ]
    #  # now filter the comb_size again!
    #  plot_table_thresholded = plot_table_thresholded[, comb_size(plot_table_thresholded) >= th]
    #  
    #  print(rownames(plot_table_thresholded))
    #}
    #print(dim(plot_table_thresholded))
    #print(rownames(plot_table_thresholded))
    
    # sort color vector according to plot_df
    #colour_list = unlist(colours[[split_prop]])[rownames(plot_table_thresholded)]
    
    
    # find which col belongs to which brain region only if the size of all combinations is exactly 1, then change the order of rows for the plot
    if ((length(unique(comb_size(plot_table_thresholded))) == 1) && (unique(comb_size(plot_table_thresholded)) == 1)) {
      # transform combination list into matrix
      mat = do.call(rbind, type.convert(strsplit(comb_name(plot_table_thresholded, readable = FALSE), ""), as.is = TRUE))
      # set colnames according to upset plot object. order is by implementation eqaul
      colnames(mat) = set_name(plot_table_thresholded)
      
      # sort according to image
      tmp = names(sort(colSums(mat), decreasing = TRUE))
      mat = mat[,tmp]
      
      mat_orig = list_to_matrix(sig_neg_correlated_genes_tissue)
      mat_orig = mat_orig[, colnames(mat)]
      
      # make list of row vectors
      target = split(mat, row(mat))
      orig = split(mat_orig, row(mat_orig))
      
      # find the ordering
      target_order = match(target, orig)
      #mat_orig[target_order, ]
      
      sorted_colnames = rownames(mat_orig[target_order, ])
      
      # sort plot rows and rownames alphabetically as given in the colours list
      ordered_colnames = names(colours[[split_prop]])[names(colours[[split_prop]]) %in% sorted_colnames]
      
      # extract the keys (number for reordering)
      target_col_order = match(ordered_colnames, sorted_colnames)
    } else {
      target_col_order = order(comb_size(plot_table_thresholded), decreasing = TRUE)
    }
    
    if (length(set_size(plot_table_thresholded)) >= 15) {
      gene_font_size = 6
    } else {
      gene_font_size = plot_props$font_size
    }
    
    object_height = 9
    object_width = 6
    
    upset_plot = ComplexHeatmap::UpSet(plot_table_thresholded, pt_size = unit(1, "mm"), comb_order = target_col_order, lwd = 1, #pt_size = unit(1, "mm"), lwd = 1, bg_col = c("#ECECEC", "#CBCBCB")
                                       height = unit(object_height, "cm"), width = unit(object_width, "cm"),
                                       row_names_side = "left", row_names_gp = gpar(fontsize = gene_font_size, fontfamily=plot_props$font_family),
                                       #top_annotation = NULL,
                                       right_annotation = rowAnnotation("Num. of\nbrain reg.\n(with sig.\nneg. corr.)" = anno_barplot(set_size(plot_table_thresholded), 
                                                                                                                                      add_numbers=TRUE,
                                                                                                                                      numbers_rot = 0,                                                                                                                                       border = FALSE, 
                                                                                                                                      #gp = gpar(fill = colour_list, col = "white", fontsize = plot_props$font_size, fontfamily=plot_props$font_family),
                                                                                                                                      gp = gpar(fill = "black", col = "white", fontsize = 4, fontfamily=plot_props$font_family),
                                                                                                                                      labels_gp = gpar(col = "black", fontsize = 6),
                                                                                                                                      numbers_gp = gpar(col = "black", fontsize = 6), 
                                                                                                                                      width = unit(1.25, "cm"),
                                                                                                                                      axis_param = list(labels_rot = 0),
                                       ),
                                       annotation_name_side = "top", annotation_name_rot = 0, annotation_name_gp = gpar(fontsize = plot_props$font_size, fontfamily=plot_props$font_family)
                                       ),
                                       
                                       top_annotation = HeatmapAnnotation("Brain region\noverlap" = anno_barplot(comb_size(plot_table_thresholded),
                                                                                                               #add_numbers=TRUE,
                                                                                                               #numbers_rot=0,
                                                                                                               ylim = c(0, max(comb_size(plot_table_thresholded))*1.1),
                                                                                                               border = FALSE,
                                                                                                               gp = gpar(fill = "black", fontsize = plot_props$font_size, fontfamily=plot_props$font_family),
                                                                                                               labels_gp = gpar(col = "black", fontsize = 6), 
                                                                                                               height = unit(1.25, "cm"),
                                                                                                               axis_param = list(labels_rot = 0)
                                       ),
                                       annotation_name_side = "left", annotation_name_gp = gpar(fontsize = plot_props$font_size, fontfamily=plot_props$font_family)
                                       ),
                                       #column_title = sprintf("%s", plot_title)
    )
    #co = column_order(upset_plot)
    #ro = row_order(upset_plot)
    
    #upset_plot
    
    # if the size of all combinations is exactly 1, then display the brain regions and remove right bar plot
    if ((length(unique(comb_size(plot_table_thresholded))) == 1) && (unique(comb_size(plot_table_thresholded)) == 1)) {
      # for some strange reason, we have to create the names in the wrong order
      col_name_square_plot = col_name_squares(sorted_colnames, 1.75, object_width, xticks_names, split_prop, colours, plot_props)
      
      # then, we have to change the order of the heatmap object because it must be the same as the order of the upset_plot object (which we alreadz sorted)
      col_name_square_plot@column_order = target_col_order
      
      # Get the order of the rows as they appear in the heatmap
      row_order = row_order(upset_plot)
      # Retrieve the row names in that order
      #ordered_row_names = rownames(plot_table_thresholded)[row_order]
      # row_name_square_plot = row_name_squares_adj(ordered_row_names, object_height, 1.5, xticks_names, "genes", gene_colours, plot_props)
      row_name_square_plot = row_name_squares_adj(rownames(plot_table_thresholded), object_height, 1.5, xticks_names, "genes", gene_colours, plot_props)
      row_name_square_plot@row_order = row_order
      
      # for some strange reason, we have to create the names in the wrong order
      #row_name_square_plot = row_name_squares(sorted_rownames, object_height, 1.5, xticks_names, split_prop, colours, plot_props)
      
      # then, we have to change the order of the heatmap object because it must be the same as the order of the upset_plot object (which we alreadz sorted)
      #row_name_square_plot@row_order = upset_plot@row_order
      
      # remove bottom annotaion
      upset_plot@top_annotation = NULL
      
      # upset_plot_with_rownames = row_name_square_plot + upset_plot
      upset_plot = col_name_square_plot %v% upset_plot

    }
    
    plot_height = 12
    plot_width = plot_props$image_width
    png(sprintf("%s_intersection_th=%s_%sx%s.png", file_name_upset, th, plot_height, plot_width), width = plot_width, height = plot_height, units = plot_props$image_units, res = plot_props$dpi)
    upset_plot = draw(upset_plot, padding = unit(c(0.1, 0, 0.1, 0), "cm"));  # bottom, left, top and right margins
    #decorate_annotation("Brain region\noverlap", {grid.text(cs[co], x = 1:nc, y = unit(cs[co], "native") + unit(1, "mm"), gp = gpar(fontsize = 6), rot = 90, vjust = 0.35, hjust = 0.15, just = "bottom", default.units = "native")})
    #decorate_annotation("Num. of brain regions", {grid.text(ss[rev(ro)], x = unit(ss[rev(ro)], "native") + unit(1, "mm"), y = 1:length(ro), gp = gpar(fontsize = 6), vjust = 0.35, hjust = 0.15, just = "left", default.units = "native")})
    dev.off()
    svglite(sprintf("%s_intersection_th=%s_%sx%s.svg", file_name_upset, th, plot_height, plot_width), width = plot_width / 2.54, height = plot_height / 2.54)
    upset_plot = draw(upset_plot, padding = unit(c(0.1, 0, 0.1, 0), "cm"));  # bottom, left, top and right margins
    #decorate_annotation("Brain region\noverlap", {grid.text(cs[co], x = 1:nc, y = unit(cs[co], "native") + unit(1, "mm"), gp = gpar(fontsize = 6), rot = 90, vjust = 0.35, hjust = 0.15, just = "bottom", default.units = "native")})
    #decorate_annotation("Num. of brain regions", {grid.text(ss[rev(ro)], x = unit(ss[rev(ro)], "native") + unit(1, "mm"), y = 1:length(ro), gp = gpar(fontsize = 6), vjust = 0.35, hjust = 0.15, just = "left", default.units = "native")})
    dev.off()
    # save rownames separatly
    if ((length(unique(comb_size(plot_table_thresholded))) == 1) && (unique(comb_size(plot_table_thresholded)) == 1)) {
      png(sprintf("%s_intersection_th=%s_%sx%s_rownames.png", file_name_upset, th, plot_height, plot_width), width = plot_width, height = plot_height, units = plot_props$image_units, res = plot_props$dpi)
      draw(row_name_square_plot, padding = unit(c(0.1, 0, 0.1, 0), "cm"));  # bottom, left, top and right margins
      #decorate_annotation("Brain region\noverlap", {grid.text(cs[co], x = 1:nc, y = unit(cs[co], "native") + unit(1, "mm"), gp = gpar(fontsize = 6), rot = 90, vjust = 0.35, hjust = 0.15, just = "bottom", default.units = "native")})
      #decorate_annotation("Num. of brain regions", {grid.text(ss[rev(ro)], x = unit(ss[rev(ro)], "native") + unit(1, "mm"), y = 1:length(ro), gp = gpar(fontsize = 6), vjust = 0.35, hjust = 0.15, just = "left", default.units = "native")})
      dev.off()
      svglite(sprintf("%s_intersection_th=%s_%sx%s_rownames.svg", file_name_upset, th, plot_height, plot_width), width = plot_width / 2.54, height = plot_height / 2.54)
      draw(row_name_square_plot, padding = unit(c(0.1, 1, 0.1, -1), "cm"));  # bottom, left, top and right margins
      #decorate_annotation("Brain region\noverlap", {grid.text(cs[co], x = 1:nc, y = unit(cs[co], "native") + unit(1, "mm"), gp = gpar(fontsize = 6), rot = 90, vjust = 0.35, hjust = 0.15, just = "bottom", default.units = "native")})
      #decorate_annotation("Num. of brain regions", {grid.text(ss[rev(ro)], x = unit(ss[rev(ro)], "native") + unit(1, "mm"), y = 1:length(ro), gp = gpar(fontsize = 6), vjust = 0.35, hjust = 0.15, just = "left", default.units = "native")})
      dev.off()
    }
    
  }
}


#---------------------------------- Load data ----------------------------------
data_input = snakemake@params$data_input
feature = snakemake@params$feature
gene_list_cas_folder_path = snakemake@params$gene_list_cas_folder_path
gene_list_target_genes_folder_path = snakemake@params$gene_list_target_genes_folder_path
split_prop = snakemake@params$split_prop
corr_params = snakemake@params$corr_params
gene_colours_file_path = snakemake@params$colour_genes
xticks_names = snakemake@params$xticks_names
colours = snakemake@params$colors
plot_props = snakemake@params$plots_props
results_folder = snakemake@params$results_folder

# load files
annot = fread(sprintf("%s/annotation_%s_%s.csv", data_input$annotation_data_folder, data_input$data_sub_set, data_input$feature_filtering), sep='\t', colClasses=c(ID="character"))
gene_colours = fread(gene_colours_file_path, sep='\t') 

input_folder = sprintf("%s/%s_%s/results_%s/matrices/%s_target_genes_corr", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, feature)

#output_folder_venn = sprintf("%s/%s_%s/results_%s/figures/%s_target_genes_corr/venn_plots", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, feature)
#dir.create(output_folder_venn, recursive=TRUE)
output_folder_upset = sprintf("%s/%s_%s/results_%s/figures/%s_target_genes_corr/upset_plots", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, feature)
dir.create(output_folder_upset, recursive=TRUE)
output_folder_tab = sprintf("%s/%s_%s/results_%s/matrices/%s_target_genes_corr/upset_plots", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, feature)
dir.create(output_folder_tab, recursive=TRUE)


#------------------------------------ Script ----------------------------------- 
for (target_genes_source in c("target_genes_all", "target_genes_without_weak")) {
  for (method in corr_params$method) {
    #dir.create(sprintf("%s/%s/%s", output_folder_venn, target_genes_source, method), recursive=TRUE)
    dir.create(sprintf("%s/%s/%s", output_folder_upset, target_genes_source, method), recursive=TRUE)
    dir.create(sprintf("%s/%s/%s", output_folder_tab, target_genes_source, method), recursive=TRUE)
    
    for (corr_th in corr_params$corr_th) {
      
      sig_neg_correlated_genes = c()
      for (ti in unique(annot[[split_prop]])) {
        if (!file.exists(sprintf("%s/%s/%s/correlation_%s.csv", input_folder, target_genes_source, method, ti))) {
          print(ti)
          next
        }
        # load correlation data
        correlation = as.data.frame(fread(sprintf("%s/%s/%s/correlation_%s.csv", input_folder, target_genes_source, method, ti), sep='\t', header = T))
        correlation$rna = c()
        padj = as.data.frame(fread(sprintf("%s/%s/%s/padj_%s.csv", input_folder, target_genes_source, method, ti), sep='\t', header = T))
        padj$rna = c()
        
        # if padj is not sorted like correlation we do so
        padj = padj[colnames(correlation)]
        
        # filter for sig. neg. correlated genes
        keep_genes = ((correlation <= -corr_th) & (padj < corr_params$sig_lvl))
        keep_genes[is.na(keep_genes)] = FALSE
        sig_neg_correlated_genes[[ti]] = colnames(correlation[,keep_genes, drop = FALSE])
      }
      
      # Get all unique elements across the vectors
      unique_elements = unique(unlist(sig_neg_correlated_genes))
      
      # Create a new list where each unique element is the name and the values are the vector names containing that element
      sig_neg_correlated_genes_tissue = setNames(lapply(unique_elements, function(x) {
        names(sig_neg_correlated_genes)[sapply(sig_neg_correlated_genes, function(vec) x %in% vec)]
      }), unique_elements)
      
      file_name_upset = sprintf("%s/%s/%s/upset_%s_corr_m=%s_corr_th=%s_rot", output_folder_upset, target_genes_source, method, feature, method, corr_th)
      file_name_tab = sprintf("%s/%s/%s/upset_%s_corr_m=%s_corr_th=%s", output_folder_tab, target_genes_source, method, feature, method, corr_th)
      create_upset_plot_rot(sig_neg_correlated_genes_tissue, split_prop, xticks_names, colours, plot_props, file_name_upset, file_name_tab)
      
      file_name_upset = sprintf("%s/%s/%s/upset_%s_corr_m=%s_corr_th=%s", output_folder_upset, target_genes_source, method, feature, method, corr_th)
      create_upset_plot(sig_neg_correlated_genes_tissue, split_prop, xticks_names, gene_colours, colours, plot_props, file_name_upset)
      
    }
  }
}


#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])


