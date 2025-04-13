suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(svglite))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(latex2exp))
suppressPackageStartupMessages(library(openxlsx))


source("./scripts/iwanthue.R")

#snakemake = readRDS("snakemake_mirna_number_per_tissue_bar_plot.rds")

set.seed(snakemake@params$parameters_porps$set_seed)

print("mirna_number_per_tissue_bar_plot") 


#---------------------------------- Load data ---------------------------------- 
data_input = snakemake@params$data_input
tissue = snakemake@params$tissue
detection_rate_list = snakemake@params$detection_rate_list
min_detect_raw = snakemake@params$min_detect_raw
xticks_names = snakemake@params$xticks_names
plots_props = snakemake@params$plots_props
colours = snakemake@params$colors
results_folder = snakemake@params$results_folder

# save rdata
#saveRDS(snakemake, file = "snakemake_mirna_number_per_tissue_bar_plot.rds")
#stop()

expr = fread(sprintf("%s_%s/%s_%s_quantification_%s.detection_rate_%s_per_%s.csv", data_input$raw_data_folder, data_input$data_sub_set, data_input$rna_class, data_input$feature_filtering, data_input$norm, data_input$detection_rate, data_input$detection_group), sep='\t', header=T)
annot = fread(sprintf("%s/annotation_%s_%s.csv", data_input$annotation_data_folder, data_input$data_sub_set, data_input$feature_filtering), sep='\t', colClasses=c(ID="character")) 
detect = fread(sprintf("%s_%s/%s_%s_detection_matrix_min_detect=%s.csv", data_input$raw_data_folder, data_input$data_sub_set, data_input$rna_class, data_input$feature_filtering, min_detect_raw), sep='\t', header=T)


# --------------------------------- Script -------------------------------------
sample_cols = colnames(expr)[2:ncol(expr)]

for (detection_rate in detection_rate_list) {
  
  # output folder
  output_folder_bar = sprintf("%s/%s_%s/results_%s/figures/expressed_per_%s/bar_plots", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, tissue)
  dir.create(output_folder_bar, recursive=TRUE)
  output_folder_upset = sprintf("%s/%s_%s/results_%s/figures/expressed_per_%s/upset_plots", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, tissue)
  dir.create(output_folder_upset, recursive=TRUE)
  output_folder_tbl = sprintf("%s/%s_%s/results_%s/matrices/expressed_per_%s", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, tissue)
  dir.create(output_folder_tbl, recursive=TRUE)
  
  to_keep_list = c()
  number_of_mirnas = c()
  for (ti in unique(annot[[tissue]])) {
    to_keep = rep(F, nrow(detect))
    sub = (annot[annot[[tissue]] == ti,][[data_input$identifier_column]])
    to_keep = to_keep | (rowMeans(detect[,..sub]) >= (as.numeric(gsub("p", "", detection_rate)) / 100))
    
    features_to_keep = detect[to_keep][[data_input$rna_class]]
    
    tmp = (expr[[data_input$rna_class]] %in% features_to_keep)
    expr_to_keep = expr[tmp]
    
    to_keep_list[[ti]] = expr_to_keep[[data_input$rna_class]]
    
    write.table(to_keep_list[[ti]], sprintf("%s/expressed_%s_in_%s=%s.csv", output_folder_tbl, data_input$rna_class, tissue, ti), sep='\t', row.names = F, col.names = F)
    write.xlsx(to_keep_list[[ti]], sprintf("%s/expressed_%s_in_%s=%s.xlsx", output_folder_tbl, data_input$rna_class, tissue, ti), colNames = FALSE, rowNames = FALSE, append = FALSE)
    
    number_of_mirnas = append(number_of_mirnas, dim(expr_to_keep)[1])
  }
  
  number_of_mirnas_df = data.frame(number_of_mirnas)
  number_of_mirnas_df$tissue = unique(annot[[tissue]])
  number_of_mirnas_df$tissue_readable = xticks_names[[tissue]][unique(annot[[tissue]])]
  
  
  # order decreasing by column "counts" (only values where column alignment == aligned)
  # get the numeric order
  ord = order(number_of_mirnas_df$number_of_mirnas, decreasing = TRUE)
  
  # get the names
  x_order = number_of_mirnas_df[ord,]$tissue_readable
  
  # bar plot
  #p = ggplot(number_of_mirnas_df, aes(x=reorder(tissue_readable, -number_of_mirnas), y=number_of_mirnas, fill=tissue)) + 
   p = ggplot(number_of_mirnas_df, aes(x=factor(tissue_readable, levels=x_order), y=number_of_mirnas, fill=tissue)) + 
    geom_bar(stat="identity") +
    scale_fill_manual(values=colours[[tissue]]) +
    scale_y_continuous(expand = c(0,0)) +
    xlab("") +
    ylab(sprintf("Number of %ss\n(≥%s counts for %s%% of the samples)", data_input$rna_class, min_detect_raw, detection_rate)) +
    theme_classic() +
    theme(legend.position="none",
          text = element_text(family = plots_props$font_family, size = plots_props$font_size),
          axis.text = element_text(angle = 0, family = plots_props$font_family, size = plots_props$font_size),
          axis.text.x = element_text(family = plots_props$font_family, size = plots_props$font_size, angle = 90, vjust = 0.5, hjust=1), 
          axis.title = element_text(family = plots_props$font_family, size = plots_props$font_size),
          plot.title = element_text(family = plots_props$font_family, size = plots_props$font_size_header),
          legend.text = element_text(family = plots_props$font_family, size = plots_props$font_size),
          plot.margin = unit(c(0.5,0.25,-0.25,0.25), "cm"), #t, r, b, l
    )
  
  ggsave(file = sprintf("%s/number_of_expressed_miRNAs_per_%s_th=%s_counts_detection_rate=%s.svg", output_folder_bar, tissue, min_detect_raw, detection_rate), plot = p, width = plots_props$image_width, height = plots_props$image_height, unit = plots_props$image_units, dpi = plots_props$dpi)
  ggsave(file = sprintf("%s/number_of_expressed_miRNAs_per_%s_th=%s_counts_detection_rate=%s.png", output_folder_bar, tissue, min_detect_raw, detection_rate), plot = p, width = plots_props$image_width, height = plots_props$image_height, unit = plots_props$image_units, dpi = plots_props$dpi)
  
  # upset plot
  mode = "distinct"
  intersection_th_list = c(1,5)
  
  for (th in intersection_th_list) {
    plot_table = make_comb_mat(to_keep_list, mode=mode)
    plot_table_thresholded = plot_table[comb_size(plot_table) >= th]
    
    cs = comb_size(plot_table_thresholded)
    ss = set_size(plot_table_thresholded)
    nc = ncol(plot_table_thresholded)
    
    # sort color vector according to plot_df
    colour_list = unlist(colours[[tissue]])[rownames(plot_table_thresholded)]
    
    rownames(plot_table_thresholded) = xticks_names[[tissue]][rownames(plot_table_thresholded)]
    
    upset_plot = ComplexHeatmap::UpSet(plot_table_thresholded, pt_size = unit(1, "mm"), lwd = 1, comb_order = order(comb_size(plot_table_thresholded), decreasing = TRUE), # bg_col = c("#ECECEC", "#CBCBCB")
                               height = unit(3.5, "cm"), width = unit(4.5, "cm"),
                               row_names_side = "left", row_names_gp = gpar(fontsize = plots_props$font_size, fontfamily=plots_props$font_family), 
                               top_annotation = HeatmapAnnotation("miRNA\noverlap" = anno_barplot(comb_size(plot_table_thresholded),
                                                                                                  #add_numbers=TRUE,
                                                                                                  ylim = c(0, max(comb_size(plot_table_thresholded))*1.1), 
                                                                                                  border = FALSE, 
                                                                                                  gp = gpar(fill = "black", fontsize = 6, fontfamily=plots_props$font_family),
                                                                                                  labels_gp = gpar(col = "black", fontsize = 4), 
                                                                                                  height = unit(1, "cm")
                                                                                                  ), 
                                                                  annotation_name_side = "left", annotation_name_rot = 90, annotation_name_gp = gpar(fontsize = plots_props$font_size, fontfamily=plots_props$font_family)
                                                                  ),
                               
                               right_annotation = rowAnnotation("Num. of expr. miRNAs\nper brain region" = anno_barplot(set_size(plot_table_thresholded), 
                                                                                                                             border = FALSE, 
                                                                                                                             gp = gpar(fill = colour_list, col = "white", fontsize = plots_props$font_size, fontfamily=plots_props$font_family),
                                                                                                                             labels_gp = gpar(col = "black", fontsize = 6), 
                                                                                                                             width = unit(2.5, "cm"),
                                                                                                                             axis_param = list(labels_rot = 0),
                                                                                                                             ),
                                                                annotation_name_side = "bottom", annotation_name_gp = gpar(fontsize = plots_props$font_size, fontfamily=plots_props$font_family)
                                                                ),
                               #column_title = sprintf("%s", plot_title)
                               )
    co = column_order(upset_plot)
    ro = row_order(upset_plot)
    
    png(sprintf("%s/number_of_expressed_miRNAs_per_%s_th=%s_counts_detection_rate=%s_intersection_th=%s_9x6.png", output_folder_upset, tissue, min_detect_raw, detection_rate, th), width = plots_props$image_width, height = plots_props$image_height, units = plots_props$image_units, res = plots_props$dpi)
    upset_plot = draw(upset_plot, padding = unit(c(-0.1, -0.25, 0.25, 0.25), "cm"));  # bottom, left, top and right margins
    decorate_annotation("miRNA\noverlap", {grid.text(cs[co], x = 1:nc, y = unit(cs[co], "native") + unit(1, "mm"), gp = gpar(fontsize = 6), rot = 90, vjust = 0.35, hjust = 0.15, just = "bottom", default.units = "native")})
    dev.off()
    svglite(sprintf("%s/number_of_expressed_miRNAs_per_%s_th=%s_counts_detection_rate=%s_intersection_th=%s_9x6.svg", output_folder_upset, tissue, min_detect_raw, detection_rate, th), width = plots_props$image_width / 2.54, height = plots_props$image_height / 2.54)
    upset_plot = draw(upset_plot, padding = unit(c(-0.1, -0.25, 0.25, 0.25), "cm"));  # bottom, left, top and right margins
    decorate_annotation("miRNA\noverlap", {grid.text(cs[co], x = 1:nc, y = unit(cs[co], "native") + unit(1, "mm"), gp = gpar(fontsize = 6), rot = 90, vjust = 0.35, hjust = 0.15, just = "bottom", default.units = "native")})
    dev.off()
    
    upset_plot = ComplexHeatmap::UpSet(plot_table_thresholded, pt_size = unit(1, "mm"), lwd = 1, comb_order = order(comb_size(plot_table_thresholded), decreasing = TRUE), # bg_col = c("#ECECEC", "#CBCBCB")
                                       height = unit(3.5, "cm"), width = unit(2, "cm"),
                                       row_names_side = "left", row_names_gp = gpar(fontsize = plots_props$font_size, fontfamily=plots_props$font_family), 
                                       top_annotation = HeatmapAnnotation("miRNA\noverlap" = anno_barplot(comb_size(plot_table_thresholded),
                                                                                                          #add_numbers=TRUE,
                                                                                                          ylim = c(0, max(comb_size(plot_table_thresholded))*1.1), 
                                                                                                          border = FALSE, 
                                                                                                          gp = gpar(fill = "black", fontsize = plots_props$font_size, fontfamily=plots_props$font_family),
                                                                                                          labels_gp = gpar(col = "black", fontsize = 4), 
                                                                                                          height = unit(1, "cm")
                                       ), 
                                       annotation_name_side = "left", annotation_name_rot = 90, annotation_name_gp = gpar(fontsize = plots_props$font_size, fontfamily=plots_props$font_family)
                                       ),
                                       
                                       right_annotation = rowAnnotation("Num. of expr. miRNAs\nper brain region" = anno_barplot(set_size(plot_table_thresholded), 
                                                                                                                                border = FALSE, 
                                                                                                                                gp = gpar(fill = colour_list, col = "white", fontsize = plots_props$font_size, fontfamily=plots_props$font_family),
                                                                                                                                labels_gp = gpar(col = "black", fontsize = 6), 
                                                                                                                                width = unit(2.5, "cm"),
                                                                                                                                axis_param = list(labels_rot = 0),
                                       ),
                                       annotation_name_side = "bottom", annotation_name_gp = gpar(fontsize = plots_props$font_size, fontfamily=plots_props$font_family)
                                       ),
                                       #column_title = sprintf("%s", plot_title)
    )
    co = column_order(upset_plot)
    ro = row_order(upset_plot)
    
    png(sprintf("%s/number_of_expressed_miRNAs_per_%s_th=%s_counts_detection_rate=%s_intersection_th=%s_6x6.png", output_folder_upset, tissue, min_detect_raw, detection_rate, th), width = 2/3 * plots_props$image_width, height = plots_props$image_height, units = plots_props$image_units, res = plots_props$dpi)
    upset_plot = draw(upset_plot, padding = unit(c(-0.1, -0.25, 0.25, 0.25), "cm"));  # bottom, left, top and right margins
    decorate_annotation("miRNA\noverlap", {grid.text(cs[co], x = 1:nc, y = unit(cs[co], "native") + unit(1, "mm"), gp = gpar(fontsize = 6), rot = 90, vjust = 0.35, hjust = 0.15, just = "bottom", default.units = "native")})
    dev.off()
    svglite(sprintf("%s/number_of_expressed_miRNAs_per_%s_th=%s_counts_detection_rate=%s_intersection_th=%s_6x6.svg", output_folder_upset, tissue, min_detect_raw, detection_rate, th), width = 2/3 * plots_props$image_width / 2.54, height = plots_props$image_height / 2.54)
    upset_plot = draw(upset_plot, padding = unit(c(-0.1, -0.25, 0.25, 0.25), "cm"));  # bottom, left, top and right margins
    decorate_annotation("miRNA\noverlap", {grid.text(cs[co], x = 1:nc, y = unit(cs[co], "native") + unit(1, "mm"), gp = gpar(fontsize = 6), rot = 90, vjust = 0.35, hjust = 0.15, just = "bottom", default.units = "native")})
    dev.off()
  }
}


#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])
