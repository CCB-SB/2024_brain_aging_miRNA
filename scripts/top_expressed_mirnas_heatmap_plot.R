suppressMessages(library(data.table))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(viridisLite))
suppressPackageStartupMessages(library(svglite))
suppressPackageStartupMessages(library(circlize))
suppressMessages(library(openxlsx))
suppressMessages(library(stringi))


# save rdata
#saveRDS(snakemake, file = "snakemake_top_expressed_mirnas_heatmap_plot.rds")
#stop()

#snakemake = readRDS("snakemake_top_expressed_mirnas_heatmap_plot.rds")

data.table::setDTthreads(1)
set.seed(snakemake@params$parameters_props$set_seed)

print("top_expressed_mirnas_heatmap_plot")


#---------------------------------- Functions ----------------------------------
plot_heatmap = function(df, row_names, i, identifier_column, prop, scientific_numbers, legend_title, colors, xticks_names, plot_props, width, height, filename_fig){
  #show_rownames = nrow(df) < 50
  show_rownames = TRUE
  show_legend = TRUE
  
  df = as.matrix(df)
  rownames(df) = row_names
  
  annotation_col = data.frame(dummy=1:ncol(df))
  #annotation_col[[prop]] = c(as.character(annot[[prop]][match(colnames(df), annot[[prop]])]))
  annotation_col[[prop]] = c(as.character(annot[[prop]][match(colnames(df), annot[[identifier_column]])]))
  rownames(annotation_col) = colnames(df)
  annotation_col$dummy = NULL
  
  annotation_colors = list()
  ucolors <- unlist(colors[[prop]])
  if(is.vector(ucolors) && all(annotation_col[[prop]] %in% names(ucolors))) {
    annotation_colors[[prop]] = ucolors
  }
  
  column_ha = HeatmapAnnotation(df = annotation_col, col = annotation_colors, gp = gpar(col = "white"), simple_anno_size = unit(0.15, "cm"), show_annotation_name = FALSE, show_legend = FALSE,
                                annotation_name_side = "left", annotation_name_gp = gpar(fontsize = plot_props$font_size, fontface = "bold", fontfamily=plot_props$font_family))
  
  if (length(unique(as.numeric(df))) <= 2) {
    #legend_ticks = seq(from = 0, to = max(t(df)), length.out = 5)
    maxi = max(df)
    mini = min(df)
    #col_fun = colorRamp2(c(0, maxi), c(viridis(100)[1],viridis(100)[100]))
    col_fun = colorRamp2(c(0, maxi), c("#DFDFB9","#525200"))
    #col_fun = colorRamp2(c(0,maxi), hcl_palette = "Greens", reverse = TRUE)
    col_fun_legend = col_fun
  } else {
    maxi = max(df)
    mini = min(df)
    #col_fun = colorRamp2(c(0,maxi), hcl_palette = "Greens", reverse = TRUE)
    col_fun = colorRamp2(c(mini,maxi), c("#DFDFB9","#525200"), reverse = TRUE)
    col_fun_legend = col_fun
  }
  
  if (show_legend == FALSE) {
    height_heatmap = height + 2
  } else {
    height_heatmap = height
  }
  
  # remove mmu-miR
  #rownames(df) = gsub("mmu-miR-", "", rownames(df))
  #rownames(df) = gsub("mmu-let-", "", rownames(df))
  # remove mmu-
  rownames(df) = gsub("mmu-", "", rownames(df))

  #df = df[rev(rownames(df)),]
  
  # features rows
  p = ComplexHeatmap::Heatmap(df, col = col_fun,
                              #rect_gp = gpar(type = "none"), 
                              #height = unit(0.005, "cm") * nrow(df), width = unit(0.3, "cm") * ncol(df),
                              height = unit(height_heatmap -3 , "cm"), width = unit(width - 6, "cm"),
                              # row settings
                              show_row_names = show_rownames,
                              row_names_side = "left", row_names_gp = gpar(fontsize = 4, fontfamily=plot_props$font_family),
                              #cluster_rows = FALSE,
                              cluster_rows = function(m) hclust(dist(m), method = "complete"), clustering_distance_rows = "euclidean", 
                              row_dend_reorder = TRUE,
                              show_row_dend = TRUE,
                              # column settings
                              show_column_names = TRUE,
                              cluster_columns = FALSE,
                              #cluster_columns = function(m) hclust(dist(m), method = "complete"), clustering_distance_columns = "euclidean", 
                              column_dend_reorder = FALSE,
                              #show_column_dend = FALSE,
                              column_names_side = "top", column_names_gp = gpar(fontsize = plot_props$font_size, fontfamily=plot_props$font_family),
                              # other settings
                              show_heatmap_legend = FALSE,
                              top_annotation = column_ha
  )
  
  if (length(unique(as.numeric(df))) <= 2) {
    # discrete legend
    legend_ticks = c(0, maxi)
    legend_tick_locations = c(0, 1)
    legend_tick_labels = c("0" = sprintf("|zscore| < %s", zscore_th), "1" = sprintf("|zscore| ≥ %s", zscore_th))
    legend = ComplexHeatmap::Legend(title = legend_title, at = legend_tick_locations, labels=legend_tick_labels, direction="horizontal",
                                    legend_gp = gpar(fill = col_fun_legend(legend_tick_locations)),
                                    title_gp = gpar(fontsize = plot_props$font_size, fontface = "bold", fontfamily=plot_props$font_family), labels_gp = gpar(fontsize = plot_props$font_size, fontfamily=plot_props$font_family))
    
  } else {
    # continuous (colorbar) legend
    #legend_ticks = round(seq(from = min(t(df)), to = max(t(df))))#, length.out = 5)). # this produces too many tiks in the legend!

    if (max(t(df)) -  min(t(df)) < 5) {
      length_out = ceiling(max(t(df)) - min(t(df)))
    } else {
      length_out = 5
    }

    legend_ticks = round(seq(from = min(t(df)), to = max(t(df)), length.out = length_out))
    legend_tick_locations = round(legend_ticks, digits = 2)
    
    if (scientific_numbers == TRUE) {
      legend_ticks = formatC(legend_ticks, format="e", digits = 2)
      legend = ComplexHeatmap::Legend(title = legend_title, col_fun = col_fun_legend, at = legend_tick_locations, labels = legend_ticks, direction="horizontal", legend_width = unit(5, plot_props$image_units),
                                      title_gp = gpar(fontsize = plot_props$font_size, fontface = "bold", fontfamily=plot_props$font_family), labels_gp = gpar(fontsize = 6, fontfamily=plot_props$font_family))
    } else {
      legend = ComplexHeatmap::Legend(title = legend_title, col_fun = col_fun_legend, at = legend_tick_locations, labels = legend_ticks, direction="horizontal", legend_width = unit(3, plot_props$image_units),
                                      title_gp = gpar(fontsize = plot_props$font_size, fontface = "bold", fontfamily=plot_props$font_family), labels_gp = gpar(fontsize = 6, fontfamily=plot_props$font_family))
    }
  }
  
  if (i == 25) {
    image_height_scaling_factor = 1
  } else if (i == 50) {
    image_height_scaling_factor = 1
  } else {
    image_height_scaling_factor = 2
  }
  
  png(sprintf("%s_%ss.png", filename, data_input$rna_class), width = 2/3 * width, height = image_height_scaling_factor * height, units = plot_props$image_units, res = plot_props$dpi)
  if (show_legend == TRUE) {
    draw(p, show_annotation_legend = FALSE, padding = unit(c(1.5, 0.25, 0, 0.25), "cm")) #top, right, bottom, left
    draw(legend, x = unit(1.3, "cm"), y = unit(0.02, "cm"), just = c("left", "bottom"))
  } else {
    draw(p, show_annotation_legend = FALSE, padding = unit(c(0.25, 0.25, 0, 0.25), "cm"))
  }
  dev.off()
  svglite(sprintf("%s_%ss.svg", filename, data_input$rna_class), width = 2/3 * width / 2.54, height = image_height_scaling_factor * height / 2.54)
  if (show_legend == TRUE) {
    draw(p, show_annotation_legend = FALSE, padding = unit(c(1.5, 0.25, 0, 0.25), "cm"))
    draw(legend, x = unit(1.3, "cm"), y = unit(0.02, "cm"), just = c("left", "bottom"))
  } else {
    draw(p, show_annotation_legend = FALSE, padding = unit(c(0.25, 0.25, 0, 0.25), "cm"))
  }
  
  dev.off()
}


#---------------------------------- Load data ----------------------------------
data_input = snakemake@params$data_input
identifier_column = data_input$identifier_column
method = snakemake@params$method
top_list = snakemake@params$top_list
sample_order = snakemake@params$sample_order
properties = snakemake@params$properties
xticks_names = snakemake@params$xticks_names
colors = snakemake@params$colors
plot_props = snakemake@params$plots_props
results_folder = snakemake@params$results_folder
 
expr = fread(sprintf("%s_%s/%s_%s_quantification_%s.detection_rate_%s_per_%s.csv", data_input$raw_data_folder, data_input$data_sub_set, data_input$rna_class, data_input$feature_filtering, data_input$norm, data_input$detection_rate, data_input$detection_group), sep='\t', header=T)
expr_log10 = fread(sprintf("%s_%s/%s_%s_quantification_%s_log10.detection_rate_%s_per_%s.csv", data_input$raw_data_folder, data_input$data_sub_set, data_input$rna_class, data_input$feature_filtering, data_input$norm, data_input$detection_rate, data_input$detection_group), sep='\t', header=T)
annot = fread(sprintf("%s/annotation_%s_%s.csv", data_input$annotation_data_folder, data_input$data_sub_set, data_input$feature_filtering), sep='\t', colClasses=c(ID="character")) 

output_folder = sprintf("%s/%s_%s/results_%s/figures/top_expressed/heatmap_complex", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set)
dir.create(output_folder, recursive=TRUE)


#------------------------------------ Script ----------------------------------- 
for(i in top_list) {
  print(i)
  for(prop in properties) {
    top_expressed = fread(sprintf("%s/%s_%s/results_%s/matrices/most_expressed/top_%s.csv", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, i))
    
    # expr
    scientific_numbers = TRUE
    
    expr_filtered = expr[expr[[data_input$rna_class]] %in% top_expressed$miRNA,]
    
    row_names = expr_filtered[[data_input$rna_class]]
    expr_filtered[[data_input$rna_class]] = c()
    
    # sort colnames
    expr_filtered = expr_filtered[, ..sample_order]
    
    if (grepl(" ", xticks_names$categories[[data_input$detection_group]]) == TRUE) {
      detection_group = strsplit(xticks_names$categories[[data_input$detection_group]], " ")[[1]][1]
    } else {
      detection_group = data_input$detection_group
    }
    legend_title = sprintf("Expr (%s, detection rate = %s%% per %s)", gsub("_norm", "", data_input$norm), gsub("p", "", data_input$detection_rate), detection_group)
    filename = sprintf("%s/top_%s", output_folder, i)
  
    plot_heatmap(expr_filtered, row_names, i, identifier_column, prop, scientific_numbers, legend_title, colors, xticks_names, plot_props, plot_props$image_width, plot_props$image_height, filename)
    
    # expr log10
    scientific_numbers = FALSE
    
    expr_log10_filtered = expr_log10[expr_log10[[data_input$rna_class]] %in% top_expressed$miRNA,]
    
    row_names = expr_log10_filtered[[data_input$rna_class]]
    expr_log10_filtered[[data_input$rna_class]] = c()
    
    # sort colnames
    expr_log10_filtered = expr_log10_filtered[, ..sample_order]
    
    if (grepl(" ", xticks_names$categories[[data_input$detection_group]]) == TRUE) {
      detection_group = strsplit(xticks_names$categories[[data_input$detection_group]], " ")[[1]][1]
    }
    legend_title = sprintf("Expr (%s, log10,\ndetection rate = %s%% per %s)", gsub("_norm", "", data_input$norm), gsub("p", "", data_input$detection_rate), detection_group)
    filename = sprintf("%s/top_%s_log10", output_folder, i)
    
    plot_heatmap(expr_log10_filtered, row_names, i, identifier_column, prop, scientific_numbers, legend_title, colors, xticks_names, plot_props, plot_props$image_width, plot_props$image_height, filename)
  }
}

#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])

