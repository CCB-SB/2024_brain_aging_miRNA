suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(svglite))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(Hmisc))
#suppressPackageStartupMessages(library(org.Mm.eg.db))

#snakemake = readRDS("snakemake_correlation_mirna_mrna_heatmap_plot.rds")

set.seed(snakemake@params$parameters_props$set_seed)

print("correlation_mirna_mrna_heatmap_plot")


#---------------------------------- Functions ----------------------------------
plot_heatmap_all = function(df, pvalues_df, plots_props, m, time, data_input, colors, filename, output_folder_fig, output_folder_tab) {
  df = t(df)
  df[abs(df) <= th] = 0
  df[df < -th] = -0.5
  df[df > th] = 0.5
  
  color_bar_name = sprintf("%s correlation (%s)", time, m) 
  
  #cmap = colorRamp2(c(-1, 0, 1), hcl_palette = "Green-Brown", reverse = TRUE)
  #col_val = 0.75
  #col_fun = c(cmap(-col_val), cmap(0), cmap(col_val))
  #col_legend = c(cmap(-col_val), cmap(col_val))
  col_fun = c(colors$direction$neg, colors$direction$neg, colors$direction$pos)
  col_legend = c(colors$direction$neg, colors$direction$pos)
  names(col_legend) = c(sprintf("negative (R ≤ -%s)", th), sprintf("positive (R ≥ %s)", th))
  
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (is.na(t(pvalues_df)[i, j])) {
    }
    else if (t(pvalues_df)[i, j] < 0.05){ 
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
  
  colnames(df) = row_names
  #print(annotation_row)
  #print(annotation_colors)
  column_ha = rowAnnotation(df = annotation_row, col = annotation_colors, gp = gpar(row = "white"), simple_anno_size = unit(0.2, "cm"), 
                            annotation_name_side = "top", annotation_name_gp = gpar(fontsize = plots_props$font_size, fontface = "bold", fontfamily=plots_props$font_family),
                            show_legend=FALSE)
  
  p = ComplexHeatmap::Heatmap(df, col = col_fun, cell_fun = cell_fun,
                              #rect_gp = gpar(col = "white", lwd = .5),
                              # height = unit(0.3, "cm") * nrow(t(df)), width = unit(0.3, "cm") * ncol(t(df)),
                              #height = unit(plots_props$image_height  - 2, "cm"), width = unit(plots_props$image_width - 4, "cm"),
                              height = unit(3, "cm"), width = unit(6.25*3, "cm"),
                              show_column_names = FALSE,
                              row_names_side = "left", row_names_gp = gpar(fontsize = 6, fontfamily=plots_props$font_family), 
                              #column_names_side = "bottom", column_names_gp = gpar(fontsize = plots_props$font_size, fontfamily=plots_props$font_family),
                              cluster_columns = function(m) hclust(dist(m), method = "complete"), show_column_dend = TRUE, clustering_distance_columns = "euclidean", 
                              #cluster_rows = function(m) hclust(dist(m), method = "complete"), show_row_dend = TRUE, clustering_distance_rows = "euclidean", 
                              column_dend_reorder = TRUE, #cluster_columns = TRUE,
                              row_dend_reorder = FALSE, cluster_rows = FALSE,
                              show_heatmap_legend = FALSE,
                              left_annotation = column_ha
  )
  #legend_ticks = seq(from = min(t(df)), to = max(t(df)), length.out = 5)
  #legend_tick_labels = round(legend_ticks, digits = 2)
  
  #legend = ComplexHeatmap::Legend(title = color_bar_name, col_fun = col_fun,
  #                                at = legend_tick_labels, direction="horizontal", #col_fun = circlize::colorRamp2(legend_ticks, viridis(5))
  #                                title_gp = gpar(fontsize = plots_props$font_size, fontface = "bold", fontfamily=plots_props$font_family), labels_gp = gpar(fontsize = plots_props$font_size, fontfamily=plots_props$font_family),
  #                                legend_width=unit(4, "cm"))
  
  legend = ComplexHeatmap::Legend(title = color_bar_name, labels = names(col_legend), legend_gp = gpar(fill = col_legend),
                                  direction="horizontal",
                                  title_gp = gpar(fontsize = plots_props$font_size, fontface = "bold", fontfamily=plots_props$font_family), labels_gp = gpar(fontsize = plots_props$font_size, fontfamily=plots_props$font_family),
                                  legend_width=unit(4, "cm"))
  
  #p = pheatmap(t(df), color = viridis(100), border_color = "white", show_colnames = TRUE, show_rownames = TRUE,
  #             fontsize = plots_props$font_size, fontsize_row = 6, fontsize_col = plots_props$font_size, 
  #             cellwidth = 8, cellheight = 8,
  #             cluster_cols = FALSE, clustering_distance_cols = "euclidean", cluster_rows = TRUE, clustering_distance_rows = "euclidean", treeheight_row = 0
  #)
  
  png(sprintf("%s_corr_th=%s.png", filename, th), width = plots_props$image_width*2.5, height = plots_props$image_height*1.25, units = plots_props$image_units, res = plots_props$dpi)
  draw(p, padding = unit(c(0.25, 0.25, -1.5, 0.25), "cm"))  # bottom, left, top and right margins
  draw(legend, x = unit(.5, "cm"), y = unit(.1, "cm"), just = c("left", "bottom"))
  dev.off()
  svglite(sprintf("%s_corr_th=%s.svg", filename, th), width = plots_props$image_width*2.5 / 2.54, height = plots_props$image_height*1.25 / 2.54)
  draw(p, padding = unit(c(0.25, 0.25, -1.5, 0.25), "cm"))  # bottom, left, top and right margins
  draw(legend, x = unit(.5, "cm"), y = unit(.1, "cm"), just = c("left", "bottom"))
  dev.off()
}

plot_heatmap_th = function(df, pvalues_df, plots_props, m, colors, th, filename, output_folder_fig, output_folder_tab) {

  color_bar_name = sprintf("correlation (%s)", m) 
    
  if (i == 1) {
    col_fun = c(colors$direction$light, colors$direction$pos)
    col_legend = c(colors$direction$pos)
    names(col_legend) = c(sprintf("positive (R ≥ %s)", th))
  } else {
    col_fun = c(colors$direction$neg, colors$direction$light)
    col_legend = c(colors$direction$neg)
    names(col_legend) = c(sprintf("negative (R ≤ -%s)", th))
  }
  
  df_mod = as.data.frame(df_mod)
  rownames(df_mod) = row_names_mod
  num_cols = ncol(df_mod)
  num_rows = nrow(df_mod)
  
  p = ComplexHeatmap::Heatmap(df_mod, col = col_fun,
                              #cell_fun = cell_fun,
                              #rect_gp = gpar(col = "white", lwd = .5),
                              #height = unit(0.3, "cm") * nrow(t(df)), width = unit(0.3, "cm") * ncol(t(df)),
                              #height = unit(plots_props$image_height  - 2, "cm"), width = unit(plots_props$image_width - 4, "cm"),
                              height = unit(0.3*num_cols, "cm"), width = unit(0.75*num_rows, "cm"),
                              #show_column_names = FALSE,
                              show_row_names = TRUE,
                              row_names_side = "left", row_names_gp = gpar(fontsize = 6, fontfamily=plots_props$font_family), 
                              column_names_side = "bottom", column_names_gp = gpar(fontsize = 6, fontfamily=plots_props$font_family),
                              cluster_columns = function(m) hclust(dist(m), method = "complete"), show_column_dend = FALSE, clustering_distance_columns = "euclidean", 
                              cluster_rows = function(m) hclust(dist(m), method = "complete"), show_row_dend = FALSE, clustering_distance_rows = "euclidean", 
                              column_dend_reorder = TRUE, #cluster_columns = FALSE,
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
  
  legend = ComplexHeatmap::Legend(title = color_bar_name, labels = names(col_legend), legend_gp = gpar(fill = col_legend),
                                  direction="horizontal",
                                  title_gp = gpar(fontsize = plots_props$font_size, fontface = "bold", fontfamily=plots_props$font_family), labels_gp = gpar(fontsize = plots_props$font_size, fontfamily=plots_props$font_family),
                                  legend_width=unit(4, "cm"))
  
  #p = pheatmap(t(df), color = viridis(100), border_color = "white", show_colnames = TRUE, show_rownames = TRUE,
  #             fontsize = plots_props$font_size, fontsize_row = 6, fontsize_col = plots_props$font_size, 
  #             cellwidth = 8, cellheight = 8,
  #             cluster_cols = FALSE, clustering_distance_cols = "euclidean", cluster_rows = TRUE, clustering_distance_rows = "euclidean", treeheight_row = 0
  #)
  p
  
  if (i == 1) {
    file_name_mod = sprintf("%s_pos", file_name)
  } else {
    file_name_mod = sprintf("%s_neg", file_name)
  }
  
  png(sprintf("%s/%s.png", output_folder_fig, file_name_mod), width = plots_props$image_width, height = plots_props$image_height, units = plots_props$image_units, res = plots_props$dpi)
  draw(p, padding = unit(c(0.25, 0.25, 0.25, 0.25), "cm"))  # bottom, left, top and right margins
  draw(legend, x = unit(.5, "cm"), y = unit(.1, "cm"), just = c("left", "bottom"))
  dev.off()
  svglite(sprintf("%s/%s.svg", output_folder_fig, file_name_mod), width = plots_props$image_width / 2.54, height = plots_props$image_height / 2.54)
  draw(p, padding = unit(c(0.25, 0.25, 0.25, 0.25), "cm"))  # bottom, left, top and right margins
  draw(legend, x = unit(.5, "cm"), y = unit(.1, "cm"), just = c("left", "bottom"))
  dev.off()
}


#---------------------------------- Load data ---------------------------------- 
data_input = snakemake@params$data_input
method = snakemake@params$method
adjustment = snakemake@params$adjustment
corr_th = snakemake@params$corr_th
plots_props = snakemake@params$plots_props
colors = snakemake@params$colors
results_folder = snakemake@params$results_folder

# mirtarbase = fread(snakemake@input$mirtarbase)

# save rdata
#saveRDS(snakemake, file = "snakemake_correlation_mirna_mrna_heatmap_plot.rds")


#------------------------------------ Script ----------------------------------- 
for (m in method) {
    
  corr_table = fread(sprintf("%s/%s_%s/results_%s/matrices/correlation_mirna_mrna/correlation_%s.csv", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, m), sep='\t', header=T)
  pvalues_adj_table = fread(sprintf("%s/%s_%s/results_%s/matrices/correlation_mirna_mrna/padj_%s.csv", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, m), sep='\t', header=T)
  
  output_folder_fig = sprintf("%s/%s_%s/results_%s/figures/corr_mirna_mrna_plots/heatmap_complex/%s", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, m)
  dir.create(output_folder_fig, recursive=TRUE)
  output_folder_tab = sprintf("%s/%s_%s/results_%s/matrices/correlation_mirna_mrna/%s", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, m)
  dir.create(output_folder_tab, recursive=TRUE)
  
  # all miRNAs
  #file_name = sprintf("corr_mirna_mrna_method=%s", m)
  #plot_heatmap_all(correlation_df, corr_pvalues_df, plots_props, m, time, data_input, colors, file_name, output_folder_fig, output_folder_tab)
  
  # only miRNAs and mRNAs with |corr| >= th
  for (th in corr_th) {
    print(th)
    
    row_names = df$rna
    if ((length(unique(df$rna == pvalues_df$rna)) == 1) & (unique(df$rna == pvalues_df$rna) == TRUE)) {
      print(unique(df$rna == pvalues_df$rna))
      
      df$rna = c()
      pvalues_df$rna = c()
      
      df_pos = df
      df_pos[df_pos < th] = 0
      df_pos[df_pos >= th] = 0.5
      
      df_neg = df
      df_neg[df_neg > -th] = 0
      df_neg[df_neg <= -th] = -0.5
      
      i = 2
      for (df_mod in list(df_neg)) { #df_pos, 
        
        if (i == 1) {
          print("pos")
        } else {
          print("neg")
        }
        
        tmp = colSums(df_mod != 0) > 0
        df_mod = df_mod[, ..tmp]
        pvalues_df_mod = pvalues_df[, ..tmp]
        tmp = colSums(is.na(df_mod)) != nrow(df_mod)
        df_mod = df_mod[, ..tmp]
        pvalues_df_mod = pvalues_df_mod[, ..tmp]
        
        #print(dim(df_mod))
        #print(dim(pvalues_df_mod))
        
        # filter for genes with info from mirTarBase
        
        
        # filter correlation
        filter_cols = colSums(df_mod) != 0
        filter_rows = rowSums(df_mod) != 0
        df_mod = df_mod[filter_rows, ..filter_cols]
        pvalues_df_mod = pvalues_df_mod[filter_rows, ..filter_cols]
        row_names_mod = row_names[filter_rows]
        
        #print("after correlation th filtering")
        #print(dim(df_mod))
        #print(dim(pvalues_df_mod))
        
        # filter significance
        sig_th = 0.01
        sig_mask = pvalues_df_mod < sig_th
        filter_cols = colSums(sig_mask) != 0
        filter_rows = rowSums(sig_mask) != 0
        df_mod = df_mod[filter_rows, ..filter_cols]
        pvalues_df_mod = pvalues_df_mod[filter_rows, ..filter_cols]
        row_names_mod = row_names_mod[filter_rows]
        
        #print("after significance filtering")
        #print(dim(df_mod))
        #print(dim(pvalues_df_mod))
        
        
        #df_label = matrix("|", nrow(df), ncol(df))
        #df_label[,tmp] = ""
      file_name = sprintf("corr_mirna_mrna_method=%s_miRNA_corr_th=%s", m, th)
      plot_heatmap_th(corr_table, pvalues_adj_table, plots_props, m, colors, th, file_name, output_folder_fig, output_folder_tab)
      }
    }
  }
}


#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])

