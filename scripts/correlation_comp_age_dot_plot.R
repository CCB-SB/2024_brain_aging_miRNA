suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(svglite))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(openxlsx))

#snakemake = readRDS("snakemake_correlation_comp_age_dot_plot.rds")

set.seed(snakemake@params$parameters_props$set_seed)

print("correlation_comp_age_dot_plot")


#---------------------------------- Functions ----------------------------------
plot_heatmap_top_features_dots = function(tbl, row_names, tissue, colors, plots_props, file_name) {
  tbl = as.data.frame(tbl)
  rownames(tbl) = row_names
  
  cell_fun = function(j, i, x, y, w, h, col) {
    radius = 1
    if(tbl[i, j] == 1) {
      ti = colnames(tbl)[j]
      grid.circle(x, y, unit(radius, "mm"), gp = gpar(fill = colors[[tissue]][[ti]], col = colors[[tissue]][[ti]]))
    } else if(tbl[i, j] == 0) {
      grid.circle(x, y, unit(radius, "mm"), gp = gpar(fill = "white", col = "white"))
    }
  }
  
  p = ComplexHeatmap::Heatmap(tbl, rect_gp = gpar(type = "none"), cell_fun = cell_fun,
                              #rect_gp = gpar(col = "white", lwd = .5),
                              # height = unit(0.3, "cm") * nrow(t(df)), width = unit(0.3, "cm") * ncol(t(df)),
                              #height = unit(plots_props$image_height  - 2, "cm"), width = unit(plots_props$image_width - 4, "cm"),
                              height = unit(2.5, "cm"), width = unit(5, "cm"),
                              show_column_names = FALSE,
                              show_row_names = TRUE,
                              row_names_side = "left", row_names_gp = gpar(fontsize = plots_props$font_size, fontfamily=plots_props$font_family), 
                              #column_names_side = "bottom", column_names_gp = gpar(fontsize = plots_props$font_size, fontfamily=plots_props$font_family),
                              #cluster_columns = function(m) hclust(dist(m), method = "complete"), show_column_dend = TRUE, clustering_distance_columns = "euclidean", 
                              #cluster_rows = function(m) hclust(dist(m), method = "complete"), show_row_dend = TRUE, clustering_distance_rows = "euclidean", 
                              column_dend_reorder = FALSE, cluster_columns = FALSE,
                              row_dend_reorder = FALSE, cluster_rows = FALSE,
                              show_heatmap_legend = FALSE,
  )
  
  #legend_ticks = seq(from = -maximal_num_tissues, to = maximal_num_tissues, length.out = 5)
  #legend_tick_labels = round(legend_ticks, digits = 2)
  
  #legend = ComplexHeatmap::Legend(col_fun = col_fun, title = color_bar_name, direction="horizontal",
  #                                at = legend_ticks, labels = legend_tick_labels)
  
  png(sprintf("%s.png", file_name), width = plots_props$image_width, height = plots_props$image_height / 2, units = plots_props$image_units, res = plots_props$dpi)
  draw(p, padding = unit(c(0.25, 0.25, 0.25, 0.25), "cm"))  # bottom, left, top and right margins
  #draw(legend, x = unit(.5, "cm"), y = unit(.1, "cm"), just = c("left", "bottom"))
  dev.off()
  svglite(sprintf("%s.svg", file_name), width = plots_props$image_width / 2.54, height = (plots_props$image_height / 2) / 2.54)
  draw(p, padding = unit(c(0.25, 0.25, 0.25, 0.25), "cm"))  # bottom, left, top and right margins
  #draw(legend, x = unit(.5, "cm"), y = unit(.1, "cm"), just = c("left", "bottom"))
  dev.off()
}

plot_heatmap_top_features_dots_merge = function(tbl, row_names, up_down_length, tissue, colors, plots_props, file_name) {
  tbl = as.matrix(tbl)
  rownames(tbl) = row_names
  
  # vector with 1 for up and -1 for down
  row_split = c(rep(1, up_down_length[1]), rep(2, up_down_length[2]))
  
  cell_fun = function(j, i, x, y, w, h, col) {
    radius = 1
    if(tbl[i, j] == 1) {
      ti = colnames(tbl)[j]
      grid.circle(x, y, unit(radius, "mm"), gp = gpar(fill = colors[[tissue]][[ti]], col = colors[[tissue]][[ti]]))
    } else if(tbl[i, j] == 0) {
      grid.circle(x, y, unit(radius, "mm"), gp = gpar(fill = "white", col = "white"))
    }
  }
  
  # decide which labels to show on the left side of the dot plot
  # if up_down_length[1] == 0 -> no positives
  # if up_down_length[2] == 0 -> no negatives
  if (up_down_length[1] == 0) {
    labels = c("Neg.")
  } else if (up_down_length[2] == 0) {
    labels = c("Pos.")
  } else {
    labels = c("Pos.", "Neg.")
  }

  left_annotation = ComplexHeatmap::rowAnnotation(foo = anno_block(gp = gpar(fill = "white"),
                                                                   labels = labels,
                                                                   labels_gp = gpar(col = "black", fontsize = 10)))
  
  p = ComplexHeatmap::Heatmap(tbl, rect_gp = gpar(type = "none"), cell_fun = cell_fun,
                              #rect_gp = gpar(col = "white", lwd = .5),
                              # height = unit(0.3, "cm") * nrow(t(df)), width = unit(0.3, "cm") * ncol(t(df)),
                              #height = unit(plots_props$image_height  - 2, "cm"), width = unit(plots_props$image_width - 4, "cm"),
                              height = unit(2.5, "cm"), width = unit(5, "cm"),
                              show_column_names = FALSE,
                              show_row_names = TRUE,
                              row_names_side = "right", row_names_gp = gpar(fontsize = plots_props$font_size, fontfamily=plots_props$font_family), 
                              #column_names_side = "bottom", column_names_gp = gpar(fontsize = plots_props$font_size, fontfamily=plots_props$font_family),
                              #cluster_columns = function(m) hclust(dist(m), method = "complete"), show_column_dend = TRUE, clustering_distance_columns = "euclidean", 
                              #cluster_rows = function(m) hclust(dist(m), method = "complete"), show_row_dend = TRUE, clustering_distance_rows = "euclidean", 
                              column_dend_reorder = FALSE, cluster_columns = FALSE,
                              row_dend_reorder = FALSE, cluster_rows = FALSE,
                              show_heatmap_legend = FALSE,
                              row_split = row_split,
                              row_title = NULL,
                              left_annotation = left_annotation
  )
  
  #legend_ticks = seq(from = -maximal_num_tissues, to = maximal_num_tissues, length.out = 5)
  #legend_tick_labels = round(legend_ticks, digits = 2)
  
  #legend = ComplexHeatmap::Legend(col_fun = col_fun, title = color_bar_name, direction="horizontal",
  #                                at = legend_ticks, labels = legend_tick_labels)
  
  png(sprintf("%s.png", file_name), width = plots_props$image_width, height = plots_props$image_height / 2, units = plots_props$image_units, res = plots_props$dpi)
  draw(p, padding = unit(c(0.25, 0.25, 0.25, 0.25), "cm"))  # bottom, left, top and right margins
  #draw(legend, x = unit(.5, "cm"), y = unit(.1, "cm"), just = c("left", "bottom"))
  dev.off()
  svglite(sprintf("%s.svg", file_name), width = plots_props$image_width / 2.54, height = (plots_props$image_height / 2) / 2.54)
  draw(p, padding = unit(c(0.25, 0.25, 0.25, 0.25), "cm"))  # bottom, left, top and right margins
  #draw(legend, x = unit(.5, "cm"), y = unit(.1, "cm"), just = c("left", "bottom"))
  dev.off()
}


#---------------------------------- Load data ---------------------------------- 
data_input = snakemake@params$data_input
#tissues = snakemake@params$tissues 
method = snakemake@params$method
props = snakemake@params$props
plots_props = snakemake@params$plots_props
colors = snakemake@params$colors
results_folder = snakemake@params$results_folder

corr_list = snakemake@params$thresholds$corr_th
number_of_corr_tissues_thresholds = snakemake@params$thresholds$num_corr_tissues_th

# save rdata
#saveRDS(snakemake, file = "snakemake_correlation_comp_age_dot_plot.rds")

output_folder_fig = sprintf("%s/%s_%s/results_%s/figures/corr_plots/heatmap_complex/global_%ss", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, data_input$rna_class)
dir.create(output_folder_fig, recursive=TRUE)
input_folder_tab = sprintf("%s/%s_%s/results_%s/matrices/aggregation_corr", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set)


#------------------------------------ Script ----------------------------------- 
for (m in method) {
  for (prop in props) {
    
    tissue = prop[1]
    time = prop[2]
    for (th in corr_list) {
      for (num_corr_tissues in number_of_corr_tissues_thresholds) {
        # load data
        file_name = sprintf("corr_method=%s_global_sig_pos_for_%s_%s_%s_corr_th=%s", m, num_corr_tissues, tissue, data_input$rna_class, th)
        file_path = sprintf("%s/%s.csv", input_folder_tab, file_name)
        if (file.exists(file_path)) {
          global_pos = as.data.frame(fread(file_path, sep = "\t"))
          # check if there are NA in the table, then it was empty
          if (any(is.na(global_pos$V1))) {
            global_pos = data.frame()
          } else {
            global_pos[[data_input$rna_class]] = global_pos$V1
            global_pos$V1 = c()
            # reverse colmnn order such that features is in front
            global_pos = global_pos[, rev(colnames(global_pos))]
          }
        } else {
          global_pos = data.frame(column1 = numeric())  
        }
        
        file_name = sprintf("corr_method=%s_global_sig_neg_for_%s_%s_%s_corr_th=%s", m, num_corr_tissues, tissue, data_input$rna_class, th)
        file_path = sprintf("%s/%s.csv", input_folder_tab, file_name)
        if (file.exists(file_path)) {
          global_neg = as.data.frame(fread(file_path, sep = "\t"))
          # check if there are NA in the table, then it was empty
          if (any(is.na(global_neg$V1))) {
            global_neg = data.frame()
          } else {
            global_neg[[data_input$rna_class]] = global_neg$V1
            global_neg$V1 = c()
            # reverse colmnn order such that features is in front
            global_neg = global_neg[, rev(colnames(global_neg))]
          }
        } else {
          global_neg = data.frame(column1 = numeric())  
        }
        
        # check if up and down re empty
        if ((nrow(global_pos) == 0) & (nrow(global_neg) == 0)) {
          print(sprintf("pos and neg empty for method=%s, tissue=%s, time=%s, corr_th=%s num_corr_tissues=%s", m, tissue, time, th, num_corr_tissues))
          next
        }
        
        # make plots for each
        # up
        # sort by rowsum
        if (nrow(global_pos) <= 1) {
          global_pos_sorted = global_pos
        } else {
          global_pos_sorted = global_pos[order(rowSums(global_pos[, -1]), decreasing=T),]
        }
        # down
        # sort by rowsum
        if (nrow(global_neg) <= 1) {
          global_neg_sorted = global_neg
        } else {
          global_neg_sorted = global_neg[order(rowSums(global_neg[, -1]), decreasing=F),]
        }
        
        # merge up and down
        uni = union(colnames(global_pos_sorted), colnames(global_neg_sorted))
        inter = intersect(colnames(global_pos_sorted), colnames(global_neg_sorted))
        complement = setdiff(uni, inter)

        #print(uni)
        #print(inter)
        #print(complement)
        
        # columns missing in the respective tables
        missing_pos = complement[!(complement %in% colnames(global_pos_sorted))]
        missing_neg = complement[!(complement %in% colnames(global_neg_sorted))]
        
        # add columns
        # some dirty hack because R data.frames cannot set column values on an empty df
        if (nrow(global_pos_sorted) != 0) {
          # normal case
          global_pos_sorted[, missing_pos] = 0
        } else {
          # case it is empty
          # add dummy row
          global_pos_sorted[1, ] = 0
          # add columns
          global_pos_sorted[, missing_pos] = 0
          # remove dummy row
          global_pos_sorted = global_pos_sorted[c(), , drop=FALSE]
        }
        # same for neg
        if (nrow(global_neg_sorted) != 0) {
          # normal case
          global_neg_sorted[, missing_neg] = 0
        } else {
          # case it is empty
          # add dummy row
          global_neg_sorted[1, ] = 0
          # add columns
          global_neg_sorted[, missing_neg] = 0
          # remove dummy row
          global_neg_sorted = global_neg_sorted[c(), , drop=FALSE]
        }        
        # change to column order of down_reg
        col_order = colnames(global_pos_sorted)
        # fix all types to data.frame
        global_neg_sorted = global_neg_sorted[, col_order]
        
        # merge to one table
        global_merge_sorted = rbind(global_pos_sorted, global_neg_sorted)
        up_down_length = c(dim(global_pos_sorted)[1], dim(global_neg_sorted)[1])
        #print("Pos:")
        #print(global_pos_sorted)
        #print("Neg:")
        #print(global_neg_sorted)
        #print("Merge:")
        #print(global_merge_sorted)
        
        # create plots
        #up
        if (nrow(global_pos_sorted) != 0) {
          file_name = sprintf("corr_method=%s_global_%s_over_%s_corr_th=%s_for_%s_%s_sig_pos_dots", m, data_input$rna_class, time, th, num_corr_tissues, tissue)
          file_name_full = sprintf("%s/%s", output_folder_fig, file_name)
          plot_heatmap_top_features_dots(global_pos_sorted[, -1], global_pos_sorted[[data_input$rna_class]], tissue, colors, plots_props, file_name_full) 
        } else {
          print(sprintf("global_pos_sorted empty for corr_th=%s num_corr_tissues=%s", th, num_corr_tissues))
        }
        
        #down
        if (nrow(global_neg_sorted) != 0) {
          file_name = sprintf("corr_method=%s_global_%s_over_%s_corr_th=%s_for_%s_%s_sig_neg_dots", m, data_input$rna_class, time, th, num_corr_tissues, tissue)
          file_name_full = sprintf("%s/%s", output_folder_fig, file_name)
          plot_heatmap_top_features_dots(global_neg_sorted[, -1], global_neg_sorted[[data_input$rna_class]], tissue, colors, plots_props, file_name_full) 
        } else {
          print(sprintf("global_neg_sorted empty for corr_th=%s num_corr_tissues=%s", th, num_corr_tissues))
        }
        
        #merge
        #the empty case was catched after loading the matrices
        file_name = sprintf("corr_method=%s_global_%s_over_%s_corr_th=%s_for_%s_%s_sig_dots", m, data_input$rna_class, time, th, num_corr_tissues, tissue)
        file_name_full = sprintf("%s/%s", output_folder_fig, file_name)
        plot_heatmap_top_features_dots_merge(global_merge_sorted[, -1], global_merge_sorted[[data_input$rna_class]], up_down_length, tissue, colors, plots_props, file_name_full)
        print("---------------------------------------------------------")
      }
    }
  }
}

#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])

