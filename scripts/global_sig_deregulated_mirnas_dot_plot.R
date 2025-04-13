suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(svglite))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(openxlsx))

#snakemake = readRDS("global_sig_deregulated_mirnas_dot_plot")

set.seed(snakemake@params$parameters_props$set_seed)

print("global_significant_deregulated")

#---------------------------------- Functions ----------------------------------
plot_heatmap_top_features_dots = function(tbl, row_names, tissue, colors, plots_props, file_name) {
  tbl = as.data.frame(tbl)
  rownames(tbl) = row_names
  
  cell_fun = function(j, i, x, y, w, h, col) {
    if(tbl[i, j] == 1) {
      ti = colnames(tbl)[j]
      grid.circle(x, y, unit(2, "mm"), gp = gpar(fill = colors[[tissue]][[ti]], col = colors[[tissue]][[ti]]))
    } else if(tbl[i, j] == 0) {
      grid.circle(x, y, unit(2, "mm"), gp = gpar(fill = "white", col = "white"))
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
  # if up_down_length[1] == 0 -> no up regulated
  # if up_down_length[2] == 0 -> no down regulated
  if (up_down_length[1] == 0) {
    labels = c("Down")
  } else if (up_down_length[2] == 0) {
    labels = c("Up")
  } else {
    labels = c("Up", "Down")
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
  
  # 6 = plots_props$image_height / 2
  # 8 = plots_props$image_height / 1.5
  
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
props = snakemake@params$props
plots_props = snakemake@params$plots_props
colors = snakemake@params$colors
number_of_deregulated_tissues_thresholds = snakemake@params$thresholds$num_dereg_tissues_th
number_of_deregulated_combinations_thresholds = snakemake@params$thresholds$num_dereg_combinations_th
ID = data_input$identifier_column
results_folder = snakemake@params$results_folder

# save rdata
#saveRDS(snakemake, file = "global_sig_deregulated_mirnas_dot_plot")

input_folder_tab = sprintf("%s/%s_%s/results_%s/matrices/aggregation_diff_exp", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set)

# create output folder
output_folder_fig = sprintf("%s/%s_%s/results_%s/figures/diff_exp/heatmap_complex/global_%ss", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, data_input$rna_class)
dir.create(output_folder_fig, recursive=TRUE)


#------------------------------------ Script -----------------------------------
#------------ sig ------------
for (prop in props) {
  
  tissue = prop[1]
  time = prop[2]
  
  for (num_dereg_combinations in number_of_deregulated_combinations_thresholds){
    for (num_dereg_tissues in number_of_deregulated_tissues_thresholds) {
      # load data
      file_name = sprintf("diff_exp_global_sig_up_%s_at_least_in_%s_%s_comps_for_%s_%s",data_input$rna_class, num_dereg_combinations, time, num_dereg_tissues, tissue)
      file_path = sprintf("%s/%s.csv", input_folder_tab, file_name)
      if (file.exists(file_path)) {
        global_up = as.data.frame(fread(file_path, sep = "\t"))
        # check if there are NA in the table, then it was empty
        if (any(is.na(global_up$V1))) {
          global_up = data.frame()
        } else {
          global_up$feature = global_up$V1
          global_up$V1 = c()
          # reverse colmnn order such that feature is in front
          global_up = global_up[, rev(colnames(global_up))]
        }
      } else {
        global_up = data.frame(column1 = numeric())  
      }

      file_name = sprintf("diff_exp_global_sig_down_%s_at_least_in_%s_%s_comps_for_%s_%s",data_input$rna_class, num_dereg_combinations, time, num_dereg_tissues, tissue)
      file_path = sprintf("%s/%s.csv", input_folder_tab, file_name)
      if (file.exists(file_path)) {
        global_down = as.data.frame(fread(file_path, sep = "\t"))
        # check if there are NA in the table, then it was empty
        if (any(is.na(global_down$V1))) {
          global_down = data.frame()
        } else {
          global_down$feature = global_down$V1
          global_down$V1 = c()
          # reverse colmnn order such that feature is in front
          global_down = global_down[, rev(colnames(global_down))]
        } 
      } else {
        global_down = data.frame(column1 = numeric())  
      }

      # check if up and down re empty
      if ((nrow(global_up) == 0) & ((nrow(global_down) == 0))) {
        print(sprintf("sig. up and sig. down empty for num_dereg_combinations=%s num_dereg_tissues=%s", num_dereg_combinations, num_dereg_tissues))
        next
      }
      
      # make plots for each
      # up
      # sort by rowsum
      if (dim(global_up)[1] <= 1) {
        global_up_sorted = global_up
      } else if (dim(global_up)[2] <= 2) {
        global_up_sorted = global_up[order(global_up[, -1], decreasing=F),]
      } else {
        global_up_sorted = global_up[order(rowSums(global_up[, -1]), decreasing=T),]
      }
      
      # down
      # sort by rowsum
      if (dim(global_down)[1] <= 1) {
        global_down_sorted = global_down
      } else if (dim(global_down)[2] <= 2) {
        global_down_sorted = global_down[order(global_down[, -1], decreasing=F),]
      } else {
        global_down_sorted = global_down[order(rowSums(global_down[, -1]), decreasing=F),]
      }
      
      # merge up and down
      uni = union(colnames(global_up_sorted), colnames(global_down_sorted))
      inter = intersect(colnames(global_up_sorted), colnames(global_down_sorted))
      complement = setdiff(uni, inter)
      
      # columns missing in the respective tables
      missing_up = complement[!(complement %in% colnames(global_up_sorted))]
      missing_down = complement[!(complement %in% colnames(global_down_sorted))]
      
      # add columns
      #####
      # original
      #global_up_sorted[, missing_up] = 0
      #global_down_sorted[, missing_down] = 0
      #####
      # some dirty hack because R data.frames cannot set column values on an empty df
      if (nrow(global_up_sorted) != 0) {
        # normal case
        global_up_sorted[, missing_up] = 0
      } else {
        # case it is empty
        # add dummy row
        global_up_sorted[1, ] = 0
        # add columns
        global_up_sorted[, missing_up] = 0
        # remove dummy row
        global_up_sorted = global_up_sorted[c(), , drop=FALSE]
      }
      # same for down
      if (nrow(global_down_sorted) != 0) {
        # normal case
        global_down_sorted[, missing_down] = 0
      } else {
        # case it is empty
        # add dummy row
        global_down_sorted[1, ] = 0
        # add columns
        global_down_sorted[, missing_down] = 0
        # remove dummy row
        global_down_sorted = global_down_sorted[c(), , drop=FALSE]
      }

      # change to column order of down_reg
      col_order = colnames(global_up_sorted)
      global_down_sorted = global_down_sorted[, col_order]
      
      # merge to one table
      global_merge_sorted = rbind(global_up_sorted, global_down_sorted)
      up_down_length = c(dim(global_up_sorted)[1], dim(global_down_sorted)[1])
      #print("Up:")
      #print(global_up_sorted)
      #print("Down:")
      #print(global_down_sorted)
      #print("Merge:")
      #print(global_merge_sorted)
      #print("---------------------------------------------------------")
      
      # create plots
      #up
      if (nrow(global_up_sorted) != 0) {
        file_name = sprintf("diff_exp_global_sig_up_%s_at_least_in_%s_%s_comps_for_%s_%s", data_input$rna_class, num_dereg_combinations, time, num_dereg_tissues, tissue)
        file_name_full = sprintf("%s/%s", output_folder_fig, file_name)
        plot_heatmap_top_features_dots(global_up_sorted[, -1, drop=FALSE], global_up_sorted$feature, tissue, colors, plots_props, file_name_full) 
      } else {
        print(sprintf("sig. global_up_sorted empty for num_dereg_combinations=%s num_dereg_tissues=%s", num_dereg_combinations, num_dereg_tissues))
      }
      
      #down
      if (nrow(global_down_sorted) != 0) {
        file_name = sprintf("diff_exp_global_sig_down_%s_at_least_in_%s_%s_comps_for_%s_%s", data_input$rna_class, num_dereg_combinations, time, num_dereg_tissues, tissue)
        file_name_full = sprintf("%s/%s", output_folder_fig, file_name)
        plot_heatmap_top_features_dots(global_down_sorted[, -1, drop=FALSE], global_down_sorted$feature, tissue, colors, plots_props, file_name_full) 
      } else {
        print(sprintf("sig. global_down_sorted empty for num_dereg_combinations=%s num_dereg_tissues=%s", num_dereg_combinations, num_dereg_tissues))
      }
      
      #merge
      #the empty case was catched after loading the matrices
      file_name = sprintf("diff_exp_global_sig_de_reg_%s_at_least_in_%s_%s_comps_for_%s_%s", data_input$rna_class, num_dereg_combinations, time, num_dereg_tissues, tissue)
      file_name_full = sprintf("%s/%s", output_folder_fig, file_name)
      plot_heatmap_top_features_dots_merge(global_merge_sorted[, -1, drop=FALSE], global_merge_sorted$feature, up_down_length, tissue, colors, plots_props, file_name_full)
    
    }
  }
}

#------------ not sig ------------
for (prop in props) {
  
  tissue = prop[1]
  time = prop[2]
  
  for (num_dereg_combinations in number_of_deregulated_combinations_thresholds){
    for (num_dereg_tissues in number_of_deregulated_tissues_thresholds) {     
      
      # load data
      file_name = sprintf("diff_exp_global_up_%s_at_least_in_%s_%s_comps_for_%s_%s", data_input$rna_class, num_dereg_combinations, time, num_dereg_tissues, tissue)
      global_up = as.data.frame(fread(sprintf("%s/%s.csv", input_folder_tab, file_name), sep = "\t"))
      # check if there are NA in the table, then it was empty
      if (any(is.na(global_up$V1))) {
        global_up = data.frame()
      } else {
        global_up$feature = global_up$V1
        global_up$V1 = c()
        # reverse colmnn order such that feature is in front
        global_up = global_up[, rev(colnames(global_up))]
      }
      
      file_name = sprintf("diff_exp_global_down_%s_at_least_in_%s_%s_comps_for_%s_%s", data_input$rna_class, num_dereg_combinations, time, num_dereg_tissues, tissue)
      global_down = as.data.frame(fread(sprintf("%s/%s.csv", input_folder_tab, file_name), sep = "\t"))
      # check if there are NA in the table, then it was empty
      if (any(is.na(global_down$V1))) {
        global_down = data.frame()
      } else {
        global_down$feature = global_down$V1
        global_down$V1 = c()
        # reverse colmnn order such that feature is in front
        global_down = global_down[, rev(colnames(global_down))]
      }
      
      # check if up and down re empty
      if ((nrow(global_up) == 0) & ((nrow(global_down) == 0))) {
        print(sprintf("up and down empty for num_dereg_combinations=%s num_dereg_tissues=%s", num_dereg_combinations, num_dereg_tissues))
        next
      }
      
      # make plots for each
      # up
      # sort by rowsum
      if (dim(global_up)[1] <= 1) {
        global_up_sorted = global_up
      } else {
        global_up_sorted = global_up[order(rowSums(global_up[, -1]), decreasing=T),]
      }
      
      # down
      # sort by rowsum
      if (dim(global_down)[1] <= 1) {
        global_down_sorted = global_down
      } else {
        global_down_sorted = global_down[order(rowSums(global_down[, -1]), decreasing=F),]
      }
      
      # merge up and down
      uni = union(colnames(global_up_sorted), colnames(global_down_sorted))
      inter = intersect(colnames(global_up_sorted), colnames(global_down_sorted))
      complement = setdiff(uni, inter)
      
      # columns missing in the respective tables
      missing_up = complement[!(complement %in% colnames(global_up_sorted))]
      missing_down = complement[!(complement %in% colnames(global_down_sorted))]
      
      # add columns
      #####
      # original
      #global_up_sorted[, missing_up] = 0
      #global_down_sorted[, missing_down] = 0
      #####
      # some dirty hack because R data.frames cannot set column values on an empty df
      if (nrow(global_up_sorted) != 0) {
        # normal case
        global_up_sorted[, missing_up] = 0
      } else {
        # case it is empty
        # add dummy row
        global_up_sorted[1, ] = 0
        # add columns
        global_up_sorted[, missing_up] = 0
        # remove dummy row
        global_up_sorted = global_up_sorted[c(), , drop=FALSE]
      }
      # same for down
      if (nrow(global_down_sorted) != 0) {
        # normal case
        global_down_sorted[, missing_down] = 0
      } else {
        # case it is empty
        # add dummy row
        global_down_sorted[1, ] = 0
        # add columns
        global_down_sorted[, missing_down] = 0
        # remove dummy row
        global_down_sorted = global_down_sorted[c(), , drop=FALSE]
      }
      
      # change to column order of down_reg
      col_order = colnames(global_up_sorted)
      global_down_sorted = global_down_sorted[, col_order]
      
      # merge to one table
      global_merge_sorted = rbind(global_up_sorted, global_down_sorted)
      up_down_length = c(dim(global_up_sorted)[1], dim(global_down_sorted)[1])
      #print("Up:")
      #print(global_up_sorted)
      #print("Down:")
      #print(global_down_sorted)
      #print("Merge:")
      #print(global_merge_sorted)
      #print("---------------------------------------------------------")
      
      # create plots
      #up
      if (nrow(global_up_sorted) != 0) {
        file_name = sprintf("diff_exp_global_up_%s_at_least_in_%s_%s_comps_for_%s_%s", data_input$rna_class, num_dereg_combinations, time, num_dereg_tissues, tissue)
        file_name_full = sprintf("%s/%s", output_folder_fig, file_name)
        plot_heatmap_top_features_dots(global_up_sorted[, -1, drop=FALSE], global_up_sorted$feature, tissue, colors, plots_props, file_name_full) 
      } else {
        print(sprintf("global_up_sorted empty for num_dereg_combinations=%s num_dereg_tissues=%s", num_dereg_combinations, num_dereg_tissues))
      }
      
      #down
      if (nrow(global_down_sorted) != 0) {
        file_name = sprintf("diff_exp_global_down_%s_at_least_in_%s_%s_comps_for_%s_%s", data_input$rna_class, num_dereg_combinations, time, num_dereg_tissues, tissue)
        file_name_full = sprintf("%s/%s", output_folder_fig, file_name)
        plot_heatmap_top_features_dots(global_down_sorted[, -1, drop=FALSE], global_down_sorted$feature, tissue, colors, plots_props, file_name_full) 
      } else {
        print(sprintf("global_down_sorted empty for num_dereg_combinations=%s num_dereg_tissues=%s", num_dereg_combinations, num_dereg_tissues))
      }
      
      #merge
      #the empty case was catched after loading the matrices
      file_name = sprintf("diff_exp_global_de_reg_%s_at_least_in_%s_%s_comps_for_%s_%s", data_input$rna_class, num_dereg_combinations, time, num_dereg_tissues, tissue)
      file_name_full = sprintf("%s/%s", output_folder_fig, file_name)
      plot_heatmap_top_features_dots_merge(global_merge_sorted[, -1, drop=FALSE], global_merge_sorted$feature, up_down_length, tissue, colors, plots_props, file_name_full)
    }
  }
}



#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])