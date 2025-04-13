options(bitmapType='cairo')
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(cowplot))
#suppressPackageStartupMessages(library(ComplexHeatmap))
#suppressPackageStartupMessages(library(viridisLite))
#suppressPackageStartupMessages(library(proxy))
suppressPackageStartupMessages(library(Mfuzz))
suppressPackageStartupMessages(library(data.table))
#suppressPackageStartupMessages(library(clValid))
suppressPackageStartupMessages(library(pbapply))
#suppressPackageStartupMessages(library(fcvalid))
suppressPackageStartupMessages(library(gridBase))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(svglite))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(openxlsx))

#snakemake = readRDS("snakemake_mfuzz_ts_analysis_age_development_clusters.rds")

data.table::setDTthreads(1)
set.seed(snakemake@params$parameters_props$set_seed)

print("mfuzz_ts_analysis_age_development_clusters")


#---------------------------------- Functions ----------------------------------
line_plot_center_lines = function(cluster_result, feature, clusters, num_of_clusters, tissue_th, plot_category, xticks_names, data_input, plot_props, output_folder_path_fig, file_name) {
  
  font_size_correction = 2.845
  
  # get colours based on tissues of the feature contained in the given tissues
  feature_elem = grep(feature, names(cluster_result$cluster), value = TRUE)
  cluster_feature_info = cluster_result$cluster[feature_elem]
  
  mask = (cluster_feature_info %in% clusters)
  cluster_feature_info_filtered = cluster_feature_info[mask]
  
  names(cluster_feature_info_filtered) = str_split_fixed(names(cluster_feature_info_filtered), "_", 2)[,2]
  
  # get center lines
  mask = (rownames(cluster_result$centers) %in% clusters)
  cluster_result_centerlines_group = cluster_result$centers[mask,]
  
  cluster_counts = table(cluster_feature_info_filtered)
  single_occ = names(cluster_counts[cluster_counts == 1]) 
  for (cl in single_occ) {
    tissue = names(cluster_feature_info_filtered)[cluster_feature_info_filtered == cl]
    rownames(cluster_result_centerlines_group)[rownames(cluster_result_centerlines_group) == cl] = sprintf("%s_%s", rownames(cluster_result_centerlines_group)[rownames(cluster_result_centerlines_group) == cl], tissue)
  }
  multiple_occ = names(cluster_counts[cluster_counts > 1]) 
  for (cl in multiple_occ) {
    tissues = names(cluster_feature_info_filtered)[cluster_feature_info_filtered == cl]
    for (ti in tissues) {
      actual_row = cluster_result_centerlines_group[rownames(cluster_result_centerlines_group) == cl,]
      cluster_result_centerlines_group = rbind(cluster_result_centerlines_group, actual_row)
      rownames(cluster_result_centerlines_group)[nrow(cluster_result_centerlines_group)] = sprintf("%s_%s", cl, ti)
    }
    cluster_result_centerlines_group = cluster_result_centerlines_group[rownames(cluster_result_centerlines_group) != cl,]
  }
  
  cluster_center_lines = melt(cluster_result_centerlines_group)
  colnames(cluster_center_lines) = c("CLUSTER", plot_category[2], "values")
  
  cluster_center_lines$tissue = str_split_fixed(cluster_center_lines$CLUSTER, "_", 2)[,2]
  cluster_center_lines$cluster_number = str_split_fixed(cluster_center_lines$CLUSTER, "_", 2)[,1]
  
  cluster_center_lines$time = factor(cluster_center_lines[, plot_category[2]], levels = sort(unique(as.numeric(cluster_center_lines[, plot_category[2]]))))
  
  line_width_list = c()
  for (time in unique(cluster_center_lines$age)) {
    cluster_center_lines_time = cluster_center_lines[cluster_center_lines$age == time,]
    cluster_number_list = c()
    for (i in 1: nrow(cluster_center_lines_time)) {
      actual_row = cluster_center_lines_time[i,]
      if (!(actual_row$cluster_number %in% cluster_number_list)) {
        cluster_number_list = append(cluster_number_list, actual_row$cluster_number)
        line_width = unname(cluster_counts[names(cluster_counts) == actual_row$cluster_number])
      } 
      line_width_list = append(line_width_list, line_width)
      line_width = line_width - 1
      }
  }
  
  cluster_center_lines$line_width = as.factor(line_width_list)
  
  legend_name = sprintf("%s for", feature)
  mask = rownames(cluster_result$membership) %in% feature_elem
  mem_feature_elem = cluster_result$membership[mask,]
  mask = (colnames(mem_feature_elem) %in% clusters)
  mem_feature_elem_filtered = mem_feature_elem[,mask]
  
  legend_labels = c()
  for (i in 1:length(cluster_center_lines$tissue)) {
    ti = cluster_center_lines$tissue[i]
    cl = cluster_center_lines$cluster_number[i]
    row = sprintf("%s_%s", feature, ti)
    mem = mem_feature_elem_filtered[row, cl]
    mem = format(round(mem * 100, 1), nsmall = 1)
    label_str = sprintf("%s in Cl. %s (mem.: %s%%)", ti, cl, mem)
    legend_labels = append(legend_labels, label_str)
  }
  cluster_center_lines$legend_labels = legend_labels
  
  # labels = setNames(cluster_center_lines$legend_labels, cluster_center_lines$tissue)
  labels = list()
  for (i in 1:length(cluster_center_lines$legend_labels)) {
    key = cluster_center_lines$tissue[i]
    value = cluster_center_lines$legend_labels[i]
    labels[[key]] = value
  }
  
  mask = (names(unlist(colours[[plot_category[1]]][unique(annot[[plot_category[1]]])])[]) %in% unique(cluster_center_lines$tissue))
  colours_filtered =unlist(colours[[plot_category[1]]][unique(annot[[plot_category[1]]])])[mask]
  
  number_of_overlapping_elements = unique(cluster_center_lines$line_width)
  line_width = seq(0, by = 3, length.out = length(number_of_overlapping_elements))
  line_width[1] = 1
  
  legend_side = "bottom"
  p = ggplot(cluster_center_lines, aes(x = time, y = values, group = CLUSTER)) + 
    geom_line(aes(color=tissue, size=line_width), alpha = 1) +
    #geom_text_repel(aes(label = Var1), data = data_plot[data_plot$Var2 == sort(unique(data_plot$Var2))[1],], color = "black", size = 3) +
    #scale_x_continuous(breaks = sort(unique(data_plot$time))) +#, labels = 1:10) 
    scale_size_manual(values=line_width, guide = "none") +
    scale_color_manual(values=colours_filtered, labels = labels) +
    xlab(sprintf("%s", xticks_names$categories[[sprintf("%s_capital", plot_category[2])]])) + 
    ylab(sprintf("Standardised expression\nvalues (%s)", str_split_fixed(data_input$norm, "_", 2)[1])) +
    scale_x_discrete(expand = c(0, 0.1)) +
    #geom_text(x = xpos, y = ypos, label = num_of_lines, hjust = 1, vjust = 1, family = plots_props$font_family, size = plots_props$font_size/font_size_correction) +
    theme_classic() +
    theme(legend.position=legend_side, 
          text = element_text(family = plot_props$font_family, size = plot_props$font_size),
          axis.text = element_text(family = plot_props$font_family, size = plot_props$font_size),
          axis.title = element_text(family = plot_props$font_family, size = plot_props$font_size),
          plot.title = element_text(family = plot_props$font_family, size = plot_props$font_size_header),
          legend.text = element_text(family = plot_props$font_family, size = 6),
          legend.title = element_text(family = plot_props$font_family, size = plot_props$font_size),
          legend.spacing.x = unit(0, "cm"),
          legend.spacing.y = unit(0, 'cm'),
          legend.key.size = unit(0.3, "cm"),
          #legend.box.margin = margin(-0.8,0,0,-4, "cm"),
          legend.margin = margin(-0.25, 0.25, 1, -1, "cm"), #t = -10, r = 5, b = -5, 
          plot.margin = margin(0.25, 0.25, -1, 0.25, "cm") #top, right, bottom, and left 
    ) +
  guides(color = guide_legend(title=legend_name, nrow = 4, title.position = "top", byrow=TRUE))
  #p
  
  # save bar plot
  width = plot_props$image_width
  height = plot_props$image_height
  ggsave(sprintf("%s/%s_legend_side=%s_%sx%s.png", output_folder_path_fig, file_name, legend_side, width, height), p, dpi=plot_props$dpi, width=width, height=height, units = plot_props$image_units)
  ggsave(sprintf("%s/%s_legend_side=%s_%sx%s.svg", output_folder_path_fig, file_name, legend_side, width, height), p, width=width, height=height, units =plot_props$image_units, limitsize = FALSE, scale = 1)
  
  width = 2/3 * plot_props$image_width
  height = plot_props$image_height
  ggsave(sprintf("%s/%s_legend_side=%s_%sx%s.png", output_folder_path_fig, file_name, legend_side, width, height), p, dpi=plot_props$dpi, width=width, height=height, units = plot_props$image_units)
  ggsave(sprintf("%s/%s_legend_side=%s_%sx%s.svg", output_folder_path_fig, file_name, legend_side, width, height), p, width=width, height=height, units =plot_props$image_units, limitsize = FALSE, scale = 1)
  
  legend_side = "right"
  p = ggplot(cluster_center_lines, aes(x = time, y = values, group = CLUSTER)) + 
    geom_line(aes(color=tissue, size=line_width), alpha = 1) +
    #geom_text_repel(aes(label = Var1), data = data_plot[data_plot$Var2 == sort(unique(data_plot$Var2))[1],], color = "black", size = 3) +
    #scale_x_continuous(breaks = sort(unique(data_plot$time))) +#, labels = 1:10) 
    scale_size_manual(values=line_width, guide = "none") +
    scale_color_manual(values=colours_filtered, labels = labels) +
    xlab(sprintf("%s", xticks_names$categories[[sprintf("%s_capital", plot_category[2])]])) + 
    ylab(sprintf("Standardised expression\nvalues (%s)", str_split_fixed(data_input$norm, "_", 2)[1])) +
    scale_x_discrete(expand = c(0, 0.1)) +
    #geom_text(x = xpos, y = ypos, label = num_of_lines, hjust = 1, vjust = 1, family = plots_props$font_family, size = plots_props$font_size/font_size_correction) +
    theme_classic() +
    theme(legend.position=legend_side, 
          text = element_text(family = plot_props$font_family, size = plot_props$font_size),
          axis.text = element_text(family = plot_props$font_family, size = plot_props$font_size),
          axis.title = element_text(family = plot_props$font_family, size = plot_props$font_size),
          plot.title = element_text(family = plot_props$font_family, size = plot_props$font_size_header),
          legend.text = element_text(family = plot_props$font_family, size = 6),
          legend.title = element_text(family = plot_props$font_family, size = plot_props$font_size),
          #legend.spacing.x = unit(0, "cm"),
          #legend.spacing.y = unit(0, 'cm'),
          legend.key.size = unit(0.3, "cm"),
          #legend.box.margin = margin(-0.8,0,0,-4, "cm"),
          legend.margin = margin(0.75, 0.25, 0.25, 0.25, "cm"), #t = -10, r = 5, b = -5, 
          #plot.margin = margin(0.75, 0.25, 0.1, 0.25, "cm") #top, right, bottom, and left 
    ) +
    guides(color = guide_legend(title=legend_name, byrow=TRUE))
  #p
  
  # save bar plot
  width = plot_props$image_width
  height = 2/3 * plot_props$image_height
  ggsave(sprintf("%s/%s_legend_side=%s_%sx%s.png", output_folder_path_fig, file_name, legend_side, width, height), p, dpi=plot_props$dpi, width=width, height=height, units = plot_props$image_units)
  ggsave(sprintf("%s/%s_legend_side=%s_%sx%s.svg", output_folder_path_fig, file_name, legend_side, width, height), p, width=width, height=height, units =plot_props$image_units, limitsize = FALSE, scale = 1)

  width = plot_props$image_width
  height = plot_props$image_height
  ggsave(sprintf("%s/%s_legend_side=%s_%sx%s.png", output_folder_path_fig, file_name, legend_side, width, height), p, dpi=plot_props$dpi, width=width, height=height, units = plot_props$image_units)
  ggsave(sprintf("%s/%s_legend_side=%s_%sx%s.svg", output_folder_path_fig, file_name, legend_side, width, height), p, width=width, height=height, units =plot_props$image_units, limitsize = FALSE, scale = 1)
}


#---------------------------------- Load data ----------------------------------
data_input = snakemake@params$data_input
cluster_props = snakemake@params$cluster_props
plot_category = snakemake@params$plot_category
diff_exp_log = snakemake@params$diff_log
thresholds = snakemake@params$thresholds
cluster_groups = snakemake@params$cluster_groups
colours = snakemake@params$colors
plot_props = snakemake@params$plots_props
results_folder = snakemake@params$results_folder
xticks_names = snakemake@params$xticks_names

# save rdata
#saveRDS(snakemake, file = "snakemake_mfuzz_ts_analysis_age_development_clusters.rds")

expr = fread(sprintf("%s_%s/%s_%s_quantification_%s.detection_rate_%s_per_%s.csv", data_input$raw_data_folder, data_input$data_sub_set, data_input$rna_class, data_input$feature_filtering, data_input$norm, data_input$detection_rate, data_input$detection_group), sep='\t', header=T)
annot = fread(sprintf("%s/annotation_%s_%s.csv", data_input$annotation_data_folder, data_input$data_sub_set, data_input$feature_filtering), sep='\t', colClasses=c(ID="character")) 
diff_exp = fread(sprintf("%s/%s_%s/results_%s/matrices/diff_exp/diff_exp%s.csv", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, diff_exp_log), sep='\t') 

expr_rownames = expr[[data_input$rna_class]]
expr[[data_input$rna_class]] = c()

input_folder_path_tab = sprintf("%s/%s_%s/results_%s/matrices/mfuzz_clustering/k=%s", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, cluster_props$num_of_clusters)

cluster_result = readRDS(file.path(snakemake@input$cluster_results_folder_path, sprintf("cluster_result_k=%s.rds", cluster_props$num_of_clusters)))
acore_dt = fread(file.path(snakemake@input$cluster_results_folder_path, sprintf("cluster_members_k=%s.csv", cluster_props$num_of_clusters)))

output_folder_path_fig = sprintf("%s/%s_%s/results_%s/figures/mfuzz_clustering/k=%s", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, cluster_props$num_of_clusters)
dir.create(output_folder_path_fig, recursive=TRUE)
#output_folder_path_tab = sprintf("%s/%s_%s/results_%s/matrices/mfuzz_clustering/k=%s", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, cluster_props$num_of_clusters)
#dir.create(output_folder_path_tab, recursive=TRUE)

  
#------------------------------------ Script ----------------------------------- 
# get medians from the diff exp table
median_table_df = data.frame()
for (tissue in unique(annot[[plot_category[1]]])){

  #print(tissue)
  
  median_table = data.table(feature=sprintf("%s_%s", diff_exp[[data_input$rna_class]], tissue))
  for (tp in paste(sort(unique(annot[[plot_category[2]]])))) {
    median_table[[tp]] = 2^(as.numeric(diff_exp[[sprintf("median_%s__%s_%s=%s", plot_category[2], tp, plot_category[1], tissue)]]))
  }
  
  median_table_df = rbind(median_table_df, data.frame(median_table, check.names = FALSE))
}


for (group in cluster_groups) {
  feature = group$feature
  clusters = group$cluster

  # make plot showing the centerlines for the given clusters 
  output_folder_path_fig_all = sprintf("%s/all/line_plots/age_development", output_folder_path_fig)
  dir.create(output_folder_path_fig_all, recursive=TRUE)
  file_name = sprintf("%s_center_lines_for_clusters_%s", feature, paste(clusters, collapse = "_"))
  line_plot_center_lines(cluster_result, feature, clusters, cluster_props$num_of_clusters, th, plot_category, xticks_names, data_input, plot_props, output_folder_path_fig_all, file_name)
  
  # make plot showing the medians per timepoint for the given tissue and feature shape indicating the cluster membership
  file_name = sprintf("%s_medians_per_timepoint_for_clusters_%s", feature, paste(clusters, collapse = "_"))
  #line_plot_medians_per_timepoint(cluster_result, cluster_info_tissue_df, cluster_props$num_of_clusters, th, plot_category, xticks_names, data_input, plot_props, output_folder_path_fig_all, file_name) 
  
  for (membership_plot_th in thresholds$membership_plot_ths) {
    # make plot showing the centerlines for the given clusters 
    output_folder_path_fig_th = sprintf("%s/%s/line_plots/age_development", output_folder_path_fig, membership_plot_th)
    dir.create(output_folder_path_fig_th, recursive=TRUE)
    file_name = sprintf("%s_center_lines_for_clusters_%s", feature, paste(clusters, collapse = "_"))
    line_plot_center_lines(cluster_result, feature, clusters, cluster_props$num_of_clusters, th, plot_category, xticks_names, data_input, plot_props, output_folder_path_fig_th, file_name)
  
    # make plot showing the medians per timepoint for the given tissue and feature shape indicating the cluster membership
    file_name = sprintf("%s_medians_per_timepoint_for_clusters_%s", feature, paste(clusters, collapse = "_"))
    #line_plot_medians_per_timepoint(cluster_result, cluster_info_tissue_df, cluster_props$num_of_clusters, th, plot_category, xticks_names, data_input, plot_props, output_folder_path_fig_th, file_name) 
  }
}


#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])
