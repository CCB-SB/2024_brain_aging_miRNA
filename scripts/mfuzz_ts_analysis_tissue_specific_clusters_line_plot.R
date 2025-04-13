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

source("scripts/mfuzz_ggplot.R")

#snakemake = readRDS("snakemake_mfuzz_ts_analysis_tissue_specific_clusters_line_plot.rds")

data.table::setDTthreads(1)
set.seed(snakemake@params$parameters_props$set_seed)

print("mfuzz_ts_analysis_tissue_specific_clusters_line_plot")


#---------------------------------- Functions ----------------------------------
line_plot_for_tissue_specific_clusters = function(cluster_result, cluster_info_tissue_df, num_of_clusters, tissue_th, plot_category, xticks_names, data_input, plot_props, output_folder_path_fig, file_name) {
  
  font_size_correction = 2.845
  
  important_tissues = c()
  for (cluster in 1:num_of_clusters) {
    mask = ((cluster_info_tissue_df$CLUSTER == cluster) & (cluster_info_tissue_df$relative_features >= tissue_th))
    cluster_info_tissue_df_cluster = cluster_info_tissue_df[mask,]
    important_tissues = rbind(important_tissues, cluster_info_tissue_df_cluster)
  }
  
  cluster_center_lines = melt(cluster_result$centers)
  colnames(cluster_center_lines) = c("CLUSTER", plot_category[2], "values")
  
  plot_df = merge(cluster_center_lines, important_tissues, by="CLUSTER")
  
  plot_df$time = factor(plot_df[, plot_category[2]], levels = sort(unique(as.numeric(plot_df[, plot_category[2]]))))
  
  plot_df$tissue_cluster = sprintf("%s_%s", plot_df$feature, plot_df$CLUSTER)
  
  p = ggplot(plot_df, aes(x = time, y = values, group = tissue_cluster)) + 
    geom_line(aes(color=feature)) +
    #geom_text_repel(aes(label = Var1), data = data_plot[data_plot$Var2 == sort(unique(data_plot$Var2))[1],], color = "black", size = 3) +
    #scale_x_continuous(breaks = sort(unique(data_plot$time))) +#, labels = 1:10) 
    scale_color_manual(values=unlist(colours[[plot_category[1]]][unique(annot[[plot_category[1]]])])) +
    xlab(sprintf("%s", xticks_names$categories[[plot_category[2]]])) + 
    ylab(sprintf("Standardised expression\nvalues (%s)", data_input$norm)) +
    #geom_text(x = xpos, y = ypos, label = num_of_lines, hjust = 1, vjust = 1, family = plots_props$font_family, size = plots_props$font_size/font_size_correction) +
    theme_classic() +
    theme(#plot.margin = margin(0.1, 0, -0.25, 0, "cm"), #top, right, bottom, and left 
          legend.position="none", 
          #legend.margin = margin(t = -10, r = 5, b = -5, l = -5),
          text = element_text(family = plot_props$font_family, size = plot_props$font_size),
          axis.text = element_text(family = plot_props$font_family, size = plot_props$font_size),
          axis.title = element_text(family = plot_props$font_family, size = plot_props$font_size),
          plot.title = element_text(family = plot_props$font_family, size = plot_props$font_size_header),
          legend.text = element_text(family = plot_props$font_family, size = plot_props$font_size),
          legend.title = element_text(family = plot_props$font_family, size = plot_props$font_size)
    ) +
    guides(color = guide_legend(nrow = 3))
  #p
  
  # save bar plot
  ggsave(sprintf("%s/%s.png", output_folder_path_fig, file_name), p, dpi=plot_props$dpi, width=plot_props$image_width, height=plot_props$image_height, units = plot_props$image_units)
  ggsave(sprintf("%s/%s.svg", output_folder_path_fig, file_name), p, width=plot_props$image_width, height=plot_props$image_height, units =plot_props$image_units, limitsize = FALSE, scale = 1)
  
  return(important_tissues)
}

line_plot_for_tissue_specific_clusters_adj = function(cluster_result, cluster_info_tissue_df, num_of_clusters, tissue_th, occ_th, plot_category, xticks_names, data_input, plot_props, output_folder_path_fig, file_name) {
  
  font_size_correction = 2.845
  
  colours_brain_regions = colours[[plot_category[1]]]
  colours_brain_regions$low = "#DFDFB9"
  colours_brain_regions_adj = unlist(colours_brain_regions[c(unique(annot[[plot_category[1]]]), "low")])
  
  important_tissues = c()
  for (cluster in 1:num_of_clusters) {
    mask = ((cluster_info_tissue_df$CLUSTER == cluster) & (cluster_info_tissue_df$relative_features >= tissue_th))
    cluster_info_tissue_df_cluster = cluster_info_tissue_df[mask,]
    important_tissues = rbind(important_tissues, cluster_info_tissue_df_cluster)
  }
  
  important_tissues$tissue_occ = 0
  for (ti in unique(important_tissues$feature)) {
    tmp = (important_tissues$feature == ti)
    important_tissues[tmp]$tissue_occ = dim(important_tissues[tmp])[1]
  }
  
  cluster_center_lines = melt(cluster_result$centers)
  colnames(cluster_center_lines) = c("CLUSTER", plot_category[2], "values")
  
  plot_df = merge(cluster_center_lines, important_tissues, by="CLUSTER")
  
  plot_df$time = factor(plot_df[, plot_category[2]], levels = sort(unique(as.numeric(plot_df[, plot_category[2]]))))
  
  plot_df$tissue_cluster = sprintf("%s_%s", plot_df$feature, plot_df$CLUSTER)
  
  plot_df$colouring = "low"
  plot_df[(plot_df$tissue_occ >= occ_th),]$colouring = plot_df[(plot_df$tissue_occ >= occ_th),]$feature

  plot_df$thickness = 0
  plot_df[(plot_df$tissue_occ >= occ_th),]$thickness = 1

  p = ggplot(plot_df, aes(x = time, y = values, group = tissue_cluster)) + 
    #geom_line(aes(color=colouring, size = thickness)) +
    # split the plotting in a way that we have the coloured thicker lines in the foreground
    geom_line(data = subset(plot_df, thickness == 0), aes(color=colouring, size = thickness)) +
    geom_line(data = subset(plot_df, thickness == 1), aes(color=colouring, size = thickness)) +
    #geom_text_repel(aes(label = Var1), data = data_plot[data_plot$Var2 == sort(unique(data_plot$Var2))[1],], color = "black", size = 3) +
    #scale_x_continuous(breaks = sort(unique(data_plot$time))) +#, labels = 1:10) 
    scale_color_manual(values=colours_brain_regions_adj) +
    scale_size_continuous(range = c(0.25, 0.75)) + 
    xlab(sprintf("%s", xticks_names$categories[[sprintf("%s_capital", plot_category[2])]])) + 
    ylab(sprintf("Std. expression\nvalues (%s)", str_split_fixed(data_input$norm, "_", 2)[1])) +
    #geom_text(x = xpos, y = ypos, label = num_of_lines, hjust = 1, vjust = 1, family = plots_props$font_family, size = plots_props$font_size/font_size_correction) +
    scale_x_discrete(expand = c(0, 0.1)) +
    theme_classic() +
    theme(#plot.margin = margin(0.1, 0, -0.25, 0, "cm"), #top, right, bottom, and left 
      legend.position="none", 
      #legend.margin = margin(t = -10, r = 5, b = -5, l = -5),
      text = element_text(family = plot_props$font_family, size = plot_props$font_size),
      axis.text = element_text(family = plot_props$font_family, size = plot_props$font_size),
      axis.title = element_text(family = plot_props$font_family, size = plot_props$font_size),
      plot.title = element_text(family = plot_props$font_family, size = plot_props$font_size_header),
      legend.text = element_text(family = plot_props$font_family, size = plot_props$font_size),
      legend.title = element_text(family = plot_props$font_family, size = plot_props$font_size)
    ) +
    guides(color = guide_legend(nrow = 3))
  #p
  
  # save bar plot
  plot_height = 2/3 * plot_props$image_height
  plot_width = 2/3 * plot_props$image_width
  ggsave(sprintf("%s/%s.png", output_folder_path_fig, file_name), p, dpi=plot_props$dpi, width=plot_width, height=plot_height, units = plot_props$image_units)
  ggsave(sprintf("%s/%s.svg", output_folder_path_fig, file_name), p, width=plot_width, height=plot_height, units =plot_props$image_units, limitsize = FALSE, scale = 1)
}

line_plot_for_tissue_specific_clusters_manual_selection_adj = function(cluster_result, cluster_info_tissue_df, num_of_clusters, tissue_th, plot_category, xticks_names, data_input, plot_props, output_folder_path_fig, file_name) {
  
  font_size_correction = 2.845
  
  colours_brain_regions = colours[[plot_category[1]]]
  colours_brain_regions$low = "#DFDFB9"
  colours_brain_regions_adj = unlist(colours_brain_regions[c(unique(annot[[plot_category[1]]]), "low")])
  
  important_tissues = c()
  for (cluster in 1:num_of_clusters) {
    mask = ((cluster_info_tissue_df$CLUSTER == cluster) & (cluster_info_tissue_df$relative_features >= tissue_th))
    cluster_info_tissue_df_cluster = cluster_info_tissue_df[mask,]
    important_tissues = rbind(important_tissues, cluster_info_tissue_df_cluster)
  }
  
  important_tissues$tissue_occ = 0
  for (ti in unique(important_tissues$feature)) {
    tmp = (important_tissues$feature == ti)
    important_tissues[tmp]$tissue_occ = dim(important_tissues[tmp])[1]
  }
  
  cluster_center_lines = melt(cluster_result$centers)
  colnames(cluster_center_lines) = c("CLUSTER", plot_category[2], "values")
  
  plot_df = merge(cluster_center_lines, important_tissues, by="CLUSTER")
  
  plot_df$time = factor(plot_df[, plot_category[2]], levels = sort(unique(as.numeric(plot_df[, plot_category[2]]))))
  
  plot_df$tissue_cluster = sprintf("%s_%s", plot_df$feature, plot_df$CLUSTER)
  
  plot_df$colouring = "low"
  plot_df[plot_df$feature == "plx",]$colouring = plot_df[plot_df$feature == "plx",]$feature
  plot_df[plot_df$feature == "pon",]$colouring = plot_df[plot_df$feature == "pon",]$feature
  
  plot_df$thickness = 0
  plot_df[((plot_df$feature == "pon") | (plot_df$feature == "plx")),]$thickness = 1
  
  p = ggplot(plot_df, aes(x = time, y = values, group = tissue_cluster)) + 
    #geom_line(aes(color=colouring, size = thickness)) +
    # split the plotting in a way that we have the coloured thicker lines in the foreground
    geom_line(data = subset(plot_df, thickness == 0), aes(color=colouring, size = thickness)) +
    geom_line(data = subset(plot_df, thickness == 1), aes(color=colouring, size = thickness)) +
    #geom_text_repel(aes(label = Var1), data = data_plot[data_plot$Var2 == sort(unique(data_plot$Var2))[1],], color = "black", size = 3) +
    #scale_x_continuous(breaks = sort(unique(data_plot$time))) +#, labels = 1:10) 
    scale_color_manual(values=colours_brain_regions_adj) +
    scale_size_continuous(range = c(0.25, 0.75)) + 
    xlab(sprintf("%s", xticks_names$categories[[sprintf("%s_capital", plot_category[2])]])) + 
    ylab(sprintf("Std. expression\nvalues (%s)", str_split_fixed(data_input$norm, "_", 2)[1])) +
    #geom_text(x = xpos, y = ypos, label = num_of_lines, hjust = 1, vjust = 1, family = plots_props$font_family, size = plots_props$font_size/font_size_correction) +
    scale_x_discrete(expand = c(0, 0.1)) +
    theme_classic() +
    theme(#plot.margin = margin(0.1, 0, -0.25, 0, "cm"), #top, right, bottom, and left 
      legend.position="none", 
      #legend.margin = margin(t = -10, r = 5, b = -5, l = -5),
      text = element_text(family = plot_props$font_family, size = plot_props$font_size),
      axis.text = element_text(family = plot_props$font_family, size = plot_props$font_size),
      axis.title = element_text(family = plot_props$font_family, size = plot_props$font_size),
      plot.title = element_text(family = plot_props$font_family, size = plot_props$font_size_header),
      legend.text = element_text(family = plot_props$font_family, size = plot_props$font_size),
      legend.title = element_text(family = plot_props$font_family, size = plot_props$font_size)
    ) +
    guides(color = guide_legend(nrow = 3))
  #p
  
  # save bar plot
  plot_height = 2/3 * plot_props$image_height
  plot_width = 2/3 * plot_props$image_width
  ggsave(sprintf("%s/%s.png", output_folder_path_fig, file_name), p, dpi=plot_props$dpi, width=plot_width, height=plot_height, units = plot_props$image_units)
  ggsave(sprintf("%s/%s.svg", output_folder_path_fig, file_name), p, width=plot_width, height=plot_height, units =plot_props$image_units, limitsize = FALSE, scale = 1)
}


#---------------------------------- Load data ----------------------------------
data_input = snakemake@params$data_input
cluster_props = snakemake@params$cluster_props
plot_category = snakemake@params$plot_category
diff_exp_log = snakemake@params$diff_log
thresholds = snakemake@params$thresholds
colours = snakemake@params$colors
plot_props = snakemake@params$plots_props
results_folder = snakemake@params$results_folder
xticks_names = snakemake@params$xticks_names

# save rdata
#saveRDS(snakemake, file = "snakemake_mfuzz_ts_analysis_tissue_specific_clusters_line_plot.rds")

expr = fread(sprintf("%s_%s/%s_%s_quantification_%s.detection_rate_%s_per_%s.csv", data_input$raw_data_folder, data_input$data_sub_set, data_input$rna_class, data_input$feature_filtering, data_input$norm, data_input$detection_rate, data_input$detection_group), sep='\t', header=T)
annot = fread(sprintf("%s/annotation_%s_%s.csv", data_input$annotation_data_folder, data_input$data_sub_set, data_input$feature_filtering), sep='\t', colClasses=c(ID="character")) 
diff_exp = fread(sprintf("%s/%s_%s/results_%s/matrices/diff_exp/diff_exp%s.csv", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, diff_exp_log), sep='\t') 

feature_list = fread(snakemake@input$feature_list, sep='\t')  # , colClasses=c(ID="character")

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

for (membership_plot_th in thresholds$membership_plot_ths) {
  median_table_rownames = median_table_df$feature
  for (th in thresholds$tissue_th) {
    # tissue containment in each cluster
    # read info table
    file_name = sprintf("cluster_overview_%s", plot_category[1])
    cluster_info_tissue_df = fread(sprintf("%s/%s/%s.csv", input_folder_path_tab, membership_plot_th, file_name), sep='\t') 
    
    # make plot showing the centerlines for the clusters with a tissue occ higher than a given threshold
    output_folder_path_fig_th = sprintf("%s/%s/line_plots", output_folder_path_fig, membership_plot_th)
    dir.create(output_folder_path_fig_th, recursive=TRUE)
    file_name = sprintf("tissue_specific_clusters_more_%s_from_one_%s_center_lines", th, plot_category[1])
    specific_tissues = line_plot_for_tissue_specific_clusters(cluster_result, cluster_info_tissue_df, cluster_props$num_of_clusters, th, plot_category, xticks_names, data_input, plot_props, output_folder_path_fig_th, file_name)
    
    # adjust by only colouring the tissues with an occurrence of more than occ_th in this plot
    occ_th = 3
    file_name = sprintf("tissue_specific_clusters_more_%s_from_one_%s_center_lines_colouring_th=%s", th, plot_category[1], occ_th)
    line_plot_for_tissue_specific_clusters_adj(cluster_result, cluster_info_tissue_df, cluster_props$num_of_clusters, th, occ_th, plot_category, xticks_names, data_input, plot_props, output_folder_path_fig_th, file_name)
    
    file_name = sprintf("tissue_specific_clusters_more_%s_from_one_%s_center_lines_colouring_adj", th, plot_category[1])
    line_plot_for_tissue_specific_clusters_manual_selection_adj(cluster_result, cluster_info_tissue_df, cluster_props$num_of_clusters, th, plot_category, xticks_names, data_input, plot_props, output_folder_path_fig_th, file_name)
    
    #cluster_line_plots_with_tissue_specific_info(median_table_df, cluster_result, acore_dt, th, specific_tissues, data_input, output_folder_path_fig_th, output_folder_path_tab, cluster_props$num_of_clusters, plot_props)
    
    # feature occurances in each cluster
    # read info table
    #file_name = sprintf("cluster_overview.membership_at_least_%s.%s_occ", membership_plot_th, data_input$rna_class)
    #cluster_info_feature_df = fread(sprintf("%s/%s/%s.csv", input_folder_path_tab, membership_plot_th, file_name), sep='\t') 
  }
}


#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])
