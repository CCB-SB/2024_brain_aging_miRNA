suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(svglite))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(ggpubr))


#snakemake = readRDS("snakemake_mfuzz_ts_analysis_scatter_plot.rds")

set.seed(snakemake@params$parameters_props$set_seed)

print("mfuzz_ts_analysis_scatter_plot")


create_overview_scatter_plot = function(plot_df, tissue_th, feature_th, x_axis_name, y_axis_name, file_name, data_input, plots_props, output_folder_fig_scatter) {
  # generate vector for colours
  colour <- as.character(plot_df$colour)
  names(colour) <- as.character(plot_df$colour)
  
  maxi = max(plot_df$fs_spec_value_max_abs)
  
  p = ggplot(plot_df, aes(x=ts_spec_value_max_rel*100, y=fs_spec_value_max_abs)) + 
    geom_hline(yintercept = feature_th, color="lightgrey", linetype="dashed", size = 0.25) + 
    geom_vline(xintercept = tissue_th*100, color="lightgrey", linetype="dashed", size = 0.25) + 
    geom_point(aes(colour=colour, size=cardinality), alpha = 0.8) +
    scale_size_continuous(range = c(0.25,4)) +
    scale_y_continuous(breaks = seq(0, maxi, by = 2)) +
    xlab(x_axis_name) + ylab(y_axis_name) +
    scale_color_manual(name="", values=colour, guide = "none") + 
    theme_classic() +
    theme(legend.position = "bottom", 
          #plot.background = element_rect(fill="transparent"),
          #panel.background = element_rect(fill="transparent"),
          text = element_text(family = plots_props$font_family, size = plots_props$font_size), 
          axis.text = element_text(family = plots_props$font_family, size = plots_props$font_size), 
          axis.title = element_text(family = plots_props$font_family, size = plots_props$font_size), 
          plot.title = element_text(family = plots_props$font_family, size = plots_props$font_size), 
          legend.text = element_text(family = plots_props$font_family, size = plots_props$font_size),
          legend.title = element_text(family = plots_props$font_family, size = plots_props$font_size, face = "bold"),
          plot.margin = unit(c(0.25,0.25,0.25,0.25), "cm"), #t, r, b, l
          legend.spacing.y = unit(0, "cm"),
          legend.margin = margin(0,0,0,0),
          legend.box.margin = margin(-0.2,0,-0.2,0, "cm"),
          legend.spacing.x = unit(0, "cm")
    ) +
    guides(size = guide_legend(title = "Cardinality", nrow = 1, title.position = "top")   # Shape legend in one row
    )
  p_label = p + geom_text(aes(label = CLUSTER_name), nudge_x = 0, nudge_y = -0.6,
                          size = 6/11.04*3.88, family = plots_props$font_family, color = "black")
  #geom_text_repel(aes(label = CLUSTER_name),
  #  max.overlaps = Inf,
  #  box.padding = 0.1,
  #  direction = "y",
  #  force_pull = 2,
  #  force = 15,
  #  point.padding = 0.3,
  #  min.segment.length = 0, segment.size = 0.25,
  #  seed = snakemake@params$parameters_porps$set_seed,
  #  size = 6/11.04*3.88, family = plots_props$font_family, color = "black")
  #p_label
  
  height = plots_props$image_height
  width = plots_props$image_width
  ggsave(sprintf("%s/cluster_overview_tissue_occ_th=%s_%s_occ_th=%s_%sx%s.png", output_folder_fig_scatter, tissue_th, data_input$rna_class, feature_th, width, height), p_label, dpi=plots_props$dpi, width = width, height = height, units = plots_props$image_units)
  ggsave(sprintf("%s/cluster_overview_tissue_occ_th=%s_%s_occ_th=%s_%sx%s.svg", output_folder_fig_scatter, tissue_th, data_input$rna_class, feature_th, width, height), p_label, width = width, height = height, units =plots_props$image_units, limitsize = FALSE, scale = 1)
  
  height = 5
  width = plots_props$image_width
  ggsave(sprintf("%s/cluster_overview_tissue_occ_th=%s_%s_occ_th=%s_%sx%s.png", output_folder_fig_scatter, tissue_th, data_input$rna_class, feature_th, width, height), p_label, dpi=plots_props$dpi, width = width, height = height, units = plots_props$image_units)
  ggsave(sprintf("%s/cluster_overview_tissue_occ_th=%s_%s_occ_th=%s_%sx%s.svg", output_folder_fig_scatter, tissue_th, data_input$rna_class, feature_th, width, height), p_label, width = width, height = height, units =plots_props$image_units, limitsize = FALSE, scale = 1)
}

create_scatter_plot = function(plot_df_filtered, tissue_spec_filtered, feature_spec_filtered, tissue, tissue_th, feature_th, x_axis_name_left, x_axis_name_right, colours, file_name, plots_props, output_folder_fig_scatter) {
  cluster_sorted = sort(plot_df_filtered$CLUSTER)
  start_step = -1
  step_size = 2
  levels_with_dividers = c(step_size*length(cluster_sorted) + start_step, rbind(sprintf("Cluster %s", as.character(rev(cluster_sorted))), rev(seq(start_step, step_size*length(cluster_sorted) - step_size, step_size))))
  half_yticks = seq(start_step, step_size*length(cluster_sorted) + abs(start_step), step_size)
  
  if ((length(half_yticks)-1)%%2 == 1) {
    unodd_compensation = c(1)
  } else {
    unodd_compensation = c()
  }
  rects_left = data.frame(ystart = half_yticks[-length(half_yticks)], yend = half_yticks[-1], col = append(rep(c(1,0), floor((length(half_yticks)-1)/2)), unodd_compensation))
  rects_right = data.frame(ystart = half_yticks[-length(half_yticks)], yend = half_yticks[-1], col = rev(append(rep(c(1,0), floor((length(half_yticks)-1)/2)), unodd_compensation)))
  
  plot_df_filtered$CLUSTER_factor = factor(plot_df_filtered$CLUSTER, levels = cluster_sorted)
  
  #visibility in plot
  fs_th = 2
  fs_absolute_features_show = c()
  for (cl in unique(feature_spec_filtered$CLUSTER)) {
    feature_spec_filtered_cl = feature_spec_filtered[feature_spec_filtered$CLUSTER == cl,]
    for (row in 1:dim(feature_spec_filtered_cl)[1]) {
      feature_spec_filtered_cl_row = feature_spec_filtered_cl[row]
      if (feature_spec_filtered_cl_row$absolute_features > fs_th) {
        fs_absolute_features_show = append(fs_absolute_features_show, 1)
      } else {
        fs_absolute_features_show = append(fs_absolute_features_show, 0)
      }
    }
  }
  
  ts_th = 0.1
  ts_relative_features_show = c()
  for (cl in unique(tissue_spec_filtered$CLUSTER)) {
    tissue_spec_filtered_cl = tissue_spec_filtered[tissue_spec_filtered$CLUSTER == cl,]
    for (row in 1:dim(tissue_spec_filtered_cl)[1]) {
      tissue_spec_filtered_cl_row = tissue_spec_filtered_cl[row]
      if (tissue_spec_filtered_cl_row$relative_features > ts_th) {
        ts_relative_features_show = append(ts_relative_features_show, 1)
      } else {
        ts_relative_features_show = append(ts_relative_features_show, 0)
      }
    }
  }
  
  feature_spec_filtered$absolute_features_show = feature_spec_filtered$absolute_features * fs_absolute_features_show
  tissue_spec_filtered$relative_features_show = tissue_spec_filtered$relative_features * ts_relative_features_show
  
  feature_spec_filtered$CLUSTER_name = sprintf("Cluster %s", feature_spec_filtered$CLUSTER)
  feature_spec_filtered$CLUSTER_name_factor = factor(feature_spec_filtered$CLUSTER_name, levels = sprintf("Cluster %s", cluster_sorted))
  tissue_spec_filtered$CLUSTER_name = sprintf("Cluster %s", tissue_spec_filtered$CLUSTER)
  tissue_spec_filtered$CLUSTER_name_factor = factor(tissue_spec_filtered$CLUSTER_name, levels = sprintf("Cluster %s", cluster_sorted))
  
  x_shift = -0.4
  
  # plot for axis labels
  y_labels_df = data.frame(x=0, y=rev(unique(tissue_spec_filtered$CLUSTER_name_factor)))
  
  scatter_y_label_plot = ggplot(y_labels_df, aes(x=x, y=y)) +
    geom_text(aes(label=y), size=plots_props$font_size * 5 / 14 * 2 / 3) +
    scale_y_discrete(breaks = rev(feature_spec_filtered$CLUSTER_name_factor), limits = levels_with_dividers, expand=c(0, 0)) +
    theme_void() +
    theme(plot.margin = unit(c(0, -0.4, 0.1, x_shift), "cm"),)
  
  # plot ratio scaling
  # for a taller narrow plot, multiply ratio_scaling > 1
  ratio_scaling = 1.5
  
  # left plot tissue spec
  # border
  tissue_spec_filtered$relative_features_show_border = ifelse(tissue_spec_filtered$relative_features >= tissue_th, 1, 0)
  
  maxi_left = max(tissue_spec_filtered$relative_features)
  
  # fix the ratio of the plot area
  x_min_left = 0
  ratio_left = (maxi_left - x_min_left) * 100 / length(levels_with_dividers) * ratio_scaling
  
  scatter_plot_left = ggplot(tissue_spec_filtered, aes(x=relative_features_show*100, y=CLUSTER_name_factor)) +
    #geom_rect(data=rects_left, aes(x=NULL, y=NULL, xmin=-Inf, xmax=Inf, ymin=ystart, ymax=yend, fill=col), show.legend=FALSE) + 
    #scale_fill_manual(values=c("#E5E5E5", "#C6C6C6")) +
    geom_vline(xintercept = tissue_th*100, color="lightgrey", linetype="dashed", size = 0.25) + 
    # points with border
    geom_point(data=tissue_spec_filtered[tissue_spec_filtered$relative_features_show_border==1, ], aes(fill=feature), color="black", shape=21, size=1.5) +
    # points without border
    geom_point(data=tissue_spec_filtered[tissue_spec_filtered$relative_features_show_border!=1, ], aes(color=feature), size=1.5) +
    # included coloured background according to the tissue coloures
    #geom_rect(data=rects, aes(x=NULL, y=NULL, xmin=-Inf, xmax=Inf, ymin=ystart, ymax=yend, fill=col), alpha = 0.4, show.legend=FALSE) + 
    ylab("") +
    xlab(x_axis_name_left) +
    coord_fixed(ratio = ratio_left) +
    # only show breaks for roman numbers
    scale_y_discrete(breaks = rev(tissue_spec_filtered$CLUSTER_name_factor), limits = levels_with_dividers, expand=c(0, 0)) +
    #scale_x_continuous(breaks = append(1,seq(0, maxi_left*100, by = 20)[-1]), limits = c(1, maxi_left*100)) +
    scale_x_reverse(lim=c(100, 1), breaks=c(100, 75, 50, 25, 10)) +
    #xlim(100,1) +
    scale_fill_manual(values=unlist(colours[[tissue]])) +
    scale_color_manual(values=unlist(colours[[tissue]])) +
    theme_classic() +
    #coord_cartesian(ylim = c(1, 4)) +  # Limit y-axis to the range of factors
    theme(plot.background = element_rect(fill='transparent', color=NA),
          legend.position="none", 
          text = element_text(family = plots_props$font_family, size = plots_props$font_size),
          axis.title = element_text(family = plots_props$font_family, size = plots_props$font_size),
          axis.text = element_text(family = plots_props$font_family, size = 6),
          axis.line = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(family = plots_props$font_family, size = 6), #angle = 90, vjust = 0.5, hjust=1, 
          #axis.ticks.x = element_blank(),
          plot.title = element_text(family = plots_props$font_family, size = plots_props$font_size_header),
          legend.text = element_text(family = plots_props$font_family, size = plots_props$font_size),
          legend.title = element_text(family = plots_props$font_family, size = plots_props$font_size),
          #plot.margin = unit(c(-0.25,1,-0.25,0.25), "cm")
          plot.margin = unit(c(0,0,0.1,x_shift), "cm"),
    )
  for (hline_pos in half_yticks[-c(1, length(half_yticks))]) {
    # durchgehende Linie
    scatter_plot_left = scatter_plot_left + geom_hline(yintercept=hline_pos, size = 0.3, color = "#727272")
  } 
  #scatter_plot_left_label = scatter_plot_left + geom_text(aes(label = feature), nudge_x = 0.008, nudge_y = -0.4,
  #                                                          size = 6/11.04*3.88, family = plots_props$font_family, color = "black")
  #scatter_plot_left
  
  # right plot feature spec
  # font colour
  feature_spec_filtered$absolute_features_show_border = ifelse(feature_spec_filtered$absolute_features >= feature_th, 1, 0)
  feature_spec_filtered$absolute_features_show_border = as.factor(feature_spec_filtered$absolute_features_show_border)
  
  maxi_right = max(feature_spec_filtered$absolute_features)
  
  # fix the ratio of the plot area
  x_min_right = fs_th + 1 - 0.5
  if (x_min_right >= maxi_right) {
    x_min_right = 0.5
  }
  ratio_right = (maxi_right - x_min_right) / length(levels_with_dividers) * ratio_scaling
  
  feature_spec_filtered$feature_name = gsub("mmu-miR-", "", feature_spec_filtered$feature)
  feature_spec_filtered$feature_name = gsub("mmu-let-", "", feature_spec_filtered$feature_name)
  
  scatter_plot_right = ggplot(feature_spec_filtered, aes(x=absolute_features_show, y=CLUSTER_name_factor)) +
    #geom_rect(data=rects_right, aes(x=NULL, y=NULL, xmin=-Inf, xmax=Inf, ymin=ystart, ymax=yend, fill=col), show.legend=FALSE) + 
    #scale_fill_manual(values=c("#E5E5E5", "#C6C6C6")) +
    geom_vline(xintercept = feature_th, color="lightgrey", linetype="dashed", size = 0.25) + 
    geom_point(aes(color=absolute_features_show_border), size = 1.5) +
    # included coloured background according to the tissue coloures
    #geom_rect(data=rects, aes(x=NULL, y=NULL, xmin=-Inf, xmax=Inf, ymin=ystart, ymax=yend, fill=col), alpha = 0.4, show.legend=FALSE) + 
    ylab("") +
    xlab(x_axis_name_right) +
    # only show breaks for roman numbers
    scale_y_discrete(breaks = rev(feature_spec_filtered$CLUSTER_name_factor), limits = levels_with_dividers, expand=c(0, 0)) +
    scale_x_continuous(breaks = seq(x_min_right + 0.5, maxi_right, by = 2), limits = c(x_min_right, maxi_right)) +
    scale_color_manual(values=c("0"="#DFDFB9", "1"="#525200")) +
    coord_fixed(ratio=ratio_right) +
    theme_classic() +
    #coord_cartesian(ylim = c(1, 4)) +  # Limit y-axis to the range of factors
    theme(plot.background = element_rect(fill='transparent', color=NA),
          legend.position="none", 
          text = element_text(family = plots_props$font_family, size = plots_props$font_size),
          axis.title = element_text(family = plots_props$font_family, size = plots_props$font_size),
          axis.text = element_text(family = plots_props$font_family, size = 6),
          axis.line = element_blank(),
          axis.text.y = element_blank(),#element_text(hjust = 0.5),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(family = plots_props$font_family, size = 6), #angle = 90, vjust = 0.5, hjust=1, 
          #axis.ticks.x = element_blank(),
          plot.title = element_text(family = plots_props$font_family, size = plots_props$font_size_header),
          legend.text = element_text(family = plots_props$font_family, size = plots_props$font_size),
          legend.title = element_text(family = plots_props$font_family, size = plots_props$font_size),
          #plot.margin = unit(c(-0.25,0.25,-0.25,-1.25), "cm")
          plot.margin = unit(c(0,0,0.1,x_shift), "cm"),
    )
  for (hline_pos in half_yticks[-c(1, length(half_yticks))]) {
    # durchgehende Linie
    scatter_plot_right = scatter_plot_right + geom_hline(yintercept=hline_pos, size = 0.3, color = "#727272")
  } 
  scatter_plot_right_label = scatter_plot_right +
    # only label above th
    geom_text(aes(label = feature_name),
              data=feature_spec_filtered[feature_spec_filtered$absolute_features_show_border==1, ],
              nudge_x = -0.2, nudge_y = 0, hjust=1,
              size = 4/11.04*3.88, family = plots_props$font_family, color = "black", 
              #position = position_jitter(width = 0, height = 0.25, seed = snakemake@params$parameters_props$set_seed)
              )
  
  # combine the two plots
  combined_plot = ggarrange(scatter_plot_left,
                            scatter_y_label_plot,
                            scatter_plot_right,
                            nrow=1,
                            align="h",
                            widths=c(1, 0.25, 1)) 
  combined_plot_label = ggarrange(scatter_plot_left,
                                  scatter_y_label_plot,
                                  scatter_plot_right_label,
                                  nrow=1,
                                  align="h",
                                  widths=c(1, 0.25, 1)) 
  
  height = plots_props$image_height
  width = plots_props$image_width
  ggsave(sprintf("%s/%s_tissue_occ_th=%s_%s_occ_th=%s_%sx%s.png", output_folder_fig_scatter, file_name, tissue_th, data_input$rna_class, feature_th, width, height), combined_plot, dpi=plots_props$dpi, width=width, height=height, units = plots_props$image_units)
  ggsave(sprintf("%s/%s_tissue_occ_th=%s_%s_occ_th=%s_%sx%s.svg", output_folder_fig_scatter, file_name, tissue_th, data_input$rna_class, feature_th, width, height), combined_plot, width=width, height=height, units = plots_props$image_units, limitsize = FALSE, scale = 1)

  ggsave(sprintf("%s/%s_tissue_occ_th=%s_%s_occ_th=%s_%sx%s_label.png", output_folder_fig_scatter, file_name, tissue_th, data_input$rna_class, feature_th, width, height), combined_plot_label, dpi=plots_props$dpi, width=width, height=height, units = plots_props$image_units)
  ggsave(sprintf("%s/%s_tissue_occ_th=%s_%s_occ_th=%s_%sx%s_label.svg", output_folder_fig_scatter, file_name, tissue_th, data_input$rna_class, feature_th, width, height), combined_plot_label, width=width, height=height, units = plots_props$image_units, limitsize = FALSE, scale = 1)
  
  height = 4
  width = plots_props$image_width
  ggsave(sprintf("%s/%s_tissue_occ_th=%s_%s_occ_th=%s_%sx%s.png", output_folder_fig_scatter, file_name, tissue_th, data_input$rna_class, feature_th, width, height), combined_plot, dpi=plots_props$dpi, width=width, height=height, units = plots_props$image_units)
  ggsave(sprintf("%s/%s_tissue_occ_th=%s_%s_occ_th=%s_%sx%s.svg", output_folder_fig_scatter, file_name, tissue_th, data_input$rna_class, feature_th, width, height), combined_plot, width=width, height=height, units = plots_props$image_units, limitsize = FALSE, scale = 1)
  
  ggsave(sprintf("%s/%s_tissue_occ_th=%s_%s_occ_th=%s_%sx%s_label.png", output_folder_fig_scatter, file_name, tissue_th, data_input$rna_class, feature_th, width, height), combined_plot_label, dpi=plots_props$dpi, width=width, height=height, units = plots_props$image_units)
  ggsave(sprintf("%s/%s_tissue_occ_th=%s_%s_occ_th=%s_%sx%s_label.svg", output_folder_fig_scatter, file_name, tissue_th, data_input$rna_class, feature_th, width, height), combined_plot_label, width=width, height=height, units = plots_props$image_units, limitsize = FALSE, scale = 1)
  
  #height = plots_props$image_height
  #width = 3/2 * plots_props$image_width 
  #ggsave(sprintf("%s/%s_tissue_occ_th=%s_%s_occ_th=%s_%sx%s.png", output_folder_fig_scatter, file_name, tissue_th, data_input$rna_class, feature_th, width, height), combined_plot, dpi=plots_props$dpi, width=width, height=height, units = plots_props$image_units)
  #ggsave(sprintf("%s/%s_tissue_occ_th=%s_%s_occ_th=%s_%sx%s.svg", output_folder_fig_scatter, file_name, tissue_th, data_input$rna_class, feature_th, width, height), combined_plot, width=width, height=height, units = plots_props$image_units, limitsize = FALSE, scale = 1)
  
  #ggsave(sprintf("%s/%s_tissue_occ_th=%s_%s_occ_th=%s_%sx%s_label.png", output_folder_fig_scatter, file_name, tissue_th, data_input$rna_class, feature_th, width, height), combined_plot_label, dpi=plots_props$dpi, width=width, height=height, units = plots_props$image_units)
  #ggsave(sprintf("%s/%s_tissue_occ_th=%s_%s_occ_th=%s_%sx%s_label.svg", output_folder_fig_scatter, file_name, tissue_th, data_input$rna_class, feature_th, width, height), combined_plot_label, width=width, height=height, units = plots_props$image_units, limitsize = FALSE, scale = 1)
}


#---------------------------------- Load data ---------------------------------- 
data_input = snakemake@params$data_input
cluster_props = snakemake@params$cluster_props
tissue = snakemake@params$plot_category[1]
time = snakemake@params$plot_category[2]
thresholds = snakemake@params$thresholds
xticks_names = snakemake@params$xticks_names
colours = snakemake@params$colors
plots_props = snakemake@params$plots_props
results_folder = snakemake@params$results_folder

# save rdata
#saveRDS(snakemake, file = "snakemake_mfuzz_ts_analysis_scatter_plot.rds")


# --------------------------------- Script -------------------------------------
for (membership_plot_th in c("all", thresholds$membership_plot_ths)) {
  for (tissue_th in thresholds$tissue_ths) {
    for (feature_th in thresholds$rna_occ_ths) { 
      # specify inout folder
      input_folder_path = sprintf("%s/%s_%s/results_%s/matrices/mfuzz_clustering/k=%s", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, cluster_props$num_of_clusters)
      print(membership_plot_th)
      tissue_spec = fread(sprintf("%s/%s/cluster_overview_brain_region.csv", input_folder_path, membership_plot_th), sep='\t', header=T)
      feature_spec = fread(sprintf("%s/%s/cluster_overview_%s_occ.csv", input_folder_path, membership_plot_th, data_input$rna_class), sep='\t', header=T)
      membership = fread(sprintf("%s/cluster_members_k=%s.csv", snakemake@input$cluster_results_folder_path, cluster_props$num_of_clusters), sep=',', header=T) 
      membership = membership[membership$MEM.SHIP >= membership_plot_th,]
      
      # create output folder
      output_folder_fig_scatter = sprintf("%s/%s_%s/results_%s/figures/mfuzz_clustering/k=%s/%s/scatter_plots", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, cluster_props$num_of_clusters, membership_plot_th)
      dir.create(output_folder_fig_scatter, recursive=TRUE)
      
      #ts_string = sprintf("%s specific", xticks_names$categories[[sprintf("%s_capital", tissue)]]) #tissue_th
      #fs_string = sprintf("%s specific", data_input$rna_class) #feature_th
      #ts_fs_string = sprintf("%s and %s specific", xticks_names$categories[[sprintf("%s_capital", tissue)]], data_input$rna_class)
      #rest_string = "Not specific"
      
      colors = list()
      colors[["ts"]] = unlist(colours[[tissue]])
      colors[["fs"]] = "#525200"
      colors[["ts_fs"]] = unlist(colours[[tissue]])
      colors[["rest"]] = "#DFDFB9" 
      
      shape_size = c()
      shape_size["ts"] = 1
      shape_size["fs"] = 1
      shape_size["ts_fs"] = 1
      shape_size["rest"] = 0
      
      ts_spec_value_max_abs = c()
      ts_spec_value_max_rel = c()
      ts_spec_tissue_max = c()
      fs_spec_value_max_abs = c()
      fs_spec_value_max_rel = c()
      fs_spec_tissue_max = c()
      cardinality = c()
      for (cl in 1:cluster_props$num_of_clusters) {
        # tissue specific
        tissue_spec_cl = tissue_spec[tissue_spec$CLUSTER == cl,]
        x = tissue_spec_cl$absolute_features
        maxi = max(x)
        if (length(which(x == maxi)) > 2) {
          print(sprintf("For cluster %s there are multiple %s (%s) with the highest occurance (%s absolute %s)", cl, tissue, paste(names(which(x == maxi)), collapse=", "), maxi, data_input$rna_class))
        }
        ts_spec_value_max_abs = append(ts_spec_value_max_abs, maxi)
        # realtive value is 0 if absolute == 0
        if (maxi != 0) {
          ts_spec_value_max_rel = append(ts_spec_value_max_rel, tissue_spec_cl[tissue_spec_cl$absolute_features == maxi,]$relative_features[1])
        } else {
          ts_spec_value_max_rel = append(ts_spec_value_max_rel, 0)
        }
        ts_spec_tissue_max = append(ts_spec_tissue_max, tissue_spec$feature[which(x == maxi)[1]])
        
        # feature specific
        feature_spec_cl = feature_spec[feature_spec$CLUSTER == cl,]
        x = feature_spec_cl$absolute_features
        maxi = max(x)
        if (length(which(x == maxi)) > 2) {
          print(sprintf("For cluster %s there are multiple features (%s) with the highest occurance (%s absolute %s)", cl, paste(feature_spec_cl$feature[which(x == maxi)], collapse=", "), maxi, data_input$rna_class))
        }
        fs_spec_value_max_abs = append(fs_spec_value_max_abs, maxi)
        if (maxi != 0) {
          fs_spec_value_max_rel = append(fs_spec_value_max_rel, feature_spec_cl[feature_spec_cl$absolute_features == maxi,]$relative_features[1])
        } else {
          fs_spec_value_max_rel = append(fs_spec_value_max_rel, 0)
        }
        fs_spec_tissue_max = append(fs_spec_tissue_max, feature_spec$feature[which(x == maxi)[1]])
        
        # cardinality of cluster for membership threshold
        tmp = membership[membership$CLUSTER == cl,]
        cardinality = append(cardinality, dim(tmp[tmp$MEM.SHIP >= membership_plot_th,])[1])
      }
      
      plot_df = data.frame(CLUSTER = seq(1, cluster_props$num_of_clusters, 1), 
                           ts_spec_value_max_abs = ts_spec_value_max_abs,
                           ts_spec_value_max_rel = ts_spec_value_max_rel,
                           ts_spec_tissue_max = ts_spec_tissue_max,
                           fs_spec_value_max_abs = fs_spec_value_max_abs,
                           fs_spec_value_max_rel = fs_spec_value_max_rel,
                           fs_spec_tissue_max = fs_spec_tissue_max,
                           cardinality = cardinality)
      
      # make it a datatable to make := work
      plot_df = as.data.table(plot_df)
      
      #plot_df$colour = colors["rest"]
      #plot_df[ts_spec_value_rel >= tissue_th, colour:=colors["ts"]]
      #plot_df[fs_spec_value >= feature_th, colour:=colors["fs"]]
      #plot_df[(ts_spec_value_rel >= tissue_th) & (fs_spec_value >= feature_th), colour:=colors["ts_fs"]]
      
      plot_df$colour = colors[["rest"]]
      tmp = ts_spec_value_max_rel >= tissue_th
      plot_df[tmp, ]$colour = colors[["ts"]][plot_df[tmp, ]$ts_spec_tissue_max]
      tmp = fs_spec_value_max_abs >= feature_th
      plot_df[tmp, ]$colour = colors[["fs"]]
      tmp = (ts_spec_value_max_rel >= tissue_th) & (fs_spec_value_max_abs >= feature_th)
      plot_df[tmp, ]$colour = colors[["ts_fs"]][plot_df[tmp, ]$ts_spec_tissue]
      
      plot_df$shape_size = shape_size["rest"]
      plot_df[ts_spec_value_max_rel >= tissue_th, ]$shape_size = shape_size["ts"]
      plot_df[fs_spec_value_max_abs >= feature_th, ]$shape_size =shape_size["fs"]
      plot_df[(ts_spec_value_max_rel >= tissue_th) & (fs_spec_value_max_abs >= feature_th), ]$shape_size = shape_size["ts_fs"]
      plot_df$shape_size_plot = plot_df$shape_size * plot_df$cardinality
      
      plot_df$CLUSTER_name = plot_df$CLUSTER * plot_df$shape_size
      plot_df$CLUSTER_name = gsub(0, "", plot_df$CLUSTER_name)
      
      x_axis_name = sprintf("Highest %s occ. (%%)", xticks_names$categories[[tissue]])
      y_axis_name = sprintf("Highest %s occ. (count)", data_input$rna_class)
      
      create_overview_scatter_plot(plot_df, tissue_th, feature_th, x_axis_name, y_axis_name, file_name, data_input, plots_props, output_folder_fig_scatter)
        
      plot_df_filtered = plot_df[plot_df$shape_size == 1,]
      
      tissue_spec_filtered = tissue_spec[tissue_spec$CLUSTER %in% plot_df_filtered$CLUSTER, ]
      feature_spec_filtered = feature_spec[feature_spec$CLUSTER %in% plot_df_filtered$CLUSTER, ]
      
      x_axis_name_left = sprintf("%s occ. (%%)", xticks_names$categories[[sprintf("%s_capital",tissue)]])
      x_axis_name_right = sprintf("%s occ. (count)", data_input$rna_class)
      file_name = "cluster_analysis"
      create_scatter_plot(plot_df_filtered, tissue_spec_filtered, feature_spec_filtered, tissue, tissue_th, feature_th, x_axis_name_left, x_axis_name_right, colours, file_name, plots_props, output_folder_fig_scatter)
    }
  }
}


#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])

