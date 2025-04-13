suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(svglite))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(stringr))


#snakemake = readRDS("snakemake_mfuzz_ts_analysis_membership_box_plot.rds")

set.seed(snakemake@params$parameters_props$set_seed)

print("mfuzz_ts_analysis_membership_box_plot")


#---------------------------------- Load data ---------------------------------- 
data_input = snakemake@params$data_input
cluster_props = snakemake@params$cluster_props
thresholds = snakemake@params$thresholds
interesting_cluster = snakemake@params$interesting_cluster
colours = snakemake@params$colors
plots_props = snakemake@params$plots_props
results_folder = snakemake@params$results_folder

# save rdata
#saveRDS(snakemake, file = "snakemake_mfuzz_ts_analysis_membership_box_plot.rds")

input_folder_path = sprintf("%s/%s_%s/results_%s/matrices/mfuzz_clustering/k=%s", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, cluster_props$num_of_clusters)

membership = fread(sprintf("%s/cluster_members_k=%s.csv", snakemake@input$cluster_results_folder_path, cluster_props$num_of_clusters), sep=',', header=T) 
#cluster_result = readRDS(file.path(input_folder_path, sprintf("cluster_result_k=%s.rds", cluster_props$num_of_clusters)))

output_folder_fig = sprintf("%s/%s_%s/results_%s/figures/mfuzz_clustering/k=%s/all/box_plots", results_folder, data_input$rna_class, data_input$detection_rate, data_input$data_sub_set, cluster_props$num_of_clusters)
dir.create(output_folder_fig, recursive=TRUE)


# --------------------------------- Script -------------------------------------

#plot_df = melt(cluster_result$membership)
plot_df = membership

plot_df$tissue = str_split_fixed(plot_df$NAME, "_", 2)[,2]

plot_df$interesting = rep(0, dim(plot_df)[1])
tmp = plot_df$CLUSTER %in% interesting_cluster
plot_df[tmp,]$interesting = 1
#plot_df$interesting = factor(plot_df$interesting)

x_axis_name = "Cluster"
y_axis_name = "Membership"
 
scatter_plot = ggplot(plot_df, aes(x=CLUSTER, y=MEM.SHIP, colour = factor(interesting))) + 
  geom_point(size = 0.1) +
  scale_color_manual(values=c("black","red")) + 
  #geom_hline(yintercept = 0.3, color="grey", linetype="dashed", size = 0.5) + 
  #geom_hline(yintercept = 0.25, color="grey", linetype="dashed", size = 0.5) + 
  geom_hline(yintercept = thresholds$membership, color="blue", linetype="dashed", size = 0.5) + 
  #geom_hline(yintercept = 0.1, color="grey", linetype="dashed", size = 0.5) + 
  xlab(x_axis_name) + ylab(y_axis_name) +
  theme_classic() +
  theme(legend.position = "none", 
        #plot.background = element_rect(fill="transparent"),
        #panel.background = element_rect(fill="transparent"),
        text = element_text(family = plots_props$font_family, size = plots_props$font_size), 
        axis.text = element_text(family = plots_props$font_family, size = plots_props$font_size), 
        axis.title = element_text(family = plots_props$font_family, size = plots_props$font_size), 
        plot.title = element_text(family = plots_props$font_family, size = plots_props$font_size), 
        legend.text = element_text(family = plots_props$font_family, size = plots_props$font_size),
        legend.title = element_text(family = plots_props$font_family, size = plots_props$font_size, face = "bold"),
        plot.margin = margin(b=-0.4, unit = "cm")) #+ 
#scatter_plot

bar_plot = ggplot(plot_df, aes(x = CLUSTER, y = MEM.SHIP, group = CLUSTER, color = factor(interesting))) +
  geom_boxplot(outlier.size=0.5) +
  scale_color_manual(values=c("black", "red")) +
  #geom_boxplot(width = 1.75, color = "black", outlier.colour="grey", outlier.shape=16, outlier.size=1, notch=FALSE, fatten=TRUE) + 
  #geom_boxplot(color="black", width = 1.75, fatten=FALSE, coef = 0, outlier.shape = NA, outlier.alpha = 0) +
  geom_hline(yintercept = thresholds$membership, color="blue", linetype="dashed", size = 0.5) +
  #geom_jitter(shape=16, position=position_jitter(0.2), alpha =0.1) +
  xlab(x_axis_name) + ylab(y_axis_name) +
  theme_classic() +
  theme(
    #plot.margin = unit(c(0,0.1,0,0.1), plots_props$image_units), #top, right, bottom, left
    plot.background = element_rect(fill='transparent', color=NA),
    legend.position="none", 
    text = element_text(family = plots_props$font_family, size = plots_props$font_size),
    axis.text = element_text(family = plots_props$font_family, size = plots_props$font_size),
    axis.title = element_text(family = plots_props$font_family, size = plots_props$font_size),
    plot.title = element_text(family = plots_props$font_family, size = plots_props$font_size_header),
    legend.text = element_text(family = plots_props$font_family, size = plots_props$font_size),
    legend.title = element_text(family = plots_props$font_family, size = plots_props$font_size)
  )
#bar_plot

ggsave(sprintf("%s/distribution_membership.png", output_folder_fig), bar_plot, dpi=plots_props$dpi, width = 2*plots_props$image_width, height = plots_props$image_height, units = plots_props$image_units)
ggsave(sprintf("%s/distribution_membership.svg", output_folder_fig), bar_plot, width = 2*plots_props$image_width, height = plots_props$image_height, units =plots_props$image_units, limitsize = FALSE, scale = 1)


#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])

