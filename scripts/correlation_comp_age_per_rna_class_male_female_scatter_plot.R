suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggridges))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(svglite))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(openxlsx))
#library(corrplot)
suppressPackageStartupMessages(library(Hmisc))

#snakemake = readRDS("snakemake_correlation_comp_age_per_rna_class_male_female_scatter_plot.rds")

set.seed(snakemake@params$parameters_props$set_seed)

print("correlation_comp_age_per_rna_class_male_female_scatter_plot")


#---------------------------------- Functions ----------------------------------


#---------------------------------- Load data ---------------------------------- 
data_input = snakemake@params$data_input
props = snakemake@params$prop_attribute
method = snakemake@params$corr_methods
corr_ths = snakemake@params$corr_th
adjustment = snakemake@params$p_val_adj
adj_p_value = snakemake@params$adj_p_value
rna_class = snakemake@params$rna_class
ID = data_input$identifier_column
male_folder_path = data_input$data_sub_set_male
female_folder_path = data_input$data_sub_set_male
xticks_names = snakemake@params$xticks_names
plots_props = snakemake@params$plots_props
colors = snakemake@params$colors
plot_size = snakemake@params$plot_size
results_folder = snakemake@params$results_folder

# save rdata
#saveRDS(snakemake, file = "snakemake_correlation_comp_age_per_rna_class_male_female_scatter_plot.rds")
#stop()

annot = fread(sprintf("%s/annotation_%s_%s.csv", data_input$annotation_data_folder, data_input$data_sub_set, data_input$feature_filtering), sep='\t', colClasses=c(ID="character")) 


#------------------------------------ Script ----------------------------------- 
for (m in method) {
  #for (prop in props) {
  split_group = props[1]
  corr_variable = props[2]
  
  for (sex_key in c("male", "female")) {
    output_folder_mele = sprintf("%s/%s_%s/results_%s/figures/corr_plots/age/ridgeline_plot", results_folder, male_folder_path, data_input$detection_rate, data_input$data_sub_set)

    # read tables
    fread(correlation_df, sprintf("%s/corr_method=%s_%ss_with_%s_per_%s.csv", output_folder_tab, m, rna_class, corr_variable, split_group), sep = "\t", row.names = TRUE)
    fread(corr_pvalues_df, sprintf("%s/corr_method=%s_adj_pvalues_%ss_with_%s_per_%s.csv", output_folder_tab, m, rna_class, corr_variable, split_group), sep = "\t", row.names = TRUE)

    corr_values_per_split_group = correlation_df[,2:ncol(correlation_df)]
    
    corr_values = unname(unlist(correlation_df[,2:ncol(correlation_df)]))

    corr_values_melted = melt(corr_values)
    colnames(corr_values_melted) = c("value", "variable")
    
    mask = (names(unlist(colors$rna_classes)) %in% unique(corr_values_melted$variable))
    colours_rna_class_filtered = unlist(colors$rna_classes)[mask]
    
    median_values = aggregate(value ~ variable, data = corr_values_melted, FUN = median, na.rm = TRUE)
    median_values = median_values[order(-median_values$value), ]  # Sort in descending order
    
    corr_values_melted$variable_fac = factor(corr_values_melted$variable, levels = median_values$variable)
    
    #corr_values_melted$variable_fac = factor(corr_values_melted$variable, levels = rev(names(colors$rna_classes)))
    
    if (m == "spearman") {
      x_axis_title = sprintf("Spearman corr. (%s) per %s", corr_variable, xticks_names$categories[[split_group]])
    } else {
      x_axis_title = sprintf("Pearson corr. (%s) per %s", corr_variable, xticks_names$categories[[split_group]])
    }
    
    ridgeline_plot = ggplot(corr_values_melted, aes(x = value, y = variable_fac, fill = variable_fac)) +
      #geom_density_ridges(scale = 1, alpha = 0.7, size = 0.3, show.legend = FALSE) +  # Ridgeline plot
      geom_density_ridges_gradient(quantile_lines = TRUE, quantiles = 2, size = 0.15, show.legend = FALSE) + 
      geom_vline(xintercept = c(-corr_th, 0, corr_th), linetype = "dashed", color = "grey", size = 0.3) +  # Threshold lines
      scale_y_discrete(expand = expansion(mult = c(0, 0))) +  # Space at the bottom for x-axis
      #scale_x_continuous(limits = c(-1,1), breaks = c(-1, -0.5, 0, 0.5, 1)) +  # Space at the bottom for x-axis
      scale_fill_manual(values = colours_rna_class_filtered) +
      xlab(x_axis_title) +
      theme_classic() +
      theme(
        #legend.position = "none", 
        #legend.title = element_blank(), #element_text(family = plots_props$font_family, size = plots_props$font_size),
        #legend.text = element_text(family = plots_props$font_family, size = plots_props$font_size, margin = margin(0, 0, 0, 0.3, "pt")),
        text = element_text(family = plots_props$font_family, size = plots_props$font_size),
        #axis.text = element_text(angle = 0, family = plots_props$font_family, size = plots_props$font_size),
        #axis.title = element_text(family = plots_props$font_family, size = plots_props$font_size),
        axis.text.y = element_text(family = plots_props$font_family, size = plots_props$font_size, vjust = -0.5, margin = margin(r = -10)),  # Labels for matrix names
        axis.title.y = element_blank(),  # Remove y-axis title
        axis.text.x = element_text(family = plots_props$font_family, size = plots_props$font_size, margin = margin(t = -2)),
        axis.title.x = element_text(family = plots_props$font_family, size = plots_props$font_size, hjust = 1),  # Remove x-axis title
        axis.ticks.y = element_blank(), 
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        #strip.text.y.left = element_text(family = plots_props$font_family, angle = 0, size=plots_props$font_size),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        #panel.spacing = unit(0, "pt"),
        #plot.background = element_rect(fill = "transparent"),
      ) +
      coord_cartesian(clip = "off")  # Prevents cutting off labels
    #ridgeline_plot
  
    # save
    width = plot_size$image_width
    height = plot_size$image_height
    ggsave(sprintf("%s/corr_method=%s_corr_th=%s_features_with_%s_per_%s_%sx%s.png", output_folder_fig, m, corr_th, corr_variable, split_group, height, width), ridgeline_plot, dpi=plots_props$dpi, width=width, height=height, units = plots_props$image_units)
    ggsave(sprintf("%s/corr_method=%s_corr_th=%s_features_with_%s_per_%s_%sx%s.svg", output_folder_fig, m, corr_th, corr_variable, split_group, height, width), ridgeline_plot, width = width, height = height, unit = plots_props$image_units, dpi = plots_props$dpi)

    aggregated_corr_combined = do.call(rbind, aggregated_corr)
    fwrite(aggregated_corr_combined, sprintf("%s/number_of_corr_corr_method=%s_corr_th=%s_%ss_with_%s_per_%s.csv", output_folder_tab, m, corr_th, rna_cl, corr_variable, split_group), sep = "\t", row.names = FALSE)
    write.xlsx(aggregated_corr_combined, sprintf("%s/number_of_corr_corr_method=%s_corr_th=%s_%ss_with_%s_per_%s.xlsx", output_folder_tab, m, corr_th, rna_cl, corr_variable, split_group), colNames = TRUE, rowNames = FALSE, append = FALSE)
  }
  print(sprintf("Done: %s", m))
}


#---------------------------------- Save step ---------------------------------- 
# creates an empty file in the step folder
file.create(snakemake@output[[1]])
